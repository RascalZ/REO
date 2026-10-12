"""Bind formal model/statistical protocol, code and cohort-v2 scientific inputs.

No original drug-response tables, CTRP response/index or expression-array bytes
are opened. Inherited data descriptors remain bound through the cohort-v2 lock;
the execution runner verifies actual development inputs before fitting.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def local(path):
    path = Path(path)
    path = path.resolve() if path.is_absolute() else (ROOT / path).resolve()
    require(path.is_relative_to(ROOT), "Formal scientific dependency must remain repository-relative")
    return path


def scientific_contract(model_path, statistics_path, protocol_path):
    model, statistics = load(model_path), load(statistics_path)
    panel = pd.read_csv(ROOT / "results/cohort_v2/main_drugs.csv", dtype=str, keep_default_na=False)
    manifest = pd.read_csv(ROOT / "results/cohort_v2/development_sample_manifest.csv", dtype=str, keep_default_na=False)
    expected_drugs = panel.gdsc_drug_id.astype(int).tolist()
    expected_models = ["mean", "tissue", "expression", "expression_matched", "reo", "expression_rf", "reo_rf"]
    require(len(expected_drugs) == len(set(expected_drugs)) == 8 and model["drug_ids"] == statistics["drug_ids_in_frozen_rank"] == expected_drugs, "Fixed eight-drug panel/rank differs across scientific configurations")
    require(model["models"] == statistics["models"] == expected_models, "Formal model families changed or duplicated")
    require(model["cohort_version"] == statistics["cohort_version"] == "v2", "Formal analysis must bind cancer cohort v2")
    for key, number in (("outer_folds", 5), ("outer_repeats", 3), ("inner_folds", 3)):
        require(model[key] == statistics[key] == number, "Formal nested-CV dimensions differ: " + key)
    require(len(model["repeat_seeds"]) == len(set(model["repeat_seeds"])) == 3, "Repeats require three distinct frozen seeds")
    require(set(manifest.role) == {"development"} and not manifest.duplicated(["gdsc_drug_id", "canonical_model_id"]).any(), "Frozen formal development manifest has duplicate/nondevelopment rows")
    require(model["elastic_net"]["regularization_mode"] == "alpha_max_relative", "Formal relative regularization path not declared")
    require(all(0 < a <= 1 for a in model["elastic_net"]["alpha_ratios"]) and 1.0 in model["elastic_net"]["alpha_ratios"], "Relative path lacks valid theoretical null endpoint")
    require(model["random_forest"]["n_jobs"] == 1 and model["resources"]["blas_threads"] == 1 and model["resources"]["parallel_jobs"] in (1, 2), "Frozen worker/thread budget invalid")
    require(statistics["primary_comparison"] == {"candidate": "reo", "reference": "expression"}, "Primary representation comparison changed")
    require(statistics["secondary_comparison"] == {"candidate": "reo_rf", "reference": "expression_rf"}, "Secondary RF comparison changed")
    require(statistics["multiplicity"]["primary_family_size"] == statistics["multiplicity"]["secondary_family_size"] == 9, "Family size must retain eight fixed drugs plus macro")
    require(statistics["bootstrap"]["draws"] == 10000 and statistics["bootstrap"]["minimum_valid_draw_fraction"] == 0.99, "Bootstrap precision/validity contract changed")
    require(statistics["no_confirmatory_p_values_from_bootstrap"] and statistics["no_equivalence_claim_from_non_significance"], "Unsupported inference claim permitted")
    require(protocol_path.exists() and len(protocol_path.read_text(encoding="utf-8")) > 1000, "Formal analysis protocol missing or incomplete")
    candidate_counts = {"mean": 1, "tissue": len(model["ridge_alphas"]),
                        "expression": len(model["elastic_net"]["alpha_ratios"]) * len(model["elastic_net"]["l1_ratios"]),
                        "expression_matched": len(model["elastic_net"]["alpha_ratios"]) * len(model["elastic_net"]["l1_ratios"]),
                        "reo": len(model["elastic_net"]["alpha_ratios"]) * len(model["elastic_net"]["l1_ratios"]),
                        "expression_rf": len(model["random_forest"]["max_features"]) * len(model["random_forest"]["min_samples_leaf"]),
                        "reo_rf": len(model["random_forest"]["max_features"]) * len(model["random_forest"]["min_samples_leaf"])}
    partitions = len(expected_drugs) * model["outer_repeats"] * model["outer_folds"]
    fits = partitions * sum(1 if kind == "mean" else count * model["inner_folds"] + 1 for kind, count in candidate_counts.items())
    tasks = partitions * len(expected_models)
    predictions = len(manifest) * len(expected_models) * model["outer_repeats"]
    require((tasks, fits, predictions) == (840, 18120, 104307), "Formal job/fit/OOF budget differs from declared full analysis")
    cohort_path = ROOT / "config/cohort_freeze_v2.json"
    cohort_lock_path = ROOT / "config/cohort_freeze_v2.lock.json"
    cohort, cohort_lock = load(cohort_path), load(cohort_lock_path)
    require(sha(cohort_path) == cohort_lock["freeze_sha256"], "Cohort v2 differs from its validated lock")
    require(sha(ROOT / "config/cohort_freeze_v1.json") == cohort["parent_cohort_freeze_sha256"], "Historical cohort-v1 linkage changed")
    require(sha(ROOT / "config/cohort_freeze_v1.lock.json") == cohort["parent_lock_sha256"], "Historical cohort-v1 lock linkage changed")
    recorded = {r["path"].replace("\\", "/"): r["sha256"] for r in cohort["files"]}
    for path in (ROOT / "results/cohort_v2/main_drugs.csv", ROOT / "results/cohort_v2/development_sample_manifest.csv"):
        require(sha(path) == recorded[path.relative_to(ROOT).as_posix()], "Frozen development identity/panel input changed")
    return {"tasks": tasks, "fits": fits, "oof_predictions": predictions,
            "drug_model_observations": len(manifest), "unique_development_groups": manifest.biological_group_id.nunique(),
            "candidate_counts": candidate_counts}, cohort_lock


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-config", default="config/formal_model_v1.json")
    parser.add_argument("--statistics-config", default="config/formal_statistics_v1.json")
    parser.add_argument("--protocol", default="docs/FORMAL_ANALYSIS_PROTOCOL.md")
    parser.add_argument("--output-lock", default="config/formal_analysis_v1.lock.json")
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--lock", action="store_true", help="Create the immutable pre-run formal scientific lock")
    action.add_argument("--verify-only", action="store_true", help="Verify an existing immutable lock without rewriting it")
    args = parser.parse_args()
    model_path, statistics_path, protocol_path, lock_path = map(local, (args.model_config, args.statistics_config, args.protocol, args.output_lock))
    contract, cohort_lock = scientific_contract(model_path, statistics_path, protocol_path)
    execution = [model_path, statistics_path, protocol_path, ROOT / "src/run_formal_evaluation.py",
                 ROOT / "src/analyze_formal_results.py", ROOT / "src/run_development_pilot.py", ROOT / "src/project_paths.py",
                 *sorted((ROOT / "src/modeling").glob("*.py"))]
    scientific = execution + [ROOT / "config/cohort_freeze_v2.json", ROOT / "config/cohort_freeze_v2.lock.json",
                             ROOT / "results/cohort_v2/main_drugs.csv", ROOT / "results/cohort_v2/development_sample_manifest.csv",
                             ROOT / "data_manifest/original_source_sha256.csv"]
    files = {path.relative_to(ROOT).as_posix(): sha(path) for path in scientific}
    require(not any(name.startswith(("reports/", "runs/")) or "RESEARCH_PLAN" in name for name in files), "Mutable research/report output entered formal freeze")
    if args.verify_only or lock_path.exists():
        require(lock_path.exists(), "Formal lock missing")
        previous = load(lock_path)
        require(previous["files"] == files and previous["contract"] == contract and previous["cohort_freeze_sha256"] == cohort_lock["freeze_sha256"], "Scientific formal lock changed; use a new version/run rather than overwrite")
        print(json.dumps({"status": "verified", "lock": lock_path.relative_to(ROOT).as_posix(), "lock_sha256": sha(lock_path), "contract": contract}, ensure_ascii=False))
        return 0
    if not args.lock:
        print(json.dumps({"status": "ready_for_lock", "contract": contract}, ensure_ascii=False))
        return 0
    result = {"schema_version": "1.0", "status": "frozen_pre_formal_analysis", "created_utc": datetime.now(timezone.utc).isoformat(),
              "model_config": model_path.relative_to(ROOT).as_posix(), "statistics_config": statistics_path.relative_to(ROOT).as_posix(),
              "protocol": protocol_path.relative_to(ROOT).as_posix(), "cohort_version": "v2", "cohort_freeze_sha256": cohort_lock["freeze_sha256"],
              "cohort_lock_sha256": sha(ROOT / "config/cohort_freeze_v2.lock.json"), "files": files,
              "execution_files": [path.relative_to(ROOT).as_posix() for path in execution], "contract": contract,
              "external_response_accessed": False, "original_response_or_expression_array_opened_by_freezer": False,
              "inherited_scientific_inputs": "Bound through the immutable cohort-v2 lock and verified by the development runner before use.",
              "prior_pilot_results_seen": True, "mutable_plan_and_reports_excluded": True}
    with lock_path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": "locked", "lock": lock_path.relative_to(ROOT).as_posix(), "lock_sha256": sha(lock_path), "contract": contract}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(json.dumps({"status": "failed", "error_type": type(exc).__name__, "error": str(exc).replace(str(ROOT), "<repository>")}, ensure_ascii=False))
        sys.exit(1)
