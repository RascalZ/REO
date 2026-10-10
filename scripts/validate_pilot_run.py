"""Independently audit pilot artifacts; never open source response/expression data.

Usage: .venv/Scripts/python.exe scripts/validate_pilot_run.py --run-dir runs/NAME
This validator uses run artifacts, the frozen development identity manifest,
cohort lock metadata, and Git source blobs. It does not import the model runner.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def artifact_label(path):
    path = Path(path).resolve()
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.name


def validate(run, revision):
    run = Path(run).resolve()
    meta = read_json(run / "run_manifest.json")
    config = meta["configuration"]
    drugs = {str(d) for d in config["drug_ids"]}
    models = set(config["models"])
    outer_count, inner_count = config["outer_folds"], config["inner_folds"]
    require(config["outer_repeats"] == 1, "This pilot validator requires one outer repeat")
    require(meta["status"] == "completed", "Pilot run is not completed")

    freeze_path = ROOT / "config/cohort_freeze_v1.json"
    require(sha(freeze_path) == meta["cohort_lock"]["freeze_sha256"], "Current cohort description differs from run lock")
    freeze = read_json(freeze_path)
    frozen_development = next(r for r in freeze["files"]
                              if str(r["path"]).replace("\\", "/") == "results/development_sample_manifest.csv")
    development_path = ROOT / "results/development_sample_manifest.csv"
    require(sha(development_path) == frozen_development["sha256"], "Development identity assignments changed")
    # No response values are present in this identity-only phase-1 manifest.
    development = pd.read_csv(development_path, dtype=str, keep_default_na=False)
    require(set(development.role) == {"development"}, "Unexpected role in development manifest")
    global_identity = development.drop_duplicates("canonical_model_id").set_index("canonical_model_id")
    require(not development.groupby("canonical_model_id").biological_group_id.nunique().gt(1).any(), "Inconsistent biological-group identities")
    require(not development.groupby("canonical_model_id").expression_row_index.nunique().gt(1).any(), "Inconsistent expression-row identities")
    selected = development[development.gdsc_drug_id.isin(drugs)].copy()
    require(set(selected.gdsc_drug_id) == drugs, "Selected drugs missing from development manifest")
    require(not selected.duplicated(["gdsc_drug_id", "canonical_model_id"]).any(), "Duplicated development drug/model identities")

    folds = pd.read_csv(run / "global_fold_manifest.csv", dtype={"canonical_model_id": str, "biological_group_id": str, "tissue": str})
    require(not folds.duplicated(["outer_fold", "inner_fold", "role", "canonical_model_id"]).any(), "Duplicated fold memberships")
    require(set(folds.outer_fold) == set(range(outer_count)), "Unexpected outer fold indices")
    for row in folds.itertuples(index=False):
        require(row.canonical_model_id in global_identity.index, "Fold identity outside development pool")
        identity = global_identity.loc[row.canonical_model_id]
        require(row.biological_group_id == identity.biological_group_id and row.tissue == identity.tissue, "Fold identity metadata mismatch")

    membership = {}
    for (outer, inner, role), frame in folds.groupby(["outer_fold", "inner_fold", "role"]):
        membership[(int(outer), int(inner), role)] = set(frame.canonical_model_id)
    expected_membership_keys = {(outer, -1, role) for outer in range(outer_count) for role in ("train", "test")}
    expected_membership_keys |= {(outer, inner, role) for outer in range(outer_count)
                                 for inner in range(inner_count) for role in ("train", "validation")}
    require(set(membership) == expected_membership_keys, "Unexpected fold roles or indices")
    all_ids = set(global_identity.index)
    outer_test_coverage = Counter()
    group_for = global_identity.biological_group_id.to_dict()
    def groups_for(ids):
        return {group_for[i] for i in ids}
    for outer in range(outer_count):
        train = membership[(outer, -1, "train")]
        test = membership[(outer, -1, "test")]
        require(not train & test and train | test == all_ids, "Incomplete or overlapping outer split")
        require(not groups_for(train) & groups_for(test), "Biological group crosses outer train/test")
        outer_test_coverage.update(test)
        inner_validation_coverage = Counter()
        for inner in range(inner_count):
            inner_train = membership[(outer, inner, "train")]
            validation = membership[(outer, inner, "validation")]
            require(not inner_train & validation and inner_train | validation == train, "Incomplete or overlapping inner split")
            require(not groups_for(inner_train) & groups_for(validation), "Biological group crosses inner train/validation")
            require(not (inner_train | validation) & test, "Outer test appears in inner folds")
            inner_validation_coverage.update(validation)
        require(set(inner_validation_coverage) == train and set(inner_validation_coverage.values()) == {1}, "Inner validation coverage is not exactly once")
    require(set(outer_test_coverage) == all_ids and set(outer_test_coverage.values()) == {1}, "Outer test coverage is not exactly once")

    predictions = pd.read_csv(run / "oof_predictions.csv", dtype={"gdsc_drug_id": str, "canonical_model_id": str, "biological_group_id": str})
    require(set(predictions.gdsc_drug_id) == drugs and set(predictions.model) == models, "Unexpected prediction drugs/models")
    require(not predictions.duplicated(["gdsc_drug_id", "model", "canonical_model_id"]).any(), "More than one OOF prediction per drug/model/identity")
    require(np.isfinite(predictions[["observed_auc", "predicted_auc"]].to_numpy()).all(), "Non-finite OOF values")
    require(predictions.observed_auc.between(0, 1).all(), "OOF observed AUC outside [0,1]")
    require(not predictions.groupby(["gdsc_drug_id", "canonical_model_id"]).observed_auc.nunique().gt(1).any(), "Models used different observed responses")
    selected_lookup = selected.set_index(["gdsc_drug_id", "canonical_model_id"])
    counts = {}
    for (drug, model), frame in predictions.groupby(["gdsc_drug_id", "model"]):
        expected_ids = set(selected[selected.gdsc_drug_id.eq(drug)].canonical_model_id)
        require(set(frame.canonical_model_id) == expected_ids and len(frame) == len(expected_ids), "Missing or extra OOF identities")
        counts.setdefault(drug, {})[model] = len(frame)
        for row in frame.itertuples(index=False):
            identity = selected_lookup.loc[(drug, row.canonical_model_id)]
            require(row.biological_group_id == identity.biological_group_id, "Prediction biological group mismatch")
            require(int(row.expression_row_index) == int(identity.expression_row_index) and row.tissue == identity.tissue, "Prediction expression identity mismatch")
            require(row.canonical_model_id in membership[(row.outer_fold, -1, "test")], "OOF prediction does not belong to outer test")
    expected_total = sum(len(selected[selected.gdsc_drug_id.eq(drug)]) for drug in drugs) * len(models)
    require(len(predictions) == expected_total == meta["n_oof_predictions"], "OOF total count mismatch")

    observed = predictions[predictions.model.eq("mean")].set_index(["gdsc_drug_id", "expression_row_index"]).observed_auc
    events = [json.loads(line) for line in (run / "fit_events.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    choices = read_json(run / "chosen_parameters.json")
    choices_lookup = {(str(c["gdsc_drug_id"]), c["model"], c["outer_fold"]): c for c in choices}
    require(len(choices_lookup) == len(choices) == len(drugs) * len(models) * outer_count, "Chosen-parameter records incomplete/duplicated")
    def grid(model):
        if model == "mean":
            return [{}]
        if model == "tissue":
            return [{"alpha": a} for a in config["ridge_alphas"]]
        return [{"alpha": a, "l1_ratio": r} for a, r in itertools.product(config["elastic_net"]["alphas"], config["elastic_net"]["l1_ratios"])]
    event_keys = Counter()
    feature_fit_count = 0
    for event in events:
        drug, model, outer = str(event["gdsc_drug_id"]), event["model"], int(event["outer_fold"])
        require(drug in drugs and model in models and outer in range(outer_count), "Unexpected fit context")
        require(event["model_kind"] == model, "Fit model kind differs from context")
        drug_rows = selected[selected.gdsc_drug_id.eq(drug)]
        drug_ids = set(drug_rows.canonical_model_id)
        id_for_expression = dict(zip(drug_rows.expression_row_index.astype(int), drug_rows.canonical_model_id))
        fit_samples = event["fit_sample_ids"]
        require(len(fit_samples) == len(set(fit_samples)), "Duplicated fit sample IDs")
        require(set(fit_samples) <= set(id_for_expression), "Fit sample outside drug development cohort")
        fitted_ids = {id_for_expression[i] for i in fit_samples}
        stage, candidate, inner = event["stage"], int(event["candidate"]), int(event["inner_fold"])
        key = (drug, model, outer, stage, candidate, inner)
        event_keys.update([key])
        if stage == "inner_training":
            require(model != "mean" and 0 <= inner < inner_count and 0 <= candidate < len(grid(model)), "Invalid inner-fit indices")
            expected_ids = membership[(outer, inner, "train")] & drug_ids
            forbidden_ids = (membership[(outer, inner, "validation")] | membership[(outer, -1, "test")]) & drug_ids
            require(event["params"] == grid(model)[candidate], "Candidate parameters do not match fixed grid")
        else:
            require(stage == "outer_training_refit" and candidate == -1 and inner == -1, "Unknown fitting stage")
            expected_ids = membership[(outer, -1, "train")] & drug_ids
            forbidden_ids = membership[(outer, -1, "test")] & drug_ids
            require(event["params"] == choices_lookup[(drug, model, outer)]["chosen_params"], "Outer refit parameters differ from inner choice")
        require(fitted_ids == expected_ids, "Actual fit membership differs from intended training partition")
        require(not fitted_ids & forbidden_ids and not groups_for(fitted_ids) & groups_for(forbidden_ids), "Fit contains current validation/test identity/group")
        responses = np.array([observed.loc[(drug, sample)] for sample in fit_samples], dtype=float)
        require(np.isclose(event["target_mean"], responses.mean(), atol=1e-12, rtol=0), "Target mean was not computed from actual training responses")
        expected_sd = responses.std() if responses.std() != 0 else 1.0
        require(np.isclose(event["target_sd"], expected_sd, atol=1e-12, rtol=0), "Target scale was not computed from actual training responses")
        if model in {"expression", "expression_matched", "reo"} and not event["fallback_reason"]:
            maximum_genes = config["feature"]["expression_genes"] if model == "expression" else config["feature"]["reo_genes"]
            require(0 < event["n_selected_genes"] <= maximum_genes and event["n_features"] > 0, "Feature builder was not recorded for this training fit")
            feature_fit_count += 1
    expected_fit_count = 0
    for drug, model, outer in itertools.product(drugs, models, range(outer_count)):
        require(event_keys[(drug, model, outer, "outer_training_refit", -1, -1)] == 1, "Outer training refit is not exactly once")
        expected_fit_count += 1
        if model != "mean":
            for candidate, inner in itertools.product(range(len(grid(model))), range(inner_count)):
                require(event_keys[(drug, model, outer, "inner_training", candidate, inner)] == 1, "Candidate/inner partition was not freshly fitted exactly once")
                expected_fit_count += 1
            choice = choices_lookup[(drug, model, outer)]
            scores = choice["inner_scores"]
            require(len(scores) == len(grid(model)), "Inner candidate scores incomplete")
            best = min(scores, key=lambda s: s["pooled_inner_rmse"])
            require(choice["chosen_params"] == best["params"] and choice["best_inner_rmse"] == best["pooled_inner_rmse"], "Chosen parameters do not minimize logged inner RMSE")
    require(len(events) == expected_fit_count == meta["n_fits"] and set(event_keys.values()) == {1}, "Fit total or unique-context count mismatch")

    summary = pd.read_csv(run / "summary_metrics.csv", dtype={"gdsc_drug_id": str})
    require(len(summary) == len(drugs) * len(models), "Summary metrics incomplete")
    for row in summary.itertuples(index=False):
        frame = predictions[predictions.gdsc_drug_id.eq(row.gdsc_drug_id) & predictions.model.eq(row.model)]
        y, pred = frame.observed_auc.to_numpy(), frame.predicted_auc.to_numpy()
        mse = np.mean((y - pred) ** 2)
        expected = {"mae": np.mean(np.abs(y - pred)), "rmse": np.sqrt(mse), "r2": 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)}
        require(row.n == len(frame), "Summary n differs from OOF records")
        for name, value in expected.items():
            require(np.isclose(getattr(row, name), value, atol=1e-12, rtol=0), "Summary metric differs from OOF records: " + name)
        if row.model == "mean":
            require(pd.isna(row.spearman) and row.spearman_unavailable_reason == "fold_constant_training_mean_not_an_individual_ranking_model", "Mean baseline given misleading pooled ranking score")
            for outer, part in frame.groupby("outer_fold"):
                require(part.predicted_auc.nunique() == 1, "Mean baseline varies within a test fold")
        else:
            require(np.isclose(row.spearman, spearmanr(y, pred).statistic, atol=1e-12, rtol=0), "Summary Spearman differs from OOF records")

    source_checks = []
    full_revision = subprocess.check_output(["git", "rev-parse", revision], cwd=ROOT, text=True).strip()
    for name, recorded_sha in meta["source_hashes"].items():
        git_path = name.replace("\\", "/")
        blob = subprocess.check_output(["git", "show", f"{revision}:{git_path}"], cwd=ROOT)
        revision_sha = hashlib.sha256(blob).hexdigest()
        current_sha = sha(ROOT / git_path)
        require(revision_sha == recorded_sha, "Recorded run source does not match immutable Git revision: " + git_path)
        source_checks.append({"path": git_path, "recorded_sha256": recorded_sha, "revision_matches": True,
                              "current_working_copy_matches": current_sha == recorded_sha})
    require(meta["external_response_accessed"] is False and not meta["blocked_external_access_attempts"], "Run reports external access or blocked attempts")
    require(not any("ctrp" in path.lower() for path in meta["input_files_accessed"]), "Access log includes an external-response path")

    return {
        "status": "passed", "validated_utc": datetime.now(timezone.utc).isoformat(),
        "run_name": config["run_name"], "run_directory": artifact_label(run),
        "validation_scope": "Recorded actual fit memberships, fold/group boundaries, OOF completeness, train-only target scaling, fixed-grid choice, OOF metric reconstruction and immutable source provenance",
        "inputs_read": ["run artifacts", "frozen development identity manifest", "cohort-lock metadata", "Git source blobs"],
        "original_GDSC_or_CTRP_response_read_by_validator": False,
        "source_revision": full_revision, "source_checks": source_checks,
        "counts": {"development_pool_identities": len(global_identity), "drug_model_oof_counts": counts,
                   "actual_fits": len(events), "expected_fits": expected_fit_count,
                   "feature_builder_fits_with_training_membership_verified": feature_fit_count,
                   "oof_predictions": len(predictions), "outer_folds": outer_count, "inner_folds": inner_count},
        "checks": {"cohort_lock_and_development_identity_hash": True, "all_group_and_fold_boundaries": True,
                   "actual_inner_fit_membership_and_exactly_once_refits": True, "outer_test_not_used_for_fitting": True,
                   "each_identity_exactly_one_oof_prediction_per_drug_model": True,
                   "target_scaling_from_actual_training_responses": True, "fixed_grid_and_logged_inner_rmse_selection": True,
                   "oof_metrics_reconstructed": True, "mean_baseline_ranking_not_interpreted": True,
                   "recorded_sources_match_immutable_run_revision": True, "run_external_access_log_clear": True},
        "limitations": ["Fit-boundary verification depends on actual-fit audit records; synthetic tests independently inspect builder inputs and perturb held-out labels/distributions.",
                        "Inner OOF predictions were not exported; parameter choice is checked against logged candidate scores, while pooled outer metrics are reconstructed from individual OOF records.",
                        "External access checks inspect the run's access log and guard tests; the validator does not reopen external responses."]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--source-revision", default="629b8f5", help="Immutable revision whose bytes were used by this pilot")
    args = parser.parse_args()
    run = args.run_dir.resolve()
    try:
        report = validate(run, args.source_revision)
    except Exception as exc:
        error = str(exc).replace(str(run), artifact_label(run)).replace(str(ROOT), "<repository>")
        report = {"status": "failed", "validated_utc": datetime.now(timezone.utc).isoformat(),
                  "error_type": type(exc).__name__, "error": error,
                  "original_GDSC_or_CTRP_response_read_by_validator": False}
    (run / "pilot_validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "report": artifact_label(run / "pilot_validation.json"), "counts": report.get("counts"), "error": report.get("error")}, ensure_ascii=False))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
