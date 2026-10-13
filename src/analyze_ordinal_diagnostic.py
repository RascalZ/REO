"""Frozen descriptive statistics for the four-case ordinal diagnostic pilot.

Reads only completed development OOF/audits. Observed-label intervals share
global biological-group bootstrap counts; all nine null runs remain point-only.
No model fits, permutation P-values or source-platform comparison to v1 occur.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from analyze_formal_results import (METRICS, clean, dump, generate_global_bootstrap_weights,
    install_response_source_guard, interval_summary, sha, weighted_metrics)

MACRO = "macro_equal_four_diagnostic_cases"


def planned_contrasts(config):
    contrasts = []
    for rule in config["selection_rules"]:
        for candidate, reference in [("rank", "expression"), ("pair_full", "rank"), ("pair_full", "expression")]:
            contrasts.append((rule + ":" + candidate + "-" + reference, rule + "_" + candidate, rule + "_" + reference))
    for representation in config["representations"]:
        contrasts.append((representation + ":association-MAD", "response_association_" + representation, "mad_" + representation))
    for model in config["models"]:
        if model != "tissue_null":
            contrasts.append((model + "-tissue_null", model, "tissue_null"))
    if len(contrasts) != 15:
        raise ValueError("Expected the frozen 15 contrasts")
    return contrasts


def statistics_config(config):
    return {"bootstrap": {"minimum_valid_draw_fraction": config["statistics"]["minimum_valid_fraction"], "unadjusted_level": .95}}


def analyze_predictions(data, global_manifest, config, groups=None, counts=None):
    """Return public point/CI/null aggregates plus private shared count arrays."""
    required = {"drug_id", "model", "repeat", "outer_fold", "permutation_id", "canonical_model_id",
                "biological_group_id", "tissue", "observed_auc", "predicted_auc"}
    if not required <= set(data):
        raise ValueError("Missing OOF columns: " + str(required - set(data)))
    data = data.copy()
    for key in ("drug_id", "repeat", "outer_fold", "permutation_id"):
        data[key] = data[key].astype(int)
    if data.duplicated(["drug_id", "model", "permutation_id", "canonical_model_id"]).any():
        raise ValueError("Duplicate OOF identities")
    if set(data.model) != set(config["models"]) or set(data.drug_id) != set(config["drug_ids"]):
        raise ValueError("Unexpected/incomplete model or case set")
    if set(data["repeat"]) != {0} or set(data.permutation_id) != set(range(config["permutation"]["count"] + 1)):
        raise ValueError("Unexpected repeat/permutation coverage")
    if set(data.outer_fold) != set(range(config["outer_folds"])):
        raise ValueError("Unexpected outer folds")
    if data.groupby(["drug_id", "canonical_model_id"]).outer_fold.nunique().ne(1).any():
        raise ValueError("Outer role changed across model/permutation")
    if groups is None or counts is None:
        groups, counts = generate_global_bootstrap_weights(global_manifest.biological_group_id,
            config["statistics"]["bootstrap_draws"], config["statistics"]["bootstrap_seed"])
    expected_groups = np.asarray(sorted(set(global_manifest.biological_group_id.astype(str))), dtype=str)
    if not np.array_equal(groups, expected_groups):
        raise ValueError("Shared global bootstrap group order mismatch")
    if counts.shape != (config["statistics"]["bootstrap_draws"], len(groups)):
        raise ValueError("Wrong bootstrap shape")
    if not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all() or not (counts.sum(1) == len(groups)).all():
        raise ValueError("Invalid shared global counts")
    positions = {group: i for i, group in enumerate(groups)}
    points, bootstrap, molecular_deviations, ns = {}, {}, {}, {}
    actual_case_groups = int(data[data.permutation_id.eq(0)].biological_group_id.nunique())
    for drug in config["drug_ids"]:
        source = global_manifest[global_manifest.gdsc_drug_id.astype(int).eq(drug)].sort_values("canonical_model_id")
        if source.empty or source.canonical_model_id.duplicated().any() or source.biological_group_id.duplicated().any():
            raise ValueError("Unique model per biological group required in each drug")
        ids = source.canonical_model_id.tolist()
        reference_groups = source.biological_group_id.to_numpy(str)
        reference_tissues = source.tissue.to_numpy(str)
        ns[drug] = len(ids)
        index = np.asarray([positions[g] for g in reference_groups])
        for perm in range(config["permutation"]["count"] + 1):
            subset = data[data.drug_id.eq(drug) & data.permutation_id.eq(perm)]
            if subset.groupby("canonical_model_id").observed_auc.nunique(dropna=False).ne(1).any():
                raise ValueError("Assigned outcomes differ across arms")
            null_records = subset[subset.model.eq("tissue_null")].set_index("canonical_model_id").reindex(ids)
            null_prediction = null_records.predicted_auc.to_numpy(float)
            y = null_records.observed_auc.to_numpy(float)
            for model in config["models"]:
                record = subset[subset.model.eq(model)].set_index("canonical_model_id")
                if set(record.index) != set(ids):
                    raise ValueError("Incomplete OOF case/model/permutation; no case dropping")
                record = record.reindex(ids)
                if not np.array_equal(record.biological_group_id.to_numpy(str), reference_groups) or not np.array_equal(record.tissue.to_numpy(str), reference_tissues):
                    raise ValueError("Frozen identity/tissue mismatch")
                if not np.array_equal(record.observed_auc.to_numpy(float), y):
                    raise ValueError("Assigned outcomes changed")
                p = record.predicted_auc.to_numpy(float)
                if not np.isfinite(y).all() or not np.isfinite(p).all():
                    raise ValueError("Nonfinite response/prediction requires preserved run failure")
                points[(drug, perm, model)] = {metric: value[0] for metric, value in weighted_metrics(y, p, np.ones((1, len(y)), dtype=int)).items()}
                molecular_deviations[(drug, perm, model)] = float(np.sqrt(np.mean((p - null_prediction) ** 2)))
                if perm == 0:
                    samples = {metric: np.full(len(counts), np.nan) for metric in METRICS}
                    for start in range(0, len(counts), 250):
                        stop = min(start + 250, len(counts))
                        values = weighted_metrics(y, p, counts[start:stop, index])
                        for metric in METRICS:
                            samples[metric][start:stop] = values[metric]
                    bootstrap[(drug, model)] = samples
    for perm in range(config["permutation"]["count"] + 1):
        for model in config["models"]:
            points[(MACRO, perm, model)] = {metric: np.mean([points[(drug, perm, model)][metric] for drug in config["drug_ids"]]) for metric in METRICS}
            # Average of four case-specific RMS prediction changes, not pooled RMS.
            molecular_deviations[(MACRO, perm, model)] = np.mean([molecular_deviations[(drug, perm, model)] for drug in config["drug_ids"]])
    for model in config["models"]:
        bootstrap[(MACRO, model)] = {metric: np.mean(np.stack([bootstrap[(drug, model)][metric] for drug in config["drug_ids"]]), axis=0) for metric in METRICS}
    point_rows, model_ci, effect_ci, null_effects = [], [], [], []
    stats_config = statistics_config(config)
    for drug in [*config["drug_ids"], MACRO]:
        for perm in range(config["permutation"]["count"] + 1):
            for model in config["models"]:
                point_rows.append({"gdsc_drug_id": drug, "permutation_id": perm, "model": model,
                    "n_unique_biological_groups": actual_case_groups if drug == MACRO else ns[drug],
                    "n_bootstrap_global_groups": len(groups),
                    "prediction_rms_difference_from_tissue_null": molecular_deviations[(drug, perm, model)],
                    **points[(drug, perm, model)], "scope": "observed_exploratory" if perm == 0 else "full_pipeline_restricted_permutation_control"})
            for metric in METRICS:
                sign = -1 if metric in ("rmse", "mae") else 1
                if perm == 0:
                    for model in config["models"]:
                        point = points[(drug, perm, model)][metric]
                        model_ci.append({"gdsc_drug_id": drug, "model": model, "metric": metric,
                            **interval_summary(point, bootstrap[(drug, model)][metric], stats_config, reason="undefined_constant_response_or_prediction")})
                for name, candidate, reference in planned_contrasts(config):
                    point = sign * (points[(drug, perm, candidate)][metric] - points[(drug, perm, reference)][metric])
                    common = {"gdsc_drug_id": drug, "contrast": name, "candidate": candidate, "reference": reference,
                        "metric": metric, "positive_favors_candidate": True}
                    if perm == 0:
                        samples = sign * (bootstrap[(drug, candidate)][metric] - bootstrap[(drug, reference)][metric])
                        effect_ci.append({**common, **interval_summary(point, samples, stats_config, reason="undefined_paired_metric"),
                            "analysis_status": "unadjusted_exploratory_conditional_OOF_not_full_train_refit"})
                    else:
                        null_effects.append({**common, "permutation_id": perm, "point": point,
                            "scope": "point_only_full_pipeline_null_no_P_value"})
    null_frame = pd.DataFrame(null_effects)
    null_summary = []
    for key, frame in null_frame.groupby(["gdsc_drug_id", "contrast", "candidate", "reference", "metric"], sort=False):
        finite = np.isfinite(frame.point)
        # Do not silently summarize a smaller-than-frozen collection of draws.
        values = frame.point.to_numpy(float)
        null_summary.append(dict(zip(["gdsc_drug_id", "contrast", "candidate", "reference", "metric"], key)) |
            {"null_runs_requested": config["permutation"]["count"], "null_runs_recorded": len(frame),
             "null_runs_valid": int(finite.sum()), "null_min": np.min(values), "null_median": np.median(values), "null_max": np.max(values),
             "summary_status": "complete" if finite.all() else "undefined_constituent_null_retained"})
    return {"point_metrics": pd.DataFrame(point_rows), "model_intervals": pd.DataFrame(model_ci),
        "paired_effect_intervals": pd.DataFrame(effect_ci), "null_contrasts": null_frame,
        "null_contrast_summary": pd.DataFrame(null_summary), "bootstrap_group_order": groups, "bootstrap_counts": counts}


def aggregate_audits(path, config):
    states, resources, chosen, features = {}, [], [], defaultdict(list)
    shared_gene_hash = {}
    seen = set()
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            task, audit = record["task"], record["outer_fit_audit"]
            key = (int(task["drug_id"]), int(task["permutation_id"]), task["model"])
            fullkey = (*key, int(task["outer_fold"]))
            if fullkey in seen:
                raise ValueError("Duplicate outer-fit audit")
            seen.add(fullkey)
            state = states.setdefault(key, {"gdsc_drug_id": key[0], "permutation_id": key[1], "model": key[2], "outer_fits": 0,
                "zero_molecular_fits": 0, "nonzero_molecular_fits": 0, "warning_outer_fits": 0, "fallback_outer_fits": 0,
                "active_feature_min": None, "active_feature_max": None, "active_feature_total": 0, "outer_unknown_tissue_samples": 0})
            state["outer_fits"] += 1
            state["warning_outer_fits"] += int(bool(audit.get("warnings")))
            state["fallback_outer_fits"] += int(bool(audit.get("fallback_reason")))
            state["outer_unknown_tissue_samples"] += record.get("outer_unknown_tissue_samples", 0)
            nonzero = audit.get("n_nonzero_molecular_coeffs")
            if key[2] != "tissue_null":
                if nonzero is None or len(audit["active_molecular_features"]) != nonzero:
                    raise ValueError("Missing/inconsistent active coefficient audit")
                state["zero_molecular_fits"] += int(nonzero == 0)
                state["nonzero_molecular_fits"] += int(nonzero > 0)
                state["active_feature_total"] += nonzero
                state["active_feature_min"] = nonzero if state["active_feature_min"] is None else min(state["active_feature_min"], nonzero)
                state["active_feature_max"] = nonzero if state["active_feature_max"] is None else max(state["active_feature_max"], nonzero)
                gene_key = (key[0], key[1], task["outer_fold"], audit["selection_rule"])
                digest = audit["selected_gene_ids_sha256"]
                if gene_key in shared_gene_hash and shared_gene_hash[gene_key] != digest:
                    raise ValueError("Candidate genes differ across representations in the same outer fit")
                shared_gene_hash[gene_key] = digest
                active = dict(zip(audit["active_molecular_features"], audit["active_molecular_coefficients_standardized"]))
                features[key].append((int(task["outer_fold"]), set(audit["selected_gene_ids"]), active))
            chosen.append({"gdsc_drug_id": key[0], "permutation_id": key[1], "model": key[2], "outer_fold": task["outer_fold"],
                "chosen_params": json.dumps(record["chosen_params"], sort_keys=True), "best_inner_rmse": record["best_inner_rmse"],
                "selected_genes": audit.get("n_selected_genes"), "candidate_features": audit.get("n_candidate_molecular_features"),
                "active_features": nonzero, "outer_solver_iterations": audit.get("n_iter"),
                "warnings": json.dumps(audit.get("warnings", []), ensure_ascii=False), "fallback_reason": audit.get("fallback_reason", "")})
            resource = record["resources"]
            resources.append({"gdsc_drug_id": key[0], "permutation_id": key[1], "model": key[2], "outer_fold": task["outer_fold"],
                "task_wall_seconds": resource["wall_seconds"], "task_cpu_seconds": resource.get("cpu_seconds"),
                "sampled_peak_rss_bytes": resource["sampled_peak_rss_bytes"],
                "windows_lifetime_peak_working_set_bytes": resource.get("windows_lifetime_peak_working_set_bytes"),
                "outer_fit_seconds": audit["fit_seconds"], "outer_prediction_seconds": audit["last_predict_seconds"]})
    expected = {(drug, perm, model, fold) for drug in config["drug_ids"] for perm in range(config["permutation"]["count"] + 1)
                for model in config["models"] for fold in range(config["outer_folds"])}
    if seen != expected:
        raise ValueError("Outer-fit audit coverage incomplete; preserve failures and stop analysis")
    stability = []
    for key, folds in features.items():
        for left, right in itertools.combinations(sorted(folds, key=lambda x: x[0]), 2):
            fa, ga, aa = left
            fb, gb, ab = right
            gene_union, active_union = ga | gb, set(aa) | set(ab)
            shared = set(aa) & set(ab)
            stability.append({"gdsc_drug_id": key[0], "permutation_id": key[1], "model": key[2], "fold_a": fa, "fold_b": fb,
                "candidate_gene_jaccard": len(ga & gb) / len(gene_union) if gene_union else np.nan,
                "active_feature_jaccard": len(shared) / len(active_union) if active_union else np.nan,
                "same_sign_active_jaccard": sum(np.sign(aa[name]) == np.sign(ab[name]) for name in shared) / len(active_union) if active_union else np.nan,
                "empty_active_union": not bool(active_union), "scope": "three_overlapping_training_folds_descriptive_not_independent_repeats"})
    resource_frame = pd.DataFrame(resources)
    summary = resource_frame.groupby(["gdsc_drug_id", "permutation_id", "model"], sort=False).agg(
        outer_tasks=("outer_fold", "size"), task_wall_seconds_sum=("task_wall_seconds", "sum"),
        task_cpu_seconds_sum=("task_cpu_seconds", "sum"), sampled_peak_rss_bytes_max=("sampled_peak_rss_bytes", "max"),
        windows_lifetime_peak_working_set_bytes_max=("windows_lifetime_peak_working_set_bytes", "max"),
        outer_fit_seconds_sum=("outer_fit_seconds", "sum"), outer_prediction_seconds_sum=("outer_prediction_seconds", "sum")).reset_index()
    return {"molecular_and_failure_summary": pd.DataFrame(states.values()), "chosen_parameters": pd.DataFrame(chosen),
            "active_feature_stability": pd.DataFrame(stability), "resource_tasks": resource_frame, "resource_summary": summary}


def aggregate_fit_events(path, config):
    states = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            event = json.loads(line)
            key = (int(event["drug_id"]), int(event["permutation_id"]), event["model"])
            state = states.setdefault(key, {"gdsc_drug_id": key[0], "permutation_id": key[1], "model": key[2],
                "actual_fit_events": 0, "all_fit_warning_events": 0, "all_fit_fallback_events": 0,
                "fit_seconds_sum": 0., "prediction_seconds_sum": 0.})
            state["actual_fit_events"] += 1
            state["all_fit_warning_events"] += int(bool(event.get("warnings")))
            state["all_fit_fallback_events"] += int(bool(event.get("fallback_reason")))
            state["fit_seconds_sum"] += event["fit_seconds"]
            state["prediction_seconds_sum"] += event["predict_seconds"]
    frame = pd.DataFrame(states.values())
    for row in frame.itertuples():
        expected = config["outer_folds"] * (1 if row.model == "tissue_null" else config["inner_folds"] * len(config["elastic_net"]["alpha_ratios"]) * len(config["elastic_net"]["l1_ratios"]) + 1)
        if row.actual_fit_events != expected:
            raise ValueError("Fit-event count differs from frozen grid")
    if len(frame) != len(config["drug_ids"]) * len(config["models"]) * (config["permutation"]["count"] + 1):
        raise ValueError("Missing fit-event task group")
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/ordinal_diagnostic_v1.json")
    parser.add_argument("--run-dir")
    parser.add_argument("--output-dir", default="results/ordinal_diagnostic")
    args = parser.parse_args()
    config_path = (ROOT / args.config).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    run = (ROOT / (args.run_dir or "runs/" + config["run_name"])).resolve()
    output = (ROOT / args.output_dir).resolve()
    for path in (config_path, run, output):
        path.relative_to(ROOT)
    if output.exists() or (run / "statistics_private").exists():
        raise FileExistsError("Existing diagnostic statistics are retained; use a new version")
    blocked = install_response_source_guard()
    started = time.perf_counter()
    metadata_path = run / "ordinal_metadata.json"
    if not metadata_path.exists():
        metadata_path = run / "formal_metadata.json"
    meta = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected_tasks = len(config["drug_ids"]) * len(config["models"]) * config["outer_folds"] * (config["permutation"]["count"] + 1)
    if meta.get("completed") is not True or meta.get("completed_tasks") != expected_tasks or meta.get("external_response_accessed") is not False:
        raise RuntimeError("Completed sealed diagnostic run required")
    if meta.get("config_sha256") != sha(config_path) or meta.get("data_guard", {}).get("blocked_attempts") != 0:
        raise RuntimeError("Changed config or nonzero data-guard attempts")
    lock_path = ROOT / "config/ordinal_diagnostic_v1.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if meta.get("freeze_sha256") != sha(lock_path):
        raise RuntimeError("Diagnostic lock mismatch")
    for filename, checksum in lock["files"].items():
        if sha(ROOT / filename) != checksum:
            raise RuntimeError("Frozen diagnostic file changed: " + filename)
    for name in ("oof_predictions.csv", "outer_fit_audits.jsonl", "fit_events.jsonl"):
        if meta.get("output_hashes", {}).get(name) != sha(run / name):
            raise RuntimeError("Bound output changed: " + name)
    manifest_path = ROOT / "results/cohort_v2/development_sample_manifest.csv"
    manifest = pd.read_csv(manifest_path, dtype={"canonical_model_id": str, "biological_group_id": str})
    data = pd.read_csv(run / "oof_predictions.csv", dtype={"canonical_model_id": str, "biological_group_id": str}, float_precision="round_trip")
    with threadpool_limits(limits=1):
        result = analyze_predictions(data, manifest, config)
    result.update(aggregate_audits(run / "outer_fit_audits.jsonl", config))
    result["all_fit_events_summary"] = aggregate_fit_events(run / "fit_events.jsonl", config)
    coverage_path = run / "permutation_coverage.csv"
    if not coverage_path.exists():
        raise FileNotFoundError("Required aggregate permutation coverage missing")
    bound_coverage = meta.get("output_hashes", {}).get("permutation_coverage.csv", meta.get("prepared_input_hashes", {}).get("permutation_coverage.csv"))
    if bound_coverage != sha(coverage_path):
        raise RuntimeError("Permutation-coverage hash not bound to completed run")
    result["permutation_coverage"] = pd.read_csv(coverage_path)
    output.mkdir(parents=True, exist_ok=False)
    private = run / "statistics_private"
    private.mkdir(exist_ok=False)
    groups, counts = result.pop("bootstrap_group_order"), result.pop("bootstrap_counts")
    pd.DataFrame({"bootstrap_group_position": np.arange(len(groups)), "biological_group_id": groups}).to_csv(private / "bootstrap_group_manifest.csv", index=False, lineterminator="\n")
    np.savez_compressed(private / "bootstrap_group_multiplicities.npz", multiplicities=counts)
    for name, frame in result.items():
        frame.to_csv(output / (name + ".csv"), index=False, encoding="utf-8-sig", lineterminator="\n")
    memory = psutil.Process().memory_info()
    provenance = {"status": "completed", "created_utc": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": time.perf_counter() - started, "process_lifetime_peak_rss_bytes": getattr(memory, "peak_wset", memory.rss),
        "source_run": run.relative_to(ROOT).as_posix(), "execution_commit": meta.get("execution_commit"),
        "analyzer_sha256": sha(Path(__file__)), "config_sha256": sha(config_path), "lock_sha256": sha(lock_path),
        "metadata_sha256": sha(metadata_path), "OOF_source_sha256": sha(run / "oof_predictions.csv"),
        "development_manifest_sha256": sha(manifest_path), "n_global_biological_groups": len(groups),
        "bootstrap_draws": len(counts), "bootstrap_seed": config["statistics"]["bootstrap_seed"],
        "bootstrap_count_matrix_sha256": hashlib.sha256(counts.tobytes(order="C")).hexdigest(),
        "bootstrap_group_order_sha256": hashlib.sha256("\0".join(groups).encode()).hexdigest(),
        "expected_outer_tasks": expected_tasks, "actual_fit_events": int(result["all_fit_events_summary"].actual_fit_events.sum()),
        "oof_rows": len(data), "no_repeat_N_inflation": True, "all_null_runs_retained": True,
        "bootstrap_shared_across_models_drugs": True, "new_model_fits": 0,
        "confidence_scope": config["statistics"]["confidence_scope"], "no_permutation_p_values": True,
        "external_response_accessed": False, "blocked_source_access_attempts": blocked,
        "output_sha256": {name + ".csv": sha(output / (name + ".csv")) for name in result}}
    dump(output / "analysis_provenance.json", provenance)
    print(json.dumps(clean(provenance), ensure_ascii=False))


if __name__ == "__main__":
    main()
