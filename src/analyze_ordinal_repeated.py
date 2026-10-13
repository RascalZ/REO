"""Frozen statistics for the fixed-eight repeated ordinal evaluation.

Every metric is computed on a complete OOF repeat before repeat averaging.
One global biological-group count matrix is shared by drugs, models and repeats.
Intervals condition on recorded OOF; they do not cover retraining uncertainty.
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
from analyze_formal_results import (METRICS, clean, generate_global_bootstrap_weights,
    install_response_source_guard, interval_summary, sha, weighted_metrics)

MACRO = "macro_equal_fixed_eight_drugs"
RULES = ("mad", "response_association")
REPRESENTATIONS = ("expression", "rank", "pair_full")
MOLECULAR200 = tuple(rule + "_" + representation for rule in RULES for representation in REPRESENTATIONS)


def write_json(path, value):
    Path(path).write_bytes((json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode("utf-8"))


def planned_contrasts(config):
    """Twenty-five fixed contrasts: 22 paired comparisons and three interactions."""
    output = [{"contrast": "primary:association_pair-strong_expression", "coefficients": {
        "response_association_pair_full": 1., "expression_strong": -1.}, "family": "primary_fixed9"}]
    for rule in RULES:
        for candidate, reference in (("pair_full", "expression"), ("rank", "expression"), ("pair_full", "rank")):
            output.append({"contrast": rule + ":" + candidate + "-" + reference,
                "coefficients": {rule + "_" + candidate: 1., rule + "_" + reference: -1.},
                "family": "representation_fixed18" if (candidate, reference) == ("pair_full", "expression") else "descriptive"})
    for representation in REPRESENTATIONS:
        output.append({"contrast": representation + ":association-MAD", "coefficients": {
            "response_association_" + representation: 1., "mad_" + representation: -1.}, "family": "descriptive"})
    for model in MOLECULAR200:
        output.append({"contrast": model + "-tissue_null", "coefficients": {model: 1., "tissue_null": -1.}, "family": "descriptive"})
    for model in MOLECULAR200:
        if model != "response_association_pair_full":
            output.append({"contrast": model + "-expression_strong", "coefficients": {
                model: 1., "expression_strong": -1.}, "family": "descriptive"})
    output.append({"contrast": "expression_strong-tissue_null", "coefficients": {
        "expression_strong": 1., "tissue_null": -1.}, "family": "descriptive"})
    for candidate, reference in (("rank", "expression"), ("pair_full", "rank"), ("pair_full", "expression")):
        output.append({"contrast": "interaction:" + candidate + "-" + reference,
            "coefficients": {"response_association_" + candidate: 1., "response_association_" + reference: -1.,
                             "mad_" + candidate: -1., "mad_" + reference: 1.}, "family": "descriptive"})
    if len(output) != 25 or any(not set(row["coefficients"]) <= set(config["models"]) for row in output):
        raise ValueError("Frozen contrast definitions disagree with model set")
    return output


def contrast_values(model_values, coefficients, metric):
    """Linear paired contrast; undefined terms propagate, never nanmean."""
    sign = -1. if metric in ("rmse", "mae") else 1.
    return sign * sum(coefficient * np.asarray(model_values[model], dtype=float) for model, coefficient in coefficients.items())


def summary_config(config):
    s = config["statistics"]
    return {"bootstrap": {"minimum_valid_draw_fraction": s["minimum_valid_fraction"], "unadjusted_level": .95},
            "multiplicity": {"family_alpha": s.get("family_alpha", .05)}}


def repeat_summary(values):
    values = np.asarray(values, dtype=float)
    return {"repeat_values": json.dumps(clean(values.tolist()), allow_nan=False), "repeat_mean": np.mean(values),
        "repeat_min": np.min(values), "repeat_max": np.max(values), "repeat_range": np.ptp(values),
        "repeat_sample_sd": np.std(values, ddof=1), "positive_repeats": int((values > 0).sum()),
        "negative_repeats": int((values < 0).sum()), "zero_repeats": int((values == 0).sum()),
        "undefined_repeats": int((~np.isfinite(values)).sum()),
        "scope": "descriptive_repeated_split_variation_not_independent_study_SE"}


def analyze_predictions(data, manifest, config, groups=None, counts=None):
    """Pure-data statistical API; no training or access to external responses."""
    required = {"drug_id", "model", "repeat", "outer_fold", "canonical_model_id", "biological_group_id", "tissue", "observed_auc", "predicted_auc"}
    if not required <= set(data):
        raise ValueError("Required OOF columns missing")
    data = data.copy()
    for field in ("drug_id", "repeat", "outer_fold"):
        data[field] = data[field].astype(int)
    drugs, models = config["drug_ids"], config["models"]
    repeats = list(range(config["outer_repeats"]))
    if set(data.drug_id) != set(drugs) or set(data.model) != set(models) or set(data["repeat"]) != set(repeats):
        raise ValueError("Missing/unplanned drug, model or repeat")
    if "permutation_id" in data and not data.permutation_id.eq(0).all():
        raise ValueError("This repeated run contains observed labels only")
    if data.duplicated(["drug_id", "model", "repeat", "canonical_model_id"]).any():
        raise ValueError("Duplicate complete-repeat OOF identity")
    if set(data.outer_fold) != set(range(config["outer_folds"])):
        raise ValueError("Unplanned outer-fold values")
    if data.groupby(["drug_id", "repeat", "canonical_model_id"]).outer_fold.nunique().ne(1).any():
        raise ValueError("Outer roles differ across paired models")
    if groups is None or counts is None:
        groups, counts = generate_global_bootstrap_weights(manifest.biological_group_id,
            config["statistics"]["bootstrap_draws"], config["statistics"]["bootstrap_seed"])
    groups, counts = np.asarray(groups, dtype=str), np.asarray(counts)
    if not np.array_equal(groups, np.asarray(sorted(set(manifest.biological_group_id.astype(str))), dtype=str)):
        raise ValueError("Shared bootstrap group order mismatch")
    if counts.shape != (config["statistics"]["bootstrap_draws"], len(groups)) or not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all() or not (counts.sum(1) == len(groups)).all():
        raise ValueError("Invalid global group bootstrap counts")
    group_position = {group: i for i, group in enumerate(groups)}
    points, boot, ns, deviations = {}, {}, {}, {}
    repeat_rows, coverage_rows = [], []
    for drug in drugs:
        source = manifest[manifest.gdsc_drug_id.astype(int).eq(drug)].sort_values("canonical_model_id")
        if source.empty or source.canonical_model_id.duplicated().any() or source.biological_group_id.duplicated().any():
            raise ValueError("One model per frozen drug-biological-group required")
        ids = source.canonical_model_id.tolist()
        ns[drug] = source.biological_group_id.nunique()
        positions = np.array([group_position[g] for g in source.biological_group_id])
        subset = data[data.drug_id.eq(drug)]
        if subset.groupby("canonical_model_id").observed_auc.nunique(dropna=False).ne(1).any():
            raise ValueError("Observed response changed across models/repeats")
        y = subset.groupby("canonical_model_id").observed_auc.first().reindex(ids).to_numpy(float)
        null_predictions = {}
        for repeat in repeats:
            reference = subset[subset.model.eq("tissue_null") & subset["repeat"].eq(repeat)].set_index("canonical_model_id")
            if set(reference.index) != set(ids):
                raise ValueError("Incomplete tissue-null repeat")
            null_predictions[repeat] = reference.reindex(ids).predicted_auc.to_numpy(float)
        for model in models:
            point = np.full((len(repeats), len(METRICS)), np.nan)
            # Keep only one model's repeat x draw array, then retain its mean.
            sample = np.full((len(repeats), len(counts), len(METRICS)), np.nan)
            deviation = np.full(len(repeats), np.nan)
            for repeat in repeats:
                record = subset[subset.model.eq(model) & subset["repeat"].eq(repeat)].set_index("canonical_model_id")
                if set(record.index) != set(ids) or set(record.outer_fold) != set(range(config["outer_folds"])):
                    raise ValueError("Incomplete paired OOF repeat; no model/case/repeat deletion")
                record = record.reindex(ids)
                if not np.array_equal(record.biological_group_id.to_numpy(str), source.biological_group_id.to_numpy(str)) or not np.array_equal(record.tissue.to_numpy(str), source.tissue.to_numpy(str)):
                    raise ValueError("Frozen OOF identity or tissue mismatch")
                if not np.array_equal(record.observed_auc.to_numpy(float), y):
                    raise ValueError("OOF labels differ from paired source")
                prediction = record.predicted_auc.to_numpy(float)
                if not np.isfinite(prediction).all() or not np.isfinite(y).all():
                    raise ValueError("Nonfinite outcome/prediction; preserve failed run before inference")
                current = weighted_metrics(y, prediction, np.ones((1, len(ids)), dtype=int))
                point[repeat] = [current[metric][0] for metric in METRICS]
                deviation[repeat] = float(np.sqrt(np.mean((prediction - null_predictions[repeat]) ** 2)))
                common = {"gdsc_drug_id": drug, "model": model, "repeat": repeat,
                    "n_unique_models": len(ids), "n_unique_biological_groups": ns[drug], "outer_folds": config["outer_folds"]}
                repeat_rows.append({**common, **{metric: point[repeat, i] for i, metric in enumerate(METRICS)},
                    "prediction_rms_difference_from_tissue_null": deviation[repeat],
                    "spearman_unavailable_reason": "constant_response_or_prediction" if not np.isfinite(point[repeat, 0]) else ""})
                coverage_rows.append({**common, "expected_oof_rows": len(ids), "observed_oof_rows": len(record), "status": "complete"})
                for start in range(0, len(counts), 250):
                    stop = min(start + 250, len(counts))
                    current = weighted_metrics(y, prediction, counts[start:stop, positions])
                    for i, metric in enumerate(METRICS):
                        sample[repeat, start:stop, i] = current[metric]
            points[(drug, model)], boot[(drug, model)] = point, np.mean(sample, axis=0)
            deviations[(drug, model)] = deviation
    ns[MACRO] = data.biological_group_id.nunique()
    for model in models:
        points[(MACRO, model)] = np.mean(np.stack([points[(drug, model)] for drug in drugs]), axis=0)
        boot[(MACRO, model)] = np.mean(np.stack([boot[(drug, model)] for drug in drugs]), axis=0)
        deviations[(MACRO, model)] = np.mean(np.stack([deviations[(drug, model)] for drug in drugs]), axis=0)
    model_rows, effects, stability, rms_rows = [], [], [], []
    sconfig = summary_config(config)
    for drug in [*drugs, MACRO]:
        for model in models:
            rms_rows.append({"gdsc_drug_id": drug, "model": model,
                "n_unique_biological_groups": ns[drug], **repeat_summary(deviations[(drug, model)])})
            for i, metric in enumerate(METRICS):
                values = points[(drug, model)][:, i]
                model_rows.append({"gdsc_drug_id": drug, "model": model, "metric": metric,
                    "n_unique_biological_groups": ns[drug], "n_bootstrap_global_groups": len(groups),
                    **interval_summary(np.mean(values), boot[(drug, model)][:, i], sconfig, reason="undefined_constituent_repeat_or_drug")})
                stability.append({"gdsc_drug_id": drug, "comparison": "model:" + model, "metric": metric, **repeat_summary(values)})
        for contrast in planned_contrasts(config):
            for i, metric in enumerate(METRICS):
                repeated = contrast_values({model: points[(drug, model)][:, i] for model in models}, contrast["coefficients"], metric)
                resampled = contrast_values({model: boot[(drug, model)][:, i] for model in models}, contrast["coefficients"], metric)
                family = contrast["family"] if metric == "spearman" else "descriptive"
                size = config["statistics"].get("primary_family_size", 9) if family == "primary_fixed9" else config["statistics"].get("representation_family_size", 18) if family == "representation_fixed18" else None
                coefficients = contrast["coefficients"]
                candidate = next((model for model, coefficient in coefficients.items() if coefficient == 1.), "") if len(coefficients) == 2 else "interaction"
                reference = next((model for model, coefficient in coefficients.items() if coefficient == -1.), "") if len(coefficients) == 2 else "difference_of_differences"
                effects.append({"gdsc_drug_id": drug, "contrast": contrast["contrast"], "candidate": candidate, "reference": reference,
                    "metric": metric, "coefficients": json.dumps(coefficients, sort_keys=True), "family_name": family,
                    "positive_favors_candidate_or_positive_interaction": True,
                    "effect_direction": "reference_minus_candidate" if metric in ("rmse", "mae") else "candidate_minus_reference",
                    **interval_summary(np.mean(repeated), resampled, sconfig, family_size=size, reason="undefined_paired_repeat_or_drug"),
                    "analysis_scope": "development_inspired_conditional_OOF_not_independent_confirmation_or_noninferiority"})
                stability.append({"gdsc_drug_id": drug, "comparison": contrast["contrast"], "metric": metric, **repeat_summary(repeated)})
    return {"repeat_metrics": pd.DataFrame(repeat_rows), "coverage_summary": pd.DataFrame(coverage_rows),
        "model_intervals": pd.DataFrame(model_rows), "paired_effect_intervals": pd.DataFrame(effects),
        "repeat_stability": pd.DataFrame(stability), "prediction_deviation_from_tissue": pd.DataFrame(rms_rows),
        "bootstrap_group_order": groups, "bootstrap_counts": counts}


def aggregate_audits(path, config):
    states, choices, resource_rows = {}, [], []
    features, shared_genes, seen = defaultdict(list), {}, set()
    for line in Path(path).open(encoding="utf-8"):
        if not line.strip():
            continue
        record = json.loads(line)
        task, audit, resource = record["task"], record["outer_fit_audit"], record["resources"]
        drug, model, repeat, fold = int(task["drug_id"]), task["model"], int(task["repeat"]), int(task["outer_fold"])
        key, context = (drug, model, repeat), (drug, model, repeat, fold)
        if context in seen:
            raise ValueError("Duplicate completed outer audit")
        seen.add(context)
        state = states.setdefault(key, {"gdsc_drug_id": drug, "model": model, "repeat": repeat, "outer_fits": 0,
            "zero_molecular_fits": 0, "nonzero_molecular_fits": 0, "warning_outer_fits": 0, "fallback_outer_fits": 0,
            "active_feature_min": None, "active_feature_max": None, "active_feature_total": 0, "outer_unknown_tissue_samples": 0})
        state["outer_fits"] += 1
        state["warning_outer_fits"] += int(bool(audit.get("warnings")))
        state["fallback_outer_fits"] += int(bool(audit.get("fallback_reason")))
        state["outer_unknown_tissue_samples"] += record.get("outer_unknown_tissue_samples", 0)
        active_count = audit.get("n_nonzero_molecular_coeffs")
        if model != "tissue_null":
            active_names, active_values = audit["active_molecular_features"], audit["active_molecular_coefficients_standardized"]
            if active_count is None or len(active_names) != active_count or len(active_values) != active_count or len(set(active_names)) != active_count or any(not np.isfinite(value) or value == 0 for value in active_values):
                raise ValueError("Invalid active coefficient audit")
            state["zero_molecular_fits"] += int(active_count == 0)
            state["nonzero_molecular_fits"] += int(active_count > 0)
            state["active_feature_total"] += active_count
            state["active_feature_min"] = active_count if state["active_feature_min"] is None else min(state["active_feature_min"], active_count)
            state["active_feature_max"] = active_count if state["active_feature_max"] is None else max(state["active_feature_max"], active_count)
            if model in MOLECULAR200:
                gene_key = (drug, repeat, fold, audit["selection_rule"])
                digest = audit["selected_gene_ids_sha256"]
                if gene_key in shared_genes and shared_genes[gene_key] != digest:
                    raise ValueError("Three representations differ in outer candidate gene set")
                shared_genes[gene_key] = digest
            features[(drug, model)].append((repeat, fold, set(audit["selected_gene_ids"]),
                dict(zip(active_names, active_values)), set(audit["fit_sample_ids"])))
        choices.append({"gdsc_drug_id": drug, "model": model, "repeat": repeat, "outer_fold": fold,
            "chosen_params": json.dumps(record["chosen_params"], sort_keys=True), "best_inner_rmse": record["best_inner_rmse"],
            "chosen_grid_boundary": record["chosen_grid_boundary"],
            "chosen_grid_boundary_details": json.dumps(record["chosen_grid_boundary_details"], sort_keys=True),
            "selected_genes": audit.get("n_selected_genes"), "candidate_features": audit.get("n_candidate_molecular_features"),
            "active_features": active_count, "outer_solver_iterations": audit.get("n_iter"), "actual_alpha": audit.get("actual_alpha"),
            "alpha_max": audit.get("alpha_max"), "warnings": json.dumps(audit.get("warnings", []), ensure_ascii=False),
            "fallback_reason": audit.get("fallback_reason", "")})
        resource_rows.append({"gdsc_drug_id": drug, "model": model, "repeat": repeat, "outer_fold": fold,
            "task_wall_seconds": resource["wall_seconds"], "fit_work_wall_seconds": resource.get("fit_work_wall_seconds"),
            "initialization_wall_seconds": resource.get("initialization_wall_seconds"), "task_cpu_seconds": resource.get("cpu_seconds"),
            "sampled_peak_rss_bytes": resource["sampled_peak_rss_bytes"],
            "windows_lifetime_peak_working_set_bytes": resource.get("windows_lifetime_peak_working_set_bytes"),
            "outer_fit_seconds": audit["fit_seconds"], "outer_prediction_seconds": audit["last_predict_seconds"]})
    expected = set(itertools.product(config["drug_ids"], config["models"], range(config["outer_repeats"]), range(config["outer_folds"])))
    if seen != expected:
        raise ValueError("Complete repeated outer-audit coverage required")
    stability_rows = []
    for (drug, model), records in features.items():
        for left, right in itertools.combinations(sorted(records, key=lambda x: x[:2]), 2):
            ra, fa, ga, aa, ta = left
            rb, fb, gb, ab, tb = right
            gene_union, active_union = ga | gb, set(aa) | set(ab)
            common = set(aa) & set(ab)
            stability_rows.append({"gdsc_drug_id": drug, "model": model, "repeat_a": ra, "outer_fold_a": fa,
                "repeat_b": rb, "outer_fold_b": fb,
                "comparison_scope": "same_repeat_different_fold" if ra == rb else "different_repeats",
                "training_sample_jaccard": len(ta & tb) / len(ta | tb),
                "candidate_gene_jaccard": len(ga & gb) / len(gene_union) if gene_union else np.nan,
                "active_feature_jaccard": len(common) / len(active_union) if active_union else np.nan,
                "same_sign_active_jaccard": sum(np.sign(aa[name]) == np.sign(ab[name]) for name in common) / len(active_union) if active_union else np.nan,
                "empty_active_union": not bool(active_union),
                "scope": "overlapping_training_sets_descriptive_not_independent_validation"})
    resource_frame, chosen_frame = pd.DataFrame(resource_rows), pd.DataFrame(choices)
    resource_summary = resource_frame.groupby(["gdsc_drug_id", "model", "repeat"], sort=False).agg(
        outer_tasks=("outer_fold", "size"), task_wall_seconds_sum=("task_wall_seconds", "sum"),
        task_cpu_seconds_sum=("task_cpu_seconds", "sum"), fit_work_wall_seconds_sum=("fit_work_wall_seconds", lambda x: x.sum(min_count=1)),
        initialization_wall_seconds_sum=("initialization_wall_seconds", lambda x: x.sum(min_count=1)),
        sampled_peak_rss_bytes_max=("sampled_peak_rss_bytes", "max"),
        windows_lifetime_peak_working_set_bytes_max=("windows_lifetime_peak_working_set_bytes", "max"),
        outer_fit_seconds_sum=("outer_fit_seconds", "sum"), outer_prediction_seconds_sum=("outer_prediction_seconds", "sum")).reset_index()
    boundary_rows = []
    for (drug, model), frame in chosen_frame.groupby(["gdsc_drug_id", "model"], sort=False):
        for params, selected in frame.groupby("chosen_params", sort=False):
            boundary_rows.append({"gdsc_drug_id": drug, "model": model, "chosen_params": params, "outer_fits": len(selected),
                "selected_boundary_outer_fits": int(selected.chosen_grid_boundary.sum()),
                "outer_solver_iterations_max": selected.outer_solver_iterations.max(),
                "zero_molecular_outer_fits": int(selected.active_features.eq(0).sum()) if model != "tissue_null" else 0,
                "scope": "null_endpoint_is_predefined_not_automatically_a_grid_defect"})
    return {"molecular_and_failure_summary": pd.DataFrame(states.values()), "chosen_parameters": chosen_frame,
        "parameter_boundary_summary": pd.DataFrame(boundary_rows), "feature_stability": pd.DataFrame(stability_rows),
        "resource_tasks": resource_frame, "resource_summary": resource_summary}


def aggregate_fit_events(path, config):
    states = {}
    for line in Path(path).open(encoding="utf-8"):
        if not line.strip():
            continue
        event = json.loads(line)
        key = (int(event["drug_id"]), event["model"], int(event["repeat"]))
        state = states.setdefault(key, {"gdsc_drug_id": key[0], "model": key[1], "repeat": key[2], "actual_fit_events": 0,
            "all_fit_warning_events": 0, "all_fit_fallback_events": 0, "fit_seconds_sum": 0., "prediction_seconds_sum": 0.})
        state["actual_fit_events"] += 1
        state["all_fit_warning_events"] += int(bool(event.get("warnings")))
        state["all_fit_fallback_events"] += int(bool(event.get("fallback_reason")))
        state["fit_seconds_sum"] += event["fit_seconds"]
        state["prediction_seconds_sum"] += event["predict_seconds"]
    expected = set(itertools.product(config["drug_ids"], config["models"], range(config["outer_repeats"])))
    if set(states) != expected:
        raise ValueError("Missing/unplanned actual-fit group")
    for key, row in states.items():
        expected_count = config["outer_folds"] * (1 if key[1] == "tissue_null" else config["inner_folds"] * len(config["elastic_net"]["alpha_ratios"]) * len(config["elastic_net"]["l1_ratios"]) + 1)
        if row["actual_fit_events"] != expected_count:
            raise ValueError("Actual fits disagree with frozen per-model grid")
    return pd.DataFrame(states.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/ordinal_repeated_v1.json")
    parser.add_argument("--run-dir")
    parser.add_argument("--output-dir", default="results/ordinal_repeated")
    args = parser.parse_args()
    config_path = (ROOT / args.config).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    run = (ROOT / (args.run_dir or "runs/" + config["run_name"])).resolve()
    output = (ROOT / args.output_dir).resolve()
    for path in (config_path, run, output):
        path.relative_to(ROOT)
    if output.exists() or (run / "statistics_private").exists():
        raise FileExistsError("Existing repeated statistics retained; use an explicit new output version")
    blocked = install_response_source_guard()
    started = time.perf_counter()
    metadata_path = run / "ordinal_repeated_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    lock_path = ROOT / "config/ordinal_repeated_v1.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if metadata.get("completed") is not True or metadata.get("status") != "completed" or metadata.get("completed_tasks") != lock["contract"]["tasks"]:
        raise RuntimeError("Completed repeated evaluation required; pause/partial results are not analyzed")
    if metadata.get("external_response_accessed") is not False or metadata.get("data_guard", {}).get("blocked_attempts") != 0:
        raise RuntimeError("External-response seal or guard failed")
    if metadata.get("config_sha256") != sha(config_path) or metadata.get("freeze_sha256") != sha(lock_path):
        raise RuntimeError("Repeated config/lock binding changed")
    for name, checksum in lock["files"].items():
        if sha(ROOT / name) != checksum:
            raise RuntimeError("Frozen file changed: " + name)
    validation_path = run / "ordinal_repeated_validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    validation_lock = validation.get("repeated_lock_sha256", validation.get("diagnostic_lock_sha256"))
    if validation.get("status") != "passed" or validation_lock != sha(lock_path):
        raise RuntimeError("Independent completed-run validation required")
    for name in ("oof_predictions.csv", "outer_fit_audits.jsonl", "fit_events.jsonl"):
        if metadata["output_hashes"].get(name) != sha(run / name):
            raise RuntimeError("Run output changed: " + name)
    manifest_path = ROOT / "results/cohort_v2/development_sample_manifest.csv"
    manifest = pd.read_csv(manifest_path, dtype={"canonical_model_id": str, "biological_group_id": str})
    data = pd.read_csv(run / "oof_predictions.csv", dtype={"canonical_model_id": str, "biological_group_id": str}, float_precision="round_trip")
    with threadpool_limits(limits=1):
        result = analyze_predictions(data, manifest, config)
    result.update(aggregate_audits(run / "outer_fit_audits.jsonl", config))
    result["all_fit_events_summary"] = aggregate_fit_events(run / "fit_events.jsonl", config)
    fits = int(result["all_fit_events_summary"].actual_fit_events.sum())
    if fits != lock["contract"]["fits"] or len(data) != lock["contract"]["oof_predictions"]:
        raise RuntimeError("Complete fit/OOF contract mismatch")
    overhead = validation.get("operational_overhead")
    if not isinstance(overhead, dict) or overhead.get("successful_pipeline_fits") != fits:
        raise RuntimeError("Operational audit missing or inconsistent")
    output.mkdir(parents=True, exist_ok=False)
    private = run / "statistics_private"
    private.mkdir(exist_ok=False)
    groups, counts = result.pop("bootstrap_group_order"), result.pop("bootstrap_counts")
    pd.DataFrame({"bootstrap_group_position": np.arange(len(groups)), "biological_group_id": groups}).to_csv(private / "bootstrap_group_manifest.csv", index=False, lineterminator="\n")
    np.savez_compressed(private / "bootstrap_group_multiplicities.npz", multiplicities=counts)
    for name, frame in result.items():
        frame.to_csv(output / (name + ".csv"), index=False, encoding="utf-8-sig", lineterminator="\n")
    write_json(output / "operational_overhead.json", overhead)
    memory = psutil.Process().memory_info()
    provenance = {"status": "completed", "created_utc": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": time.perf_counter() - started, "process_lifetime_peak_rss_bytes": getattr(memory, "peak_wset", memory.rss),
        "source_run": run.relative_to(ROOT).as_posix(), "execution_commit": metadata.get("execution_commit"),
        "analyzer_sha256": sha(Path(__file__)), "config_sha256": sha(config_path), "lock_sha256": sha(lock_path),
        "metadata_sha256": sha(metadata_path), "validation_sha256": sha(validation_path),
        "OOF_source_sha256": sha(run / "oof_predictions.csv"), "development_manifest_sha256": sha(manifest_path),
        "n_global_biological_groups": len(groups), "n_observed_biological_groups": data.biological_group_id.nunique(),
        "bootstrap_draws": len(counts), "bootstrap_seed": config["statistics"]["bootstrap_seed"],
        "bootstrap_count_matrix_sha256": hashlib.sha256(counts.tobytes(order="C")).hexdigest(),
        "bootstrap_group_order_sha256": hashlib.sha256("\0".join(groups).encode()).hexdigest(),
        "expected_outer_tasks": lock["contract"]["tasks"], "actual_fit_events": fits, "oof_rows": len(data),
        "outer_repeats": config["outer_repeats"], "no_repeat_N_inflation": True,
        "bootstrap_shared_across_models_drugs_repeats": True, "new_model_fits": 0,
        "confidence_scope": config["statistics"]["confidence_scope"], "noninferiority_or_equivalence_tested": False,
        "primary_family": "9 Spearman effects: association pair versus strong expression; separate from 18 representation effects",
        "representation_family": "18 Spearman effects: pair-expression within each rule at eight drugs plus macro",
        "no_joint_family_guarantee_or_confirmatory_P_values": True,
        "new_permutation_runs": 0, "external_response_accessed": False, "blocked_source_access_attempts": blocked,
        "output_sha256": {**{name + ".csv": sha(output / (name + ".csv")) for name in result},
                          "operational_overhead.json": sha(output / "operational_overhead.json")}}
    write_json(output / "analysis_provenance.json", provenance)
    print(json.dumps(clean(provenance), ensure_ascii=False))


if __name__ == "__main__":
    main()
