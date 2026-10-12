"""Fixed-protocol paired inference for formal development OOF records.

All ranks/metrics are recomputed under shared global BIOGROUP multiplicities.
Three repeats remain three complete prediction sets, never three times the N.
Intervals condition on observed OOF records; no fit/tuning is rerun here.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
METRICS = ("spearman", "rmse", "mae", "r2")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def dump(path, value):
    Path(path).write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def generate_global_bootstrap_weights(group_ids, draws, seed):
    """Return lexicographic global group order and one shared count matrix."""
    groups = np.asarray(sorted(set(map(str, group_ids))), dtype=str)
    if len(groups) < 1:
        raise ValueError("Global bootstrap needs at least one biological group")
    rng = np.random.Generator(np.random.PCG64(seed))
    counts = rng.multinomial(len(groups), np.full(len(groups), 1.0 / len(groups)), size=draws)
    dtype = np.uint16 if len(groups) <= np.iinfo(np.uint16).max else np.uint32
    return groups, counts.astype(dtype)


def weighted_midrank(values, weights):
    """Exact ranks of duplicated observations represented by integer weights.

    A zero-weight source row has no effect. Exact tied source values share the
    average rank of all bootstrap copies in their entire value group.
    """
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.ndim != 1 or weights.ndim != 2 or weights.shape[1] != len(values):
        raise ValueError("Expected values[N] and multiplicities[B,N]")
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or (weights < 0).any():
        raise ValueError("Ranks require finite values and nonnegative multiplicities")
    if not np.equal(weights, np.floor(weights)).all():
        raise ValueError("Bootstrap multiplicities must be integer counts")
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    start = np.r_[0, np.flatnonzero(sorted_values[1:] != sorted_values[:-1]) + 1]
    value_group = np.cumsum(np.r_[True, sorted_values[1:] != sorted_values[:-1]]) - 1
    sums = np.add.reduceat(weights[:, order], start, axis=1)
    mid = np.cumsum(sums, axis=1) - 0.5 * sums + 0.5
    ranks_sorted = mid[:, value_group]
    ranks = np.empty_like(ranks_sorted)
    ranks[:, order] = ranks_sorted
    return ranks


def weighted_metrics(y, prediction, weights, ranking_enabled=True):
    """Compute exact duplicate-row bootstrap metrics without duplicating rows."""
    y = np.asarray(y, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    if y.ndim != 1 or prediction.shape != y.shape or w.ndim != 2 or w.shape[1] != len(y):
        raise ValueError("Invalid paired metric shapes")
    result = {name: np.full(len(w), np.nan) for name in METRICS}
    if not np.isfinite(y).all() or not np.isfinite(prediction).all():
        return result
    n = w.sum(axis=1)
    nonempty = n > 0
    error = prediction - y
    squared = w @ np.square(error)
    with np.errstate(divide="ignore", invalid="ignore"):
        result["rmse"][nonempty] = np.sqrt(squared[nonempty] / n[nonempty])
        result["mae"][nonempty] = (w @ np.abs(error))[nonempty] / n[nonempty]
    ry = weighted_midrank(y, w)
    center = (n + 1.0) / 2.0
    cy = ry - center[:, None]
    rank_y_ss = np.sum(w * cy * cy, axis=1)
    variable_y = rank_y_ss > 0
    # Center first to reduce catastrophic cancellation for near-constant AUC.
    y_centered = y - y.mean()
    with np.errstate(divide="ignore", invalid="ignore"):
        yc_mean = (w @ y_centered) / n
        ss_y = w @ np.square(y_centered) - n * np.square(yc_mean)
        r2_valid = nonempty & variable_y & (ss_y > 0)
        result["r2"][r2_valid] = 1.0 - squared[r2_valid] / ss_y[r2_valid]
    if ranking_enabled:
        rp = weighted_midrank(prediction, w)
        cp = rp - center[:, None]
        rank_p_ss = np.sum(w * cp * cp, axis=1)
        valid = (n >= 2) & variable_y & (rank_p_ss > 0)
        cross = np.sum(w * cy * cp, axis=1)
        result["spearman"][valid] = cross[valid] / np.sqrt(rank_y_ss[valid] * rank_p_ss[valid])
    return result


def interval_summary(point, samples, configuration, family_size=None, reason=""):
    """Explicit validity accounting; never turn bootstrap tails into P-values."""
    samples = np.asarray(samples, dtype=np.float64)
    valid = np.isfinite(samples)
    fraction = float(valid.mean()) if len(samples) else 0.0
    point_defined = bool(np.isfinite(point))
    ok = point_defined and fraction >= configuration["bootstrap"]["minimum_valid_draw_fraction"]
    output = {"point": point, "bootstrap_requested": len(samples), "bootstrap_valid": int(valid.sum()),
        "bootstrap_invalid": int((~valid).sum()), "bootstrap_valid_fraction": fraction,
        "conditional_95_lower": np.nan, "conditional_95_upper": np.nan,
        "family_adjusted_conditional_lower": np.nan, "family_adjusted_conditional_upper": np.nan,
        "family_size": family_size, "ci_status": "reported" if ok else "withheld_undefined_point" if not point_defined else "withheld_insufficient_valid_draws",
        "undefined_reason": reason if not point_defined else "" if valid.all() else "One or more fixed-repeat/paired metrics undefined in these resamples; retained invalid-count diagnostic.",
        "interval_scope": "conditional_on_observed_OOF_records_not_full_train_refit_or_guaranteed_algorithmic_CV_coverage"}
    if ok:
        alpha = 1.0 - configuration["bootstrap"]["unadjusted_level"]
        output["conditional_95_lower"], output["conditional_95_upper"] = np.quantile(samples[valid], [alpha / 2.0, 1.0 - alpha / 2.0], method="linear")
        if family_size is not None:
            tail = configuration["multiplicity"]["family_alpha"] / (2 * family_size)
            output["family_adjusted_conditional_lower"], output["family_adjusted_conditional_upper"] = np.quantile(samples[valid], [tail, 1 - tail], method="linear")
    return output


def analyze_oof(oof, manifest, configuration, weights=None):
    """Pure-data API for independent synthetic tests and the formal CLI.

    Returns public aggregate frames plus shared bootstrap weights/group order.
    Input OOF needs drug/model/repeat/outer_fold/canonical/group/y/prediction.
    """
    data = oof.copy()
    manifest = manifest.copy()
    if "drug_id" in data and "gdsc_drug_id" not in data:
        data = data.rename(columns={"drug_id": "gdsc_drug_id"})
    required = {"gdsc_drug_id", "model", "repeat", "outer_fold", "canonical_model_id", "biological_group_id", "observed_auc", "predicted_auc"}
    if not required <= set(data.columns):
        raise ValueError("Missing formal OOF fields: " + ",".join(sorted(required - set(data.columns))))
    for frame in (data, manifest):
        frame["gdsc_drug_id"] = frame.gdsc_drug_id.astype(int)
        frame["canonical_model_id"] = frame.canonical_model_id.astype(str)
        frame["biological_group_id"] = frame.biological_group_id.astype(str)
    data["repeat"] = data["repeat"].astype(int)
    repeats = range(configuration["outer_repeats"])
    drugs = configuration["drug_ids_in_frozen_rank"]
    models = configuration["models"]
    if data.duplicated(["gdsc_drug_id", "model", "repeat", "canonical_model_id"]).any():
        raise ValueError("Duplicate OOF identity in a complete repeat")
    if not set(data.model) <= set(models) or not set(data.gdsc_drug_id) <= set(drugs) or not set(data["repeat"]) <= set(repeats):
        raise ValueError("Unplanned model/drug/repeat in formal OOF input")
    if weights is None:
        groups, counts = generate_global_bootstrap_weights(manifest.biological_group_id, configuration["bootstrap"]["draws"], configuration["bootstrap"]["seed"])
    else:
        groups = np.asarray(sorted(set(manifest.biological_group_id)), dtype=str)
        counts = np.asarray(weights)
    if counts.ndim != 2 or counts.shape[1] != len(groups) or counts.shape[0] != configuration["bootstrap"]["draws"]:
        raise ValueError("Shared bootstrap count matrix has the wrong shape")
    if (counts < 0).any() or not np.equal(counts, np.floor(counts)).all() or not (counts.sum(axis=1) == len(groups)).all():
        raise ValueError("Invalid shared global group multiplicities")
    group_position = {group: i for i, group in enumerate(groups)}
    by_model = {}
    points = {}
    repeat_rows = []
    status_reasons = {}
    chunk_size = configuration["bootstrap"]["computational_block_draws"]
    b = len(counts)
    for drug in drugs:
        source = manifest[manifest.gdsc_drug_id.eq(drug)].sort_values("canonical_model_id")
        if source.empty or source.canonical_model_id.duplicated().any():
            raise ValueError("Expected source manifest must contain unique per-drug canonical identities")
        ids = source.canonical_model_id.tolist()
        source_groups = source.set_index("canonical_model_id").biological_group_id
        subset = data[data.gdsc_drug_id.eq(drug)]
        if not set(subset.canonical_model_id) <= set(ids):
            raise ValueError("OOF contains a model absent from frozen drug manifest")
        if not subset.biological_group_id.eq(subset.canonical_model_id.map(source_groups)).all():
            raise ValueError("OOF biological group disagrees with frozen identity")
        invariant = subset.groupby("canonical_model_id").observed_auc.nunique(dropna=True)
        if (invariant > 1).any():
            raise ValueError("One frozen drug-model outcome changed between predictions/repeats")
        y = subset.groupby("canonical_model_id").observed_auc.first().reindex(ids).to_numpy(dtype=np.float64)
        positions = np.asarray([group_position[group] for group in source.biological_group_id])
        for model in models:
            repeat_predictions = []
            point_matrix = np.full((len(repeats), len(METRICS)), np.nan)
            boot_matrix = np.full((len(repeats), b, len(METRICS)), np.nan)
            reasons = []
            for repeat in repeats:
                record = subset[subset.model.eq(model) & subset["repeat"].eq(repeat)].set_index("canonical_model_id")
                p = record.predicted_auc.reindex(ids).to_numpy(dtype=np.float64)
                missing_rows = len(ids) - len(record)
                reason = "incomplete_repeat" if missing_rows else "nonfinite_response_or_prediction" if not np.isfinite(y).all() or not np.isfinite(p).all() else ""
                values = weighted_metrics(y, p, np.ones((1, len(ids)), dtype=int), ranking_enabled=model != "mean")
                for j, metric in enumerate(METRICS):
                    point_matrix[repeat, j] = values[metric][0]
                rank_reason = reason
                if not rank_reason:
                    if model == "mean":
                        rank_reason = "training_mean_has_no_individual_ranking"
                    elif np.isfinite(y).all() and np.ptp(y) == 0:
                        rank_reason = "constant_response"
                    elif np.isfinite(p).all() and np.ptp(p) == 0:
                        rank_reason = "constant_prediction"
                repeat_rows.append({"gdsc_drug_id": drug, "model": model, "repeat": repeat,
                    "n_expected_models": len(ids), "n_unique_biological_groups": source.biological_group_id.nunique(), "missing_prediction_rows": missing_rows,
                    "repeat_status": "complete" if not reason else reason, "spearman_unavailable_reason": rank_reason if not np.isfinite(point_matrix[repeat, 0]) else "",
                    **{name: point_matrix[repeat, j] for j, name in enumerate(METRICS)}})
                reasons.append(rank_reason)
                repeat_predictions.append(p)
            for start in range(0, b, chunk_size):
                stop = min(start + chunk_size, b)
                shared = counts[start:stop, positions]
                for repeat, p in enumerate(repeat_predictions):
                    values = weighted_metrics(y, p, shared, ranking_enabled=model != "mean")
                    for j, metric in enumerate(METRICS):
                        boot_matrix[repeat, start:stop, j] = values[metric]
            points[(drug, model)] = point_matrix
            by_model[(drug, model)] = boot_matrix
            status_reasons[(drug, model)] = "|".join(sorted(set(reason for reason in reasons if reason)))
    model_rows = []
    stability_rows = []
    for drug in drugs:
        n_groups = manifest[manifest.gdsc_drug_id.eq(drug)].biological_group_id.nunique()
        for model in models:
            p = points[(drug, model)]
            bootstrap_mean = np.mean(by_model[(drug, model)], axis=0)
            for j, metric in enumerate(METRICS):
                point = np.mean(p[:, j])
                reason = status_reasons[(drug, model)] or "undefined_fixed_repeat_metric"
                model_rows.append({"gdsc_drug_id": drug, "model": model, "metric": metric, "n_unique_biological_groups": n_groups,
                    **interval_summary(point, bootstrap_mean[:, j], configuration, reason=reason)})
                stability_rows.append({"gdsc_drug_id": drug, "comparison": "model:" + model, "metric": metric,
                    "repeat_values": json.dumps(clean(p[:, j].tolist()), allow_nan=False), "repeat_mean": point,
                    "repeat_min": np.min(p[:, j]), "repeat_max": np.max(p[:, j]), "repeat_range": np.ptp(p[:, j]), "repeat_sample_sd": np.std(p[:, j], ddof=1)})
    for model in models:
        for j, metric in enumerate(METRICS):
            macro_repeats = np.mean(np.stack([points[(drug, model)][:, j] for drug in drugs]), axis=0)
            macro_boot = np.mean(np.stack([np.mean(by_model[(drug, model)][:, :, j], axis=0) for drug in drugs]), axis=0)
            bad_drugs = [drug for drug in drugs if not np.isfinite(points[(drug, model)][:, j]).all()]
            reason = "undefined_constituent_drugs:" + ",".join(map(str, bad_drugs)) if bad_drugs else ""
            model_rows.append({"gdsc_drug_id": "macro_equal_fixed_drugs", "model": model, "metric": metric,
                "n_unique_biological_groups": len(groups), **interval_summary(np.mean(macro_repeats), macro_boot, configuration, reason=reason)})
            stability_rows.append({"gdsc_drug_id": "macro_equal_fixed_drugs", "comparison": "model:" + model, "metric": metric,
                "repeat_values": json.dumps(clean(macro_repeats.tolist()), allow_nan=False), "repeat_mean": np.mean(macro_repeats),
                "repeat_min": np.min(macro_repeats), "repeat_max": np.max(macro_repeats), "repeat_range": np.ptp(macro_repeats), "repeat_sample_sd": np.std(macro_repeats, ddof=1)})
    comparisons = [dict(name="primary_EN", **configuration["primary_comparison"]), dict(name="secondary_RF", **configuration["secondary_comparison"])] + configuration["auxiliary_comparisons"]
    effect_rows = []
    for comparison in comparisons:
        name, candidate, reference = comparison["name"], comparison["candidate"], comparison["reference"]
        effect_points = {}
        effect_boot = {}
        effect_repeat_values = {}
        for drug in drugs:
            a, r = points[(drug, candidate)], points[(drug, reference)]
            ba, br = by_model[(drug, candidate)], by_model[(drug, reference)]
            for j, metric in enumerate(METRICS):
                sign = -1 if metric in ["rmse", "mae"] else 1
                values = sign * (a[:, j] - r[:, j])
                boot_values = sign * (ba[:, :, j] - br[:, :, j])
                effect_points[(drug, metric)] = np.mean(values)
                effect_boot[(drug, metric)] = np.mean(boot_values, axis=0)
                effect_repeat_values[(drug, metric)] = values
                stability_rows.append({"gdsc_drug_id": drug, "comparison": name, "metric": metric,
                    "repeat_values": json.dumps(clean(values.tolist()), allow_nan=False), "repeat_mean": np.mean(values),
                    "repeat_min": np.min(values), "repeat_max": np.max(values), "repeat_range": np.ptp(values), "repeat_sample_sd": np.std(values, ddof=1)})
            with np.errstate(divide="ignore", invalid="ignore"):
                ratios = (r[:, 1] - a[:, 1]) / r[:, 1]
                boot_ratios = (br[:, :, 1] - ba[:, :, 1]) / br[:, :, 1]
            ratios[~np.isfinite(ratios)] = np.nan
            boot_ratios[~np.isfinite(boot_ratios)] = np.nan
            effect_points[(drug, "relative_rmse_improvement")] = np.mean(ratios)
            effect_boot[(drug, "relative_rmse_improvement")] = np.mean(boot_ratios, axis=0)
            effect_repeat_values[(drug, "relative_rmse_improvement")] = ratios
            stability_rows.append({"gdsc_drug_id": drug, "comparison": name, "metric": "relative_rmse_improvement",
                "repeat_values": json.dumps(clean(ratios.tolist()), allow_nan=False), "repeat_mean": np.mean(ratios),
                "repeat_min": np.min(ratios), "repeat_max": np.max(ratios), "repeat_range": np.ptp(ratios), "repeat_sample_sd": np.std(ratios, ddof=1)})
        for metric in (*METRICS, "relative_rmse_improvement"):
            for drug in [*drugs, "macro_equal_fixed_drugs"]:
                family = configuration["multiplicity"]["primary_family_size"] if name == "primary_EN" and metric == "spearman" else configuration["multiplicity"]["secondary_family_size"] if name == "secondary_RF" and metric == "spearman" else None
                if isinstance(drug, int):
                    point, boot = effect_points[(drug, metric)], effect_boot[(drug, metric)]
                    unknown = [drug] if not np.isfinite(point) else []
                    n_groups = manifest[manifest.gdsc_drug_id.eq(drug)].biological_group_id.nunique()
                else:
                    point = np.mean([effect_points[(d, metric)] for d in drugs])
                    boot = np.mean(np.stack([effect_boot[(d, metric)] for d in drugs]), axis=0)
                    unknown = [d for d in drugs if not np.isfinite(effect_points[(d, metric)])]
                    n_groups = len(groups)
                    macro_repeats = np.mean(np.stack([effect_repeat_values[(d, metric)] for d in drugs]), axis=0)
                    stability_rows.append({"gdsc_drug_id": drug, "comparison": name, "metric": metric,
                        "repeat_values": json.dumps(clean(macro_repeats.tolist()), allow_nan=False), "repeat_mean": np.mean(macro_repeats),
                        "repeat_min": np.min(macro_repeats), "repeat_max": np.max(macro_repeats), "repeat_range": np.ptp(macro_repeats), "repeat_sample_sd": np.std(macro_repeats, ddof=1)})
                reason = "undefined_constituent_drugs:" + ",".join(map(str, unknown)) if unknown else ""
                effect_rows.append({"gdsc_drug_id": drug, "comparison": name, "candidate": candidate, "reference": reference, "metric": metric,
                    "n_unique_biological_groups": n_groups, "fixed_drug_count": len(drugs) if isinstance(drug, str) else 1,
                    **interval_summary(point, boot, configuration, family_size=family, reason=reason)})
    return {"repeat_metrics": pd.DataFrame(repeat_rows), "model_intervals": pd.DataFrame(model_rows),
        "paired_effect_intervals": pd.DataFrame(effect_rows), "repeat_stability": pd.DataFrame(stability_rows),
        "bootstrap_group_order": groups, "bootstrap_counts": counts,
        "undefined_repeat_reasons": dict(Counter(row["spearman_unavailable_reason"] for row in repeat_rows if row["spearman_unavailable_reason"]))}


def install_response_source_guard():
    attempts = []
    def hook(event, args):
        if event != "open" or not args or not isinstance(args[0], (str, bytes, Path)):
            return
        path = Path(os.fsdecode(args[0])).resolve()
        forbidden = path.is_relative_to(ROOT / "data_raw") or path.suffix == ".npy" or "ctrp" in path.name.casefold() or path.name.startswith("GDSC2_fitted_dose_response")
        if forbidden:
            attempts.append(path.name)
            raise PermissionError("Formal OOF analyzer does not read raw responses, CTRP sources or expression arrays")
    sys.addaudithook(hook)
    return attempts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--statistics-config", default="config/formal_statistics_v1.json")
    parser.add_argument("--output-dir", default="results/formal_development")
    args = parser.parse_args()
    def local(value):
        p = Path(value)
        p = p.resolve() if p.is_absolute() else (ROOT / p).resolve()
        p.relative_to(ROOT)
        return p
    run, configuration_path, output = map(local, [args.run_dir, args.statistics_config, args.output_dir])
    blocked = install_response_source_guard()
    configuration = json.loads(configuration_path.read_text(encoding="utf-8"))
    protocol_path = ROOT / "docs/FORMAL_ANALYSIS_PROTOCOL.md"
    meta = json.loads((run / "formal_metadata.json").read_text(encoding="utf-8"))
    if meta.get("completed") is not True or meta.get("completed_tasks") != meta.get("expected_tasks"):
        raise RuntimeError("Completed formal development run required; no partial-result selection")
    assert meta.get("external_response_accessed") is False, "Formal metadata must explicitly record no external response access"
    assert meta.get("data_guard", {}).get("blocked_attempts") == 0, "Explicit zero external source access attempts required"
    assert meta.get("expected_tasks") == len(configuration["drug_ids_in_frozen_rank"]) * len(configuration["models"]) * configuration["outer_repeats"] * configuration["outer_folds"]
    assert meta.get("statistics_sha256") == sha(configuration_path), "Formal run must bind the frozen statistics config hash"
    assert meta.get("protocol_sha256") == sha(protocol_path), "Formal run must bind the pre-run protocol hash"
    cohort_lock = json.loads((ROOT / "config/cohort_freeze_v2.lock.json").read_text(encoding="utf-8"))
    assert meta.get("cohort_sha256") == cohort_lock["freeze_sha256"]
    assert sha(ROOT / "config/cohort_freeze_v2.json") == cohort_lock["freeze_sha256"]
    formal_lock_path = ROOT / "config/formal_analysis_v1.lock.json"
    assert meta.get("freeze_sha256") == sha(formal_lock_path), "Formal statistical/model lock no longer matches run"
    formal_lock = json.loads(formal_lock_path.read_text(encoding="utf-8"))
    for path in [configuration_path, protocol_path, Path(__file__).resolve()]:
        name = path.relative_to(ROOT).as_posix()
        assert formal_lock["files"].get(name) == sha(path), "Frozen analysis implementation or protocol changed: " + name
    manifest_path = ROOT / "results/cohort_v2/development_sample_manifest.csv"
    assert formal_lock["files"].get(manifest_path.relative_to(ROOT).as_posix()) == sha(manifest_path), "Frozen development identity manifest changed"
    manifest = pd.read_csv(manifest_path, dtype={"canonical_model_id": str, "biological_group_id": str})
    oof_path = run / "oof_predictions.csv"
    audits_path = run / "outer_fit_audits.jsonl"
    assert meta.get("output_hashes", {}).get("oof_predictions.csv") == sha(oof_path), "Formal OOF output hash mismatch"
    assert meta.get("output_hashes", {}).get("outer_fit_audits.jsonl") == sha(audits_path), "Formal audit output hash mismatch"
    oof = pd.read_csv(oof_path, dtype={"canonical_model_id": str, "biological_group_id": str})
    started = time.perf_counter()
    with threadpool_limits(limits=1):
        result = analyze_oof(oof, manifest, configuration)
    output.mkdir(parents=True, exist_ok=False)
    for name in ["repeat_metrics", "model_intervals", "paired_effect_intervals", "repeat_stability"]:
        result[name].to_csv(output / (name + ".csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame({"bootstrap_group_position": np.arange(len(result["bootstrap_group_order"])), "biological_group_id": result["bootstrap_group_order"]}).to_csv(output / "bootstrap_group_manifest.csv", index=False, encoding="utf-8-sig")
    np.savez_compressed(output / "bootstrap_group_multiplicities.npz", multiplicities=result["bootstrap_counts"])
    if audits_path.exists():
        aggregate = {}
        # Stream because each task may carry a large complete feature audit.
        with audits_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                event = json.loads(line)
                task = event["task"]
                row = event.get("outer_fit_audit", event)
                key = (int(task["drug_id"]), task["model"], int(task["repeat"]))
                state = aggregate.setdefault(key, {"gdsc_drug_id": key[0], "model": key[1], "repeat": key[2], "outer_fits": 0, "zero_molecular_coefficients_outer_fits": 0, "warning_outer_fits": 0, "fallback_outer_fits": 0})
                state["outer_fits"] += 1
                state["zero_molecular_coefficients_outer_fits"] += int(row.get("n_nonzero_molecular_coeffs") == 0)
                state["warning_outer_fits"] += int(bool(row.get("warnings")))
                state["fallback_outer_fits"] += int(bool(row.get("fallback_reason")))
        pd.DataFrame(aggregate.values()).to_csv(output / "molecular_and_failure_audit_summary.csv", index=False, encoding="utf-8-sig")
    counts = result["bootstrap_counts"]
    public = {"status": "completed", "statistics_protocol_version": configuration["protocol_version"],
        "created_utc": datetime.now(timezone.utc).isoformat(), "wall_seconds": time.perf_counter() - started,
        "source_run": run.relative_to(ROOT).as_posix(), "code_execution_commit": meta.get("execution_commit"),
        "analyzer_sha256": sha(Path(__file__)), "statistics_config_sha256": sha(configuration_path), "statistics_protocol_sha256": sha(protocol_path),
        "formal_analysis_lock_sha256": meta["freeze_sha256"], "cohort_freeze_sha256": meta["cohort_sha256"],
        "OOF_source_sha256": sha(oof_path), "development_manifest_sha256": sha(manifest_path),
        "outer_fit_audits_source_sha256": sha(audits_path),
        "n_global_biological_groups": len(result["bootstrap_group_order"]), "bootstrap_draws": len(counts),
        "bootstrap_seed": configuration["bootstrap"]["seed"], "bootstrap_generator": "PCG64.multinomial",
        "bootstrap_count_matrix_sha256": hashlib.sha256(counts.tobytes(order="C")).hexdigest(),
        "bootstrap_group_order_sha256": hashlib.sha256("\0".join(result["bootstrap_group_order"]).encode()).hexdigest(),
        "numpy_version": np.__version__, "pandas_version": pd.__version__, "no_repeat_N_inflation": True,
        "bootstrap_shared_across_models_repeats_drugs": True, "undefined_repeat_reasons": result["undefined_repeat_reasons"],
        "external_response_accessed": False, "blocked_source_access_attempts": blocked,
        "interval_scope": configuration["bootstrap"]["interval_scope"], "multiplicity": configuration["multiplicity"],
        "failure_task_policy": configuration["undefined_policy"], "no_confirmatory_P_or_equivalence_claims": True}
    dump(output / "analysis_provenance.json", public)
    print(json.dumps({"status": "completed", "output": output.relative_to(ROOT).as_posix(), "wall_seconds": public["wall_seconds"], "bootstrap_draws": len(counts)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
