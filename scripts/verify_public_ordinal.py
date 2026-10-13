"""Verify published ordinal diagnostic hashes and aggregate consistency.

Standard-library only; never opens private OOF, checkpoints, response data or
live metadata. This verifies the public snapshot, not the private statistical
calculation or the truth of scientific claims. Mutable reports are excluded.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
MACRO = "macro_equal_four_diagnostic_cases"
METRICS = ("spearman", "rmse", "mae", "r2")
FILES = {
    "point_metrics.csv", "model_intervals.csv", "paired_effect_intervals.csv",
    "null_contrasts.csv", "null_contrast_summary.csv", "molecular_and_failure_summary.csv",
    "chosen_parameters.csv", "active_feature_stability.csv", "resource_tasks.csv",
    "resource_summary.csv", "all_fit_events_summary.csv", "permutation_coverage.csv",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def number(value):
    return float(value) if value not in ("", None) else math.nan


def integer(value):
    result = number(value)
    require(math.isfinite(result) and result.is_integer(), "Expected an integer count")
    return int(result)


def close(actual, expected, message):
    actual, expected = number(actual), number(expected)
    require((math.isnan(actual) and math.isnan(expected)) or math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12), message)


def table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames is not None and len(reader.fieldnames) == len(set(reader.fieldnames)), "Missing/duplicate CSV header: " + path.name)
        require(not any(name in reader.fieldnames for name in ("canonical_model_id", "biological_group_id", "observed_auc", "predicted_auc")), "Individual-level column in public aggregate: " + path.name)
        result = list(reader)
    require(result and all(None not in row and None not in row.values() for row in result), "Empty or malformed CSV: " + path.name)
    return result


def index_rows(rows, columns, expected, label):
    indexed = {tuple(row[name] for name in columns): row for row in rows}
    require(len(indexed) == len(rows), "Duplicate aggregate key: " + label)
    require(set(indexed) == expected, "Missing/unplanned aggregate coverage: " + label)
    return indexed


def contrasts(config):
    result = []
    for rule in config["selection_rules"]:
        for candidate, reference in (("rank", "expression"), ("pair_full", "rank"), ("pair_full", "expression")):
            result.append((rule + ":" + candidate + "-" + reference, rule + "_" + candidate, rule + "_" + reference))
    for representation in config["representations"]:
        result.append((representation + ":association-MAD", "response_association_" + representation, "mad_" + representation))
    result.extend((model + "-tissue_null", model, "tissue_null") for model in config["models"] if model != "tissue_null")
    require(len(result) == 15, "Frozen contrast count changed")
    return result


def check_intervals(rows, draws, minimum_fraction):
    for row in rows:
        valid, invalid = integer(row["bootstrap_valid"]), integer(row["bootstrap_invalid"])
        require(integer(row["bootstrap_requested"]) == draws and 0 <= valid <= draws and valid + invalid == draws, "Bootstrap count inconsistency")
        close(row["bootstrap_valid_fraction"], valid / draws, "Bootstrap valid fraction mismatch")
        point, lower, upper = [number(row[key]) for key in ("point", "conditional_95_lower", "conditional_95_upper")]
        status = row["ci_status"]
        if status == "reported":
            require(all(map(math.isfinite, (point, lower, upper))) and lower <= upper and valid / draws >= minimum_fraction, "Invalid reported conditional interval")
        else:
            require(status in ("withheld_undefined_point", "withheld_insufficient_valid_draws") and math.isnan(lower) and math.isnan(upper), "Invalid withheld interval")
            require((not math.isfinite(point)) if status == "withheld_undefined_point" else valid / draws < minimum_fraction, "Withheld interval reason inconsistent")
        require(row["family_size"] == "" and row["family_adjusted_conditional_lower"] == "" and row["family_adjusted_conditional_upper"] == "", "Exploratory diagnostic unexpectedly contains family-adjusted intervals")
        require(row["interval_scope"] == "conditional_on_observed_OOF_records_not_full_train_refit_or_guaranteed_algorithmic_CV_coverage", "Interval scope changed")


def verify(results_directory, root=ROOT):
    results_directory, root = Path(results_directory).resolve(), Path(root).resolve()
    # Fail explicitly when training/statistics have not yet produced the final
    # public snapshot. Missing outputs must never count as a passing CI check.
    provenance_path = results_directory / "analysis_provenance.json"
    require(provenance_path.is_file(), "Completed public ordinal provenance is missing")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    require(provenance.get("status") == "completed", "Public analysis is incomplete")
    config_path, lock_path = root / "config/ordinal_diagnostic_v1.json", root / "config/ordinal_diagnostic_v1.lock.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    for field, path in (("config_sha256", config_path), ("lock_sha256", lock_path), ("analyzer_sha256", root / "src/analyze_ordinal_diagnostic.py")):
        require(provenance.get(field) == sha(path), "Public frozen source binding mismatch: " + field)
    require(lock["files"]["src/analyze_ordinal_diagnostic.py"] == provenance["analyzer_sha256"], "Analyzer does not match scientific freeze")
    require(provenance.get("source_run") == "runs/" + config["run_name"], "Unplanned source run")
    require(provenance.get("external_response_accessed") is False and provenance.get("blocked_source_access_attempts") == [], "External response access guard is not clear")
    require(provenance.get("no_permutation_p_values") is True and provenance.get("new_model_fits") == 0, "Statistics-only/no-P-value boundary changed")
    for flag in ("no_repeat_N_inflation", "all_null_runs_retained", "bootstrap_shared_across_models_drugs"):
        require(provenance.get(flag) is True, "Missing analysis boundary: " + flag)
    require(provenance["confidence_scope"] == config["statistics"]["confidence_scope"], "Confidence scope differs from frozen config")
    require(set(provenance.get("output_sha256", {})) == FILES, "Public output manifest is incomplete or contains an unplanned artifact")
    frames = {}
    for filename in sorted(FILES):
        path = results_directory / filename
        require(path.is_file(), "Missing public aggregate: " + filename)
        require(sha(path) == provenance["output_sha256"][filename], "Public file hash mismatch: " + filename)
        frames[filename] = table(path)
    drugs = list(map(str, config["drug_ids"]))
    scopes, models = [*drugs, MACRO], config["models"]
    permutations = list(map(str, range(config["permutation"]["count"] + 1)))
    nulls, folds = permutations[1:], list(map(str, range(config["outer_folds"])))
    molecular = [model for model in models if model != "tissue_null"]
    contract = lock["contract"]
    for field, expected in (("expected_outer_tasks", contract["tasks"]), ("actual_fit_events", contract["fits"]),
                            ("oof_rows", contract["oof_predictions"]), ("n_global_biological_groups", contract["global_development_groups"]),
                            ("bootstrap_draws", config["statistics"]["bootstrap_draws"]), ("bootstrap_seed", config["statistics"]["bootstrap_seed"])):
        require(provenance.get(field) == expected, "Provenance count/seed differs from frozen contract: " + field)
    point = index_rows(frames["point_metrics.csv"], ("gdsc_drug_id", "permutation_id", "model"), set(itertools.product(scopes, permutations, models)), "point metrics")
    for key, row in point.items():
        drug, perm, model = key
        require(integer(row["n_bootstrap_global_groups"]) == contract["global_development_groups"], "Global group count changed")
        n = integer(row["n_unique_biological_groups"])
        require(n == 620 if drug != MACRO else 620 <= n <= contract["global_development_groups"], "Case independent-group count invalid")
        for metric in METRICS:
            value = number(row[metric])
            require(not math.isinf(value), "Infinite point metric")
            if math.isfinite(value):
                require((-1.000000001 <= value <= 1.000000001) if metric == "spearman" else value >= 0 if metric in ("rmse", "mae") else value <= 1.000000001, "Metric outside its mathematical range")
        require(number(row["prediction_rms_difference_from_tissue_null"]) >= 0, "Negative prediction RMS difference")
        if model == "tissue_null":
            close(row["prediction_rms_difference_from_tissue_null"], 0, "Tissue-null reference differs from itself")
        if drug == MACRO:
            for metric in (*METRICS, "prediction_rms_difference_from_tissue_null"):
                close(row[metric], sum(number(point[(d, perm, model)][metric]) for d in drugs) / len(drugs), "Four-case macro does not equal fixed-case arithmetic mean")
    model_ci = index_rows(frames["model_intervals.csv"], ("gdsc_drug_id", "model", "metric"), set(itertools.product(scopes, models, METRICS)), "model intervals")
    for (drug, model, metric), row in model_ci.items():
        close(row["point"], point[(drug, "0", model)][metric], "Model CI point differs from observed OOF aggregate")
    planned = contrasts(config)
    comparison_names = [item[0] for item in planned]
    effects = index_rows(frames["paired_effect_intervals.csv"], ("gdsc_drug_id", "contrast", "metric"), set(itertools.product(scopes, comparison_names, METRICS)), "paired intervals")
    null_effects = index_rows(frames["null_contrasts.csv"], ("gdsc_drug_id", "permutation_id", "contrast", "metric"), set(itertools.product(scopes, nulls, comparison_names, METRICS)), "all null effects")
    for scope in scopes:
        for name, candidate, reference in planned:
            for metric in METRICS:
                sign = -1 if metric in ("rmse", "mae") else 1
                for perm in permutations:
                    row = effects[(scope, name, metric)] if perm == "0" else null_effects[(scope, perm, name, metric)]
                    require(row["candidate"] == candidate and row["reference"] == reference and row["positive_favors_candidate"] == "True", "Paired contrast identity/direction changed")
                    close(row["point"], sign * (number(point[(scope, perm, candidate)][metric]) - number(point[(scope, perm, reference)][metric])), "Paired effect disagrees with published point metrics")
    check_intervals(list(model_ci.values()) + list(effects.values()), config["statistics"]["bootstrap_draws"], config["statistics"]["minimum_valid_fraction"])
    summaries = index_rows(frames["null_contrast_summary.csv"], ("gdsc_drug_id", "contrast", "metric"), set(itertools.product(scopes, comparison_names, METRICS)), "null summaries")
    for (scope, contrast, metric), row in summaries.items():
        values = [number(null_effects[(scope, perm, contrast, metric)]["point"]) for perm in nulls]
        valid = sum(math.isfinite(value) for value in values)
        require(integer(row["null_runs_requested"]) == len(nulls) == integer(row["null_runs_recorded"]) and integer(row["null_runs_valid"]) == valid, "Null run omission/count mismatch")
        for name, fn in (("null_min", min), ("null_median", statistics.median), ("null_max", max)):
            close(row[name], fn(values) if valid == len(values) else math.nan, "Null range/median changed or undefined runs dropped")
    grouped = set(itertools.product(drugs, permutations, models))
    task_keys = set(itertools.product(drugs, permutations, models, folds))
    audit = index_rows(frames["molecular_and_failure_summary.csv"], ("gdsc_drug_id", "permutation_id", "model"), grouped, "molecular audit")
    params = index_rows(frames["chosen_parameters.csv"], ("gdsc_drug_id", "permutation_id", "model", "outer_fold"), task_keys, "chosen parameters")
    fit_events = index_rows(frames["all_fit_events_summary.csv"], ("gdsc_drug_id", "permutation_id", "model"), grouped, "all fit events")
    resource_tasks = index_rows(frames["resource_tasks.csv"], ("gdsc_drug_id", "permutation_id", "model", "outer_fold"), task_keys, "task resources")
    resources = index_rows(frames["resource_summary.csv"], ("gdsc_drug_id", "permutation_id", "model"), grouped, "resource summaries")
    total_fits = 0
    for key in grouped:
        a, f, r = audit[key], fit_events[key], resources[key]
        model = key[2]
        require(integer(a["outer_fits"]) == integer(r["outer_tasks"]) == len(folds), "Outer-fit summary count changed")
        expected_fits = len(folds) * (1 if model == "tissue_null" else config["inner_folds"] * len(config["elastic_net"]["alpha_ratios"]) * len(config["elastic_net"]["l1_ratios"]) + 1)
        require(integer(f["actual_fit_events"]) == expected_fits, "Actual fit count differs from frozen per-model budget")
        total_fits += expected_fits
        chosen = [params[(*key, fold)] for fold in folds]
        if model != "tissue_null":
            active = [integer(row["active_features"]) for row in chosen]
            require(integer(a["zero_molecular_fits"]) == sum(n == 0 for n in active) and integer(a["nonzero_molecular_fits"]) == sum(n > 0 for n in active), "Zero/nonzero outer-fit summary mismatch")
            require(integer(a["active_feature_total"]) == sum(active) and integer(a["active_feature_min"]) == min(active) and integer(a["active_feature_max"]) == max(active), "Active-feature summary mismatch")
            for row in chosen:
                if row["selected_genes"] == "" or row["candidate_features"] == "":
                    require(bool(row["fallback_reason"]) and integer(row["active_features"]) == 0, "Missing feature counts without a retained zero-molecular fallback")
                    genes, features = 0, 0
                else:
                    genes, features = integer(row["selected_genes"]), integer(row["candidate_features"])
                require(0 <= genes <= config["feature"]["n_genes"] and 0 <= integer(row["active_features"]) <= features, "Gene/active-feature budget invalid")
                require(features == (genes * (genes - 1) // 2 if model.endswith("pair_full") else genes), "Representation feature count does not match selected genes")
                selected_params = json.loads(row["chosen_params"])
                require(selected_params["alpha_ratio"] in config["elastic_net"]["alpha_ratios"] and selected_params["l1_ratio"] in config["elastic_net"]["l1_ratios"], "Unplanned selected EN parameter")
        require(integer(a["warning_outer_fits"]) == sum(bool(json.loads(row["warnings"])) for row in chosen), "Outer-warning summary mismatch")
        require(integer(a["fallback_outer_fits"]) == sum(bool(row["fallback_reason"]) for row in chosen), "Outer-fallback summary mismatch")
        for name in ("all_fit_warning_events", "all_fit_fallback_events"):
            require(0 <= integer(f[name]) <= expected_fits, "Invalid fit warning/fallback count")
        tasks = [resource_tasks[(*key, fold)] for fold in folds]
        for source, target, fn in (("task_wall_seconds", "task_wall_seconds_sum", sum), ("task_cpu_seconds", "task_cpu_seconds_sum", sum),
            ("outer_fit_seconds", "outer_fit_seconds_sum", sum), ("outer_prediction_seconds", "outer_prediction_seconds_sum", sum),
            ("sampled_peak_rss_bytes", "sampled_peak_rss_bytes_max", max), ("windows_lifetime_peak_working_set_bytes", "windows_lifetime_peak_working_set_bytes_max", max)):
            values = [number(t[source]) for t in tasks]
            require(all(math.isnan(v) or (math.isfinite(v) and v >= 0) for v in values), "Invalid resource measurement")
            if all(math.isfinite(v) for v in values):
                close(r[target], fn(values), "Resource aggregate does not match task values")
    require(total_fits == provenance["actual_fit_events"], "Total fit count mismatch")
    pair_keys = {(d, p, model, a, b) for d, p, model in itertools.product(drugs, permutations, molecular) for a, b in itertools.combinations(folds, 2)}
    stability = index_rows(frames["active_feature_stability.csv"], ("gdsc_drug_id", "permutation_id", "model", "fold_a", "fold_b"), pair_keys, "descriptive feature stability")
    for row in stability.values():
        for field in ("candidate_gene_jaccard", "active_feature_jaccard", "same_sign_active_jaccard"):
            value = number(row[field])
            require(math.isnan(value) or 0 <= value <= 1, "Jaccard outside [0,1]")
        if row["empty_active_union"] == "True":
            require(math.isnan(number(row["active_feature_jaccard"])) and math.isnan(number(row["same_sign_active_jaccard"])), "Empty active sets must not masquerade as perfect stability")
    coverage = index_rows(frames["permutation_coverage.csv"], ("drug_id", "permutation_id"), set(itertools.product(drugs, permutations)), "permutation coverage")
    for (drug, perm), row in coverage.items():
        n, singleton, identities, unchanged = [integer(row[field]) for field in ("n_groups", "singleton_groups", "identity_donors", "unchanged_response_count")]
        require(n == 620 and 0 <= singleton <= identities <= unchanged <= n, "Invalid permutation identity/response coverage")
        close(row["changed_response_fraction"], (n - unchanged) / n, "Permutation change fraction mismatch")
        if perm == "0":
            require(identities == unchanged == n, "Observed labels unexpectedly permuted")
    require(sum(integer(coverage[(drug, "0")]["n_groups"]) for drug in drugs) * len(models) * len(permutations) == provenance["oof_rows"], "OOF row contract inconsistent with public coverage")
    return {"status": "passed", "files_verified": len(FILES), "outer_tasks": len(task_keys), "successful_pipeline_fits": total_fits,
        "observed_paired_intervals": len(effects), "null_contrast_points": len(null_effects),
        "private_data_opened": False, "new_statistics_or_fits": False,
        "scope": "Public hashes and aggregate identities/counts only; private OOF/rank/bootstrap/checkpoint and scientific validity require separate full validation. Mutable reports, figures and operational-overhead summaries are excluded."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=ROOT / "results/ordinal_diagnostic")
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.results_dir), ensure_ascii=False))
        return 0
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        print(json.dumps({"status": "failed", "error": str(exc), "private_data_opened": False}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    sys.exit(main())
