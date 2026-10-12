"""Independently audit pilot artifacts; never open source response/expression data.

Usage: .venv/Scripts/python.exe scripts/validate_pilot_run.py --run-dir runs/NAME
This validator uses run artifacts, the frozen development identity manifest,
cohort lock metadata, and Git source blobs. It does not import the model runner.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
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


def names_sha(values):
    return hashlib.sha256("\0".join(values).encode()).hexdigest()


def source_provenance(meta, revision):
    revision = meta.get("code_execution_commit") or revision or "629b8f5"
    full_revision = subprocess.check_output(["git", "rev-parse", revision], cwd=ROOT, text=True).strip()
    source_checks = []
    for name, recorded_sha in meta["source_hashes"].items():
        git_path = name.replace("\\", "/")
        blob = subprocess.check_output(["git", "show", f"{revision}:{git_path}"], cwd=ROOT)
        require(hashlib.sha256(blob).hexdigest() == recorded_sha, "Recorded run source does not match immutable Git revision: " + git_path)
        source_checks.append({"path": git_path, "recorded_sha256": recorded_sha, "revision_matches": True,
                              "current_working_copy_matches": sha(ROOT / git_path) == recorded_sha})
    require(meta["external_response_accessed"] is False and not meta["blocked_external_access_attempts"], "Run reports external access or blocked attempts")
    require(not any("ctrp" in path.lower() for path in meta["input_files_accessed"]), "Access log includes an external-response path")
    return full_revision, source_checks


def validate_relative_path_event(event, config):
    if event["model"] not in {"expression", "expression_matched", "reo"}:
        return
    if config["elastic_net"].get("regularization_mode", "absolute") != "alpha_max_relative":
        return
    if event["alpha_max"] is None:
        require(event["fallback_reason"].startswith(("empty_features:", "constant_training_response")), "Missing alpha_max without documented fallback")
        return
    maximum, actual, ratio = event["alpha_max"], event["actual_alpha"], event["params"]["alpha_ratio"]
    require(np.isfinite(maximum) and maximum >= 0 and 0 < ratio <= 1, "Invalid training alpha_max or path ratio")
    require(np.isclose(actual, maximum * ratio, atol=1e-14, rtol=1e-12) and event["path_ratio"] == ratio, "Actual alpha differs from current-training relative path")
    require(isinstance(event["n_nonzero_molecular_coeffs"], int) and 0 <= event["n_nonzero_molecular_coeffs"] <= event["n_features"], "Invalid molecular coefficient count")
    if maximum == 0:
        require(event["n_nonzero_molecular_coeffs"] == 0 and event["n_iter"] == 0, "Zero alpha_max did not retain explicit zero molecular solution")
    if ratio == 1.0:
        require(event["n_nonzero_molecular_coeffs"] == 0, "Strongest relative path endpoint is not a zero molecular solution")
        if maximum > 0 and "null_endpoint" in event:
            require(event["null_endpoint"] is True and event["solver_reason"] == "theoretical_alpha_max_null_endpoint" and not event["fallback_reason"], "Analytic null endpoint incorrectly treated as fallback or non-null fit")


def validate_outer_feature_audits(run, config, outer_events):
    for (drug, model, outer), event in outer_events.items():
        if event["n_selected_genes"] is None:
            continue
        payload = read_json(run / "features" / f"{drug}_{model}_outer{outer}.json")
        genes, indices, names = payload["gene_ids"], payload["gene_indices"], payload["feature_names"]
        require(len(genes) == len(indices) == event["n_selected_genes"] and len(set(genes)) == len(genes) and len(set(indices)) == len(indices), "Outer feature gene identities are duplicated/misaligned")
        require(payload["training_samples"] == len(event["fit_sample_ids"]), "Outer feature audit training count differs from actual fit")
        require(names_sha(genes) == event["selected_gene_ids_sha256"] and names_sha(names) == event["built_feature_names_sha256"], "Outer feature lists do not match actual-fit schema hashes")
        identity_for_index = dict(zip(indices, genes))
        pairs = payload["pair_indices_original_input_columns"]
        if payload["kind"] == "reo":
            require(genes == sorted(genes) and len(pairs) == len(names), "REO gene orientation or pair alignment is inconsistent")
            require(len({tuple(pair) for pair in pairs}) == len(pairs), "Duplicated candidate REO pairs")
            for pair, name in zip(pairs, names):
                require(len(pair) == 2 and pair[0] in identity_for_index and pair[1] in identity_for_index, "REO pair references a gene outside selected identities")
                first, second = identity_for_index[pair[0]], identity_for_index[pair[1]]
                require(first < second and name == f"REO:{first}>{second}", "REO stable name cannot be traced to canonical pair indices")
        else:
            require(names == ["log2TPM:" + gene for gene in genes] and pairs == [], "Expression schema differs from stable gene identities")
        if model.endswith("_rf"):
            forest = payload["random_forest"]
            require(forest["n_estimators"] == config["random_forest"]["n_estimators"] and forest["n_jobs"] == 1, "Forest tree/parallel budget differs from fixed config")
            require(forest["max_features"] == event["params"]["max_features"] and forest["min_samples_leaf"] == event["params"]["min_samples_leaf"], "Forest audit hyperparameters differ from actual chosen fit")
            require(forest["molecular_features_used"] == event["n_features"] == len(forest["molecular_feature_names"]), "Forest molecular feature count mismatch")
            if model == "reo_rf":
                selection = forest["pair_selection"]
                selected = selection["selected_feature_indices"]
                selected_names = selection["selected_feature_names"]
                require(len(selected) == len(set(selected)) <= config["random_forest"]["reo_max_pairs"] and all(0 <= index < len(names) for index in selected), "RF pair-selection indices exceed candidate schema or budget")
                require(selected_names == [names[index] for index in selected] == forest["molecular_feature_names"], "RF selected names do not match frozen candidate columns")
                require(selection["selected_pairs_original_input_gene_indices"] == [pairs[index] for index in selected], "RF selected pairs do not trace to original gene indices")
                require(names_sha(selected_names) == event["rf_selected_pair_names_sha256"], "RF selected pair schema differs from actual-fit hash")
                associations = selection["selected_signed_associations"]
                require(len(associations) == len(selected) and np.isfinite(associations).all() and all(abs(value) <= 1 for value in associations), "RF association audit invalid")
                require(selection["training_samples"] == len(event["fit_sample_ids"]) and selection["association"] == config["random_forest"]["association"], "RF supervised-selector training scope/association differs from protocol")
            else:
                require(forest["molecular_feature_names"] == names, "Expression RF changed molecular columns after training gene selection")
        else:
            require(event["n_features"] == len(names), "EN feature count differs from stable schema")


def validate_timing_run(run, meta, selected, membership, global_identity, revision):
    events = [json.loads(line) for line in (run / "fit_events.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    config = meta["configuration"]
    drugs = {str(d) for d in config["drug_ids"]}
    kinds = set(config["models"]) & {"expression", "reo", "expression_rf", "reo_rf"}
    require(meta["n_oof_predictions"] == 0 and not (run / "oof_predictions.csv").exists() and not (run / "chosen_parameters.json").exists(), "Timing probe contains prediction/selection outputs")
    observed = pd.read_csv(run / "development_responses.csv", dtype={"gdsc_drug_id": str})
    response = observed.set_index(["gdsc_drug_id", "expression_row_index"]).response_auc
    contexts = Counter()
    within_rms = []
    for event in events:
        drug, model = str(event["gdsc_drug_id"]), event["model"]
        require(drug in drugs and model in kinds and event["stage"] == "training_only_timing_probe" and event["outer_fold"] == 0 and event["inner_fold"] == 0 and event["candidate"] == -1, "Unexpected timing-probe context")
        identities = selected[selected.gdsc_drug_id.eq(drug)]
        lookup = dict(zip(identities.expression_row_index.astype(int), identities.canonical_model_id))
        ids = event["fit_sample_ids"]
        require(len(ids) == len(set(ids)) and set(ids) <= set(lookup), "Timing probe has invalid sample IDs")
        actual = {lookup[i] for i in ids}
        require(actual == membership[(0, 0, "train")] & set(identities.canonical_model_id), "Timing fit crosses its current training partition")
        require(event["prediction_calls"] == 0 and event["predict_seconds"] == 0, "Timing probe performed held-out predictions")
        y = np.array([response.loc[(drug, sample)] for sample in ids], dtype=float)
        require(np.isclose(event["target_mean"], y.mean(), atol=1e-12, rtol=0) and np.isclose(event["target_sd"], y.std() if y.std() else 1., atol=1e-12, rtol=0), "Timing target scaling not training-only")
        validate_relative_path_event(event, config)
        categories = np.array([global_identity.loc[lookup[sample], "tissue"] for sample in ids])
        residual_squares = sum(np.sum((y[categories == tissue] - y[categories == tissue].mean()) ** 2)
                               for tissue in np.unique(categories) if np.ptp(y[categories == tissue]) > 0)
        within_rms.append(float(np.sqrt(residual_squares / len(y))))
        if model == "reo_rf":
            require(0 <= event["n_features"] <= config["random_forest"]["reo_max_pairs"], "RF timing pair budget exceeded")
        contexts.update([(drug, model)])
    require(set(contexts) == set(itertools.product(drugs, kinds)) and set(contexts.values()) == {1} and len(events) == meta["n_fits"], "Timing fits incomplete/duplicated")
    full_revision, sources = source_provenance(meta, revision)
    return {"status": "passed", "validated_utc": datetime.now(timezone.utc).isoformat(), "run_mode": "training_only_timing",
            "run_directory": artifact_label(run), "source_revision": full_revision, "source_checks": sources,
            "counts": {"development_pool_identities": len(global_identity), "timing_fits": len(events), "oof_predictions": 0},
            "within_tissue_response_variation": {"minimum_positive_fit_residual_rms": min((value for value in within_rms if value > 0), default=None), "fits_with_exactly_zero_within_tissue_response": sum(value == 0 for value in within_rms)},
            "original_GDSC_or_CTRP_response_read_by_validator": False,
            "checks": {"training_membership": True, "no_tuning_or_heldout_predictions": True, "current_training_relative_path": True,
                       "recorded_sources_match_immutable_run_revision": True, "run_external_access_log_clear": True}}


def validate(run, revision):
    run = Path(run).resolve()
    meta = read_json(run / "run_manifest.json")
    config = meta["configuration"]
    drugs = {str(d) for d in config["drug_ids"]}
    models = set(config["models"])
    outer_count, inner_count = config["outer_folds"], config["inner_folds"]
    require(config["outer_repeats"] == 1, "This pilot validator requires one outer repeat")
    require(meta["status"] in {"completed", "completed_timing_only"}, "Pilot run is not completed")

    cohort_version = config.get("cohort_version", "v1")
    require(cohort_version in {"v1", "v2"}, "Unsupported cohort version")
    cohort_directory = ROOT / ("results/cohort_v2" if cohort_version == "v2" else "results")
    freeze_path = ROOT / f"config/cohort_freeze_{cohort_version}.json"
    require(sha(freeze_path) == meta["cohort_lock"]["freeze_sha256"], "Current cohort description differs from run lock")
    freeze = read_json(freeze_path)
    development_path = cohort_directory / "development_sample_manifest.csv"
    frozen_development = next(r for r in freeze["files"]
                              if str(r["path"]).replace("\\", "/") == development_path.relative_to(ROOT).as_posix())
    if cohort_version == "v2":
        require(sha(ROOT / "config/cohort_freeze_v1.json") == freeze["parent_cohort_freeze_sha256"] == meta["cohort_lock"]["parent_cohort_freeze_sha256"], "V2 run parent freeze binding changed")
        require(sha(ROOT / "config/cohort_freeze_v1.lock.json") == freeze["parent_lock_sha256"] == meta["cohort_lock"]["parent_lock_sha256"], "V2 run parent lock binding changed")
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
    if meta["status"] == "completed_timing_only":
        return validate_timing_run(run, meta, selected, membership, global_identity, revision)

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
        if model.endswith("_rf"):
            rf = config["random_forest"]
            return [{"max_features": f, "min_samples_leaf": leaf}
                    for f, leaf in itertools.product(rf["max_features"], rf["min_samples_leaf"])]
        if config["elastic_net"].get("regularization_mode", "absolute") == "alpha_max_relative":
            return [{"alpha_ratio": a, "l1_ratio": r}
                    for a, r in itertools.product(config["elastic_net"]["alpha_ratios"], config["elastic_net"]["l1_ratios"])]
        return [{"alpha": a, "l1_ratio": r} for a, r in itertools.product(config["elastic_net"]["alphas"], config["elastic_net"]["l1_ratios"])]
    event_keys = Counter()
    feature_fit_count = 0
    within_tissue_fit_rms = []
    within_tissue_positive_group_rms = []
    within_tissue_positive_ranges = []
    wholly_tissue_constant_fits = []
    context_hashes = defaultdict(lambda: defaultdict(set))
    context_thresholds = defaultdict(list)
    outer_events = {}
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
        categories = np.array([global_identity.loc[id_for_expression[sample], "tissue"] for sample in fit_samples])
        within_squares = 0.0
        for tissue in np.unique(categories):
            values = responses[categories == tissue]
            span = float(np.ptp(values))
            if span == 0:
                continue
            residual = values - values.mean()
            within_squares += float(np.sum(residual ** 2))
            within_tissue_positive_group_rms.append(float(np.sqrt(np.mean(residual ** 2))))
            within_tissue_positive_ranges.append(span)
        within_tissue_fit_rms.append(float(np.sqrt(within_squares / len(responses))))
        if within_squares == 0:
            wholly_tissue_constant_fits.append({"drug": drug, "model": model, "stage": stage, "outer_fold": outer, "inner_fold": inner})
            if cohort_version == "v2" and model == "reo_rf" and event["n_selected_genes"] is not None:
                require(event["n_features"] == 0, "Exactly tissue-constant training response caused spurious RF pair selection")
        if cohort_version == "v2":
            require(event["prediction_calls"] == 1 and event["predict_seconds"] >= 0 and event["fit_seconds"] >= 0, "Missing actual-fit prediction/resource audit")
            validate_relative_path_event(event, config)
            if model == "reo_rf":
                require(0 <= event["n_features"] <= config["random_forest"]["reo_max_pairs"], "RF pair budget exceeded")
            context = (drug, model, outer, stage, inner)
            for field in ("selected_gene_ids_sha256", "built_feature_names_sha256", "rf_selected_pair_names_sha256"):
                if event.get(field) is not None:
                    context_hashes[context][field].add(event[field])
            if model in {"expression", "expression_matched", "reo"} and event.get("alpha_max") is not None:
                context_thresholds[context].append(event["alpha_max"] * event["params"]["l1_ratio"])
        if model in {"expression", "expression_matched", "reo", "expression_rf", "reo_rf"} and event["n_selected_genes"] is not None:
            maximum_genes = (config["feature"]["expression_genes"] if model in {"expression", "expression_rf"}
                             else config["feature"].get("matched_expression_genes", config["feature"]["reo_genes"]) if model == "expression_matched"
                             else config["feature"]["reo_genes"])
            require(0 < event["n_selected_genes"] <= maximum_genes and event["n_features"] >= 0, "Feature builder was not recorded for this training fit")
            feature_fit_count += 1
        if stage == "outer_training_refit":
            outer_events[(drug, model, outer)] = event
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
    for hashes in context_hashes.values():
        require(all(len(values) == 1 for values in hashes.values()), "Training gene/pair schema changed across hyperparameters within the same partition")
    for values in context_thresholds.values():
        require(np.allclose(values, values[0], atol=1e-12, rtol=1e-10), "Current-training alpha_max*l1_ratio changed across path candidates")
    if cohort_version == "v2":
        validate_outer_feature_audits(run, config, outer_events)

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

    full_revision, source_checks = source_provenance(meta, revision)

    return {
        "status": "passed", "validated_utc": datetime.now(timezone.utc).isoformat(),
        "run_name": config["run_name"], "run_directory": artifact_label(run), "cohort_version": cohort_version,
        "validation_scope": "Recorded actual fit memberships, fold/group boundaries, OOF completeness, train-only target scaling, fixed-grid choice, OOF metric reconstruction and immutable source provenance",
        "inputs_read": ["run artifacts", "frozen development identity manifest", "cohort-lock metadata", "Git source blobs"],
        "original_GDSC_or_CTRP_response_read_by_validator": False,
        "source_revision": full_revision, "source_checks": source_checks,
        "counts": {"development_pool_identities": len(global_identity), "drug_model_oof_counts": counts,
                   "actual_fits": len(events), "expected_fits": expected_fit_count,
                   "feature_builder_fits_with_training_membership_verified": feature_fit_count,
                   "oof_predictions": len(predictions), "outer_folds": outer_count, "inner_folds": inner_count},
        "within_tissue_response_variation": {"fits_checked": len(events), "fits_with_all_tissues_exactly_constant": len(wholly_tissue_constant_fits),
                     "minimum_positive_fit_residual_rms": min((value for value in within_tissue_fit_rms if value > 0), default=None),
                     "minimum_positive_tissue_group_residual_rms": min(within_tissue_positive_group_rms, default=None),
                     "minimum_positive_within_tissue_response_range": min(within_tissue_positive_ranges, default=None),
                     "exact_constant_edge_activated_in_this_run": bool(wholly_tissue_constant_fits)},
        "checks": {"cohort_lock_and_development_identity_hash": True, "all_group_and_fold_boundaries": True,
                   "actual_inner_fit_membership_and_exactly_once_refits": True, "outer_test_not_used_for_fitting": True,
                   "each_identity_exactly_one_oof_prediction_per_drug_model": True,
                   "target_scaling_from_actual_training_responses": True, "fixed_grid_and_logged_inner_rmse_selection": True,
                   "oof_metrics_reconstructed": True, "mean_baseline_ranking_not_interpreted": True,
                   "recorded_sources_match_immutable_run_revision": True, "run_external_access_log_clear": True,
                   "relative_path_and_rf_schema_checks": True if cohort_version == "v2" else "not_applicable"},
        "limitations": ["Fit-boundary verification depends on actual-fit audit records; synthetic tests independently inspect builder inputs and perturb held-out labels/distributions.",
                        "Inner OOF predictions were not exported; parameter choice is checked against logged candidate scores, while pooled outer metrics are reconstructed from individual OOF records.",
                        "External access checks inspect the run's access log and guard tests; the validator does not reopen external responses."]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--source-revision", default=None, help="Fallback immutable execution revision; run_manifest.code_execution_commit takes priority")
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
