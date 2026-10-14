"""Independent split and output-contract validation for tissue LOTO v1.3."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import run_tissue_loto_v1_3 as loto


ROOT = loto.ROOT


def validate_preflight(run: Path, config: dict, model_config: dict, lock_sha: str) -> dict:
    meta = loto.json_load(run / "loto_metadata.json")
    if meta["status"] != "preflight_completed" or meta["freeze_sha256"] != lock_sha:
        raise ValueError("No valid frozen LOTO preflight exists")
    for name, expected in meta["prepared_input_hashes"].items():
        if loto.sha(run / name) != expected:
            raise ValueError("Prepared input hash mismatch: " + name)

    response = pd.read_csv(
        run / "development_responses.csv",
        dtype={"gdsc_drug_id": str, "canonical_model_id": str, "biological_group_id": str, "tissue": str},
        keep_default_na=False, float_precision="round_trip",
    )
    plan = loto.read_jsonl(run / "task_plan.jsonl")
    tasks = loto.task_list(config, model_config)
    if len(plan) != len(tasks) or len(tasks) != 1472:
        raise ValueError("LOTO task plan has the wrong size")
    assignment_cache = {}
    warning_counts = Counter()
    warning_contexts = 0
    for i, (planned, task) in enumerate(zip(plan, tasks)):
        if planned["task"] != task:
            raise ValueError(f"Task ordering mismatch at index {i}")
        data = response.loc[response.gdsc_drug_id.eq(str(task["drug_id"]))].sort_values("canonical_model_id").reset_index(drop=True)
        split_task = {key: task[key] for key in ("drug_id", "heldout_tissue", "tissue_index", "inner_seed")}
        train, test, inner, details = loto.compute_split(data, split_task, model_config)
        assignment = details.pop("validation_assignment")
        if task["model"] == "tissue_null":
            warning_counts.update(details["splitter_warnings"])
            warning_contexts += int(bool(details["splitter_warnings"]))
        if details != planned["split"]:
            raise ValueError(f"Independent split reconstruction mismatch at task {i}")
        assignment_cache[(int(task["drug_id"]), task["heldout_tissue"])] = assignment
        if set(data.iloc[test].tissue) != {task["heldout_tissue"]}:
            raise ValueError("Outer test rows include a tissue other than the declared holdout")
        if task["heldout_tissue"] in set(data.iloc[train].tissue):
            raise ValueError("Held-out tissue appears in outer training rows")
        if len(inner) != 3:
            raise ValueError("An outer training set does not have three inner folds")

    inner = pd.read_csv(run / "loto_inner_fold_manifest.csv", dtype={"canonical_model_id": str, "biological_group_id": str, "tissue": str})
    if inner.duplicated(["drug_id", "heldout_tissue", "canonical_model_id"]).any():
        raise ValueError("A training sample is duplicated in the inner split manifest")
    expected_rows = 0
    for drug in config["drug_ids"]:
        data = response.loc[response.gdsc_drug_id.eq(str(drug))]
        for tissue in config["tissue_order"]:
            expected_rows += int(data.tissue.ne(tissue).sum())
            subset = inner.loc[inner.drug_id.eq(drug) & inner.heldout_tissue.eq(tissue)]
            if len(subset) != int(data.tissue.ne(tissue).sum()):
                raise ValueError(f"Inner validation manifest coverage mismatch for {drug}/{tissue}")
            counts = subset.validation_fold.value_counts().to_dict()
            if sorted(counts) != [0, 1, 2]:
                raise ValueError(f"Inner split is missing a validation fold for {drug}/{tissue}")
            if set(subset.tissue) != set(data.loc[data.tissue.ne(tissue), "tissue"]):
                raise ValueError("Inner manifest contains unexpected tissues")
            observed_assignment = dict(zip(subset.canonical_model_id.astype(str), subset.validation_fold.astype(int)))
            if observed_assignment != assignment_cache[(int(drug), tissue)]:
                raise ValueError(f"Inner-fold membership mismatch for {drug}/{tissue}")
    if len(inner) != expected_rows:
        raise ValueError("Inner split manifest row count mismatch")
    preflight_manifest = loto.json_load(run / "preflight_manifest.json")
    declared_split_summary = preflight_manifest["split_summary"]
    declared_warnings = declared_split_summary.get("inner_splitter_warning_message_counts", {})
    if declared_warnings != dict(warning_counts):
        raise ValueError("Recorded inner splitter warning counts differ from independent reconstruction")
    summary = {
        "status": "preflight_validated_no_model_fits",
        "tasks": len(tasks),
        "outer_drug_tissue_splits": 184,
        "inner_manifest_rows": len(inner),
        "all_23_tissues_present": True,
        "training_excludes_heldout_tissue": True,
        "grouped_inner_fold_coverage_complete": True,
        "inner_splitter_warning_contexts": warning_contexts,
        "inner_splitter_warning_message_counts": dict(warning_counts),
        "expected_fits": config["expected_contract"]["fits"],
        "expected_oof_predictions": config["expected_contract"]["oof_predictions"],
    }
    return summary


def validate_run(run: Path, config: dict | None = None, model_config: dict | None = None,
                 lock: dict | None = None, lock_sha: str | None = None) -> dict:
    if config is None:
        lock, lock_sha, _, _ = loto.load_execution_lock()
        config = loto.json_load(loto.EXEC_CONFIG_PATH)
        model_config = loto.json_load(loto.MODEL_CONFIG_PATH)
    meta = loto.json_load(run / "loto_metadata.json")
    if meta["status"] != "completed" or not meta["completed"]:
        raise ValueError("LOTO model run is not complete")
    if meta["freeze_sha256"] != lock_sha:
        raise ValueError("LOTO run was produced under a different lock")
    expected = config["expected_contract"]
    observed = meta["output_counts"]
    if (observed["tasks"], observed["fit_events"], observed["oof_predictions"]) != (
        expected["tasks"], expected["fits"], expected["oof_predictions"]
    ):
        raise ValueError("LOTO run metadata does not match the fixed workload")
    if meta.get("controller_data_guard", {}).get("external_response_accessed") is not False or meta.get("controller_data_guard", {}).get("blocked_attempts"):
        raise ValueError("Controller data guard reported external response access")

    for name, expected in meta["prepared_input_hashes"].items():
        if loto.sha(run / name) != expected:
            raise ValueError("Prepared input hash mismatch: " + name)
    for name, expected in meta["output_hashes"].items():
        if loto.sha(run / name) != expected:
            raise ValueError("Aggregated output hash mismatch: " + name)

    response = pd.read_csv(
        run / "development_responses.csv",
        dtype={"gdsc_drug_id": str, "canonical_model_id": str, "biological_group_id": str, "tissue": str},
        keep_default_na=False, float_precision="round_trip",
    )
    tasks = loto.task_list(config, model_config)
    checkpoints = {}
    zero_solutions = Counter()
    warning_count = 0
    fallback_count = 0
    checkpoint_hashes = meta["checkpoint_hashes"]
    for task in tasks:
        path = run / "checkpoints" / loto.checkpoint_name(task)
        if not loto.validate_checkpoint(path, task, meta["freeze_sha256"]):
            raise ValueError("Missing or invalid checkpoint: " + path.name)
        if checkpoint_hashes.get(path.name) != loto.sha(path):
            raise ValueError("Checkpoint hash mismatch: " + path.name)
        obj = loto.json_load(path)
        expected_data = response.loc[response.gdsc_drug_id.eq(str(task["drug_id"]))].sort_values("canonical_model_id").reset_index(drop=True)
        train, test, inner, split = loto.compute_split(
            expected_data,
            {key: task[key] for key in ("drug_id", "heldout_tissue", "tissue_index", "inner_seed")},
            model_config,
        )
        expected_test = sorted(expected_data.iloc[test].canonical_model_id.astype(str).tolist())
        expected_train_expr = sorted(expected_data.expression_row_index.astype(int).to_numpy()[train].tolist())
        expected_inner_expr = [
            sorted(expected_data.expression_row_index.astype(int).to_numpy()[train[np.asarray(inner_train, dtype=int)]].tolist())
            for inner_train, inner_valid in inner
        ]
        if obj.get("expected_test_ids") != expected_test:
            raise ValueError("Checkpoint test identity mismatch: " + path.name)
        if obj.get("expected_training_expression_ids") != expected_train_expr:
            raise ValueError("Checkpoint training expression membership mismatch: " + path.name)
        if obj.get("inner_training_expression_ids") != expected_inner_expr:
            raise ValueError("Checkpoint inner training membership mismatch: " + path.name)
        if obj.get("inner_validation_assignment_sha256") != split["inner_validation_assignment_sha256"]:
            raise ValueError("Checkpoint inner validation split mismatch: " + path.name)
        assigned_y = expected_data.response_auc.to_numpy(dtype=np.float64)
        if obj.get("assigned_response_sha256") != hashlib.sha256(assigned_y.tobytes()).hexdigest():
            raise ValueError("Checkpoint response vector differs from its frozen source: " + path.name)
        if obj.get("data_guard", {}).get("external_response_accessed") is not False or obj.get("data_guard", {}).get("blocked_attempts"):
            raise ValueError("External-data guard reported access/block attempts: " + path.name)
        if obj.get("status") != "completed" or obj.get("error"):
            raise ValueError("Failed task is present in completed LOTO run")
        if len(obj["oof_rows"]) != len(expected_test):
            raise ValueError("Checkpoint OOF row count mismatch: " + path.name)
        if int(obj["result"]["outer_unknown_tissue_samples"]) != len(expected_test):
            raise ValueError("Unknown-tissue sample count differs from held-out test size")
        if task["model"] == "tissue_null" and obj["result"]["metrics"]["spearman"] is not None:
            raise ValueError("Constant held-out tissue-null Spearman must be NA")
        expected_fit_count = 1 if task["model"] == "tissue_null" else 31
        if len(obj["fit_events"]) != expected_fit_count:
            raise ValueError("Per-task fit count mismatch: " + path.name)
        if task["model"] != "tissue_null":
            inner_events = [event for event in obj["fit_events"] if event.get("stage") == "inner_training"]
            outer_events = [event for event in obj["fit_events"] if event.get("stage") == "outer_training_refit"]
            inner_coordinates = Counter((int(event["candidate"]), int(event["inner_fold"])) for event in inner_events)
            expected_coordinates = Counter((candidate, fold) for candidate in range(10) for fold in range(3))
            if inner_coordinates != expected_coordinates or len(outer_events) != 1:
                raise ValueError("Inner tuning grid or outer refit event coverage mismatch: " + path.name)
        test_expr_ids = set(expected_data.expression_row_index.astype(int).to_numpy()[test].tolist())
        for event in obj["fit_events"]:
            fit_ids = set(map(int, event.get("fit_sample_ids", [])))
            if fit_ids & test_expr_ids:
                raise ValueError("An inner/outer fit includes held-out sample IDs")
            warning_count += int(bool(event.get("fit_warning") or event.get("warning")))
            fallback_count += int(bool(event.get("fallback_reason")))
        outer = obj["result"].get("outer_fit_audit", {})
        if task["model"] != "tissue_null":
            molecular = outer.get("n_nonzero_molecular_coeffs")
            if molecular is not None and int(molecular) == 0:
                zero_solutions[task["model"]] += 1
        key = (task["drug_id"], task["heldout_tissue"], task["model"])
        if key in checkpoints:
            raise ValueError("Duplicate LOTO checkpoint mapping")
        checkpoints[key] = obj

    oof = pd.read_csv(
        run / "oof_predictions.csv",
        dtype={"canonical_model_id": str, "biological_group_id": str, "tissue": str,
               "heldout_tissue": str, "drug_name": str, "model": str},
        keep_default_na=False, float_precision="round_trip",
    )
    if len(oof) != config["expected_contract"]["oof_predictions"]:
        raise ValueError("LOTO OOF table row count mismatch")
    if "inner_seed" not in oof.columns:
        raise ValueError("LOTO OOF table is missing the locked inner_seed task identity")
    expected_seed = oof.heldout_tissue.map(config["inner_seed_resolution"]["seed_by_tissue"])
    if expected_seed.isna().any() or not np.array_equal(
        oof.inner_seed.to_numpy(dtype=np.int64), expected_seed.to_numpy(dtype=np.int64)
    ):
        raise ValueError("LOTO OOF inner_seed differs from its locked held-out tissue seed")
    if oof.duplicated(["drug_id", "model", "canonical_model_id"]).any():
        raise ValueError("A sample has duplicate predictions for a drug/model")
    if (oof.tissue != oof.heldout_tissue).any():
        raise ValueError("An OOF row is not assigned to the tissue held out for prediction")
    if not np.isfinite(oof[["observed_auc", "predicted_auc"]].to_numpy(dtype=float)).all():
        raise ValueError("Nonfinite LOTO OOF value")

    expected_keys = set()
    for task in tasks:
        data = response.loc[response.gdsc_drug_id.eq(str(task["drug_id"]))]
        expected_keys.update((int(task["drug_id"]), task["model"], str(model_id))
                             for model_id in data.loc[data.tissue.eq(task["heldout_tissue"]), "canonical_model_id"])
    observed_keys = set(zip(oof.drug_id.astype(int), oof.model.astype(str), oof.canonical_model_id.astype(str)))
    if observed_keys != expected_keys:
        raise ValueError("LOTO OOF identity set does not match all frozen drug/model/tissue cases")
    oof_for_merge = oof.assign(_drug_key=oof.drug_id.astype(int).astype(str))
    source_for_merge = response[["gdsc_drug_id", "canonical_model_id", "biological_group_id", "tissue", "response_auc"]]
    merged = oof_for_merge.merge(
        source_for_merge, left_on=["_drug_key", "canonical_model_id"],
        right_on=["gdsc_drug_id", "canonical_model_id"], how="left", validate="many_to_one",
    )
    if merged.response_auc.isna().any():
        raise ValueError("An OOF sample has no corresponding frozen source response")
    if not np.array_equal(merged.observed_auc.to_numpy(dtype=np.float64), merged.response_auc.to_numpy(dtype=np.float64)):
        raise ValueError("Observed OOF labels differ from the frozen source responses")
    for field in ("biological_group_id", "tissue"):
        if not (merged[field + "_x"].astype(str) == merged[field + "_y"].astype(str)).all():
            raise ValueError(f"Observed OOF {field} differs from the frozen source identity")

    shared = oof.groupby(["drug_id", "heldout_tissue", "canonical_model_id"], sort=False).agg(
        model_count=("model", "nunique"), y_count=("observed_auc", "nunique"),
        tissue_count=("tissue", "nunique"), group_count=("biological_group_id", "nunique"),
    )
    if len(shared) != int(oof[["drug_id", "heldout_tissue", "canonical_model_id"]].drop_duplicates().shape[0]):
        raise ValueError("LOTO shared sample grouping failed")
    if not shared.eq({"model_count": 8, "y_count": 1, "tissue_count": 1, "group_count": 1}).all().all():
        raise ValueError("LOTO model arms do not share identities, labels, tissues and groups")

    fit_events = loto.read_jsonl(run / "fit_events.jsonl")
    if len(fit_events) != config["expected_contract"]["fits"]:
        raise ValueError("Aggregated fit event count mismatch")
    audit_rows = loto.read_jsonl(run / "outer_fit_audits.jsonl")
    if len(audit_rows) != config["expected_contract"]["tasks"]:
        raise ValueError("Outer-fit audit row count mismatch")
    audit_keys = {
        (int(row["task"]["drug_id"]), row["task"]["heldout_tissue"], row["task"]["model"])
        for row in audit_rows
    }
    expected_audit_keys = {(int(task["drug_id"]), task["heldout_tissue"], task["model"]) for task in tasks}
    if audit_keys != expected_audit_keys:
        raise ValueError("Outer-fit audit identities differ from the locked task list")
    if len(checkpoint_hashes) != config["expected_contract"]["tasks"]:
        raise ValueError("Checkpoint registry is incomplete")
    fit_count_by_task = Counter((int(row["drug_id"]), row["heldout_tissue"], row["model"]) for row in fit_events)
    expected_fit_counts = Counter({
        (int(task["drug_id"]), task["heldout_tissue"], task["model"]):
        (1 if task["model"] == "tissue_null" else 31) for task in tasks
    })
    if fit_count_by_task != expected_fit_counts:
        raise ValueError("Fit events are not complete for every locked outer task")

    summary = {
        "status": "completed_contract_validated",
        "tasks": len(checkpoints),
        "fits": len(fit_events),
        "oof_predictions": len(oof),
        "shared_drug_tissue_sample_sets": len(shared),
        "unknown_tissue_samples_equal_test_samples": True,
        "training_includes_heldout_tissue": False,
        "checkpoint_hashes_verified": len(checkpoint_hashes),
        "zero_molecular_outer_fits_by_model": dict(zero_solutions),
        "fit_warning_events": warning_count,
        "fit_fallback_events": fallback_count,
        "external_response_accessed": False,
        "blocked_external_access_attempts": 0,
    }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--preflight", action="store_true", help="Validate no-fit split preflight only")
    args = parser.parse_args()
    lock, lock_sha, _, _ = loto.load_execution_lock()
    config = loto.json_load(loto.EXEC_CONFIG_PATH)
    model_config = loto.json_load(loto.MODEL_CONFIG_PATH)
    run = (args.run_dir or ROOT / "runs" / config["run_name"]).resolve()
    run.relative_to(ROOT / "runs")
    if args.preflight:
        result = validate_preflight(run, config, model_config, lock_sha)
    else:
        result = validate_run(run, config, model_config, lock, lock_sha)
        loto.atomic_json(run / "validation_manifest.json", result)
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
