"""Execute the frozen 23-tissue leave-one-tissue-out sensitivity run.

The runner never reads CTRP inputs. It performs a no-fit split preflight first,
then one fresh worker process per locked tissue × drug × model task.
"""
from __future__ import annotations

from ordinal_repeated_runtime import set_preimport_threads
EARLY_THREADS = set_preimport_threads()

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import threading
import time
import traceback
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_info

import run_formal_evaluation as base
from ordinal_repeated_models import model_spec, evaluate_repeated_outer_fold
from ordinal_repeated_runtime import atomic_json, sha, system_memory, process_tree_memory, pressure_reasons
from run_development_pilot import execution_source_snapshot, install_data_guard, verify_development_inputs


ROOT = base.ROOT
EXEC_CONFIG_PATH = ROOT / "config/loto_execution_v1_2.json"
LOCK_PATH = ROOT / "config/loto_execution_v1_2.lock.json"
MODEL_CONFIG_PATH = ROOT / "config/ordinal_repeated_v1.json"
SOURCE_RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"
PREPARED_FILES = (
    "development_responses.csv",
    "global_fold_manifest.csv",
    "loto_inner_fold_manifest.csv",
    "task_plan.jsonl",
    "preflight_manifest.json",
)
SOURCE_FILES = (
    "src/run_tissue_loto_v1_2.py",
    "src/validate_tissue_loto_v1_2.py",
    "src/analyze_tissue_loto_v1_2.py",
    "src/run_formal_evaluation.py",
    "src/run_development_pilot.py",
    "src/project_paths.py",
    "src/ordinal_repeated_runtime.py",
    "src/ordinal_repeated_models.py",
    "src/ordinal_diagnostic_models.py",
    "src/modeling/nested_cv.py",
    "src/modeling/features.py",
    "src/modeling/reo_rf.py",
)
OOF_FIELDS = [
    "drug_id", "drug_name", "model", "repeat", "outer_fold", "canonical_model_id",
    "biological_group_id", "tissue", "observed_auc", "predicted_auc", "permutation_id",
    "selection_rule", "representation", "n_genes", "heldout_tissue", "tissue_index",
]


def json_load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_jsonl(path: Path, rows: list[dict]) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, allow_nan=False, sort_keys=True) + "\n")
    os.replace(temp, path)


def read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def hash_ids(values) -> str:
    payload = "\n".join(sorted(map(str, values))).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def hash_inner_assignment(mapping: dict[str, int]) -> str:
    payload = "\n".join(f"{key}\t{mapping[key]}" for key in sorted(mapping)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_execution_lock() -> tuple[dict, str, str, dict]:
    lock = json_load(LOCK_PATH)
    config = json_load(EXEC_CONFIG_PATH)
    model_config = json_load(MODEL_CONFIG_PATH)
    if sha(EXEC_CONFIG_PATH) != lock["configuration_sha256"]:
        raise ValueError("LOTO execution config differs from its lock")
    for name, expected in lock["locked_file_hashes"].items():
        if sha(ROOT / name) != expected:
            raise ValueError("Frozen LOTO input changed: " + name)
    if lock["execution_commit"] != lock["source_snapshot"]["git_commit"]:
        raise ValueError("Frozen LOTO execution commit is inconsistent")
    os.environ["REO_EXECUTION_COMMIT"] = lock["execution_commit"]
    revision, source_hashes = execution_source_snapshot([ROOT / name for name in SOURCE_FILES])
    if revision != lock["execution_commit"] or source_hashes != lock["source_snapshot"]["source_hashes"]:
        raise ValueError("Runner source snapshot differs from the frozen LOTO lock")
    if lock["p5e_lock_sha256"] != sha(ROOT / config["source_p5e_lock"]):
        raise ValueError("P5e protocol lock changed")
    if lock["model_config_sha256"] != sha(MODEL_CONFIG_PATH):
        raise ValueError("Fixed-eight model config changed")
    return lock, sha(LOCK_PATH), revision, source_hashes


def task_list(config: dict, model_config: dict) -> list[dict]:
    tasks = []
    seeds = config["inner_seed_resolution"]["seed_by_tissue"]
    for tissue_index, tissue in enumerate(config["tissue_order"]):
        for drug_id in config["drug_ids"]:
            for model in config["models"]:
                tasks.append(
                    dict(
                        drug_id=int(drug_id), repeat=0, outer_fold=tissue_index,
                        tissue_index=tissue_index, heldout_tissue=tissue,
                        inner_seed=int(seeds[tissue]), model=model, permutation_id=0,
                        **model_spec(model, model_config),
                    )
                )
    return tasks


def checkpoint_name(task: dict) -> str:
    return f"t{task['tissue_index']:02d}_d{task['drug_id']}_{task['model']}.json"


def validate_checkpoint(path: Path, task: dict, freeze_sha: str) -> bool:
    try:
        return base.validate_resume_checkpoint(path, task, freeze_sha)
    except Exception:
        return False


def compute_split(data: pd.DataFrame, task: dict, model_config: dict) -> tuple[np.ndarray, np.ndarray, list, dict]:
    tissue = data.tissue.astype(str).to_numpy()
    heldout = np.flatnonzero(tissue == task["heldout_tissue"])
    train = np.flatnonzero(tissue != task["heldout_tissue"])
    if not len(heldout) or not len(train):
        raise ValueError(f"Empty LOTO train/test partition: {task['drug_id']}/{task['heldout_tissue']}")
    groups = data.biological_group_id.astype(str).to_numpy()
    if set(groups[train]) & set(groups[heldout]):
        raise ValueError("A biological group crosses the held-out tissue boundary")
    train_frame = data.iloc[train]
    captured: list[str] = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        inner = base.make_grouped_splits(
            train_frame.biological_group_id.astype(str).to_numpy(),
            train_frame.tissue.astype(str).to_numpy(),
            int(model_config["inner_folds"]), int(task["inner_seed"]),
        )
    captured.extend(str(item.message) for item in caught)
    validation_assignment: dict[str, int] = {}
    for fold, (inner_train, inner_valid) in enumerate(inner):
        inner_train = np.asarray(inner_train, dtype=int)
        inner_valid = np.asarray(inner_valid, dtype=int)
        if set(groups[train][inner_train]) & set(groups[train][inner_valid]):
            raise ValueError("An inner training/validation fold shares biological groups")
        if set(inner_train) | set(inner_valid) != set(range(len(train))) or set(inner_train) & set(inner_valid):
            raise ValueError("An inner fold does not partition the outer training samples")
        for i in inner_valid:
            model_id = str(train_frame.iloc[int(i)].canonical_model_id)
            if model_id in validation_assignment:
                raise ValueError("A training sample appears in multiple inner validation folds")
            validation_assignment[model_id] = fold
    if len(validation_assignment) != len(train):
        raise ValueError("Inner validation folds do not cover all outer training samples exactly once")
    indices = data.expression_row_index.astype(int).to_numpy()
    details = {
        "n_train": int(len(train)),
        "n_test": int(len(heldout)),
        "n_train_groups": int(len(set(groups[train]))),
        "n_test_groups": int(len(set(groups[heldout]))),
        "train_ids_sha256": hash_ids(data.iloc[train].canonical_model_id),
        "test_ids_sha256": hash_ids(data.iloc[heldout].canonical_model_id),
        "train_expression_indices_sha256": hash_ids(indices[train]),
        "test_expression_indices_sha256": hash_ids(indices[heldout]),
        "inner_validation_assignment_sha256": hash_inner_assignment(validation_assignment),
        "inner_validation_counts": [
            sum(value == fold for value in validation_assignment.values())
            for fold in range(int(model_config["inner_folds"]))
        ],
        "splitter_warnings": captured,
    }
    return train, heldout, inner, {**details, "validation_assignment": validation_assignment}


def preflight(config: dict, model_config: dict, run: Path, freeze_sha: str,
              revision: str, source_hashes: dict, verification: dict) -> dict:
    if run.exists():
        raise FileExistsError(f"LOTO run path already exists: {run}; do not overwrite prior evidence")
    run.mkdir(parents=True)
    (run / "checkpoints").mkdir()
    source_response = SOURCE_RUN / "development_responses.csv"
    shutil.copyfile(source_response, run / "development_responses.csv")
    responses = pd.read_csv(
        run / "development_responses.csv",
        dtype={"gdsc_drug_id": str, "canonical_model_id": str, "biological_group_id": str, "tissue": str},
        keep_default_na=False, float_precision="round_trip",
    )
    if len(responses) != 4967 or sorted(responses.gdsc_drug_id.astype(int).unique()) != sorted(config["drug_ids"]):
        raise ValueError("Prepared source response inventory does not match the frozen 8-drug panel")
    if not np.isfinite(responses.response_auc.to_numpy(dtype=np.float64)).all():
        raise ValueError("Nonfinite development response")
    if responses.duplicated(["gdsc_drug_id", "canonical_model_id"]).any():
        raise ValueError("Duplicate drug × model response row")
    p5e_coverage = json_load(ROOT / "config/p5e_protocol_v1.lock.json")["coverage_freeze"]["drug_tissue_group_counts"]
    for drug in config["drug_ids"]:
        drug_frame = responses.loc[responses.gdsc_drug_id.eq(str(drug))]
        if drug_frame.biological_group_id.nunique() != len(drug_frame):
            raise ValueError(f"Drug {drug} has duplicate rows for an independent biological group")
        observed_counts = drug_frame.groupby("tissue").biological_group_id.nunique().to_dict()
        frozen_counts = p5e_coverage[str(drug)]
        if observed_counts != frozen_counts:
            raise ValueError(f"Drug/tissue group counts differ from frozen P5e coverage for drug {drug}")
    observed_tissues = sorted(responses.tissue.astype(str).unique().tolist())
    if observed_tissues != config["tissue_order"]:
        raise ValueError(f"Frozen tissue order/inventory mismatch: {observed_tissues}")

    identity = responses.drop_duplicates("canonical_model_id").sort_values("canonical_model_id")
    if identity.canonical_model_id.nunique() != 626 or identity.biological_group_id.nunique() != 626:
        raise ValueError("The frozen biological group inventory changed")
    outer_rows = []
    for tissue_index, heldout_tissue in enumerate(config["tissue_order"]):
        for row in identity.itertuples(index=False):
            outer_rows.append(
                dict(
                    canonical_model_id=row.canonical_model_id,
                    biological_group_id=row.biological_group_id,
                    tissue=row.tissue,
                    repeat=0,
                    outer_fold=tissue_index,
                    inner_fold=-1,
                    role="test" if row.tissue == heldout_tissue else "train",
                    heldout_tissue=heldout_tissue,
                )
            )
    pd.DataFrame(outer_rows).to_csv(run / "global_fold_manifest.csv", index=False, encoding="utf-8-sig", lineterminator="\n")

    tasks = task_list(config, model_config)
    if len(tasks) != config["expected_contract"]["tasks"]:
        raise ValueError("LOTO task budget mismatch")
    plan_rows = []
    inner_rows = []
    split_cache: dict[tuple[int, str], dict] = {}
    coverage_counts = {}
    for tissue_index, heldout_tissue in enumerate(config["tissue_order"]):
        for drug_id in config["drug_ids"]:
            subset = responses.loc[responses.gdsc_drug_id.eq(str(drug_id))].sort_values("canonical_model_id").reset_index(drop=True)
            task_seed = config["inner_seed_resolution"]["seed_by_tissue"][heldout_tissue]
            split_task = dict(drug_id=int(drug_id), heldout_tissue=heldout_tissue, tissue_index=tissue_index, inner_seed=task_seed)
            train, test, inner, details = compute_split(subset, split_task, model_config)
            assignment = details.pop("validation_assignment")
            split_cache[(int(drug_id), heldout_tissue)] = details
            coverage_counts[f"{drug_id}|{heldout_tissue}"] = details["n_test"]
            for model_id, validation_fold in sorted(assignment.items()):
                row = subset.loc[subset.canonical_model_id.eq(model_id)].iloc[0]
                inner_rows.append(
                    dict(
                        drug_id=int(drug_id), tissue_index=tissue_index,
                        heldout_tissue=heldout_tissue, canonical_model_id=model_id,
                        biological_group_id=row.biological_group_id, tissue=row.tissue,
                        validation_fold=int(validation_fold), inner_seed=int(task_seed),
                    )
                )
    pd.DataFrame(inner_rows).to_csv(run / "loto_inner_fold_manifest.csv", index=False, encoding="utf-8-sig", lineterminator="\n")
    for task in tasks:
        details = split_cache[(task["drug_id"], task["heldout_tissue"])]
        plan_rows.append({"task": task, "split": details})
    write_jsonl(run / "task_plan.jsonl", plan_rows)

    expected_fits = sum(1 if task["model"] == "tissue_null" else 31 for task in tasks)
    expected_oof = sum(row["split"]["n_test"] for row in plan_rows)
    contract = config["expected_contract"]
    if (expected_fits, expected_oof) != (contract["fits"], contract["oof_predictions"]):
        raise ValueError(f"LOTO fit/OOF workload mismatch: {expected_fits}/{expected_oof}")

    source_hashes_in = {
        "source_development_responses_sha256": sha(source_response),
        "prepared_development_responses_sha256": sha(run / "development_responses.csv"),
    }
    if source_hashes_in["source_development_responses_sha256"] != source_hashes_in["prepared_development_responses_sha256"]:
        raise ValueError("Prepared response copy is not byte-identical to the frozen source")
    prepared_hashes = {name: sha(run / name) for name in PREPARED_FILES if name != "preflight_manifest.json"}
    split_summary = {
        "tasks": len(tasks), "fits": expected_fits, "oof_predictions": expected_oof,
        "outer_drug_tissue_splits": len(coverage_counts),
        "outer_tissue_sample_counts": coverage_counts,
        "inner_splitter_warning_contexts": sum(bool(value["splitter_warnings"]) for value in split_cache.values()),
        "inner_splitter_warning_count": sum(len(value["splitter_warnings"]) for value in split_cache.values()),
        "inner_splitter_warning_message_counts": dict(Counter(
            message for value in split_cache.values() for message in value["splitter_warnings"]
        )),
        "total_inner_manifest_rows": len(inner_rows),
        "unknown_tissue_equals_test_samples": True,
        "training_includes_heldout_tissue": False,
    }
    manifest = {
        "status": "preflight_completed_no_model_fits",
        "freeze_sha256": freeze_sha,
        "execution_commit": revision,
        "locked_source_hashes": source_hashes,
        "prepared_input_hashes": prepared_hashes,
        "source_response_hashes": source_hashes_in,
        "development_input_verification": verification,
        "seed_by_tissue": config["inner_seed_resolution"]["seed_by_tissue"],
        "split_summary": split_summary,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(run / "preflight_manifest.json", manifest)
    prepared_hashes["preflight_manifest.json"] = sha(run / "preflight_manifest.json")
    context = {
        "freeze_sha256": freeze_sha,
        "execution_commit": revision,
        "source_hashes": source_hashes,
        "prepared_input_hashes": prepared_hashes,
        "execution_config": config,
        "model_config": model_config,
        "task_count": len(tasks),
    }
    atomic_json(run / "worker_context.json", context)
    context_hash = sha(run / "worker_context.json")
    atomic_json(
        run / "loto_metadata.json",
        {
            "status": "preflight_completed",
            "completed": False,
            "run_name": config["run_name"],
            "freeze_sha256": freeze_sha,
            "execution_commit": revision,
            "worker_context_sha256": context_hash,
            "prepared_input_hashes": prepared_hashes,
            "source_hashes": source_hashes,
            "expected_tasks": contract["tasks"],
            "expected_fit_events": contract["fits"],
            "expected_oof_predictions": contract["oof_predictions"],
            "preflight": split_summary,
            "sessions": [],
            "checkpoint_hashes": {},
            "external_response_accessed": False,
            "development_input_verification": verification,
        },
    )
    return {**split_summary, "run_dir": str(run), "worker_context_sha256": context_hash}


def load_worker_context(run: Path, expected_sha: str) -> dict:
    context_path = run / "worker_context.json"
    if sha(context_path) != expected_sha:
        raise ValueError("Worker context checksum mismatch")
    context = json_load(context_path)
    for name, expected in context["prepared_input_hashes"].items():
        if sha(run / name) != expected:
            raise ValueError("Prepared LOTO input changed: " + name)
    for name, expected in context["source_hashes"].items():
        if sha(ROOT / name) != expected:
            raise ValueError("Frozen LOTO source changed: " + name)
    return context


def execute_task(task_index: int, run: Path, context: dict) -> bool:
    config = context["execution_config"]
    model_config = context["model_config"]
    plan_rows = read_jsonl(run / "task_plan.jsonl")
    tasks = [item["task"] for item in plan_rows]
    plan = plan_rows[task_index]
    task = tasks[task_index]
    started = time.perf_counter()
    process = psutil.Process()
    cpu0 = process.cpu_times()
    samples = []
    access: set = set()
    blocked: list = []
    stop = threading.Event()
    obj = dict(
        schema_version="1.0", task=task, freeze_sha256=context["freeze_sha256"],
        status="failed", oof_rows=[], result=None, fit_events=[], error=None,
        expected_split=plan["split"],
        execution_commit=context["execution_commit"], source_hashes=context["source_hashes"],
    )

    def sample():
        while not stop.is_set():
            try:
                mem = process.memory_info()
                samples.append(dict(rss=int(mem.rss), private=int(getattr(mem, "private", mem.vms)), **system_memory()))
            except Exception:
                pass
            stop.wait(0.1)

    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    try:
        base.worker_init(str(run), model_config, context["freeze_sha256"])
        worker = base.WORKER
        access, blocked = worker["access"], worker["blocked"]
        worker["responses"] = pd.read_csv(
            run / "development_responses.csv",
            dtype={"gdsc_drug_id": str, "canonical_model_id": str, "biological_group_id": str, "tissue": str},
            keep_default_na=False, float_precision="round_trip",
        )
        data, x = base.get_drug(task["drug_id"])
        y = data.response_auc.to_numpy(dtype=np.float64)
        train, test, inner, details = compute_split(data, task, model_config)
        comparable = {key: value for key, value in details.items() if key != "validation_assignment"}
        if comparable != plan["split"]:
            raise ValueError("Worker recomputed split differs from preflight split plan")
        if len(inner) != int(config["inner_folds"]):
            raise ValueError("Unexpected inner-fold count")
        test_ids = sorted(data.iloc[test].canonical_model_id.astype(str).tolist())
        training_expression_ids = sorted(x[train, -1].astype(int).tolist())
        inner_training_expression_ids = [
            sorted(x[train[np.asarray(inner_train, dtype=int)], -1].astype(int).tolist())
            for inner_train, inner_valid in inner
        ]
        obj.update(
            expected_test_ids=test_ids,
            expected_training_expression_ids=training_expression_ids,
            inner_training_expression_ids=inner_training_expression_ids,
            assigned_response_sha256=hashlib.sha256(y.tobytes()).hexdigest(),
            inner_validation_assignment_sha256=details["inner_validation_assignment_sha256"],
        )
        events = obj["fit_events"]
        result = evaluate_repeated_outer_fold(
            x, y, data.biological_group_id.astype(str).to_numpy(), train, test, inner,
            task["model"], model_config, worker["stable"], audit_sink=events.append,
        )
        predictions = result.pop("predictions")
        if len(predictions) != len(test) or not np.isfinite(predictions).all():
            raise ValueError("Incomplete or nonfinite held-out predictions")
        if int(result["outer_unknown_tissue_samples"]) != len(test):
            raise ValueError("Every held-out sample must be an unknown tissue to the fitted model")
        if blocked:
            raise ValueError("A blocked external-data access attempt was recorded")
        for index, predicted in zip(test, predictions):
            row = data.iloc[int(index)]
            obj["oof_rows"].append(
                dict(
                    drug_id=task["drug_id"], drug_name=row.drug_name, model=task["model"],
                    repeat=0, outer_fold=task["outer_fold"], canonical_model_id=row.canonical_model_id,
                    biological_group_id=row.biological_group_id, tissue=row.tissue,
                    observed_auc=float(y[index]), predicted_auc=float(predicted), permutation_id=0,
                    selection_rule=task["selection_rule"], representation=task["representation"],
                    n_genes=task["n_genes"], heldout_tissue=task["heldout_tissue"],
                    tissue_index=task["tissue_index"],
                )
            )
        result.pop("fit_events", None)
        obj.update(status="completed", result=result)
    except BaseException:
        obj["error"] = traceback.format_exc()
    finally:
        stop.set()
        sampler.join()
        cpu = process.cpu_times()
        mem = process.memory_info()
        if not samples:
            samples.append(dict(rss=int(mem.rss), private=int(getattr(mem, "private", mem.vms)), **system_memory()))
        obj["resources"] = dict(
            wall_seconds=time.perf_counter() - started,
            cpu_seconds=cpu.user + cpu.system - cpu0.user - cpu0.system,
            sampled_peak_rss_bytes=max(item["rss"] for item in samples),
            sampled_peak_private_bytes=max(item["private"] for item in samples),
            minimum_available_ram_bytes=min(item["available_ram_bytes"] for item in samples),
            minimum_commit_free_bytes=min(item["system_commit_free_bytes"] for item in samples),
            windows_lifetime_peak_working_set_bytes=getattr(mem, "peak_wset", None),
            pid=process.pid, process_create_time=process.create_time(), controller_pid=process.ppid(),
            task_process_isolated=True, preimport_threads_applied=EARLY_THREADS,
            threadpools=threadpool_info(),
        )
        obj["data_guard"] = dict(external_response_accessed=False, blocked_attempts=len(blocked), accessed_path_count=len(access))
        atomic_json(run / "checkpoints" / checkpoint_name(task), obj)
    return obj["status"] == "completed"


def aggregate(run: Path, tasks: list[dict], freeze_sha: str) -> dict:
    counts = dict(tasks=0, fit_events=0, oof_predictions=0)
    with (run / "oof_predictions.csv").open("w", newline="", encoding="utf-8-sig") as oof, \
         (run / "fit_events.jsonl").open("w", encoding="utf-8", newline="\n") as events_file, \
         (run / "outer_fit_audits.jsonl").open("w", encoding="utf-8", newline="\n") as audits_file:
        writer = csv.DictWriter(oof, fieldnames=OOF_FIELDS)
        writer.writeheader()
        for task in tasks:
            path = run / "checkpoints" / checkpoint_name(task)
            if not validate_checkpoint(path, task, freeze_sha):
                raise ValueError("Failed or tampered LOTO task cannot be omitted: " + path.name)
            obj = json_load(path)
            counts["tasks"] += 1
            writer.writerows(obj["oof_rows"])
            counts["oof_predictions"] += len(obj["oof_rows"])
            for event in obj["fit_events"]:
                events_file.write(json.dumps({**task, **event}, ensure_ascii=False, allow_nan=False) + "\n")
                counts["fit_events"] += 1
            audits_file.write(
                json.dumps(
                    dict(task=task, **obj["result"], resources=obj["resources"], data_guard=obj["data_guard"]),
                    ensure_ascii=False, allow_nan=False,
                ) + "\n"
            )
    return counts


def worker_command(run: Path, task_index: int, context_sha: str) -> list[str]:
    return [sys.executable, "-u", str(Path(__file__).resolve()), "--run-dir", str(run),
            "--worker-task", str(task_index), "--context-sha", context_sha]


def run_controller(run: Path, resume: bool, lock: dict, freeze_sha: str, revision: str,
                   config: dict, model_config: dict) -> int:
    metadata_path = run / "loto_metadata.json"
    meta = json_load(metadata_path)
    if meta["freeze_sha256"] != freeze_sha or meta["execution_commit"] != revision:
        raise ValueError("LOTO resume lock or source revision mismatch")
    if not resume and meta["status"] != "preflight_completed":
        raise ValueError("Initial execution requires a completed no-fit preflight")
    if resume and meta["status"] not in {"resource_paused", "failed"}:
        raise ValueError("Only a paused/failed LOTO run can be resumed")
    for name, expected in meta["prepared_input_hashes"].items():
        if sha(run / name) != expected:
            raise ValueError("Prepared LOTO input changed: " + name)
    tasks = task_list(config, model_config)
    freeze_sha = meta["freeze_sha256"]
    context_sha = meta["worker_context_sha256"]
    registered = meta["checkpoint_hashes"]
    pending = []
    for index, task in enumerate(tasks):
        path = run / "checkpoints" / checkpoint_name(task)
        if path.exists() and validate_checkpoint(path, task, freeze_sha):
            if path.name in registered and registered[path.name] != sha(path):
                raise ValueError("Registered LOTO checkpoint changed: " + path.name)
            registered[path.name] = sha(path)
            continue
        if path.exists():
            failed_dir = run / "failed_attempts"
            failed_dir.mkdir(exist_ok=True)
            dest = failed_dir / f"{path.stem}_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"
            os.replace(path, dest)
            registered.pop(path.name, None)
        pending.append((index, task))

    session = dict(
        started_utc=datetime.now(timezone.utc).isoformat(), pending_at_start=len(pending),
        completed_tasks_this_session=0, status="running", max_inflight_observed=1,
        resource_wait_seconds=0.0, child_wall_seconds_sum=0.0,
    )
    meta["sessions"].append(session)
    meta.update(status="running", completed=False, completed_tasks=len(tasks) - len(pending))
    atomic_json(metadata_path, meta)
    started = time.perf_counter()
    peak = private_peak = 0
    minimum = system_memory()["available_ram_bytes"]
    commit_minimum = system_memory()["system_commit_free_bytes"]
    exit_code = 0
    child = None
    session_no = len(meta["sessions"])
    resource_path = run / f"resource_samples_session{session_no}.csv"
    fields = ["elapsed_seconds", "process_tree_rss_bytes", "process_tree_private_bytes", "process_tree_count",
              "available_ram_bytes", "system_commit_free_bytes", "system_commit_bytes", "system_commit_limit_bytes", "commit_measurement", "inflight"]
    try:
        with resource_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()

            def record(inflight: int) -> dict:
                nonlocal peak, private_peak, minimum, commit_minimum
                memory = process_tree_memory()
                peak = max(peak, memory["process_tree_rss_bytes"])
                private_peak = max(private_peak, memory["process_tree_private_bytes"])
                minimum = min(minimum, memory["available_ram_bytes"])
                commit_minimum = min(commit_minimum, memory["system_commit_free_bytes"])
                writer.writerow({"elapsed_seconds": time.perf_counter() - started, "inflight": inflight, **memory})
                handle.flush()
                return memory

            for index, task in pending:
                wait_start = time.perf_counter()
                while True:
                    memory = record(0)
                    reasons = pressure_reasons(memory, config["resources"])
                    if not reasons:
                        break
                    waited = time.perf_counter() - wait_start
                    print(json.dumps(dict(event="resource_wait", reasons=reasons, wait_seconds=waited, **memory)), flush=True)
                    if waited >= config["resources"]["pressure_max_wait_seconds"]:
                        session["resource_wait_seconds"] += waited
                        session["status"] = "resource_paused"
                        meta["status"] = "resource_paused"
                        exit_code = 75
                        print(json.dumps(dict(event="resource_paused", completed_tasks=meta["completed_tasks"], pending_task=task)), flush=True)
                        break
                    time.sleep(config["resources"]["pressure_poll_seconds"])
                if exit_code:
                    break
                session["resource_wait_seconds"] += time.perf_counter() - wait_start
                child_start = time.perf_counter()
                env = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8", REO_EXECUTION_COMMIT=revision)
                for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
                    env[key] = "1"
                child = subprocess.Popen(worker_command(run, index, context_sha), cwd=ROOT, env=env)
                last_heartbeat = -60.0
                while child.poll() is None:
                    record(1)
                    elapsed = time.perf_counter() - started
                    if elapsed - last_heartbeat >= 60:
                        print(json.dumps(dict(event="loto_heartbeat", completed_tasks=meta["completed_tasks"], total_tasks=len(tasks), elapsed_seconds=round(elapsed, 1), inflight=1)), flush=True)
                        last_heartbeat = elapsed
                    time.sleep(1)
                child_wall = time.perf_counter() - child_start
                session["child_wall_seconds_sum"] += child_wall
                path = run / "checkpoints" / checkpoint_name(task)
                if not path.exists():
                    atomic_json(path, dict(schema_version="1.0", task=task, freeze_sha256=freeze_sha,
                        status="failed", oof_rows=[], result=None, fit_events=[],
                        error=f"Worker exited {child.returncode} without a complete checkpoint",
                        resources=dict(wall_seconds=child_wall, sampled_peak_rss_bytes=0),
                        data_guard=dict(external_response_accessed=False, blocked_attempts=0)))
                obj = json_load(path)
                registered[path.name] = sha(path)
                if child.returncode or not validate_checkpoint(path, task, freeze_sha):
                    print(json.dumps(dict(event="loto_task_failed", task=task, error=obj.get("error"), returncode=child.returncode)), flush=True)
                    raise RuntimeError(obj.get("error") or "LOTO worker exited nonzero")
                session["completed_tasks_this_session"] += 1
                meta["completed_tasks"] = len(tasks) - len(pending) + session["completed_tasks_this_session"]
                print(json.dumps(dict(event="loto_task_completed", task=task,
                    completed_tasks=meta["completed_tasks"], total_tasks=len(tasks),
                    wall_seconds=obj["resources"]["wall_seconds"], fit_count=len(obj["fit_events"]),
                    sampled_peak_rss_bytes=obj["resources"]["sampled_peak_rss_bytes"])), flush=True)
                del obj
                atomic_json(metadata_path, meta)
                child = None
        if not exit_code:
            counts = aggregate(run, tasks, freeze_sha)
            expected = config["expected_contract"]
            if (counts["tasks"], counts["fit_events"], counts["oof_predictions"]) != (expected["tasks"], expected["fits"], expected["oof_predictions"]):
                raise ValueError(f"LOTO complete workload contract mismatch: {counts}")
            meta.update(
                completed=True, status="completed", output_counts=counts,
                output_hashes={name: sha(run / name) for name in ("oof_predictions.csv", "fit_events.jsonl", "outer_fit_audits.jsonl")},
                python_version=platform.python_version(),
                package_versions={name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scikit-learn", "scipy", "psutil")},
            )
            session["status"] = "completed"
    except BaseException:
        if child is not None and child.poll() is None:
            child.wait()
        meta.update(completed=False, status="failed", failure=traceback.format_exc())
        session["status"] = "failed"
        raise
    finally:
        session.update(
            wall_seconds=time.perf_counter() - started,
            finished_utc=datetime.now(timezone.utc).isoformat(),
            sampled_process_tree_peak_rss_bytes=peak,
            sampled_process_tree_peak_private_bytes=private_peak,
            minimum_available_ram_bytes=minimum,
            minimum_commit_free_bytes=commit_minimum,
        )
        atomic_json(metadata_path, meta)
    if not exit_code:
        print(json.dumps(dict(event="loto_run_completed", counts=meta["output_counts"], wall_seconds=session["wall_seconds"])), flush=True)
    return exit_code


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--preflight-only", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--resume", action="store_true")
    parser.add_argument("--worker-task", type=int)
    parser.add_argument("--context-sha")
    args = parser.parse_args()
    if args.worker_task is not None:
        if not args.run_dir or not args.context_sha:
            raise ValueError("Worker requires --run-dir and --context-sha")
        run = args.run_dir.resolve()
        run.relative_to(ROOT / "runs")
        context = load_worker_context(run, args.context_sha)
        tasks = read_jsonl(run / "task_plan.jsonl")
        if not 0 <= args.worker_task < len(tasks):
            raise IndexError("Worker task index out of bounds")
        return 0 if execute_task(args.worker_task, run, context) else 1

    lock, freeze_sha, revision, source_hashes = load_execution_lock()
    config = json_load(EXEC_CONFIG_PATH)
    model_config = json_load(MODEL_CONFIG_PATH)
    run = (args.run_dir or ROOT / "runs" / config["run_name"]).resolve()
    run.relative_to(ROOT / "runs")
    if not (args.preflight_only or args.run or args.resume):
        raise ValueError("Choose exactly one mode: --preflight-only, --run, or --resume")
    access, blocked = install_data_guard()
    if blocked:
        raise ValueError("A blocked external-data access attempt occurred before LOTO setup")
    verification = verify_development_inputs("v2")
    if args.preflight_only:
        summary = preflight(config, model_config, run, freeze_sha, revision, source_hashes, verification)
        print(json.dumps({"event": "loto_preflight_completed", **summary}, ensure_ascii=False))
        return 0
    meta_path = run / "loto_metadata.json"
    if not meta_path.exists():
        raise FileNotFoundError("A successful no-fit LOTO preflight is required before fitting")
    meta = json_load(meta_path)
    if not meta.get("development_input_verification"):
        raise ValueError("LOTO run lacks development input verification")
    meta["controller_data_guard"] = dict(external_response_accessed=False, blocked_attempts=len(blocked), accessed_path_count=len(access))
    atomic_json(meta_path, meta)
    return run_controller(run, args.resume, lock, freeze_sha, revision, config, model_config)


if __name__ == "__main__":
    sys.exit(main())
