"""Create the second P6 lock after final models and response-free preflight.

The CTRPv2 response file and response index are intentionally never opened by
this program. Its expected response SHA is copied from already-frozen curation
metadata; source verification happens only inside the later hard-gated executor.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from ordinal_repeated_runtime import sha  # noqa: E402
from freeze_external_validation_v1_1 import (  # noqa: E402
    EXPECTED_DRUGS, EXPECTED_MODELS, build_lock, install_pre_external_guard,
)

FIRST_LOCK = ROOT / "config/external_validation_v1_1.lock.json"
SECOND_LOCK = ROOT / "config/external_execution_v1_1.lock.json"
FINAL_RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-full-development-final-v1_1"
EXTERNAL_RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-external-v1_1"
RAW_RESPONSE_NAME = "data_raw/ctrp/v20.data.curves_post_qc.txt"
SOURCE_FIELD = "area_under_curve"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def verify_first_lock() -> tuple[dict, str]:
    first = _json(FIRST_LOCK)
    if first.get("status") != "frozen_pre_final_fit" or first.get("external_response_accessed") is not False:
        raise ValueError("First P6 lock is not valid for final fitting")
    current = build_lock(first["execution_commit"])
    for key in ("status", "execution_commit", "execution_files", "files", "contract", "environment"):
        if first.get(key) != current[key]:
            raise ValueError("First P6 lock changed or its bound source/input differs: " + key)
    return first, sha(FIRST_LOCK)


def model_hashes_from_manifest(path: Path) -> tuple[dict[str, str], int]:
    """Independently verify every planned drug/arm artifact and fit budget."""
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]
    expected = {(drug, arm) for drug in EXPECTED_DRUGS for arm in EXPECTED_MODELS}
    seen: set[tuple[int, str]] = set()
    hashes: dict[str, str] = {}
    fits = 0
    for row in rows:
        task = row.get("task", {})
        key = (task.get("drug_id"), task.get("model"))
        if key not in expected or key in seen:
            raise ValueError("Final-model manifest has missing, duplicate or unexpected task")
        seen.add(key)
        expected_fits = 1 if key[1] == "tissue_null" else 31
        if row.get("fit_count") != expected_fits or row.get("fresh_process_reload_verified") is not True:
            raise ValueError("Final-model task lacks expected fit/reload evidence")
        if row.get("artifact_path") != f"models/d{key[0]}_{key[1]}.joblib":
            raise ValueError("Final-model artifact path differs from fixed task identity")
        artifact = (FINAL_RUN / row["artifact_path"]).resolve()
        artifact.relative_to((FINAL_RUN / "models").resolve())
        name = artifact.relative_to(ROOT.resolve()).as_posix()
        if not artifact.is_file() or row.get("artifact_sha256") != sha(artifact):
            raise ValueError("Final-model artifact missing or SHA mismatch: " + name)
        hashes[name] = row["artifact_sha256"]
        fits += expected_fits
    if seen != expected or len(hashes) != 64 or fits != 1744:
        raise ValueError("Incomplete fixed 64-artifact/1,744-fit final-model panel")
    return hashes, fits


def _prediction_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or ())
        for forbidden in ("observed_auc", "area_under_curve", "response_auc", "ctrp_auc"):
            if forbidden in columns:
                raise ValueError("Response column leaked into response-free prediction manifest")
        required = {"drug_id", "model", "canonical_model_id", "status", "predicted_auc"}
        if not required <= columns:
            raise ValueError("Response-free prediction manifest lacks required columns")
        counts: dict[tuple[int, str], int] = {}
        for row in reader:
            key = (int(row["drug_id"]), row["model"])
            counts[key] = counts.get(key, 0) + 1
    if set(counts) != {(drug, arm) for drug in EXPECTED_DRUGS for arm in EXPECTED_MODELS}:
        raise ValueError("Response-free prediction manifest omits a fixed arm")
    if sum(counts.values()) != 9480:
        raise ValueError("Response-free prediction row count differs from 1,185 x 8")
    return sum(counts.values())


def build_execution_lock() -> dict:
    first, first_sha = verify_first_lock()
    final_meta_path = FINAL_RUN / "final_model_metadata.json"
    final_manifest_path = FINAL_RUN / "model_manifest.jsonl"
    final_summary_path = FINAL_RUN / "final_model_summary.json"
    final_audit_path = FINAL_RUN / "final_model_independent_audit.json"
    preflight_path = EXTERNAL_RUN / "preflight_metadata.json"
    prediction_path = EXTERNAL_RUN / "prediction_manifest.csv"
    final = _json(final_meta_path)
    summary = _json(final_summary_path)
    audit = _json(final_audit_path)
    preflight = _json(preflight_path)
    model_hashes, fit_count = model_hashes_from_manifest(final_manifest_path)
    if (final.get("status"), final.get("completed"), final.get("completed_tasks")) != ("completed", True, 64):
        raise ValueError("Full-development final run is incomplete")
    if final.get("freeze_sha256") != first_sha or final.get("external_response_accessed") is not False:
        raise ValueError("Full-development final run first-lock/response contract failed")
    if summary.get("artifacts") != 64 or summary.get("fits") != fit_count or summary.get("model_manifest_sha256") != sha(final_manifest_path):
        raise ValueError("Final-model aggregate differs from verified model manifest")
    if (audit.get("status"), audit.get("artifacts"), audit.get("fit_events"),
        audit.get("external_response_accessed")) != ("completed_contract_validated", 64, 1744, False):
        raise ValueError("Independent final-model audit gate has not passed")
    if audit.get("pre_fit_lock_sha256") != first_sha or audit.get("model_manifest_sha256") != sha(final_manifest_path):
        raise ValueError("Independent final-model audit source binding differs")
    if preflight.get("status") != "ready_no_response" or preflight.get("external_response_accessed") is not False:
        raise ValueError("External prediction preflight is not response-free and ready")
    if (preflight.get("first_lock_sha256"), preflight.get("final_model_metadata_sha256"),
        preflight.get("final_model_manifest_sha256")) != (first_sha, sha(final_meta_path), sha(final_manifest_path)):
        raise ValueError("External prediction preflight final-model binding differs")
    if preflight.get("model_artifact_hashes") != model_hashes:
        raise ValueError("External prediction preflight artifact hashes differ")
    prediction_rows = _prediction_rows(prediction_path)
    if preflight.get("prediction_rows") != prediction_rows or preflight.get("prediction_manifest_sha256") != sha(prediction_path):
        raise ValueError("External preflight prediction manifest hash/count differs")
    for name in ("src/run_external_validation_v1.py", "src/external_validation_v1.py"):
        if preflight.get("source_hashes", {}).get(name) != first["files"][name]:
            raise ValueError("External preflight scoring source differs from first lock: " + name)
    if preflight.get("external_manifest_sha256") != first["files"]["results/cohort_v2/external_primary_sample_manifest.csv"]:
        raise ValueError("External preflight roster differs from first lock")
    if preflight.get("global_split_manifest_sha256") != first["files"]["results/cohort_v2/global_split_manifest.csv"]:
        raise ValueError("External preflight global split differs from first lock")
    if preflight.get("expression_tpm_sha256") != first["files"]["data_processed/expression_tpm.npy"] or \
       preflight.get("expression_genes_sha256") != first["files"]["data_processed/expression_genes.tsv"]:
        raise ValueError("External preflight expression input differs from first lock")
    curation = _json(ROOT / "data_processed/ctrp_curation_summary.json")
    if curation.get("response_store") != RAW_RESPONSE_NAME or curation.get("response_field") != SOURCE_FIELD:
        raise ValueError("Frozen CTRP source-field metadata differs")
    expected_raw_sha = curation.get("raw_curve_sha256")
    if not isinstance(expected_raw_sha, str) or len(expected_raw_sha) != 64:
        raise ValueError("Frozen CTRP response source SHA is unavailable")
    source_hashes = {name: first["files"][name] for name in first["execution_files"]}
    return {
        "schema_version": "1.0",
        "status": "frozen_pre_external_response",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "first_lock_path": FIRST_LOCK.relative_to(ROOT).as_posix(),
        "first_lock_sha256": first_sha,
        "execution_commit": first["execution_commit"],
        "source_hashes": source_hashes,
        "final_model_metadata_sha256": sha(final_meta_path),
        "final_model_manifest_sha256": sha(final_manifest_path),
        "final_model_independent_audit_sha256": sha(final_audit_path),
        "preflight_metadata_sha256": sha(preflight_path),
        "prediction_manifest_sha256": sha(prediction_path),
        "model_artifact_hashes": model_hashes,
        "model_count": 64,
        "fit_events": fit_count,
        "prediction_rows": prediction_rows,
        "response_source_path": RAW_RESPONSE_NAME,
        "response_source_sha256": expected_raw_sha,
        "response_column": SOURCE_FIELD,
        "external_response_accessed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--ready", action="store_true")
    modes.add_argument("--lock", action="store_true")
    modes.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    _, blocked = install_pre_external_guard()
    current = build_execution_lock()
    if args.verify_only:
        saved = _json(SECOND_LOCK)
        for key, value in current.items():
            if key != "created_utc" and saved.get(key) != value:
                raise ValueError("Second P6 lock differs from current frozen artifacts: " + key)
        status = "verified"
    elif args.lock:
        with SECOND_LOCK.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(current, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
        status = "locked"
    else:
        status = "ready"
    if blocked:
        raise ValueError("CTRP response guard observed a forbidden access attempt")
    print(json.dumps({"status": status, "model_count": 64, "fit_events": 1744,
                      "prediction_rows": current["prediction_rows"],
                      "second_lock_sha256": sha(SECOND_LOCK) if SECOND_LOCK.exists() else None,
                      "external_response_accessed": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
