"""Freeze the P6 v1.1 guard amendment before any full-development fit.

This program reads only committed protocol/code, frozen development inputs and
response-free cohort metadata. It installs the development CTRP data guard.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ordinal_repeated_runtime import sha  # noqa: E402
from run_development_pilot import verify_development_inputs  # noqa: E402
from project_paths import original_input  # noqa: E402


EXPECTED_DRUGS = [1931, 1036, 1191, 1930, 1199, 1091, 1010, 1922]
EXPECTED_RESPONSE_SOURCE = "data_raw/ctrp/v20.data.curves_post_qc.txt"
EXPECTED_RESPONSE_FIELD = "area_under_curve"
EXPECTED_RESPONSE_SHA256 = "b97bb9004f97ff7d5910c0f66c6d26163a2ce88529eaefc5df0e51624e5ac6a9"
PREDECESSOR_LOCK = "config/external_validation_v1.lock.json"
PREDECESSOR_LOCK_SHA256 = "dcb06ab198aa52db6f264c4b33e7b6b352740b97ef15f0c73722ed33b0b46594"
EXPECTED_MODELS = ["tissue_null", "expression_strong"] + [
    rule + "_" + representation
    for rule in ("mad", "response_association")
    for representation in ("expression", "rank", "pair_full")
]
CONTRACT = {
    "final_model_tasks": 64,
    "final_fit_events": 1744,
    "development_groups": 626,
    "development_drug_model_observations": 4967,
    "external_eligible_groups": 175,
    "external_drug_model_observations": 1185,
}

EXECUTION_FILES = [
    "config/external_validation_v1_1.json",
    "docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md",
    "reports/阶段12_P6端点剂量无响应审计.md",
    "reports/阶段12d_P6_v1首锁后零拟合启动故障.md",
    "reports/阶段12e_P6_v1_1冻结前独立审查.md",
    "reports/阶段12f_P6_v1_1首锁前验收.md",
    "src/final_model_v1.py",
    "src/run_final_model_v1.py",
    "src/external_validation_v1.py",
    "src/run_external_validation_v1.py",
    "src/analyze_formal_results.py",
    "scripts/freeze_external_validation_v1_1.py",
    "scripts/validate_final_model_v1.py",
    "scripts/freeze_external_execution_v1.py",
    "src/curate_ctrp.py",
    "src/curate_drugs.py",
    "src/freeze_cohorts.py",
    "src/ordinal_repeated_models.py",
    "src/ordinal_diagnostic_models.py",
    "src/ordinal_repeated_runtime.py",
    "src/run_formal_evaluation.py",
    "src/run_development_pilot.py",
    "src/project_paths.py",
    "tests/test_final_model_v1.py",
    "tests/test_external_validation_v1.py",
    "tests/test_run_external_validation_v1.py",
    "tests/test_external_freeze_v1.py",
    "tests/test_external_freeze_v1_1.py",
    "tests/test_validate_final_model_v1.py",
    "tests/test_external_execution_freeze_v1.py",
    "config/ordinal_repeated_v1.json",
]
EXECUTION_FILES += [
    p.relative_to(ROOT).as_posix() for p in sorted((ROOT / "src/modeling").glob("*.py"))
]

DIRECT_INPUTS = [
    PREDECESSOR_LOCK,
    "config/external_validation_v1.json",
    "docs/EXTERNAL_VALIDATION_PROTOCOL_V1.md",
    "config/ordinal_repeated_v1.lock.json",
    "config/cohort_freeze_v2.json",
    "config/cohort_freeze_v2.lock.json",
    "config/cohort_rules.json",
    "results/cohort_v2/development_sample_manifest.csv",
    "results/cohort_v2/external_primary_sample_manifest.csv",
    "results/cohort_v2/main_drugs.csv",
    "results/cohort_v2/global_split_manifest.csv",
    "results/cohort_v2/public_cohort_snapshot.json",
    "data_processed/drug_identity_map.csv",
    "data_processed/ctrp_curation_summary.json",
    "data_raw/ctrp/v20._COLUMNS.txt",
    "data_raw/ctrp/v20._README.txt",
    "data_raw/ctrp/v20.meta.per_compound.txt",
    "data_processed/expression_tpm.npy",
    "data_processed/expression_genes.tsv",
    "data_processed/expression_samples.tsv",
    "runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1/formal_metadata.json",
    "runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1/development_responses.csv",
]


def install_pre_external_guard() -> tuple[set[str], list[str]]:
    """Permit source metadata, but forbid CTRP responses and their row index."""
    ctrp_root = (ROOT / "data_raw/ctrp").resolve()
    metadata_allowlist = {
        (ROOT / name).resolve() for name in (
            "data_raw/ctrp/v20._COLUMNS.txt",
            "data_raw/ctrp/v20._README.txt",
            "data_raw/ctrp/v20.meta.per_compound.txt",
        )
    }
    forbidden_processed = {
        (ROOT / "data_processed" / name).resolve() for name in (
            "ctrp_curve_index.csv", "ctrp_pair_availability.csv", "ctrp_drug_availability.csv"
        )
    }
    legacy = original_input("legacy_ctrp_response").resolve()
    accessed: set[str] = set()
    blocked: list[str] = []

    def hook(event, args):
        if event != "open" or not args or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        try:
            path = Path(os.fsdecode(args[0])).resolve()
        except (ValueError, OSError):
            return
        if (path == legacy or path in forbidden_processed or
            (path.is_relative_to(ctrp_root) and path not in metadata_allowlist)):
            blocked.append(str(path))
            raise PermissionError("CTRP response/index access is forbidden before second P6 lock")
        if path.is_relative_to(ROOT / "data_processed") or path in metadata_allowlist:
            accessed.add(str(path))

    sys.addaudithook(hook)
    return accessed, blocked


def _committed_snapshot(names: list[str], revision: str | None = None) -> tuple[str, dict[str, str]]:
    rev = subprocess.check_output(
        ["git", "rev-parse", revision or "HEAD"], cwd=ROOT, text=True
    ).strip()
    found: dict[str, str] = {}
    for name in names:
        raw = subprocess.check_output(["git", "show", rev + ":" + name], cwd=ROOT)
        current = sha(ROOT / name)
        if hashlib.sha256(raw).hexdigest() != current:
            raise ValueError("Scientific source differs from committed revision: " + name)
        found[name] = current
    return rev, found


def _verify_configuration(config: dict) -> None:
    if config.get("protocol_status") != "frozen_pre_final_fit":
        raise ValueError("P6 protocol/config remain draft; do not freeze or fit")
    if config.get("pre_fit_gates_passed") is not True:
        raise ValueError("P6 pre-fit gates have not passed")
    if config.get("external_endpoint", {}).get("column_and_dose_contract_status") != "verified_pre_fit_without_external_response":
        raise ValueError("CTRP field/dose contract remains unverified")
    endpoint = config.get("external_endpoint", {})
    if (endpoint.get("source_field"), endpoint.get("raw_source_path"),
        endpoint.get("raw_source_sha256_from_prior_curation")) != (
            EXPECTED_RESPONSE_FIELD, EXPECTED_RESPONSE_SOURCE, EXPECTED_RESPONSE_SHA256
        ):
        raise ValueError("CTRP endpoint source/field/checksum contract changed")
    if config.get("external_response_access") != "forbidden_until_all_freeze_gates_pass":
        raise ValueError("Premature external-response access requested")
    if config.get("drug_ids") != EXPECTED_DRUGS or config.get("models") != EXPECTED_MODELS:
        raise ValueError("Frozen eight-drug/eight-arm panel differs from development")
    training = config.get("final_training", {})
    if (training.get("inner_folds"), training.get("inner_seed")) != (3, 2026101501):
        raise ValueError("Full-development inner split differs from the P6 decision")
    statistics = config.get("statistics", {})
    bootstrap = statistics.get("bootstrap", {})
    if (bootstrap.get("draws"), bootstrap.get("seed")) != (20000, 2026101502):
        raise ValueError("External bootstrap budget/seed differs from the P6 decision")
    if bootstrap.get("minimum_valid_fraction") != 0.99 or bootstrap.get("minimum_valid_draws") != 19800:
        raise ValueError("External bootstrap valid fraction differs from the P6 decision")
    if statistics.get("primary_comparison", {}).get("family_size") != 9:
        raise ValueError("Primary contrast family changed")
    if statistics.get("same_candidate_representation_family", {}).get("family_size") != 18:
        raise ValueError("Representation contrast family changed")
    if training.get("expected_tasks") not in (None, CONTRACT["final_model_tasks"]):
        raise ValueError("Final model task budget differs")
    if training.get("expected_fit_events") != CONTRACT["final_fit_events"]:
        raise ValueError("Final model fit budget differs")
    if training.get("expected_tasks") != CONTRACT["final_model_tasks"]:
        raise ValueError("Final model task budget differs")
    if (config.get("pre_fit_lock_path"), config.get("pre_external_execution_lock_path")) != (
        "config/external_validation_v1_1.lock.json", "config/external_execution_v1_1.lock.json"
    ):
        raise ValueError("P6 two-lock locations changed")


def _environment() -> dict:
    return {
        "python": platform.python_version(),
        "packages": {name: importlib.metadata.version(name) for name in
                     ("numpy", "pandas", "scipy", "scikit-learn", "joblib", "threadpoolctl", "psutil")},
    }


def verify_approval_status(review: str, acceptance: str) -> None:
    """Require one exact top-level status line from each review document."""
    review_states = [line.strip() for line in review.splitlines()
                     if line.strip().startswith("审查状态：")]
    acceptance_states = [line.strip() for line in acceptance.splitlines()
                         if line.strip().startswith("验收状态：")]
    if review_states != ["审查状态：可首锁"] or acceptance_states != ["验收状态：通过"]:
        raise ValueError("Independent review or pre-fit acceptance gate remains open")


def build_lock(revision: str | None = None) -> dict:
    config = json.loads((ROOT / "config/external_validation_v1_1.json").read_text(encoding="utf-8"))
    _verify_configuration(config)
    review = (ROOT / "reports/阶段12e_P6_v1_1冻结前独立审查.md").read_text(encoding="utf-8")
    acceptance = (ROOT / "reports/阶段12f_P6_v1_1首锁前验收.md").read_text(encoding="utf-8")
    verify_approval_status(review, acceptance)
    predecessor = json.loads((ROOT / PREDECESSOR_LOCK).read_text(encoding="utf-8"))
    if (sha(ROOT / PREDECESSOR_LOCK) != PREDECESSOR_LOCK_SHA256 or
            predecessor.get("status") != "frozen_pre_final_fit" or
            predecessor.get("external_response_accessed") is not False or
            predecessor.get("contract") != CONTRACT):
        raise ValueError("Preserved P6 v1 predecessor lock changed")
    verify_development_inputs("v2")
    parent = json.loads((ROOT / "config/ordinal_repeated_v1.lock.json").read_text(encoding="utf-8"))
    for category in ("files", "v1_files", "diagnostic_files"):
        for name, expected in parent[category].items():
            if sha(ROOT / name) != expected:
                raise ValueError("Historical scientific input changed: " + name)
    cohort = json.loads((ROOT / "results/cohort_v2/public_cohort_snapshot.json").read_text(encoding="utf-8"))
    if cohort["status"] != "validated_and_locked":
        raise ValueError("Cancer-model cohort v2 is not validated")
    curation = json.loads((ROOT / "data_processed/ctrp_curation_summary.json").read_text(encoding="utf-8"))
    if (curation.get("response_store"), curation.get("response_field"), curation.get("raw_curve_sha256")) != (
        EXPECTED_RESPONSE_SOURCE, EXPECTED_RESPONSE_FIELD, EXPECTED_RESPONSE_SHA256
    ):
        raise ValueError("Prior response-free curation source/field/checksum differs")
    counts = cohort["cohort_counts"]
    for key, expected in [
        ("development_unique_groups", 626),
        ("external_unique_groups", 175),
        ("development_drug_model_observations", 4967),
        ("external_drug_model_observations", 1185),
    ]:
        if counts.get(key) != expected:
            raise ValueError("Frozen cohort count changed: " + key)
    git_revision, source_hashes = _committed_snapshot(EXECUTION_FILES, revision)
    files = {**source_hashes, **{name: sha(ROOT / name) for name in DIRECT_INPUTS}}
    return {
        "schema_version": "1.0",
        "status": "frozen_pre_final_fit",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": "docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md",
        "configuration": "config/external_validation_v1_1.json",
        "predecessor_lock_sha256": PREDECESSOR_LOCK_SHA256,
        "execution_commit": git_revision,
        "execution_files": EXECUTION_FILES,
        "files": files,
        "historical_repeated_lock_sha256": sha(ROOT / "config/ordinal_repeated_v1.lock.json"),
        "cohort_lock_sha256": sha(ROOT / "config/cohort_freeze_v2.lock.json"),
        "contract": CONTRACT,
        "environment": _environment(),
        "external_response_accessed": False,
        "historical_outcomes_seen": True,
        "mutable_plan_reports_validators_excluded": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ready", action="store_true", help="Verify source readiness without writing a lock")
    group.add_argument("--lock", action="store_true", help="Write a new immutable lock after all gates pass")
    group.add_argument("--verify-only", action="store_true", help="Verify the existing lock")
    args = parser.parse_args()
    access, blocked = install_pre_external_guard()
    path = ROOT / "config/external_validation_v1_1.lock.json"
    if args.verify_only:
        previous = json.loads(path.read_text(encoding="utf-8"))
        current = build_lock(previous["execution_commit"])
        for key in ("status", "protocol", "configuration", "predecessor_lock_sha256", "execution_commit", "execution_files", "files", "contract", "environment"):
            if previous.get(key) != current[key]:
                raise ValueError("P6 lock verification failed: " + key)
        status = "verified"
    else:
        current = build_lock()
        if args.lock:
            with path.open("x", encoding="utf-8", newline="\n") as handle:
                json.dump(current, handle, ensure_ascii=False, indent=2, allow_nan=False)
                handle.write("\n")
            status = "locked"
        else:
            status = "ready"
    if blocked:
        raise ValueError("External response guard blocked an attempted access")
    print(json.dumps({
        "status": status,
        "lock_sha256": sha(path) if path.exists() else None,
        "execution_commit": current["execution_commit"],
        "contract": current["contract"],
        "external_response_accessed": False,
        "development_data_paths_accessed": len(access),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
