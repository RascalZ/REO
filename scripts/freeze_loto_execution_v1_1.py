"""Pre-fit freeze and provenance check for the LOTO execution amendment v1.1."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from project_paths import REPOSITORY as ROOT
from run_development_pilot import execution_source_snapshot, install_data_guard, sha, verify_development_inputs
from modeling.nested_cv import candidate_grid


LOCK_PATH = ROOT / "config/loto_execution_v1_1.lock.json"
CONFIG_PATH = ROOT / "config/loto_execution_v1_1.json"
MODEL_CONFIG = ROOT / "config/ordinal_repeated_v1.json"
SOURCE_RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"
SOURCE_FILES = (
    "src/run_tissue_loto_v1_1.py",
    "src/validate_tissue_loto_v1_1.py",
    "src/analyze_tissue_loto_v1_1.py",
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
LOCKED_FILES = (
    "docs/LOTO_EXECUTION_AMENDMENT_V1_1.md",
    "docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md",
    "config/loto_execution_v1_1.json",
    "config/tissue_sensitivity_loto_v1.json",
    "config/p5e_protocol_v1.lock.json",
    "config/ordinal_repeated_v1.json",
    "config/ordinal_repeated_v1.lock.json",
    "config/synthetic_mechanism_execution_v1_1.lock.json",
    "config/synthetic_mechanisms_v1_1.json",
    "docs/SYNTHETIC_MECHANISM_EXECUTION_V1_1.md",
    "scripts/freeze_loto_execution_v1_1.py",
    "scripts/audit_synthetic_selection_consistency.py",
    "src/run_tissue_loto_v1_1.py",
    "src/validate_tissue_loto_v1_1.py",
    "src/analyze_tissue_loto_v1_1.py",
    "reports/阶段9_合成机制控制.md",
    "reports/阶段9b_候选基因集一致性审计.md",
    "results/synthetic_mechanisms/synthetic_mechanism_manifest.json",
    "results/synthetic_mechanisms/representation_gene_set_audit.csv",
)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build_lock() -> dict:
    config = read_json(CONFIG_PATH)
    plan = read_json(ROOT / "config/tissue_sensitivity_loto_v1.json")
    model_config = read_json(MODEL_CONFIG)
    p5e_lock = read_json(ROOT / "config/p5e_protocol_v1.lock.json")
    repeated_lock = read_json(ROOT / "config/ordinal_repeated_v1.lock.json")

    if config["inner_seed_resolution"]["seed_base"] != plan["inner_seed_base"]:
        raise ValueError("LOTO amendment must use the existing machine-configured seed base")
    if config["tissue_order"] != sorted(plan["tissues"]):
        raise ValueError("LOTO amendment tissue order differs from the frozen tissue inventory")
    if config["drug_ids"] != plan["drug_ids"] or config["models"] != plan["models"]:
        raise ValueError("LOTO amendment drug/model panel differs from P5e v1")
    expected_seeds = {tissue: int(plan["inner_seed_base"]) + i for i, tissue in enumerate(config["tissue_order"])}
    if config["inner_seed_resolution"]["seed_by_tissue"] != expected_seeds:
        raise ValueError("LOTO tissue seed map does not follow the resolved deterministic rule")
    if config["expected_contract"] != plan["expected_contract"]:
        raise ValueError("LOTO execution contract differs from the frozen P5e contract")
    if len(candidate_grid("expression", model_config)) != 10 or int(model_config["inner_folds"]) != 3:
        raise ValueError("Fixed-eight inner grid/fold contract changed")

    original_protocol = (ROOT / "docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md").read_text(encoding="utf-8")
    if "2026101801 + sorted_tissue_index" not in original_protocol:
        raise ValueError("Expected original prose seed declaration was not found; review amendment applicability")
    if int(plan["inner_seed_base"]) != 2026101401:
        raise ValueError("Unexpected source machine-configured LOTO seed base")
    if p5e_lock["inputs"].get("docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md") != sha(ROOT / "docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md"):
        raise ValueError("Original P5e protocol no longer matches its frozen lock")
    if p5e_lock["inputs"].get("config/tissue_sensitivity_loto_v1.json") != sha(ROOT / "config/tissue_sensitivity_loto_v1.json"):
        raise ValueError("Original LOTO config no longer matches its frozen lock")
    if repeated_lock.get("files", {}).get("config/ordinal_repeated_v1.json") != sha(MODEL_CONFIG):
        raise ValueError("Fixed-eight source model configuration changed")
    if p5e_lock["inputs"].get("config/ordinal_repeated_v1.lock.json") != sha(ROOT / "config/ordinal_repeated_v1.lock.json"):
        raise ValueError("Fixed-eight source model lock changed")

    synthetic = read_json(ROOT / "results/synthetic_mechanisms/synthetic_mechanism_manifest.json")
    if synthetic["status"] != "completed" or synthetic["counts"] != {
        "tasks": 960, "fits": 26160, "oof_predictions": 76800,
        "cohort_summaries": 320, "scenario_model_summaries": 32,
        "all_checkpoints_hash_verified": True,
    }:
        raise ValueError("Synthetic mechanism run did not satisfy its completed contract")
    if synthetic["execution_lock_sha256"] != sha(ROOT / "config/synthetic_mechanism_execution_v1_1.lock.json"):
        raise ValueError("Synthetic predecessor execution lock changed")
    if synthetic["external_response_accessed"] or synthetic["ctrp_accessed"]:
        raise ValueError("Synthetic control accessed external responses")
    for name, expected in synthetic["outputs"].items():
        if sha(ROOT / name) != expected:
            raise ValueError("Synthetic predecessor output hash mismatch: " + name)
    with (ROOT / "results/synthetic_mechanisms/representation_gene_set_audit.csv").open(encoding="utf-8-sig", newline="") as handle:
        audit_rows = list(csv.DictReader(handle))
    if len(audit_rows) != 240 or not all(row["candidate_gene_sets_identical"] == "True" for row in audit_rows):
        raise ValueError("Synthetic candidate-set audit is incomplete or has mismatched representations")

    access, blocked = install_data_guard()
    development_verification = verify_development_inputs("v2")
    if blocked:
        raise ValueError("Blocked external-data access was attempted during LOTO pre-freeze")
    source_meta = read_json(SOURCE_RUN / "ordinal_repeated_metadata.json")
    if source_meta.get("status") != "completed" or source_meta.get("external_response_accessed") is not False:
        raise ValueError("Fixed-eight development source run is not complete or violates the data gate")
    source_response = SOURCE_RUN / "development_responses.csv"
    if sha(source_response) != p5e_lock["inputs"].get(
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/development_responses.csv"
    ):
        raise ValueError("Development response source differs from the P5e frozen input")

    execution_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    os.environ["REO_EXECUTION_COMMIT"] = execution_commit
    revision, source_hashes = execution_source_snapshot([ROOT / name for name in SOURCE_FILES])
    if revision != execution_commit:
        raise ValueError("Execution source snapshot is not the current committed revision")

    locked_hashes = {name: sha(ROOT / name) for name in LOCKED_FILES}
    prerequisite_hashes = {
        "source_run_metadata_sha256": sha(SOURCE_RUN / "ordinal_repeated_metadata.json"),
        "source_development_responses_sha256": sha(source_response),
        "source_oof_sha256": sha(SOURCE_RUN / "oof_predictions.csv"),
        "source_development_manifest_sha256": sha(ROOT / "results/cohort_v2/development_sample_manifest.csv"),
        "development_input_freeze_sha256": development_verification["freeze_sha256"],
        "development_freeze_lock_sha256": sha(ROOT / "config/cohort_freeze_v2.lock.json"),
    }
    return {
        "schema_version": "1.0",
        "status": "frozen_before_loto_model_fits",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_version": config["protocol_version"],
        "execution_commit": execution_commit,
        "source_snapshot": {"git_commit": revision, "source_hashes": source_hashes},
        "configuration_sha256": sha(CONFIG_PATH),
        "model_config_sha256": sha(MODEL_CONFIG),
        "p5e_lock_sha256": sha(ROOT / config["source_p5e_lock"]),
        "source_model_lock_sha256": sha(ROOT / config["source_model_lock"]),
        "locked_file_hashes": locked_hashes,
        "prerequisite_hashes": prerequisite_hashes,
        "development_input_verification": development_verification,
        "seed_resolution": {
            "original_prose_base": 2026101801,
            "frozen_machine_config_base": 2026101401,
            "chosen_base": config["inner_seed_resolution"]["seed_base"],
            "reason": "Follow the frozen executable config field; preserve the original protocol/lock and document this pre-fit clarification in v1.1.",
            "seed_by_tissue": config["inner_seed_resolution"]["seed_by_tissue"],
        },
        "synthetic_predecessor": {
            "manifest_sha256": sha(ROOT / "results/synthetic_mechanisms/synthetic_mechanism_manifest.json"),
            "candidate_audit_sha256": sha(ROOT / "results/synthetic_mechanisms/representation_gene_set_audit.csv"),
            "tasks": 960, "fits": 26160, "oof_predictions": 76800,
        },
        "expected_contract": config["expected_contract"],
        "external_response_accessed": False,
        "ctrp_accessed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", action="store_true", help="Create the immutable pre-fit LOTO lock")
    args = parser.parse_args()
    current = build_lock()
    if LOCK_PATH.exists():
        previous = read_json(LOCK_PATH)
        for key, value in current.items():
            if key != "created_utc" and previous.get(key) != value:
                raise ValueError(f"Existing LOTO execution lock differs at {key}")
        print(json.dumps({"status": "verified", "lock_sha256": sha(LOCK_PATH)}, ensure_ascii=False))
        return
    if not args.lock:
        print(json.dumps({"status": "ready_to_lock", "execution_commit": current["execution_commit"],
                          "seed_by_tissue": current["seed_resolution"]["seed_by_tissue"],
                          "expected_contract": current["expected_contract"]}, ensure_ascii=False))
        return
    LOCK_PATH.write_text(json.dumps(current, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "locked", "lock_sha256": sha(LOCK_PATH)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
