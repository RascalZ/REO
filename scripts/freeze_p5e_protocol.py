"""Freeze the P5e synthetic, stratified-OOF, and LOTO protocol inputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"
SOURCE_LOCK = ROOT / "config/ordinal_repeated_v1.lock.json"
LOCK_PATH = ROOT / "config/p5e_protocol_v1.lock.json"
DRUG_IDS = [1931, 1036, 1191, 1930, 1199, 1091, 1010, 1922]
MODELS = [
    "tissue_null", "expression_strong", "mad_expression", "mad_rank", "mad_pair_full",
    "response_association_expression", "response_association_rank", "response_association_pair_full",
]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def tissue_coverage(response_csv: Path) -> dict:
    counts: dict[tuple[str, str], set[str]] = {}
    with response_csv.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            drug = row["gdsc_drug_id"]
            tissue = row["tissue"]
            counts.setdefault((drug, tissue), set()).add(row["biological_group_id"])
    by_drug = {str(drug): {} for drug in DRUG_IDS}
    for (drug, tissue), groups in counts.items():
        if int(drug) in DRUG_IDS:
            by_drug[drug][tissue] = len(groups)
    common_tissues = sorted(
        tissue for tissue in set.intersection(*(set(by_drug[str(d)]) for d in DRUG_IDS))
        if min(by_drug[str(d)][tissue] for d in DRUG_IDS) >= 20
    )
    return {"drug_tissue_group_counts": by_drug, "primary_display_tissues": common_tissues}


def build() -> dict:
    protocol_path = ROOT / "docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md"
    loto_path = ROOT / "config/tissue_sensitivity_loto_v1.json"
    synthetic_path = ROOT / "config/synthetic_mechanisms_v1.json"
    runner_path = ROOT / "scripts/freeze_p5e_protocol.py"
    protocol = protocol_path.read_text(encoding="utf-8")
    loto = load_json(loto_path)
    synthetic = load_json(synthetic_path)
    source_meta_path = SOURCE_RUN / "ordinal_repeated_metadata.json"
    source_oof = SOURCE_RUN / "oof_predictions.csv"
    source_response = SOURCE_RUN / "development_responses.csv"
    source_meta = load_json(source_meta_path)

    assert "此方案先于本阶段新增分层统计和任何模拟/LOTO模型拟合冻结" in protocol
    assert loto["drug_ids"] == DRUG_IDS and loto["models"] == MODELS
    assert len(loto["tissues"]) == 23 and loto["outer_folds"] == 23 and loto["inner_folds"] == 3
    assert loto["expected_contract"] == {
        "tasks": 1472, "fits": 40112, "oof_predictions": 39736,
        "unknown_tissue_samples_equal_test_samples": True,
        "training_includes_heldout_tissue": False,
    }
    assert synthetic["scenario_order"] == [
        "tissue_only", "weak_molecular", "common_sample_amplitude", "predefined_pair_order"
    ]
    assert synthetic["expected_contract"] == {
        "cohorts": 40, "tasks": 960, "fits": 26160,
        "per_sample_model_oof_predictions": 76800,
        "new_real_development_model_fits": 0, "synthetic_model_fits": 26160,
    }
    assert source_meta["status"] == "completed"
    assert source_meta["external_response_accessed"] is False
    assert source_meta["completed_tasks"] == 960
    assert source_meta["expected_fit_events"] == 26160
    assert source_meta["expected_oof_predictions"] == 119208
    assert sha(source_oof) == source_meta["output_hashes"]["oof_predictions.csv"]
    assert sha(SOURCE_LOCK) == "ef8f588b99c42cb54db73b495127bc41882fa3e5b450e62167c6b46c53d583e7"

    coverage = tissue_coverage(source_response)
    assert coverage["primary_display_tissues"] == loto["expected_primary_display_tissues"]
    assert len(coverage["drug_tissue_group_counts"]) == 8

    inputs = [
        "docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md",
        "config/tissue_sensitivity_loto_v1.json",
        "config/synthetic_mechanisms_v1.json",
        "scripts/freeze_p5e_protocol.py",
        "config/ordinal_repeated_v1.json",
        "config/ordinal_repeated_v1.lock.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/ordinal_repeated_metadata.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/oof_predictions.csv",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/development_responses.csv",
        "results/cohort_v2/development_sample_manifest.csv",
    ]
    return {
        "schema_version": "1.0",
        "status": "frozen_before_p5e_derived_analysis_and_model_fits",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_version": "P5e-v1",
        "source_run_execution_commit": source_meta["execution_commit"],
        "inputs": {name: sha(ROOT / name) for name in inputs},
        "source_run": source_meta["run_name"],
        "source_run_lock_sha256": sha(SOURCE_LOCK),
        "source_run_freeze_sha256": source_meta["freeze_sha256"],
        "source_oof_sha256": sha(source_oof),
        "external_response_accessed": False,
        "external_response_access": "forbidden until internal analysis/model export gates are complete",
        "coverage_freeze": coverage,
        "loto_contract": loto["expected_contract"],
        "synthetic_contract": synthetic["expected_contract"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", action="store_true", help="Create the immutable protocol lock")
    args = parser.parse_args()
    current = build()
    if LOCK_PATH.exists():
        previous = load_json(LOCK_PATH)
        for key, value in current.items():
            if key != "created_utc" and previous.get(key) != value:
                raise ValueError(f"Frozen P5e lock mismatch at {key}")
        print(json.dumps({"status": "verified", "lock_sha256": sha(LOCK_PATH)}, ensure_ascii=False))
        return
    if not args.lock:
        print(json.dumps({"status": "ready_for_lock", "primary_display_tissues": coverage_summary()}, ensure_ascii=False))
        return
    with LOCK_PATH.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(current, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "locked", "lock_sha256": sha(LOCK_PATH)}, ensure_ascii=False))


def coverage_summary() -> list[str]:
    return tissue_coverage(SOURCE_RUN / "development_responses.csv")["primary_display_tissues"]


if __name__ == "__main__":
    main()
