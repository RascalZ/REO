"""Freeze the corrected tissue-null-residual contrasts before aggregation."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK_PATH = ROOT / "config/tissue_residual_increment_v1_1.lock.json"
RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"
RESULTS = ROOT / "results/tissue_sensitivity"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build() -> dict:
    config_path = ROOT / "config/tissue_oof_residual_increment_v1_1.json"
    correction_path = ROOT / "docs/P5E_TISSUE_RHO_CLARIFICATION_V1_1.md"
    p5e_lock = ROOT / "config/p5e_protocol_v1.lock.json"
    source_lock = ROOT / "config/tissue_oof_analysis_v1.lock.json"
    source_metrics = RESULTS / "oof_tissue_metrics.csv"
    source_drug_macro = RESULTS / "tissue_macro_by_repeat.csv"
    source_manifest = RESULTS / "tissue_sensitivity_manifest.json"
    p5e = read_json(p5e_lock)
    old_analysis = read_json(source_lock)
    old_manifest = read_json(source_manifest)
    config = read_json(config_path)
    if p5e["external_response_accessed"] is not False or old_manifest["external_response_accessed"] is not False:
        raise ValueError("CTRP/external access is prohibited")
    if old_manifest["new_model_fits"] != 0 or len(config["primary_contrasts"]) != 4:
        raise ValueError("Unexpected source or contrast configuration")
    if old_analysis["expected_contract"]["metric_rows"] != 4416:
        raise ValueError("Unexpected tissue OOF grid")
    if sha(source_metrics) != old_manifest["output_hashes"]["results/tissue_sensitivity/oof_tissue_metrics.csv"]:
        raise ValueError("Stage8 tissue metric input has changed")
    if sha(source_drug_macro) != old_manifest["output_hashes"]["results/tissue_sensitivity/tissue_macro_by_repeat.csv"]:
        raise ValueError("Stage8 per-drug macro input has changed")
    files = [
        "scripts/freeze_p5e_tissue_residual_increment.py",
        "src/analyze_tissue_residual_increment.py",
        "docs/P5E_TISSUE_RHO_CLARIFICATION_V1_1.md",
        "config/tissue_oof_residual_increment_v1_1.json",
        "config/p5e_protocol_v1.lock.json",
        "config/tissue_oof_analysis_v1.lock.json",
        "config/tissue_oof_fixed8_macro_v1.lock.json",
        "results/tissue_sensitivity/oof_tissue_metrics.csv",
        "results/tissue_sensitivity/tissue_macro_by_repeat.csv",
        "results/tissue_sensitivity/tissue_sensitivity_manifest.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/ordinal_repeated_metadata.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/oof_predictions.csv",
    ]
    return {
        "schema_version": "1.0", "status": "frozen_before_residual_increment_scoring",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_addendum_sha256": sha(correction_path),
        "configuration_sha256": sha(config_path),
        "analysis_script_sha256": sha(ROOT / "src/analyze_tissue_residual_increment.py"),
        "protocol_lock_sha256": sha(p5e_lock),
        "source_analysis_lock_sha256": sha(source_lock),
        "source_tissue_metrics_sha256": sha(source_metrics),
        "source_stage8_manifest_sha256": sha(source_manifest),
        "inputs": {name: sha(ROOT / name) for name in files},
        "expected_contract": {"tissue_rows": 2208, "drug_macro_rows": 96, "fixed8_macro_rows": 12,
                              "n_drugs": 8, "n_repeats": 3, "n_tissues": 23,
                              "primary_tissues": config["tissue_sets"]["primary10"]},
        "new_model_fits": 0,
        "external_response_accessed": False,
        "method_correction": "Use tissue-null residual-scale rho for molecular increment; retain prior absolute-rho deltas as diagnostics",
        "ci_or_p_value_added": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", action="store_true")
    args = parser.parse_args()
    current = build()
    if LOCK_PATH.exists():
        old = read_json(LOCK_PATH)
        for key, value in current.items():
            if key != "created_utc" and old.get(key) != value:
                raise ValueError(f"Residual increment lock mismatch at {key}")
        print(json.dumps({"status": "verified", "lock_sha256": sha(LOCK_PATH)}))
        return
    if not args.lock:
        print(json.dumps({"status": "ready_for_lock", "expected_contract": current["expected_contract"]}))
        return
    with LOCK_PATH.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(current, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "locked", "lock_sha256": sha(LOCK_PATH)}))


if __name__ == "__main__":
    main()
