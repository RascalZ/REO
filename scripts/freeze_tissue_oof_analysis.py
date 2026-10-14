"""Freeze the descriptive P5e tissue-stratified OOF computation."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"
LOCK = ROOT / "config/p5e_protocol_v1.lock.json"
OUT = ROOT / "config/tissue_oof_analysis_v1.lock.json"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build() -> dict:
    protocol_lock = read_json(LOCK)
    meta = read_json(RUN / "ordinal_repeated_metadata.json")
    if protocol_lock["status"] != "frozen_before_p5e_derived_analysis_and_model_fits":
        raise ValueError("P5e protocol is not frozen")
    if sha(LOCK) != "deb15cd16258b85c86f1605a3099b85319b33a69f086c8ffe912ba12543e0e1b":
        raise ValueError("Unexpected P5e protocol lock")
    for relative, expected in protocol_lock["inputs"].items():
        if sha(ROOT / relative) != expected:
            raise ValueError("P5e protocol input changed: " + relative)
    if meta["status"] != "completed" or meta["external_response_accessed"] is not False:
        raise ValueError("Fixed eight-drug source OOF is not a valid sealed-development input")
    oof = RUN / "oof_predictions.csv"
    if sha(oof) != protocol_lock["source_oof_sha256"]:
        raise ValueError("Frozen OOF source SHA mismatch")
    inputs = [
        "scripts/freeze_tissue_oof_analysis.py",
        "src/analyze_tissue_sensitivity.py",
        "config/p5e_protocol_v1.lock.json",
        "config/tissue_sensitivity_loto_v1.json",
        "config/ordinal_repeated_v1.lock.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/ordinal_repeated_metadata.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/oof_predictions.csv",
    ]
    return {
        "schema_version": "1.0",
        "status": "frozen_before_tissue_oof_scoring",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_lock_sha256": sha(LOCK),
        "analysis_script_sha256": sha(ROOT / "src/analyze_tissue_sensitivity.py"),
        "inputs": {name: sha(ROOT / name) for name in inputs},
        "expected_contract": {
            "oof_rows": 119208, "metric_rows": 4416, "macro_rows": 192,
            "drug_ids": [1931, 1036, 1191, 1930, 1199, 1091, 1010, 1922],
            "repeats": [0, 1, 2], "models": 8, "tissues": 23,
            "display_tissues": protocol_lock["coverage_freeze"]["primary_display_tissues"],
            "new_model_fits": 0,
        },
        "scoring_rules": {
            "within_tissue_residual_rho": "Spearman(observed - tissue_null, model_prediction - tissue_null); NA for constant/undefined values",
            "rmse_gain": "reference RMSE minus candidate RMSE; paired on identical frozen OOF rows",
            "model_delta_rho": "candidate absolute within-tissue Spearman minus expression_strong on shared OOF rows",
            "repeat_handling": "score each repeat separately; never pool OOF repeats or treat as new samples",
            "uncertainty": "no new confidence intervals or p-values; descriptive only",
            "external_response_access": "forbidden",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", action="store_true")
    args = parser.parse_args()
    current = build()
    if OUT.exists():
        existing = read_json(OUT)
        for key, value in current.items():
            if key != "created_utc" and existing.get(key) != value:
                raise ValueError(f"Tissue OOF analysis lock mismatch at {key}")
        print(json.dumps({"status": "verified", "lock_sha256": sha(OUT)}))
        return
    if not args.lock:
        print(json.dumps({"status": "ready_for_lock", "expected": current["expected_contract"]}))
        return
    with OUT.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(current, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "locked", "lock_sha256": sha(OUT)}))


if __name__ == "__main__":
    main()
