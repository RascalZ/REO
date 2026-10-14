"""Freeze the protocol-specified equal-drug summary of tissue OOF metrics."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK_PATH = ROOT / "config/tissue_oof_fixed8_macro_v1.lock.json"
SOURCE_DIR = ROOT / "results/tissue_sensitivity"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build() -> dict:
    protocol = ROOT / "config/p5e_protocol_v1.lock.json"
    analysis = ROOT / "config/tissue_oof_analysis_v1.lock.json"
    macro_csv = SOURCE_DIR / "tissue_macro_by_repeat.csv"
    stage8_manifest = SOURCE_DIR / "tissue_sensitivity_manifest.json"
    protocol_obj = read_json(protocol)
    analysis_obj = read_json(analysis)
    manifest = read_json(stage8_manifest)
    if protocol_obj["external_response_accessed"] is not False or manifest["external_response_accessed"] is not False:
        raise ValueError("External response access is forbidden")
    if analysis_obj["expected_contract"]["macro_rows"] != 192:
        raise ValueError("Unexpected frozen per-drug tissue summary contract")
    if manifest["new_model_fits"] != 0 or sha(macro_csv) != manifest["output_hashes"][
        "results/tissue_sensitivity/tissue_macro_by_repeat.csv"
    ]:
        raise ValueError("Stage 8 aggregate source changed")
    inputs = [
        "scripts/freeze_tissue_oof_fixed8_macro.py",
        "src/summarize_tissue_oof_fixed8.py",
        "config/p5e_protocol_v1.lock.json",
        "config/tissue_oof_analysis_v1.lock.json",
        "results/tissue_sensitivity/oof_tissue_metrics.csv",
        "results/tissue_sensitivity/tissue_macro_by_repeat.csv",
        "results/tissue_sensitivity/tissue_sensitivity_manifest.json",
    ]
    return {
        "schema_version": "1.0", "status": "frozen_before_fixed8_macro_aggregation",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_lock_sha256": sha(protocol), "source_analysis_lock_sha256": sha(analysis),
        "inputs": {name: sha(ROOT / name) for name in inputs},
        "expected_contract": {"rows": 24, "repeats": [0, 1, 2], "models": 8, "drugs": 8, "tissues_per_drug": 10},
        "aggregation_rule": "Equal-weight the eight precomputed per-drug means after equal-weighting the fixed ten primary tissues within drug and repeat; strict metrics require all eight drugs and all ten tissues; report available denominator without replacing strict NA",
        "repeat_rule": "Keep repeat values separate; any three-repeat mean/range is descriptive and not an inferential interval",
        "external_response_access": "forbidden",
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
                raise ValueError(f"Fixed8 tissue macro lock mismatch at {key}")
        print(json.dumps({"status": "verified", "lock_sha256": sha(LOCK_PATH)}))
        return
    if not args.lock:
        print(json.dumps({"status": "ready_for_lock", "expected": current["expected_contract"]}))
        return
    with LOCK_PATH.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(current, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "locked", "lock_sha256": sha(LOCK_PATH)}))


if __name__ == "__main__":
    main()
