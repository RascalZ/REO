"""Freeze the status-aligned v1.3 runner for null-residual tissue contrasts."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK_PATH = ROOT / "config/tissue_residual_increment_v1_3.lock.json"
RESULTS = ROOT / "results/tissue_sensitivity"
RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def build() -> dict:
    config = ROOT / "config/tissue_oof_residual_increment_v1_1.json"
    correction = ROOT / "docs/P5E_TISSUE_RHO_CLARIFICATION_V1_1.md"
    prior_lock = ROOT / "config/tissue_residual_increment_v1_1.lock.json"
    prior_runner_lock = ROOT / "config/tissue_residual_increment_v1_2.lock.json"
    prior_lock_obj = read_json(prior_lock)
    prior_runner_lock_obj = read_json(prior_runner_lock)
    stage8_manifest = read_json(RESULTS / "tissue_sensitivity_manifest.json")
    p5e = read_json(ROOT / "config/p5e_protocol_v1.lock.json")
    if sha(prior_lock) != "bc6c591216b5d13218b93cbf0f149ae9c8d7f483341371d97b5e8c86eff47dd8":
        raise ValueError("Original v1.1 method lock changed")
    if sha(prior_runner_lock) != "8514a695a02d4e9d835d3cbfdd8e7003964da9d97fb35109e8f280b6c6bc6ac5":
        raise ValueError("v1.2 runner integration lock changed")
    if prior_lock_obj["external_response_accessed"] or p5e["external_response_accessed"]:
        raise ValueError("CTRP/external response access is prohibited")
    if stage8_manifest["new_model_fits"] != 0:
        raise ValueError("This summary must not include new model fits")
    for relative, expected in p5e["inputs"].items():
        if sha(ROOT / relative) != expected:
            raise ValueError("Frozen P5e input changed: " + relative)
    files = [
        "scripts/freeze_p5e_tissue_residual_increment_v1_3.py",
        "src/analyze_tissue_residual_increment_v1_3.py",
        "scripts/freeze_p5e_tissue_residual_increment_v1_2.py",
        "src/analyze_tissue_residual_increment_v1_2.py",
        "src/analyze_tissue_residual_increment.py",
        "docs/P5E_TISSUE_RHO_CLARIFICATION_V1_1.md",
        "config/tissue_oof_residual_increment_v1_1.json",
        "config/tissue_residual_increment_v1_1.lock.json",
        "config/tissue_residual_increment_v1_2.lock.json",
        "config/p5e_protocol_v1.lock.json",
        "config/tissue_oof_analysis_v1.lock.json",
        "results/tissue_sensitivity/oof_tissue_metrics.csv",
        "results/tissue_sensitivity/tissue_macro_by_repeat.csv",
        "results/tissue_sensitivity/tissue_sensitivity_manifest.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/ordinal_repeated_metadata.json",
        "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/oof_predictions.csv",
    ]
    return {
        "schema_version": "1.0", "status": "frozen_before_residual_increment_scoring_v1_3",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parent_method_lock_sha256": sha(prior_lock),
        "parent_runner_lock_sha256": sha(prior_runner_lock),
        "protocol_lock_sha256": sha(ROOT / "config/p5e_protocol_v1.lock.json"),
        "configuration_sha256": sha(config), "protocol_addendum_sha256": sha(correction),
        "source_tissue_metrics_sha256": sha(RESULTS / "oof_tissue_metrics.csv"),
        "inputs": {name: sha(ROOT / name) for name in files},
        "expected_contract": {"tissue_rows": 2208, "drug_macro_rows": 96, "fixed8_macro_rows": 12,
                              "n_drugs": 8, "n_repeats": 3, "n_tissues": 23,
                              "primary_tissues": read_json(config)["tissue_sets"]["primary10"]},
        "new_model_fits": 0, "external_response_accessed": False,
        "correction_scope": "Report field and available-cell denominator correction plus frozen-status binding; no scientific metric/formula changes",
        "preserved_implementation_attempts": [
            {"version": "v1.1", "failure": "Report generation stopped after aggregate CSV writes due to nonexistent report field",
             "model_fits": 0, "external_access": False},
            {"version": "v1.2", "failure": "Stopped at pre-score lock-status check because the runner and lock used different status strings",
             "model_fits": 0, "external_access": False},
        ],
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
                raise ValueError(f"v1.3 residual increment lock mismatch at {key}")
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
