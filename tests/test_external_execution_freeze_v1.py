"""Synthetic-only checks for the second P6 lock; no CTRP response access."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import freeze_external_execution_v1 as gate  # noqa: E402
from ordinal_repeated_runtime import sha  # noqa: E402


class SecondLockSyntheticTests(unittest.TestCase):
    def test_exact_64_artifact_inventory_and_fit_budget(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            final = root / "runs/final"
            (final / "models").mkdir(parents=True)
            rows = []
            for drug in gate.EXPECTED_DRUGS:
                for arm in gate.EXPECTED_MODELS:
                    relative = f"models/d{drug}_{arm}.joblib"
                    artifact = final / relative
                    artifact.write_bytes(f"synthetic:{drug}:{arm}".encode())
                    rows.append({"task": {"drug_id": drug, "model": arm},
                                 "fit_count": 1 if arm == "tissue_null" else 31,
                                 "fresh_process_reload_verified": True,
                                 "artifact_path": relative,
                                 "artifact_sha256": sha(artifact)})
            manifest = final / "model_manifest.jsonl"
            manifest.write_text("\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8")
            with patch.object(gate, "ROOT", root), patch.object(gate, "FINAL_RUN", final):
                hashes, fits = gate.model_hashes_from_manifest(manifest)
                self.assertEqual((len(hashes), fits), (64, 1744))
                rows[0]["fit_count"] = 2
                manifest.write_text("\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "fit/reload"):
                    gate.model_hashes_from_manifest(manifest)
                rows[0]["fit_count"] = 1
                manifest.write_text("\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8")
                (final / rows[0]["artifact_path"]).write_bytes(b"corrupt")
                with self.assertRaisesRegex(ValueError, "SHA mismatch"):
                    gate.model_hashes_from_manifest(manifest)

    def test_no_response_prediction_inventory_is_exact(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "prediction_manifest.csv"
            fields = ["drug_id", "model", "canonical_model_id", "status", "predicted_auc"]
            def write(extra_field=None, drop_last=False):
                with path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields + ([extra_field] if extra_field else []))
                    writer.writeheader()
                    rows_written = 0
                    for index, drug in enumerate(gate.EXPECTED_DRUGS):
                        n = 149 if index == 0 else 148
                        for arm in gate.EXPECTED_MODELS:
                            for sample in range(n):
                                if drop_last and index == 7 and sample == n - 1 and arm == gate.EXPECTED_MODELS[-1]:
                                    continue
                                row = {"drug_id": drug, "model": arm,
                                       "canonical_model_id": f"S{sample}",
                                       "status": "predicted", "predicted_auc": "0.5"}
                                if extra_field:
                                    row[extra_field] = "0.5"
                                writer.writerow(row)
                                rows_written += 1
                return rows_written
            self.assertEqual(write(), 9480)
            self.assertEqual(gate._prediction_rows(path), 9480)
            write(extra_field="area_under_curve")
            with self.assertRaisesRegex(ValueError, "Response column leaked"):
                gate._prediction_rows(path)
            write(drop_last=True)
            with self.assertRaisesRegex(ValueError, "row count"):
                gate._prediction_rows(path)


if __name__ == "__main__":
    unittest.main()
