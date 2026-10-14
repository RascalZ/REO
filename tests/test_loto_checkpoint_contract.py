"""No-fit regression tests for LOTO task identity in checkpoints and OOF rows."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import run_formal_evaluation as checkpoint_validator
import run_tissue_loto_v1_3 as loto


class LotoCheckpointIdentity(unittest.TestCase):
    def setUp(self):
        self.task = {
            "drug_id": 1931,
            "repeat": 0,
            "outer_fold": 0,
            "tissue_index": 0,
            "heldout_tissue": "Biliary Tract",
            "inner_seed": 2026101401,
            "model": "tissue_null",
            "permutation_id": 0,
            "selection_rule": "none",
            "representation": "tissue_null",
            "n_genes": 0,
        }
        self.row = {
            **self.task,
            "drug_name": "synthetic",
            "canonical_model_id": "model-001",
            "biological_group_id": "group-001",
            "tissue": "Biliary Tract",
            "observed_auc": 0.5,
            "predicted_auc": 0.5,
        }
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name) / "checkpoint.json"

    def tearDown(self):
        self.temp.cleanup()

    def write_checkpoint(self, row):
        checkpoint_validator.atomic_json(self.path, {
            "schema_version": "1.0",
            "task": self.task,
            "freeze_sha256": "frozen-test-lock",
            "status": "completed",
            "oof_rows": [row],
            "expected_test_ids": ["model-001"],
        })

    def test_inner_seed_is_required_and_preserved_in_oof_contract(self):
        self.assertIn("inner_seed", loto.OOF_FIELDS)
        row = dict(self.row)
        row.pop("inner_seed")
        self.write_checkpoint(row)
        with self.assertRaises(KeyError):
            checkpoint_validator.validate_resume_checkpoint(
                self.path, self.task, "frozen-test-lock"
            )

        self.write_checkpoint(self.row)
        self.assertTrue(checkpoint_validator.validate_resume_checkpoint(
            self.path, self.task, "frozen-test-lock"
        ))

    def test_wrong_inner_seed_is_rejected(self):
        row = dict(self.row, inner_seed=2026101801)
        self.write_checkpoint(row)
        with self.assertRaisesRegex(ValueError, "Checkpoint row context invalid"):
            checkpoint_validator.validate_resume_checkpoint(
                self.path, self.task, "frozen-test-lock"
            )


if __name__ == "__main__":
    unittest.main()
