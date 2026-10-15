"""Synthetic/config-only P6 freeze checks; never open CTRP response or index."""
from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from freeze_external_validation_v1 import (  # noqa: E402
    CONTRACT, DIRECT_INPUTS, EXECUTION_FILES, _verify_configuration,
)
import freeze_external_validation_v1 as freeze_module  # noqa: E402


class FreezeConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.configuration = json.loads((ROOT / "config/external_validation_v1.json").read_text(encoding="utf-8"))
        cls.draft = copy.deepcopy(cls.configuration)
        cls.draft["protocol_status"] = "draft_pre_final_fit"
        cls.draft["pre_fit_gates_passed"] = False

    def frozen(self):
        config = copy.deepcopy(self.configuration)
        config["protocol_status"] = "frozen_pre_final_fit"
        config["pre_fit_gates_passed"] = True
        config["external_endpoint"]["column_and_dose_contract_status"] = "verified_pre_fit_without_external_response"
        return config

    def test_draft_refuses_final_fitting(self):
        with self.assertRaisesRegex(ValueError, "remain draft"):
            _verify_configuration(self.draft)
        for field, value in (
            ("pre_fit_gates_passed", False),
            ("external_response_access", "allowed"),
        ):
            changed = self.frozen()
            changed[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                _verify_configuration(changed)
        changed = self.frozen()
        changed["external_endpoint"]["column_and_dose_contract_status"] = "pending"
        with self.assertRaisesRegex(ValueError, "unverified"):
            _verify_configuration(changed)

    def test_frozen_panel_budget_statistics_and_two_locks(self):
        config = self.frozen()
        _verify_configuration(config)
        variants = (
            ("drug_ids", list(reversed(config["drug_ids"]))),
            ("models", config["models"][:-1]),
        )
        for field, value in variants:
            changed = copy.deepcopy(config)
            changed[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                _verify_configuration(changed)
        for path, value in (
            (("final_training", "inner_seed"), 1),
            (("final_training", "expected_fit_events"), 1736),
            (("statistics", "bootstrap", "draws"), 1000),
            (("statistics", "bootstrap", "minimum_valid_fraction"), .95),
            (("statistics", "primary_comparison", "family_size"), 8),
            (("external_endpoint", "source_field"), "ic50"),
            (("pre_external_execution_lock_path",), "same_as_pre_fit.lock"),
        ):
            changed = copy.deepcopy(config)
            node = changed
            for key in path[:-1]:
                node = node[key]
            node[path[-1]] = value
            with self.subTest(path=path), self.assertRaises(ValueError):
                _verify_configuration(changed)
        self.assertEqual(CONTRACT["final_model_tasks"], 64)
        self.assertEqual(CONTRACT["final_fit_events"], 1744)

    def test_source_contract_contains_evidence_and_no_ctrp_response(self):
        for name in ("reports/阶段12_P6端点剂量无响应审计.md", "src/external_validation_v1.py",
                     "tests/test_final_model_v1.py", "tests/test_external_validation_v1.py"):
            self.assertIn(name, EXECUTION_FILES)
        joined = "\n".join(EXECUTION_FILES + DIRECT_INPUTS).lower()
        for forbidden in ("curves_post_qc.txt", "ctrp_curve_index.csv", "ctrp--drug--cellline--auc"):
            self.assertNotIn(forbidden, joined)

    def test_review_status_requires_one_exact_approval_line(self):
        approval = "审查状态：可首锁\n"
        acceptance = "验收状态：通过\n"
        freeze_module.verify_approval_status(approval, acceptance)
        blocked = "审查状态：阻断缺陷待修复\n未来才可写审查状态：可首锁。\n"
        with self.assertRaisesRegex(ValueError, "gate remains open"):
            freeze_module.verify_approval_status(blocked, acceptance)
        with self.assertRaisesRegex(ValueError, "gate remains open"):
            freeze_module.verify_approval_status(approval + approval, acceptance)
        with self.assertRaisesRegex(ValueError, "gate remains open"):
            freeze_module.verify_approval_status(approval, "验收状态：未通过\n" + acceptance)

    def test_pre_external_guard_allows_only_synthetic_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            ctrp = root / "data_raw/ctrp"
            ctrp.mkdir(parents=True)
            metadata = ctrp / "v20._COLUMNS.txt"
            response = ctrp / "v20.data.curves_post_qc.txt"
            metadata.write_text("synthetic columns", encoding="utf-8")
            response.write_text("synthetic response", encoding="utf-8")
            with patch.object(freeze_module, "ROOT", root), \
                    patch.object(freeze_module, "original_input", return_value=root / "legacy.csv"):
                accessed, blocked = freeze_module.install_pre_external_guard()
            self.assertEqual(metadata.read_text(encoding="utf-8"), "synthetic columns")
            self.assertIn(str(metadata.resolve()), accessed)
            with self.assertRaisesRegex(PermissionError, "forbidden"):
                response.read_text(encoding="utf-8")
            self.assertEqual(blocked, [str(response.resolve())])


if __name__ == "__main__":
    unittest.main()
