"""P6 v1.1 guard-amendment checks; synthetic paths and old-lock bytes only."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import freeze_external_validation_v1_1 as amended  # noqa: E402
import freeze_external_execution_v1 as execution  # noqa: E402
import run_external_validation_v1 as external  # noqa: E402
import run_final_model_v1 as final  # noqa: E402
import validate_final_model_v1 as validator  # noqa: E402
from ordinal_repeated_runtime import sha  # noqa: E402


class GuardAmendmentTests(unittest.TestCase):
    def test_predecessor_lock_is_preserved_and_all_v1_1_paths_agree(self):
        old = ROOT / amended.PREDECESSOR_LOCK
        self.assertEqual(sha(old), amended.PREDECESSOR_LOCK_SHA256)
        self.assertFalse(json.loads(old.read_text(encoding="utf-8"))["external_response_accessed"])
        config = json.loads((ROOT / "config/external_validation_v1_1.json").read_text(encoding="utf-8"))
        amended._verify_configuration(config)
        self.assertEqual(config["pre_fit_lock_path"], "config/external_validation_v1_1.lock.json")
        self.assertEqual(config["pre_external_execution_lock_path"], "config/external_execution_v1_1.lock.json")
        self.assertEqual(final.DEFAULT_LOCK, external.FIRST_LOCK)
        self.assertEqual(final.DEFAULT_LOCK, validator.DEFAULT_LOCK)
        self.assertEqual(final.DEFAULT_LOCK, execution.FIRST_LOCK)
        self.assertEqual(external.SECOND_LOCK, execution.SECOND_LOCK)
        self.assertEqual(final.DEFAULT_EXTERNAL_CONFIG, external.EXTERNAL_CONFIG)
        self.assertEqual(final.DEFAULT_EXTERNAL_CONFIG, validator.DEFAULT_EXTERNAL_CONFIG)
        self.assertEqual(final.DEFAULT_RUN, external.FINAL_RUN)
        self.assertEqual(final.DEFAULT_RUN, execution.FINAL_RUN)
        self.assertEqual(external.EXTERNAL_RUN, execution.EXTERNAL_RUN)
        for name in ("config/external_validation_v1_1.json",
                     "docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md",
                     "reports/阶段12d_P6_v1首锁后零拟合启动故障.md",
                     "reports/阶段12e_P6_v1_1冻结前独立审查.md",
                     "reports/阶段12f_P6_v1_1首锁前验收.md",
                     "scripts/freeze_external_validation_v1_1.py"):
            self.assertIn(name, amended.EXECUTION_FILES)
        self.assertIn(amended.PREDECESSOR_LOCK, amended.DIRECT_INPUTS)
        self.assertNotIn("data_raw/ctrp/v20.data.curves_post_qc.txt",
                         amended.EXECUTION_FILES + amended.DIRECT_INPUTS)

    def test_new_guard_allows_only_pure_ctrp_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            raw = root / "data_raw/ctrp"
            processed = root / "data_processed"
            raw.mkdir(parents=True)
            processed.mkdir()
            allowed = [raw / name for name in
                       ("v20._COLUMNS.txt", "v20._README.txt", "v20.meta.per_compound.txt")]
            forbidden = [raw / "v20.data.curves_post_qc.txt",
                         processed / "ctrp_curve_index.csv",
                         root / "legacy_ctrp_response.csv"]
            for path in allowed + forbidden:
                path.write_text("synthetic", encoding="utf-8")
            with patch.object(amended, "ROOT", root), \
                    patch.object(amended, "original_input", return_value=forbidden[-1]):
                accessed, blocked = amended.install_pre_external_guard()
            for path in allowed:
                self.assertEqual(path.read_text(encoding="utf-8"), "synthetic")
                self.assertIn(str(path.resolve()), accessed)
            for path in forbidden:
                with self.assertRaises(PermissionError):
                    path.read_text(encoding="utf-8")
                self.assertIn(str(path.resolve()), blocked)


if __name__ == "__main__":
    unittest.main()
