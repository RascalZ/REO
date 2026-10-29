"""Regression checks for the Stage 7 repeat-0-only display-count defect."""

from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from analyze_development_diagnostics import aggregate_zero_counts


def fixture():
    return [
        {"gdsc_drug_id": "1931", "model": "pair", "repeat": str(repeat),
         "outer_fits": "5", "zero_molecular_fits": str(count)}
        for repeat, count in enumerate((3, 2, 2))
    ]


class DevelopmentZeroCountAggregationTests(unittest.TestCase):
    def test_counts_all_repeats_instead_of_repeat_zero(self):
        self.assertEqual(aggregate_zero_counts(fixture()), {(1931, "pair"): 7})

    def test_missing_repeat_does_not_silently_shorten_denominator(self):
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            aggregate_zero_counts(fixture()[:2])

    def test_duplicate_repeat_is_not_double_counted(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            aggregate_zero_counts(fixture() + [fixture()[0]])

    def test_counts_cannot_exceed_saved_outer_fits(self):
        rows = fixture()
        rows[0]["zero_molecular_fits"] = "6"
        with self.assertRaisesRegex(ValueError, "Invalid"):
            aggregate_zero_counts(rows)


if __name__ == "__main__":
    unittest.main()
