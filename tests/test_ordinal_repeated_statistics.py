"""Synthetic invariants for repeat-first metrics and shared paired inference."""
import copy
import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from analyze_ordinal_repeated import MACRO, analyze_predictions, contrast_values, planned_contrasts
from analyze_formal_results import generate_global_bootstrap_weights


def fixture():
    config = json.loads((ROOT / "config/ordinal_repeated_v1.json").read_text(encoding="utf-8"))
    config["statistics"]["bootstrap_draws"] = 80
    config["statistics"]["minimum_valid_fraction"] = .9
    manifest, rows = [], []
    for drug_index, drug in enumerate(config["drug_ids"]):
        n = 12 if drug_index == 0 else 11
        y = np.arange(n, dtype=float) / n
        for i in range(n):
            identity = {"canonical_model_id": f"ID{i:02d}", "biological_group_id": f"G{i:02d}", "tissue": "A" if i < 6 else "B"}
            manifest.append({"gdsc_drug_id": drug, **identity})
            for repeat in range(3):
                for mi, model in enumerate(config["models"]):
                    if model == "tissue_null":
                        prediction = .2 if i < 6 else .7
                    elif model == "expression_strong":
                        prediction = [y, y[::-1], np.roll(y, 3)][repeat][i]
                    else:
                        prediction = y[i] + .23 * np.sin(i + repeat + mi)
                    rows.append({"drug_id": drug, "model": model, "repeat": repeat,
                        "outer_fold": (i + repeat) % 5, "observed_auc": y[i], "predicted_auc": prediction, **identity})
    # A no-coverage source group still belongs to the global bootstrap universe.
    manifest.append({"gdsc_drug_id": 9999, "canonical_model_id": "ID99", "biological_group_id": "G99", "tissue": "C"})
    manifest = pd.DataFrame(manifest)
    groups, counts = generate_global_bootstrap_weights(manifest.biological_group_id, 80, 84)
    return pd.DataFrame(rows), manifest, config, groups, counts


class OrdinalRepeatedStatisticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.manifest, cls.config, cls.groups, cls.counts = fixture()
        cls.result = analyze_predictions(cls.data, cls.manifest, cls.config, cls.groups, cls.counts)

    def test_twenty_five_unique_contrasts_and_two_separate_families(self):
        contrasts = planned_contrasts(self.config)
        self.assertEqual(len(contrasts), 25)
        self.assertEqual(len({c["contrast"] for c in contrasts}), 25)
        self.assertEqual(sum(len(c["coefficients"]) == 2 for c in contrasts), 22)
        self.assertEqual(sum(len(c["coefficients"]) == 4 for c in contrasts), 3)
        self.assertEqual(sum(c["family"] == "primary_fixed9" for c in contrasts), 1)
        self.assertEqual(sum(c["family"] == "representation_fixed18" for c in contrasts), 2)
        frame = self.result["paired_effect_intervals"]
        self.assertEqual(len(frame), 9 * 25 * 4)
        self.assertEqual(frame.family_size.eq(9).sum(), 9)
        self.assertEqual(frame.family_size.eq(18).sum(), 18)
        self.assertTrue(frame[frame.metric.ne("spearman")].family_size.isna().all())

    def test_metric_is_computed_before_repeat_averaging(self):
        drug = self.config["drug_ids"][0]
        records = self.data[self.data.drug_id.eq(drug) & self.data.model.eq("expression_strong")]
        expected = np.mean([spearmanr(part.observed_auc, part.predicted_auc).statistic for _, part in records.groupby("repeat")])
        averaged = records.groupby("canonical_model_id").agg(y=("observed_auc", "first"), p=("predicted_auc", "mean"))
        wrong = spearmanr(averaged.y, averaged.p).statistic
        row = self.result["model_intervals"].query("gdsc_drug_id == @drug and model == 'expression_strong' and metric == 'spearman'").iloc[0]
        self.assertAlmostEqual(row.point, expected)
        self.assertGreater(abs(row.point - wrong), .05)

    def test_shared_global_counts_and_actual_independent_n(self):
        np.testing.assert_array_equal(self.result["bootstrap_counts"], self.counts)
        self.assertEqual(len(self.result["bootstrap_group_order"]), 13)
        macro = self.result["model_intervals"][self.result["model_intervals"].gdsc_drug_id.eq(MACRO)]
        self.assertTrue(macro.n_bootstrap_global_groups.eq(13).all())
        self.assertTrue(macro.n_unique_biological_groups.eq(12).all())
        for model in self.config["models"]:
            for metric in ("spearman", "rmse", "mae", "r2"):
                subset = self.result["model_intervals"]
                subset = subset[subset.model.eq(model) & subset.metric.eq(metric)]
                self.assertAlmostEqual(subset[subset.gdsc_drug_id.eq(MACRO)].iloc[0].point,
                                       subset[subset.gdsc_drug_id.ne(MACRO)].point.mean())

    def test_rmse_direction_and_paired_difference_of_differences(self):
        values = {"response_association_rank": np.array([.4, .7]), "response_association_expression": np.array([.1, .3]),
                  "mad_rank": np.array([.3, .2]), "mad_expression": np.array([.2, .1])}
        coefficients = {"response_association_rank": 1., "response_association_expression": -1., "mad_rank": -1., "mad_expression": 1.}
        expected = (values["response_association_rank"] - values["response_association_expression"]) - (values["mad_rank"] - values["mad_expression"])
        np.testing.assert_allclose(contrast_values(values, coefficients, "spearman"), expected)
        np.testing.assert_allclose(contrast_values(values, coefficients, "rmse"), -expected)
        effects = self.result["paired_effect_intervals"]
        for drug in [*self.config["drug_ids"], MACRO]:
            for metric in ("spearman", "rmse", "mae", "r2"):
                block = effects[effects.gdsc_drug_id.eq(drug) & effects.metric.eq(metric)].set_index("contrast")
                self.assertAlmostEqual(block.loc["interaction:pair_full-expression", "point"],
                    block.loc["interaction:pair_full-rank", "point"] + block.loc["interaction:rank-expression", "point"])

    def test_undefined_repeat_propagates_to_drug_and_macro(self):
        data = self.data.copy()
        drug = self.config["drug_ids"][0]
        data.loc[data.drug_id.eq(drug) & data.model.eq("mad_rank") & data["repeat"].eq(1), "predicted_auc"] = .4
        result = analyze_predictions(data, self.manifest, self.config, self.groups, self.counts)
        rows = result["model_intervals"]
        rows = rows[rows.gdsc_drug_id.isin([drug, MACRO]) & rows.model.eq("mad_rank") & rows.metric.eq("spearman")]
        self.assertEqual(len(rows), 2)
        self.assertTrue(rows.point.isna().all())
        self.assertTrue(rows.ci_status.eq("withheld_undefined_point").all())
        self.assertTrue(rows.conditional_95_lower.isna().all())

    def test_missing_rows_and_mismatched_groups_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            analyze_predictions(self.data.iloc[1:], self.manifest, self.config, self.groups, self.counts)
        changed = self.data.copy()
        changed.loc[0, "biological_group_id"] = "WRONG"
        with self.assertRaisesRegex(ValueError, "identity or tissue mismatch"):
            analyze_predictions(changed, self.manifest, self.config, self.groups, self.counts)


if __name__ == "__main__":
    unittest.main()
