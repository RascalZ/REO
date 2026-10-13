"""Tests for shared-resample OOF path decomposition, with no fitted models."""
import copy
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from analyze_oof_decomposition import analyze_decomposition, decompose_arrays
from analyze_formal_results import generate_global_bootstrap_weights, weighted_metrics


def fixture():
    drugs = [1, 2]
    rows, manifest = [], []
    for drug in drugs:
        n = 7 if drug == 1 else 6
        for i in range(n):
            group = "G" + str(i // 2)  # Multiple observations stay grouped.
            manifest.append({"gdsc_drug_id": drug, "drug_name": "D" + str(drug),
                "canonical_model_id": "ID" + str(i), "biological_group_id": group})
            for repeat in range(3):
                for k, model in enumerate(["expression", "expression_matched", "reo"]):
                    y = [0.1, 0.2, 0.2, 0.4, 0.65, 0.8, 0.85][i]
                    p = y + .13 * np.sin((i + 1) * (k + 1) + repeat) + .02 * k
                    rows.append({"drug_id": drug, "model": model, "repeat": repeat,
                        "outer_fold": i % 2, "canonical_model_id": "ID" + str(i),
                        "biological_group_id": group, "observed_auc": y, "predicted_auc": p})
    config = {"outer_repeats": 3, "drug_ids_in_frozen_rank": drugs,
        "bootstrap": {"draws": 200, "seed": 62, "computational_block_draws": 41,
                      "minimum_valid_draw_fraction": .9, "unadjusted_level": .95}}
    groups, counts = generate_global_bootstrap_weights([m["biological_group_id"] for m in manifest], 200, 62)
    return pd.DataFrame(rows), pd.DataFrame(manifest), config, groups, counts


class DecompositionTests(unittest.TestCase):
    def test_identity_every_draw_and_favorable_error_direction(self):
        rng = np.random.default_rng(42)
        a, b, c = rng.random((3, 3, 1000))
        rho, audit = decompose_arrays(a, b, c)
        rmse, _ = decompose_arrays(a, b, c, -1)
        np.testing.assert_allclose(rho["total"], rho["candidate_range"] + rho["representation_and_fitting"], atol=1e-15)
        np.testing.assert_allclose(rmse["total"], -rho["total"])
        self.assertEqual(audit["n_finite_coordinates"], 3000)

    def test_missing_metric_is_recorded_not_silently_removed(self):
        result, audit = decompose_arrays([1, np.nan], [2, 3], [4, 5])
        self.assertTrue(np.isnan(result["candidate_range"][1]))
        self.assertEqual(audit["n_invalid_coordinates_retained"], 1)

    def test_paired_interval_is_not_endpoint_subtraction(self):
        # Correlated source statistics give a constant paired change. Subtracting
        # marginal interval endpoints would incorrectly produce a broad interval.
        a = np.linspace(-1, 1, 1000)
        b, c = a + .2, a + .7
        effects, _ = decompose_arrays(a, b, c)
        ci = np.quantile(effects["candidate_range"], [.025, .975])
        marginal_subtraction = np.quantile(b, [.025, .975]) - np.quantile(a, [.975, .025])
        np.testing.assert_allclose(ci, [.2, .2])
        self.assertGreater(np.ptp(marginal_subtraction), 3)

    def test_grouped_weights_match_explicit_duplicate_spearman(self):
        y = np.array([1., 1., 3., 4., 5.])
        p = np.array([4., 3., 3., 1., 2.])
        w = np.array([[2, 2, 0, 1, 1], [0, 0, 2, 3, 2]])
        observed = weighted_metrics(y, p, w)
        for i, weights in enumerate(w):
            np.testing.assert_allclose(observed["spearman"][i], spearmanr(np.repeat(y, weights), np.repeat(p, weights)).statistic)
            np.testing.assert_allclose(observed["rmse"][i], np.sqrt(np.mean(np.repeat((y - p) ** 2, weights))))

    def test_complete_outputs_macro_and_identity(self):
        data, manifest, config, groups, counts = fixture()
        result = analyze_decomposition(data, manifest, config, groups, counts)
        rows = result["paired_decomposition_intervals"]
        self.assertEqual(len(rows), 18)
        self.assertLess(result["identity_audit"].max_abs_identity_error.max(), 1e-12)
        for metric in ["spearman", "rmse"]:
            for component in ["candidate_range", "representation_and_fitting", "total"]:
                selected = rows[rows.metric.eq(metric) & rows.component.eq(component)]
                macro = selected[selected.gdsc_drug_id.eq("macro_equal_fixed_drugs")].iloc[0]
                drugs = selected[selected.gdsc_drug_id.isin([1, 2])]
                np.testing.assert_allclose(macro.point, drugs.point.mean(), atol=1e-14)

    def test_missing_or_mismatched_id_and_unpaired_group_order_fail(self):
        data, manifest, config, groups, counts = fixture()
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            analyze_decomposition(data.iloc[1:], manifest, config, groups, counts)
        altered = data.copy()
        altered.loc[0, "biological_group_id"] = "WRONG"
        with self.assertRaisesRegex(ValueError, "group identity"):
            analyze_decomposition(altered, manifest, config, groups, counts)
        with self.assertRaisesRegex(ValueError, "group order"):
            analyze_decomposition(data, manifest, config, groups[::-1], counts)


if __name__ == "__main__":
    unittest.main()
