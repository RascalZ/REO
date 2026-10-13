"""Synthetic tests of frozen ordinal pilot statistics (no model fitting)."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from analyze_ordinal_diagnostic import MACRO, aggregate_audits, analyze_predictions, planned_contrasts
from analyze_formal_results import generate_global_bootstrap_weights


def fixture():
    config = json.loads((ROOT / "config/ordinal_diagnostic_v1.json").read_text(encoding="utf-8"))
    config["statistics"]["bootstrap_draws"] = 100
    config["statistics"]["minimum_valid_fraction"] = .9
    config["permutation"]["count"] = 2
    data, manifest = [], []
    for drug in config["drug_ids"]:
        for i in range(8):
            identity = {"canonical_model_id": f"ID{i}", "biological_group_id": f"G{i}", "tissue": "A" if i < 4 else "B"}
            manifest.append({"gdsc_drug_id": drug, **identity})
            for perm in range(3):
                y = .1 + .1 * (i if perm == 0 else (i + perm) % 8)
                for mi, model in enumerate(config["models"]):
                    # identical observed pairs are an intentional zero-CI control
                    prediction = (.2 if i < 4 else .6) if mi == 0 else y + .13 * np.sin(i + mi + perm)
                    if model == "mad_rank":
                        prediction = y + .13 * np.sin(i + 1 + perm)
                    data.append({"drug_id": drug, "model": model, "permutation_id": perm, "repeat": 0,
                        "outer_fold": i % 3, "observed_auc": y, "predicted_auc": prediction, **identity})
    # Extra frozen groups without coverage in these diagnostic cases remain in
    # the global resampling universe, as in the 626-group actual cohort.
    manifest.append({"gdsc_drug_id": 9999, "canonical_model_id": "ID9", "biological_group_id": "G9", "tissue": "C"})
    manifest = pd.DataFrame(manifest)
    groups, counts = generate_global_bootstrap_weights(manifest.biological_group_id, 100, 1)
    return pd.DataFrame(data), manifest, config, groups, counts


class OrdinalStatisticsTests(unittest.TestCase):
    def test_fifteen_contrasts_directions_and_coverage(self):
        data, manifest, config, groups, counts = fixture()
        self.assertEqual(len(planned_contrasts(config)), 15)
        result = analyze_predictions(data, manifest, config, groups, counts)
        self.assertEqual(len(result["point_metrics"]), 5 * 3 * 7)
        self.assertEqual(len(result["paired_effect_intervals"]), 5 * 15 * 4)
        self.assertEqual(len(result["null_contrasts"]), 5 * 2 * 15 * 4)
        self.assertTrue(result["null_contrasts"].permutation_id.isin([1, 2]).all())
        effects = result["paired_effect_intervals"]
        identical = effects[effects.contrast.eq("mad:rank-expression")]
        np.testing.assert_array_equal(identical.point, 0.)
        np.testing.assert_array_equal(identical.conditional_95_lower, 0.)
        np.testing.assert_array_equal(identical.conditional_95_upper, 0.)
        points = result["point_metrics"]
        for metric in ("spearman", "rmse", "mae", "r2"):
            row = effects[effects.gdsc_drug_id.eq(MACRO) & effects.contrast.eq("mad:pair_full-expression") & effects.metric.eq(metric)].iloc[0]
            selected = points[points.gdsc_drug_id.eq(MACRO) & points.permutation_id.eq(0)].set_index("model")
            difference = selected.loc["mad_pair_full", metric] - selected.loc["mad_expression", metric]
            self.assertAlmostEqual(row.point, -difference if metric in ("rmse", "mae") else difference)

    def test_shared_global_weights_kept_even_for_noncase_groups(self):
        data, manifest, config, groups, counts = fixture()
        result = analyze_predictions(data, manifest, config, groups, counts)
        self.assertEqual(len(result["bootstrap_group_order"]), 9)
        np.testing.assert_array_equal(result["bootstrap_counts"], counts)
        point = result["point_metrics"]
        self.assertTrue(point[point.gdsc_drug_id.eq(MACRO)].n_unique_biological_groups.eq(8).all())
        self.assertTrue(point.n_bootstrap_global_groups.eq(9).all())
        for model in config["models"]:
            values = point[point.permutation_id.eq(0) & point.model.eq(model)]
            self.assertAlmostEqual(values[values.gdsc_drug_id.eq(MACRO)].iloc[0].rmse,
                                   values[values.gdsc_drug_id.ne(MACRO)].rmse.mean())

    def test_invalid_null_metric_retained_without_smaller_run_summary(self):
        data, manifest, config, groups, counts = fixture()
        data.loc[data.permutation_id.eq(1) & data.model.eq("mad_expression"), "predicted_auc"] = .5
        result = analyze_predictions(data, manifest, config, groups, counts)
        summary = result["null_contrast_summary"]
        row = summary[summary.gdsc_drug_id.eq(MACRO) & summary.contrast.eq("mad_expression-tissue_null") & summary.metric.eq("spearman")].iloc[0]
        self.assertEqual(row.null_runs_recorded, 2)
        self.assertEqual(row.null_runs_valid, 1)
        self.assertTrue(np.isnan(row.null_median))

    def test_missing_or_changed_assigned_outcomes_and_folds_fail(self):
        data, manifest, config, groups, counts = fixture()
        with self.assertRaisesRegex(ValueError, "Incomplete OOF"):
            analyze_predictions(data.iloc[1:], manifest, config, groups, counts)
        changed = data.copy()
        changed.loc[0, "observed_auc"] += .1
        with self.assertRaisesRegex(ValueError, "Assigned outcomes"):
            analyze_predictions(changed, manifest, config, groups, counts)
        changed = data.copy()
        changed.loc[0, "outer_fold"] = 2
        with self.assertRaisesRegex(ValueError, "Outer role"):
            analyze_predictions(changed, manifest, config, groups, counts)

    def test_empty_active_sets_are_not_perfect_stability(self):
        _, _, config, _, _ = fixture()
        records = []
        for drug in config["drug_ids"]:
            for perm in range(3):
                for model in config["models"]:
                    rule = "none" if model == "tissue_null" else "response_association" if model.startswith("response_association") else "mad"
                    for fold in range(3):
                        records.append({"task": {"drug_id": drug, "permutation_id": perm, "model": model, "outer_fold": fold},
                            "chosen_params": {}, "best_inner_rmse": .1,
                            "outer_fit_audit": {"n_nonzero_molecular_coeffs": None if model == "tissue_null" else 0,
                                "active_molecular_features": [], "active_molecular_coefficients_standardized": [],
                                "selected_gene_ids": ["A", "B"], "selected_gene_ids_sha256": "same", "selection_rule": rule,
                                "fit_seconds": 1., "last_predict_seconds": .1},
                            "resources": {"wall_seconds": 2., "sampled_peak_rss_bytes": 1000}})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "audit.jsonl"
            path.write_text("\n".join(json.dumps(row) for row in records), encoding="utf-8")
            result = aggregate_audits(path, config)
            self.assertTrue(result["active_feature_stability"].active_feature_jaccard.isna().all())
            self.assertTrue(result["active_feature_stability"].empty_active_union.all())
            records[4]["outer_fit_audit"]["selected_gene_ids_sha256"] = "different"
            path.write_text("\n".join(json.dumps(row) for row in records), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Candidate genes differ"):
                aggregate_audits(path, config)


if __name__ == "__main__":
    unittest.main()
