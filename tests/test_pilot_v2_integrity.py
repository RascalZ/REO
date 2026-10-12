"""Synthetic acceptance of inner-fit EN paths and supervised REO/RF selection.

No real drug-response or expression data is read. The fixtures deliberately
change held-out information while fixing legitimate training observations.
"""
from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
from sklearn.ensemble import RandomForestRegressor

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from modeling.features import FoldFeatureBuilder
from modeling.nested_cv import TrainOnlyRegressor, evaluate_outer_fold, make_grouped_splits
from modeling.reo_rf import TrainOnlyReoPairSelector


def residualize_within_tissue(values, tissue):
    """Independent fixture calculation, rather than using model helper code."""
    result = np.asarray(values, dtype=np.float64).copy()
    for category in np.unique(tissue):
        rows = tissue == category
        result[rows] -= result[rows].mean(axis=0)
    return result


class RelativePathTrainingBoundaryTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(77123)
        self.groups = np.repeat(np.arange(24), 2).astype(str)
        self.tissue = np.repeat(np.arange(24) % 3, 2)
        self.ids = tuple(f"G{i}" for i in range(6))
        genes = rng.lognormal(1.3, 1.0, size=(len(self.groups), len(self.ids))).astype(np.float64)
        genes[:, 0] = 100.0 + np.arange(len(self.groups))
        self.x = np.column_stack((genes, self.tissue, np.arange(len(self.groups)))).astype(np.float64)
        self.y = 0.3 * np.log2(genes[:, 1] + 1.0) + 0.5 * self.tissue + rng.normal(0, 0.07, size=len(self.groups))
        self.config = {
            "feature": {"reo_genes": 6, "expression_genes": 6, "min_tpm": 1.0,
                        "min_detection_fraction": 0.1, "max_pair_tie_fraction": 0.9},
            "elastic_net": {"regularization_mode": "alpha_max_relative", "alpha_ratios": [1.0, 0.1],
                            "l1_ratios": [0.5, 1.0], "max_iter": 2000, "tol": 1e-5},
            "ridge_alphas": [0.1, 1.0],
            "random_forest": {"n_estimators": 20, "max_features": ["sqrt", 1.0],
                              "min_samples_leaf": [1, 2], "n_jobs": 1,
                              "reo_max_pairs": 3, "association": "within_tissue_pearson", "seed": 731},
        }
        self.train, self.test = make_grouped_splits(self.groups, self.tissue, 3, 732)[0]
        self.inner = make_grouped_splits(self.groups[self.train], self.tissue[self.train], 2, 733)

    def evaluate(self, kind, x=None, y=None):
        return evaluate_outer_fold(self.x if x is None else x, self.y if y is None else y,
                                   self.groups, self.train, self.test, self.inner,
                                   kind, self.config, self.ids)

    def test_alpha_max_is_independently_reconstructed_from_current_training_rows(self):
        for l1 in (0.5, 1.0):
            with self.subTest(l1_ratio=l1):
                model = TrainOnlyRegressor("expression", self.config, self.ids,
                                           alpha_ratio=0.1, l1_ratio=l1)
                model.fit(self.x[self.train], self.y[self.train])
                raw = np.log2(self.x[self.train, :-2] + 1.0)
                standardized = (raw - raw.mean(axis=0)) / raw.std(axis=0)
                response = self.y[self.train]
                response = (response - response.mean()) / response.std()
                zr = residualize_within_tissue(standardized, self.tissue[self.train])
                yr = residualize_within_tissue(response, self.tissue[self.train])
                expected = np.max(np.abs(zr.T @ yr)) / (len(self.train) * l1)
                audit = model.audit()
                self.assertAlmostEqual(audit["alpha_max"], expected, places=12)
                self.assertAlmostEqual(audit["actual_alpha"], expected * 0.1, places=12)

    def test_heldout_labels_and_expression_cannot_change_paths_tuning_or_pair_selection(self):
        changed_y = self.y.copy()
        changed_y[self.test] = 10000.0 + np.arange(len(self.test)) * 57.0
        changed_x = self.x.copy()
        changed_x[self.test, :-2] = np.array([1e12, 0., 1e9, 1e5, 1., 100.])
        for kind in ("expression", "expression_matched", "reo", "expression_rf", "reo_rf"):
            with self.subTest(kind=kind):
                original = self.evaluate(kind)
                label_counterfactual = self.evaluate(kind, y=changed_y)
                expression_counterfactual = self.evaluate(kind, x=changed_x)
                self.assertEqual(original["chosen_params"], label_counterfactual["chosen_params"])
                self.assertEqual(original["chosen_params"], expression_counterfactual["chosen_params"])
                np.testing.assert_array_equal(original["predictions"], label_counterfactual["predictions"])
                for counterfactual in (label_counterfactual, expression_counterfactual):
                    self.assertEqual(original["feature_audit"], counterfactual["feature_audit"])
                    if kind in {"expression", "expression_matched", "reo"}:
                        for first, second in zip(original["fit_events"], counterfactual["fit_events"]):
                            self.assertEqual(first["alpha_max"], second["alpha_max"])
                            self.assertEqual(first["actual_alpha"], second["actual_alpha"])

    def test_zero_alpha_max_retains_unpenalized_tissue_effects(self):
        tissue = np.arange(len(self.x)) % 2
        x = self.x.copy()
        x[:, -2] = tissue
        y = np.where(tissue == 0, -1.0, 1.0)
        model = TrainOnlyRegressor("expression", self.config, self.ids,
                                   alpha_ratio=0.1, l1_ratio=0.5)
        model.fit(x, y)
        self.assertEqual(model.audit()["alpha_max"], 0.0)
        self.assertEqual(model.audit()["n_nonzero_molecular_coeffs"], 0)
        np.testing.assert_allclose(model.predict(x), y, atol=1e-12, rtol=0)
        unknown = x[:3].copy()
        unknown[:, -2] = 99
        np.testing.assert_allclose(model.predict(unknown), y.mean(), atol=1e-12, rtol=0)

    def test_alpha_max_endpoint_has_exact_zero_coefficients_and_tissue_constant_predictions(self):
        for kind in ("expression", "reo"):
            for l1 in (0.5, 1.0):
                with self.subTest(kind=kind, l1_ratio=l1):
                    model = TrainOnlyRegressor(kind, self.config, self.ids, alpha_ratio=1.0, l1_ratio=l1)
                    model.fit(self.x[self.train], self.y[self.train])
                    audit = model.audit()
                    self.assertTrue(audit["null_endpoint"])
                    self.assertEqual(audit["solver_reason"], "theoretical_alpha_max_null_endpoint")
                    self.assertEqual(audit["fallback_reason"], "")
                    np.testing.assert_array_equal(model.regressor_.coef_, np.zeros_like(model.regressor_.coef_))
                    prediction = model.predict(self.x[self.test])
                    for category in np.unique(self.tissue[self.test]):
                        values = prediction[self.tissue[self.test] == category]
                        np.testing.assert_array_equal(values, np.full(len(values), values[0]))
                    probe = self.x[self.test].copy()
                    probe[:, :-2] = np.array([1e12, 0., 1e9, 1e5, 1., 100.])
                    np.testing.assert_array_equal(model.predict(probe), prediction)

    def test_rf_gene_builder_never_receives_tissue_or_audit_columns(self):
        original_fit = FoldFeatureBuilder.fit
        actual_rows = []
        def record(estimator, x, y=None):
            array = np.asarray(x)
            self.assertEqual(array.shape[1], len(self.ids))
            actual_rows.append(tuple(sorted((array[:, 0] - 100).astype(int).tolist())))
            return original_fit(estimator, x, y)
        with patch.object(FoldFeatureBuilder, "fit", record):
            self.evaluate("reo_rf")
        expected = [tuple(sorted(self.train[inner_train].tolist())) for inner_train, _ in self.inner]
        full_train = tuple(sorted(self.train.tolist()))
        self.assertEqual(actual_rows.count(full_train), 1)
        for rows in actual_rows:
            self.assertIn(rows, expected + [full_train])
            self.assertFalse(set(rows) & set(self.test.tolist()))

    def test_rf_pair_selector_receives_only_the_current_training_encoding_and_labels(self):
        original_fit = TrainOnlyReoPairSelector.fit
        calls = []
        def record(selector, z, y, *, tissue_design, feature_names):
            calls.append((np.array(z, copy=True), np.array(y, copy=True),
                          np.array(tissue_design, copy=True), tuple(feature_names)))
            return original_fit(selector, z, y, tissue_design=tissue_design, feature_names=feature_names)
        with patch.object(TrainOnlyReoPairSelector, "fit", record):
            result = self.evaluate("reo_rf")
        self.assertEqual(len(calls), len(result["fit_events"]))
        for (z, y, d, names), event in zip(calls, result["fit_events"]):
            rows = np.array(event["fit_sample_ids"], dtype=int)
            builder = FoldFeatureBuilder(kind="reo", n_genes=6, min_tpm=1., min_detection_fraction=0.1,
                                         max_pair_tie_fraction=0.9, stable_gene_ids=self.ids).fit(self.x[rows, :-2])
            np.testing.assert_array_equal(z, builder.transform(self.x[rows, :-2]))
            self.assertEqual(names, tuple(builder.get_feature_names_out()))
            expected_y = (self.y[rows] - self.y[rows].mean()) / self.y[rows].std()
            np.testing.assert_allclose(y, expected_y, atol=1e-12, rtol=0)
            expected_tissue = self.tissue[rows]
            self.assertEqual(d.shape, (len(rows), len(np.unique(expected_tissue))))
            for category, column in zip(np.unique(expected_tissue), d.T):
                np.testing.assert_array_equal(column, expected_tissue == category)

    def test_forest_receives_original_molecular_values_plus_training_one_hot(self):
        original_fit = RandomForestRegressor.fit
        for kind in ("reo_rf", "expression_rf"):
            with self.subTest(kind=kind):
                calls = []
                def record(forest, x, y, *args, **kwargs):
                    calls.append((np.array(x, copy=True), np.array(y, copy=True)))
                    return original_fit(forest, x, y, *args, **kwargs)
                model = TrainOnlyRegressor(kind, self.config, self.ids, max_features="sqrt", min_samples_leaf=1)
                with patch.object(RandomForestRegressor, "fit", record):
                    model.fit(self.x[self.train], self.y[self.train])
                self.assertEqual(len(calls), 1)
                molecular = model.builder_.transform(self.x[self.train, :-2])
                if kind == "reo_rf":
                    molecular = molecular[:, model.pair_selector_.selected_indices_]
                    self.assertTrue(set(np.unique(molecular)) <= {-1, 0, 1})
                one_hot = np.column_stack([self.tissue[self.train] == category
                                           for category in np.unique(self.tissue[self.train])]).astype(float)
                np.testing.assert_array_equal(calls[0][0], np.column_stack((molecular, one_hot)))
                expected_y = (self.y[self.train] - self.y[self.train].mean()) / self.y[self.train].std()
                np.testing.assert_allclose(calls[0][1], expected_y, atol=1e-12, rtol=0)


class SupervisedPairSelectionTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(826)
        self.tissue = np.repeat(np.arange(3), 16)
        self.d = np.eye(3)[self.tissue]
        self.z = rng.integers(-1, 2, size=(48, 8)).astype(np.int8)
        # A pair can vary globally while being constant within every tissue;
        # that apparent tissue-only association must be excluded.
        self.z[:, 7] = self.tissue - 1
        self.names = tuple(f"REO:G{i}>G{i + 10}" for i in range(8))
        self.y = 0.8 * self.z[:, 0] + 0.2 * self.z[:, 1] + 3.0 * self.tissue + rng.normal(0, 0.02, size=48)

    def fit(self, z=None, y=None, names=None):
        return TrainOnlyReoPairSelector(max_pairs=3, association="within_tissue_pearson").fit(
            self.z if z is None else z, self.y if y is None else y,
            tissue_design=self.d, feature_names=self.names if names is None else names)

    def test_pair_ranking_matches_independent_within_tissue_pearson(self):
        fitted = self.fit()
        zr = residualize_within_tissue(self.z, self.tissue)
        yr = residualize_within_tissue(self.y, self.tissue)
        associations = []
        for column in range(zr.shape[1]):
            if np.std(zr[:, column]) == 0:
                continue
            value = np.dot(zr[:, column], yr) / np.sqrt(np.sum(zr[:, column] ** 2) * np.sum(yr ** 2))
            associations.append((column, value))
        expected = sorted(associations, key=lambda item: (-abs(item[1]), self.names[item[0]]))[:3]
        np.testing.assert_array_equal(fitted.selected_indices_, [i for i, _ in expected])
        np.testing.assert_allclose(fitted.signed_associations_, [v for _, v in expected], atol=1e-12, rtol=0)
        self.assertNotIn(7, fitted.selected_indices_)
        np.testing.assert_array_equal(fitted.transform(self.z), self.z[:, fitted.selected_indices_])

    def test_transform_cannot_refit_supervised_pair_selection(self):
        fitted = self.fit()
        selected = fitted.selected_indices_.copy()
        names = tuple(fitted.selected_feature_names_)
        audit = copy.deepcopy(fitted.audit_metadata_)
        heldout = -self.z.copy()
        fitted.transform(heldout)
        np.testing.assert_array_equal(fitted.selected_indices_, selected)
        self.assertEqual(tuple(fitted.selected_feature_names_), names)
        self.assertEqual(fitted.audit_metadata_, audit)
        np.testing.assert_array_equal(fitted.transform(self.z), self.z[:, selected])

    def test_training_response_changes_select_the_corresponding_informative_pair(self):
        first = self.fit(y=self.z[:, 0].astype(float) + self.tissue * 10.)
        second = self.fit(y=self.z[:, 1].astype(float) + self.tissue * 10.)
        self.assertEqual(first.selected_indices_[0], 0)
        self.assertEqual(second.selected_indices_[0], 1)

    def test_zero_within_tissue_response_retains_no_molecular_pairs(self):
        fitted = self.fit(y=self.tissue.astype(float))
        self.assertEqual(len(fitted.selected_indices_), 0)
        self.assertEqual(fitted.transform(self.z).shape, (48, 0))

    def test_reordered_pair_columns_preserve_stable_feature_selection(self):
        first = self.fit()
        order = np.array([7, 2, 0, 6, 4, 1, 5, 3])
        second = self.fit(z=self.z[:, order], names=tuple(self.names[i] for i in order))
        self.assertEqual(tuple(first.selected_feature_names_), tuple(second.selected_feature_names_))
        np.testing.assert_array_equal(first.transform(self.z), second.transform(self.z[:, order]))


if __name__ == "__main__":
    unittest.main()
