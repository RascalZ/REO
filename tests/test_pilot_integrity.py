"""Synthetic, counterfactual acceptance tests for the development-only pilot.

These tests deliberately do not load expression or response values from the
project cohorts. Changes to held-out information must not affect fitted
training artifacts; repeated biological groups must never straddle folds.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from openpyxl import Workbook
from sklearn.exceptions import NotFittedError

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from modeling.features import FoldFeatureBuilder
from modeling.nested_cv import TrainOnlyRegressor, evaluate_outer_fold, make_grouped_splits
from run_development_pilot import read_selected_responses


class FoldFeatureIntegrityTests(unittest.TestCase):
    def setUp(self):
        # Scrambled gene columns, multiple exact ties, varying pair direction,
        # a completely undetected gene, and a detected but constant gene.
        self.ids = ("ENSG000003", "ENSG000001", "ENSG000002", "ENSG000005", "ENSG000004")
        self.x = np.array([
            [0, 0, 3, 0, 10],
            [1, 2, 3, 0, 10],
            [4, 4, 6, 0, 10],
            [2, 1, 4, 0, 10],
            [3, 5, 3, 0, 10],
            [8, 7, 1, 0, 10],
        ], dtype=np.float64)

    def builder(self, kind, **overrides):
        kwargs = dict(kind=kind, n_genes=3, min_tpm=1.0,
                      min_detection_fraction=0.1, max_pair_tie_fraction=0.9,
                      stable_gene_ids=self.ids)
        kwargs.update(overrides)
        return FoldFeatureBuilder(**kwargs)

    def test_reo_has_exact_three_value_encoding_and_canonical_direction(self):
        fitted = self.builder("reo").fit(self.x)
        actual = fitted.transform(self.x)
        self.assertEqual(actual.dtype, np.dtype("int8"))
        self.assertEqual(set(np.unique(actual)), {-1, 0, 1})
        self.assertEqual(len(fitted.pair_indices_), 3)
        for position, (left, right) in enumerate(fitted.pair_indices_):
            self.assertLess(self.ids[left], self.ids[right])
            expected = np.sign(self.x[:, left] - self.x[:, right]).astype(np.int8)
            np.testing.assert_array_equal(actual[:, position], expected)
            np.testing.assert_array_equal(actual[self.x[:, left] == self.x[:, right], position], 0)

    def test_column_reordering_preserves_gene_identity_and_pair_direction(self):
        original = self.builder("reo").fit(self.x)
        order = np.array([2, 4, 0, 3, 1])
        reordered_ids = tuple(self.ids[i] for i in order)
        reordered = self.builder("reo", stable_gene_ids=reordered_ids).fit(self.x[:, order])
        self.assertEqual(list(original.get_feature_names_out()), list(reordered.get_feature_names_out()))
        np.testing.assert_array_equal(original.transform(self.x), reordered.transform(self.x[:, order]))

    def test_heldout_distribution_cannot_change_training_artifacts(self):
        heldout_a = self.x[[0, 1]].copy()
        # The originally undetected and constant columns now have extreme
        # held-out variance. A global or transform-time selector would retain
        # them, displacing legitimate training features.
        heldout_b = np.array([[1e9, 1e9, 1e9, 0, 1e9], [0, 0, 0, 1e9, 0]], dtype=np.float64)
        for kind in ("expression", "reo"):
            with self.subTest(kind=kind):
                fitted = self.builder(kind).fit(self.x)
                gene_ids_before = tuple(fitted.gene_ids_)
                names_before = tuple(fitted.get_feature_names_out())
                stats_before = copy.deepcopy(fitted.training_statistics_)
                training_before = fitted.transform(self.x).copy()
                fitted.transform(heldout_a)
                fitted.transform(heldout_b)
                self.assertEqual(tuple(fitted.gene_ids_), gene_ids_before)
                self.assertEqual(tuple(fitted.get_feature_names_out()), names_before)
                self.assertEqual(fitted.training_statistics_, stats_before)
                np.testing.assert_array_equal(fitted.transform(self.x), training_before)
                self.assertNotIn("ENSG000004", fitted.gene_ids_)
                self.assertNotIn("ENSG000005", fitted.gene_ids_)

    def test_feature_selection_is_response_label_independent(self):
        for kind in ("expression", "reo"):
            with self.subTest(kind=kind):
                first = self.builder(kind).fit(self.x, np.arange(len(self.x), dtype=float))
                second = self.builder(kind).fit(self.x, np.array([999., -10., 2., 60., -1., 8.]))
                self.assertEqual(list(first.get_feature_names_out()), list(second.get_feature_names_out()))
                np.testing.assert_array_equal(first.transform(self.x), second.transform(self.x))

    def test_expression_uses_fold_selected_log2_tpm(self):
        fitted = self.builder("expression").fit(self.x)
        np.testing.assert_array_equal(fitted.transform(self.x), np.log2(self.x[:, fitted.gene_indices_] + 1.0))

    def test_training_constant_and_high_tie_pairs_are_not_reintroduced(self):
        ids = ("G1", "G2", "G3")
        constant_pair = np.array([[1, 2, 3], [2, 3, 1], [3, 4, 2],
                                  [4, 5, 8], [5, 6, 7], [6, 7, 0]], dtype=np.float64)
        high_tie_pair = constant_pair.copy()
        high_tie_pair[:, 1] = high_tie_pair[:, 0]
        high_tie_pair[-1, 1] = 7.0
        for fixture, tie_limit in ((constant_pair, 0.9), (high_tie_pair, 0.5)):
            with self.subTest(tie_limit=tie_limit):
                fitted = self.builder("reo", stable_gene_ids=ids,
                                      max_pair_tie_fraction=tie_limit).fit(fixture)
                names = list(fitted.get_feature_names_out())
                self.assertNotIn("REO:G1>G2", names)
                # The held-out direction reverses and the previously high-tie
                # pair becomes fully variable. Its exclusion remains frozen.
                heldout = np.array([[20., 1., 10.], [1., 20., 2.]], dtype=np.float64)
                transformed = fitted.transform(heldout)
                self.assertEqual(transformed.shape[1], len(names))
                self.assertEqual(list(fitted.get_feature_names_out()), names)

    def test_no_eligible_features_fail_explicitly(self):
        for kind in ("expression", "reo"):
            with self.subTest(kind=kind):
                with self.assertRaises(ValueError):
                    self.builder(kind).fit(np.zeros_like(self.x))

    def test_invalid_expression_or_gene_identity_is_rejected(self):
        for bad_value in (-1.0, np.nan, np.inf):
            invalid = self.x.copy()
            invalid[0, 0] = bad_value
            with self.subTest(value=str(bad_value)):
                with self.assertRaises(ValueError):
                    self.builder("reo").fit(invalid)
        duplicate_ids = (self.ids[0], self.ids[0], *self.ids[2:])
        with self.assertRaises(ValueError):
            self.builder("reo", stable_gene_ids=duplicate_ids).fit(self.x)

    def test_failed_refit_cannot_reuse_an_earlier_fold_feature_set(self):
        fitted = self.builder("reo").fit(self.x)
        invalid = self.x.copy()
        invalid[0, 0] = -1.0
        with self.assertRaises(ValueError):
            fitted.fit(invalid)
        with self.assertRaises(NotFittedError):
            fitted.transform(self.x)


class NestedBoundaryIntegrityTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(919)
        self.groups = np.repeat(np.arange(12), 2).astype(str)
        self.tissue = np.repeat(np.arange(12) % 3, 2)
        self.sample_ids = np.arange(len(self.groups))
        genes = rng.lognormal(1.5, 1.0, size=(len(self.groups), 4)).astype(np.float64)
        # A unique first-gene value makes actual selector training rows
        # identifiable without using the audit ID columns in model features.
        genes[:, 0] = 100.0 + self.sample_ids
        self.x = np.column_stack((genes, self.tissue, self.sample_ids)).astype(np.float64)
        self.y = (0.2 * np.log2(genes[:, 1] + 1.0) + 0.1 * self.tissue
                  + rng.normal(0, 0.03, size=len(self.groups)))
        self.ids = ("G1", "G2", "G3", "G4")
        self.config = {
            "feature": {"reo_genes": 4, "expression_genes": 4,
                        "min_tpm": 1.0, "min_detection_fraction": 0.1,
                        "max_pair_tie_fraction": 0.9},
            "elastic_net": {"alphas": [0.1, 0.01], "l1_ratios": [0.5],
                            "max_iter": 1000, "tol": 0.001},
            "ridge_alphas": [0.1, 1.0],
        }
        self.outer_splits = make_grouped_splits(self.groups, self.tissue, 3, 173)
        self.train, self.test = self.outer_splits[0]
        self.inner = make_grouped_splits(self.groups[self.train], self.tissue[self.train], 2, 174)

    def evaluate(self, kind, y=None):
        return evaluate_outer_fold(self.x, self.y if y is None else y,
                                   self.groups, self.train, self.test, self.inner,
                                   kind, self.config, self.ids)

    def test_grouped_folds_are_disjoint_complete_and_reproducible(self):
        groups = np.repeat(np.arange(12), np.tile([1, 2, 3, 4], 3)).astype(str)
        tissue = np.concatenate([np.full(count, group % 3)
                                 for group, count in enumerate(np.tile([1, 2, 3, 4], 3))])
        first = make_grouped_splits(groups, tissue, 3, 391)
        second = make_grouped_splits(groups, tissue, 3, 391)
        outer_tests = []
        for (train, test), (again_train, again_test) in zip(first, second):
            np.testing.assert_array_equal(train, again_train)
            np.testing.assert_array_equal(test, again_test)
            self.assertFalse(set(groups[train]) & set(groups[test]))
            self.assertEqual(set(train) | set(test), set(range(len(groups))))
            self.assertFalse(set(train) & set(test))
            outer_tests.extend(test.tolist())
            for inner_train, inner_test in make_grouped_splits(groups[train], tissue[train], 2, 392):
                self.assertFalse(set(groups[train][inner_train]) & set(groups[train][inner_test]))
        self.assertEqual(sorted(outer_tests), list(range(len(groups))))

    def test_outer_test_labels_cannot_change_tuning_or_predictions(self):
        changed = self.y.copy()
        changed[self.test] = 5000.0 + np.arange(len(self.test)) * 300.0
        for kind in ("mean", "tissue", "expression", "reo", "expression_matched"):
            with self.subTest(kind=kind):
                original = self.evaluate(kind)
                counterfactual = self.evaluate(kind, changed)
                self.assertEqual(original["chosen_params"], counterfactual["chosen_params"])
                np.testing.assert_array_equal(original["predictions"], counterfactual["predictions"])

    def test_each_inner_fit_rebuilds_features_from_its_own_training_rows(self):
        original_fit = FoldFeatureBuilder.fit
        for kind in ("expression", "reo", "expression_matched"):
            with self.subTest(kind=kind):
                actual_fits = []

                def recording_fit(estimator, x, y=None):
                    array = np.asarray(x)
                    self.assertEqual(array.shape[1], len(self.ids),
                                     "tissue/audit columns entered the feature selector")
                    row_ids = tuple(sorted((array[:, 0] - 100.0).astype(int).tolist()))
                    actual_fits.append(row_ids)
                    return original_fit(estimator, x, y)

                with patch.object(FoldFeatureBuilder, "fit", recording_fit):
                    self.evaluate(kind)
                full_training = tuple(sorted(self.train.tolist()))
                expected_inner = [tuple(sorted(self.train[inner_train].tolist()))
                                  for inner_train, _ in self.inner]
                for observed in actual_fits:
                    self.assertIn(observed, expected_inner + [full_training])
                    self.assertFalse(set(observed) & set(self.test.tolist()))
                self.assertEqual(actual_fits.count(full_training), 1,
                                 "outer-training selector may only be fitted for the final refit")
                for membership in expected_inner:
                    self.assertGreaterEqual(actual_fits.count(membership), 1,
                                            "each inner training partition needs its own fitted selector")

    def test_unseen_tissue_baseline_uses_training_global_mean(self):
        model = TrainOnlyRegressor("tissue", self.config, self.ids, alpha=1.0)
        model.fit(self.x[self.train], self.y[self.train])
        unseen = self.x[self.test].copy()
        unseen[:, -2] = 999.0
        np.testing.assert_array_equal(model.predict(unseen),
                                      np.full(len(unseen), self.y[self.train].mean()))
        self.assertNotIn(999.0, model.training_tissue_codes_)

    def test_unseen_tissue_in_gene_model_retains_gene_component(self):
        model = TrainOnlyRegressor("expression", self.config, self.ids,
                                   alpha=0.01, l1_ratio=0.5)
        model.fit(self.x[self.train], self.y[self.train])
        unseen = self.x[self.test].copy()
        unseen[:, -2] = 999.0
        different_unknown_label = unseen.copy()
        different_unknown_label[:, -2] = 1000.0
        np.testing.assert_array_equal(model.predict(unseen), model.predict(different_unknown_label))
        self.assertGreater(np.ptp(model.predict(unseen)), 1e-8,
                           "unknown tissue incorrectly discarded the fitted gene component")
        self.assertNotIn(999.0, model.training_tissue_codes_)
        self.assertNotIn(1000.0, model.training_tissue_codes_)

    def test_heldout_expression_cannot_change_training_feature_audit(self):
        shifted = self.x.copy()
        shifted[self.test, :4] = np.array([1e12, 0., 1e9, 1e5])
        for kind in ("expression", "reo", "expression_matched"):
            with self.subTest(kind=kind):
                original = self.evaluate(kind)
                counterfactual = evaluate_outer_fold(shifted, self.y, self.groups,
                                                      self.train, self.test, self.inner,
                                                      kind, self.config, self.ids)
                self.assertEqual(original["chosen_params"], counterfactual["chosen_params"])
                self.assertEqual(original["feature_audit"], counterfactual["feature_audit"])

    def test_audit_sample_ids_cannot_affect_predictions(self):
        alternate_ids = self.x.copy()
        alternate_ids[:, -1] = 1e8 + self.sample_ids * 17.0
        for kind in ("tissue", "expression", "reo"):
            with self.subTest(kind=kind):
                original = self.evaluate(kind)
                counterfactual = evaluate_outer_fold(alternate_ids, self.y, self.groups,
                                                      self.train, self.test, self.inner,
                                                      kind, self.config, self.ids)
                self.assertEqual(original["chosen_params"], counterfactual["chosen_params"])
                np.testing.assert_array_equal(original["predictions"], counterfactual["predictions"])


class DevelopmentAccessGuardTests(unittest.TestCase):
    def test_guard_blocks_external_response_reads_and_allows_development_inputs(self):
        # Audit hooks cannot be removed, so this test runs in a fresh Python
        # process and uses temporary mock files. No real cohort data is read.
        script = r'''
import json, os, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import run_development_pilot as runner
root = Path(sys.argv[2]).resolve()
runner.ROOT = root
mock_originals = {"gdsc2": root / "original" / "gdsc.xlsx",
                  "legacy_ctrp_response": root / "original" / "legacy_ctrp.tsv"}
runner.original_input = lambda name: mock_originals[name]
forbidden = [root / "data_raw" / "ctrp" / "curves.txt",
             mock_originals["legacy_ctrp_response"],
             root / "data_processed" / "ctrp_curve_index.csv",
             root / "data_processed" / "ctrp_pair_availability.csv",
             root / "data_processed" / "ctrp_drug_availability.csv"]
allowed = [root / "data_processed" / "expression_genes.tsv",
           mock_originals["gdsc2"], root / "results" / "development_sample_manifest.csv"]
for path in forbidden + allowed:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("mock data only\n", encoding="utf-8")
access, blocked = runner.install_data_guard()
for path in forbidden:
    try:
        path.read_bytes()
    except PermissionError:
        pass
    else:
        raise AssertionError("External response read was not blocked: " + str(path))
try:
    descriptor = os.open(os.fsencode(forbidden[0]), os.O_RDONLY)
except PermissionError:
    pass
else:
    os.close(descriptor)
    raise AssertionError("Byte-path os.open external response read was not blocked")
for path in allowed:
    assert path.read_text(encoding="utf-8") == "mock data only\n"
assert str(allowed[0]) in access and str(allowed[1]) in access
assert len(blocked) == 6
print(json.dumps({"blocked_attempts": len(blocked), "allowed_development_reads": len(allowed)}))
'''
        with tempfile.TemporaryDirectory(prefix="reo_pilot_guard_") as temp:
            result = subprocess.run([sys.executable, "-c", script, str(REPO / "src"), temp],
                                    text=True, capture_output=True, timeout=30, check=False)
        self.assertEqual(result.returncode, 0, msg=result.stderr or result.stdout)
        report = json.loads(result.stdout)
        self.assertEqual(report["blocked_attempts"], 6)
        self.assertEqual(report["allowed_development_reads"], 3)


class FrozenResponseRowTests(unittest.TestCase):
    @staticmethod
    def write_fixture(path, selected_response=0.52):
        workbook = Workbook()
        sheet = workbook.active
        sheet.append(["DRUG_ID", "NLME_CURVE_ID", "MIN_CONC", "MAX_CONC", "AUC"])
        # Deliberately unusable unselected response: only the frozen row can
        # enter the pilot output or endpoint eligibility checks.
        sheet.append([999, 200, 0.1, 10.0, "not a numeric response"])
        sheet.append([1036, 201, 0.01, 10.0, selected_response])
        workbook.save(path)
        workbook.close()

    @staticmethod
    def frozen_row(**changes):
        row = {"gdsc_source_row_number": "3", "gdsc_drug_id": "1036",
               "nlme_curve_id": "201", "min_conc": "0.01", "max_conc": "10.0",
               "canonical_model_id": "ACH-SYNTHETIC", "role": "development"}
        return pd.DataFrame([{**row, **changes}])

    def test_only_frozen_response_rows_are_extracted(self):
        with tempfile.TemporaryDirectory(prefix="reo_response_fixture_") as temp:
            source = Path(temp) / "source.xlsx"
            output = Path(temp) / "development.csv"
            self.write_fixture(source)
            result = read_selected_responses(self.frozen_row(), source, output)
            self.assertEqual(len(result), 1)
            self.assertEqual(result.iloc[0].gdsc_drug_id, "1036")
            self.assertEqual(result.iloc[0].response_auc, 0.52)
            self.assertEqual(result.iloc[0].role, "development")

    def test_wrong_drug_curve_or_dose_window_is_rejected(self):
        for change in ({"gdsc_drug_id": "999"}, {"nlme_curve_id": "999"},
                       {"min_conc": "0.1"}, {"max_conc": "1.0"}):
            with self.subTest(change=change):
                with tempfile.TemporaryDirectory(prefix="reo_wrong_row_fixture_") as temp:
                    source = Path(temp) / "source.xlsx"
                    self.write_fixture(source)
                    with self.assertRaises(AssertionError):
                        read_selected_responses(self.frozen_row(**change), source,
                                                Path(temp) / "development.csv")

    def test_duplicate_or_absent_frozen_source_rows_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix="reo_missing_row_fixture_") as temp:
            source = Path(temp) / "source.xlsx"
            self.write_fixture(source)
            duplicate = pd.concat([self.frozen_row(), self.frozen_row()], ignore_index=True)
            with self.assertRaises(AssertionError):
                read_selected_responses(duplicate, source, Path(temp) / "duplicate.csv")
            with self.assertRaises(AssertionError):
                read_selected_responses(self.frozen_row(gdsc_source_row_number="999"), source,
                                        Path(temp) / "missing.csv")

    def test_out_of_domain_development_auc_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="reo_bad_auc_fixture_") as temp:
            source = Path(temp) / "source.xlsx"
            self.write_fixture(source, selected_response=1.5)
            with self.assertRaises(AssertionError):
                read_selected_responses(self.frozen_row(), source, Path(temp) / "development.csv")


if __name__ == "__main__":
    unittest.main()
