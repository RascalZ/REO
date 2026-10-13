"""Numerical regression tests for the versioned ordinal diagnostic module."""
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from ordinal_diagnostic_models import DiagnosticFeatureBuilder, DiagnosticRegressor, evaluate_diagnostic_outer_fold


CONFIG = {'feature': {'diagnostic_genes': 4, 'min_tpm': 1.0, 'min_detection_fraction': 0.1},
          'elastic_net': {'regularization_mode': 'alpha_max_relative', 'alpha_ratios': [1.0, 0.3],
                          'l1_ratios': [0.5, 1.0], 'max_iter': 4000, 'tol': 0.001, 'selection': 'cyclic'}}


class OrdinalFeatureTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(112)
        self.x = np.exp(rng.normal(1.0, 0.8, (24, 7))).astype(np.float64)
        self.ids = tuple(f'G{i:02d}' for i in range(7))
        self.tissue = np.repeat([0., 1.], 12)
        self.y = 0.3+0.03*np.log2(self.x[:, 1]+1)+self.tissue*0.25+rng.normal(0, .005, 24)
        self.design = np.column_stack((self.x, self.tissue, np.arange(24)))

    def test_full_pairs_and_average_ranks_are_invertible_with_ties(self):
        x = np.asarray([[2, 2, 6, 1], [3, 9, 9, 3], [8, 1, 1, 8], [1, 6, 3, 6]], dtype=np.float64)
        kwargs = dict(n_genes=4, stable_gene_ids=('A', 'B', 'C', 'D'), min_tpm=0, min_detection_fraction=0)
        pair = DiagnosticFeatureBuilder(kind='pair_full', **kwargs).fit(x)
        rank = DiagnosticFeatureBuilder(kind='rank', **kwargs).fit(x)
        p, r = pair.transform(x), rank.transform(x)
        np.testing.assert_array_equal(rank.gene_indices_, pair.gene_indices_)
        reconstructed = np.zeros((len(x), 4, 4), dtype=np.int8)
        for k, (i, j) in enumerate(pair.pair_indices_):
            reconstructed[:, i, j] = p[:, k]
            reconstructed[:, j, i] = -p[:, k]
        np.testing.assert_array_equal(2.5+reconstructed.sum(axis=2)/2, r)
        for k, (i, j) in enumerate(pair.pair_indices_):
            np.testing.assert_array_equal(np.sign(r[:, i]-r[:, j]), p[:, k])

    def test_gene_selection_is_shared_and_application_cannot_refit_it(self):
        for rule in ('mad', 'response_association'):
            built = [DiagnosticFeatureBuilder(kind=kind, selection_rule=rule, n_genes=4, stable_gene_ids=self.ids).fit(self.x, self.y, self.tissue) for kind in ('expression', 'rank', 'pair_full')]
            self.assertEqual(built[0].gene_ids_, built[1].gene_ids_)
            self.assertEqual(built[0].gene_ids_, built[2].gene_ids_)
            before = built[1].get_audit_metadata(include_feature_names=True)
            external = self.x[:3].copy()*1e9
            built[1].transform(external)
            self.assertEqual(before, built[1].get_audit_metadata(include_feature_names=True))
            np.testing.assert_array_equal(built[1].transform(self.x), rankdata(self.x[:, built[1].gene_indices_], axis=1, method='average'))

    def test_analytical_null_predictions_match_policy_for_every_representation(self):
        null = DiagnosticRegressor('tissue_null', config=CONFIG, stable_gene_ids=self.ids).fit(self.design, self.y)
        xt = self.design[:4].copy()
        xt[:, -2] = [0., 1., 9., 9.]
        expected = null.predict(xt)
        for rule in ('mad', 'response_association'):
            for kind in ('expression', 'rank', 'pair_full'):
                model = DiagnosticRegressor(kind, selection_rule=rule, config=CONFIG, stable_gene_ids=self.ids, alpha_ratio=1., l1_ratio=.5).fit(self.design, self.y)
                np.testing.assert_array_equal(model.predict(xt), expected)
                self.assertEqual(model.n_nonzero_molecular_coeffs_, 0)
                self.assertEqual(model.audit()['prediction_sample_ids'], [0, 1, 2, 3])
                self.assertEqual(len(model.audit()['fit_response_sha256']), 64)

    def test_outer_outcomes_do_not_change_fit_or_prediction(self):
        train, test = np.arange(18), np.arange(18, 24)
        inner = [(np.setdiff1d(np.arange(18), val), val) for val in np.array_split(np.arange(18), 3)]
        common = dict(X=self.design, groups=np.arange(24), train_idx=train, test_idx=test,
                      inner_splits=inner, model_kind='rank', selection_rule='response_association', config=CONFIG, stable_gene_ids=self.ids)
        first = evaluate_diagnostic_outer_fold(y=self.y, **common)
        changed = self.y.copy()
        changed[test] += 1000
        second = evaluate_diagnostic_outer_fold(y=changed, **common)
        np.testing.assert_array_equal(first['predictions'], second['predictions'])
        self.assertEqual(first['chosen_params'], second['chosen_params'])
        self.assertEqual(first['feature_audit'], second['feature_audit'])
        self.assertEqual(len(first['fit_events']), 13)
        for event in first['fit_events']:
            self.assertFalse(set(event['fit_sample_ids']) & set(test))

    def test_no_mad_reselection_after_response_association(self):
        x = np.column_stack((np.arange(12)+100., np.arange(12)%3+1., np.arange(12)%2+3.)).astype(np.float64)
        y = x[:, 1].copy()
        builder = DiagnosticFeatureBuilder(kind='expression', selection_rule='response_association', n_genes=1, stable_gene_ids=('A', 'B', 'C')).fit(x, y, np.zeros(12))
        self.assertEqual(builder.gene_ids_, ('B',))


if __name__ == '__main__':
    unittest.main()
