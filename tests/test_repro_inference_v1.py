"""Contract tests using handwritten state only: no fit or real observations."""
from copy import deepcopy
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from repro_inference_v1 import (checked_response_free_path, load_bundle,
                                load_json_response_free, predict_bundle)


class PortableInferenceTests(unittest.TestCase):
    def setUp(self):
        # Keep temporary files under the workspace on restricted Windows hosts.
        (ROOT/'.pytest_cache').mkdir(exist_ok=True)
        self.bundle = load_bundle(ROOT/'examples/repro_inference_v1/toy_rank_model.json')
        self.sample = load_json_response_free(ROOT/'examples/repro_inference_v1/toy_input.json')

    def one_row(self, tpm=None, tissue='TOY_TISSUE_A'):
        return dict(schema_version='reo-inference-input-1.0', sample_ids=['handmade'],
                    biological_group_ids=['handmade_group'], tissues=[tissue],
                    gene_ids=self.bundle['input_gene_ids'],
                    raw_tpm=[[1., 1., 0., 3., 100.] if tpm is None else tpm])

    def null_bundle(self, declared=False, kind='tissue_null'):
        bundle = deepcopy(self.bundle)
        bundle.update(model_kind='tissue_null' if not declared else 'rank', fitted_kind=kind,
                      declared_molecular_arm=declared, candidate_gene_indices=[], candidate_gene_ids=[],
                      feature_names=[], pair_indices=[], feature_mean=[], feature_scale=[],
                      feature_tissue_means=[], coefficients=[], n_nonzero_molecular_coeffs=0)
        return bundle

    def test_hand_calculated_midrank_fwl_and_unseen(self):
        result = predict_bundle(self.bundle, self.sample)
        np.testing.assert_allclose([result['predictions'][i] for i in (0, 1, 2, 6)],
                                   [.395, .68, .5, .575], atol=1e-15, rtol=0)
        self.assertEqual(result['unknown_tissue'], [False, False, True, False, None, False, True])
        self.assertEqual(result['statuses'], ['predicted']*3+['abstained']*3+['predicted'])
        self.assertEqual(result['reasons'][3:6], ['missing_molecular_input', 'missing_or_invalid_tissue', 'negative_molecular_input'])

    def test_rank_needs_inactive_candidates_not_only_active_gene(self):
        # Only first coefficient is active, but other candidates change its rank.
        a = predict_bundle(self.bundle, self.one_row([1., 1., 0., 3., 100.]))['predictions'][0]
        b = predict_bundle(self.bundle, self.one_row([1., 2., 2., 3., 100.]))['predictions'][0]
        self.assertAlmostEqual(a-b, .075)

    def test_zero_molecular_and_true_null_have_different_gateway(self):
        bundle = deepcopy(self.bundle)
        bundle.update(coefficients=[0.]*4, n_nonzero_molecular_coeffs=0)
        sample = self.one_row([1., 1., 0., 3., -1.])
        zero = predict_bundle(bundle, sample)
        self.assertEqual(zero['predictions'], [None])
        self.assertTrue(zero['zero_molecular_solution'])
        null = predict_bundle(self.null_bundle(), sample)
        self.assertEqual(null['predictions'], [.4])
        self.assertIsNone(null['zero_molecular_solution'])
        sample['tissues'] = ['TOY_NEW_TISSUE']
        self.assertEqual(predict_bundle(self.null_bundle(), sample)['predictions'], [.5])

    def test_zero_valid_tpm_predicts_tissue_policy(self):
        bundle = deepcopy(self.bundle)
        bundle.update(coefficients=[0.]*4, n_nonzero_molecular_coeffs=0)
        self.assertEqual(predict_bundle(bundle, self.one_row())['predictions'], [.4])

    def test_empty_feature_fallback_is_still_declared_molecular(self):
        sample = self.one_row()
        sample['raw_tpm'] = [None]
        result = predict_bundle(self.null_bundle(declared=True), sample)
        self.assertEqual(result['statuses'], ['abstained'])
        self.assertEqual(result['reasons'], ['missing_molecular_input'])

    def test_mean_fallback_remains_mean_not_tissue_effect(self):
        result = predict_bundle(self.null_bundle(declared=True, kind='mean'), self.one_row())
        self.assertEqual(result['predictions'], [.5])

    def test_expression_encoding_scale_then_projection(self):
        bundle = deepcopy(self.bundle)
        bundle.update(model_kind='expression', fitted_kind='expression',
                      feature_names=['log2TPM:'+x for x in bundle['candidate_gene_ids']],
                      feature_mean=[0.]*4, feature_scale=[2., 1., 1., 1.])
        # log2(1+1)/2 = .5; then subtract tissue projection .1.
        self.assertAlmostEqual(predict_bundle(bundle, self.one_row())['predictions'][0], .42)

    def test_pair_orientation_and_exact_ties(self):
        bundle = deepcopy(self.bundle)
        pairs = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        bundle.update(model_kind='pair_full', fitted_kind='pair_full', pair_indices=pairs,
                      feature_names=['REO:'+bundle['input_gene_ids'][a]+'>'+bundle['input_gene_ids'][b] for a, b in pairs],
                      feature_mean=[0.]*6, feature_scale=[1.]*6, feature_tissue_means=[[0.]*6]*2,
                      coefficients=[1., 2., 3., 0., 0., 0.], n_nonzero_molecular_coeffs=3)
        # The first three pair codes are [0,+1,-1].
        self.assertAlmostEqual(predict_bundle(bundle, self.one_row())['predictions'][0], .2)
        # A strict nextafter change must break the first exact tie, without epsilon.
        bundle['coefficients'] = [1., 0., 0., 0., 0., 0.]
        bundle['n_nonzero_molecular_coeffs'] = 1
        changed = self.one_row([float(np.nextafter(1., 2.)), 1., 0., 3., 100.])
        self.assertAlmostEqual(predict_bundle(bundle, changed)['predictions'][0], .6)

    def test_labels_and_schema_order_are_rejected(self):
        sample = self.one_row()
        sample['response_auc'] = [.7]
        with self.assertRaises(ValueError): predict_bundle(self.bundle, sample)
        sample = self.one_row()
        sample['gene_ids'] = sample['gene_ids'][::-1]
        with self.assertRaises(ValueError): predict_bundle(self.null_bundle(), sample)

    def test_structural_and_duplicate_identity_failures(self):
        for field in ('sample_ids', 'biological_group_ids'):
            sample = deepcopy(self.sample)
            sample[field][1] = sample[field][0]
            with self.assertRaises(ValueError): predict_bundle(self.bundle, sample)
        sample = self.one_row()
        sample['raw_tpm'] = [[1., 2.]]
        with self.assertRaises(ValueError): predict_bundle(self.bundle, sample)
        sample = self.one_row()
        sample['raw_tpm'][0][0] = [1.]
        with self.assertRaises(ValueError): predict_bundle(self.null_bundle(), sample)

    def test_invalid_identity_or_tissue_abstains_for_null_too(self):
        for field, reason in [('sample_ids', 'missing_or_invalid_identity'),
                              ('biological_group_ids', 'missing_or_invalid_identity'),
                              ('tissues', 'missing_or_invalid_tissue')]:
            sample = self.one_row()
            sample[field] = [None]
            self.assertEqual(predict_bundle(self.null_bundle(), sample)['reasons'], [reason])

    def test_missing_nonnumeric_nonfinite_and_negative_are_retained(self):
        for value, reason in [(None, 'missing_molecular_input'), ('bad', 'nonnumeric_molecular_input'),
                              (float('nan'), 'nonfinite_molecular_input'), (float('inf'), 'nonfinite_molecular_input'),
                              (-1., 'negative_molecular_input'), (True, 'nonnumeric_molecular_input')]:
            sample = self.one_row()
            sample['raw_tpm'][0][-1] = value  # not a selected gene
            result = predict_bundle(self.bundle, sample)
            self.assertEqual(result['predictions'], [None])
            self.assertEqual(result['reasons'], [reason])

    def test_invalid_bundle_dimensions_and_scales_rejected(self):
        for key, value in [('feature_scale', [0.]*4), ('coefficients', [0.]*3),
                           ('feature_tissue_means', [[0.]*4]), ('target_sd', 0.),
                           ('n_nonzero_molecular_coeffs', 0), ('training_tissue_codes', [1., 0.])]:
            bundle = deepcopy(self.bundle)
            bundle[key] = value
            with self.assertRaises(ValueError): predict_bundle(bundle, self.one_row())

    def test_json_gzip_and_duplicate_or_nonstandard_json_rejection(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'.pytest_cache') as folder:
            self.assertTrue(Path(folder).resolve().is_relative_to(ROOT))
            path = Path(folder)/'toy.json.gz'
            with gzip.open(path, 'wt', encoding='utf-8') as handle: json.dump(self.bundle, handle)
            self.assertEqual(load_bundle(path)['coefficients'], self.bundle['coefficients'])
            bad = Path(folder)/'bad.json'
            for content in ('{"a":1,"a":2}', '{"a":NaN}'):
                bad.write_text(content, encoding='utf-8')
                with self.assertRaises(ValueError): load_json_response_free(bad)

    def test_sealed_paths_blocked_without_opening(self):
        for path in ('data_raw/ctrp/does_not_exist.json', 'data_processed/ctrp_curve_index.json',
                     'runs/fake-external/none.json', 'legacy/CTRP--Drug--Cellline--AUC.txt'):
            with self.assertRaises(PermissionError): checked_response_free_path(ROOT/path)

    def test_cli_creates_new_file_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'.pytest_cache') as folder:
            self.assertTrue(Path(folder).resolve().is_relative_to(ROOT))
            output = Path(folder)/'new.json'
            command = [sys.executable, str(ROOT/'scripts/predict_repro_v1.py'),
                       '--model', str(ROOT/'examples/repro_inference_v1/toy_rank_model.json'),
                       '--input', str(ROOT/'examples/repro_inference_v1/toy_input.json'), '--output', str(output)]
            first = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            original = output.read_bytes()
            second = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(second.returncode, 2)
            self.assertEqual(output.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
