"""Adversarial synthetic tests for the independent final-model auditor."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))

from final_model_v1 import make_global_inner_manifest
from ordinal_repeated_runtime import sha
from run_final_model_v1 import execute_fit_task, execute_verify_task, save_checkpoint
import validate_final_model_v1 as auditor


class IndependentFinalModelAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run = self.root/'runs/synthetic-final'
        (self.run/'checkpoints').mkdir(parents=True)
        (self.run/'models').mkdir()
        (self.root/'data_processed').mkdir()
        rng = np.random.default_rng(132)
        n = 30
        self.stable = tuple(f'ENSG_SYN_{i:02d}' for i in range(6))
        self.raw = rng.lognormal(mean=1., size=(n, 6)).astype(np.float64)
        np.save(self.root/'data_processed/expression_tpm.npy', self.raw)
        pd.DataFrame(dict(ensembl_id=self.stable,
                          gene_type=['protein_coding']*6)).to_csv(
            self.root/'data_processed/expression_genes.tsv', sep='\t', index=False)
        tissue = np.tile(['A', 'B', 'C'], 10)
        self.data = pd.DataFrame(dict(gdsc_drug_id=['1']*n,
            canonical_model_id=[f'SYN-{i:03d}' for i in range(n)],
            biological_group_id=[f'GROUP-{i:03d}' for i in range(n)],
            tissue=tissue, expression_row_index=np.arange(n),
            response_auc=np.round(.4+.04*np.log2(self.raw[:, 0]+1)
                +rng.normal(0, .01, n), 10)))
        self.data.to_csv(self.run/'development_responses.csv', index=False, float_format='%.17g')
        self.data = pd.read_csv(self.run/'development_responses.csv',
            dtype={'gdsc_drug_id': str, 'canonical_model_id': str,
                   'biological_group_id': str, 'tissue': str},
            keep_default_na=False, float_precision='round_trip')
        self.y = self.data.response_auc.to_numpy(dtype=np.float64)
        self.X = np.column_stack((self.raw,
            self.data.tissue.map({'A': 0, 'B': 1, 'C': 2}),
            np.arange(n))).astype(np.float64)
        self.manifest, _ = make_global_inner_manifest(
            self.data[['canonical_model_id', 'biological_group_id', 'tissue']],
            n_folds=3, seed=2026101501)
        self.manifest.to_csv(self.run/'global_inner_fold_manifest.csv', index=False)
        self.ordinal = dict(inner_folds=3, selection_rules=['mad', 'response_association'],
            representations=['expression', 'rank', 'pair_full'],
            feature=dict(n_genes=5, strong_expression_genes=6,
                         min_tpm=1., min_detection_fraction=.1),
            elastic_net=dict(regularization_mode='alpha_max_relative',
                alpha_ratios=[1., .3, .1, .03, .01], l1_ratios=[.5, 1.],
                max_iter=4000, tol=.001, selection='cyclic'),
            resources=dict(rss_sampling_interval_seconds=.01))
        self.task = dict(drug_id=1, model='response_association_pair_full',
                         selection_rule='response_association',
                         representation='pair_full', n_genes=5)
        self.context = dict(freeze_sha256='synthetic-frozen',
            ordinal_config=self.ordinal,
            external_config={'final_training': {'inner_folds': 3}},
            tissue_code={'A': 0, 'B': 1, 'C': 2},
            development_file_hashes={
                'data_processed/expression_genes.tsv': sha(self.root/'data_processed/expression_genes.tsv'),
                'data_processed/expression_tpm.npy': sha(self.root/'data_processed/expression_tpm.npy')})
        with patch('run_final_model_v1.load_drug', return_value=(
                self.data, self.X, self.y, self.stable)):
            self.assertTrue(execute_fit_task(self.run, self.task, self.context))
            self.assertTrue(execute_verify_task(self.run, self.task, self.context))
        self.checkpoint_path = self.run/'checkpoints/d1_response_association_pair_full.json'
        self.checkpoint = json.loads(self.checkpoint_path.read_text(encoding='utf-8'))
        self.grid = [dict(alpha_ratio=ratio, l1_ratio=l1)
            for ratio in self.ordinal['elastic_net']['alpha_ratios']
            for l1 in self.ordinal['elastic_net']['l1_ratios']]

    def test_independent_checkpoint_membership_grid_and_artifact_contract(self):
        obj = auditor.verify_checkpoint(self.checkpoint_path, self.task,
            'synthetic-frozen', self.data, self.manifest, self.grid, self.run)
        self.assertEqual(obj['fit_count'], 31)
        with patch.object(auditor, 'ROOT', self.root):
            result = auditor.verify_artifact_child(self.run, self.task, obj, self.context)
        self.assertEqual(result['artifact_sha256'], obj['artifact_sha256'])
        self.assertEqual(result['sentinel_prediction_sha256'],
                         obj['result']['reload_sentinel_prediction_sha256'])

    def test_relabelled_inner_fit_is_rejected_even_with_resigned_checkpoint(self):
        corrupted = dict(self.checkpoint)
        corrupted['fit_events'] = [dict(event) for event in corrupted['fit_events']]
        corrupted['fit_events'][0]['fit_sample_ids'] = list(
            corrupted['fit_events'][0]['fit_sample_ids'])
        corrupted['fit_events'][0]['fit_sample_ids'][0] = 999999
        save_checkpoint(self.checkpoint_path, corrupted)
        with self.assertRaisesRegex(ValueError, 'journal/checkpoint mismatch'):
            auditor.verify_checkpoint(self.checkpoint_path, self.task,
                'synthetic-frozen', self.data, self.manifest, self.grid, self.run)

    def test_journal_byte_corruption_is_rejected_before_artifact_loading(self):
        journal = self.run/self.checkpoint['fit_event_journal_path']
        with journal.open('ab') as handle:
            handle.write(b'{"torn":')
        with self.assertRaisesRegex(ValueError, 'journal path/SHA mismatch'):
            auditor.verify_checkpoint(self.checkpoint_path, self.task,
                'synthetic-frozen', self.data, self.manifest, self.grid, self.run)

    def test_artifact_byte_corruption_is_rejected_before_loading(self):
        artifact = self.run/self.checkpoint['artifact_path']
        with artifact.open('ab') as handle:
            handle.write(b'corruption')
        with self.assertRaisesRegex(ValueError, 'SHA mismatch'):
            auditor.verify_checkpoint(self.checkpoint_path, self.task,
                'synthetic-frozen', self.data, self.manifest, self.grid, self.run)

    def test_complete_626_group_manifest_validation_rejects_overlap(self):
        names = [f'SYN-{i:03d}' for i in range(626)]
        ids = pd.DataFrame(dict(canonical_model_id=names,
            biological_group_id=[f'G-{i:03d}' for i in range(626)],
            tissue=['A']*250+['B']*250+['C']*126))
        manifest, _ = make_global_inner_manifest(ids, n_folds=3, seed=2026101501)
        auditor.verify_global_manifest(ids, manifest)
        bad = manifest.copy()
        bad.loc[0, 'biological_group_id'] = bad.loc[1, 'biological_group_id']
        with self.assertRaisesRegex(ValueError, 'manifest identity metadata mismatch'):
            auditor.verify_global_manifest(ids, bad)


if __name__ == '__main__':
    unittest.main()
