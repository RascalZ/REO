"""Independent synthetic checks for the versioned ordinal diagnostic.

No original expression, drug-response or external-validation files are read.
These checks use deliberately adverse ties, tissue confounding, unknown tissue,
and changed held-out outcomes rather than mirroring the model implementation.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ordinal_diagnostic_models import (DiagnosticFeatureBuilder, DiagnosticRegressor,
                                        evaluate_diagnostic_outer_fold)
from run_ordinal_diagnostic import (make_permutations, task_list, atomic_json, checkpoint_name,
                                    validate_checkpoint, validate_prepared_checkpoint, aggregate)
from run_formal_evaluation import make_global_split_records, local_split_indices
sys.path.insert(0, str(ROOT / 'scripts'))
from validate_ordinal_diagnostic import verify_permutations, verify_folds, check_payload
import validate_ordinal_diagnostic as independent_validator


def config(genes=4):
    return {"feature": {"diagnostic_genes": genes, "min_tpm": 1.0,
                        "min_detection_fraction": 0.1},
            "elastic_net": {"regularization_mode": "alpha_max_relative",
                            "alpha_ratios": [1.0, 0.1], "l1_ratios": [0.5],
                            "max_iter": 4000, "tol": 0.001, "selection": "cyclic"}}


class OrdinalRepresentationTests(unittest.TestCase):
    def test_complete_pairs_and_gene_identified_midranks_roundtrip_with_ties(self):
        raw = np.array([[0, 0, 4, 9], [6, 2, 6, 0], [1, 4, 7, 10],
                        [2, 5, 8, 11], [3, 6, 9, 12], [4, 7, 10, 13]], dtype=np.float64)
        ids = ("D", "B", "A", "C")
        rank = DiagnosticFeatureBuilder("rank", "mad", 4, stable_gene_ids=ids).fit(raw)
        pair = DiagnosticFeatureBuilder("pair_full", "mad", 4, stable_gene_ids=ids).fit(raw)
        self.assertEqual(rank.gene_ids_, pair.gene_ids_)
        z, ranks = pair.transform(raw), rank.transform(raw)
        left, right = np.triu_indices(4, k=1)
        signed = np.zeros((len(raw), 4, 4), dtype=np.float64)
        signed[:, left, right] = z
        signed[:, right, left] = -z
        np.testing.assert_array_equal(ranks, 2.5 + 0.5 * signed.sum(axis=2))
        np.testing.assert_array_equal(z, np.sign(ranks[:, left]-ranks[:, right]))
        self.assertEqual(z.shape[1], 6)
        self.assertEqual(pair.audit_metadata_["constant_pairs_removed"], 0)

    def test_rank_uses_selected_genes_not_full_transcriptome_or_gene_standardization(self):
        raw = np.array([[1, 40, 10000], [3, 42, 10000], [9, 43, 10000],
                        [27, 41, 10000], [81, 44, 10000]], dtype=np.float64)
        model = DiagnosticFeatureBuilder("rank", "mad", 2, stable_gene_ids=("A", "B", "C")).fit(raw)
        expected = rankdata(raw[:, model.gene_indices_], axis=1, method="average")
        np.testing.assert_array_equal(model.transform(raw), expected)
        changed = raw.copy()
        changed[:, 2] = [0, 1, 2, 3, 4]
        np.testing.assert_array_equal(model.transform(changed), expected)
        self.assertEqual(model.gene_ids_, ("A", "B"))

    def test_exact_float64_pair_order_is_not_rounded_to_ties(self):
        raw = np.array([[1, np.nextafter(1., 2.), 5], [2, 3, 7], [4, 2, 8],
                        [7, 5, 9]], dtype=np.float64)
        pair = DiagnosticFeatureBuilder("pair_full", "mad", 3, stable_gene_ids=("A", "B", "C")).fit(raw)
        self.assertEqual(pair.transform(raw)[0, 0], -1)

    def test_all_rules_share_identical_ordered_genes_across_representations(self):
        rng = np.random.default_rng(919)
        raw = rng.lognormal(size=(36, 9)).astype(np.float64)
        y, tissue = rng.normal(size=36), np.repeat([0, 1, 2], 12)
        ids = tuple(f"G{i}" for i in range(8, -1, -1))
        for rule in ("mad", "response_association"):
            models = [DiagnosticFeatureBuilder(kind, rule, 4, stable_gene_ids=ids).fit(raw, y, tissue)
                      for kind in ("expression", "rank", "pair_full")]
            self.assertEqual(len({m.gene_ids_ for m in models}), 1)
            self.assertEqual(models[0].gene_ids_, tuple(sorted(models[0].gene_ids_)))
            self.assertEqual(len({m.audit_metadata_["selected_gene_ids_sha256"] for m in models}), 1)

    def test_response_association_removes_pure_tissue_confounding(self):
        # Gene T perfectly identifies tissue but has no within-tissue variation.
        tissue = np.repeat([0., 1.], 6)
        signal = np.tile(np.array([-3., -2., -1., 1., 2., 3.]), 2)
        y = 100*tissue+signal
        raw = np.column_stack((np.exp2(2+4*tissue)-1,
                               np.exp2(5+signal/4)-1,
                               np.exp2(4+np.tile([1., -1., 0., 0., -1., 1.], 2))-1)).astype(np.float64)
        model = DiagnosticFeatureBuilder("expression", "response_association", 1,
                                         stable_gene_ids=("T", "S", "N")).fit(raw, y, tissue)
        self.assertEqual(model.gene_ids_, ("S",))
        self.assertAlmostEqual(model.selected_scores_[0], 1., places=12)

    def test_response_rule_does_not_reintroduce_positive_mad_filter(self):
        raw = np.array([[0., 0.]]*8+[[4., 0.], [8., 4.]], dtype=np.float64)
        # Both genes have zero global MAD, but have within-tissue variance.
        y, tissue = np.arange(10, dtype=float), np.zeros(10)
        model = DiagnosticFeatureBuilder("expression", "response_association", 2,
                                         stable_gene_ids=("A", "B")).fit(raw, y, tissue)
        self.assertEqual(model.gene_ids_, ("A", "B"))

    def test_complete_pair_arm_retains_training_constant_and_high_tie_pairs(self):
        raw = np.column_stack((np.arange(1., 21.), np.arange(1., 21.),
                               np.arange(30., 50.))).astype(np.float64)
        model = DiagnosticFeatureBuilder("pair_full", "mad", 3,
                                         stable_gene_ids=("A", "B", "C")).fit(raw)
        self.assertEqual(model.n_output_features_, 3)
        self.assertEqual(model.audit_metadata_["training_constant_pairs_retained"], 3)
        self.assertEqual(model.audit_metadata_["training_high_tie_pairs_retained"], 1)
        np.testing.assert_array_equal(model.transform(raw)[0], [0, -1, -1])


class OrdinalPredictionBoundaryTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(311)
        self.raw = rng.lognormal(mean=2., size=(36, 6)).astype(np.float64)
        self.tissue = np.tile([0., 1., 2.], 12)
        self.x = np.column_stack((self.raw, self.tissue, np.arange(36, dtype=float)))
        self.y = 0.2*self.tissue+0.06*np.log2(self.raw[:, 0]+1)+rng.normal(0, 0.01, 36)
        self.ids = tuple(f"G{i}" for i in range(6))
        self.cfg = config()

    def test_alpha_max_zero_molecular_predictions_equal_same_policy_tissue_null(self):
        test = self.x[:5].copy()
        test[-1, -2] = 99
        null = DiagnosticRegressor("tissue_null", config=self.cfg, stable_gene_ids=self.ids).fit(self.x, self.y)
        expected = null.predict(test)
        self.assertAlmostEqual(expected[-1], self.y.mean(), places=15)
        for rule in ("mad", "response_association"):
            for kind in ("expression", "rank", "pair_full"):
                model = DiagnosticRegressor(kind, rule, self.cfg, self.ids, alpha_ratio=1.).fit(self.x, self.y)
                np.testing.assert_array_equal(model.predict(test), expected)
                self.assertEqual(model.n_nonzero_molecular_coeffs_, 0)

    def test_empty_features_keep_tissue_reference_and_constant_target_is_retained(self):
        x = self.x.copy()
        x[:, :-2] = 0.
        null = DiagnosticRegressor("tissue_null", config=self.cfg, stable_gene_ids=self.ids).fit(x, self.y)
        for kind in ("expression", "rank", "pair_full"):
            model = DiagnosticRegressor(kind, "response_association", self.cfg, self.ids, alpha_ratio=0.1).fit(x, self.y)
            self.assertEqual(model.fallback_reason_, "empty_features_same_policy_tissue_null")
            np.testing.assert_array_equal(model.predict(x), null.predict(x))
        model = DiagnosticRegressor("expression", "mad", self.cfg, self.ids).fit(self.x, np.full(36, .7))
        self.assertEqual(model.fallback_reason_, "constant_training_response")
        np.testing.assert_array_equal(model.predict(self.x), np.full(36, .7))

    def test_held_out_outcomes_cannot_change_selected_genes_hpo_or_predictions(self):
        train, test = np.arange(27), np.arange(27, 36)
        inner = [(np.flatnonzero(np.arange(27) % 3 != i), np.flatnonzero(np.arange(27) % 3 == i)) for i in range(3)]
        groups = np.arange(36)
        original = evaluate_diagnostic_outer_fold(self.x, self.y, groups, train, test, inner,
                                                   "rank", "response_association", self.cfg, self.ids)
        changed_y = self.y.copy()
        changed_y[test] = np.linspace(90., -90., len(test))
        changed = evaluate_diagnostic_outer_fold(self.x, changed_y, groups, train, test, inner,
                                                  "rank", "response_association", self.cfg, self.ids)
        np.testing.assert_array_equal(original["predictions"], changed["predictions"])
        self.assertEqual(original["chosen_params"], changed["chosen_params"])
        self.assertEqual(original["inner_scores"], changed["inner_scores"])
        self.assertEqual(original["feature_audit"], changed["feature_audit"])
        self.assertEqual(len(original["fit_events"]), 7)
        for event in original["fit_events"]:
            self.assertFalse(set(event["fit_sample_ids"]) & set(test))

    def test_nonzero_active_coefficients_have_traceable_names_and_original_scale(self):
        model = DiagnosticRegressor("expression", "response_association", self.cfg, self.ids,
                                    alpha_ratio=0.01).fit(self.x, self.y)
        audit = model.feature_audit()
        records = audit["active_molecular_coefficients"]
        self.assertGreater(len(records), 0)
        self.assertEqual(len(records), model.n_nonzero_molecular_coeffs_)
        for item in records:
            index = audit["feature_names"].index(item["feature_name"])
            coef = item["coefficient_standardized_feature_standardized_response"]
            self.assertEqual(coef, model.regressor_.coef_[index])
            self.assertAlmostEqual(item["coefficient_encoded_feature_original_response"],
                                   model.target_sd_*coef/model.scaler_.scale_[index], places=14)
            self.assertEqual(item["direction"], "positive" if coef > 0 else "negative")


class OrdinalPermutationAndCheckpointTests(unittest.TestCase):
    def setUp(self):
        self.cfg = config()
        self.cfg.update(drug_ids=[1], outer_folds=3, outer_repeats=1, inner_folds=3,
                        repeat_seeds=[3201], selection_rules=['mad', 'response_association'],
                        representations=['expression', 'rank', 'pair_full'],
                        models=['tissue_null']+[r+'_'+k for r in ['mad', 'response_association']
                                                for k in ['expression', 'rank', 'pair_full']],
                        permutation={'count': 2, 'seed': 5307})
        self.response = pd.DataFrame(dict(gdsc_drug_id=['1']*18,
            drug_name=['Synthetic compound']*18, canonical_model_id=[f'SYN-{i:02d}' for i in range(18)],
            biological_group_id=[f'GROUP-{i:02d}' for i in range(18)],
            tissue=['A']*9+['B']*9, expression_row_index=list(range(18)),
            response_auc=np.linspace(.1, .9, 18)))
        self.response.loc[0, 'response_auc'] = np.nextafter(.1, 1.)

    def test_permutation_is_reproducible_bijection_within_tissue_and_exact_label_hex(self):
        mapping, coverage = make_permutations(self.response, self.cfg)
        shuffled, _ = make_permutations(self.response.sample(frac=1, random_state=99), self.cfg)
        pd.testing.assert_frame_equal(mapping, shuffled)
        assignments, independent = verify_permutations(self.response, mapping, self.cfg)
        self.assertEqual(len(assignments), 3)
        pd.testing.assert_frame_equal(coverage[independent.columns], independent)
        zero = mapping[mapping.permutation_id.eq(0)]
        self.assertEqual(zero.target_canonical_model_id.tolist(), zero.donor_canonical_model_id.tolist())
        self.assertEqual(float.fromhex(zero.iloc[0].assigned_auc_hex), self.response.iloc[0].response_auc)
        for (_, _), block in mapping.groupby(['permutation_id', 'tissue']):
            self.assertEqual(set(block.target_canonical_model_id), set(block.donor_canonical_model_id))

    def test_independent_permutation_verifier_rejects_changed_label_or_donor(self):
        mapping, _ = make_permutations(self.response, self.cfg)
        altered = mapping.copy()
        altered.loc[0, 'assigned_auc_hex'] = float(.99).hex()
        with self.assertRaisesRegex(ValueError, 'Assigned response'):
            verify_permutations(self.response, altered, self.cfg)
        altered = mapping.copy()
        altered.loc[0, 'donor_canonical_model_id'] = 'SYN-10'
        with self.assertRaisesRegex(ValueError, 'donor'):
            verify_permutations(self.response, altered, self.cfg)

    def test_singleton_is_preserved_and_reported_not_dropped(self):
        data = self.response.copy()
        data.loc[0, 'tissue'] = 'Singleton'
        mapping, coverage = make_permutations(data, self.cfg)
        block = mapping[mapping.target_canonical_model_id.eq('SYN-00')]
        self.assertTrue(block.donor_canonical_model_id.eq('SYN-00').all())
        self.assertTrue(coverage.singleton_groups.eq(1).all())
        verify_permutations(data, mapping, self.cfg)

    def test_formal_diagnostic_cartesian_budget_is_840_tasks_and_22440_fits(self):
        cfg = json.loads((ROOT/'config/ordinal_diagnostic_v1.json').read_text(encoding='utf-8'))
        tasks = task_list(cfg)
        self.assertEqual(len(tasks), 840)
        self.assertEqual(len({checkpoint_name(t) for t in tasks}), 840)
        candidates = len(cfg['elastic_net']['alpha_ratios'])*len(cfg['elastic_net']['l1_ratios'])
        fits = sum(1 if t['model'] == 'tissue_null' else candidates*cfg['inner_folds']+1 for t in tasks)
        self.assertEqual(fits, 22440)
        self.assertEqual(sum(t['permutation_id'] == 0 for t in tasks), 84)

    def test_tissue_null_executes_once_without_inner_tuning(self):
        rng = np.random.default_rng(644)
        x = np.column_stack((rng.lognormal(size=(18, 4)), np.repeat([0., 1.], 9), np.arange(18)))
        y = np.linspace(.1, .9, 18)
        train, test = np.arange(12), np.arange(12, 18)
        inner = [(np.flatnonzero(np.arange(12) % 3 != i), np.flatnonzero(np.arange(12) % 3 == i)) for i in range(3)]
        result = evaluate_diagnostic_outer_fold(x, y, np.arange(18), train, test, inner,
                                                'tissue_null', 'none', self.cfg, ('A', 'B', 'C', 'D'))
        self.assertEqual(len(result['fit_events']), 1)
        self.assertEqual(result['inner_scores'], [])
        self.assertIsNone(result['best_inner_rmse'])

    def test_checkpoint_checksum_and_bound_prepared_labels_are_both_required(self):
        folds, _ = make_global_split_records(self.response, self.cfg)
        verify_folds(self.response, folds, self.cfg)
        task = task_list(self.cfg)[0]
        train, test, inner = local_split_indices(self.response, folds, task, self.cfg['inner_folds'])
        mapping, _ = make_permutations(self.response, self.cfg)
        lookup = {(int(d), int(p)): dict(zip(b.target_canonical_model_id, b.assigned_auc_hex))
                  for (d, p), b in mapping.groupby(['drug_id', 'permutation_id'])}
        yy = self.response.response_auc.to_numpy(dtype=np.float64)
        rows = [{**task, 'drug_name': row.drug_name, 'canonical_model_id': row.canonical_model_id,
                 'biological_group_id': row.biological_group_id, 'tissue': row.tissue,
                 'observed_auc': float(row.response_auc), 'predicted_auc': .5}
                for _, row in self.response.iloc[test].iterrows()]
        obj = dict(task=task, freeze_sha256='f'*64, status='completed', oof_rows=rows,
                   expected_test_ids=sorted(self.response.iloc[test].canonical_model_id),
                   expected_training_expression_ids=sorted(train.tolist()),
                   outer_test_expression_ids=sorted(test.tolist()),
                   inner_training_expression_ids=[sorted(train[it].tolist()) for it, iv in inner],
                   assigned_response_sha256=hashlib.sha256(yy.tobytes()).hexdigest())
        with tempfile.TemporaryDirectory(prefix='reo_ordinal_checkpoint_') as folder:
            path = Path(folder)/checkpoint_name(task)
            atomic_json(path, obj)
            self.assertTrue(validate_checkpoint(path, task, 'f'*64))
            validate_prepared_checkpoint(path, task, self.response, folds, lookup, self.cfg)
            check_payload(path)
            altered = json.loads(path.read_text(encoding='utf-8'))
            altered['oof_rows'][0]['predicted_auc'] = .2
            path.write_text(json.dumps(altered), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                validate_checkpoint(path, task, 'f'*64)
            altered = copy.deepcopy(obj)
            altered['oof_rows'][0]['observed_auc'] = .123456
            atomic_json(path, altered)
            # A valid checksum is not proof of correct labels.
            self.assertTrue(validate_checkpoint(path, task, 'f'*64))
            with self.assertRaisesRegex(ValueError, 'response/identity'):
                validate_prepared_checkpoint(path, task, self.response, folds, lookup, self.cfg)

    def test_completed_synthetic_run_passes_independent_end_to_end_validator(self):
        # Full artifact-schema integration, synthetic data only. Retain the
        # frozen cohort cardinality so the production validator's group check
        # is exercised; use four genes and one alpha to keep runtime small.
        n = 626
        rng = np.random.default_rng(451)
        cfg = copy.deepcopy(self.cfg)
        cfg['permutation']['count'] = 0
        cfg['feature']['n_genes'] = 4
        cfg['elastic_net']['alpha_ratios'] = [1.]
        tissue = np.arange(n) % 3
        raw = rng.lognormal(mean=2., size=(n, 4)).astype(np.float64)
        x = np.column_stack((raw, tissue.astype(float), np.arange(n, dtype=float)))
        y = .1*tissue+rng.uniform(.1, .5, size=n)
        response = pd.DataFrame(dict(gdsc_drug_id=['1']*n, drug_name=['Synthetic']*n,
            canonical_model_id=[f'S{i:04d}' for i in range(n)],
            biological_group_id=[f'G{i:04d}' for i in range(n)], tissue=[f'T{i}' for i in tissue],
            expression_row_index=np.arange(n), response_auc=y))
        folds, _ = make_global_split_records(response, cfg)
        mapping, coverage = make_permutations(response, cfg)
        tasks = task_list(cfg)
        with tempfile.TemporaryDirectory(prefix='reo_ordinal_e2e_') as folder:
            root = Path(folder)
            run = root/'runs'/'synthetic'
            for path in (root/'config', root/'results'/'cohort_v2', root/'data_processed', run/'checkpoints'):
                path.mkdir(parents=True, exist_ok=True)
            cfg_path = root/'config'/'ordinal_diagnostic_v1.json'
            cfg_path.write_text(json.dumps(cfg), encoding='utf-8')
            (root/'config'/'cohort_freeze_v2.lock.json').write_text(json.dumps({'freeze_sha256': 'c'*64}), encoding='utf-8')
            manifest = response.copy()
            manifest['role'] = 'development'
            manifest.to_csv(root/'results'/'cohort_v2'/'development_sample_manifest.csv', index=False)
            pd.DataFrame({'ensembl_id': ['A', 'B', 'C', 'D'], 'gene_type': ['protein_coding']*4}).to_csv(
                root/'data_processed'/'expression_genes.tsv', index=False, sep='\t')
            for name, frame in [('development_responses.csv', response), ('global_fold_manifest.csv', folds),
                                ('permutation_manifest.csv', mapping), ('permutation_coverage.csv', coverage)]:
                frame.to_csv(run/name, index=False)
            lock = {'files': {'config/ordinal_diagnostic_v1.json': independent_validator.sha(cfg_path)},
                    'v1_files': {}, 'execution_files': [],
                    'contract': {'tasks': 21, 'fits': 75, 'oof_predictions': n*7, 'drug_model_observations': n}}
            lock_path = root/'config'/'ordinal_diagnostic_v1.lock.json'
            lock_path.write_text(json.dumps(lock), encoding='utf-8')
            frozen = independent_validator.sha(lock_path)
            registered = {}
            for task in tasks:
                train, test, inner = local_split_indices(response, folds, task, cfg['inner_folds'])
                result = evaluate_diagnostic_outer_fold(x, y, response.biological_group_id.to_numpy(), train, test,
                    inner, task['representation'], task['selection_rule'], cfg, ('A', 'B', 'C', 'D'))
                pred = result.pop('predictions')
                events = result.pop('fit_events')
                rows = [{**task, 'drug_name': 'Synthetic', 'canonical_model_id': response.iloc[index].canonical_model_id,
                         'biological_group_id': response.iloc[index].biological_group_id,
                         'tissue': response.iloc[index].tissue, 'observed_auc': float(y[index]),
                         'predicted_auc': float(prediction)} for index, prediction in zip(test, pred)]
                obj = dict(task=task, freeze_sha256=frozen, status='completed', error=None, oof_rows=rows,
                    fit_events=events, result=result, expected_test_ids=sorted(response.iloc[test].canonical_model_id),
                    expected_training_expression_ids=sorted(train.tolist()), outer_test_expression_ids=sorted(test.tolist()),
                    inner_training_expression_ids=[sorted(train[it].tolist()) for it, iv in inner],
                    assigned_response_sha256=hashlib.sha256(y.tobytes()).hexdigest(),
                    resources={'wall_seconds': .1, 'cpu_seconds': .1, 'sampled_peak_rss_bytes': 1024,
                               'threadpools': []}, data_guard={'external_response_accessed': False, 'blocked_attempts': 0})
                path = run/'checkpoints'/checkpoint_name(task)
                atomic_json(path, obj)
                registered[path.name] = independent_validator.sha(path)
            counts = aggregate(run, tasks, frozen)
            prepared = ['development_responses.csv', 'global_fold_manifest.csv', 'permutation_manifest.csv', 'permutation_coverage.csv']
            outputs = ['oof_predictions.csv', 'fit_events.jsonl', 'outer_fit_audits.jsonl', 'permutation_coverage.csv']
            meta = dict(configuration=cfg, config_sha256=independent_validator.sha(cfg_path), freeze_sha256=frozen,
                cohort_sha256='c'*64, execution_commit='synthetic-no-git-source', source_hashes={},
                prepared_input_hashes={name: independent_validator.sha(run/name) for name in prepared},
                external_response_accessed=False, completed=True, status='completed',
                data_guard={'external_response_accessed': False, 'blocked_attempts': 0},
                checkpoint_hashes=registered, expected_tasks=21, completed_tasks=21,
                expected_fit_events=75, expected_oof_predictions=n*7, output_counts=counts,
                output_hashes={name: independent_validator.sha(run/name) for name in outputs}, preserved_failed_attempts=0)
            (run/'ordinal_metadata.json').write_text(json.dumps(meta), encoding='utf-8')
            with patch.object(independent_validator, 'ROOT', root):
                report = independent_validator.validate(run)
            self.assertEqual(report['status'], 'passed')
            self.assertEqual(report['counts']['actual_fits'], 75)
            self.assertEqual(report['counts']['zero_molecular_outer_fits'], 18)


class OrdinalOperationalOverheadTests(unittest.TestCase):
    """Very small synthetic file fixtures: no training or prediction occurs."""
    def fixture(self, root):
        root = Path(root)
        (root/'failed_attempts').mkdir()
        (root/'checkpoints').mkdir()
        task = dict(drug_id=1, repeat=0, outer_fold=0, model='response_association_pair_full',
                    permutation_id=0, selection_rule='response_association', representation='pair_full')
        events = []
        for candidate, inner in [(0, 0), (0, 1), (0, 2), (1, 0)]:
            events.append(dict(stage='inner_training', candidate=candidate, inner_fold=inner,
                params={'alpha_ratio': 1. if candidate == 0 else .1, 'l1_ratio': .5},
                prediction_calls=1, fit_sample_ids=[2, 3], fit_response_sha256='y'*64,
                selected_gene_ids_sha256='g'*64))
        failed = dict(task=task, freeze_sha256='f'*64, status='failed', error='Traceback\n    model.fit(x, y)\nMemoryError',
                      oof_rows=[], fit_events=events,
                      resources={'wall_seconds': 2.8, 'cpu_seconds': 2.7, 'sampled_peak_rss_bytes': 1000})
        initial = root/'failed_attempts'/'pending.json'
        atomic_json(initial, failed)
        failed_sha = independent_validator.sha(initial)
        initial.rename(root/'failed_attempts'/('failed_'+failed_sha+'.json'))
        successful = root/'checkpoints'/'old_success.json'
        successful.write_bytes(b'synthetic-success-retained')
        sessions = [dict(status='failed', started_utc='2026-10-13T00:00:00+00:00',
                         finished_utc='2026-10-13T00:00:11+00:00', wall_seconds=10.,
                         completed_tasks_this_session=6, sampled_process_tree_peak_rss_bytes=1000,
                         minimum_available_ram_bytes=2000),
                    dict(status='completed', started_utc='2026-10-13T00:00:31+00:00',
                         finished_utc='2026-10-13T00:00:42+00:00', wall_seconds=10.,
                         completed_tasks_this_session=834, sampled_process_tree_peak_rss_bytes=900,
                         minimum_available_ram_bytes=2500)]
        meta = dict(freeze_sha256='f'*64, execution_commit='synthetic-revision', configuration=config(),
                    prepared_input_hashes={}, sessions=sessions)
        snapshot = dict(metadata=meta, checkpoints={
            successful.name: {'status': 'completed', 'sha256': independent_validator.sha(successful)},
            'original_failed.json': {'status': 'failed', 'sha256': failed_sha}})
        (root/'recovery_before_fixture.json').write_text(json.dumps(snapshot), encoding='utf-8')
        audit = dict(prediction_calls=0, prediction_sample_ids=[], model_kind='pair_full',
                     selection_rule='response_association', fit_sample_ids=[2, 3],
                     fit_response_sha256='y'*64, selected_gene_ids_sha256='g'*64)
        probe = dict(status='passed', freeze_sha256='f'*64, execution_commit='synthetic-revision', task=task,
            outer_or_validation_predictions_computed=False, data_guard={'blocked_attempts': 0, 'external_response_accessed': False},
            checkpoint_sha256_before=failed_sha, checkpoint_sha256_after=failed_sha, checkpoint_unchanged=True,
            thread_environment={'OPENBLAS_NUM_THREADS': '1'}, threadpools=[{'num_threads': 1}],
            fit_audit=audit, inner_fold=1, candidate=1, n_training_samples=2,
            training_expression_ids_sha256=hashlib.sha256(np.asarray([2, 3], dtype='<i8').tobytes()).hexdigest(),
            parameters={'alpha_ratio': .1, 'l1_ratio': .5}, wall_seconds=.5,
            peak_rss_bytes=800, peak_private_bytes=600, minimum_commit_free_bytes=4000)
        probe_path = root/'probe_failed_inner_fit_fixture.json'
        probe_path.write_text(json.dumps(probe), encoding='utf-8')
        return meta, task, probe_path, successful

    def test_failed_partial_fits_and_resource_probe_are_not_successful_design_fits(self):
        with tempfile.TemporaryDirectory(prefix='reo_operational_fixture_') as folder:
            meta, task, _, _ = self.fixture(folder)
            report = independent_validator.operational_overhead(folder, meta, 22440, 100., [task])
            self.assertEqual(report['successful_pipeline_fits'], 22440)
            self.assertEqual(report['failed_attempt_completed_fit_events'], 4)
            self.assertEqual(report['resource_probe_fit_events'], 1)
            self.assertEqual(report['completed_fit_events_total'], 22445)
            self.assertEqual(report['failed_incomplete_fit_attempts'], 1)
            self.assertEqual(report['attempted_fit_total_lower_bound'], 22446)
            self.assertEqual(report['controller_session_wall_seconds_sum'], 20.)
            self.assertEqual(report['elapsed_first_start_to_final_finish_seconds'], 42.)
            self.assertEqual(report['inter_session_gap_seconds'], 20.)
            self.assertEqual(report['successful_checkpoint_hashes_preserved'], 1)
            self.assertEqual(report['failed_checkpoint_hashes_preserved'], 1)

    def test_resource_probe_cannot_include_heldout_prediction(self):
        with tempfile.TemporaryDirectory(prefix='reo_operational_prediction_') as folder:
            meta, task, path, _ = self.fixture(folder)
            obj = json.loads(path.read_text(encoding='utf-8'))
            obj['fit_audit']['prediction_calls'] = 1
            path.write_text(json.dumps(obj), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'issued a prediction'):
                independent_validator.operational_overhead(folder, meta, 22440, 100., [task])

    def test_pre_recovery_successful_checkpoint_must_remain_byte_identical(self):
        with tempfile.TemporaryDirectory(prefix='reo_operational_preservation_') as folder:
            meta, task, _, path = self.fixture(folder)
            path.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'replaced or changed'):
                independent_validator.operational_overhead(folder, meta, 22440, 100., [task])


class OrdinalStreamingAuditTests(unittest.TestCase):
    def rows(self):
        return [dict(drug_id=1, permutation_id=0, repeat=0, outer_fold=0, model='rank',
                     canonical_model_id='S0', biological_group_id='G0', tissue='A',
                     observed_auc=.1, predicted_auc=float(np.nextafter(.1, 1.))),
                dict(drug_id=1, permutation_id=0, repeat=0, outer_fold=0, model='rank',
                     canonical_model_id='S1', biological_group_id='G1', tissue='B',
                     observed_auc=.2, predicted_auc=.3)]

    def write(self, path, rows):
        with Path(path).open('w', encoding='utf-8-sig', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def test_streaming_preserves_exact_float64_roundtrip_and_entire_row_schema(self):
        rows = self.rows()
        with tempfile.TemporaryDirectory(prefix='reo_stream_fixture_') as folder:
            path = Path(folder)/'oof.csv'
            self.write(path, rows)
            stream = independent_validator.aggregate_oof_stream(path)
            independent_validator.compare_oof_block(stream, rows[:1])
            independent_validator.compare_oof_block(stream, rows[1:])
            self.assertIsNone(next(stream, None))

    def test_stream_rejects_one_ulp_prediction_change_or_context_reordering(self):
        rows = self.rows()
        examples = []
        changed = copy.deepcopy(rows)
        changed[0]['predicted_auc'] = .1
        examples.append(changed)
        examples.append(list(reversed(rows)))
        changed = copy.deepcopy(rows)
        changed[0]['drug_id'] = 2
        examples.append(changed)
        with tempfile.TemporaryDirectory(prefix='reo_stream_tamper_') as folder:
            path = Path(folder)/'oof.csv'
            for example in examples:
                self.write(path, example)
                stream = independent_validator.aggregate_oof_stream(path)
                try:
                    with self.assertRaises(ValueError):
                        independent_validator.compare_oof_block(stream, rows)
                finally:
                    stream.close()

    def test_stream_rejects_missing_rows_and_exposes_extra_rows(self):
        rows = self.rows()
        with tempfile.TemporaryDirectory(prefix='reo_stream_coverage_') as folder:
            path = Path(folder)/'oof.csv'
            self.write(path, rows[:1])
            stream = independent_validator.aggregate_oof_stream(path)
            with self.assertRaisesRegex(ValueError, 'missing row'):
                independent_validator.compare_oof_block(stream, rows)
            self.write(path, rows+[rows[0]])
            stream = independent_validator.aggregate_oof_stream(path)
            independent_validator.compare_oof_block(stream, rows)
            self.assertIsNotNone(next(stream, None))
            stream.close()

    def test_zero_prediction_hash_is_id_order_independent_and_float_bit_exact(self):
        rows = self.rows()
        digest = independent_validator.prediction_sha(rows)
        self.assertEqual(digest, independent_validator.prediction_sha(list(reversed(rows))))
        changed = copy.deepcopy(rows)
        changed[0]['predicted_auc'] = .1
        self.assertNotEqual(digest, independent_validator.prediction_sha(changed))
        changed[0]['predicted_auc'] = 0.
        positive_zero = independent_validator.prediction_sha(changed)
        changed[0]['predicted_auc'] = -0.
        self.assertNotEqual(positive_zero, independent_validator.prediction_sha(changed))
        with self.assertRaisesRegex(ValueError, 'Duplicate identity'):
            independent_validator.prediction_sha(rows+[rows[0]])


if __name__ == "__main__":
    unittest.main()
