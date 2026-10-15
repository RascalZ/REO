"""Synthetic final-fit checks; no project response or CTRP file is opened."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

from final_model_v1 import (fit_full_development_model, local_inner_splits,
                            make_global_inner_manifest, prediction_sha256)
from external_validation_v1 import names_sha256, validate_model_contract
from run_final_model_v1 import (make_reload_sentinels, read_checkpoint, save_checkpoint,
                                task_list, execute_fit_task, execute_verify_task, _aggregate,
                                _archive_stale_temp, validate_prepared_checkpoint,
                                fit_journal_path,
                                prepare_resource_samples,
                                recover_interrupted_task_if_needed)
from ordinal_repeated_runtime import sha
import run_final_model_v1 as final_runner


def small_configuration():
    return dict(drug_ids=[1], models=['tissue_null', 'response_association_pair_full'],
        selection_rules=['mad', 'response_association'],
        representations=['expression', 'rank', 'pair_full'], inner_folds=3,
        feature=dict(n_genes=5, strong_expression_genes=6, min_tpm=1.,
                     min_detection_fraction=.1),
        elastic_net=dict(regularization_mode='alpha_max_relative',
            alpha_ratios=[1., .3, .1, .03, .01], l1_ratios=[.5, 1.],
            max_iter=4000, tol=.001, selection='cyclic'))


class FinalModelSyntheticTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(12345)
        n = 30
        self.ids = pd.DataFrame(dict(canonical_model_id=[f'SYN-{j:03d}' for j in range(n)],
            biological_group_id=[f'SYN-GROUP-{j:03d}' for j in range(n)],
            tissue=np.tile(['A', 'B', 'C'], 10)))
        self.manifest, self.warnings = make_global_inner_manifest(self.ids, n_folds=3, seed=2026101501)
        self.stable = tuple(f'ENSG_SYN_{j:02d}' for j in range(6))
        raw = rng.lognormal(mean=1., sigma=.7, size=(n, 6)).astype(np.float64)
        tissue_codes = self.ids.tissue.map({'A': 0, 'B': 1, 'C': 2}).to_numpy()
        self.X = np.column_stack((raw, tissue_codes, np.arange(n))).astype(np.float64)
        self.y = (.35 + .04*np.log2(raw[:, 0]+1) - .02*np.log2(raw[:, 2]+1)
                  + .02*tissue_codes + rng.normal(0, .01, n))
        self.cfg = small_configuration()

    def test_fixed_real_panel_has_64_tasks_and_1744_fits(self):
        self.assertEqual(final_runner.DEFAULT_LOCK.name, 'external_validation_v1_1.lock.json')
        self.assertEqual(final_runner.DEFAULT_EXTERNAL_CONFIG.name, 'external_validation_v1_1.json')
        self.assertTrue(final_runner.DEFAULT_RUN.name.endswith('-final-v1_1'))
        ordinal = json.loads((ROOT/'config/ordinal_repeated_v1.json').read_text(encoding='utf-8'))
        external = json.loads(final_runner.DEFAULT_EXTERNAL_CONFIG.read_text(encoding='utf-8'))
        tasks = task_list(external, ordinal)
        self.assertEqual(len(tasks), 64)
        self.assertEqual(len({(t['drug_id'], t['model']) for t in tasks}), 64)
        self.assertEqual(sum(1 if t['model'] == 'tissue_null' else 31 for t in tasks), 1744)

    def test_one_global_grouped_split_and_per_drug_intersection(self):
        self.assertEqual(len(self.manifest), 90)
        self.assertEqual(len(self.manifest.loc[self.manifest.role.eq('validation')]), 30)
        data = self.ids.iloc[np.flatnonzero(np.arange(30) % 5 != 0)].reset_index(drop=True)
        inner = local_inner_splits(data, self.manifest, 3)
        self.assertEqual(sum(len(iv) for _, iv in inner), len(data))
        for it, iv in inner:
            self.assertFalse(set(data.iloc[it].biological_group_id) &
                             set(data.iloc[iv].biological_group_id))
            self.assertEqual(set(it) | set(iv), set(range(len(data))))

    def test_pair_full_training_has_30_inner_fits_and_one_full_refit(self):
        inner = local_inner_splits(self.ids, self.manifest, 3)
        model, result = fit_full_development_model(self.X, self.y,
            self.ids.biological_group_id.to_numpy(), inner,
            'response_association_pair_full', self.cfg, self.stable)
        self.assertEqual(result['fit_count'], 31)
        self.assertEqual(len(result['inner_scores']), 10)
        self.assertEqual(result['chosen_params'], min(result['inner_scores'],
            key=lambda row: row['pooled_inner_rmse'])['params'])
        for event in result['fit_events']:
            if event['stage'] == 'inner_training':
                it, iv = inner[event['inner_fold']]
                self.assertEqual(set(event['fit_sample_ids']), set(self.X[it, -1].astype(int)))
                self.assertFalse(set(event['fit_sample_ids']) & set(self.X[iv, -1].astype(int)))
            else:
                self.assertEqual(event['stage'], 'full_development_refit')
                self.assertEqual(set(event['fit_sample_ids']), set(range(30)))
        self.assertEqual(result['feature_audit']['output_features'], 10)
        self.assertEqual(result['development_prediction_sha256'],
                         prediction_sha256(model.predict(self.X)))

    def test_tissue_null_has_one_fit_and_unknown_tissue_falls_back_to_global_mean(self):
        inner = local_inner_splits(self.ids, self.manifest, 3)
        model, result = fit_full_development_model(self.X, self.y,
            self.ids.biological_group_id.to_numpy(), inner,
            'tissue_null', self.cfg, self.stable)
        self.assertEqual(result['fit_count'], 1)
        self.assertEqual(result['inner_scores'], [])
        self.assertIsNone(result['best_inner_rmse'])
        sentinels = make_reload_sentinels(self.X, {'A': 0, 'B': 1, 'C': 2})
        self.assertEqual(sentinels[1, -2], 3.)
        self.assertTrue(np.all(sentinels[2, :-2] == 0.))
        self.assertAlmostEqual(float(model.predict(sentinels)[1]), float(self.y.mean()), places=14)

    def test_duplicate_or_overlapping_groups_are_rejected(self):
        bad = self.ids.copy()
        bad.loc[1, 'biological_group_id'] = bad.loc[0, 'biological_group_id']
        with self.assertRaisesRegex(ValueError, 'one identity per biological group'):
            make_global_inner_manifest(bad, n_folds=3, seed=2026101501)
        inner = local_inner_splits(self.ids, self.manifest, 3)
        bad_inner = list(inner)
        it, iv = bad_inner[0]
        bad_inner[0] = (np.r_[it, iv[0]], iv)
        with self.assertRaisesRegex(ValueError, 'partition'):
            fit_full_development_model(self.X, self.y,
                self.ids.biological_group_id.to_numpy(), bad_inner,
                'tissue_null', self.cfg, self.stable)

    def test_trusted_joblib_reload_in_fresh_python_process_matches_sentinels(self):
        inner = local_inner_splits(self.ids, self.manifest, 3)
        model, _ = fit_full_development_model(self.X, self.y,
            self.ids.biological_group_id.to_numpy(), inner,
            'tissue_null', self.cfg, self.stable)
        sentinels = make_reload_sentinels(self.X, {'A': 0, 'B': 1, 'C': 2})
        expected = model.predict(sentinels)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            joblib.dump(model, path/'model.joblib', compress=3)
            np.save(path/'sentinels.npy', sentinels)
            code = ("import joblib,numpy as np,sys;"
                    "sys.path.insert(0,sys.argv[1]);"
                    "m=joblib.load(sys.argv[2]);"
                    "x=np.load(sys.argv[3]);"
                    "np.save(sys.argv[4],m.predict(x))")
            subprocess.run([sys.executable, '-c', code, str(ROOT/'src'),
                str(path/'model.joblib'), str(path/'sentinels.npy'), str(path/'pred.npy')],
                check=True, capture_output=True, text=True)
            np.testing.assert_allclose(np.load(path/'pred.npy'), expected, rtol=0, atol=1e-12)

    def test_checkpoint_contract_detects_model_corruption(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            (out/'checkpoints').mkdir(); (out/'models').mkdir()
            artifact = out/'models/d1_tissue_null.joblib'
            artifact.write_bytes(b'synthetic-only-artifact')
            task = dict(drug_id=1, model='tissue_null', selection_rule='none',
                        representation='tissue_null', n_genes=0)
            journal = fit_journal_path(out, task)
            journal.parent.mkdir()
            journal.write_text('{"stage":"full_development_refit"}\n', encoding='utf-8')
            obj = dict(schema_version='1.0', task=task, freeze_sha256='frozen',
                status='completed', fit_count=1, fit_events=[{'stage': 'full_development_refit'}],
                artifact_path='models/d1_tissue_null.joblib', artifact_sha256=sha(artifact),
                fit_event_journal_path='fit_event_journals/'+journal.name,
                fit_event_journal_sha256=sha(journal), fit_event_journal_parse_error=None,
                fresh_process_reload_verified=True)
            checkpoint = out/'checkpoints/d1_tissue_null.json'
            save_checkpoint(checkpoint, obj)
            self.assertEqual(read_checkpoint(checkpoint, task, 'frozen')['status'], 'completed')
            artifact.write_bytes(b'corrupted')
            with self.assertRaisesRegex(ValueError, 'artifact SHA mismatch'):
                read_checkpoint(checkpoint, task, 'frozen')

    def test_synthetic_worker_checkpoint_and_reload_preserve_complete_model(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            (out/'checkpoints').mkdir(); (out/'models').mkdir()
            self.manifest.to_csv(out/'global_inner_fold_manifest.csv', index=False)
            data = self.ids.copy()
            data['response_auc'] = self.y
            data['gdsc_drug_id'] = '1'
            data['expression_row_index'] = np.arange(len(data))
            context = dict(freeze_sha256='synthetic-frozen', ordinal_config={**self.cfg,
                    'resources': {'rss_sampling_interval_seconds': .01}},
                external_config={'final_training': {'inner_folds': 3}},
                tissue_code={'A': 0, 'B': 1, 'C': 2},
                development_file_hashes={
                    'data_processed/expression_genes.tsv': 'synthetic-gene-hash',
                    'data_processed/expression_tpm.npy': 'synthetic-tpm-hash'})
            task = dict(drug_id=1, model='response_association_pair_full',
                        selection_rule='response_association', representation='pair_full', n_genes=5)
            with patch('run_final_model_v1.load_drug', return_value=(data, self.X, self.y, self.stable)):
                fitted = execute_fit_task(out, task, context)
                checkpoint = out/'checkpoints/d1_response_association_pair_full.json'
                self.assertTrue(fitted, json.loads(checkpoint.read_text(encoding='utf-8')).get('error'))
                pending = read_checkpoint(checkpoint, task, 'synthetic-frozen')
                self.assertEqual(pending['status'], 'fitted_pending_verification')
                self.assertEqual(pending['fit_count'], 31)
                journal = fit_journal_path(out, task)
                self.assertEqual(pending['fit_event_journal_sha256'], sha(journal))
                self.assertEqual(len(journal.read_text(encoding='utf-8').splitlines()), 31)
                self.assertTrue(execute_verify_task(out, task, context))
            completed = read_checkpoint(checkpoint, task, 'synthetic-frozen')
            self.assertEqual(completed['status'], 'completed')
            self.assertEqual(completed['reload_verification']['development_n'], 30)
            self.assertEqual(completed['reload_verification']['sentinel_policy'],
                             'known tissue, unseen tissue, all-zero exact-tie sample')
            self.assertTrue(completed['fresh_process_reload_verified'])
            self.assertTrue(completed['artifact_sha256'])
            artifact = joblib.load(out/completed['artifact_path'])
            trained = artifact['model']
            validate_model_contract(trained, representation='pair_full',
                stable_gene_ids=artifact['stable_gene_ids'],
                training_tissue_codes=trained.training_tissue_codes_,
                expected_gene_sha256=names_sha256(self.stable),
                expected_selected_gene_ids=completed['result']['feature_audit']['gene_ids'],
                expected_feature_sha256=completed['result']['full_fit_audit']['built_feature_names_sha256'])
            validate_prepared_checkpoint(completed, task, data, self.manifest, 3)
            changed = data.copy()
            changed.loc[0, 'response_auc'] += .001
            with self.assertRaisesRegex(ValueError, 'response differs'):
                validate_prepared_checkpoint(completed, task, changed, self.manifest, 3)

    def test_failed_synthetic_fit_is_preserved_as_failed_checkpoint(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            (out/'checkpoints').mkdir(); (out/'models').mkdir()
            context = dict(freeze_sha256='synthetic-frozen', ordinal_config={**self.cfg,
                    'resources': {'rss_sampling_interval_seconds': .01}})
            task = dict(drug_id=1, model='tissue_null', selection_rule='none',
                        representation='tissue_null', n_genes=0)
            with patch('run_final_model_v1.load_drug', side_effect=RuntimeError('synthetic failure')):
                self.assertFalse(execute_fit_task(out, task, context))
            failure = read_checkpoint(out/'checkpoints/d1_tissue_null.json', task,
                                      'synthetic-frozen', verify_artifact=False)
            self.assertEqual(failure['status'], 'failed')
            self.assertIn('synthetic failure', failure['error'])
            self.assertEqual(failure['fit_count'], 0)
            self.assertEqual(failure['fit_event_journal_sha256'],
                             sha(fit_journal_path(out, task)))

    def test_mid_fit_failure_preserves_fsynced_events_and_refuses_implicit_refit(self):
        from ordinal_diagnostic_models import DiagnosticRegressor
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            (out/'checkpoints').mkdir(); (out/'models').mkdir()
            self.manifest.to_csv(out/'global_inner_fold_manifest.csv', index=False)
            data = self.ids.copy()
            data['response_auc'] = self.y
            data['gdsc_drug_id'] = '1'
            data['expression_row_index'] = np.arange(len(data))
            context = dict(freeze_sha256='synthetic-frozen',
                ordinal_config={**self.cfg, 'resources': {'rss_sampling_interval_seconds': .01}},
                external_config={'final_training': {'inner_folds': 3}},
                tissue_code={'A': 0, 'B': 1, 'C': 2})
            task = dict(drug_id=1, model='response_association_pair_full',
                selection_rule='response_association', representation='pair_full', n_genes=5)
            original_fit = DiagnosticRegressor.fit
            calls = []
            def fail_on_fifth(model, X, y):
                calls.append(1)
                if len(calls) == 5:
                    raise RuntimeError('injected failure before fifth fit')
                return original_fit(model, X, y)
            with patch('run_final_model_v1.load_drug', return_value=(data, self.X, self.y, self.stable)), \
                    patch.object(DiagnosticRegressor, 'fit', fail_on_fifth):
                self.assertFalse(execute_fit_task(out, task, context))
            path = out/'checkpoints/d1_response_association_pair_full.json'
            failed = read_checkpoint(path, task, 'synthetic-frozen')
            journal = fit_journal_path(out, task)
            journal_before = journal.read_bytes()
            self.assertEqual(failed['status'], 'failed')
            self.assertEqual(failed['fit_count'], 4)
            self.assertEqual(len(failed['fit_events']), 4)
            self.assertEqual(len(journal_before.splitlines()), 4)
            self.assertEqual(failed['fit_event_journal_sha256'], sha(journal))
            self.assertIn('injected failure before fifth fit', failed['error'])
            self.assertEqual([e['candidate'] for e in failed['fit_events']], [0, 0, 0, 1])
            with self.assertRaisesRegex(ValueError, 'Refusing to overwrite'):
                execute_fit_task(out, task, context)
            self.assertEqual(journal.read_bytes(), journal_before)
            journal.write_bytes(journal_before+b'{"torn":')
            with self.assertRaisesRegex(ValueError, 'journal path/SHA mismatch'):
                read_checkpoint(path, task, 'synthetic-frozen')

    def test_crash_without_checkpoint_recovers_partial_journal_as_failed(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            (out/'checkpoints').mkdir(); (out/'models').mkdir()
            task = dict(drug_id=1, model='tissue_null', selection_rule='none',
                        representation='tissue_null', n_genes=0)
            journal = fit_journal_path(out, task)
            journal.parent.mkdir()
            journal.write_text('{"stage":"full_development_refit"}\n', encoding='utf-8')
            partial = journal.read_bytes()
            artifact = out/'models/d1_tissue_null.joblib'
            artifact.write_bytes(b'preserved-orphan')
            with self.assertRaisesRegex(ValueError, 'preserved as failed'):
                recover_interrupted_task_if_needed(out, task, 'synthetic-frozen')
            obj = read_checkpoint(out/'checkpoints/d1_tissue_null.json',
                task, 'synthetic-frozen')
            self.assertEqual(obj['fit_count'], 1)
            self.assertEqual(obj['fit_events'], [{'stage': 'full_development_refit'}])
            self.assertEqual(obj['artifact_sha256'], sha(artifact))
            self.assertEqual(journal.read_bytes(), partial)
            self.assertEqual(obj['status'], 'failed')

    def test_worker_freeze_rejects_context_mutation_after_lock_verification(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root/'data_processed').mkdir()
            genes = root/'data_processed/expression_genes.tsv'
            tpm = root/'data_processed/expression_tpm.npy'
            genes.write_bytes(b'synthetic gene list')
            tpm.write_bytes(b'synthetic tpm matrix')
            lock = root/'external.lock.json'; lock.write_text('{}', encoding='utf-8')
            external = root/'external.json'; external.write_text('{"a":1}', encoding='utf-8')
            ordinal = root/'ordinal.json'; ordinal.write_text('{"b":2}', encoding='utf-8')
            environment = {'python': 'synthetic', 'packages': {}}
            context = dict(freeze_sha256=sha(lock), execution_commit='synthetic-commit',
                source_hashes={'src/synthetic.py': 'synthetic-hash'}, environment=environment,
                external_config={'a': 1}, ordinal_config={'b': 2},
                development_file_hashes={
                    'data_processed/expression_genes.tsv': sha(genes),
                    'data_processed/expression_tpm.npy': sha(tpm)},
                prepared_input_hashes={'development_responses.csv': 'synthetic-response-sha'})
            patchers = [patch.object(final_runner, 'ROOT', root),
                patch.object(final_runner, 'DEFAULT_LOCK', lock),
                patch.object(final_runner, 'DEFAULT_EXTERNAL_CONFIG', external),
                patch.object(final_runner, 'DEFAULT_ORDINAL_CONFIG', ordinal),
                patch.object(final_runner, 'verify_freeze', return_value=(
                    {'environment': environment}, 'synthetic-commit',
                    {'src/synthetic.py': 'synthetic-hash'})),
                patch.object(final_runner, '_response_source_verified',
                    return_value=(root/'response.csv', 'synthetic-response-sha', 'meta-sha'))]
            for patcher in patchers:
                patcher.start()
                self.addCleanup(patcher.stop)
            final_runner.verify_worker_freeze(context)
            changed = {**context, 'freeze_sha256': 'different-lock'}
            with self.assertRaisesRegex(ValueError, 'context differs'):
                final_runner.verify_worker_freeze(changed)
            changed = {**context, 'external_config': {'a': 999}}
            with self.assertRaisesRegex(ValueError, 'config differs'):
                final_runner.verify_worker_freeze(changed)

    def test_direct_worker_cli_rechecks_first_lock_before_any_fit(self):
        with tempfile.TemporaryDirectory() as folder:
            for flag in ('--worker-task', '--verify-task'):
                with self.subTest(flag=flag), \
                        patch.object(sys, 'argv', ['run_final_model_v1.py', '--run-dir', folder,
                                                  flag, '0', '--context-sha', 'synthetic']), \
                        patch.object(final_runner, '_checked_run_dir', return_value=Path(folder)), \
                        patch.object(final_runner, 'install_pre_external_guard', return_value=(set(), [])), \
                        patch.object(final_runner, 'load_context', return_value={'freeze_sha256': 'old'}), \
                        patch.object(final_runner, 'verify_worker_freeze',
                                     side_effect=ValueError('pre-fit lock rejected')) as gate, \
                        patch.object(final_runner, 'execute_fit_task') as fit, \
                        patch.object(final_runner, 'execute_verify_task') as verify:
                    with self.assertRaisesRegex(ValueError, 'pre-fit lock rejected'):
                        final_runner.main()
                    gate.assert_called_once()
                    fit.assert_not_called()
                    verify.assert_not_called()

    def test_narrow_prefit_guard_allows_ctrp_metadata_but_blocks_response_and_index(self):
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            (base/'data_raw/ctrp').mkdir(parents=True)
            (base/'data_processed').mkdir()
            (base/'data_raw/ctrp/v20._COLUMNS.txt').write_text('metadata', encoding='utf-8')
            (base/'data_raw/ctrp/v20.data.curves_post_qc.txt').write_text('response', encoding='utf-8')
            (base/'data_processed/ctrp_curve_index.csv').write_text('index', encoding='utf-8')
            code = '''import json,sys
from pathlib import Path
sys.path.insert(0,sys.argv[2])
import freeze_external_validation_v1_1 as freeze
freeze.ROOT=Path(sys.argv[1])
access,blocked=freeze.install_pre_external_guard()
base=freeze.ROOT
metadata=(base/'data_raw/ctrp/v20._COLUMNS.txt').read_text(encoding='utf-8')
result={}
for key,path in [('response',base/'data_raw/ctrp/v20.data.curves_post_qc.txt'),
                 ('index',base/'data_processed/ctrp_curve_index.csv')]:
    try:
        path.read_text(encoding='utf-8')
        result[key]='allowed'
    except PermissionError:
        result[key]='blocked'
print(json.dumps(dict(metadata=metadata,result=result,access=sorted(access),blocked=blocked)))'''
            child = subprocess.run([sys.executable, '-c', code, str(base), str(ROOT/'scripts')],
                cwd=ROOT, check=True, capture_output=True, text=True)
            outcome = json.loads(child.stdout.strip())
            self.assertEqual(outcome['metadata'], 'metadata')
            self.assertEqual(outcome['result'], {'response': 'blocked', 'index': 'blocked'})
            self.assertEqual(len(outcome['blocked']), 2)
            self.assertIn(str(base/'data_raw/ctrp/v20._COLUMNS.txt'), outcome['access'])

    def test_cli_rejects_custom_protocol_paths_before_worker_or_controller(self):
        with tempfile.TemporaryDirectory() as folder:
            for flag in ('--config', '--ordinal-config', '--lock'):
                for action in ((), ('--worker-task', '0'), ('--verify-task', '0')):
                    with self.subTest(flag=flag, action=action), \
                            patch.object(sys, 'argv', ['run_final_model_v1.py',
                                '--run-dir', folder, flag, str(Path(folder)/'custom.json'),
                                *action]), \
                            patch.object(final_runner, '_checked_run_dir') as checked_run, \
                            patch.object(final_runner, 'install_pre_external_guard') as guard, \
                            patch.object(final_runner, 'verify_freeze') as freeze, \
                            patch.object(final_runner, 'execute_fit_task') as fit, \
                            patch.object(final_runner, 'execute_verify_task') as verify:
                        with self.assertRaisesRegex(ValueError, 'fixed pre-fit protocol path'):
                            final_runner.main()
                        checked_run.assert_not_called()
                        guard.assert_not_called()
                        freeze.assert_not_called()
                        fit.assert_not_called()
                        verify.assert_not_called()

    def test_resource_sample_name_collision_archives_prior_partial(self):
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            prior = out/'resource_samples_session2.csv'
            prior.write_bytes(b'old partial resource sample\n')
            target = prepare_resource_samples(out, 2)
            self.assertEqual(target, prior)
            self.assertFalse(target.exists())
            archived = list((out/'interrupted_partial_files').glob(
                'resource_samples_session2.csv.*'))
            self.assertEqual(len(archived), 1)
            self.assertEqual(archived[0].read_bytes(), b'old partial resource sample\n')
            with target.open('x', encoding='utf-8') as handle:
                handle.write('new session\n')
            self.assertEqual(target.read_text(encoding='utf-8'), 'new session\n')

    def test_checkpoint_only_aggregate_is_atomic_and_preserves_stale_partial(self):
        ordinal = json.loads((ROOT/'config/ordinal_repeated_v1.json').read_text(encoding='utf-8'))
        external = json.loads(final_runner.DEFAULT_EXTERNAL_CONFIG.read_text(encoding='utf-8'))
        tasks = task_list(external, ordinal)
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            (out/'checkpoints').mkdir(); (out/'models').mkdir()
            for task in tasks:
                stem = f"d{task['drug_id']}_{task['model']}"
                artifact = out/'models'/(stem+'.joblib')
                artifact.write_bytes(stem.encode('ascii'))
                nfit = 1 if task['model'] == 'tissue_null' else 31
                obj = dict(schema_version='1.0', task=task, freeze_sha256='frozen',
                    status='completed', fit_count=nfit, fit_events=[{}]*nfit,
                    artifact_path='models/'+artifact.name, artifact_sha256=sha(artifact),
                    fit_event_journal_path='fit_event_journals/'+stem+'.jsonl',
                    fresh_process_reload_verified=True,
                    result={'chosen_params': {}, 'best_inner_rmse': None,
                        'full_fit_audit': {'n_nonzero_molecular_coeffs': 0,
                            'fitted_kind': 'synthetic', 'fallback_reason': '',
                            'selected_gene_ids_sha256': None},
                        'development_prediction_sha256': 'synthetic-prediction-hash'})
                journal = fit_journal_path(out, task)
                journal.parent.mkdir(exist_ok=True)
                journal.write_text('{}\n'*nfit, encoding='utf-8')
                obj['fit_event_journal_sha256'] = sha(journal)
                obj['fit_event_journal_parse_error'] = None
                save_checkpoint(out/'checkpoints'/(stem+'.json'), obj)
            (out/'model_manifest.jsonl.tmp').write_text('old partial', encoding='utf-8')
            summary = _aggregate(out, tasks, 'frozen')
            self.assertEqual(summary['artifacts'], 64)
            self.assertEqual(summary['fits'], 1744)
            self.assertEqual(len((out/'model_manifest.jsonl').read_text(encoding='utf-8').splitlines()), 64)
            archived = list((out/'interrupted_partial_files').glob('model_manifest.jsonl.tmp.*'))
            self.assertEqual(len(archived), 1)
            self.assertEqual(archived[0].read_text(encoding='utf-8'), 'old partial')
            orphan = out/'models/orphan.joblib'
            orphan.write_bytes(b'incomplete-artifact')
            _archive_stale_temp(orphan, out)
            self.assertFalse(orphan.exists())
            self.assertEqual(len(list((out/'interrupted_partial_files').glob('orphan.joblib.*'))), 1)


if __name__ == '__main__':
    unittest.main()
