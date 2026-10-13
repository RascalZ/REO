"""Synthetic-only independent contracts for the additive eight-drug run.

No private expression/response data, old run checkpoint or CTRP file is read.
Small synthetic fits establish counterfactual training boundaries and strong
gene budgets; full artifact validation uses a temporary synthetic cohort.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import importlib
import io
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
from ordinal_repeated_models import evaluate_repeated_outer_fold, model_spec
from run_formal_evaluation import make_global_split_records, local_split_indices, atomic_json
import validate_ordinal_repeated as audit
import run_ordinal_repeated_stage as launcher
import ordinal_repeated_progress as progress


def config():
    cfg = json.loads((ROOT/'config/ordinal_repeated_v1.json').read_text(encoding='utf-8'))
    cfg['feature'].update(n_genes=4, strong_expression_genes=8)
    cfg['elastic_net'].update(alpha_ratios=[1., .03], l1_ratios=[.5])
    return cfg


def synthetic_problem(n=36, genes=12):
    rng = np.random.default_rng(81119)
    raw = rng.lognormal(mean=2, size=(n, genes)).astype(np.float64)
    tissue = np.arange(n) % 3
    x = np.column_stack((raw, tissue, np.arange(n, dtype=float)))
    y = .08*tissue + .03*np.log2(raw[:, 0]+1) + rng.uniform(.01, .05, n)
    train, test = np.arange(n-9), np.arange(n-9, n)
    inner = [(np.flatnonzero(np.arange(len(train)) % 3 != i),
              np.flatnonzero(np.arange(len(train)) % 3 == i)) for i in range(3)]
    ids = tuple(f'G{i:05d}' for i in range(genes-1, -1, -1))
    return x, y, np.arange(n), train, test, inner, ids


class RepeatedModelContracts(unittest.TestCase):
    def test_fixed_grid_has_960_tasks_26160_fits_and_observed_labels_only(self):
        cfg = json.loads((ROOT/'config/ordinal_repeated_v1.json').read_text(encoding='utf-8'))
        tasks = audit.tasks_for(cfg)
        self.assertEqual(len(tasks), 960)
        self.assertEqual(len({audit.filename(t) for t in tasks}), 960)
        self.assertEqual({t['permutation_id'] for t in tasks}, {0})
        self.assertEqual(cfg['drug_ids'], [1931, 1036, 1191, 1930, 1199, 1091, 1010, 1922])
        self.assertEqual(cfg['repeat_seeds'], [20261013, 20261014, 20261015])
        candidates = len(cfg['elastic_net']['alpha_ratios'])*len(cfg['elastic_net']['l1_ratios'])
        self.assertEqual(sum(1 if t['model'] == 'tissue_null' else candidates*cfg['inner_folds']+1 for t in tasks), 26160)
        self.assertEqual(sum([625, 622, 620, 620, 620, 620, 620, 620])*3*8, 119208)
        for task in tasks:
            self.assertEqual({k: task[k] for k in ('selection_rule', 'representation', 'n_genes')},
                             model_spec(task['model'], cfg))

    def test_runner_cartesian_tasks_and_names_match_independent_contract(self):
        runner = importlib.import_module('run_ordinal_repeated')
        cfg = json.loads((ROOT/'config/ordinal_repeated_v1.json').read_text(encoding='utf-8'))
        tasks = audit.tasks_for(cfg)
        self.assertEqual(runner.task_list(cfg), tasks)
        self.assertEqual([runner.checkpoint_name(t) for t in tasks], [audit.filename(t) for t in tasks])

    def test_strong_reference_actually_builds_2000_genes_and_matched_is_subset(self):
        problem = synthetic_problem(n=30, genes=2011)
        cfg = config()
        cfg['feature'].update(n_genes=200, strong_expression_genes=2000)
        cfg['elastic_net']['alpha_ratios'] = [1.]
        original = copy.deepcopy(cfg)
        results = {}
        for arm in ('expression_strong', 'mad_expression'):
            results[arm] = evaluate_repeated_outer_fold(*problem[:6], arm, cfg, problem[6])
        strong = results['expression_strong']['feature_audit']
        matched = results['mad_expression']['feature_audit']
        self.assertEqual(len(strong['gene_ids']), 2000)
        self.assertEqual(len(matched['gene_ids']), 200)
        self.assertTrue(set(matched['gene_ids']) <= set(strong['gene_ids']))
        self.assertEqual(strong['gene_ids'], sorted(strong['gene_ids']))
        self.assertEqual(cfg, original, 'Adapter must not mutate shared/frozen configuration')
        self.assertEqual({event['candidate_gene_budget'] for event in results['expression_strong']['fit_events']}, {2000})

    def test_all_arms_retain_train_only_states_when_outer_labels_or_expression_change(self):
        x, y, groups, train, test, inner, genes = synthetic_problem()
        cfg = config()
        altered_y = y.copy()
        altered_y[test] = np.linspace(80., -80., len(test))
        altered_x = x.copy()
        altered_x[test, :-2] = altered_x[test, :-2][:, ::-1]*1000
        for model in cfg['models']:
            baseline = evaluate_repeated_outer_fold(x, y, groups, train, test, inner, model, cfg, genes)
            label = evaluate_repeated_outer_fold(x, altered_y, groups, train, test, inner, model, cfg, genes)
            expression = evaluate_repeated_outer_fold(altered_x, y, groups, train, test, inner, model, cfg, genes)
            np.testing.assert_array_equal(label['predictions'], baseline['predictions'])
            for changed in (label, expression):
                self.assertEqual(changed['chosen_params'], baseline['chosen_params'])
                self.assertEqual(changed['inner_scores'], baseline['inner_scores'])
                self.assertEqual(changed['feature_audit'], baseline['feature_audit'])
                for event in changed['fit_events']:
                    self.assertFalse(set(event['fit_sample_ids']) & set(test))
                    self.assertFalse(set(event['fit_sample_ids']) & set(event['prediction_sample_ids']))

    def test_shared_gene_hashes_and_exact_zero_predictions_match_tissue_null(self):
        problem = synthetic_problem()
        cfg = config()
        cfg['elastic_net']['alpha_ratios'] = [1.]
        results = {m: evaluate_repeated_outer_fold(*problem[:6], m, cfg, problem[6]) for m in cfg['models']}
        reference = results['tissue_null']['predictions']
        for name, result in results.items():
            np.testing.assert_array_equal(result['predictions'], reference)
            self.assertEqual(result['outer_fit_audit']['n_nonzero_molecular_coeffs'], 0)
        for rule in cfg['selection_rules']:
            result = [results[rule+'_'+rep] for rep in cfg['representations']]
            for events in zip(*(r['fit_events'] for r in result)):
                self.assertEqual(len({e['selected_gene_ids_sha256'] for e in events}), 1)
        self.assertEqual(len(results['tissue_null']['fit_events']), 1)
        self.assertEqual(results['tissue_null']['inner_scores'], [])

    def test_source_adapter_leaves_historical_numerical_sources_unchanged(self):
        files = [ROOT/'src/ordinal_diagnostic_models.py', ROOT/'src/modeling/nested_cv.py',
                 ROOT/'config/ordinal_diagnostic_v1.json', ROOT/'config/formal_model_v1.json']
        before = {p: p.read_bytes() for p in files}
        p = synthetic_problem()
        evaluate_repeated_outer_fold(*p[:6], 'expression_strong', config(), p[6])
        self.assertEqual(before, {p: p.read_bytes() for p in files})


class RepeatedIndependentAuditContracts(unittest.TestCase):
    def test_checkpoint_checksum_does_not_replace_prepared_identity_and_label_validation(self):
        runner = importlib.import_module('run_ordinal_repeated')
        cfg = config()
        cfg.update(drug_ids=[1], outer_folds=3, outer_repeats=1, repeat_seeds=[99])
        n = 30
        y = np.linspace(.01, .9, n)
        y[0] = np.nextafter(.1, 1.)
        data = pd.DataFrame(dict(gdsc_drug_id=['1']*n, drug_name=['Synthetic']*n,
            canonical_model_id=[f'S{i:03}' for i in range(n)], biological_group_id=[f'B{i:03}' for i in range(n)],
            tissue=['T'+str(i%3) for i in range(n)], expression_row_index=np.arange(n), response_auc=y))
        folds, _ = make_global_split_records(data, cfg)
        task = audit.tasks_for(cfg)[0]
        train, test, inner = local_split_indices(data, folds, task, 3)
        rows = [{**task, 'drug_name': 'Synthetic', 'canonical_model_id': data.iloc[i].canonical_model_id,
                 'biological_group_id': data.iloc[i].biological_group_id, 'tissue': data.iloc[i].tissue,
                 'observed_auc': float(y[i]), 'predicted_auc': .5} for i in test]
        obj = dict(task=task, freeze_sha256='f'*64, status='completed', oof_rows=rows,
                   expected_test_ids=sorted(data.iloc[test].canonical_model_id),
                   expected_training_expression_ids=sorted(train.tolist()), outer_test_expression_ids=sorted(test.tolist()),
                   inner_training_expression_ids=[sorted(train[it].tolist()) for it, iv in inner],
                   assigned_response_sha256=audit.y_sha(y))
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/audit.filename(task)
            atomic_json(path, obj)
            self.assertTrue(runner.validate_checkpoint(path, task, 'f'*64))
            runner.verify_prepared_checkpoint(path, task, data, folds, cfg)
            corrupted = json.loads(path.read_text(encoding='utf-8'))
            corrupted['oof_rows'][0]['observed_auc'] = .987654
            path.write_text(json.dumps(corrupted), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                runner.validate_checkpoint(path, task, 'f'*64)
            altered = copy.deepcopy(obj)
            altered['oof_rows'][0]['observed_auc'] = .987654
            atomic_json(path, altered)
            self.assertTrue(runner.validate_checkpoint(path, task, 'f'*64))
            with self.assertRaisesRegex(ValueError, 'outcome/identity'):
                runner.verify_prepared_checkpoint(path, task, data, folds, cfg)

    def test_worker_context_binds_exact_configuration_sources_and_prepared_scope(self):
        meta = dict(configuration={'model': 'test'}, freeze_sha256='f'*64, execution_commit='synthetic',
                    source_hashes={'src/synthetic.py': 'a'*64}, prepared_input_hashes={'development_responses.csv': 'b'*64})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'worker_context.json'
            context = dict(configuration=meta['configuration'], freeze_sha256=meta['freeze_sha256'],
                           execution_commit=meta['execution_commit'], source_hashes=meta['source_hashes'],
                           prepared_input_hashes=dict(meta['prepared_input_hashes']))
            path.write_text(json.dumps(context), encoding='utf-8')
            meta['prepared_input_hashes']['worker_context.json'] = audit.sha(path)
            audit.verify_worker_context(folder, meta)
            context['source_hashes'] = {'src/synthetic.py': 'c'*64}
            path.write_text(json.dumps(context), encoding='utf-8')
            meta['prepared_input_hashes']['worker_context.json'] = audit.sha(path)
            with self.assertRaisesRegex(ValueError, 'source or lock'):
                audit.verify_worker_context(folder, meta)

    def test_fresh_child_rechecks_prepared_and_scientific_bytes_before_loading_model(self):
        runner = importlib.import_module('run_ordinal_repeated')
        with tempfile.TemporaryDirectory() as folder:
            root, run = Path(folder), Path(folder)/'run'
            run.mkdir()
            source = root/'synthetic_source.py'
            source.write_text('VALUE = 1\n', encoding='utf-8')
            prepared = run/'development_responses.csv'
            prepared.write_text('response_auc\n0.10000000000000002\n', encoding='utf-8')
            context = dict(configuration={}, freeze_sha256='f'*64, execution_commit='synthetic',
                           prepared_input_hashes={prepared.name: audit.sha(prepared)},
                           source_hashes={source.name: audit.sha(source)})
            path = run/'worker_context.json'
            path.write_text(json.dumps(context), encoding='utf-8')
            digest = audit.sha(path)
            with patch.object(runner, 'ROOT', root):
                self.assertEqual(runner.load_worker_context(run, digest), context)
                source.write_text('VALUE = 2\n', encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'scientific source changed'):
                    runner.load_worker_context(run, digest)
                source.write_text('VALUE = 1\n', encoding='utf-8')
                prepared.write_text('response_auc\n0.1\n', encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'prepared input changed'):
                    runner.load_worker_context(run, digest)

    def test_worker_loads_nontrivial_response_decimals_with_float64_roundtrip(self):
        runner = importlib.import_module('run_ordinal_repeated')
        cfg = config()
        cfg.update(drug_ids=[1], outer_folds=3, outer_repeats=1, repeat_seeds=[17])
        x, y, groups, _, _, _, genes = synthetic_problem(n=30)
        y[0] = np.nextafter(.1, 1.)
        data = pd.DataFrame(dict(gdsc_drug_id=['1']*len(y), drug_name=['Synthetic']*len(y),
            canonical_model_id=[f'M{i:03}' for i in range(len(y))], biological_group_id=[f'G{i:03}' for i in range(len(y))],
            tissue=['T'+str(int(v)) for v in x[:, -2]], expression_row_index=np.arange(len(y)), response_auc=y))
        folds, _ = make_global_split_records(data, cfg)
        task = audit.tasks_for(cfg)[0]
        state = dict(access=[], blocked=[], stable=genes, folds=folds, config=cfg)
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder)
            (run/'checkpoints').mkdir()
            data.to_csv(run/'development_responses.csv', index=False)
            original_default = pd.read_csv(run/'development_responses.csv').response_auc.to_numpy()
            self.assertNotEqual(audit.y_sha(original_default), audit.y_sha(y),
                                'Synthetic response must expose default-parser rounding')
            state['responses'] = pd.read_csv(run/'development_responses.csv')
            with patch.object(runner.base, 'WORKER', state), patch.object(runner.base, 'worker_init'), \
                 patch.object(runner.base, 'get_drug', side_effect=lambda drug: (state['responses'], x)):
                success = runner.execute_task(task, run, cfg, 'f'*64,
                                              dict(execution_commit='synthetic', source_hashes={}))
            self.assertTrue(success)
            checkpoint = audit.check_payload(run/'checkpoints'/audit.filename(task))
            self.assertEqual(checkpoint['assigned_response_sha256'], audit.y_sha(y))
            runner.verify_prepared_checkpoint(run/'checkpoints'/audit.filename(task), task, data, folds, cfg)

    def test_fresh_process_identity_allows_pid_recycling_but_not_process_reuse(self):
        seen = set()
        resource = dict(pid=21, controller_pid=11, process_create_time=1234.,
                        task_process_isolated=True, preimport_threads_applied=True)
        audit.verify_worker_isolation(resource, seen)
        with self.assertRaisesRegex(ValueError, 'reused'):
            audit.verify_worker_isolation(resource, seen)
        audit.verify_worker_isolation(dict(resource, process_create_time=1235.), seen)
        with self.assertRaisesRegex(ValueError, 'pre-import'):
            audit.verify_worker_isolation(dict(resource, pid=23, preimport_threads_applied=False), seen)
        with self.assertRaisesRegex(ValueError, 'controller'):
            audit.verify_worker_isolation(dict(resource, pid=11), seen)

    def test_mad_budget_comparison_rejects_missing_or_non_nested_genes(self):
        audit.check_mad_budget_sets({'fold': {200: frozenset(['A']), 2000: frozenset(['A', 'B'])}})
        with self.assertRaisesRegex(ValueError, 'contained'):
            audit.check_mad_budget_sets({'fold': {200: frozenset(['C']), 2000: frozenset(['A', 'B'])}})
        with self.assertRaisesRegex(ValueError, 'Missing'):
            audit.check_mad_budget_sets({'fold': {200: frozenset(['A'])}})

    def test_streaming_checks_exact_float_bits_typed_context_and_fit_payload(self):
        row = dict(drug_id=1, repeat=0, model='x', observed_auc=.1, predicted_auc=-0.)
        encoded = {k: str(v) for k, v in row.items()}
        audit.compare_oof_block(iter([encoded]), [row])
        for key, value in [('repeat', '1'), ('predicted_auc', '0.0'), ('observed_auc', str(np.nextafter(.1, 1.)))]:
            with self.assertRaises(ValueError):
                audit.compare_oof_block(iter([{**encoded, key: value}]), [row])
        with self.assertRaisesRegex(ValueError, 'typed fields'):
            audit.compare_json_record(iter([dict(row, repeat=0.)]), row)

    def test_global_repeat_and_inner_fold_annotations_are_independently_reconstructed(self):
        cfg = config()
        n = 60
        data = pd.DataFrame(dict(canonical_model_id=[f'M{i:03}' for i in range(n)],
                                 biological_group_id=[f'G{i:03}' for i in range(n)], tissue=['T'+str(i%3) for i in range(n)]))
        folds, _ = make_global_split_records(data, cfg)
        audit.verify_folds(data, folds, cfg)
        changed = folds.copy()
        changed.loc[0, 'biological_group_id'] = 'WRONG'
        with self.assertRaisesRegex(ValueError, 'annotation'):
            audit.verify_folds(data, changed, cfg)
        with self.assertRaisesRegex(ValueError, 'membership'):
            audit.verify_folds(data, folds.iloc[1:], cfg)

    def test_safe_pause_stops_launcher_before_audit_and_keeps_thread_environment(self):
        calls = []
        def run(command, **kwargs):
            calls.append((command, kwargs))
            return type('Result', (), {'returncode': 75})()
        with tempfile.TemporaryDirectory() as folder, patch.object(launcher.subprocess, 'run', side_effect=run), patch('sys.stdout', new=io.StringIO()):
            code = launcher.execute_stages([('fit', ['fit.py']), ('validate', ['qa.py'])], 'new-scientific-sha', folder)
            self.assertEqual(code, 75)
            self.assertEqual(len(calls), 1)
            self.assertEqual(len(list(Path(folder).glob('*.log'))), 1)
            self.assertTrue(all(calls[0][1]['env'][name] == '1' for name in launcher.THREAD_VARIABLES))

    def test_progress_safe_pause_does_not_count_missing_tasks_as_failure(self):
        events = [dict(event='ordinal_repeated_start', total_tasks=960, pending_tasks=959),
                  dict(event='ordinal_repeated_task_finished', status='completed', task={'id': 2}, wall_seconds=.2),
                  dict(event='resource_paused', reason='commit headroom')]
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'fit.log'
            path.write_text('\n'.join(json.dumps(e) for e in events)+'\n{partial', encoding='utf-8')
            value = progress.progress_from_log(path, 960)
        self.assertEqual(value['completed_tasks'], 2)
        self.assertEqual(value['failed_tasks_this_session'], 0)
        self.assertEqual(value['status'], 'resource_paused')

    def test_paused_session_is_retained_as_overhead_without_inventing_failure_fits(self):
        meta = dict(freeze_sha256='f'*64, execution_commit='synthetic', prepared_input_hashes={}, sessions=[
            dict(status='resource_paused', started_utc='2026-10-13T01:00:00+00:00',
                 finished_utc='2026-10-13T01:00:10+00:00', wall_seconds=10., completed_tasks_this_session=0),
            dict(status='completed', started_utc='2026-10-13T01:00:30+00:00',
                 finished_utc='2026-10-13T01:00:40+00:00', wall_seconds=10., completed_tasks_this_session=1)])
        with tempfile.TemporaryDirectory() as folder:
            value = audit.operational_overhead(folder, meta, 31, 9., [])
        self.assertEqual(value['resource_paused_session_count'], 1)
        self.assertEqual(value['failed_session_count'], 0)
        self.assertEqual(value['completed_fit_events_total'], 31)
        self.assertEqual(value['inter_session_gap_seconds'], 20.)

    def test_failed_attempt_and_pre_recovery_success_hashes_are_preserved(self):
        task = audit.tasks_for(config())[0]
        meta = dict(freeze_sha256='f'*64, execution_commit='synthetic', prepared_input_hashes={}, sessions=[])
        failed = dict(task=task, freeze_sha256='f'*64, status='failed', error='Traceback: model.fit( failed',
                      oof_rows=[], fit_events=[{'prediction_calls': 1}]*4, resources={'wall_seconds': 2.})
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder)
            (run/'failed_attempts').mkdir(); (run/'checkpoints').mkdir()
            path = run/'failed_attempts'/'temporary.json'
            atomic_json(path, failed)
            failed_sha = audit.sha(path)
            path.rename(path.with_name('failed_'+failed_sha+'.json'))
            success = run/'checkpoints'/'success.json'
            success.write_text('{"status":"completed"}', encoding='utf-8')
            snapshot = dict(metadata=meta, checkpoints={
                'success.json': dict(status='completed', sha256=audit.sha(success)),
                'failed.json': dict(status='failed', sha256=failed_sha)})
            (run/'recovery_before_synthetic.json').write_text(json.dumps(snapshot), encoding='utf-8')
            result = audit.operational_overhead(run, meta, 31, 1., [task])
            self.assertEqual(result['completed_fit_events_total'], 35)
            self.assertEqual(result['attempted_fit_total_lower_bound'], 36)
            self.assertEqual(result['successful_checkpoint_hashes_preserved'], 1)
            self.assertEqual(result['failed_checkpoint_hashes_preserved'], 1)
            success.write_text('{"status":"changed"}', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'Successful checkpoint changed'):
                audit.operational_overhead(run, meta, 31, 1., [task])

    def test_synthetic_complete_artifacts_pass_streaming_audit_and_edited_oof_fails(self):
        cfg = config()
        cfg.update(drug_ids=[1], outer_folds=3, outer_repeats=1, repeat_seeds=[8121])
        cfg['elastic_net']['alpha_ratios'] = [1.]
        n, p = 626, 8
        rng = np.random.default_rng(4141)
        tissue = np.arange(n) % 3
        x = np.column_stack((rng.lognormal(size=(n, p)), tissue, np.arange(n)))
        y = .1*tissue+rng.uniform(.1, .5, n)
        genes = tuple(f'G{i}' for i in range(p))
        response = pd.DataFrame(dict(gdsc_drug_id=['1']*n, drug_name=['Synthetic']*n,
            canonical_model_id=[f'M{i:04}' for i in range(n)], biological_group_id=[f'B{i:04}' for i in range(n)],
            tissue=['T'+str(t) for t in tissue], expression_row_index=np.arange(n), response_auc=y))
        folds, _ = make_global_split_records(response, cfg)
        tasks = audit.tasks_for(cfg)
        with tempfile.TemporaryDirectory(prefix='reo_repeated_synthetic_') as folder:
            root, run = Path(folder), Path(folder)/'run'
            for directory in (root/'config', root/'results/cohort_v2', root/'data_processed', run/'checkpoints'):
                directory.mkdir(parents=True, exist_ok=True)
            config_path = root/'config/ordinal_repeated_v1.json'
            config_path.write_text(json.dumps(cfg), encoding='utf-8')
            (root/'config/cohort_freeze_v2.lock.json').write_text(json.dumps({'freeze_sha256': 'c'*64}), encoding='utf-8')
            response.assign(role='development').to_csv(root/'results/cohort_v2/development_sample_manifest.csv', index=False)
            pd.DataFrame(dict(ensembl_id=genes, gene_type=['protein_coding']*p)).to_csv(root/'data_processed/expression_genes.tsv', sep='\t', index=False)
            response.to_csv(run/'development_responses.csv', index=False)
            folds.to_csv(run/'global_fold_manifest.csv', index=False)
            lock = dict(files={'config/ordinal_repeated_v1.json': audit.sha(config_path)}, v1_files={}, diagnostic_files={}, execution_files=[],
                        contract=dict(tasks=24, fits=87, oof_predictions=n*8, drug_model_observations=n))
            lock_path = root/'config/ordinal_repeated_v1.lock.json'
            lock_path.write_text(json.dumps(lock), encoding='utf-8')
            frozen, registered = audit.sha(lock_path), {}
            with (run/'oof_predictions.csv').open('w', newline='', encoding='utf-8') as fo, (run/'fit_events.jsonl').open('w', encoding='utf-8') as ff, (run/'outer_fit_audits.jsonl').open('w', encoding='utf-8') as fa:
                writer = None
                for number, task in enumerate(tasks):
                    train, test, inner = local_split_indices(response, folds, task, 3)
                    result = evaluate_repeated_outer_fold(x, y, response.biological_group_id.to_numpy(), train, test, inner,
                                                         task['model'], cfg, genes)
                    predictions, events = result.pop('predictions'), result.pop('fit_events')
                    rows = [{**task, 'drug_name': 'Synthetic', 'canonical_model_id': response.iloc[i].canonical_model_id,
                             'biological_group_id': response.iloc[i].biological_group_id, 'tissue': response.iloc[i].tissue,
                             'observed_auc': float(y[i]), 'predicted_auc': float(pred)} for i, pred in zip(test, predictions)]
                    resources = dict(wall_seconds=.1, cpu_seconds=.1, sampled_peak_rss_bytes=1024,
                                     threadpools=[dict(num_threads=1)], pid=100+number, controller_pid=99,
                                     process_create_time=1000.+number, task_process_isolated=True, preimport_threads_applied=True,
                                     initialization_wall_seconds=.02, fit_work_wall_seconds=.07)
                    obj = dict(task=task, freeze_sha256=frozen, status='completed', error=None, oof_rows=rows,
                               execution_commit='synthetic', source_hashes={},
                               fit_events=events, result=result, resources=resources,
                               expected_test_ids=sorted(response.iloc[test].canonical_model_id),
                               expected_training_expression_ids=sorted(train.tolist()), outer_test_expression_ids=sorted(test.tolist()),
                               inner_training_expression_ids=[sorted(train[it].tolist()) for it, iv in inner],
                               assigned_response_sha256=audit.y_sha(y), data_guard={'external_response_accessed': False, 'blocked_attempts': 0})
                    path = run/'checkpoints'/audit.filename(task)
                    atomic_json(path, obj)
                    registered[path.name] = audit.sha(path)
                    if writer is None:
                        writer = csv.DictWriter(fo, fieldnames=list(rows[0])); writer.writeheader()
                    writer.writerows(rows)
                    for event in events:
                        ff.write(json.dumps({**task, **event})+'\n')
                    fa.write(json.dumps(dict(task=task, **result, resources=resources, data_guard=obj['data_guard']))+'\n')
            outputs = ['oof_predictions.csv', 'fit_events.jsonl', 'outer_fit_audits.jsonl']
            meta = dict(configuration=cfg, config_sha256=audit.sha(config_path), freeze_sha256=frozen,
                        cohort_sha256='c'*64, execution_commit='synthetic', source_hashes={}, external_response_accessed=False,
                        prepared_input_hashes={name: audit.sha(run/name) for name in ['development_responses.csv', 'global_fold_manifest.csv']},
                        completed=True, status='completed', data_guard={'external_response_accessed': False, 'blocked_attempts': 0},
                        checkpoint_hashes=registered, expected_tasks=24, completed_tasks=24, expected_fit_events=87,
                        expected_oof_predictions=n*8, output_counts={'tasks': 24, 'fit_events': 87, 'oof_predictions': n*8, 'failed': 0},
                        output_hashes={name: audit.sha(run/name) for name in outputs}, preserved_failed_attempts=0)
            context = dict(configuration=cfg, freeze_sha256=frozen, execution_commit='synthetic', source_hashes={},
                           prepared_input_hashes=dict(meta['prepared_input_hashes']))
            (run/'worker_context.json').write_text(json.dumps(context), encoding='utf-8')
            meta['prepared_input_hashes']['worker_context.json'] = audit.sha(run/'worker_context.json')
            (run/'ordinal_repeated_metadata.json').write_text(json.dumps(meta), encoding='utf-8')
            with patch.object(audit, 'ROOT', root):
                report = audit.validate(run)
                self.assertEqual(report['status'], 'passed')
                self.assertEqual(report['counts']['actual_fits'], 87)
                self.assertEqual(report['fresh_worker_processes_verified'], 24)
                text = (run/'oof_predictions.csv').read_text(encoding='utf-8')
                altered = text.splitlines()
                values = altered[1].split(',')
                values[-1] = str(np.nextafter(float(values[-1]), 1.))
                altered[1] = ','.join(values)
                (run/'oof_predictions.csv').write_text('\n'.join(altered)+'\n', encoding='utf-8')
                with self.assertRaisesRegex(ValueError, 'float64'):
                    audit.validate(run)


if __name__ == '__main__':
    unittest.main()
