"""Synthetic-only P6 orchestration tests; never open project CTRP responses."""
from __future__ import annotations

from contextlib import ExitStack
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import psutil

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

import run_external_validation_v1 as runner


DRUGS = (11, 12, 13, 14, 15, 16, 17, 18)
ARMS = ('tissue_null', 'expression_strong', 'mad_expression', 'mad_rank',
        'mad_pair_full', 'response_association_expression',
        'response_association_rank', 'response_association_pair_full')


def synthetic_scoring_tables():
    config = {
        'drug_ids': DRUGS, 'models': ARMS,
        'statistics': {
            'minimum_independent_common_groups_for_spearman': 3,
            'primary_comparison': {
                'candidate': 'response_association_pair_full',
                'reference': 'expression_strong', 'family_size': 9},
            'same_candidate_representation_family': {
                'contrasts': [
                    {'candidate': 'mad_pair_full', 'reference': 'mad_expression'},
                    {'candidate': 'response_association_pair_full',
                     'reference': 'response_association_expression'}],
                'family_size': 18},
            'bootstrap': {'draws': 50, 'seed': 1234,
                          'minimum_valid_fraction': .5}}}
    external, predictions, responses = [], [], []
    targets = [0.09, 0.32, 0.17, 0.52, 0.41, 0.64]
    noise = [.1, -.05, .07, -.02, .04, -.08]
    for drug in DRUGS:
        for j in range(6):
            identity = f'C{j}'
            group = f'G{j}'
            external.append({'gdsc_drug_id': drug, 'canonical_model_id': identity,
                             'biological_group_id': group})
            responses.append({'drug_id': drug, 'canonical_model_id': identity,
                              'biological_group_id': group, 'observed_auc': targets[j]})
            for arm in ARMS:
                fitted = (targets[j] if arm.endswith('pair_full') else
                          targets[j] + noise[j] if arm != 'tissue_null' else
                          .2 + .02 * j)
                predictions.append({'drug_id': drug, 'model': arm,
                    'canonical_model_id': identity, 'biological_group_id': group,
                    'predicted_auc': fitted, 'status': 'predicted',
                    'abstain_reason': '', 'unknown_training_tissue': j == 0,
                    'zero_molecular_solution': False})
    return config, pd.DataFrame(external), pd.DataFrame(predictions), pd.DataFrame(responses)


def synthetic_frozen_roster():
    """The real frozen cardinalities, with only invented IDs and no outcomes."""
    genes = [f'C{i:03d}' for i in range(175)]
    drug_panel = pd.DataFrame({'gdsc_drug_id': DRUGS,
                               'ctrp_compound_id': [1000 + i for i in range(8)],
                               'primary_drug_rank': list(range(1, 9))})
    rows = []
    next_source, next_experiment, duplicated = 2, 10000, 0
    for di, drug in enumerate(DRUGS):
        count = 149 if di == 7 else 148
        for j in range(count):
            index = (j + 20 * di) % 175
            canonical = genes[index]
            two = duplicated < 16
            duplicated += int(two)
            size = 2 if two else 1
            lines = tuple(range(next_source, next_source + size))
            experiments = tuple(range(next_experiment, next_experiment + size))
            next_source += size
            next_experiment += size
            rows.append({'gdsc_drug_id': str(drug),
                'ctrp_compound_id': str(1000 + di),
                'canonical_model_id': canonical,
                'biological_group_id': f'G{index:03d}',
                'ccle_id': f'CCLE{index:03d}',
                'expression_row_index': str(index),
                'tissue': 'T',
                'ctrp_experiment_ids': '|'.join(map(str, experiments)),
                'ctrp_source_row_numbers': '|'.join(map(str, lines)),
                'n_source_experiments': str(size),
                'identity_qc_statuses': 'SNP-matched-reference',
                'future_response_aggregation': 'median_of_included_experiments',
                'role': 'external_primary', 'cohort_version': 'v2'})
    external = pd.DataFrame(rows)
    development = pd.DataFrame({'canonical_model_id': ['D0'],
                                'biological_group_id': ['DG0']})
    split = pd.DataFrame([{
        'canonical_model_id': name, 'biological_group_id': f'G{i:03d}',
        'expression_row_index': str(i), 'tissue': 'T', 'ccle_id': f'CCLE{i:03d}',
        'global_role': 'CTRP_external_pool_never_development',
        'in_selected_external_primary': 'True'}
        for i, name in enumerate(genes)] + [{
        'canonical_model_id': 'D0', 'biological_group_id': 'DG0',
        'expression_row_index': '175', 'tissue': 'T', 'ccle_id': 'DEV',
        'global_role': 'GDSC_development_pool_never_external',
        'in_selected_external_primary': 'False'}])
    samples = pd.DataFrame({'expression_row_index': list(map(str, range(175))),
                            'depmap_id': genes,
                            'ccle_id': [f'CCLE{i:03d}' for i in range(175)]})
    return external, development, split, drug_panel, samples


class OrchestratorSyntheticTests(unittest.TestCase):
    def test_roster_counts_identity_and_source_pointer_contract(self):
        tables = synthetic_frozen_roster()
        report = runner.validate_frozen_roster(*tables, {'drug_ids': DRUGS, 'models': ARMS})
        self.assertEqual((report['eligible_drug_model_rows'], report['external_unique_groups'],
                          report['selected_source_rows']), (1185, 175, 1201))
        corrupted = tables[0].copy()
        corrupted.loc[1, 'ctrp_source_row_numbers'] = corrupted.loc[0, 'ctrp_source_row_numbers']
        with self.assertRaisesRegex(ValueError, 'Duplicate source row pointer'):
            runner.validate_frozen_roster(corrupted, *tables[1:],
                                          {'drug_ids': DRUGS, 'models': ARMS})
        header_pointer = tables[0].copy()
        header_pointer.loc[0, 'ctrp_source_row_numbers'] = '1|2'
        with self.assertRaisesRegex(ValueError, 'header offset'):
            runner.validate_frozen_roster(header_pointer, *tables[1:],
                                          {'drug_ids': DRUGS, 'models': ARMS})

    def test_source_pointer_median_and_bad_label_are_preserved(self):
        external = pd.DataFrame([
            {'gdsc_drug_id': 11, 'ctrp_compound_id': 500,
             'canonical_model_id': 'A', 'biological_group_id': 'G1',
             'ctrp_experiment_ids': '100|101', 'ctrp_source_row_numbers': '2|3',
             'n_source_experiments': 2},
            {'gdsc_drug_id': 12, 'ctrp_compound_id': 501,
             'canonical_model_id': 'B', 'biological_group_id': 'G2',
             'ctrp_experiment_ids': '102', 'ctrp_source_row_numbers': '4',
             'n_source_experiments': 1}])
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / 'invented_curves.tsv'
            source.write_text('experiment_id\tmaster_cpd_id\tarea_under_curve\n'
                '100\t500\t0.1\n101\t500\t0.3\n102\t501\tNaN\n', encoding='utf-8')
            labels, audit = runner.read_frozen_curve_responses(source, external)
            self.assertEqual(audit['selected_source_rows'], 3)
            self.assertEqual(audit['invalid_response_drug_model_rows'], 1)
            self.assertEqual(labels.canonical_model_id.tolist(), ['A'])
            self.assertAlmostEqual(labels.observed_auc.iloc[0], .2)
            source.write_text('experiment_id\tmaster_cpd_id\tarea_under_curve\n'
                '100\t500\t0.1\n999\t500\t0.3\n102\t501\tNaN\n', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'mismatches frozen pointer'):
                runner.read_frozen_curve_responses(source, external)
            source.write_text('experiment_id\tmaster_cpd_id\tarea_under_curve\n'
                '100\t500\t0.1\n101\t500\tNaN\n102\t501\t0.4\n', encoding='utf-8')
            labels, audit = runner.read_frozen_curve_responses(source, external)
            self.assertEqual(labels.canonical_model_id.tolist(), ['B'])
            self.assertEqual(audit['invalid_response_drug_model_rows'], 1)

    def test_missing_second_lock_refuses_before_response(self):
        with tempfile.TemporaryDirectory() as folder:
            missing = Path(folder) / 'second.lock.json'
            with patch.object(runner, 'SECOND_LOCK', missing):
                with self.assertRaisesRegex(PermissionError, 'absent'):
                    runner.verify_second_lock_before_response()
            with patch.object(runner, 'SECOND_LOCK', missing), \
                 patch.object(runner, 'read_frozen_curve_responses',
                              side_effect=AssertionError('raw reader reached')):
                with self.assertRaisesRegex(PermissionError, 'absent'):
                    runner.execute()

    def test_same_lock_replay_is_explicit_and_keeps_failed_attempt(self):
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder) / 'external'
            with patch.object(runner, 'EXTERNAL_RUN', run):
                start1, stem1 = runner._next_execution_attempt('a' * 64, resume=False)
                runner._write_json_exclusive(start1, {'second_lock_sha256': 'a' * 64})
                runner._write_json_exclusive(Path(str(stem1) + '_failed.json'),
                                             {'second_lock_sha256': 'a' * 64, 'error': 'synthetic failure'})
                with self.assertRaises(FileExistsError):
                    runner._next_execution_attempt('a' * 64, resume=False)
                with self.assertRaises(PermissionError):
                    runner._next_execution_attempt('b' * 64, resume=True)
                start2, _ = runner._next_execution_attempt('a' * 64, resume=True)
                self.assertEqual(start2.name, 'attempt_0002_started.json')
                self.assertTrue(Path(str(stem1) + '_failed.json').exists())
                runner._write_json_exclusive(start2, {'second_lock_sha256': 'a' * 64})
                with self.assertRaisesRegex(PermissionError, 'running or unaudited'):
                    runner._next_execution_attempt('a' * 64, resume=True)

    def test_unclean_attempt_needs_dead_process_same_lock_and_preserves_marker(self):
        with tempfile.TemporaryDirectory() as folder:
            run = Path(folder) / 'external'
            with patch.object(runner, 'EXTERNAL_RUN', run):
                start, stem = runner._next_execution_attempt('a' * 64, resume=False)
                process = psutil.Process(os.getpid())
                runner._write_json_exclusive(start, {
                    'second_lock_sha256': 'a' * 64, 'host': runner.platform.node(),
                    'pid': process.pid, 'process_create_time': process.create_time()})
                with self.assertRaisesRegex(PermissionError, 'still alive'):
                    runner._next_execution_attempt('a' * 64, resume_unclean=True)
                with self.assertRaises(PermissionError):
                    runner._next_execution_attempt('b' * 64, resume_unclean=True)
                with patch.object(runner, '_same_process_alive', return_value=False):
                    next_start, _ = runner._next_execution_attempt('a' * 64,
                                                                   resume_unclean=True)
                self.assertEqual(next_start.name, 'attempt_0002_started.json')
                abandoned = Path(str(stem) + '_unclean_abandoned.json')
                self.assertTrue(start.exists())
                self.assertTrue(abandoned.exists())
                self.assertEqual(json.loads(abandoned.read_text(encoding='utf-8'))
                                 ['original_started_sha256'], runner.sha(start))

    def test_missing_first_lock_stops_preflight_before_data(self):
        with tempfile.TemporaryDirectory() as folder:
            missing = Path(folder) / 'first.lock.json'
            with patch.object(runner, 'FIRST_LOCK', missing), \
                 patch.object(runner, 'install_pre_external_guard', return_value=(set(), [])) as guard, \
                 patch.object(runner, '_check_final_run',
                              side_effect=AssertionError('models should not load')):
                with self.assertRaises(FileNotFoundError):
                    runner.preflight()
            guard.assert_called_once()

    def test_response_free_preflight_writes_all_9480_fixed_rows(self):
        external, development, split, panel, samples = synthetic_frozen_roster()
        development = development.assign(gdsc_drug_id=str(DRUGS[0]), tissue='T')
        development = pd.concat([development.assign(gdsc_drug_id=str(drug))
                                 for drug in DRUGS], ignore_index=True)
        config = {'drug_ids': DRUGS, 'models': ARMS,
                  'protocol_status': 'frozen_pre_final_fit',
                  'pre_fit_gates_passed': True}
        ordinal = {'selection_rules': ['mad', 'response_association'],
                   'representations': ['expression', 'rank', 'pair_full'],
                   'feature': {'n_genes': 2, 'strong_expression_genes': 3}}
        stable = ('ENSG1', 'ENSG2', 'ENSG3')
        with tempfile.TemporaryDirectory() as folder, ExitStack() as patches:
            root = Path(folder)
            paths = {
                'ROOT': root, 'EXTERNAL_RUN': root/'runs/external',
                'FIRST_LOCK': root/'config/first.lock.json',
                'EXTERNAL_CONFIG': root/'config/external.json',
                'ORDINAL_CONFIG': root/'config/ordinal.json',
                'EXTERNAL_MANIFEST': root/'results/external.csv',
                'DEVELOPMENT_MANIFEST': root/'results/development.csv',
                'GLOBAL_SPLIT': root/'results/split.csv',
                'MAIN_DRUGS': root/'results/drugs.csv',
                'EXPRESSION_GENES': root/'data/genes.tsv',
                'EXPRESSION_SAMPLES': root/'data/samples.tsv',
                'EXPRESSION_TPM': root/'data/tpm.npy',
                'FINAL_RUN': root/'runs/final'}
            for name, path in paths.items():
                patches.enter_context(patch.object(runner, name, path))
            for path in paths.values():
                if isinstance(path, Path):
                    path.parent.mkdir(parents=True, exist_ok=True)
            (paths['FINAL_RUN']/'models').mkdir(parents=True)
            (paths['FINAL_RUN']/'checkpoints').mkdir()
            external.to_csv(paths['EXTERNAL_MANIFEST'], index=False)
            development.to_csv(paths['DEVELOPMENT_MANIFEST'], index=False)
            split.to_csv(paths['GLOBAL_SPLIT'], index=False)
            panel.to_csv(paths['MAIN_DRUGS'], index=False)
            samples.to_csv(paths['EXPRESSION_SAMPLES'], sep='\t', index=False)
            pd.DataFrame({'expression_column_index': range(3),
                          'ensembl_id': stable,
                          'gene_type': ['protein_coding'] * 3}).to_csv(
                              paths['EXPRESSION_GENES'], sep='\t', index=False)
            np.save(paths['EXPRESSION_TPM'], np.full((175, 3), 2., dtype=np.float64))
            paths['EXTERNAL_CONFIG'].write_text(json.dumps(config), encoding='utf-8')
            paths['ORDINAL_CONFIG'].write_text(json.dumps(ordinal), encoding='utf-8')
            paths['FIRST_LOCK'].write_text('{}', encoding='utf-8')
            (paths['FINAL_RUN']/'worker_context.json').write_text(
                json.dumps({'tissue_code': {'T': 0}}), encoding='utf-8')
            for name in ('final_model_metadata.json', 'model_manifest.jsonl',
                         'final_model_independent_audit.json'):
                (paths['FINAL_RUN']/name).write_text('{}', encoding='utf-8')
            first_files = {}
            for name in ('src/run_external_validation_v1.py', 'src/external_validation_v1.py'):
                source = root/name
                source.parent.mkdir(exist_ok=True)
                source.write_text('synthetic frozen source', encoding='utf-8')
                first_files[name] = runner.sha(source)
            first = {'execution_files': list(first_files), 'files': first_files}
            first_sha = runner.sha(paths['FIRST_LOCK'])
            expression_genes_sha = runner.sha(paths['EXPRESSION_GENES'])
            expression_tpm_sha = runner.sha(paths['EXPRESSION_TPM'])
            records, model_hashes, artifacts = [], {}, {}
            for drug in DRUGS:
                for arm in ARMS:
                    task = {'drug_id': drug, 'model': arm,
                            **runner.model_spec(arm, ordinal)}
                    name = f'models/d{drug}_{arm}.joblib'
                    records.append({'task': task, 'artifact_path': name,
                                    'artifact_sha256': 'a' * 64})
                    model_hashes[(paths['FINAL_RUN']/name).relative_to(root).as_posix()] = 'a' * 64
                    artifacts[(paths['FINAL_RUN']/name).resolve()] = {
                        'task': task, 'freeze_sha256': first_sha,
                        'stable_gene_ids': stable, 'tissue_code': {'T': 0},
                        'expression_genes_sha256': expression_genes_sha,
                        'expression_tpm_sha256': expression_tpm_sha,
                        'model': object()}
            patches.enter_context(patch.object(runner, 'verify_freeze',
                                               return_value=(first, 'synthetic_commit', {})))
            patches.enter_context(patch.object(runner, '_check_final_run',
                                               return_value=({}, records, model_hashes)))
            patches.enter_context(patch.object(runner, 'read_checkpoint',
                                               return_value={'result': {'feature_audit': None}}))
            patches.enter_context(patch.object(runner.joblib, 'load',
                                               side_effect=lambda path: artifacts[Path(path).resolve()]))
            def fake_predict(_model, batch, *, drug_id, arm, representation,
                             expected_selected_gene_ids, expected_feature_names,
                             artifact_sha256):
                return pd.DataFrame({
                    'drug_id': drug_id, 'model': arm,
                    'canonical_model_id': batch.identity.canonical_model_id,
                    'biological_group_id': batch.identity.biological_group_id,
                    'tissue': batch.identity.tissue,
                    'eligible_denominator_member': True,
                    'predicted_auc': .25, 'status': 'predicted',
                    'abstain_reason': '', 'unknown_training_tissue': False,
                    'zero_molecular_solution': False,
                    'artifact_sha256': artifact_sha256})
            patches.enter_context(patch.object(runner, 'predict_external_arm',
                                               side_effect=fake_predict))
            # The real development guard remains installed here. It denies any
            # project CTRP raw/index read while this synthetic preflight runs.
            meta = runner.preflight()
            output = pd.read_csv(paths['EXTERNAL_RUN']/'prediction_manifest.csv')
            self.assertEqual((meta['prediction_rows'], len(output)), (9480, 9480))
            self.assertEqual(output.groupby(['drug_id', 'model']).size().shape[0], 64)
            self.assertEqual(set(output.columns), set(runner.PREDICTION_FIELDS))
            self.assertNotIn('observed_auc', output.columns)
            self.assertNotIn('area_under_curve', output.columns)
            self.assertFalse(meta['external_response_accessed'])
            self.assertFalse(meta['external_response_index_accessed'])
            self.assertGreaterEqual(meta['resource_usage']['wall_seconds'], 0)
            self.assertGreater(meta['resource_usage']['process_rss_bytes'], 0)

    def test_synthetic_execution_records_resources_and_does_not_repeat_success(self):
        with tempfile.TemporaryDirectory() as folder, ExitStack() as patches:
            root = Path(folder)
            run = root/'runs/external'
            run.mkdir(parents=True)
            raw_name = 'data_raw/ctrp/invented_curves.tsv'
            raw = root/raw_name
            raw.parent.mkdir(parents=True)
            raw.write_text('experiment_id\tmaster_cpd_id\tarea_under_curve\n'
                           '100\t500\t0.2\n', encoding='utf-8')
            external_path = root/'results/external.csv'
            external_path.parent.mkdir()
            pd.DataFrame([{'gdsc_drug_id': 11, 'ctrp_compound_id': 500,
                           'canonical_model_id': 'A', 'biological_group_id': 'G1',
                           'ctrp_experiment_ids': '100', 'ctrp_source_row_numbers': '2',
                           'n_source_experiments': 1}]).to_csv(external_path, index=False)
            prediction_path = run/'prediction_manifest.csv'
            prediction_path.write_text('drug_id,model,canonical_model_id,biological_group_id,predicted_auc,status\n'
                                       '11,tissue_null,A,G1,0.3,predicted\n', encoding='utf-8')
            preflight = run/'preflight_metadata.json'
            preflight.write_text(json.dumps({'identity_audit': {'synthetic': True},
                                             'prediction_coverage': []}), encoding='utf-8')
            second_path = root/'config/second.lock.json'
            second_path.parent.mkdir()
            second_path.write_text('{}', encoding='utf-8')
            config_path = root/'config/external.json'
            config_path.write_text('{}', encoding='utf-8')
            second = {'response_source_path': raw_name,
                      'response_source_sha256': runner.sha(raw),
                      'response_column': 'area_under_curve',
                      'prediction_manifest_sha256': runner.sha(prediction_path),
                      'preflight_metadata_sha256': runner.sha(preflight),
                      'status': 'frozen_pre_external_response'}
            for name, value in {
                'ROOT': root, 'EXTERNAL_RUN': run, 'SECOND_LOCK': second_path,
                'RAW_RESPONSE_NAME': raw_name, 'EXTERNAL_MANIFEST': external_path,
                'EXTERNAL_CONFIG': config_path}.items():
                patches.enter_context(patch.object(runner, name, value))
            patches.enter_context(patch.object(runner, 'verify_second_lock_before_response',
                                               return_value=second))
            patches.enter_context(patch.object(runner, 'score_locked_panel',
                                               return_value={'scope': 'synthetic_only'}))
            result = runner.execute()
            self.assertEqual(result['status'], 'completed_one_time_external_evaluation')
            self.assertGreaterEqual(result['resource_usage']['wall_seconds'], 0)
            self.assertGreater(result['resource_usage']['process_rss_bytes'], 0)
            self.assertTrue((run/'execution_attempts/attempt_0001_completed.json').exists())
            self.assertTrue((run/'evaluation_summary.json').exists())
            with patch.object(runner, 'read_frozen_curve_responses',
                              side_effect=AssertionError('raw response reader reached again')):
                with self.assertRaises(FileExistsError):
                    runner.execute()

    def test_scoring_fixed_eight_shared_groups_and_strict_macro(self):
        config, external, prediction, response = synthetic_scoring_tables()
        result = runner.score_locked_panel(config, external, prediction, response,
                                           expected_global_groups=6)
        self.assertEqual(set(result['contrasts']), {
            'primary', 'representation_mad', 'representation_response_association'})
        primary = result['contrasts']['primary']['inference']
        self.assertEqual(primary['bootstrap_seed'], 1234)
        self.assertEqual(primary['shared_global_groups'], 6)
        self.assertIsNotNone(primary['strict_macro']['delta_rho'])
        self.assertEqual(result['contrasts']['primary']['family_size'], 9)
        self.assertEqual(result['contrasts']['representation_mad']['family_size'], 18)
        self.assertEqual(len(result['paired_descriptive']), 88)
        self.assertTrue(all(item['descriptive_only'] for item in result['paired_descriptive']))
        missing_one_drug = response.loc[~response.drug_id.eq(18)].copy()
        failed = runner.score_locked_panel(config, external, prediction,
                                           missing_one_drug, expected_global_groups=6)
        self.assertIsNone(failed['contrasts']['primary']['inference']['strict_macro']['delta_rho'])
        self.assertEqual(failed['contrasts']['primary']['coverage']['18']['n_response_valid_ids'], 0)
        self.assertEqual(failed['contrasts']['primary']['inference']['per_drug'][18]['undefined_reason'],
                         'insufficient_independent_groups;insufficient_independent_groups')


if __name__ == '__main__':
    unittest.main()
