"""Synthetic public aggregate snapshots; standard library, no private run access."""
from __future__ import annotations

import csv
import importlib.util
import itertools
import json
import math
from pathlib import Path
import statistics
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('verify_public_ordinal_repeated', ROOT/'scripts/verify_public_ordinal_repeated.py')
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


def write_csv(path, rows):
    with Path(path).open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows({k: '' if isinstance(v, float) and math.isnan(v) else v for k, v in row.items()} for row in rows)


def repeat_summary(values):
    finite = all(map(math.isfinite, values))
    return dict(repeat_values=json.dumps([v if math.isfinite(v) else None for v in values], allow_nan=False),
                repeat_mean=sum(values)/3, repeat_min=min(values) if finite else math.nan,
                repeat_max=max(values) if finite else math.nan,
                repeat_range=max(values)-min(values) if finite else math.nan,
                repeat_sample_sd=statistics.stdev(values) if finite else math.nan,
                positive_repeats=sum(v > 0 for v in values), negative_repeats=sum(v < 0 for v in values),
                zero_repeats=sum(v == 0 for v in values), undefined_repeats=sum(not math.isfinite(v) for v in values),
                scope=VERIFY.REPEAT_SCOPE)


def interval(point, family_size=None):
    defined = math.isfinite(point)
    return dict(point=point, bootstrap_requested=10000, bootstrap_valid=10000 if defined else 0,
                bootstrap_invalid=0 if defined else 10000, bootstrap_valid_fraction=1. if defined else 0.,
                conditional_95_lower=point-.01 if defined else '', conditional_95_upper=point+.01 if defined else '',
                family_size='' if family_size is None else family_size,
                family_adjusted_conditional_lower=point-.02 if defined and family_size is not None else '',
                family_adjusted_conditional_upper=point+.02 if defined and family_size is not None else '',
                ci_status='reported' if defined else 'withheld_undefined_point',
                undefined_reason='' if defined else 'undefined_constituent_repeat_or_drug', interval_scope=VERIFY.INTERVAL_SCOPE)


def make_snapshot(output, undefined_repeat=False, nonzero_signal=False):
    config = VERIFY.load(ROOT/'config/ordinal_repeated_v1.json')
    lock = VERIFY.load(ROOT/'config/ordinal_repeated_v1.lock.json')
    drugs, models, repeats, folds = list(map(str, config['drug_ids'])), config['models'], ['0', '1', '2'], ['0', '1', '2', '3', '4']
    scopes = drugs+[VERIFY.MACRO]
    signal_model = 'response_association_pair_full' if nonzero_signal else None
    frames = {name: [] for name in VERIFY.CSV_FILES}
    values = {}
    for drug, model, repeat in itertools.product(drugs, models, repeats):
        d, r = drugs.index(drug), int(repeat)
        metric = dict(spearman=.15+.02*d+.005*r, rmse=.09+.003*d+.001*r,
                      mae=.05+.002*d+.001*r, r2=.3+.01*d-.002*r,
                      prediction_rms_difference_from_tissue_null=0.)
        if model == signal_model:
            metric.update(spearman=metric['spearman']+.02, rmse=metric['rmse']-.003,
                          mae=metric['mae']-.001, r2=metric['r2']+.01,
                          prediction_rms_difference_from_tissue_null=.01)
        if undefined_repeat and drug == drugs[0] and repeat == '0':
            metric['spearman'] = math.nan
        values[(drug, model, repeat)] = metric
        base = dict(gdsc_drug_id=drug, model=model, repeat=repeat, n_unique_models=VERIFY.GROUP_COUNTS[drug],
                    n_unique_biological_groups=VERIFY.GROUP_COUNTS[drug], outer_folds=5)
        frames['repeat_metrics.csv'].append(dict(**base, **metric,
            spearman_unavailable_reason='' if math.isfinite(metric['spearman']) else 'constant_response_or_prediction'))
        frames['coverage_summary.csv'].append(dict(**base, expected_oof_rows=VERIFY.GROUP_COUNTS[drug],
                                                   observed_oof_rows=VERIFY.GROUP_COUNTS[drug], status='complete'))
    for model, repeat in itertools.product(models, repeats):
        values[(VERIFY.MACRO, model, repeat)] = {metric: sum(values[(drug, model, repeat)][metric] for drug in drugs)/8
            for metric in (*VERIFY.METRICS, 'prediction_rms_difference_from_tissue_null')}
    for scope, model in itertools.product(scopes, models):
        frames['prediction_deviation_from_tissue.csv'].append(dict(gdsc_drug_id=scope, model=model,
            n_unique_biological_groups=626 if scope == VERIFY.MACRO else VERIFY.GROUP_COUNTS[scope],
            **repeat_summary([values[(scope, model, r)]['prediction_rms_difference_from_tissue_null'] for r in repeats])))
        for metric in VERIFY.METRICS:
            repeated = [values[(scope, model, r)][metric] for r in repeats]
            frames['model_intervals.csv'].append(dict(gdsc_drug_id=scope, model=model, metric=metric,
                n_unique_biological_groups=626 if scope == VERIFY.MACRO else VERIFY.GROUP_COUNTS[scope], n_bootstrap_global_groups=626,
                **interval(sum(repeated)/3)))
            frames['repeat_stability.csv'].append(dict(gdsc_drug_id=scope, comparison='model:'+model, metric=metric, **repeat_summary(repeated)))
    for scope in scopes:
        for name, coefficients, family in VERIFY.contrasts():
            for metric in VERIFY.METRICS:
                actual_family = family if metric == 'spearman' else 'descriptive'
                size = 9 if actual_family == 'primary_fixed9' else 18 if actual_family == 'representation_fixed18' else None
                sign = -1 if metric in ('rmse', 'mae') else 1
                repeated = [sign*sum(c*values[(scope, m, r)][metric] for m, c in coefficients.items()) for r in repeats]
                candidate = next((m for m, c in coefficients.items() if c == 1.), '') if len(coefficients) == 2 else 'interaction'
                reference = next((m for m, c in coefficients.items() if c == -1.), '') if len(coefficients) == 2 else 'difference_of_differences'
                frames['paired_effect_intervals.csv'].append(dict(gdsc_drug_id=scope, contrast=name, metric=metric,
                    coefficients=json.dumps(coefficients, sort_keys=True), family_name=actual_family,
                    candidate=candidate, reference=reference, positive_favors_candidate_or_positive_interaction=True,
                    effect_direction='reference_minus_candidate' if sign == -1 else 'candidate_minus_reference',
                    analysis_scope=VERIFY.ANALYSIS_SCOPE, **interval(sum(repeated)/3, size)))
                frames['repeat_stability.csv'].append(dict(gdsc_drug_id=scope, comparison=name, metric=metric, **repeat_summary(repeated)))
    for drug, model, repeat in itertools.product(drugs, models, repeats):
        null = model == 'tissue_null'
        active = model == signal_model
        base = dict(gdsc_drug_id=drug, model=model, repeat=repeat)
        selected = json.dumps({} if null else {'alpha_ratio': .1 if active else 1., 'l1_ratio': .5}, sort_keys=True)
        details = {} if null else {'alpha_ratio': dict(chosen=.1 if active else 1., minimum=.01, maximum=1., at_boundary=not active, search_extent_diagnostic=True),
                                   'l1_ratio': dict(chosen=.5, minimum=.5, maximum=1., at_boundary=True, search_extent_diagnostic=False)}
        n = 2000 if model == 'expression_strong' else 200
        frames['molecular_and_failure_summary.csv'].append(dict(**base, outer_fits=5, zero_molecular_fits=0 if null or active else 5,
            nonzero_molecular_fits=5 if active else 0, warning_outer_fits=0, fallback_outer_fits=0,
            active_feature_min='' if null else int(active), active_feature_max='' if null else int(active),
            active_feature_total=5 if active else 0, outer_unknown_tissue_samples=0))
        fits = 5 if null else 155
        frames['all_fit_events_summary.csv'].append(dict(**base, actual_fit_events=fits, all_fit_warning_events=0,
            all_fit_fallback_events=0, fit_seconds_sum=.01*fits, prediction_seconds_sum=.001*fits))
        for fold in folds:
            frames['chosen_parameters.csv'].append(dict(**base, outer_fold=fold, chosen_params=selected,
                best_inner_rmse='' if null else .1, chosen_grid_boundary=not null and not active, chosen_grid_boundary_details=json.dumps(details, sort_keys=True),
                selected_genes='' if null else n, candidate_features='' if null else n*(n-1)//2 if model.endswith('pair_full') else n,
                active_features=int(active), outer_solver_iterations=int(active), actual_alpha='' if null else .01 if active else .1, alpha_max='' if null else .1,
                warnings='[]', fallback_reason=''))
            frames['resource_tasks.csv'].append(dict(**base, outer_fold=fold, task_wall_seconds=2., task_cpu_seconds=1.,
                initialization_wall_seconds=.4, fit_work_wall_seconds=1.5, sampled_peak_rss_bytes=1000,
                windows_lifetime_peak_working_set_bytes=1200, outer_fit_seconds=.01, outer_prediction_seconds=.001))
        frames['resource_summary.csv'].append(dict(**base, outer_tasks=5, task_wall_seconds_sum=10., task_cpu_seconds_sum=5.,
            initialization_wall_seconds_sum=2., fit_work_wall_seconds_sum=7.5, sampled_peak_rss_bytes_max=1000,
            windows_lifetime_peak_working_set_bytes_max=1200, outer_fit_seconds_sum=.05, outer_prediction_seconds_sum=.005))
    for drug, model in itertools.product(drugs, models):
        null = model == 'tissue_null'
        active = model == signal_model
        frames['parameter_boundary_summary.csv'].append(dict(gdsc_drug_id=drug, model=model,
            chosen_params=json.dumps({} if null else {'alpha_ratio': .1 if active else 1., 'l1_ratio': .5}, sort_keys=True), outer_fits=15,
            selected_boundary_outer_fits=0 if null or active else 15, outer_solver_iterations_max=int(active), zero_molecular_outer_fits=0 if null or active else 15,
            scope='null_endpoint_is_predefined_not_automatically_a_grid_defect'))
        if not null:
            for (ra, fa), (rb, fb) in itertools.combinations(list(itertools.product(repeats, folds)), 2):
                frames['feature_stability.csv'].append(dict(gdsc_drug_id=drug, model=model, repeat_a=ra, outer_fold_a=fa,
                    repeat_b=rb, outer_fold_b=fb, comparison_scope='same_repeat_different_fold' if ra == rb else 'different_repeats',
                    training_sample_jaccard=.7, candidate_gene_jaccard=1., active_feature_jaccard=1. if active else '', same_sign_active_jaccard=1. if active else '',
                    empty_active_union=not active, scope='overlapping_training_sets_descriptive_not_independent_validation'))
    for name, rows in frames.items():
        write_csv(output/name, rows)
    overhead = dict(successful_pipeline_fits=26160, failed_attempt_completed_fit_events=0, resource_probe_fit_events=0,
        completed_fit_events_total=26160, failed_incomplete_fit_attempts=0, attempted_fit_total_lower_bound=26160,
        successful_task_wall_seconds_sum=1920., failed_attempt_wall_seconds_sum=0., resource_probe_wall_seconds=0.,
        controller_session_wall_seconds_sum=3590., elapsed_first_start_to_final_finish_seconds=3600., inter_session_gap_seconds=0.,
        resource_paused_session_count=0, failed_session_count=0, failed_attempts=[], probes=[], recovery_snapshots=[],
        successful_checkpoint_hashes_preserved=0, failed_checkpoint_hashes_preserved=0,
        sessions=[dict(session=1, status='completed', started_utc='2026-10-13T01:00:00+00:00',
                       finished_utc='2026-10-13T02:00:00+00:00', wall_seconds=3590., completed_tasks_this_session=960)])
    (output/'operational_overhead.json').write_text(json.dumps(overhead), encoding='utf-8')
    provenance = dict(status='completed', config_sha256=VERIFY.sha(ROOT/'config/ordinal_repeated_v1.json'),
        lock_sha256=VERIFY.sha(ROOT/'config/ordinal_repeated_v1.lock.json'), analyzer_sha256=VERIFY.sha(ROOT/'src/analyze_ordinal_repeated.py'),
        execution_commit='5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43', source_run='runs/'+config['run_name'],
        external_response_accessed=False, blocked_source_access_attempts=[], new_model_fits=0, new_permutation_runs=0,
        no_repeat_N_inflation=True, bootstrap_shared_across_models_drugs_repeats=True,
        no_joint_family_guarantee_or_confirmatory_P_values=True, noninferiority_or_equivalence_tested=False,
        confidence_scope=config['statistics']['confidence_scope'], expected_outer_tasks=lock['contract']['tasks'],
        actual_fit_events=lock['contract']['fits'], oof_rows=lock['contract']['oof_predictions'], outer_repeats=3,
        n_global_biological_groups=626, n_observed_biological_groups=626, bootstrap_draws=10000,
        bootstrap_seed=config['statistics']['bootstrap_seed'], output_sha256={name: VERIFY.sha(output/name) for name in VERIFY.FILES})
    (output/'analysis_provenance.json').write_text(json.dumps(provenance), encoding='utf-8')
    return provenance


class PublicRepeatedSnapshotTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='reo_public_repeated_')
        self.addCleanup(self.temp.cleanup)
        self.output = Path(self.temp.name)

    def rehash(self, name):
        path = self.output/'analysis_provenance.json'
        provenance = VERIFY.load(path)
        provenance['output_sha256'][name] = VERIFY.sha(self.output/name)
        path.write_text(json.dumps(provenance), encoding='utf-8')

    def mutate_csv(self, name, mutate):
        rows = VERIFY.table(self.output/name)
        mutate(rows)
        write_csv(self.output/name, rows)
        self.rehash(name)

    def test_complete_public_synthetic_contract_passes_without_private_reads(self):
        make_snapshot(self.output)
        result = VERIFY.verify(self.output)
        self.assertEqual(result['files_verified'], 14)
        self.assertEqual(result['outer_tasks'], 960)
        self.assertEqual(result['successful_pipeline_fits'], 26160)
        self.assertEqual(result['paired_effect_intervals'], 900)
        self.assertEqual(result['primary_family_effects'], 9)
        self.assertEqual(result['representation_family_effects'], 18)
        self.assertEqual(result['zero_molecular_outer_fits'], 840)
        self.assertFalse(result['private_data_opened'])
        self.assertFalse(result['new_statistics_or_fits'])

    def test_undefined_repeat_is_retained_through_macro_contrast_and_intervals(self):
        make_snapshot(self.output, undefined_repeat=True)
        self.assertEqual(VERIFY.verify(self.output)['status'], 'passed')

    def test_nonzero_effects_use_correct_RMSE_direction_and_interaction_algebra(self):
        make_snapshot(self.output, nonzero_signal=True)
        result = VERIFY.verify(self.output)
        self.assertEqual(result['nonzero_molecular_outer_fits'], 120)
        effects = VERIFY.table(self.output/'paired_effect_intervals.csv')
        primary = next(r for r in effects if r['gdsc_drug_id'] == VERIFY.MACRO and
                       r['contrast'] == 'primary:association_pair-strong_expression' and r['metric'] == 'rmse')
        self.assertAlmostEqual(float(primary['point']), .003, places=12)
        interaction = next(r for r in effects if r['gdsc_drug_id'] == VERIFY.MACRO and
                           r['contrast'] == 'interaction:pair_full-expression' and r['metric'] == 'spearman')
        self.assertAlmostEqual(float(interaction['point']), .02, places=12)
        def mutate(rows):
            row = next(r for r in rows if r['contrast'] == 'primary:association_pair-strong_expression' and r['metric'] == 'rmse')
            row['point'] = str(-float(row['point']))
        self.mutate_csv('paired_effect_intervals.csv', mutate)
        with self.assertRaisesRegex(ValueError, 'effect differs'):
            VERIFY.verify(self.output)

    def test_missing_outputs_fail_instead_of_skip(self):
        with self.assertRaisesRegex(ValueError, 'provenance is missing'):
            VERIFY.verify(self.output)
        make_snapshot(self.output)
        (self.output/'operational_overhead.json').unlink()
        with self.assertRaisesRegex(ValueError, 'Missing public aggregate'):
            VERIFY.verify(self.output)

    def test_byte_tamper_in_operational_accounting_is_rejected(self):
        make_snapshot(self.output)
        path = self.output/'operational_overhead.json'
        path.write_bytes(path.read_bytes()+b'\n')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            VERIFY.verify(self.output)

    def test_rehashed_family_size_change_is_rejected(self):
        make_snapshot(self.output)
        def mutate(rows):
            next(row for row in rows if row['family_name'] == 'primary_fixed9')['family_size'] = '18'
        self.mutate_csv('paired_effect_intervals.csv', mutate)
        with self.assertRaisesRegex(ValueError, 'family size'):
            VERIFY.verify(self.output)

    def test_rehashed_interaction_coefficients_are_checked(self):
        make_snapshot(self.output)
        def mutate(rows):
            row = next(r for r in rows if r['candidate'] == 'interaction')
            coefficients = json.loads(row['coefficients'])
            coefficients[next(iter(coefficients))] *= -1
            row['coefficients'] = json.dumps(coefficients, sort_keys=True)
        self.mutate_csv('paired_effect_intervals.csv', mutate)
        with self.assertRaisesRegex(ValueError, 'coefficients'):
            VERIFY.verify(self.output)

    def test_repeated_rows_cannot_be_counted_as_independent_extra_N(self):
        make_snapshot(self.output)
        self.mutate_csv('repeat_metrics.csv', lambda rows: rows[0].update(n_unique_biological_groups='1875'))
        with self.assertRaisesRegex(ValueError, 'inflated independent N'):
            VERIFY.verify(self.output)

    def test_rehashed_missing_effect_or_duplicate_repeat_is_rejected(self):
        make_snapshot(self.output)
        self.mutate_csv('paired_effect_intervals.csv', lambda rows: rows.pop())
        with self.assertRaisesRegex(ValueError, 'coverage'):
            VERIFY.verify(self.output)
        make_snapshot(self.output)
        self.mutate_csv('repeat_metrics.csv', lambda rows: rows.append(dict(rows[0])))
        with self.assertRaisesRegex(ValueError, 'Duplicate aggregate'):
            VERIFY.verify(self.output)

    def test_strong_gene_count_and_actual_fit_budget_are_checked(self):
        make_snapshot(self.output)
        def mutate(rows):
            row = next(r for r in rows if r['model'] == 'expression_strong')
            row.update(selected_genes='2001', candidate_features='2001')
        self.mutate_csv('chosen_parameters.csv', mutate)
        with self.assertRaisesRegex(ValueError, 'gene or representation feature budget'):
            VERIFY.verify(self.output)
        make_snapshot(self.output)
        self.mutate_csv('all_fit_events_summary.csv', lambda rows: rows[0].update(actual_fit_events='6'))
        with self.assertRaisesRegex(ValueError, 'actual fit count'):
            VERIFY.verify(self.output)

    def test_all_zero_repeat_cannot_gain_a_prediction_deviation(self):
        make_snapshot(self.output)
        def mutate(rows):
            next(r for r in rows if r['model'] == 'expression_strong')['zero_molecular_fits'] = '4'
        self.mutate_csv('molecular_and_failure_summary.csv', mutate)
        with self.assertRaisesRegex(ValueError, 'Zero/nonzero'):
            VERIFY.verify(self.output)

    def test_rehashed_operational_extra_fits_cannot_disappear_or_be_invented(self):
        make_snapshot(self.output)
        path = self.output/'operational_overhead.json'
        value = VERIFY.load(path)
        value['completed_fit_events_total'] += 1
        path.write_text(json.dumps(value), encoding='utf-8')
        self.rehash(path.name)
        with self.assertRaisesRegex(ValueError, 'accounting identity'):
            VERIFY.verify(self.output)

    def test_public_individual_identifiers_and_false_NI_claim_are_rejected(self):
        make_snapshot(self.output)
        self.mutate_csv('repeat_metrics.csv', lambda rows: [r.update(canonical_model_id='not-public') for r in rows])
        with self.assertRaisesRegex(ValueError, 'Individual-level'):
            VERIFY.verify(self.output)
        make_snapshot(self.output)
        path = self.output/'analysis_provenance.json'
        provenance = VERIFY.load(path)
        provenance['noninferiority_or_equivalence_tested'] = True
        path.write_text(json.dumps(provenance), encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'noninferiority'):
            VERIFY.verify(self.output)


if __name__ == '__main__':
    unittest.main()
