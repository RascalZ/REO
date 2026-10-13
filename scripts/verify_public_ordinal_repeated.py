"""Standard-library verification of the completed public repeated-run snapshot.

Checks hashes and aggregate contracts only. Never reads private run metadata,
OOF, response tables, bootstrap identities, checkpoints or expression arrays.
Missing/incomplete outputs fail; no private statistics or model fits are run.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
MACRO = 'macro_equal_fixed_eight_drugs'
METRICS = ('spearman', 'rmse', 'mae', 'r2')
RULES = ('mad', 'response_association')
REPRESENTATIONS = ('expression', 'rank', 'pair_full')
GROUP_COUNTS = {'1931': 625, '1036': 622, '1191': 620, '1930': 620,
                '1199': 620, '1091': 620, '1010': 620, '1922': 620}
CSV_FILES = {
    'repeat_metrics.csv', 'coverage_summary.csv', 'model_intervals.csv',
    'paired_effect_intervals.csv', 'repeat_stability.csv',
    'prediction_deviation_from_tissue.csv', 'molecular_and_failure_summary.csv',
    'chosen_parameters.csv', 'parameter_boundary_summary.csv', 'feature_stability.csv',
    'resource_tasks.csv', 'resource_summary.csv', 'all_fit_events_summary.csv',
}
FILES = CSV_FILES | {'operational_overhead.json'}
INTERVAL_SCOPE = 'conditional_on_observed_OOF_records_not_full_train_refit_or_guaranteed_algorithmic_CV_coverage'
REPEAT_SCOPE = 'descriptive_repeated_split_variation_not_independent_study_SE'
ANALYSIS_SCOPE = 'development_inspired_conditional_OOF_not_independent_confirmation_or_noninferiority'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def number(value):
    require(not isinstance(value, bool), 'Boolean cannot replace a numeric count/measurement')
    return math.nan if value in ('', None) else float(value)


def integer(value):
    value = number(value)
    require(math.isfinite(value) and value.is_integer(), 'Expected a finite integer count')
    return int(value)


def boolean(value):
    require(value in ('True', 'False'), 'CSV boolean must be True or False')
    return value == 'True'


def close(actual, expected, message):
    actual, expected = number(actual), number(expected)
    require((math.isnan(actual) and math.isnan(expected)) or
            math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12), message)


def table(path):
    path = Path(path)
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)), 'Missing/duplicate CSV header: '+path.name)
        require(not set(reader.fieldnames) & {'canonical_model_id', 'biological_group_id', 'observed_auc', 'predicted_auc',
                                               'target_canonical_model_id', 'donor_canonical_model_id'},
                'Individual-level columns in public aggregate: '+path.name)
        rows = list(reader)
    require(rows and all(None not in row and None not in row.values() for row in rows), 'Empty/malformed public aggregate: '+path.name)
    return rows


def index_rows(rows, columns, expected, label):
    indexed = {tuple(row[column] for column in columns): row for row in rows}
    require(len(indexed) == len(rows), 'Duplicate aggregate key: '+label)
    require(set(indexed) == set(expected), 'Missing/unplanned aggregate coverage: '+label)
    return indexed


def contrasts():
    output = [('primary:association_pair-strong_expression',
               {'response_association_pair_full': 1., 'expression_strong': -1.}, 'primary_fixed9')]
    for rule in RULES:
        for candidate, reference in (('pair_full', 'expression'), ('rank', 'expression'), ('pair_full', 'rank')):
            output.append((rule+':'+candidate+'-'+reference,
                           {rule+'_'+candidate: 1., rule+'_'+reference: -1.},
                           'representation_fixed18' if (candidate, reference) == ('pair_full', 'expression') else 'descriptive'))
    for representation in REPRESENTATIONS:
        output.append((representation+':association-MAD',
                       {'response_association_'+representation: 1., 'mad_'+representation: -1.}, 'descriptive'))
    molecular = [rule+'_'+representation for rule in RULES for representation in REPRESENTATIONS]
    for model in molecular:
        output.append((model+'-tissue_null', {model: 1., 'tissue_null': -1.}, 'descriptive'))
    for model in molecular:
        if model != 'response_association_pair_full':
            output.append((model+'-expression_strong', {model: 1., 'expression_strong': -1.}, 'descriptive'))
    output.append(('expression_strong-tissue_null', {'expression_strong': 1., 'tissue_null': -1.}, 'descriptive'))
    for candidate, reference in (('rank', 'expression'), ('pair_full', 'rank'), ('pair_full', 'expression')):
        output.append(('interaction:'+candidate+'-'+reference,
                       {'response_association_'+candidate: 1., 'response_association_'+reference: -1.,
                        'mad_'+candidate: -1., 'mad_'+reference: 1.}, 'descriptive'))
    require(len(output) == 25 and len({c[0] for c in output}) == 25, 'Internal fixed contrast contract changed')
    return output


def check_interval(row, config, family_size=None):
    draws = config['statistics']['bootstrap_draws']
    valid, invalid = integer(row['bootstrap_valid']), integer(row['bootstrap_invalid'])
    require(integer(row['bootstrap_requested']) == draws and 0 <= valid <= draws and valid+invalid == draws,
            'Bootstrap valid/invalid counts inconsistent')
    close(row['bootstrap_valid_fraction'], valid/draws, 'Bootstrap valid fraction inconsistent')
    point, lo, hi = [number(row[k]) for k in ('point', 'conditional_95_lower', 'conditional_95_upper')]
    status, threshold = row['ci_status'], config['statistics']['minimum_valid_fraction']
    if status == 'reported':
        require(all(map(math.isfinite, (point, lo, hi))) and lo <= hi and valid/draws >= threshold,
                'Invalid reported conditional interval')
    elif status == 'withheld_undefined_point':
        require(math.isnan(point) and math.isnan(lo) and math.isnan(hi), 'Undefined point was not retained as NA')
    else:
        require(status == 'withheld_insufficient_valid_draws' and math.isfinite(point) and
                valid/draws < threshold and math.isnan(lo) and math.isnan(hi), 'Invalid withheld-interval policy')
    adjusted = [number(row[k]) for k in ('family_adjusted_conditional_lower', 'family_adjusted_conditional_upper')]
    if family_size is None:
        require(row['family_size'] == '' and all(map(math.isnan, adjusted)), 'Descriptive/model interval incorrectly family adjusted')
    else:
        require(integer(row['family_size']) == family_size, 'Frozen comparison family size changed')
        if status == 'reported':
            require(all(map(math.isfinite, adjusted)) and adjusted[0] <= lo+1e-12 and adjusted[1] >= hi-1e-12,
                    'Adjusted interval narrower than unadjusted interval')
        else:
            require(all(map(math.isnan, adjusted)), 'Withheld interval has adjusted endpoints')
    require(row['interval_scope'] == INTERVAL_SCOPE, 'Conditional uncertainty scope changed')


def check_repeat_summary(row, expected):
    actual = json.loads(row['repeat_values'])
    require(len(actual) == len(expected) == 3, 'Repeated split values do not retain all three repeats')
    for a, e in zip(actual, expected):
        close(a, e, 'Repeat value differs from complete OOF metric/contrast')
    finite = all(math.isfinite(v) for v in expected)
    wanted = {'repeat_mean': sum(expected)/3,
              'repeat_min': min(expected) if finite else math.nan,
              'repeat_max': max(expected) if finite else math.nan,
              'repeat_range': max(expected)-min(expected) if finite else math.nan,
              'repeat_sample_sd': statistics.stdev(expected) if finite else math.nan}
    for name, value in wanted.items():
        close(row[name], value, 'Three-repeat summary changed or undefined repeat dropped: '+name)
    for name, count in (('positive_repeats', sum(v > 0 for v in expected)),
                        ('negative_repeats', sum(v < 0 for v in expected)),
                        ('zero_repeats', sum(v == 0 for v in expected)),
                        ('undefined_repeats', sum(not math.isfinite(v) for v in expected))):
        require(integer(row[name]) == count, 'Repeat sign/undefined count changed')
    require(row['scope'] == REPEAT_SCOPE, 'Repeated splits incorrectly described as independent samples')


def check_operational(overhead, task_resources, contract):
    require(integer(overhead['successful_pipeline_fits']) == contract['fits'], 'Operational design fit count changed')
    failures, probes, sessions = overhead['failed_attempts'], overhead['probes'], overhead['sessions']
    require(isinstance(failures, list) and probes == [] and isinstance(sessions, list) and sessions,
            'Operational audit missing completed sessions or includes an unplanned probe')
    completed = sum(integer(f['completed_fit_events']) for f in failures)
    incomplete = sum(integer(f['incomplete_fit_attempts_confirmed_by_traceback']) for f in failures)
    require(all(integer(f['completed_fit_events']) >= 0 and integer(f['incomplete_fit_attempts_confirmed_by_traceback']) >= 0
                and number(f['wall_seconds']) >= 0 for f in failures), 'Negative failed-attempt accounting')
    require(len({f['sha256'] for f in failures}) == len(failures), 'Duplicated failed-attempt original SHA')
    require(integer(overhead['failed_attempt_completed_fit_events']) == completed and
            integer(overhead['failed_incomplete_fit_attempts']) == incomplete and
            integer(overhead['resource_probe_fit_events']) == 0, 'Failure/probe fit ledger mismatch')
    require(integer(overhead['completed_fit_events_total']) == contract['fits']+completed and
            integer(overhead['attempted_fit_total_lower_bound']) == contract['fits']+completed+incomplete,
            'Completed/attempted fit accounting identity failed')
    close(overhead['failed_attempt_wall_seconds_sum'], sum(number(f['wall_seconds']) for f in failures), 'Failed wall time omitted')
    close(overhead['resource_probe_wall_seconds'], 0, 'Unplanned resource-probe time')
    close(overhead['successful_task_wall_seconds_sum'], sum(number(row['task_wall_seconds']) for row in task_resources.values()),
          'Successful task wall sum differs from public task records')
    intervals, successful = [], 0
    for i, session in enumerate(sessions, 1):
        require(integer(session['session']) == i and session['status'] in ('completed', 'failed', 'resource_paused'), 'Live/unknown operational session')
        start, end = datetime.fromisoformat(session['started_utc']), datetime.fromisoformat(session['finished_utc'])
        require(start.tzinfo is not None and end.tzinfo is not None and end >= start, 'Invalid session timestamps')
        intervals.append((start, end))
        require(math.isfinite(number(session['wall_seconds'])) and number(session['wall_seconds']) >= 0, 'Invalid controller wall time')
        count = integer(session['completed_tasks_this_session'])
        require(0 <= count <= contract['tasks'], 'Invalid completed-session task count')
        successful += count
    require(successful <= contract['tasks'] and sessions[-1]['status'] == 'completed', 'Operational run is not completed or duplicates successful tasks')
    require(intervals == sorted(intervals) and all(intervals[i][0] >= intervals[i-1][1] for i in range(1, len(intervals))),
            'Controller sessions overlap or are out of order')
    close(overhead['controller_session_wall_seconds_sum'], sum(number(s['wall_seconds']) for s in sessions), 'Session timer sum inconsistent')
    close(overhead['elapsed_first_start_to_final_finish_seconds'], (intervals[-1][1]-intervals[0][0]).total_seconds(), 'First-to-final elapsed interval inconsistent')
    close(overhead['inter_session_gap_seconds'], sum((intervals[i][0]-intervals[i-1][1]).total_seconds() for i in range(1, len(intervals))),
          'Inter-session gap accounting inconsistent')
    require(integer(overhead['resource_paused_session_count']) == sum(s['status'] == 'resource_paused' for s in sessions) and
            integer(overhead['failed_session_count']) == sum(s['status'] == 'failed' for s in sessions), 'Pause/failure sessions conflated')
    require(0 <= integer(overhead['successful_checkpoint_hashes_preserved']) <= contract['tasks'] and
            0 <= integer(overhead['failed_checkpoint_hashes_preserved']) <= len(failures), 'Impossible preserved checkpoint count')


def verify(results_directory, root=ROOT):
    directory, root = Path(results_directory).resolve(), Path(root).resolve()
    provenance_path = directory/'analysis_provenance.json'
    require(provenance_path.is_file(), 'Completed public repeated provenance is missing')
    provenance = load(provenance_path)
    require(provenance.get('status') == 'completed', 'Public repeated analysis is incomplete')
    config_path, lock_path = root/'config/ordinal_repeated_v1.json', root/'config/ordinal_repeated_v1.lock.json'
    config, lock = load(config_path), load(lock_path)
    for field, path in [('config_sha256', config_path), ('lock_sha256', lock_path),
                        ('analyzer_sha256', root/'src/analyze_ordinal_repeated.py')]:
        require(provenance.get(field) == sha(path), 'Public frozen-source binding mismatch: '+field)
    for name in lock['execution_files']:
        # execution_files contains only public source/config/protocol files;
        # lock.files also lists private inputs, which must never be opened here.
        require(name.startswith(('src/', 'scripts/', 'config/', 'docs/')) and '..' not in Path(name).parts,
                'Unexpected nonpublic scientific source path')
        require(sha(root/name) == lock['files'][name], 'Frozen public scientific source changed: '+name)
    require(provenance['execution_commit'] == '5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43', 'Unplanned scientific execution commit')
    require(provenance['source_run'] == 'runs/'+config['run_name'], 'Unplanned source run')
    for flag in ('no_repeat_N_inflation', 'bootstrap_shared_across_models_drugs_repeats', 'no_joint_family_guarantee_or_confirmatory_P_values'):
        require(provenance.get(flag) is True, 'Missing scientific reporting boundary: '+flag)
    require(provenance.get('external_response_accessed') is False and provenance.get('blocked_source_access_attempts') == [] and
            provenance.get('noninferiority_or_equivalence_tested') is False and
            integer(provenance.get('new_model_fits')) == integer(provenance.get('new_permutation_runs')) == 0,
            'External seal, no-new-fit/null or no-noninferiority boundary changed')
    require(provenance['confidence_scope'] == config['statistics']['confidence_scope'], 'Uncertainty scope differs from frozen config')
    require(set(provenance['output_sha256']) == FILES, 'Missing/unplanned public output hash manifest')
    frames = {}
    for filename in sorted(FILES):
        path = directory/filename
        require(path.is_file(), 'Missing public aggregate: '+filename)
        require(sha(path) == provenance['output_sha256'][filename], 'Public aggregate hash mismatch: '+filename)
        if filename in CSV_FILES:
            frames[filename] = table(path)
    contract = lock['contract']
    require(config['expected_contract'] == contract and contract['tasks'] == 960 and contract['fits'] == 26160 and
            contract['oof_predictions'] == 119208 and contract['drug_model_observations'] == sum(GROUP_COUNTS.values()),
            'Frozen repeated evaluation contract changed')
    drugs, models = list(map(str, config['drug_ids'])), config['models']
    require(drugs == list(GROUP_COUNTS) and config['outer_repeats'] == 3 and config['outer_folds'] == 5 and config['permutation']['count'] == 0,
            'Fixed panel/repeated-fold/observed-label design changed')
    scopes, repeats, folds = [*drugs, MACRO], ['0', '1', '2'], ['0', '1', '2', '3', '4']
    for field, expected in [('expected_outer_tasks', 960), ('actual_fit_events', 26160), ('oof_rows', 119208),
                            ('outer_repeats', 3), ('n_global_biological_groups', 626), ('n_observed_biological_groups', 626),
                            ('bootstrap_draws', config['statistics']['bootstrap_draws']), ('bootstrap_seed', config['statistics']['bootstrap_seed'])]:
        require(integer(provenance[field]) == expected, 'Provenance count/seed differs from contract: '+field)
    grouped = set(itertools.product(drugs, models, repeats))
    task_keys = set(itertools.product(drugs, models, repeats, folds))
    points = index_rows(frames['repeat_metrics.csv'], ('gdsc_drug_id', 'model', 'repeat'), grouped, 'complete-repeat metrics')
    coverage = index_rows(frames['coverage_summary.csv'], ('gdsc_drug_id', 'model', 'repeat'), grouped, 'repeat coverage')
    for key, row in points.items():
        drug, model, repeat = key
        n = GROUP_COUNTS[drug]
        cov = coverage[key]
        for record in (row, cov):
            require(integer(record['n_unique_models']) == integer(record['n_unique_biological_groups']) == n and
                    integer(record['outer_folds']) == 5, 'Repeated rows inflated independent N or changed fold count')
        require(integer(cov['expected_oof_rows']) == integer(cov['observed_oof_rows']) == n and cov['status'] == 'complete',
                'Incomplete complete-repeat OOF coverage')
        for metric in METRICS:
            value = number(row[metric])
            require(not math.isinf(value), 'Infinite metric')
            if math.isfinite(value):
                require(-1.000000001 <= value <= 1.000000001 if metric == 'spearman' else
                        value >= 0 if metric in ('rmse', 'mae') else value <= 1.000000001, 'Metric outside mathematical range')
        require(bool(row['spearman_unavailable_reason']) == math.isnan(number(row['spearman'])), 'Undefined Spearman reason omitted/invented')
        deviation = number(row['prediction_rms_difference_from_tissue_null'])
        require(math.isfinite(deviation) and deviation >= 0, 'Invalid prediction deviation')
        if model == 'tissue_null':
            close(deviation, 0, 'Tissue reference differs from itself')
    require(sum(integer(c['observed_oof_rows']) for c in coverage.values()) == contract['oof_predictions'], 'OOF row total differs from coverage')
    def repeat_values(drug, model, metric):
        if drug == MACRO:
            return [sum(number(points[(d, model, repeat)][metric]) for d in drugs)/8 for repeat in repeats]
        return [number(points[(drug, model, repeat)][metric]) for repeat in repeats]
    model_ci = index_rows(frames['model_intervals.csv'], ('gdsc_drug_id', 'model', 'metric'),
                          itertools.product(scopes, models, METRICS), 'model intervals')
    planned = contrasts()
    effects = index_rows(frames['paired_effect_intervals.csv'], ('gdsc_drug_id', 'contrast', 'metric'),
                         itertools.product(scopes, [c[0] for c in planned], METRICS), 'paired and interaction intervals')
    stability = index_rows(frames['repeat_stability.csv'], ('gdsc_drug_id', 'comparison', 'metric'),
                           itertools.product(scopes, ['model:'+m for m in models]+[c[0] for c in planned], METRICS), 'repeat variation')
    deviations = index_rows(frames['prediction_deviation_from_tissue.csv'], ('gdsc_drug_id', 'model'),
                            itertools.product(scopes, models), 'same-policy prediction deviation')
    for (drug, model, metric), row in model_ci.items():
        values = repeat_values(drug, model, metric)
        close(row['point'], sum(values)/3, 'Model point not mean of all three complete-repeat metrics')
        require(integer(row['n_unique_biological_groups']) == (626 if drug == MACRO else GROUP_COUNTS[drug]) and
                integer(row['n_bootstrap_global_groups']) == 626, 'Model interval independent N inflated')
        check_interval(row, config)
        check_repeat_summary(stability[(drug, 'model:'+model, metric)], values)
    for (drug, model), row in deviations.items():
        require(integer(row['n_unique_biological_groups']) == (626 if drug == MACRO else GROUP_COUNTS[drug]), 'Prediction-deviation independent N changed')
        check_repeat_summary(row, repeat_values(drug, model, 'prediction_rms_difference_from_tissue_null'))
    family_counts = Counter()
    for drug in scopes:
        for name, coefficients, family in planned:
            for metric in METRICS:
                row = effects[(drug, name, metric)]
                actual_family = family if metric == 'spearman' else 'descriptive'
                size = 9 if actual_family == 'primary_fixed9' else 18 if actual_family == 'representation_fixed18' else None
                family_counts[actual_family] += 1
                require(row['family_name'] == actual_family and json.loads(row['coefficients']) == coefficients,
                        'Frozen contrast coefficients or comparison family changed')
                require(boolean(row['positive_favors_candidate_or_positive_interaction']) and row['analysis_scope'] == ANALYSIS_SCOPE,
                        'Effect direction or development-inspiration boundary changed')
                expected_direction = 'reference_minus_candidate' if metric in ('rmse', 'mae') else 'candidate_minus_reference'
                require(row['effect_direction'] == expected_direction, 'Error-metric paired direction reversed')
                candidate = next((m for m, c in coefficients.items() if c == 1.), '') if len(coefficients) == 2 else 'interaction'
                reference = next((m for m, c in coefficients.items() if c == -1.), '') if len(coefficients) == 2 else 'difference_of_differences'
                require(row['candidate'] == candidate and row['reference'] == reference, 'Paired/interaction identity changed')
                sign = -1 if metric in ('rmse', 'mae') else 1
                values = [sign*sum(c*repeat_values(drug, m, metric)[r] for m, c in coefficients.items()) for r in range(3)]
                close(row['point'], sum(values)/3, 'Paired/interaction effect differs from complete-repeat metrics')
                check_interval(row, config, size)
                check_repeat_summary(stability[(drug, name, metric)], values)
    require(len(effects) == 900 and family_counts == {'primary_fixed9': 9, 'representation_fixed18': 18, 'descriptive': 873},
            'Nine/18 separate comparison families or 900-effect budget changed')
    audits = index_rows(frames['molecular_and_failure_summary.csv'], ('gdsc_drug_id', 'model', 'repeat'), grouped, 'molecular audit')
    params = index_rows(frames['chosen_parameters.csv'], ('gdsc_drug_id', 'model', 'repeat', 'outer_fold'), task_keys, 'chosen parameters')
    events = index_rows(frames['all_fit_events_summary.csv'], ('gdsc_drug_id', 'model', 'repeat'), grouped, 'fit events')
    resources = index_rows(frames['resource_summary.csv'], ('gdsc_drug_id', 'model', 'repeat'), grouped, 'resource summaries')
    resource_tasks = index_rows(frames['resource_tasks.csv'], ('gdsc_drug_id', 'model', 'repeat', 'outer_fold'), task_keys, 'task resources')
    total_fits = zeros = nonzeros = 0
    boundaries = {}
    for key in grouped:
        drug, model, repeat = key
        summary, fit, resource = audits[key], events[key], resources[key]
        chosen = [params[(*key, fold)] for fold in folds]
        null = model == 'tissue_null'
        expected_fits = 5 if null else 155
        require(integer(summary['outer_fits']) == integer(resource['outer_tasks']) == 5 and
                integer(fit['actual_fit_events']) == expected_fits, 'Outer or actual fit count differs from frozen arm budget')
        total_fits += expected_fits
        active = [integer(r['active_features']) for r in chosen]
        require(all(a >= 0 for a in active), 'Negative active-feature count')
        zero, nonzero = (0, 0) if null else (sum(a == 0 for a in active), sum(a > 0 for a in active))
        require(integer(summary['zero_molecular_fits']) == zero and integer(summary['nonzero_molecular_fits']) == nonzero and
                integer(summary['active_feature_total']) == (0 if null else sum(active)), 'Zero/nonzero molecular counts inconsistent')
        zeros += zero; nonzeros += nonzero
        if not null:
            require(integer(summary['active_feature_min']) == min(active) and integer(summary['active_feature_max']) == max(active), 'Active-feature bounds inconsistent')
        else:
            require(not any(active) and summary['active_feature_min'] == summary['active_feature_max'] == '', 'Tissue reference counted as a molecular fit')
        for row in chosen:
            selected = json.loads(row['chosen_params'])
            if null:
                require(selected == {} and row['best_inner_rmse'] == '' and row['selected_genes'] == row['candidate_features'] == '', 'Tissue reference tuned or assigned genes')
            else:
                require(set(selected) == {'alpha_ratio', 'l1_ratio'} and selected['alpha_ratio'] in config['elastic_net']['alpha_ratios'] and
                        selected['l1_ratio'] in config['elastic_net']['l1_ratios'], 'Unplanned selected EN parameters')
                require(math.isfinite(number(row['best_inner_rmse'])) and number(row['best_inner_rmse']) >= 0, 'Invalid inner score')
                if row['selected_genes'] == '' or row['candidate_features'] == '':
                    require(bool(row['fallback_reason']) and integer(row['active_features']) == 0, 'Missing counts without retained fallback')
                    genes = features = 0
                else:
                    genes, features = integer(row['selected_genes']), integer(row['candidate_features'])
                budget = 2000 if model == 'expression_strong' else 200
                require(0 <= genes <= budget and features == (genes*(genes-1)//2 if model.endswith('pair_full') else genes) and
                        0 <= integer(row['active_features']) <= features, 'Strong/matched gene or representation feature budget invalid')
                if selected['alpha_ratio'] == 1:
                    require(integer(row['active_features']) == 0, 'Analytic null endpoint has active molecular features')
            boundary = boolean(row['chosen_grid_boundary'])
            details = json.loads(row['chosen_grid_boundary_details'])
            require(boundary == any(d['at_boundary'] and d['search_extent_diagnostic'] for d in details.values()), 'Grid-boundary flag inconsistent')
            group = boundaries.setdefault((drug, model, row['chosen_params']), [])
            group.append(row)
        if zero == 5:
            close(points[key]['prediction_rms_difference_from_tissue_null'], 0, 'All-zero molecular repeat differs from same-policy tissue reference')
            for metric in METRICS:
                close(points[key][metric], points[(drug, 'tissue_null', repeat)][metric], 'Zero-molecular complete repeat metrics differ from tissue reference')
        require(integer(summary['warning_outer_fits']) == sum(bool(json.loads(r['warnings'])) for r in chosen) and
                integer(summary['fallback_outer_fits']) == sum(bool(r['fallback_reason']) for r in chosen), 'Warnings/fallbacks omitted from outer summary')
        for name, outer in [('all_fit_warning_events', 'warning_outer_fits'), ('all_fit_fallback_events', 'fallback_outer_fits')]:
            require(integer(summary[outer]) <= integer(fit[name]) <= expected_fits, 'Actual-fit warning/fallback count inconsistent')
        tasks = [resource_tasks[(*key, fold)] for fold in folds]
        for task in tasks:
            for field in ('task_wall_seconds', 'task_cpu_seconds', 'initialization_wall_seconds', 'fit_work_wall_seconds',
                          'outer_fit_seconds', 'outer_prediction_seconds', 'sampled_peak_rss_bytes'):
                value = number(task[field])
                require(math.isfinite(value) and value >= 0, 'Missing/nonfinite/negative task resource measurement')
            require(number(task['sampled_peak_rss_bytes']) > 0 and number(task['task_wall_seconds']) > 0 and
                    number(task['initialization_wall_seconds'])+number(task['fit_work_wall_seconds']) <= number(task['task_wall_seconds'])+1e-9,
                    'Task resource phases inconsistent')
        for source, target, fn in [('task_wall_seconds', 'task_wall_seconds_sum', sum), ('task_cpu_seconds', 'task_cpu_seconds_sum', sum),
            ('fit_work_wall_seconds', 'fit_work_wall_seconds_sum', sum), ('initialization_wall_seconds', 'initialization_wall_seconds_sum', sum),
            ('outer_fit_seconds', 'outer_fit_seconds_sum', sum), ('outer_prediction_seconds', 'outer_prediction_seconds_sum', sum),
            ('sampled_peak_rss_bytes', 'sampled_peak_rss_bytes_max', max),
            ('windows_lifetime_peak_working_set_bytes', 'windows_lifetime_peak_working_set_bytes_max', max)]:
            values = [number(t[source]) for t in tasks]
            require(not any(math.isinf(v) or v < 0 for v in values), 'Invalid optional resource measurement')
            finite = [v for v in values if math.isfinite(v)]
            close(resource[target], fn(finite) if finite else math.nan, 'Resource aggregate differs from all task records')
        require(number(fit['fit_seconds_sum']) >= number(resource['outer_fit_seconds_sum'])-1e-9 and
                number(fit['prediction_seconds_sum']) >= number(resource['outer_prediction_seconds_sum'])-1e-9,
                'All-fit time omits outer fit/prediction work')
    require(total_fits == 26160 and zeros+nonzeros == 840, 'Successful fits or molecular outer-fit count incomplete')
    boundary_table = index_rows(frames['parameter_boundary_summary.csv'], ('gdsc_drug_id', 'model', 'chosen_params'), boundaries, 'parameter-boundary summary')
    for key, chosen in boundaries.items():
        row = boundary_table[key]
        require(integer(row['outer_fits']) == len(chosen) and
                integer(row['selected_boundary_outer_fits']) == sum(boolean(r['chosen_grid_boundary']) for r in chosen) and
                integer(row['zero_molecular_outer_fits']) == (0 if key[1] == 'tissue_null' else sum(integer(r['active_features']) == 0 for r in chosen)),
                'Parameter-boundary counts inconsistent')
        close(row['outer_solver_iterations_max'], max(number(r['outer_solver_iterations']) for r in chosen), 'Boundary solver maximum inconsistent')
    pairs = list(itertools.combinations(list(itertools.product(repeats, folds)), 2))
    feature_keys = {(d, m, *left, *right) for d in drugs for m in models if m != 'tissue_null' for left, right in pairs}
    feature = index_rows(frames['feature_stability.csv'], ('gdsc_drug_id', 'model', 'repeat_a', 'outer_fold_a', 'repeat_b', 'outer_fold_b'), feature_keys, 'feature stability')
    for row in feature.values():
        require(row['comparison_scope'] == ('same_repeat_different_fold' if row['repeat_a'] == row['repeat_b'] else 'different_repeats') and
                row['scope'] == 'overlapping_training_sets_descriptive_not_independent_validation', 'Feature-stability dependence scope changed')
        for field in ('training_sample_jaccard', 'candidate_gene_jaccard', 'active_feature_jaccard', 'same_sign_active_jaccard'):
            value = number(row[field])
            require(math.isnan(value) or 0 <= value <= 1, 'Jaccard outside [0,1]')
        empty = boolean(row['empty_active_union'])
        require(empty == math.isnan(number(row['active_feature_jaccard'])) == math.isnan(number(row['same_sign_active_jaccard'])),
                'Empty active sets must retain undefined stability')
        if not empty:
            require(number(row['same_sign_active_jaccard']) <= number(row['active_feature_jaccard'])+1e-12, 'Signed overlap exceeds feature overlap')
    check_operational(load(directory/'operational_overhead.json'), resource_tasks, contract)
    return dict(status='passed', files_verified=len(FILES), outer_tasks=len(task_keys), successful_pipeline_fits=total_fits,
                paired_effect_intervals=len(effects), primary_family_effects=family_counts['primary_fixed9'],
                representation_family_effects=family_counts['representation_fixed18'], independent_biological_groups=626,
                zero_molecular_outer_fits=zeros, nonzero_molecular_outer_fits=nonzeros,
                private_data_opened=False, new_statistics_or_fits=False,
                scope='Public hashes, source binding and aggregate algebra/coverage only; not private OOF/bootstrap reproduction, an independent study, algorithmic uncertainty or scientific validity certification.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', type=Path, default=ROOT/'results/ordinal_repeated')
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.results_dir), ensure_ascii=False))
        return 0
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as exc:
        print(json.dumps(dict(status='failed', error=str(exc), private_data_opened=False), ensure_ascii=False))
        return 1


if __name__ == '__main__':
    sys.exit(main())
