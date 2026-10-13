"""Independent, post-process integrity audit of the frozen eight-drug repeated ordinal evaluation.

Run only after the training process has stopped: reading live Windows metadata
can hold a replacement lock. This script reads prepared development artifacts
and gene identities, never expression-array bytes or original/external response
tables. Failed/partial results remain visible and never certify completion.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import struct
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import StratifiedGroupKFold

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def names_sha(values):
    return hashlib.sha256('\0'.join(values).encode('utf-8')).hexdigest()


def y_sha(values):
    return hashlib.sha256(np.asarray(values, dtype='<f8').tobytes(order='C')).hexdigest()


def aggregate_oof_stream(path):
    """Yield one CSV row at a time, releasing it after checkpoint comparison."""
    with Path(path).open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames is not None and len(reader.fieldnames) == len(set(reader.fieldnames)),
                'Aggregate OOF header is missing or has duplicate fields')
        yield from reader


def compare_oof_block(stream, expected_rows):
    """Compare every typed field exactly, including float64 round-trip bits."""
    for expected in expected_rows:
        actual = next(stream, None)
        require(actual is not None and set(actual) == set(expected), 'Aggregate OOF missing row or changed column schema')
        for key, value in expected.items():
            require(isinstance(actual[key], str), 'Aggregate CSV field missing or repeated')
            if isinstance(value, bool):
                raise ValueError('Unexpected boolean in numeric/string OOF schema')
            if isinstance(value, int):
                require(int(actual[key]) == value, 'Aggregate OOF integer context/identity differs from checkpoint: '+key)
            elif isinstance(value, float):
                parsed = float(actual[key])
                require(np.isfinite(parsed) and struct.pack('<d', parsed) == struct.pack('<d', value),
                        'Aggregate OOF float64 value differs from checkpoint: '+key)
            elif isinstance(value, str):
                require(actual[key] == value, 'Aggregate OOF string context/identity differs from checkpoint: '+key)
            else:
                raise ValueError('Unsupported checkpoint OOF field type: '+key)


def prediction_sha(rows):
    """Hash canonical IDs and exact binary float64 predictions, in ID order."""
    ordered = sorted(rows, key=lambda r: r['canonical_model_id'])
    require(len({r['canonical_model_id'] for r in ordered}) == len(ordered), 'Duplicate identity in zero-prediction fingerprint')
    digest = hashlib.sha256()
    for row in ordered:
        value = float(row['predicted_auc'])
        require(np.isfinite(value), 'Nonfinite zero-prediction fingerprint value')
        digest.update(row['canonical_model_id'].encode('utf-8'))
        digest.update(b'\0')
        digest.update(struct.pack('<d', value))
    return digest.hexdigest()


def tasks_for(cfg):
    tasks = []
    for drug, repeat, fold, model in itertools.product(
            cfg['drug_ids'], range(cfg['outer_repeats']), range(cfg['outer_folds']), cfg['models']):
        if model == 'tissue_null':
            rule, kind, budget = 'none', 'tissue_null', 0
        elif model == 'expression_strong':
            rule, kind, budget = 'mad', 'expression', cfg['feature']['strong_expression_genes']
        else:
            rule = next(r for r in cfg['selection_rules'] if model.startswith(r+'_'))
            kind, budget = model[len(rule)+1:], cfg['feature']['n_genes']
        tasks.append(dict(drug_id=int(drug), repeat=repeat, outer_fold=fold, model=model,
                          permutation_id=0, selection_rule=rule, representation=kind, n_genes=budget))
    return tasks


def filename(task):
    return f"d{task['drug_id']}_r{task['repeat']}_f{task['outer_fold']}_{task['model']}.json"


def check_payload(path):
    obj = load(path)
    digest = obj.pop('payload_sha256', None)
    actual = hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, allow_nan=False,
                                      separators=(',', ':')).encode('utf-8')).hexdigest()
    require(digest == actual, 'Checkpoint payload checksum mismatch: '+Path(path).name)
    return obj


def verify_folds(global_ids, folds, cfg):
    require(not folds.duplicated(['repeat', 'outer_fold', 'inner_fold', 'role', 'canonical_model_id']).any(), 'Duplicate fold role')
    memberships = {tuple(k): set(v.canonical_model_id) for k, v in folds.groupby(['repeat', 'outer_fold', 'inner_fold', 'role'])}
    intended = {}
    ids = global_ids.sort_values('canonical_model_id').reset_index(drop=True)
    group_by_id = dict(zip(ids.canonical_model_id, ids.biological_group_id))
    for repeat, seed in enumerate(cfg['repeat_seeds']):
        outer = StratifiedGroupKFold(cfg['outer_folds'], shuffle=True, random_state=seed)
        for f, (train, test) in enumerate(outer.split(np.zeros((len(ids), 1)), ids.tissue, ids.biological_group_id)):
            for role, positions in [('train', train), ('test', test)]:
                intended[(repeat, f, -1, role)] = set(ids.iloc[positions].canonical_model_id)
            inner = StratifiedGroupKFold(cfg['inner_folds'], shuffle=True, random_state=seed+100+f)
            for i, (it, iv) in enumerate(inner.split(np.zeros((len(train), 1)), ids.iloc[train].tissue, ids.iloc[train].biological_group_id)):
                for role, positions in [('train', it), ('validation', iv)]:
                    intended[(repeat, f, i, role)] = set(ids.iloc[train[positions]].canonical_model_id)
            require(not {group_by_id[i] for i in intended[(repeat, f, -1, 'train')]} &
                        {group_by_id[i] for i in intended[(repeat, f, -1, 'test')]}, 'Group crosses outer boundary')
    require(memberships == intended, 'Fold membership differs from frozen seed/global cohort/inner training contract')
    for _, row in folds.iterrows():
        identity = row.canonical_model_id
        source = ids[ids.canonical_model_id.eq(identity)].iloc[0]
        require(row.biological_group_id == source.biological_group_id and row.tissue == source.tissue, 'Fold identity annotation changed')
    return memberships


def null_predictions(data, training_ids, test_ids, assignment):
    """Analytic tissue means; independent of the fitted molecular estimator."""
    train = data.loc[list(training_ids)]
    target = np.asarray([assignment[name] for name in train.index], dtype=np.float64)
    overall = float(target[0]) if np.ptp(target) == 0 else float(target.mean())
    means = {t: float(target[train.tissue.eq(t)].mean()) for t in sorted(set(train.tissue))}
    return np.asarray([means.get(data.loc[name, 'tissue'], overall) for name in test_ids])


def operational_overhead(run, meta, successful_pipeline_fits, successful_task_wall_seconds_sum,
                         expected_tasks):
    """Separate completed design fits, archived failed attempts and safe pauses."""
    run = Path(run)
    failures, failed_hashes = [], set()
    failed_completed = failed_incomplete = 0
    failed_wall = 0.0
    for path in sorted((run/'failed_attempts').glob('*.json')):
        obj, digest = check_payload(path), sha(path)
        require(obj['status'] == 'failed' and obj['error'] and obj['task'] in expected_tasks and
                obj['freeze_sha256'] == meta['freeze_sha256'], 'Archived failure belongs to a different scientific task')
        require(path.stem.endswith(digest), 'Original failure SHA not retained in archive filename')
        require(all(e.get('prediction_calls') == 1 for e in obj['fit_events']), 'Incomplete audit counted as a completed fit event')
        completed = len(obj['fit_events'])
        incomplete = int('model.fit(' in obj['error'])
        wall = float(obj['resources']['wall_seconds'])
        require(np.isfinite(wall) and wall >= 0, 'Invalid failed-attempt wall time')
        failed_completed += completed
        failed_incomplete += incomplete
        failed_wall += wall
        failed_hashes.add(digest)
        failures.append(dict(file=path.name, sha256=digest, task=obj['task'],
                             completed_fit_events=completed, incomplete_fit_attempts_confirmed_by_traceback=incomplete,
                             wall_seconds=wall))
    # This protocol schedules no resource-only fits. An unexpected probe must be
    # accounted for explicitly rather than disappearing into the design budget.
    require(not list(run.glob('probe_*.json')), 'Unplanned resource probe needs separate protocol-aware accounting')
    recovered_success, recovered_failure, snapshots = set(), set(), []
    for path in sorted(run.glob('recovery_before_*.json')):
        snapshot = load(path)
        previous = snapshot['metadata']
        require(previous['freeze_sha256'] == meta['freeze_sha256'] and
                previous['execution_commit'] == meta['execution_commit'] and
                previous['prepared_input_hashes'] == meta['prepared_input_hashes'], 'Recovery changed frozen source or prepared inputs')
        for name, entry in snapshot['checkpoints'].items():
            if entry['status'] == 'completed':
                require(sha(run/'checkpoints'/name) == entry['sha256'], 'Successful checkpoint changed during recovery')
                recovered_success.add(entry['sha256'])
            else:
                require(entry['status'] == 'failed' and entry['sha256'] in failed_hashes,
                        'Original failed checkpoint is not preserved')
                recovered_failure.add(entry['sha256'])
        snapshots.append(dict(file=path.name, sha256=sha(path), checkpoint_count=len(snapshot['checkpoints'])))
    sessions = []
    for number, session in enumerate(meta.get('sessions', []), 1):
        require(session['status'] in {'completed', 'failed', 'resource_paused'} and
                'finished_utc' in session, 'Audit requires every controller session to have exited')
        require(np.isfinite(session['wall_seconds']) and session['wall_seconds'] >= 0, 'Invalid session wall time')
        sessions.append(dict(session=number, **session))
    intervals = sorted((datetime.fromisoformat(s['started_utc']), datetime.fromisoformat(s['finished_utc'])) for s in sessions)
    require(all(end >= start for start, end in intervals), 'Controller timestamps run backwards')
    require(all(intervals[i][0] >= intervals[i-1][1] for i in range(1, len(intervals))), 'Overlapping controller sessions')
    elapsed = (intervals[-1][1]-intervals[0][0]).total_seconds() if intervals else 0.0
    gap = sum((intervals[i][0]-intervals[i-1][1]).total_seconds() for i in range(1, len(intervals)))
    total = successful_pipeline_fits + failed_completed
    return dict(successful_pipeline_fits=successful_pipeline_fits,
                failed_attempt_completed_fit_events=failed_completed, resource_probe_fit_events=0,
                completed_fit_events_total=total, failed_incomplete_fit_attempts=failed_incomplete,
                attempted_fit_total_lower_bound=total+failed_incomplete,
                successful_task_wall_seconds_sum=successful_task_wall_seconds_sum,
                failed_attempt_wall_seconds_sum=failed_wall, resource_probe_wall_seconds=0.0,
                controller_session_wall_seconds_sum=sum(s['wall_seconds'] for s in sessions),
                elapsed_first_start_to_final_finish_seconds=elapsed, inter_session_gap_seconds=gap,
                resource_paused_session_count=sum(s['status'] == 'resource_paused' for s in sessions),
                failed_session_count=sum(s['status'] == 'failed' for s in sessions),
                failed_attempts=failures, probes=[], recovery_snapshots=snapshots, sessions=sessions,
                successful_checkpoint_hashes_preserved=len(recovered_success),
                failed_checkpoint_hashes_preserved=len(recovered_failure),
                interpretation=['Safe pauses start no extra scientific fit and do not omit pending tasks.',
                                'Failed-attempt completed events are extra operational work, not independent observations.',
                                'Controller wall time already includes task work; inter-session gaps include investigation/recovery, not only idle time.'])


def jsonl_stream(path):
    with Path(path).open(encoding='utf-8') as handle:
        for line in handle:
            yield json.loads(line)


def compare_json_record(stream, expected):
    actual = next(stream, None)
    require(actual is not None and json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), 'Aggregate fit/audit JSON differs from checkpoint, including typed fields')


def check_mad_budget_sets(contexts, expected_budgets=(200, 2000), require_complete=True):
    smaller, larger = sorted(expected_budgets)
    for budgets in contexts.values():
        require(set(budgets) <= {smaller, larger}, 'Unexpected MAD candidate budget')
        if smaller in budgets and larger in budgets:
            require(budgets[smaller] <= budgets[larger], 'Matched MAD genes are not contained in same-training strong MAD genes')
        elif require_complete:
            require(len(budgets) == 2, 'Missing matched/strong MAD candidate budget')


def verify_worker_isolation(resource, process_identities):
    require(resource.get('task_process_isolated') is True and
            resource.get('preimport_threads_applied') is True,
            'Task lacks fresh-process isolation or pre-import thread policy')
    pid, controller = resource.get('pid'), resource.get('controller_pid')
    created = resource.get('process_create_time')
    require(type(pid) is int and type(controller) is int and pid > 0 and controller > 0 and pid != controller,
            'Worker process identity is not distinct from its controller')
    require(type(created) in (float, int) and np.isfinite(created) and created > 0,
            'Worker process creation timestamp is missing/invalid')
    identity = (pid, float(created))
    require(identity not in process_identities, 'Worker process reused for more than one outer task')
    process_identities.add(identity)


def verify_worker_context(run, meta):
    path = Path(run)/'worker_context.json'
    require(meta['prepared_input_hashes'].get('worker_context.json') == sha(path),
            'Immutable worker context is missing or not bound to metadata')
    context = load(path)
    require(context['configuration'] == meta['configuration'] and
            context['freeze_sha256'] == meta['freeze_sha256'] and
            context['execution_commit'] == meta['execution_commit'] and
            context['source_hashes'] == meta['source_hashes'],
            'Worker context changed frozen configuration, source or lock')
    expected = {name: digest for name, digest in meta['prepared_input_hashes'].items()
                if name != 'worker_context.json'}
    require(context['prepared_input_hashes'] == expected,
            'Worker context does not bind all and only prepared development artifacts')


def verify_feature_schema(feature, event, task, cfg, all_gene_ids):
    genes = event['selected_gene_ids']
    require(genes == sorted(genes) and len(genes) == len(set(genes)), 'Selected genes are not unique stable-ID order')
    require(set(genes) <= set(all_gene_ids), 'Selected gene outside frozen protein-coding annotation')
    if feature is None:
        require(not genes and event['n_selected_genes'] is None, 'Missing fitted feature schema')
        require(not event['active_molecular_features'], 'Active features without schema')
        return
    require(feature['gene_ids'] == genes and feature['kind'] == task['representation'] and feature['selection_rule'] == task['selection_rule'], 'Outer feature identity/rule mismatch')
    indices = feature['gene_indices']
    require(len(indices) == len(set(indices)) == len(genes) and all(0 <= i < len(all_gene_ids) for i in indices) and [all_gene_ids[i] for i in indices] == genes, 'Gene indices do not map to annotation IDs')
    require(names_sha(genes) == event['selected_gene_ids_sha256'] == feature['selected_gene_ids_sha256'], 'Gene hash mismatch')
    kind = task['representation']
    if kind == 'pair_full':
        pairs = list(itertools.combinations(range(len(genes)), 2))
        names = [f'REO:{genes[i]}>{genes[j]}' for i, j in pairs]
        require(feature['pair_indices_original_input_columns'] == [[indices[i], indices[j]] for i, j in pairs], 'Full pair set omitted/reoriented')
        require(feature['constant_pairs_removed'] == feature['high_tie_pairs_removed'] == 0 and feature['retained_pairs'] == len(names), 'Full pair representation applies an unplanned pair filter')
    else:
        prefix = 'sample_rank:' if kind == 'rank' else 'log2TPM:'
        names = [prefix+g for g in genes]
    require(feature['feature_names'] == names and names_sha(names) == event['built_feature_names_sha256'] == feature['built_feature_names_sha256'], 'Encoded feature schema mismatch')
    require(event['n_features'] == len(names) and feature['training_samples'] == len(event['fit_sample_ids']), 'Feature shape/training scope mismatch')
    require(len(genes) <= task['n_genes'] and feature['gene_shortfall'] == task['n_genes']-len(genes), 'Gene budget/shortfall not preserved')
    active = feature['active_molecular_coefficients']
    require([r['feature_name'] for r in active] == event['active_molecular_features'] and
            [r['coefficient_standardized_feature_standardized_response'] for r in active] == event['active_molecular_coefficients_standardized'], 'Active coefficient payload differs from actual-fit audit')
    for row in active:
        coef = row['coefficient_standardized_feature_standardized_response']
        scale = row['training_feature_scale']
        require(row['feature_name'] in names and np.isfinite([coef, scale]).all() and coef != 0 and scale > 0, 'Invalid active molecular coefficient')
        require(np.isclose(row['coefficient_encoded_feature_original_response'], event['target_sd']*coef/scale, atol=1e-13, rtol=1e-12), 'Original-scale coefficient cannot be traced to training scale')
        require(row['direction'] == ('positive' if coef > 0 else 'negative'), 'Active coefficient sign mismatch')


def validate(run, allow_partial=False, config_path=None, lock_path=None):
    run = Path(run).resolve()
    meta = load(run/'ordinal_repeated_metadata.json')
    cfg = meta['configuration']
    require(meta.get('status') != 'running', 'Do not audit a live controller metadata file')
    require(cfg.get('permutation', {}).get('count', 0) == 0,
            'Repeated evaluation must use observed responses only')
    config_path = Path(config_path) if config_path is not None else ROOT/'config/ordinal_repeated_v1.json'
    lock_path = Path(lock_path) if lock_path is not None else ROOT/'config/ordinal_repeated_v1.lock.json'
    lock = load(lock_path)
    require(sha(lock_path) == meta['freeze_sha256'], 'Diagnostic pre-run freeze changed')
    for category in ('files', 'v1_files', 'diagnostic_files'):
        for name, expected in lock[category].items():
            require(sha(ROOT/name) == expected, 'Frozen/preserved scientific dependency changed: '+name)
    require(meta['config_sha256'] == sha(config_path) and
            cfg == load(config_path), 'Run configuration differs from frozen configuration')
    require(meta['cohort_sha256'] == load(ROOT/'config/cohort_freeze_v2.lock.json')['freeze_sha256'], 'Cohort-v2 linkage mismatch')
    require(set(meta['source_hashes']) == set(lock['execution_files']) and
            all(lock['files'][name] == expected for name, expected in meta['source_hashes'].items()),
            'Execution-source audit omitted or changed frozen dependencies')
    for name, expected in meta['source_hashes'].items():
        blob = subprocess.check_output(['git', 'show', f"{meta['execution_commit']}:{name.replace(chr(92), '/')}"] , cwd=ROOT)
        require(hashlib.sha256(blob).hexdigest() == expected, 'Recorded execution source differs from immutable Git blob: '+name)
    for name, expected in meta['prepared_input_hashes'].items():
        require(sha(run/name) == expected, 'Prepared input changed: '+name)
    verify_worker_context(run, meta)
    require(meta['external_response_accessed'] is False, 'Metadata reports external-response access')
    if meta.get('completed'):
        require(meta['data_guard']['blocked_attempts'] == 0, 'Controller attempted sealed-data access')
    manifest = pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv', dtype=str, keep_default_na=False)
    require(set(manifest.role) == {'development'}, 'Nondevelopment identity in source manifest')
    global_ids = manifest.drop_duplicates('canonical_model_id').sort_values('canonical_model_id')
    require(len(global_ids) == global_ids.biological_group_id.nunique() == 626, 'Independent group structure differs from frozen cohort')
    response = pd.read_csv(run/'development_responses.csv', dtype={'gdsc_drug_id': str}, keep_default_na=False, float_precision='round_trip')
    require(set(response.gdsc_drug_id.astype(int)) == set(cfg['drug_ids']) and not response.duplicated(['gdsc_drug_id', 'canonical_model_id']).any(), 'Prepared response task scope/duplicates changed')
    require(len(response) == lock['contract']['drug_model_observations'], 'Prepared response count differs from freeze')
    source_prefix = 'runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1/'
    if source_prefix+'development_responses.csv' in lock['files']:
        require(meta['source_prepared_sha256'] == lock['files'][source_prefix+'development_responses.csv'] and
                meta['source_formal_metadata_sha256'] == lock['files'][source_prefix+'formal_metadata.json'],
                'Inherited formal-v1 prepared source fingerprints changed')
        parent = pd.read_csv(ROOT/(source_prefix+'development_responses.csv'), dtype={'gdsc_drug_id': str},
                             keep_default_na=False, float_precision='round_trip')
        parent = parent[parent.gdsc_drug_id.astype(int).isin(cfg['drug_ids'])]
        ordering = ['gdsc_drug_id', 'canonical_model_id']
        pd.testing.assert_frame_equal(response.sort_values(ordering).reset_index(drop=True),
                                      parent.sort_values(ordering).reset_index(drop=True), check_exact=True)
    small_manifest = manifest[manifest.gdsc_drug_id.astype(int).isin(cfg['drug_ids'])]
    for drug in cfg['drug_ids']:
        current = response[response.gdsc_drug_id.astype(int).eq(drug)].set_index('canonical_model_id')
        original = small_manifest[small_manifest.gdsc_drug_id.astype(int).eq(drug)].set_index('canonical_model_id')
        require(set(current.index) == set(original.index), 'Drug observations differ from frozen IDs')
        for field in ('biological_group_id', 'tissue', 'expression_row_index'):
            require(current[field].astype(str).sort_index().equals(original[field].astype(str).sort_index()), 'Prepared biological/tissue/expression identity changed: '+field)
    memberships = verify_folds(global_ids, pd.read_csv(run/'global_fold_manifest.csv', keep_default_na=False), cfg)
    gene_table = pd.read_csv(ROOT/'data_processed/expression_genes.tsv', sep='\t', dtype=str, keep_default_na=False)
    all_gene_ids = gene_table[gene_table.gene_type.eq('protein_coding')].ensembl_id.tolist()
    expected_tasks = tasks_for(cfg)
    tasks = {filename(t): t for t in expected_tasks}
    require(len(tasks) == meta['expected_tasks'] == lock['contract']['tasks'], 'Outer task grid changed')
    grid = [dict(alpha_ratio=a, l1_ratio=l) for a, l in itertools.product(cfg['elastic_net']['alpha_ratios'], cfg['elastic_net']['l1_ratios'])]
    require(meta['expected_fit_events'] == lock['contract']['fits'], 'Fit budget metadata changed')
    require(sum(1 if task['model'] == 'tissue_null' else len(grid)*cfg['inner_folds']+1 for task in expected_tasks)
            == meta['expected_fit_events'], 'Frozen fit budget does not match actual arm/candidate Cartesian design')
    require(len(response)*cfg['outer_repeats']*len(cfg['models']) == lock['contract']['oof_predictions']
            == meta['expected_oof_predictions'], 'Frozen OOF budget does not match complete drug/arm/repeat coverage')
    registered = meta['checkpoint_hashes']
    require(set(registered) <= set(tasks), 'Registered unknown checkpoint')
    counts = Counter()
    successful_task_wall_seconds_sum = 0.0
    shared_genes = defaultdict(set)
    mad_budget_genes = {}
    active_context = None
    matched_null, zero_predictions = {}, []
    existing_names = {p.name for p in (run/'checkpoints').glob('*.json')}
    require(existing_names <= set(tasks), 'Unexpected checkpoint filename')
    if not allow_partial:
        require(existing_names == set(tasks), 'Complete diagnostic checkpoint set is missing tasks')
        require(not list((run/'checkpoints').glob('*.tmp')), 'Unfinished atomic checkpoint artifact remains')
    aggregate_stream = None if allow_partial else aggregate_oof_stream(run/'oof_predictions.csv')
    # The frozen writer concatenates checkpoints in this Cartesian task order.
    # Streaming follows the same order, so no 119,208-row Python object list or
    # duplicate full-table sorting is needed to verify exact aggregate content.
    fit_stream = None if allow_partial else jsonl_stream(run/'fit_events.jsonl')
    audit_stream = None if allow_partial else jsonl_stream(run/'outer_fit_audits.jsonl')
    process_identities = set()
    for task in expected_tasks:
        path = run/'checkpoints'/filename(task)
        if path.name not in existing_names:
            continue
        if path.name in registered:
            require(sha(path) == registered[path.name], 'Registered checkpoint changed: '+path.name)
        else:
            require(allow_partial, 'Unregistered checkpoint in completed run')
        obj = check_payload(path)
        require(obj['task'] == task and obj['freeze_sha256'] == meta['freeze_sha256'], 'Checkpoint task/freeze mismatch')
        if obj['status'] != 'completed':
            require(obj['status'] == 'failed' and obj['error'], 'Failed checkpoint lacks failure reason')
            counts['failed_checkpoints'] += 1
            require(allow_partial, 'Failed checkpoint cannot certify complete diagnostic')
            continue
        require(obj['error'] is None and not obj['data_guard']['external_response_accessed'] and obj['data_guard']['blocked_attempts'] == 0, 'Successful task contains error/sealed-data access')
        require(obj.get('execution_commit') == meta['execution_commit'] and
                obj.get('source_hashes') == meta['source_hashes'], 'Successful worker source differs from frozen controller version')
        resource = obj['resources']
        verify_worker_isolation(resource, process_identities)
        require(np.isfinite([resource['wall_seconds'], resource['cpu_seconds'], resource['sampled_peak_rss_bytes']]).all() and
                resource['wall_seconds'] > 0 and resource['cpu_seconds'] >= 0 and resource['sampled_peak_rss_bytes'] > 0,
                'Task elapsed/CPU time or measured memory record missing/invalid')
        require(all(pool['num_threads'] == cfg.get('resources', {}).get('blas_threads', 1)
                    for pool in resource['threadpools']), 'Measured native thread count exceeds frozen worker budget')
        initialization, fit_work = resource.get('initialization_wall_seconds'), resource.get('fit_work_wall_seconds')
        require(type(initialization) in (int, float) and type(fit_work) in (int, float) and
                np.isfinite([initialization, fit_work]).all() and initialization >= 0 and fit_work >= 0 and
                initialization+fit_work <= resource['wall_seconds']+1e-9,
                'Task preparation/fitting resource phases missing or inconsistent')
        successful_task_wall_seconds_sum += float(resource['wall_seconds'])
        drug, repeat, outer = task['drug_id'], task['repeat'], task['outer_fold']
        context = (drug, repeat, outer)
        if context != active_context:
            check_mad_budget_sets(mad_budget_genes,
                                 (cfg['feature']['n_genes'], cfg['feature']['strong_expression_genes']),
                                 not allow_partial)
            mad_budget_genes = {}
            active_context = context
        data = response[response.gdsc_drug_id.astype(int).eq(drug)].sort_values('canonical_model_id').set_index('canonical_model_id')
        available = set(data.index)
        assignment = dict(zip(data.index, data.response_auc))
        def expected(inner, role):
            return memberships[(repeat, outer, inner, role)] & available
        def expr(identities):
            return {int(data.loc[i, 'expression_row_index']) for i in identities}
        train, test = expected(-1, 'train'), expected(-1, 'test')
        require(obj['expected_test_ids'] == sorted(test) and set(obj['expected_training_expression_ids']) == expr(train) and set(obj['outer_test_expression_ids']) == expr(test), 'Checkpoint declared IDs differ from bound outer split')
        require([set(v) for v in obj['inner_training_expression_ids']] == [expr(expected(i, 'train')) for i in range(cfg['inner_folds'])], 'Checkpoint inner declarations differ from bound splits')
        require(obj['assigned_response_sha256'] == y_sha([assignment[i] for i in data.index]), 'Full drug assigned response fingerprint mismatch')
        rows = obj['oof_rows']
        require(len(rows) == len(test) and {r['canonical_model_id'] for r in rows} == test, 'OOF identities duplicated/missing')
        for row in rows:
            identity = row['canonical_model_id']
            require(all(row[k] == v for k, v in task.items()) and row['tissue'] == data.loc[identity, 'tissue'] and row['biological_group_id'] == data.loc[identity, 'biological_group_id'], 'OOF task/identity annotation mismatch')
            require(row['observed_auc'] == assignment[identity] and np.isfinite([row['observed_auc'], row['predicted_auc']]).all(), 'OOF outcome is not exact frozen prepared development label')
        y = np.array([r['observed_auc'] for r in rows]); prediction = np.array([r['predicted_auc'] for r in rows])
        result = obj['result']; actual_metrics = result['metrics']
        for name, value in [('rmse', np.sqrt(np.mean((y-prediction)**2))), ('mae', np.mean(np.abs(y-prediction)))]:
            require(np.isclose(actual_metrics[name], value, atol=1e-12, rtol=0), 'Outer metric differs from recorded individual predictions')
        rho = None if np.ptp(y) == 0 or np.ptp(prediction) == 0 else float(spearmanr(y, prediction).statistic)
        require(actual_metrics['spearman'] is None if rho is None else np.isclose(actual_metrics['spearman'], rho, atol=1e-12, rtol=0), 'Outer Spearman mismatch')
        if task['model'] == 'tissue_null':
            require(result['chosen_params'] == {} and result['best_inner_rmse'] is None and result['inner_scores'] == [], 'Tissue null was tuned')
        else:
            scores = result['inner_scores']
            require(len(scores) == len(grid) and [s['candidate'] for s in scores] == list(range(len(grid))) and [s['params'] for s in scores] == grid, 'Candidate search budget/order changed')
            require(all(np.isfinite(s['pooled_inner_rmse']) for s in scores), 'Invalid inner score silently ignored')
            best = min(scores, key=lambda s: s['pooled_inner_rmse'])
            require(result['chosen_params'] == best['params'] and result['best_inner_rmse'] == best['pooled_inner_rmse'], 'HPO did not use first minimum inner RMSE')
        lookup = data.reset_index().set_index('expression_row_index')
        seen, thresholds = Counter(), defaultdict(list)
        outer_event = None
        for event in obj['fit_events']:
            stage, candidate, inner = event['stage'], event['candidate'], event['inner_fold']
            seen[(stage, candidate, inner)] += 1
            if stage == 'outer_training_refit':
                require(candidate == inner == -1 and event['params'] == result['chosen_params'], 'Outer fit context/parameters mismatch')
                fitting, predicting = train, test
                outer_event = event
            else:
                require(stage == 'inner_training' and task['model'] != 'tissue_null' and 0 <= candidate < len(grid) and 0 <= inner < cfg['inner_folds'], 'Invalid inner fit event')
                require(event['params'] == grid[candidate], 'Inner actual candidate differs from frozen grid')
                fitting, predicting = expected(inner, 'train'), expected(inner, 'validation')
            actual = event['fit_sample_ids']; predicted_ids = event['prediction_sample_ids']
            require(len(actual) == len(set(actual)) and set(actual) == expr(fitting) and not set(actual) & expr(test), 'Actual fitting used held-out or omitted training identities')
            require(len(predicted_ids) == len(set(predicted_ids)) and set(predicted_ids) == expr(predicting) and not set(predicted_ids) & set(actual), 'Actual prediction identities differ from intended holdout')
            local = lookup.loc[actual]
            target = np.array([assignment[i] for i in local.canonical_model_id], dtype=np.float64)
            require(event['fit_response_sha256'] == y_sha(target), 'Actual training labels differ from frozen prepared development targets')
            target_mean = target[0] if np.ptp(target) == 0 else target.mean()
            target_sd = 1. if np.ptp(target) == 0 else target.std()
            require(np.isclose(event['target_mean'], target_mean, atol=1e-13, rtol=0) and np.isclose(event['target_sd'], target_sd, atol=1e-13, rtol=0), 'Target standardization escaped current training scope')
            require(event['model_kind'] == task['representation'] and event['selection_rule'] == task['selection_rule'], 'Actual fit model/rule differs from assigned arm')
            require(event.get('declared_model') == task['model'] and event.get('candidate_gene_budget') == task['n_genes'],
                    'Actual fit strong/matched arm or gene budget differs from task')
            require(event['prediction_calls'] == 1 and event['predict_seconds'] >= 0 and event['fit_seconds'] >= 0, 'Actual fit/predict timing missing')
            genes = event['selected_gene_ids']
            require(genes == sorted(genes) and len(genes) == len(set(genes)) <= task['n_genes'], 'Actual-fit gene order/budget invalid')
            if event['selected_gene_ids_sha256'] is not None:
                require(names_sha(genes) == event['selected_gene_ids_sha256'], 'Actual-fit selected gene hash disagrees with IDs')
            require(len(event['active_molecular_features']) == len(event['active_molecular_coefficients_standardized']) == event['n_nonzero_molecular_coeffs'], 'Zero/nonzero count differs from retained active coefficients')
            require(all(np.isfinite(v) and v != 0 for v in event['active_molecular_coefficients_standardized']), 'Active coefficient record contains zero/nonfinite value')
            if task['model'] != 'tissue_null':
                sharing = (drug, repeat, outer, task['selection_rule'], task['n_genes'], stage, inner)
                shared_genes[sharing].add(event['selected_gene_ids_sha256'])
                if task['selection_rule'] == 'mad':
                    budgets = mad_budget_genes.setdefault((stage, inner), {})
                    if task['n_genes'] in budgets:
                        require(budgets[task['n_genes']] == frozenset(genes), 'Same-training MAD genes changed across candidates/arms')
                    budgets[task['n_genes']] = frozenset(genes)
            if event['alpha_max'] is not None:
                ratio, l1 = event['params']['alpha_ratio'], event['params']['l1_ratio']
                require(np.isclose(event['actual_alpha'], event['alpha_max']*ratio, rtol=1e-12, atol=1e-14), 'Alpha not determined by current fit relative path')
                thresholds[(stage, inner)].append(event['alpha_max']*l1)
                if ratio == 1. or event['exact_tissue_constant_target']:
                    require(event['n_nonzero_molecular_coeffs'] == 0, 'Null endpoint/exact-tissue-constant fit has nonzero molecular coefficients')
            counts['warning_fits'] += bool(event['warnings'])
            counts['fallback_fits'] += bool(event['fallback_reason'])
        intended_events = {('outer_training_refit', -1, -1)}
        if task['model'] != 'tissue_null':
            intended_events |= {('inner_training', c, i) for c in range(len(grid)) for i in range(cfg['inner_folds'])}
        require(set(seen) == intended_events and set(seen.values()) == {1}, 'Actual fits missing, repeated, or budget altered')
        require(all(np.allclose(v, v[0], atol=1e-12, rtol=1e-10) for v in thresholds.values()), 'Training alpha maximum inconsistent across relative-path candidates')
        require(all(outer_event.get(k) == v for k, v in result['outer_fit_audit'].items()), 'Returned outer audit differs from recorded final actual-fit event')
        verify_feature_schema(result['feature_audit'], outer_event, task, cfg, all_gene_ids)
        null_key = (drug, repeat, outer)
        if task['model'] == 'tissue_null' or outer_event['n_nonzero_molecular_coeffs'] == 0:
            expected_null = null_predictions(data, sorted(train), [r['canonical_model_id'] for r in rows], assignment)
            require(np.allclose(prediction, expected_null, atol=2e-14, rtol=0), 'Zero-molecular predictions differ from training-only tissue means')
            digest = prediction_sha(rows)
            if task['model'] == 'tissue_null':
                matched_null[null_key] = digest
            else:
                zero_predictions.append((null_key, digest))
                counts['zero_molecular_outer_fits'] += 1
        if aggregate_stream is not None:
            compare_oof_block(aggregate_stream, rows)
            for fit_event in obj['fit_events']:
                compare_json_record(fit_stream, {**task, **fit_event})
            compare_json_record(audit_stream, dict(task=task, **obj['result'], resources=obj['resources'], data_guard=obj['data_guard']))
        counts['completed_tasks'] += 1; counts['actual_fits'] += len(obj['fit_events']); counts['oof_rows'] += len(rows)
        # Full pair schemas and completed inner-fit events are task-local.
        del obj, result, outer_event, event, rows
    require(all(len(v) == 1 for v in shared_genes.values()), 'Same training partition/rule selects different genes across representations or HPO candidates')
    check_mad_budget_sets(mad_budget_genes,
                         (cfg['feature']['n_genes'], cfg['feature']['strong_expression_genes']),
                         not allow_partial)
    for context, digest in zero_predictions:
        if context not in matched_null:
            require(allow_partial, 'Zero-molecular fit lacks matched tissue-null arm')
        else:
            require(digest == matched_null[context], 'Zero-molecular prediction not bitwise identical to its same-policy tissue null')
    archives = list((run/'failed_attempts').glob('*.json'))
    for path in archives:
        failed = check_payload(path)
        require(failed['status'] == 'failed' and failed['error'] and failed['freeze_sha256'] == meta['freeze_sha256'] and failed['task'] in expected_tasks, 'Historical failed attempt is not valid/preserved')
    counts['preserved_failed_attempts'] = len(archives)
    if not allow_partial:
        require(meta.get('completed') is True and meta['status'] == 'completed', 'Training did not complete')
        require(counts['completed_tasks'] == meta['completed_tasks'] == meta['expected_tasks'] and
                counts['actual_fits'] == meta['expected_fit_events'] and counts['oof_rows'] == meta['expected_oof_predictions'], 'Complete task/fit/prediction budget not met')
        require(set(registered) == set(tasks), 'Completed run registry lacks complete task coverage')
        for name, expected in meta['output_hashes'].items():
            require(sha(run/name) == expected, 'Aggregated output changed: '+name)
        require(next(aggregate_stream, None) is None and next(fit_stream, None) is None and next(audit_stream, None) is None,
                'Aggregate contains extra or duplicate rows after full checkpoint coverage')
        # Complete Cartesian task coverage plus exact per-task test-ID sets and
        # independently reconstructed disjoint outer folds proves exactly one
        # OOF row per frozen drug identity/arm/permutation/repeat. Every CSV row
        # was compared to those checkpoints above, including every typed field.
        require(meta.get('preserved_failed_attempts') == len(archives), 'Historical failure count omitted from metadata')
        require(meta['output_counts'] == {'tasks': counts['completed_tasks'], 'fit_events': counts['actual_fits'],
                                         'oof_predictions': counts['oof_rows'], 'failed': 0}, 'Aggregate completion summary disagrees with independently counted checkpoints')
    overhead = operational_overhead(run, meta, counts['actual_fits'], successful_task_wall_seconds_sum, expected_tasks)
    return dict(status='partial_verified' if allow_partial else 'passed', validated_utc=datetime.now(timezone.utc).isoformat(),
                execution_commit=meta['execution_commit'], repeated_lock_sha256=meta['freeze_sha256'], diagnostic_lock_sha256=meta['freeze_sha256'],
                counts=dict(counts), shared_gene_contexts_checked=len(shared_genes), fresh_worker_processes_verified=len(process_identities),
                permutation_maps_verified=0, biological_groups=len(global_ids),
                operational_overhead=overhead,
                memory_policy='One checkpoint and its OOF block at a time; aggregate CSV streamed in frozen task order; zero-molecular prediction sets retained only as exact binary-float hashes',
                original_or_external_response_or_expression_array_read=False,
                scope='Independent frozen-source, task, actual training/prediction identity and label, same-gene, active-coefficient, zero-molecular tissue-null and aggregate coverage checks',
                limitations=['Actual inner pooled RMSE is checked for candidate selection consistency; individual inner OOF predictions are not exported for independent score reconstruction.',
                             'Conditional OOF integrity does not establish superiority, noninferiority, full train-refit uncertainty or absence of every possible bias.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'config/ordinal_repeated_v1.json')
    parser.add_argument('--lock', type=Path, default=ROOT/'config/ordinal_repeated_v1.lock.json')
    parser.add_argument('--run-dir', type=Path)
    parser.add_argument('--allow-partial', action='store_true', help='Only after an interrupted controller has exited')
    args = parser.parse_args()
    run = args.run_dir.resolve() if args.run_dir else ROOT/'runs'/load(args.config)['run_name']
    try:
        report = validate(run, args.allow_partial, args.config.resolve(), args.lock.resolve())
    except Exception as exc:
        report = dict(status='failed', validated_utc=datetime.now(timezone.utc).isoformat(),
                      error_type=type(exc).__name__, error=str(exc).replace(str(ROOT), '<repository>'))
    name = 'ordinal_repeated_validation_partial.json' if args.allow_partial else 'ordinal_repeated_validation.json'
    (run/name).write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    if report['status'] == 'passed':
        (run/'operational_overhead.json').write_text(json.dumps(report['operational_overhead'], ensure_ascii=False,
                                                            indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: report.get(k) for k in ('status', 'counts', 'error')}, ensure_ascii=False))
    return int(report['status'] == 'failed')


if __name__ == '__main__':
    sys.exit(main())
