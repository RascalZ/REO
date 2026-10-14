"""Run frozen synthetic mechanism controls with task isolation and checkpoints."""
from __future__ import annotations

from ordinal_repeated_runtime import set_preimport_threads
EARLY_THREADS = set_preimport_threads()

import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import psutil
from threadpoolctl import threadpool_info

from modeling.nested_cv import make_grouped_splits
from ordinal_repeated_models import evaluate_repeated_outer_fold, model_spec
from ordinal_repeated_runtime import atomic_json, process_tree_memory, pressure_reasons, sha, system_memory
from synthetic_mechanisms import array_sha256, atomic_savez, generate_cohort

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / 'config/synthetic_mechanisms_v1_1.json'
DEFAULT_LOCK = ROOT / 'config/synthetic_mechanism_execution_v1_1.lock.json'
DEFAULT_RUN = ROOT / 'runs/fixed8_synthetic_mechanism_controls_v1_1'
RESOURCE_PAUSE_EXIT = 75


def load_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding='utf-8'))


def to_native(value):
    if isinstance(value, dict):
        return {str(k): to_native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_native(v) for v in value]
    if isinstance(value, np.ndarray):
        return to_native(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        result = float(value)
        return result if np.isfinite(result) else None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def sorted_ids_sha(values) -> str:
    canonical = '\0'.join(sorted(str(value) for value in values))
    return hashlib.sha256(canonical.encode('utf-8')).hexdigest()


def runtime_model_config(synthetic: dict, source: dict) -> dict:
    cfg = copy.deepcopy(source)
    cfg['models'] = list(synthetic['models'])
    cfg['selection_rules'] = ['mad', 'response_association']
    cfg['representations'] = ['expression', 'rank', 'pair_full']
    cfg['outer_folds'] = int(synthetic['outer_split']['n_splits'])
    cfg['outer_repeats'] = 1
    cfg['inner_folds'] = int(synthetic['inner_split']['n_splits'])
    cfg['feature']['n_genes'] = int(synthetic['elastic_net']['feature_budget'])
    cfg['feature']['diagnostic_genes'] = int(synthetic['elastic_net']['feature_budget'])
    cfg['feature']['strong_expression_genes'] = int(synthetic['elastic_net']['expression_strong_budget'])
    cfg['elastic_net']['alpha_ratios'] = list(synthetic['elastic_net']['alpha_ratios'])
    cfg['elastic_net']['l1_ratios'] = list(synthetic['elastic_net']['l1_ratios'])
    if cfg['elastic_net']['regularization_mode'] != 'alpha_max_relative':
        raise ValueError('Synthetic model configuration must use the frozen relative-alpha path')
    if cfg['elastic_net']['alpha_ratios'] != [1.0, 0.3, 0.1, 0.03, 0.01] or cfg['elastic_net']['l1_ratios'] != [0.5, 1.0]:
        raise ValueError('Synthetic EN grid differs from the frozen contract')
    return cfg


def tasks_for(config: dict) -> list[dict]:
    tasks = []
    for scenario in config['scenario_order']:
        for replicate in range(int(config['replicates_per_scenario'])):
            cohort_key = f'{scenario}_r{replicate:02d}'
            for outer_fold in range(int(config['outer_split']['n_splits'])):
                for model in config['models']:
                    tasks.append(dict(scenario=scenario,
                        scenario_index=int(config['scenario_indices'][scenario]),
                        replicate=replicate, cohort_key=cohort_key,
                        outer_fold=outer_fold, model=model))
    if len(tasks) != int(config['expected_contract']['tasks']):
        raise ValueError('Synthetic outer-task count differs from the frozen contract')
    return tasks


def checkpoint_name(task: dict) -> str:
    return (f"{task['scenario']}_r{task['replicate']:02d}_f{task['outer_fold']}_"
            f"{task['model']}.json")


def expected_fit_count(task: dict, model_config: dict) -> int:
    if task['model'] == 'tissue_null':
        return 1
    return (int(model_config['inner_folds']) * len(model_config['elastic_net']['alpha_ratios'])
            * len(model_config['elastic_net']['l1_ratios']) + 1)


def task_checkpoint_valid(path: Path, task: dict, freeze_sha: str, model_config: dict,
                          expected_samples: list[int] | None = None) -> bool:
    try:
        obj = load_json(path)
        if obj.get('status') != 'completed' or obj.get('task') != task or obj.get('freeze_sha256') != freeze_sha:
            return False
        digest = obj.get('payload_sha256')
        without = dict(obj); without.pop('payload_sha256', None)
        calculated = hashlib.sha256(json.dumps(without, sort_keys=True, ensure_ascii=False,
            allow_nan=False, separators=(',', ':')).encode('utf-8')).hexdigest()
        if digest != calculated:
            return False
        rows = obj.get('oof_rows', [])
        sample_ids = [int(row['synthetic_sample_id']) for row in rows]
        if len(sample_ids) != len(set(sample_ids)) or len(sample_ids) != 80:
            return False
        if expected_samples is not None and sorted(sample_ids) != sorted(expected_samples):
            return False
        return int(obj.get('fit_count', -1)) == expected_fit_count(task, model_config)
    except (OSError, ValueError, KeyError, TypeError):
        return False


def expected_outer_sample_ids(run_dir: Path, task: dict) -> list[int]:
    split_path = run_dir / 'prepared' / f"splits_{task['cohort_key']}.npz"
    cohort_path = run_dir / 'prepared' / f"cohort_{task['cohort_key']}.npz"
    with np.load(split_path, allow_pickle=False) as sf, np.load(cohort_path, allow_pickle=False) as cf:
        indices = np.asarray(sf[f"outer_{task['outer_fold']}_test"], dtype=int)
        return np.asarray(cf['sample_id'])[indices].astype(int).tolist()


def compact_fit_event(event: dict) -> dict:
    out = dict(event)
    fit_ids = out.pop('fit_sample_ids', [])
    prediction_ids = out.pop('prediction_sample_ids', [])
    out['fit_sample_count'] = len(fit_ids)
    out['fit_sample_ids_sha256'] = sorted_ids_sha(fit_ids)
    out['prediction_sample_count'] = len(prediction_ids)
    out['prediction_sample_ids_sha256'] = sorted_ids_sha(prediction_ids)
    selected = out.pop('selected_gene_ids', [])
    out['selected_gene_count'] = len(selected)
    out['selected_gene_ids_sha256_audit'] = sorted_ids_sha(selected) if selected else None
    if out.get('stage') == 'outer_training_refit':
        out['selected_gene_ids_outer_refit'] = selected
    active_names = out.pop('active_molecular_features', [])
    coefficients = out.pop('active_molecular_coefficients_standardized', [])
    out['active_molecular_feature_count_audit'] = len(active_names)
    out['active_molecular_feature_names_sha256'] = sorted_ids_sha(active_names) if active_names else None
    if coefficients:
        magnitudes = np.abs(np.asarray(coefficients, dtype=np.float64))
        out['active_coefficient_abs_min'] = float(np.min(magnitudes))
        out['active_coefficient_abs_median'] = float(np.median(magnitudes))
        out['active_coefficient_abs_max'] = float(np.max(magnitudes))
    return to_native(out)


def compact_outer_audit(audit: dict | None) -> dict | None:
    if audit is None:
        return None
    keep = ('model_kind', 'fitted_kind', 'fallback_reason', 'target_mean', 'target_sd',
            'fit_seconds', 'predict_seconds', 'n_iter', 'warnings', 'actual_alpha', 'alpha_max',
            'path_ratio', 'null_endpoint', 'solver_reason', 'exact_tissue_constant_target',
            'n_nonzero_molecular_coeffs', 'dual_gap', 'n_selected_genes', 'n_features',
            'n_candidate_molecular_features', 'n_molecular_features_used', 'n_tissue_features',
            'selected_gene_ids_sha256', 'built_feature_names_sha256', 'selection_rule',
            'tissue_policy', 'fit_response_sha256')
    return to_native({key: audit.get(key) for key in keep if key in audit})


def compact_feature_audit(audit: dict | None) -> dict | None:
    if audit is None:
        return None
    genes = list(audit.get('gene_ids', []))
    return to_native({
        'kind': audit.get('kind'), 'selection_rule': audit.get('selection_rule'),
        'training_samples': audit.get('training_samples'), 'input_genes': audit.get('input_genes'),
        'requested_genes': audit.get('requested_genes'), 'selected_genes': audit.get('selected_genes'),
        'gene_shortfall': audit.get('gene_shortfall'), 'eligible_genes': audit.get('eligible_genes'),
        'output_features': audit.get('output_features'), 'selected_gene_ids_sha256': audit.get('selected_gene_ids_sha256'),
        'built_feature_names_sha256': audit.get('built_feature_names_sha256'),
        'exact_tissue_constant_selection_response': audit.get('exact_tissue_constant_selection_response'),
        'gene_ids': genes,
    })


def execute_task(task: dict, run_dir: Path, context: dict, freeze_sha: str) -> bool:
    started = time.perf_counter()
    process = psutil.Process()
    cpu0 = process.cpu_times()
    samples = []
    stop = threading.Event()
    events = []
    obj = dict(schema_version='1.0', task=task, freeze_sha256=freeze_sha,
               status='failed', fit_count=0, oof_rows=[], result=None,
               fit_events=events, error=None,
               external_response_accessed=False, real_development_data_accessed=False,
               ctrp_accessed=False, execution_commit=context['execution_commit'],
               source_hashes=context['source_hashes'])

    def sample():
        while not stop.is_set():
            try:
                memory = process.memory_info()
                samples.append(dict(rss=int(memory.rss), private=int(getattr(memory, 'private', memory.vms)),
                                    **system_memory()))
            except Exception:
                pass
            stop.wait(float(context['synthetic_config']['resources']['rss_sampling_interval_seconds']))

    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    try:
        cohort_path = run_dir / 'prepared' / f"cohort_{task['cohort_key']}.npz"
        split_path = run_dir / 'prepared' / f"splits_{task['cohort_key']}.npz"
        for path in (cohort_path, split_path):
            relative = path.relative_to(run_dir).as_posix()
            if sha(path) != context['prepared_input_hashes'][relative]:
                raise ValueError('Prepared synthetic input hash mismatch: ' + relative)
        with np.load(cohort_path, allow_pickle=False) as cohort, np.load(split_path, allow_pickle=False) as split_file:
            tpm = np.asarray(cohort['tpm'], dtype=np.float64)
            y = np.asarray(cohort['response'], dtype=np.float64)
            tissue = np.asarray(cohort['tissue'], dtype=np.int16)
            sample_id = np.asarray(cohort['sample_id'], dtype=np.int32)
            groups = np.asarray(cohort['groups']).astype(str)
            fold = int(task['outer_fold'])
            train = np.asarray(split_file[f'outer_{fold}_train'], dtype=np.int64)
            test = np.asarray(split_file[f'outer_{fold}_test'], dtype=np.int64)
            inner = [(np.asarray(split_file[f'inner_{fold}_{inner_fold}_train'], dtype=np.int64),
                      np.asarray(split_file[f'inner_{fold}_{inner_fold}_test'], dtype=np.int64))
                     for inner_fold in range(int(context['synthetic_config']['inner_split']['n_splits']))]
        if tpm.shape != (240, 220) or y.shape != (240,) or tissue.shape != (240,) or sample_id.shape != (240,):
            raise ValueError('Prepared cohort shape differs from frozen DGP contract')
        if len(np.unique(groups)) != len(groups) or len(set(groups[train]) & set(groups[test])):
            raise ValueError('Synthetic outer groups overlap or are not unique')
        if set(train) & set(test) or set(train) | set(test) != set(range(len(y))):
            raise ValueError('Outer indices do not partition the prepared cohort')
        if len(test) != 80 or len(inner) != 3:
            raise ValueError('Synthetic fold size/count differs from the frozen split contract')
        x = np.column_stack((tpm, tissue.astype(np.float64), sample_id.astype(np.float64)))
        stable_gene_ids = tuple(f'SIMGENE_{index:04d}' for index in range(1, 221))
        fit_start = time.perf_counter()
        result = evaluate_repeated_outer_fold(
            x, y, groups, train, test, inner, task['model'], context['model_config'],
            stable_gene_ids, audit_sink=events.append)
        fit_work_seconds = time.perf_counter() - fit_start
        prediction = np.asarray(result.pop('predictions'), dtype=np.float64)
        if prediction.shape != (len(test),) or not np.isfinite(prediction).all():
            raise ValueError('Outer synthetic predictions are incomplete or nonfinite')
        rows = []
        for index, predicted in zip(test, prediction):
            rows.append(dict(synthetic_sample_id=int(sample_id[index]),
                biological_group_id=str(groups[index]), tissue=int(tissue[index]),
                observed_response=float(y[index]), predicted_response=float(predicted)))
        outer_audit = result.pop('outer_fit_audit', None)
        feature_audit = result.pop('feature_audit', None)
        result.pop('fit_events', None)
        result['outer_fit_audit_compact'] = compact_outer_audit(outer_audit)
        result['outer_feature_audit_compact'] = compact_feature_audit(feature_audit)
        result['fit_count'] = len(events)
        expected = expected_fit_count(task, context['model_config'])
        if len(events) != expected:
            raise ValueError(f'Fit event count {len(events)} differs from expected {expected}')
        obj['oof_rows'] = rows
        obj['fit_events'] = [compact_fit_event(event) for event in events]
        obj['fit_count'] = len(events)
        obj['result'] = to_native(result)
        obj['fit_work_wall_seconds'] = fit_work_seconds
        obj['status'] = 'completed'
    except BaseException:
        obj['error'] = traceback.format_exc()
        # Preserve any completed fit audit when a task fails before serialization.
        try:
            obj['fit_events'] = [compact_fit_event(event) for event in events]
            obj['fit_count'] = len(events)
        except Exception:
            obj['fit_events'] = []
    finally:
        stop.set(); sampler.join()
        cpu = process.cpu_times(); memory = process.memory_info()
        if not samples:
            samples.append(dict(rss=int(memory.rss), private=int(getattr(memory, 'private', memory.vms)),
                                **system_memory()))
        obj['resources'] = dict(wall_seconds=time.perf_counter() - started,
            cpu_seconds=(cpu.user + cpu.system - cpu0.user - cpu0.system),
            sampled_peak_rss_bytes=max(item['rss'] for item in samples),
            sampled_peak_private_bytes=max(item['private'] for item in samples),
            minimum_available_ram_bytes=min(item['available_ram_bytes'] for item in samples),
            minimum_commit_free_bytes=min(item['system_commit_free_bytes'] for item in samples),
            pid=process.pid, task_process_isolated=True, preimport_threads_applied=EARLY_THREADS,
            threadpools=threadpool_info(),
            memory_scope='Fresh worker process; RSS/private and system memory sampled every 0.1 seconds, not a guaranteed peak')
        path = run_dir / 'checkpoints' / checkpoint_name(task)
        atomic_json(path, to_native(obj))
    return obj['status'] == 'completed'


def verify_sources(lock_path: Path) -> tuple[dict, str]:
    lock = load_json(lock_path)
    if lock.get('status') != 'locked_before_synthetic_fits':
        raise ValueError('Synthetic execution lock has an unexpected status')
    for name, expected in lock['inputs'].items():
        if sha(ROOT / name) != expected:
            raise ValueError('Frozen synthetic method/source changed: ' + name)
    return lock, sha(lock_path)


def git_revision() -> str:
    result = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, check=True,
                            capture_output=True, text=True)
    return result.stdout.strip()


def cohort_and_split_files(run_dir: Path, config: dict) -> tuple[dict, dict]:
    prepared = run_dir / 'prepared'
    prepared.mkdir(parents=True, exist_ok=True)
    records = {}
    checks = []
    for scenario in config['scenario_order']:
        for replicate in range(int(config['replicates_per_scenario'])):
            cohort_key = f'{scenario}_r{replicate:02d}'
            data = generate_cohort(config, scenario, replicate)
            groups = np.asarray([f'{cohort_key}_sample_{index:03d}'
                                 for index in range(1, int(config['n_samples']) + 1)], dtype='U64')
            cohort_file = prepared / f'cohort_{cohort_key}.npz'
            atomic_savez(cohort_file, tpm=data['tpm'], response=data['response'],
                         tissue=data['tissue'], sample_id=data['sample_id'], groups=groups)
            seed = 2026101400 + 100 * int(config['scenario_indices'][scenario]) + replicate
            outer_seed = seed + 1000
            outer = make_grouped_splits(groups, data['tissue'], int(config['outer_split']['n_splits']), outer_seed)
            fold_arrays = {}
            coverage = np.zeros(len(groups), dtype=np.int8)
            for fold, (train, test) in enumerate(outer):
                if set(groups[train]) & set(groups[test]):
                    raise ValueError('Generated synthetic outer split leaks groups')
                coverage[test] += 1
                fold_arrays[f'outer_{fold}_train'] = np.asarray(train, dtype=np.int32)
                fold_arrays[f'outer_{fold}_test'] = np.asarray(test, dtype=np.int32)
                inner_seed = seed + 1100 + fold
                inner = make_grouped_splits(groups[train], data['tissue'][train],
                    int(config['inner_split']['n_splits']), inner_seed)
                inner_coverage = np.zeros(len(train), dtype=np.int8)
                for inner_fold, (inner_train, inner_test) in enumerate(inner):
                    if set(groups[train][inner_train]) & set(groups[train][inner_test]):
                        raise ValueError('Generated synthetic inner split leaks groups')
                    inner_coverage[inner_test] += 1
                    fold_arrays[f'inner_{fold}_{inner_fold}_train'] = np.asarray(inner_train, dtype=np.int32)
                    fold_arrays[f'inner_{fold}_{inner_fold}_test'] = np.asarray(inner_test, dtype=np.int32)
                if not np.all(inner_coverage == 1):
                    raise ValueError('Generated inner folds do not cover outer-training rows exactly once')
            if not np.all(coverage == 1):
                raise ValueError('Generated outer folds do not cover samples exactly once')
            if any(len(fold_arrays[f'outer_{fold}_test']) != 80 for fold in range(3)):
                raise ValueError('Outer synthetic folds must each contain exactly 80 samples')
            split_file = prepared / f'splits_{cohort_key}.npz'
            atomic_savez(split_file, **fold_arrays)
            records[cohort_key] = {
                'scenario': scenario, 'scenario_index': int(config['scenario_indices'][scenario]),
                'replicate': replicate, 'seed': seed, 'outer_seed': outer_seed,
                'inner_seeds': [seed + 1100 + fold for fold in range(3)],
                'cohort_file': cohort_file.relative_to(run_dir).as_posix(),
                'cohort_file_sha256': sha(cohort_file),
                'split_file': split_file.relative_to(run_dir).as_posix(),
                'split_file_sha256': sha(split_file),
                'array_hashes': {name: array_sha256(data[name]) for name in ('tpm', 'response', 'tissue', 'sample_id')},
                'checks': data['checks'],
            }
            checks.append(data['checks'])
            print(json.dumps(dict(event='synthetic_cohort_prepared', cohort=cohort_key,
                tpm_min=data['checks']['tpm_min'], tpm_max=data['checks']['tpm_max'],
                common_amplitude_rank_invariant=data['checks']['common_amplitude_rank_invariant'],
                common_amplitude_pair_invariant=data['checks']['common_amplitude_pair_invariant'])), flush=True)
    return records, {'cohorts': len(records), 'checks': checks}


def worker_context(run_dir: Path, synthetic_config: dict, model_config: dict,
                   source_hashes: dict, cohort_records: dict, freeze_sha: str,
                   execution_commit: str) -> dict:
    prepared_hashes = {}
    for record in cohort_records.values():
        prepared_hashes[record['cohort_file']] = record['cohort_file_sha256']
        prepared_hashes[record['split_file']] = record['split_file_sha256']
    context = dict(schema_version='1.0', freeze_sha256=freeze_sha,
        execution_commit=execution_commit, source_hashes=source_hashes,
        synthetic_config=synthetic_config, model_config=model_config,
        prepared_input_hashes=prepared_hashes)
    path = run_dir / 'worker_context.json'
    atomic_json(path, context)
    return context


def verify_context(run_dir: Path, expected_sha: str, lock_path: Path) -> dict:
    path = run_dir / 'worker_context.json'
    if sha(path) != expected_sha:
        raise ValueError('Worker context hash mismatch')
    context = load_json(path)
    lock = load_json(lock_path)
    if context['freeze_sha256'] != sha(lock_path):
        raise ValueError('Worker context is bound to another execution lock')
    for name, expected in context['source_hashes'].items():
        if sha(ROOT / name) != expected:
            raise ValueError('Worker frozen source changed: ' + name)
    for name, expected in context['prepared_input_hashes'].items():
        if sha(run_dir / name) != expected:
            raise ValueError('Worker prepared synthetic input changed: ' + name)
    return context


def atomic_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='raise')
        writer.writeheader(); writer.writerows(rows)
    os.replace(temp, path)


def recovery_snapshot(run_dir: Path, metadata: dict) -> None:
    checkpoints = {}
    for path in sorted((run_dir / 'checkpoints').glob('*.json')):
        checkpoints[path.name] = {'sha256': sha(path), 'status': load_json(path).get('status')}
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    atomic_json(run_dir / f'recovery_before_{stamp}.json',
                {'metadata_status': metadata.get('status'), 'checkpoint_snapshot': checkpoints})


def child_command(run_dir: Path, task_index: int, context_sha: str, lock_path: Path) -> list[str]:
    return [sys.executable, '-u', str(Path(__file__).resolve()), '--run-dir', str(run_dir),
            '--worker-task', str(task_index), '--context-sha', context_sha,
            '--lock', str(lock_path)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG)
    parser.add_argument('--lock', type=Path, default=DEFAULT_LOCK)
    parser.add_argument('--run-dir', type=Path, default=DEFAULT_RUN)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--worker-task', type=int)
    parser.add_argument('--context-sha')
    args = parser.parse_args()
    run_dir = args.run_dir.resolve(); run_dir.relative_to(ROOT / 'runs')
    config_path = args.config.resolve(); lock_path = args.lock.resolve()
    if args.worker_task is not None:
        context = verify_context(run_dir, args.context_sha, lock_path)
        config = context['synthetic_config']; tasks = tasks_for(config)
        if not 0 <= args.worker_task < len(tasks):
            raise ValueError('Worker task index is outside the frozen task list')
        ok = execute_task(tasks[args.worker_task], run_dir, context, sha(lock_path))
        return 0 if ok else 1

    config = load_json(config_path)
    lock, freeze_sha = verify_sources(lock_path)
    if sha(config_path) != lock['inputs']['config/synthetic_mechanisms_v1_1.json']:
        raise ValueError('Synthetic execution configuration differs from its frozen hash')
    source_config_path = ROOT / config['source_model_config']
    source_config = load_json(source_config_path)
    model_config = runtime_model_config(config, source_config)
    tasks = tasks_for(config)
    metadata_path = run_dir / 'synthetic_mechanism_metadata.json'
    existed_before = run_dir.exists()
    if existed_before and not args.resume:
        raise ValueError('Existing run requires explicit --resume; new run cannot reuse a prior directory')
    if not existed_before and args.resume:
        raise ValueError('--resume was requested but the frozen run directory does not exist')
    run_dir.mkdir(parents=True, exist_ok=True)
    if args.resume:
        metadata = load_json(metadata_path)
        if metadata['freeze_sha256'] != freeze_sha or metadata['execution_commit'] != git_revision():
            raise ValueError('Resume execution lock or source commit mismatch')
        for name, expected in metadata['prepared_input_hashes'].items():
            if sha(run_dir / name) != expected:
                raise ValueError('Prepared synthetic inputs changed; preserving run and stopping')
        for name, expected in metadata.get('registered_checkpoint_hashes', {}).items():
            path = run_dir / 'checkpoints' / name
            if sha(path) != expected or not task_checkpoint_valid(path, metadata['tasks_by_name'][name], freeze_sha, model_config):
                raise ValueError('Registered synthetic checkpoint changed or failed validation')
        if metadata.get('status') == 'completed':
            raise ValueError('Completed synthetic run is immutable; do not resume it')
        recovery_snapshot(run_dir, metadata)
        context = verify_context(run_dir, metadata['prepared_input_hashes']['worker_context.json'], lock_path)
    else:
        if any(run_dir.iterdir()):
            raise ValueError('Run directory is nonempty without a resumeable metadata record')
        (run_dir / 'checkpoints').mkdir(exist_ok=True)
        cohort_records, cohort_summary = cohort_and_split_files(run_dir, config)
        source_hashes = lock['inputs']
        execution_commit = git_revision()
        context = worker_context(run_dir, config, model_config, source_hashes, cohort_records, freeze_sha, execution_commit)
        prepared_hashes = dict(context['prepared_input_hashes'])
        prepared_hashes['worker_context.json'] = sha(run_dir / 'worker_context.json')
        task_by_name = {checkpoint_name(task): task for task in tasks}
        metadata = dict(schema_version='1.0', run_name=config['run_name'],
            created_utc=datetime.now(timezone.utc).isoformat(), status='prepared', completed=False,
            config_sha256=sha(config_path), freeze_sha256=freeze_sha,
            execution_commit=execution_commit, source_hashes=source_hashes,
            synthetic_config=config, model_config=model_config, cohort_records=cohort_records,
            cohort_preparation_summary=cohort_summary, prepared_input_hashes=prepared_hashes,
            expected_tasks=config['expected_contract']['tasks'], expected_fit_events=config['expected_contract']['fits'],
            expected_oof_predictions=config['expected_contract']['per_sample_model_oof_predictions'],
            completed_tasks=0, registered_checkpoint_hashes={},
            tasks_by_name=task_by_name, sessions=[], external_response_accessed=False,
            real_development_data_accessed=False, ctrp_accessed=False)
        atomic_json(metadata_path, metadata)

    session = dict(started_utc=datetime.now(timezone.utc).isoformat(), status='running',
        pending_at_start=0, completed_tasks_this_session=0, max_inflight_observed=1,
        resource_wait_seconds=0.0, child_wall_seconds_sum=0.0)
    registered = metadata['registered_checkpoint_hashes']
    pending = []
    for index, task in enumerate(tasks):
        path = run_dir / 'checkpoints' / checkpoint_name(task)
        expected_samples = expected_outer_sample_ids(run_dir, task)
        if path.exists():
            if task_checkpoint_valid(path, task, freeze_sha, model_config, expected_samples):
                if path.name in registered and sha(path) != registered[path.name]:
                    raise ValueError('Registered task checkpoint SHA changed: ' + path.name)
                registered[path.name] = sha(path)
                continue
            failed_dir = run_dir / 'failed_attempts'; failed_dir.mkdir(exist_ok=True)
            destination = failed_dir / f'{path.stem}_{sha(path)}.json'
            if destination.exists():
                raise ValueError('Repeated failed checkpoint already archived; inspect before rerunning: ' + path.name)
            shutil.move(str(path), str(destination))
            registered.pop(path.name, None)
        pending.append((index, task))
    session['pending_at_start'] = len(pending)
    metadata['sessions'].append(session); metadata.update(completed=False, status='running')
    metadata['completed_tasks'] = len(tasks) - len(pending)
    atomic_json(metadata_path, metadata)
    print(json.dumps(dict(event='synthetic_run_start', pending_tasks=len(pending), total_tasks=len(tasks),
                          expected_fits=config['expected_contract']['fits'], execution_commit=metadata['execution_commit'])), flush=True)
    started = time.perf_counter(); exit_code = 0; child = None
    parent_peak_rss = parent_peak_private = 0
    minimum_available = system_memory()['available_ram_bytes']
    minimum_commit = system_memory()['system_commit_free_bytes']
    samples_path = run_dir / f"resource_samples_session{len(metadata['sessions'])}.csv"
    fields = ['elapsed_seconds', 'process_tree_rss_bytes', 'process_tree_private_bytes', 'process_tree_count',
              'available_ram_bytes', 'system_commit_free_bytes', 'system_commit_bytes',
              'system_commit_limit_bytes', 'commit_measurement', 'inflight']
    try:
        with samples_path.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
            def record(inflight: int):
                nonlocal parent_peak_rss, parent_peak_private, minimum_available, minimum_commit
                memory = process_tree_memory()
                parent_peak_rss = max(parent_peak_rss, memory['process_tree_rss_bytes'])
                parent_peak_private = max(parent_peak_private, memory['process_tree_private_bytes'])
                minimum_available = min(minimum_available, memory['available_ram_bytes'])
                minimum_commit = min(minimum_commit, memory['system_commit_free_bytes'])
                writer.writerow(dict(elapsed_seconds=time.perf_counter() - started, inflight=inflight, **memory))
                stream.flush()
                return memory
            for index, task in pending:
                wait_start = time.perf_counter()
                while True:
                    memory = record(0); reasons = pressure_reasons(memory, config['resources'])
                    if not reasons:
                        break
                    waited = time.perf_counter() - wait_start
                    print(json.dumps(dict(event='resource_wait', task=task, reasons=reasons,
                        wait_seconds=round(waited, 1), available_ram_bytes=memory['available_ram_bytes'],
                        system_commit_free_bytes=memory['system_commit_free_bytes'])), flush=True)
                    if waited >= config['resources']['pressure_max_wait_seconds']:
                        session['status'] = 'resource_paused'; metadata['status'] = 'resource_paused'; exit_code = RESOURCE_PAUSE_EXIT
                        print(json.dumps(dict(event='resource_paused', completed_tasks=metadata['completed_tasks'], pending_task=task)), flush=True)
                        break
                    time.sleep(config['resources']['pressure_poll_seconds'])
                if exit_code:
                    break
                session['resource_wait_seconds'] += time.perf_counter() - wait_start
                child_start = time.perf_counter()
                env = dict(os.environ, PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8',
                    REO_EXECUTION_COMMIT=metadata['execution_commit'])
                for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
                    env[name] = '1'
                child = subprocess.Popen(child_command(run_dir, index,
                    metadata['prepared_input_hashes']['worker_context.json'], lock_path),
                    cwd=ROOT, env=env)
                heartbeat = -60.0
                while child.poll() is None:
                    record(1)
                    elapsed = time.perf_counter() - started
                    if elapsed - heartbeat >= 60:
                        print(json.dumps(dict(event='synthetic_heartbeat', completed_tasks=metadata['completed_tasks'],
                            elapsed_seconds=round(elapsed, 1), current_task=task)), flush=True)
                        heartbeat = elapsed
                    time.sleep(1)
                child_wall = time.perf_counter() - child_start
                session['child_wall_seconds_sum'] += child_wall
                path = run_dir / 'checkpoints' / checkpoint_name(task)
                if not path.exists():
                    atomic_json(path, dict(schema_version='1.0', task=task, freeze_sha256=freeze_sha,
                        status='failed', fit_count=0, oof_rows=[], result=None, fit_events=[],
                        error=f'Worker exited {child.returncode} without writing a checkpoint',
                        resources=dict(wall_seconds=child_wall), external_response_accessed=False,
                        real_development_data_accessed=False, ctrp_accessed=False))
                obj = load_json(path)
                if child.returncode or not task_checkpoint_valid(
                        path, task, freeze_sha, model_config, expected_outer_sample_ids(run_dir, task)):
                    session['status'] = 'failed'; metadata['status'] = 'failed'
                    atomic_json(metadata_path, metadata)
                    print(json.dumps(dict(event='synthetic_task_failed', task=task,
                        exit_code=child.returncode, fit_count=obj.get('fit_count'), error=obj.get('error'))), flush=True)
                    raise RuntimeError('Synthetic task failed; checkpoint and completed fit audit preserved')
                registered[path.name] = sha(path)
                session['completed_tasks_this_session'] += 1
                metadata['completed_tasks'] += 1
                print(json.dumps(dict(event='synthetic_task_completed', task=task,
                    completed_tasks=metadata['completed_tasks'], total_tasks=len(tasks),
                    fit_count=obj['fit_count'], wall_seconds=round(obj['resources']['wall_seconds'], 2),
                    sampled_peak_rss_bytes=obj['resources']['sampled_peak_rss_bytes'])), flush=True)
                atomic_json(metadata_path, metadata)
                child = None
        if not exit_code:
            completed = sum(1 for task in tasks if task_checkpoint_valid(
                run_dir / 'checkpoints' / checkpoint_name(task), task, freeze_sha, model_config))
            fit_count = 0; oof_count = 0
            for task in tasks:
                obj = load_json(run_dir / 'checkpoints' / checkpoint_name(task))
                fit_count += int(obj['fit_count']); oof_count += len(obj['oof_rows'])
            expected = config['expected_contract']
            if (completed, fit_count, oof_count) != (expected['tasks'], expected['fits'], expected['per_sample_model_oof_predictions']):
                raise ValueError('Completed synthetic run differs from the frozen task/fit/OOF contract')
            metadata.update(completed=True, status='completed', output_counts=dict(tasks=completed, fits=fit_count, oof_predictions=oof_count),
                checkpoint_count=len(registered), registered_checkpoint_hashes=registered,
                external_response_accessed=False, real_development_data_accessed=False, ctrp_accessed=False)
            session['status'] = 'completed'
    except BaseException:
        if child is not None and child.poll() is None:
            child.wait()
        if session.get('status') == 'running':
            session['status'] = 'failed'; metadata['status'] = 'failed'
        metadata['failure'] = traceback.format_exc()
        raise
    finally:
        session.update(wall_seconds=time.perf_counter() - started,
            finished_utc=datetime.now(timezone.utc).isoformat(),
            sampled_process_tree_peak_rss_bytes=parent_peak_rss,
            sampled_process_tree_peak_private_bytes=parent_peak_private,
            minimum_available_ram_bytes=minimum_available,
            minimum_commit_free_bytes=minimum_commit)
        atomic_json(metadata_path, metadata)
    if not exit_code:
        print(json.dumps(dict(event='synthetic_run_completed', counts=metadata['output_counts'],
            wall_seconds=round(session['wall_seconds'], 1))), flush=True)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
