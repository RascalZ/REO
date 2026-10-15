"""Frozen, response-free final-model training on GDSC development data only.

This program must not run before the external-validation protocol, config,
source and lock have been committed.  It never opens CTRP data.  Each of the
fixed 64 drug/arm artifacts is fitted in one fresh process, then reloaded and
checked in another fresh process.  Interrupted aggregation is rebuilt solely
from hashed completed checkpoints; no completed model is refitted.
"""
from __future__ import annotations

from ordinal_repeated_runtime import set_preimport_threads
EARLY_THREADS = set_preimport_threads()
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import threading
import time
import traceback

import joblib
import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_info, threadpool_limits

from final_model_v1 import (fit_full_development_model, local_inner_splits,
                            make_global_inner_manifest, prediction_sha256,
                            validate_global_inner_manifest)
from ordinal_repeated_models import model_spec
from ordinal_repeated_runtime import (atomic_json, pressure_reasons,
                                      process_tree_memory, sha, system_memory)
from run_development_pilot import ROOT, install_data_guard, verify_development_inputs


SOURCE_RUN = ROOT/'runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1'
DEFAULT_RUN = ROOT/'runs/gdsc2-cancer-v2-fixed8-full-development-final-v1'
DEFAULT_EXTERNAL_CONFIG = ROOT/'config/external_validation_v1.json'
DEFAULT_ORDINAL_CONFIG = ROOT/'config/ordinal_repeated_v1.json'
DEFAULT_LOCK = ROOT/'config/external_validation_v1.lock.json'
PREPARED_NAMES = ('development_responses.csv', 'global_inner_fold_manifest.csv', 'worker_context.json')


def now_utc():
    return datetime.now(timezone.utc).isoformat()


def current_environment():
    names = ('numpy', 'pandas', 'scipy', 'scikit-learn', 'joblib', 'threadpoolctl', 'psutil')
    return dict(python=platform.python_version(),
                packages={name: importlib.metadata.version(name) for name in names})


def task_list(external_config, ordinal_config):
    drugs = list(map(int, external_config['drug_ids']))
    models = list(external_config['models'])
    if drugs != list(map(int, ordinal_config['drug_ids'])) or models != list(ordinal_config['models']):
        raise ValueError('Final 8-drug/8-arm panel differs from frozen repeated evaluation')
    if len(drugs) != len(set(drugs)) or len(models) != len(set(models)) or len(drugs) != 8 or len(models) != 8:
        raise ValueError('Final panel must contain exactly eight unique drugs and eight unique arms')
    if int(external_config['final_training']['inner_folds']) != int(ordinal_config['inner_folds']):
        raise ValueError('Final inner-fold count differs from frozen model policy')
    if int(external_config['final_training']['inner_folds']) != 3:
        raise ValueError('Final full-development model requires three inner folds')
    if int(external_config['final_training']['inner_seed']) != 2026101501:
        raise ValueError('Unexpected final-development split seed')
    return [dict(drug_id=drug, model=model, **model_spec(model, ordinal_config))
            for drug in drugs for model in models]


def task_stem(task):
    return f"d{task['drug_id']}_{task['model']}"


def _payload_digest(obj):
    data = {key: value for key, value in obj.items() if key != 'payload_sha256'}
    return hashlib.sha256(json.dumps(data, sort_keys=True, ensure_ascii=False,
                                     allow_nan=False, separators=(',', ':')).encode('utf-8')).hexdigest()


def save_checkpoint(path, obj):
    obj = dict(obj)
    obj['payload_sha256'] = _payload_digest(obj)
    atomic_json(path, obj)


def fit_journal_path(out, task):
    return Path(out)/'fit_event_journals'/(task_stem(task)+'.jsonl')


def read_fit_journal(path):
    """Return the fsynced JSONL prefix, preserving any torn tail on disk."""
    path = Path(path)
    events, parse_error = [], None
    with path.open('rb') as handle:
        for number, line in enumerate(handle, 1):
            try:
                if not line.endswith(b'\n'):
                    raise ValueError('unterminated JSONL event')
                event = json.loads(line.decode('utf-8'))
                if not isinstance(event, dict):
                    raise ValueError('fit event is not an object')
                events.append(event)
            except (UnicodeError, json.JSONDecodeError, ValueError) as exc:
                parse_error = f'line {number}: {exc}'
                if handle.read(1):
                    parse_error += '; later bytes remain unread'
                break
    return events, parse_error


def _attach_fit_journal(obj, out, task):
    path = fit_journal_path(out, task)
    if not path.is_file():
        return
    events, parse_error = read_fit_journal(path)
    obj['fit_event_journal_path'] = 'fit_event_journals/'+path.name
    obj['fit_event_journal_sha256'] = sha(path)
    obj['fit_count'] = len(events)
    obj['fit_events'] = events
    obj['fit_event_journal_parse_error'] = parse_error


def preserve_interrupted_task(out, task, freeze_sha, reason):
    """Fail closed after a worker dies before writing its atomic checkpoint."""
    path = Path(out)/'checkpoints'/(task_stem(task)+'.json')
    if path.exists():
        raise ValueError('Interrupted checkpoint already exists')
    artifact_path = Path(out)/'models'/(task_stem(task)+'.joblib')
    journal_path = fit_journal_path(out, task)
    journal_path.parent.mkdir(exist_ok=True)
    if not journal_path.exists():
        with journal_path.open('x', encoding='utf-8', newline='\n') as handle:
            handle.flush()
            os.fsync(handle.fileno())
    obj = dict(schema_version='1.0', task=task, freeze_sha256=freeze_sha,
        status='failed', fit_count=0, fit_events=[], result=None, error=reason,
        artifact_path='models/'+artifact_path.name,
        artifact_sha256=sha(artifact_path) if artifact_path.is_file() else None)
    _attach_fit_journal(obj, out, task)
    save_checkpoint(path, obj)
    return obj


def recover_interrupted_task_if_needed(out, task, freeze_sha):
    path = Path(out)/'checkpoints'/(task_stem(task)+'.json')
    if path.exists():
        return
    if fit_journal_path(out, task).exists() or \
            (Path(out)/'models'/(task_stem(task)+'.joblib')).exists():
        preserve_interrupted_task(out, task, freeze_sha,
            'Interrupted worker left a fit event journal or artifact without a checkpoint')
        raise ValueError('Interrupted final-model task preserved as failed; inspect before recovery')


def read_checkpoint(path, task, freeze_sha, *, verify_artifact=True):
    obj = json.loads(Path(path).read_text(encoding='utf-8'))
    if obj.get('payload_sha256') != _payload_digest(obj):
        raise ValueError('Final-model checkpoint payload hash mismatch: '+str(path))
    if obj.get('task') != task or obj.get('freeze_sha256') != freeze_sha:
        raise ValueError('Final-model checkpoint task/freeze mismatch: '+str(path))
    if obj['status'] not in {'fitted_pending_verification', 'completed', 'failed'}:
        raise ValueError('Unknown final-model checkpoint state')
    journal_path = Path(path).parents[1]/'fit_event_journals'/(task_stem(task)+'.jsonl')
    if obj.get('fit_event_journal_path') != 'fit_event_journals/'+journal_path.name or \
            not journal_path.is_file() or sha(journal_path) != obj.get('fit_event_journal_sha256'):
        raise ValueError('Final-model fit event journal path/SHA mismatch: '+str(path))
    events, parse_error = read_fit_journal(journal_path)
    if obj.get('fit_count') != len(events) or obj.get('fit_events') != events or \
            obj.get('fit_event_journal_parse_error') != parse_error:
        raise ValueError('Final-model fit event journal/checkpoint mismatch: '+str(path))
    if obj['status'] != 'failed' and parse_error is not None:
        raise ValueError('Completed final-model fit journal has a torn event')
    if obj['status'] in {'fitted_pending_verification', 'completed'}:
        if obj.get('artifact_path') != 'models/'+task_stem(task)+'.joblib':
            raise ValueError('Final-model artifact name/task mismatch')
        if verify_artifact:
            artifact = Path(path).parents[1]/obj['artifact_path']
            if not artifact.is_file() or sha(artifact) != obj['artifact_sha256']:
                raise ValueError('Final-model artifact SHA mismatch: '+str(artifact))
        expected = 1 if task['model'] == 'tissue_null' else 31
        if obj.get('fit_count') != expected or len(obj.get('fit_events', [])) != expected:
            raise ValueError('Final-model fit-count contract mismatch')
        if obj['status'] == 'completed' and not obj.get('fresh_process_reload_verified'):
            raise ValueError('Completed model lacks fresh-process reload verification')
    return obj


def validate_prepared_checkpoint(obj, task, response, manifest, n_folds):
    """Bind a resumed fit to the exact prepared identities, labels and folds."""
    if obj['status'] == 'failed':
        return
    data = response.loc[response.gdsc_drug_id.eq(str(task['drug_id']))].sort_values(
        'canonical_model_id', kind='stable').reset_index(drop=True)
    if obj['development_canonical_ids'] != data.canonical_model_id.tolist():
        raise ValueError('Checkpoint development canonical IDs differ from prepared input')
    if obj['development_groups'] != data.biological_group_id.tolist():
        raise ValueError('Checkpoint biological groups differ from prepared input')
    indices = data.expression_row_index.astype(int).tolist()
    if obj['development_expression_row_indices'] != indices or obj['result']['development_sample_ids'] != indices:
        raise ValueError('Checkpoint expression mapping differs from prepared input')
    y = data.response_auc.to_numpy(dtype=np.float64)
    expected_y = hashlib.sha256(np.ascontiguousarray(y, dtype='<f8').tobytes()).hexdigest()
    if obj['result']['development_response_sha256'] != expected_y:
        raise ValueError('Checkpoint response differs from prepared input')
    inner = local_inner_splits(data, manifest, n_folds)
    expected_validation = [data.iloc[iv].expression_row_index.astype(int).tolist() for _, iv in inner]
    if obj['inner_validation_expression_row_indices'] != expected_validation:
        raise ValueError('Checkpoint inner folds differ from prepared input')
    if obj['result']['full_fit_audit']['fit_sample_ids'] != indices:
        raise ValueError('Final fit audit lacks full development membership')
    for event in obj['fit_events']:
        if event['stage'] == 'full_development_refit':
            expected = indices
        elif event['stage'] == 'inner_training':
            expected = data.iloc[inner[event['inner_fold']][0]].expression_row_index.astype(int).tolist()
        else:
            raise ValueError('Unexpected final-fit event stage')
        if event['fit_sample_ids'] != expected:
            raise ValueError('Final-fit event membership differs from frozen folds')


def _checked_run_dir(path):
    out = Path(path).resolve()
    out.relative_to((ROOT/'runs').resolve())
    if out == (ROOT/'runs').resolve():
        raise ValueError('Run output must be a child of runs/')
    return out


def _blob_sha(revision, name):
    blob = subprocess.check_output(['git', 'show', revision+':'+name], cwd=ROOT)
    return hashlib.sha256(blob).hexdigest()


def verify_freeze(lock_path, external_config_path, ordinal_config_path):
    lock = json.loads(Path(lock_path).read_text(encoding='utf-8'))
    if lock.get('status') != 'frozen_pre_final_fit':
        raise ValueError('External-validation lock is not in a frozen state')
    files = lock['files']
    for name, expected in files.items():
        if sha(ROOT/name) != expected:
            raise ValueError('Frozen external-validation file changed: '+name)
    for path in (external_config_path, ordinal_config_path):
        name = Path(path).resolve().relative_to(ROOT.resolve()).as_posix()
        if files.get(name) != sha(path):
            raise ValueError('Final model config not bound by external-validation lock: '+name)
    revision = subprocess.check_output(['git', 'rev-parse', lock['execution_commit']],
                                       cwd=ROOT, text=True).strip()
    if revision != lock['execution_commit']:
        raise ValueError('Execution commit must be a full immutable Git SHA')
    source_hashes = {}
    for name in lock['execution_files']:
        if name not in files or _blob_sha(revision, name) != files[name]:
            raise ValueError('Execution source differs from frozen pre-lock Git commit: '+name)
        source_hashes[name] = files[name]
    contract = lock['contract']
    if contract.get('final_model_tasks') != 64 or contract.get('final_fit_events') != 1744:
        raise ValueError('External lock has unexpected final-fit count contract')
    if lock.get('environment') != current_environment():
        raise ValueError('Final-fit Python/package environment differs from pre-fit lock')
    return lock, revision, source_hashes


def _response_source_verified():
    meta = json.loads((SOURCE_RUN/'formal_metadata.json').read_text(encoding='utf-8'))
    source = SOURCE_RUN/'development_responses.csv'
    if not meta['completed'] or sha(source) != meta['prepared_input_hashes']['development_responses.csv']:
        raise ValueError('Historical formal development response provenance failed')
    return source, sha(source), sha(SOURCE_RUN/'formal_metadata.json')


def _archive_stale_temp(path, out):
    path = Path(path)
    if path.exists():
        folder = out/'interrupted_partial_files'
        folder.mkdir(exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        os.replace(path, folder/(path.name+'.'+stamp+'.'+sha(path)))


def _atomic_write_text(path, text, out):
    path = Path(path)
    temp = path.with_name(path.name+'.tmp')
    _archive_stale_temp(temp, out)
    with temp.open('w', encoding='utf-8', newline='\n') as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    _archive_stale_temp(path, out)
    os.replace(temp, path)


def _atomic_joblib(path, obj, out):
    path = Path(path)
    temp = path.with_name(path.name+'.tmp')
    _archive_stale_temp(temp, out)
    joblib.dump(obj, temp, compress=3)
    with temp.open('r+b') as handle:
        os.fsync(handle.fileno())
    os.replace(temp, path)


def load_context(out, context_sha):
    out = _checked_run_dir(out)
    path = out/'worker_context.json'
    if not context_sha or sha(path) != context_sha:
        raise ValueError('Final-model worker context SHA mismatch')
    context = json.loads(path.read_text(encoding='utf-8'))
    if context.get('environment') != current_environment():
        raise ValueError('Final-model worker environment differs from pre-fit lock')
    for name, expected in context['prepared_input_hashes'].items():
        if sha(out/name) != expected:
            raise ValueError('Prepared final-model input changed: '+name)
    for name, expected in context['source_hashes'].items():
        if sha(ROOT/name) != expected:
            raise ValueError('Final-model scientific source changed: '+name)
    return context


def verify_worker_freeze(context):
    """A direct worker/verifier invocation must pass the immutable pre-fit gate."""
    lock, revision, source_hashes = verify_freeze(
        DEFAULT_LOCK, DEFAULT_EXTERNAL_CONFIG, DEFAULT_ORDINAL_CONFIG)
    if context.get('freeze_sha256') != sha(DEFAULT_LOCK) or \
            context.get('execution_commit') != revision or \
            context.get('source_hashes') != source_hashes or \
            context.get('environment') != lock['environment']:
        raise ValueError('Worker context differs from frozen protocol/source/environment')
    if context.get('external_config') != json.loads(DEFAULT_EXTERNAL_CONFIG.read_text(encoding='utf-8')) or \
            context.get('ordinal_config') != json.loads(DEFAULT_ORDINAL_CONFIG.read_text(encoding='utf-8')):
        raise ValueError('Worker config differs from frozen pre-fit configuration')
    for name, expected in context['development_file_hashes'].items():
        if name not in ('data_processed/expression_genes.tsv', 'data_processed/expression_tpm.npy') or \
                sha(ROOT/name) != expected:
            raise ValueError('Worker development expression differs from frozen prepared input')
    if set(context['development_file_hashes']) != {
            'data_processed/expression_genes.tsv', 'data_processed/expression_tpm.npy'}:
        raise ValueError('Worker development expression manifest incomplete')
    source, expected_response_sha, _ = _response_source_verified()
    if context['prepared_input_hashes'].get('development_responses.csv') != expected_response_sha:
        raise ValueError('Worker prepared response differs from frozen development source')
    return lock


def load_drug(out, task, context):
    """Load only frozen GDSC development response and frozen CCLE TPM."""
    response = pd.read_csv(out/'development_responses.csv',
        dtype={'gdsc_drug_id': str, 'canonical_model_id': str,
               'biological_group_id': str, 'tissue': str},
        keep_default_na=False, float_precision='round_trip')
    data = response.loc[response.gdsc_drug_id.eq(str(task['drug_id']))].sort_values(
        'canonical_model_id', kind='stable').reset_index(drop=True)
    if data.empty or data.canonical_model_id.duplicated().any() or data.biological_group_id.duplicated().any():
        raise ValueError('Missing or duplicate full-development drug observations')
    if not set(data.tissue).issubset(context['tissue_code']):
        raise ValueError('Unknown development tissue code')
    genes = pd.read_csv(ROOT/'data_processed/expression_genes.tsv', sep='\t',
                        dtype=str, keep_default_na=False)
    gene_idx = genes.index[genes.gene_type.eq('protein_coding')].to_numpy()
    stable = tuple(genes.iloc[gene_idx].ensembl_id)
    if len(stable) != len(set(stable)) or len(stable) == 0:
        raise ValueError('Duplicate or empty stable protein-coding genes')
    matrix = np.load(ROOT/'data_processed/expression_tpm.npy', mmap_mode='r')
    if matrix.dtype != np.float64:
        raise ValueError('Raw TPM precision changed')
    row_idx = data.expression_row_index.astype(int).to_numpy()
    if len(set(row_idx)) != len(row_idx) or row_idx.min() < 0 or row_idx.max() >= matrix.shape[0]:
        raise ValueError('Invalid expression-row mapping')
    X = np.empty((len(data), len(stable)+2), dtype=np.float64)
    X[:, :-2] = matrix[np.ix_(row_idx, gene_idx)]
    X[:, -2] = data.tissue.map(context['tissue_code']).to_numpy(dtype=np.float64)
    X[:, -1] = row_idx
    y = data.response_auc.to_numpy(dtype=np.float64)
    if not np.isfinite(X).all() or (X[:, :-2] < 0).any() or not np.isfinite(y).all():
        raise ValueError('Invalid raw development TPM or AUC')
    return data, X, y, stable


def make_reload_sentinels(X, tissue_code):
    """Prediction-only probes derived without an external expression/response."""
    if len(X) < 3:
        raise ValueError('At least three development rows required for reload sentinels')
    sentinels = X[:3].copy()
    sentinels[:, -1] = np.arange(int(np.max(X[:, -1]))+1,
                                 int(np.max(X[:, -1]))+4, dtype=np.float64)
    sentinels[1, -2] = max(tissue_code.values())+1
    sentinels[2, :-2] = 0.0  # exact pair ties; no missing or imputation
    return sentinels


def execute_fit_task(out, task, context, blocked_attempts=None):
    started = time.perf_counter()
    process = psutil.Process()
    cpu0 = process.cpu_times()
    samples = []
    stop = threading.Event()
    path = out/'checkpoints'/(task_stem(task)+'.json')
    artifact_path = out/'models'/(task_stem(task)+'.joblib')
    journal_path = fit_journal_path(out, task)
    if path.exists() or artifact_path.exists() or journal_path.exists():
        raise ValueError('Refusing to overwrite an existing final-model checkpoint, artifact or event journal')
    journal_path.parent.mkdir(exist_ok=True)
    def sample():
        while not stop.is_set():
            try:
                mem = process.memory_info()
                samples.append(dict(rss=int(mem.rss), private=int(getattr(mem, 'private', mem.vms)),
                                    **system_memory()))
            except psutil.Error:
                pass
            stop.wait(context['ordinal_config']['resources']['rss_sampling_interval_seconds'])
    sampler = threading.Thread(target=sample, daemon=True)
    sampler.start()
    obj = dict(schema_version='1.0', task=task, freeze_sha256=context['freeze_sha256'],
               status='failed', fit_count=0, fit_events=[], result=None, error=None,
               artifact_path='models/'+artifact_path.name, artifact_sha256=None)
    try:
        with journal_path.open('x', encoding='utf-8', newline='\n') as journal, \
                threadpool_limits(limits=1):
            journal.flush(); os.fsync(journal.fileno())
            def record_fit_event(event):
                journal.write(json.dumps(event, ensure_ascii=False, allow_nan=False,
                    sort_keys=True, separators=(',', ':'))+'\n')
                journal.flush()
                os.fsync(journal.fileno())
            data, X, y, stable = load_drug(out, task, context)
            manifest = pd.read_csv(out/'global_inner_fold_manifest.csv',
                dtype={'canonical_model_id': str, 'biological_group_id': str, 'tissue': str},
                keep_default_na=False)
            inner = local_inner_splits(data, manifest, context['external_config']['final_training']['inner_folds'])
            obj['development_canonical_ids'] = data.canonical_model_id.tolist()
            obj['development_groups'] = data.biological_group_id.tolist()
            obj['development_expression_row_indices'] = X[:, -1].astype(int).tolist()
            obj['inner_validation_expression_row_indices'] = [X[iv, -1].astype(int).tolist() for _, iv in inner]
            model, result = fit_full_development_model(
                X, y, data.biological_group_id.to_numpy(), inner, task['model'],
                context['ordinal_config'], stable, audit_sink=record_fit_event)
            journal_events, journal_error = read_fit_journal(journal_path)
            if journal_error is not None or journal_events != result['fit_events']:
                raise ValueError('Fsynced fit event journal differs from fitted model audit')
            sentinels = make_reload_sentinels(X, context['tissue_code'])
            sentinel_predictions = model.predict(sentinels)
            if not np.isfinite(sentinel_predictions).all():
                raise ValueError('Nonfinite synthetic/development-only reload sentinel prediction')
            result['reload_sentinel_predictions'] = sentinel_predictions.tolist()
            result['reload_sentinel_prediction_sha256'] = prediction_sha256(sentinel_predictions)
            obj['fit_events'] = result.pop('fit_events')
            obj['fit_count'] = result['fit_count']
            obj['result'] = result
            artifact = dict(schema_version='1.0', task=task, freeze_sha256=context['freeze_sha256'],
                            model=model, stable_gene_ids=stable, tissue_code=context['tissue_code'],
                            expression_genes_sha256=context['development_file_hashes']['data_processed/expression_genes.tsv'],
            expression_tpm_sha256=context['development_file_hashes']['data_processed/expression_tpm.npy'],
                            development_prediction_sha256=result['development_prediction_sha256'])
            _atomic_joblib(artifact_path, artifact, out)
            obj['artifact_sha256'] = sha(artifact_path)
            if blocked_attempts:
                raise ValueError('Forbidden external-response path attempted during final fit')
            obj['status'] = 'fitted_pending_verification'
    except BaseException:
        obj['error'] = traceback.format_exc()
    finally:
        _attach_fit_journal(obj, out, task)
        if obj['status'] != 'failed' and obj.get('fit_event_journal_parse_error') is not None:
            obj['status'] = 'failed'
            obj['error'] = 'Fsynced fit event journal contains an incomplete event'
        stop.set(); sampler.join()
        mem = process.memory_info(); cpu = process.cpu_times()
        if not samples:
            samples.append(dict(rss=int(mem.rss), private=int(getattr(mem, 'private', mem.vms)), **system_memory()))
        obj['resources'] = dict(wall_seconds=time.perf_counter()-started,
            cpu_seconds=cpu.user+cpu.system-cpu0.user-cpu0.system,
            sampled_peak_rss_bytes=max(s['rss'] for s in samples),
            sampled_peak_private_bytes=max(s['private'] for s in samples),
            minimum_available_ram_bytes=min(s['available_ram_bytes'] for s in samples),
            minimum_commit_free_bytes=min(s['system_commit_free_bytes'] for s in samples),
            windows_lifetime_peak_working_set_bytes=getattr(mem, 'peak_wset', None),
            pid=process.pid, controller_pid=process.ppid(), process_create_time=process.create_time(),
            preimport_threads_applied=EARLY_THREADS, threadpools=threadpool_info())
        obj['data_guard'] = dict(external_response_accessed=False,
                                 blocked_attempts=len(blocked_attempts or ()))
        save_checkpoint(path, obj)
    print(json.dumps(dict(event='final_fit_task', task=task, status=obj['status'],
                          fits=obj['fit_count'], artifact_sha256=obj['artifact_sha256'])), flush=True)
    return obj['status'] == 'fitted_pending_verification'


def execute_verify_task(out, task, context, blocked_attempts=None):
    """Executed in a distinct Python process after the fitting worker exits."""
    path = out/'checkpoints'/(task_stem(task)+'.json')
    obj = read_checkpoint(path, task, context['freeze_sha256'])
    if obj['status'] != 'fitted_pending_verification':
        raise ValueError('Fresh-process verifier requires a pending fitted artifact')
    started = time.perf_counter()
    try:
        data, X, y, stable = load_drug(out, task, context)
        artifact = joblib.load(out/obj['artifact_path'])  # trusted local SHA-bound artifact only
        from ordinal_diagnostic_models import DiagnosticRegressor
        if not isinstance(artifact, dict) or not isinstance(artifact.get('model'), DiagnosticRegressor):
            raise ValueError('Reloaded artifact is not a complete DiagnosticRegressor')
        if (artifact.get('task') != task or artifact.get('freeze_sha256') != context['freeze_sha256'] or
                artifact.get('stable_gene_ids') != stable or artifact.get('tissue_code') != context['tissue_code'] or
                artifact.get('expression_genes_sha256') != context['development_file_hashes']['data_processed/expression_genes.tsv'] or
                artifact.get('expression_tpm_sha256') != context['development_file_hashes']['data_processed/expression_tpm.npy']):
            raise ValueError('Reloaded artifact context mismatch')
        model = artifact['model']
        if set(model.fit_sample_ids_) != set(X[:, -1].astype(int)):
            raise ValueError('Reloaded model was not trained on the full drug development set')
        if data.canonical_model_id.tolist() != obj['development_canonical_ids']:
            raise ValueError('Reloaded development identity mismatch')
        if hashlib.sha256(np.ascontiguousarray(y, dtype='<f8').tobytes()).hexdigest() != obj['result']['development_response_sha256']:
            raise ValueError('Reloaded development response mismatch')
        predictions = model.predict(X)
        if prediction_sha256(predictions) != obj['result']['development_prediction_sha256']:
            raise ValueError('Fresh-process predictions differ from fitting-process predictions')
        sentinels = make_reload_sentinels(X, context['tissue_code'])
        sentinel_predictions = model.predict(sentinels)
        expected_sentinel = np.asarray(obj['result']['reload_sentinel_predictions'], dtype=np.float64)
        if not np.allclose(sentinel_predictions, expected_sentinel, atol=1e-12, rtol=0) or \
                prediction_sha256(sentinel_predictions) != obj['result']['reload_sentinel_prediction_sha256']:
            raise ValueError('Fresh-process sentinel predictions differ')
        audit = model.audit()
        old = obj['result']['full_fit_audit']
        for key in ('fitted_kind', 'fallback_reason', 'n_nonzero_molecular_coeffs',
                    'selected_gene_ids_sha256', 'built_feature_names_sha256'):
            if audit[key] != old[key]:
                raise ValueError('Reloaded final-model audit differs: '+key)
        if blocked_attempts:
            raise ValueError('Forbidden external-response path attempted during final-model reload')
        obj['fresh_process_reload_verified'] = True
        obj['reload_verification'] = dict(verified_utc=now_utc(), wall_seconds=time.perf_counter()-started,
            prediction_sha256=prediction_sha256(predictions),
            sentinel_prediction_sha256=prediction_sha256(sentinel_predictions),
            sentinel_policy='known tissue, unseen tissue, all-zero exact-tie sample', development_n=len(X),
            artifact_sha256=obj['artifact_sha256'], verifier_pid=os.getpid())
        obj['status'] = 'completed'
    except BaseException:
        obj['status'] = 'failed'
        obj['error'] = traceback.format_exc()
    save_checkpoint(path, obj)
    print(json.dumps(dict(event='final_reload_verification', task=task, status=obj['status'],
                          prediction_sha256=obj.get('reload_verification', {}).get('prediction_sha256'))), flush=True)
    return obj['status'] == 'completed'


def _child_command(out, flag, task_index, context_sha):
    return [sys.executable, '-u', str(Path(__file__).resolve()), '--run-dir', str(out),
            flag, str(task_index), '--context-sha', context_sha]


def _recovery_snapshot(out, meta):
    records = {}
    for path in sorted((out/'checkpoints').glob('*.json')):
        records[path.name] = sha(path)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    atomic_json(out/f'recovery_before_{stamp}.json', dict(metadata=meta, checkpoint_hashes=records,
        model_artifact_hashes={p.name: sha(p) for p in sorted((out/'models').glob('*.joblib'))}))


def _close_interrupted_session(out, meta):
    if not meta.get('sessions') or meta['sessions'][-1]['status'] != 'running':
        return
    old = meta['sessions'][-1]
    path = out/f'resource_samples_session{len(meta["sessions"])}.csv'
    last_elapsed = None
    if path.exists():
        with path.open('r', encoding='utf-8', newline='') as handle:
            for row in csv.DictReader(handle):
                if row.get('elapsed_seconds'):
                    last_elapsed = float(row['elapsed_seconds'])
    old.update(status='interrupted_unknown_cause', wall_seconds_lower_bound=last_elapsed,
               wall_seconds_is_lower_bound=True,
               interrupted_detected_utc=now_utc())


def prepare_resource_samples(out, session_number):
    """Keep a pre-metadata-crash sample file instead of truncating it."""
    path = Path(out)/f'resource_samples_session{session_number}.csv'
    _archive_stale_temp(path, out)
    return path


def _aggregate(out, tasks, freeze_sha):
    rows, fit_events = [], 0
    for task in tasks:
        obj = read_checkpoint(out/'checkpoints'/(task_stem(task)+'.json'), task, freeze_sha)
        if obj['status'] != 'completed':
            raise ValueError('Cannot aggregate an incomplete final-model task')
        fit_events += obj['fit_count']
        rows.append(dict(task=task, artifact_path=obj['artifact_path'], artifact_sha256=obj['artifact_sha256'],
            checkpoint_sha256=sha(out/'checkpoints'/(task_stem(task)+'.json')),
            fit_count=obj['fit_count'], chosen_params=obj['result']['chosen_params'],
            best_inner_rmse=obj['result']['best_inner_rmse'],
            n_nonzero_molecular_coeffs=obj['result']['full_fit_audit']['n_nonzero_molecular_coeffs'],
            fitted_kind=obj['result']['full_fit_audit']['fitted_kind'],
            fallback_reason=obj['result']['full_fit_audit']['fallback_reason'],
            selected_gene_ids_sha256=obj['result']['full_fit_audit']['selected_gene_ids_sha256'],
            development_prediction_sha256=obj['result']['development_prediction_sha256'],
            fresh_process_reload_verified=True))
    if len(rows) != 64 or fit_events != 1744:
        raise ValueError('64-artifact/1,744-fit final contract failed')
    body = ''.join(json.dumps(row, ensure_ascii=False, allow_nan=False)+'\n' for row in rows)
    _atomic_write_text(out/'model_manifest.jsonl', body, out)
    summary = dict(schema_version='1.0', artifacts=64, fits=fit_events,
        zero_molecular_artifacts=sum(row['n_nonzero_molecular_coeffs'] == 0
                                     for row in rows if row['n_nonzero_molecular_coeffs'] is not None),
        failed_artifacts=0, model_manifest_sha256=sha(out/'model_manifest.jsonl'),
        external_response_accessed=False)
    _atomic_write_text(out/'final_model_summary.json',
                       json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False)+'\n', out)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=DEFAULT_EXTERNAL_CONFIG)
    parser.add_argument('--ordinal-config', type=Path, default=DEFAULT_ORDINAL_CONFIG)
    parser.add_argument('--lock', type=Path, default=DEFAULT_LOCK)
    parser.add_argument('--run-dir', type=Path, default=DEFAULT_RUN)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--worker-task', type=int)
    parser.add_argument('--verify-task', type=int)
    parser.add_argument('--context-sha')
    args = parser.parse_args()
    for supplied, frozen, label in (
            (args.config, DEFAULT_EXTERNAL_CONFIG, '--config'),
            (args.ordinal_config, DEFAULT_ORDINAL_CONFIG, '--ordinal-config'),
            (args.lock, DEFAULT_LOCK, '--lock')):
        if supplied.resolve() != frozen.resolve():
            raise ValueError(f'{label} must use the fixed pre-fit protocol path: {frozen}')
    out = _checked_run_dir(args.run_dir)
    access, blocked = install_data_guard()
    if args.worker_task is not None or args.verify_task is not None:
        context = load_context(out, args.context_sha)
        verify_worker_freeze(context)
        tasks = task_list(context['external_config'], context['ordinal_config'])
        index = args.worker_task if args.worker_task is not None else args.verify_task
        if index is None or not 0 <= index < len(tasks):
            raise ValueError('Worker task index outside frozen panel')
        success = (execute_fit_task(out, tasks[index], context, blocked) if args.worker_task is not None else
                   execute_verify_task(out, tasks[index], context, blocked))
        if blocked:
            raise ValueError('Forbidden external-response read attempt')
        return 0 if success else 1

    external_path = args.config.resolve(); ordinal_path = args.ordinal_config.resolve()
    lock_path = args.lock.resolve()
    external = json.loads(external_path.read_text(encoding='utf-8'))
    ordinal = json.loads(ordinal_path.read_text(encoding='utf-8'))
    tasks = task_list(external, ordinal)
    lock, revision, source_hashes = verify_freeze(lock_path, external_path, ordinal_path)
    freeze_sha = sha(lock_path)
    verified_development = verify_development_inputs('v2')
    source, response_source_sha, formal_meta_sha = _response_source_verified()
    if out.exists() != args.resume:
        raise ValueError('Existing final run requires --resume; a new run must not resume')
    out.mkdir(parents=True, exist_ok=True)
    (out/'checkpoints').mkdir(exist_ok=True)
    (out/'models').mkdir(exist_ok=True)
    (out/'fit_event_journals').mkdir(exist_ok=True)
    meta_path = out/'final_model_metadata.json'
    if args.resume:
        meta = json.loads(meta_path.read_text(encoding='utf-8'))
        if meta['status'] == 'completed':
            raise ValueError('Completed final-model run is immutable')
        if meta['freeze_sha256'] != freeze_sha or meta['execution_commit'] != revision:
            raise ValueError('Final-model resume freeze/source mismatch')
        for name, expected in meta['prepared_input_hashes'].items():
            if sha(out/name) != expected:
                raise ValueError('Final-model prepared input changed: '+name)
        for name, expected in meta['checkpoint_hashes'].items():
            if sha(out/'checkpoints'/name) != expected:
                raise ValueError('Registered final checkpoint changed: '+name)
        _recovery_snapshot(out, meta)
        _close_interrupted_session(out, meta)
    else:
        if lock['contract'].get('final_model_tasks') != len(tasks):
            raise ValueError('Frozen final-model task count mismatch')
        shutil.copyfile(source, out/'development_responses.csv')
        if sha(out/'development_responses.csv') != response_source_sha:
            raise ValueError('Prepared development-response copy differs from formal source')
        responses = pd.read_csv(out/'development_responses.csv',
            dtype={'gdsc_drug_id': str, 'canonical_model_id': str,
                   'biological_group_id': str, 'tissue': str},
            keep_default_na=False, float_precision='round_trip')
        if len(responses) != 4967 or set(map(int, responses.gdsc_drug_id)) != set(external['drug_ids']):
            raise ValueError('Frozen fixed-eight development-response coverage changed')
        identities = pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv',
            dtype=str, keep_default_na=False).drop_duplicates('canonical_model_id')
        if len(identities) != 626 or set(identities.role) != {'development'}:
            raise ValueError('Frozen global development cohort changed')
        manifest, split_warnings = make_global_inner_manifest(identities,
            n_folds=external['final_training']['inner_folds'], seed=external['final_training']['inner_seed'])
        manifest.to_csv(out/'global_inner_fold_manifest.csv', index=False,
                        encoding='utf-8-sig', lineterminator='\n')
        validate_global_inner_manifest(identities,
            pd.read_csv(out/'global_inner_fold_manifest.csv', keep_default_na=False),
            external['final_training']['inner_folds'])
        prepared = {name: sha(out/name) for name in PREPARED_NAMES[:-1]}
        development_hashes = {name: sha(ROOT/name) for name in
            ('data_processed/expression_genes.tsv', 'data_processed/expression_tpm.npy')}
        context = dict(schema_version='1.0', freeze_sha256=freeze_sha,
            execution_commit=revision, source_hashes=source_hashes,
            environment=lock['environment'],
            prepared_input_hashes=prepared, development_file_hashes=development_hashes,
            external_config=external, ordinal_config=ordinal,
            tissue_code={name: i for i, name in enumerate(sorted(set(responses.tissue)))})
        atomic_json(out/'worker_context.json', context)
        prepared['worker_context.json'] = sha(out/'worker_context.json')
        meta = dict(schema_version='1.0', status='prepared', completed=False,
            created_utc=now_utc(), freeze_sha256=freeze_sha, execution_commit=revision,
            source_hashes=source_hashes, external_config_sha256=sha(external_path),
            ordinal_config_sha256=sha(ordinal_path),
            development_input_verification=verified_development,
            response_source_sha256=response_source_sha,
            formal_metadata_sha256=formal_meta_sha,
            expected_tasks=64, expected_fits=1744, completed_tasks=0,
            prepared_input_hashes=prepared, checkpoint_hashes={},
            split_warnings=split_warnings, external_response_accessed=False,
            environment=current_environment(),
            sessions=[])
        atomic_json(meta_path, meta)

    pending = []
    prepared_response = pd.read_csv(out/'development_responses.csv',
        dtype={'gdsc_drug_id': str, 'canonical_model_id': str,
               'biological_group_id': str, 'tissue': str},
        keep_default_na=False, float_precision='round_trip')
    prepared_manifest = pd.read_csv(out/'global_inner_fold_manifest.csv',
        dtype={'canonical_model_id': str, 'biological_group_id': str, 'tissue': str},
        keep_default_na=False)
    for index, task in enumerate(tasks):
        path = out/'checkpoints'/(task_stem(task)+'.json')
        if not path.exists():
            # An interrupted worker with any journal/artifact cannot be
            # silently refitted: its partial work is a preserved failure.
            recover_interrupted_task_if_needed(out, task, freeze_sha)
            pending.append((index, task, 'fit'))
            continue
        obj = read_checkpoint(path, task, freeze_sha)
        validate_prepared_checkpoint(obj, task, prepared_response, prepared_manifest,
                                     external['final_training']['inner_folds'])
        meta['checkpoint_hashes'][path.name] = sha(path)
        if obj['status'] == 'completed':
            continue
        if obj['status'] == 'fitted_pending_verification':
            pending.append((index, task, 'verify'))
        else:
            raise ValueError('Failed final-model task is preserved; inspect it before an explicit recovery')
    session = dict(started_utc=now_utc(), status='running',
                   pending_at_start=len(pending), completed_tasks_this_session=0,
                   refitted_completed_models=0, parallel_jobs=1,
                   resource_wait_seconds=0.0, child_wall_seconds_sum=0.0)
    meta['sessions'].append(session)
    meta.update(status='running', completed=False, completed_tasks=len(tasks)-len(pending))
    atomic_json(meta_path, meta)
    started = time.perf_counter()
    resource_path = prepare_resource_samples(out, len(meta['sessions']))
    fields = ['elapsed_seconds', 'process_tree_rss_bytes', 'process_tree_private_bytes',
              'process_tree_count', 'available_ram_bytes', 'system_commit_free_bytes',
              'system_commit_bytes', 'system_commit_limit_bytes', 'commit_measurement', 'inflight']
    tree_peak = private_peak = 0
    min_ram = system_memory()['available_ram_bytes']
    min_commit = system_memory()['system_commit_free_bytes']
    exit_code = 0
    child = None
    try:
        with resource_path.open('x', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            def record(inflight):
                nonlocal tree_peak, private_peak, min_ram, min_commit
                memory = process_tree_memory()
                tree_peak = max(tree_peak, memory['process_tree_rss_bytes'])
                private_peak = max(private_peak, memory['process_tree_private_bytes'])
                min_ram = min(min_ram, memory['available_ram_bytes'])
                min_commit = min(min_commit, memory['system_commit_free_bytes'])
                writer.writerow(dict(elapsed_seconds=time.perf_counter()-started,
                                     inflight=inflight, **memory))
                handle.flush()
                return memory
            for index, task, action in pending:
                wait_start = time.perf_counter()
                while True:
                    memory = record(0)
                    reasons = pressure_reasons(memory, ordinal['resources'])
                    if not reasons:
                        break
                    waited = time.perf_counter()-wait_start
                    print(json.dumps(dict(event='final_model_resource_wait', task=task,
                                          reasons=reasons, wait_seconds=waited)), flush=True)
                    if waited >= ordinal['resources']['pressure_max_wait_seconds']:
                        session['resource_wait_seconds'] += waited
                        session['status'] = meta['status'] = 'resource_paused'
                        exit_code = 75
                        break
                    time.sleep(ordinal['resources']['pressure_poll_seconds'])
                if exit_code:
                    break
                session['resource_wait_seconds'] += time.perf_counter()-wait_start
                if action == 'fit':
                    commands = [('--worker-task', 'fitted_pending_verification'),
                                ('--verify-task', 'completed')]
                else:
                    commands = [('--verify-task', 'completed')]
                for flag, expected_status in commands:
                    child_start = time.perf_counter()
                    child = subprocess.Popen(_child_command(out, flag, index,
                        meta['prepared_input_hashes']['worker_context.json']), cwd=ROOT,
                        env=dict(os.environ, PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8'))
                    while child.poll() is None:
                        record(1)
                        time.sleep(1)
                    session['child_wall_seconds_sum'] += time.perf_counter()-child_start
                    checkpoint_path = out/'checkpoints'/(task_stem(task)+'.json')
                    if not checkpoint_path.exists():
                        preserve_interrupted_task(out, task, freeze_sha,
                            f'Worker exited {child.returncode} without checkpoint; '
                            'completed fit events recovered from fsynced journal')
                    obj = read_checkpoint(checkpoint_path, task, freeze_sha)
                    meta['checkpoint_hashes'][checkpoint_path.name] = sha(checkpoint_path)
                    atomic_json(meta_path, meta)
                    if child.returncode != 0 or obj['status'] != expected_status:
                        raise RuntimeError('Final worker/verification failed: '+str(task)+': '+str(obj.get('error')))
                    child = None
                session['completed_tasks_this_session'] += 1
                meta['completed_tasks'] += 1
                atomic_json(meta_path, meta)
                print(json.dumps(dict(event='final_model_completed_task', task=task,
                                      completed_tasks=meta['completed_tasks'], total_tasks=len(tasks))), flush=True)
        if not exit_code:
            if blocked:
                raise ValueError('Forbidden external-response path attempted during final-model controller run')
            # An overnight run may span unrelated filesystem changes. Recheck
            # all frozen development/source inputs before declaring completion.
            verify_development_inputs('v2')
            verify_freeze(lock_path, external_path, ordinal_path)
            summary = _aggregate(out, tasks, freeze_sha)
            meta.update(status='completed', completed=True, completed_tasks=64,
                output_summary=summary,
                output_hashes={name: sha(out/name) for name in
                    ('model_manifest.jsonl', 'final_model_summary.json')},
                external_response_accessed=False)
            session['status'] = 'completed'
    except BaseException:
        if child is not None and child.poll() is None:
            child.wait()
        meta.update(status='failed', completed=False, failure=traceback.format_exc())
        session['status'] = 'failed'
        raise
    finally:
        session.update(finished_utc=now_utc(), wall_seconds=time.perf_counter()-started,
            sampled_process_tree_peak_rss_bytes=tree_peak,
            sampled_process_tree_peak_private_bytes=private_peak,
            minimum_available_ram_bytes=min_ram, minimum_commit_free_bytes=min_commit,
            rss_scope='Parent and one child RSS may double-count shared pages; sampled private/system commit also retained')
        meta['data_guard'] = dict(external_response_accessed=False,
                                  blocked_attempts=len(blocked), accessed_path_count=len(access))
        atomic_json(meta_path, meta)
    print(json.dumps(dict(event='final_model_run_finished', status=meta['status'],
                          completed_tasks=meta['completed_tasks'], output_summary=meta.get('output_summary'))), flush=True)
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
