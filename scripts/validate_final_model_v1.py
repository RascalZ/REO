"""Independent, response-blind audit of all frozen full-development models.

The parent validates the pre-fit lock, 64 checkpoint contracts and 1,744 fit
events without importing the final-model runner.  Each trusted SHA-verified
joblib artifact is then loaded in a separate Python process for schema and
prediction verification.  No CTRP file or response index is opened.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gc
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

import joblib
import numpy as np
import pandas as pd
import psutil
from sklearn.model_selection import StratifiedGroupKFold

from external_validation_v1 import names_sha256, validate_model_contract
from run_development_pilot import verify_development_inputs
from freeze_external_validation_v1_1 import install_pre_external_guard
from ordinal_repeated_runtime import sha

DEFAULT_RUN = ROOT/'runs/gdsc2-cancer-v2-fixed8-full-development-final-v1_1'
DEFAULT_LOCK = ROOT/'config/external_validation_v1_1.lock.json'
DEFAULT_EXTERNAL_CONFIG = ROOT/'config/external_validation_v1_1.json'
DEFAULT_AUDIT_NAME = 'final_model_independent_audit.json'


def current_environment():
    names = ('numpy', 'pandas', 'scipy', 'scikit-learn', 'joblib', 'threadpoolctl', 'psutil')
    return dict(python=platform.python_version(),
                packages={name: importlib.metadata.version(name) for name in names})


def _as_run(path):
    run = Path(path).resolve()
    if not run.is_relative_to((ROOT/'runs').resolve()) or run == (ROOT/'runs').resolve():
        raise ValueError('Audit run must be a child of repository runs/')
    return run


def _digest(obj):
    content = {k: v for k, v in obj.items() if k != 'payload_sha256'}
    return hashlib.sha256(json.dumps(content, sort_keys=True, ensure_ascii=False,
        allow_nan=False, separators=(',', ':')).encode('utf-8')).hexdigest()


def _response_hash(y):
    return hashlib.sha256(np.ascontiguousarray(y, dtype='<f8').tobytes()).hexdigest()


def _prediction_hash(values):
    arr = np.ascontiguousarray(values, dtype='<f8')
    if arr.ndim != 1 or not np.isfinite(arr).all():
        raise ValueError('Nonfinite or nonvector reloaded prediction')
    return hashlib.sha256(arr.tobytes()).hexdigest()


def _source_blob_hash(revision, name):
    return hashlib.sha256(subprocess.check_output(['git', 'show', revision+':'+name], cwd=ROOT)).hexdigest()


def verify_pre_fit_lock(path):
    lock = json.loads(Path(path).read_text(encoding='utf-8'))
    if lock.get('status') != 'frozen_pre_final_fit' or lock.get('external_response_accessed') is not False:
        raise ValueError('Pre-fit lock status/access contract failed')
    revision = subprocess.check_output(['git', 'rev-parse', lock['execution_commit']],
        cwd=ROOT, text=True).strip()
    if revision != lock['execution_commit']:
        raise ValueError('Pre-fit lock has noncanonical execution commit')
    files = lock['files']
    for name, expected in files.items():
        if sha(ROOT/name) != expected:
            raise ValueError('Pre-fit file SHA changed: '+name)
    for name in lock['execution_files']:
        if name not in files or _source_blob_hash(revision, name) != files[name]:
            raise ValueError('Pre-fit source not bound to immutable Git commit: '+name)
    contract = lock['contract']
    if contract['final_model_tasks'] != 64 or contract['final_fit_events'] != 1744 or \
            contract['development_groups'] != 626 or contract['development_drug_model_observations'] != 4967:
        raise ValueError('Pre-fit lock count contract changed')
    if lock.get('environment') != current_environment():
        raise ValueError('Python/package environment differs from pre-fit lock')
    verify_development_inputs('v2')
    return lock


def expected_tasks(config, ordinal):
    if config['drug_ids'] != ordinal['drug_ids'] or config['models'] != ordinal['models']:
        raise ValueError('P6 panel differs from repeated-evaluation panel')
    if len(config['drug_ids']) != 8 or len(config['models']) != 8 or \
            len(set(config['drug_ids'])) != 8 or len(set(config['models'])) != 8:
        raise ValueError('Fixed panel lacks exactly eight unique drugs and arms')
    if config['final_training']['inner_folds'] != 3 or config['final_training']['inner_seed'] != 2026101501:
        raise ValueError('Unexpected final inner split')
    grid = [dict(alpha_ratio=ratio, l1_ratio=l1)
            for ratio in ordinal['elastic_net']['alpha_ratios']
            for l1 in ordinal['elastic_net']['l1_ratios']]
    if len(grid) != 10:
        raise ValueError('Molecular inner candidate grid changed')
    tasks = []
    for drug in config['drug_ids']:
        for model in config['models']:
            if model == 'tissue_null':
                spec = dict(selection_rule='none', representation='tissue_null', n_genes=0)
            elif model == 'expression_strong':
                spec = dict(selection_rule='mad', representation='expression',
                            n_genes=int(ordinal['feature']['strong_expression_genes']))
            else:
                rule, representation = next(((r, model[len(r)+1:])
                    for r in ordinal['selection_rules'] if model.startswith(r+'_')), (None, None))
                if rule is None or representation not in ordinal['representations']:
                    raise ValueError('Unrecognized final-model arm: '+model)
                spec = dict(selection_rule=rule, representation=representation,
                            n_genes=int(ordinal['feature']['n_genes']))
            tasks.append(dict(drug_id=int(drug), model=model, **spec))
    return tasks, grid


def verify_global_manifest(identity, manifest, n_folds=3):
    identity = identity[['canonical_model_id', 'biological_group_id', 'tissue']].drop_duplicates()
    if len(identity) != 626 or identity.canonical_model_id.duplicated().any() or \
            identity.biological_group_id.duplicated().any():
        raise ValueError('Global identity cohort is not 626 unique biological groups')
    canonical = set(identity.canonical_model_id)
    if set(manifest.inner_fold) != set(range(n_folds)) or set(manifest.role) != {'train', 'validation'}:
        raise ValueError('Global inner manifest fold/role contract failed')
    reference = identity.set_index('canonical_model_id')
    for fold in range(n_folds):
        block = manifest.loc[manifest.inner_fold.eq(fold)]
        if len(block) != 626 or set(block.canonical_model_id) != canonical or \
                block.canonical_model_id.duplicated().any():
            raise ValueError('Global inner fold not a complete disjoint partition')
        train = block.loc[block.role.eq('train')]
        valid = block.loc[block.role.eq('validation')]
        if train.empty or valid.empty or set(train.biological_group_id) & set(valid.biological_group_id):
            raise ValueError('Global inner fold has empty side or group leakage')
        for row in block.itertuples():
            original = reference.loc[row.canonical_model_id]
            if row.biological_group_id != original.biological_group_id or row.tissue != original.tissue:
                raise ValueError('Global inner manifest identity metadata mismatch')
    counts = manifest.loc[manifest.role.eq('validation')].groupby('canonical_model_id').size()
    if counts.to_dict() != {name: 1 for name in canonical}:
        raise ValueError('Global inner validation membership is not exactly once')
    ordered = identity.sort_values('canonical_model_id', kind='stable').reset_index(drop=True)
    splitter = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=2026101501)
    for fold, (_, validation) in enumerate(splitter.split(
            np.zeros((len(ordered), 1)), ordered.tissue, ordered.biological_group_id)):
        expected = set(ordered.iloc[validation].canonical_model_id)
        actual = set(manifest.loc[manifest.inner_fold.eq(fold) &
                                  manifest.role.eq('validation'), 'canonical_model_id'])
        if actual != expected:
            raise ValueError('Global inner assignment differs from fixed seed 2026101501')


def _drug_design(data, manifest):
    positions = []
    for fold in range(3):
        block = manifest.loc[manifest.inner_fold.eq(fold)]
        validation_ids = set(block.loc[block.role.eq('validation'), 'canonical_model_id'])
        valid = np.flatnonzero(data.canonical_model_id.isin(validation_ids).to_numpy())
        train = np.flatnonzero(~data.canonical_model_id.isin(validation_ids).to_numpy())
        if not len(train) or not len(valid):
            raise ValueError('Drug inner intersection is empty')
        if set(data.iloc[train].biological_group_id) & set(data.iloc[valid].biological_group_id):
            raise ValueError('Drug inner biological-group leakage')
        positions.append((train, valid))
    if sorted(np.concatenate([iv for _, iv in positions]).tolist()) != list(range(len(data))):
        raise ValueError('Drug inner OOF coverage incomplete or duplicated')
    return positions


def verify_checkpoint(path, task, lock_sha, data, manifest, grid, run):
    obj = json.loads(path.read_text(encoding='utf-8'))
    if obj.get('payload_sha256') != _digest(obj) or obj.get('task') != task or \
            obj.get('freeze_sha256') != lock_sha or obj.get('status') != 'completed' or \
            obj.get('fresh_process_reload_verified') is not True:
        raise ValueError('Final checkpoint payload/status/freeze contract failed: '+path.name)
    journal_name = f"d{task['drug_id']}_{task['model']}.jsonl"
    journal_path = run/'fit_event_journals'/journal_name
    if obj.get('fit_event_journal_path') != 'fit_event_journals/'+journal_name or \
            not journal_path.is_file() or sha(journal_path) != obj.get('fit_event_journal_sha256'):
        raise ValueError('Final fit event journal path/SHA mismatch: '+journal_name)
    journal_events = []
    with journal_path.open('rb') as handle:
        for number, line in enumerate(handle, 1):
            if not line.endswith(b'\n'):
                raise ValueError(f'Final fit event journal line {number} is incomplete')
            event = json.loads(line.decode('utf-8'))
            if not isinstance(event, dict):
                raise ValueError(f'Final fit event journal line {number} is not an object')
            journal_events.append(event)
    if obj.get('fit_event_journal_parse_error') is not None or \
            journal_events != obj.get('fit_events') or len(journal_events) != obj.get('fit_count'):
        raise ValueError('Final fit event journal/checkpoint mismatch: '+journal_name)
    artifact_name = f"d{task['drug_id']}_{task['model']}.joblib"
    if obj.get('artifact_path') != 'models/'+artifact_name:
        raise ValueError('Final checkpoint artifact name differs from task')
    artifact = run/obj['artifact_path']
    if sha(artifact) != obj['artifact_sha256']:
        raise ValueError('Final artifact SHA mismatch: '+artifact_name)
    if data.canonical_model_id.duplicated().any() or data.biological_group_id.duplicated().any():
        raise ValueError('Drug input contains duplicate model/group')
    sample_ids = data.expression_row_index.astype(int).tolist()
    y = data.response_auc.to_numpy(dtype=np.float64)
    if obj['development_canonical_ids'] != data.canonical_model_id.tolist() or \
            obj['development_groups'] != data.biological_group_id.tolist() or \
            obj['development_expression_row_indices'] != sample_ids or \
            obj['result']['development_sample_ids'] != sample_ids or \
            obj['result']['development_response_sha256'] != _response_hash(y):
        raise ValueError('Final checkpoint differs from frozen development identity/label')
    splits = _drug_design(data, manifest)
    validation_ids = [data.iloc[iv].expression_row_index.astype(int).tolist() for _, iv in splits]
    if obj['inner_validation_expression_row_indices'] != validation_ids:
        raise ValueError('Final checkpoint inner fold differs from frozen global assignment')
    expected_fits = 1 if task['model'] == 'tissue_null' else 31
    events = obj['fit_events']
    if obj['fit_count'] != expected_fits or len(events) != expected_fits or \
            obj['result']['fit_count'] != expected_fits:
        raise ValueError('Final checkpoint fit count differs from frozen task budget')
    observed_contexts = Counter()
    for event in events:
        stage = event['stage']
        if event['declared_model'] != task['model'] or event['candidate_gene_budget'] != task['n_genes'] or \
                event['model_kind'] != task['representation'] or event['selection_rule'] != task['selection_rule']:
            raise ValueError('Fit event model specification mismatch')
        if stage == 'inner_training':
            candidate = event['candidate']; fold = event['inner_fold']
            if task['model'] == 'tissue_null' or not (0 <= candidate < 10 and 0 <= fold < 3) or \
                    event['params'] != grid[candidate]:
                raise ValueError('Unexpected inner candidate/fold/grid')
            train, valid = splits[fold]
            expected = data.iloc[train].expression_row_index.astype(int).tolist()
            if set(event['fit_sample_ids']) & set(data.iloc[valid].expression_row_index.astype(int)):
                raise ValueError('Validation sample entered training fit')
            observed_contexts[(candidate, fold)] += 1
        elif stage == 'full_development_refit':
            if event['candidate'] != -1 or event['inner_fold'] != -1 or \
                    event['params'] != obj['result']['chosen_params']:
                raise ValueError('Unexpected full-development refit parameters')
            expected = sample_ids
            observed_contexts[('refit', -1)] += 1
        else:
            raise ValueError('Unexpected final-model fit stage')
        if event['fit_sample_ids'] != expected or event['fit_response_sha256'] != \
                _response_hash(data.iloc[train].response_auc.to_numpy(dtype=np.float64)
                    if stage == 'inner_training' else y):
            raise ValueError('Fit event training membership/response hash mismatch')
        if event['training_encoded_matrix_shape'][0] != len(expected):
            raise ValueError('Fit event encoded feature matrix row count mismatch')
        selected = event.get('selected_gene_ids', [])
        if len(selected) > task['n_genes'] or selected != sorted(set(selected)):
            raise ValueError('Selected gene list exceeds budget or is not unique sorted')
        if selected and event.get('selected_gene_ids_sha256') != names_sha256(selected):
            raise ValueError('Fit event selected gene hash mismatch')
    required_contexts = {('refit', -1): 1}
    if task['model'] != 'tissue_null':
        required_contexts.update({(candidate, fold): 1 for candidate in range(10) for fold in range(3)})
        scores = obj['result']['inner_scores']
        if len(scores) != 10 or any(scores[i]['candidate'] != i or scores[i]['params'] != grid[i] or
                                    not np.isfinite(scores[i]['pooled_inner_rmse']) for i in range(10)):
            raise ValueError('Final inner score table/grid mismatch')
        first_best = min(scores, key=lambda record: record['pooled_inner_rmse'])
        if obj['result']['chosen_params'] != first_best['params'] or \
                obj['result']['best_inner_rmse'] != first_best['pooled_inner_rmse']:
            raise ValueError('Final selected parameter differs from first exact RMSE minimum')
    elif obj['result']['inner_scores'] != [] or obj['result']['chosen_params'] != {} or \
            obj['result']['best_inner_rmse'] is not None:
        raise ValueError('Tissue null unexpectedly tuned')
    if observed_contexts != required_contexts:
        raise ValueError('Missing or duplicated inner/full fit context')
    if obj['result']['full_fit_audit']['fit_sample_ids'] != sample_ids or \
            obj['data_guard'].get('external_response_accessed') is not False or \
            obj['data_guard'].get('blocked_attempts', 0) != 0:
        raise ValueError('Full final fit membership or external data guard failed')
    if obj['result']['development_prediction_sha256'] != obj['reload_verification']['prediction_sha256'] or \
            obj['result']['reload_sentinel_prediction_sha256'] != obj['reload_verification']['sentinel_prediction_sha256']:
        raise ValueError('Stored fresh-process reload hashes differ')
    return obj


def verify_artifact_child(run, task, checkpoint, context):
    """Load exactly one SHA-verified trusted local artifact in this process."""
    artifact_path = run/checkpoint['artifact_path']
    if sha(artifact_path) != checkpoint['artifact_sha256']:
        raise ValueError('Artifact changed before fresh-process load')
    artifact = joblib.load(artifact_path)
    model = artifact['model']
    if artifact.get('task') != task or artifact.get('freeze_sha256') != checkpoint['freeze_sha256'] or \
            artifact.get('tissue_code') != context['tissue_code'] or \
            artifact.get('expression_genes_sha256') != context['development_file_hashes']['data_processed/expression_genes.tsv'] or \
            artifact.get('expression_tpm_sha256') != context['development_file_hashes']['data_processed/expression_tpm.npy']:
        raise ValueError('Reloaded model provenance differs from frozen development')
    genes = pd.read_csv(ROOT/'data_processed/expression_genes.tsv', sep='\t',
        dtype=str, keep_default_na=False)
    positions = genes.index[genes.gene_type.eq('protein_coding')].to_numpy()
    stable = tuple(genes.iloc[positions].ensembl_id)
    if artifact.get('stable_gene_ids') != stable:
        raise ValueError('Reloaded model raw stable gene order differs')
    response = pd.read_csv(run/'development_responses.csv',
        dtype={'gdsc_drug_id': str, 'canonical_model_id': str,
               'biological_group_id': str, 'tissue': str},
        keep_default_na=False, float_precision='round_trip')
    data = response.loc[response.gdsc_drug_id.eq(str(task['drug_id']))].sort_values(
        'canonical_model_id', kind='stable').reset_index(drop=True)
    matrix = np.load(ROOT/'data_processed/expression_tpm.npy', mmap_mode='r')
    indices = data.expression_row_index.astype(int).to_numpy()
    X = np.empty((len(data), len(stable)+2), dtype=np.float64)
    X[:, :-2] = matrix[np.ix_(indices, positions)]
    X[:, -2] = data.tissue.map(context['tissue_code']).to_numpy(dtype=np.float64)
    X[:, -1] = indices
    del matrix
    if not np.isfinite(X).all() or (X[:, :-2] < 0).any():
        raise ValueError('Reload verification encountered invalid development TPM')
    audit = checkpoint['result']['full_fit_audit']
    selected = checkpoint['result']['feature_audit']
    validate_model_contract(model, representation=task['representation'], stable_gene_ids=stable,
        training_tissue_codes=audit['training_tissue_codes'],
        expected_gene_sha256=names_sha256(stable),
        expected_selected_gene_ids=[] if selected is None else selected['gene_ids'],
        expected_feature_names=[] if selected is None else selected['feature_names'],
        expected_feature_sha256=names_sha256([]) if selected is None else audit['built_feature_names_sha256'])
    if model.fit_sample_ids_ != indices.tolist() or model.fit_response_sha256_ != \
            checkpoint['result']['development_response_sha256']:
        raise ValueError('Reloaded model training membership/label hash mismatch')
    if model.feature_audit() != selected:
        raise ValueError('Reloaded model feature/coefficient audit differs from checkpoint')
    full_prediction = model.predict(X)
    if _prediction_hash(full_prediction) != checkpoint['result']['development_prediction_sha256']:
        raise ValueError('Independent full-development predictions differ after reload')
    sentinels = X[:3].copy()
    sentinels[:, -1] = np.arange(int(np.max(X[:, -1]))+1, int(np.max(X[:, -1]))+4)
    sentinels[1, -2] = max(context['tissue_code'].values())+1
    sentinels[2, :-2] = 0.0
    predicted = model.predict(sentinels)
    expected = np.asarray(checkpoint['result']['reload_sentinel_predictions'], dtype=np.float64)
    if not np.allclose(predicted, expected, atol=1e-12, rtol=0) or \
            _prediction_hash(predicted) != checkpoint['result']['reload_sentinel_prediction_sha256']:
        raise ValueError('Independent known/unknown/tie sentinel predictions differ after reload')
    return dict(task=task, artifact_sha256=checkpoint['artifact_sha256'],
        development_prediction_sha256=_prediction_hash(full_prediction),
        sentinel_prediction_sha256=_prediction_hash(predicted),
        selected_gene_ids_sha256=audit['selected_gene_ids_sha256'],
        n_nonzero_molecular_coeffs=audit['n_nonzero_molecular_coeffs'],
        fallback_reason=audit['fallback_reason'], fresh_process_pid=os.getpid(),
        rss_bytes_at_completion=psutil.Process().memory_info().rss)


def _write_audit(path, result):
    path = Path(path)
    if path.exists():
        previous = path.with_name(path.name+'.previous.'+sha(path))
        if not previous.exists():
            os.replace(path, previous)
    temp = path.with_name(path.name+'.tmp')
    if temp.exists():
        previous = path.with_name(temp.name+'.interrupted.'+sha(temp))
        os.replace(temp, previous)
    with temp.open('w', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n'); handle.flush(); os.fsync(handle.fileno())
    os.replace(temp, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=DEFAULT_RUN)
    parser.add_argument('--lock', type=Path, default=DEFAULT_LOCK)
    parser.add_argument('--artifact-task', type=int)
    parser.add_argument('--context-sha')
    args = parser.parse_args()
    run = _as_run(args.run_dir)
    access, blocked = install_pre_external_guard()
    started = time.perf_counter()
    lock_path = args.lock.resolve()
    lock_sha = sha(lock_path)
    external = json.loads(DEFAULT_EXTERNAL_CONFIG.read_text(encoding='utf-8'))
    ordinal = json.loads((ROOT/'config/ordinal_repeated_v1.json').read_text(encoding='utf-8'))
    tasks, grid = expected_tasks(external, ordinal)
    if args.artifact_task is not None:
        if not 0 <= args.artifact_task < len(tasks):
            raise ValueError('Artifact task index outside fixed panel')
        context_path = run/'worker_context.json'
        if not args.context_sha or sha(context_path) != args.context_sha:
            raise ValueError('Artifact child context SHA differs from parent preflight')
        context = json.loads(context_path.read_text(encoding='utf-8'))
        if context['freeze_sha256'] != lock_sha:
            raise ValueError('Artifact child context/lock mismatch')
        if context.get('environment') != current_environment():
            raise ValueError('Artifact child environment differs from pre-fit lock')
        task = tasks[args.artifact_task]
        path = run/'checkpoints'/f"d{task['drug_id']}_{task['model']}.json"
        checkpoint = json.loads(path.read_text(encoding='utf-8'))
        if checkpoint.get('payload_sha256') != _digest(checkpoint) or checkpoint['status'] != 'completed':
            raise ValueError('Artifact child checkpoint contract failed')
        result = verify_artifact_child(run, task, checkpoint, context)
        if blocked:
            raise ValueError('Forbidden external-response read attempt')
        print(json.dumps(result, ensure_ascii=False, allow_nan=False), flush=True)
        return 0

    audit_path = run/DEFAULT_AUDIT_NAME
    audit = dict(schema_version='1.0', status='failed', created_utc=datetime.now(timezone.utc).isoformat(),
        run_dir=str(run), pre_fit_lock_sha256=lock_sha, external_response_accessed=False,
        expected_tasks=64, expected_fit_events=1744, checkpoint_tasks_checked=0,
        fit_events_checked=0, fresh_process_artifacts_checked=0,
        artifact_results=[], error=None)
    try:
        lock = verify_pre_fit_lock(lock_path)
        metadata_path = run/'final_model_metadata.json'
        metadata_initial_sha = sha(metadata_path)
        metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
        if metadata['status'] != 'completed' or not metadata['completed'] or \
                metadata['freeze_sha256'] != lock_sha or \
                metadata['execution_commit'] != lock['execution_commit'] or \
                metadata['completed_tasks'] != 64 or metadata['expected_fits'] != 1744 or \
                metadata.get('external_response_accessed') is not False:
            raise ValueError('Final-model run metadata is incomplete or provenance differs')
        for name, expected in metadata['prepared_input_hashes'].items():
            if sha(run/name) != expected:
                raise ValueError('Final-model prepared input SHA mismatch: '+name)
        for name, expected in metadata['output_hashes'].items():
            if sha(run/name) != expected:
                raise ValueError('Final-model aggregate output SHA mismatch: '+name)
        manifest_rows = [json.loads(line) for line in
            (run/'model_manifest.jsonl').read_text(encoding='utf-8').splitlines() if line.strip()]
        if len(manifest_rows) != 64 or [row.get('task') for row in manifest_rows] != tasks:
            raise ValueError('Final-model manifest does not contain the fixed 64 tasks in order')
        summary = json.loads((run/'final_model_summary.json').read_text(encoding='utf-8'))
        if summary.get('artifacts') != 64 or summary.get('fits') != 1744 or \
                summary.get('model_manifest_sha256') != sha(run/'model_manifest.jsonl') or \
                summary.get('external_response_accessed') is not False:
            raise ValueError('Final-model summary/manifest count and SHA contract failed')
        expected_stems = {f"d{task['drug_id']}_{task['model']}" for task in tasks}
        checkpoint_stems = {path.stem for path in (run/'checkpoints').glob('*.json')}
        artifact_stems = {path.stem for path in (run/'models').glob('*.joblib')}
        if checkpoint_stems != expected_stems or artifact_stems != expected_stems:
            raise ValueError('Final-model checkpoint/artifact directory has missing or extra task files')
        if any(path.name.endswith('.tmp') for folder in (run, run/'checkpoints', run/'models')
               for path in folder.iterdir() if path.is_file()):
            raise ValueError('Unresolved temporary final-model file remains outside incident archive')
        source = 'runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1/development_responses.csv'
        if sha(run/'development_responses.csv') != lock['files'][source]:
            raise ValueError('Final-model prepared responses differ from frozen source')
        context = json.loads((run/'worker_context.json').read_text(encoding='utf-8'))
        if context['freeze_sha256'] != lock_sha or context['external_config'] != external or \
                context['ordinal_config'] != ordinal or \
                context.get('environment') != lock['environment']:
            raise ValueError('Final-model context does not bind frozen configuration')
        identity = pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv',
            dtype=str, keep_default_na=False)
        manifest = pd.read_csv(run/'global_inner_fold_manifest.csv',
            dtype={'canonical_model_id': str, 'biological_group_id': str, 'tissue': str},
            keep_default_na=False)
        verify_global_manifest(identity, manifest)
        response = pd.read_csv(run/'development_responses.csv',
            dtype={'gdsc_drug_id': str, 'canonical_model_id': str,
                   'biological_group_id': str, 'tissue': str},
            keep_default_na=False, float_precision='round_trip')
        if len(response) != 4967 or set(map(int, response.gdsc_drug_id)) != set(external['drug_ids']):
            raise ValueError('Final-model response coverage changed')
        if response.canonical_model_id.nunique() > 626 or \
                not set(response.canonical_model_id).issubset(set(identity.canonical_model_id)):
            raise ValueError('Final-model responses do not belong to frozen development cohort')
        zero = fallback = 0
        for index, task in enumerate(tasks):
            data = response.loc[response.gdsc_drug_id.eq(str(task['drug_id']))].sort_values(
                'canonical_model_id', kind='stable').reset_index(drop=True)
            path = run/'checkpoints'/f"d{task['drug_id']}_{task['model']}.json"
            checkpoint = verify_checkpoint(path, task, lock_sha, data, manifest, grid, run)
            if metadata['checkpoint_hashes'].get(path.name) != sha(path):
                raise ValueError('Final-model registered checkpoint SHA mismatch')
            row = manifest_rows[index]
            if row.get('artifact_path') != checkpoint['artifact_path'] or \
                    row.get('artifact_sha256') != checkpoint['artifact_sha256'] or \
                    row.get('checkpoint_sha256') != sha(path) or \
                    row.get('fit_count') != checkpoint['fit_count'] or \
                    row.get('development_prediction_sha256') != checkpoint['result']['development_prediction_sha256']:
                raise ValueError('Final-model aggregate manifest row differs from checkpoint')
            audit['checkpoint_tasks_checked'] += 1
            audit['fit_events_checked'] += checkpoint['fit_count']
            core = checkpoint['result']['full_fit_audit']
            zero += int(core['n_nonzero_molecular_coeffs'] == 0)
            fallback += int(bool(core['fallback_reason']))
            result = subprocess.run([sys.executable, '-u', str(Path(__file__).resolve()),
                '--run-dir', str(run), '--lock', str(lock_path), '--artifact-task', str(index),
                '--context-sha', metadata['prepared_input_hashes']['worker_context.json']],
                cwd=ROOT, capture_output=True, text=True)
            if result.returncode != 0:
                raise ValueError('Fresh-process artifact audit failed '+path.name+': '+result.stderr[-4000:])
            payload = json.loads(result.stdout.strip().splitlines()[-1])
            if payload['task'] != task or payload['artifact_sha256'] != checkpoint['artifact_sha256']:
                raise ValueError('Fresh-process artifact audit identity/hash mismatch')
            audit['artifact_results'].append(payload)
            audit['fresh_process_artifacts_checked'] += 1
            gc.collect()
        if audit['checkpoint_tasks_checked'] != 64 or audit['fit_events_checked'] != 1744 or \
                audit['fresh_process_artifacts_checked'] != 64:
            raise ValueError('Final-model 64/1744/64 audit contract failed')
        # Close the TOCTOU window of a long 64-process audit before sealing it.
        verify_pre_fit_lock(lock_path)
        if sha(metadata_path) != metadata_initial_sha or \
                sha(run/'model_manifest.jsonl') != metadata['output_hashes']['model_manifest.jsonl']:
            raise ValueError('Final-model metadata/manifest changed during independent audit')
        for payload in audit['artifact_results']:
            task = payload['task']
            artifact = run/'models'/f"d{task['drug_id']}_{task['model']}.joblib"
            if sha(artifact) != payload['artifact_sha256']:
                raise ValueError('Final artifact changed during independent audit')
        if blocked:
            raise ValueError('Forbidden CTRP response path attempted during audit')
        audit.update(status='completed_contract_validated', artifacts=64, fit_events=1744,
            zero_molecular_artifacts=zero,
            fallback_artifacts=fallback, final_model_metadata_sha256=sha(metadata_path),
            model_manifest_sha256=sha(run/'model_manifest.jsonl'),
            source_execution_commit=lock['execution_commit'],
            development_group_count=626, development_drug_model_observations=4967,
            external_response_accessed=False,
            wall_seconds=time.perf_counter()-started,
            sampled_parent_rss_bytes=psutil.Process().memory_info().rss)
    except BaseException:
        audit['error'] = traceback.format_exc()
        audit['wall_seconds'] = time.perf_counter()-started
        _write_audit(audit_path, audit)
        raise
    _write_audit(audit_path, audit)
    print(json.dumps(dict(status=audit['status'], tasks=audit['checkpoint_tasks_checked'],
        fits=audit['fit_events_checked'], reloaded=audit['fresh_process_artifacts_checked'],
        zero_molecular_artifacts=audit['zero_molecular_artifacts'],
        external_response_accessed=False, audit_sha256=sha(audit_path)), ensure_ascii=False), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
