"""Hard-gated P6 response-free preflight and one-time CTRPv2 evaluation.

``--preflight`` verifies the frozen v2 identity roster and 64 saved models,
then writes predictions/abstentions without opening CTRP response or response
index. ``--execute`` refuses to read a response unless the *second* immutable
execution lock binds that preflight, all models, source code and the first lock.
The latter mode is intentionally never invoked during protocol preparation.
"""
from __future__ import annotations

from ordinal_repeated_runtime import set_preimport_threads
set_preimport_threads()

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gc
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback

import joblib
import numpy as np
import pandas as pd
import psutil

from external_validation_v1 import (
    bootstrap_paired_rho, make_shared_bootstrap, paired_common_table,
    paired_rho_difference, predict_external_arm, prepare_gateway,
    spearman_checked,
)
from ordinal_repeated_models import model_spec
from ordinal_repeated_runtime import sha
from run_development_pilot import ROOT
from run_final_model_v1 import read_checkpoint, task_list, task_stem, verify_freeze

sys.path.insert(0, str(ROOT/'scripts'))
from freeze_external_validation_v1_1 import install_pre_external_guard


FINAL_RUN = ROOT / 'runs/gdsc2-cancer-v2-fixed8-full-development-final-v1_1'
EXTERNAL_RUN = ROOT / 'runs/gdsc2-cancer-v2-fixed8-external-v1_1'
FIRST_LOCK = ROOT / 'config/external_validation_v1_1.lock.json'
SECOND_LOCK = ROOT / 'config/external_execution_v1_1.lock.json'
EXTERNAL_CONFIG = ROOT / 'config/external_validation_v1_1.json'
ORDINAL_CONFIG = ROOT / 'config/ordinal_repeated_v1.json'
EXTERNAL_MANIFEST = ROOT / 'results/cohort_v2/external_primary_sample_manifest.csv'
DEVELOPMENT_MANIFEST = ROOT / 'results/cohort_v2/development_sample_manifest.csv'
GLOBAL_SPLIT = ROOT / 'results/cohort_v2/global_split_manifest.csv'
MAIN_DRUGS = ROOT / 'results/cohort_v2/main_drugs.csv'
EXPRESSION_GENES = ROOT / 'data_processed/expression_genes.tsv'
EXPRESSION_SAMPLES = ROOT / 'data_processed/expression_samples.tsv'
EXPRESSION_TPM = ROOT / 'data_processed/expression_tpm.npy'
RAW_RESPONSE_NAME = 'data_raw/ctrp/v20.data.curves_post_qc.txt'
PREDICTION_FIELDS = (
    'drug_id', 'model', 'canonical_model_id', 'biological_group_id', 'tissue',
    'eligible_denominator_member', 'predicted_auc', 'status', 'abstain_reason',
    'unknown_training_tissue', 'zero_molecular_solution', 'artifact_sha256',
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def resource_snapshot(started: float) -> dict:
    info = psutil.Process(os.getpid()).memory_info()
    return {'wall_seconds': time.perf_counter() - started,
            'process_rss_bytes': int(info.rss),
            'windows_lifetime_peak_working_set_bytes':
                int(info.peak_wset) if getattr(info, 'peak_wset', None) is not None else None}


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8'))


def _same_file_hashes(files: dict[str, str], names: tuple[str, ...]) -> None:
    for name in names:
        path = (ROOT / name).resolve()
        path.relative_to(ROOT.resolve())
        if files.get(name) != sha(path):
            raise ValueError('Frozen source/input hash mismatch: ' + name)


def parse_pipe_numbers(raw: str, *, field: str) -> tuple[int, ...]:
    parts = str(raw).split('|')
    if not parts or any(not part.isdecimal() for part in parts):
        raise ValueError('Invalid frozen ' + field + ' pointer list')
    values = tuple(map(int, parts))
    if len(set(values)) != len(values):
        raise ValueError('Duplicate frozen ' + field + ' pointer')
    return values


def validate_frozen_roster(external: pd.DataFrame, development: pd.DataFrame,
                           global_split: pd.DataFrame, main_drugs: pd.DataFrame,
                           expression_samples: pd.DataFrame,
                           config: dict) -> dict:
    """No-response check of panel, group isolation and exact source pointers."""
    expected = {'gdsc_drug_id', 'ctrp_compound_id', 'canonical_model_id',
                'biological_group_id', 'ccle_id', 'expression_row_index', 'tissue',
                'ctrp_experiment_ids', 'ctrp_source_row_numbers', 'n_source_experiments',
                'identity_qc_statuses', 'future_response_aggregation', 'role',
                'cohort_version'}
    if not expected.issubset(external):
        raise ValueError('Frozen external roster has missing identity/pointer columns')
    drugs = tuple(map(int, config['drug_ids']))
    models = tuple(config['models'])
    if len(drugs) != 8 or len(set(drugs)) != 8 or len(models) != 8 or len(set(models)) != 8:
        raise ValueError('External panel is not exactly eight drugs and eight arms')
    if len(external) != 1185 or external.biological_group_id.nunique() != 175:
        raise ValueError('Frozen v2 external roster count changed')
    if external.duplicated(['gdsc_drug_id', 'canonical_model_id']).any():
        raise ValueError('Duplicate drug/canonical ID in external roster')
    if external[['canonical_model_id', 'biological_group_id']].drop_duplicates().canonical_model_id.nunique() != 175:
        raise ValueError('External canonical model universe changed')
    by_id = external.groupby('canonical_model_id', sort=False)
    if any(block.biological_group_id.nunique() != 1 or block.expression_row_index.nunique() != 1 or
           block.tissue.nunique() != 1 or block.ccle_id.nunique() != 1 for _, block in by_id):
        raise ValueError('Conflicting external canonical identity across drugs')
    if external[['canonical_model_id', 'biological_group_id']].drop_duplicates().biological_group_id.duplicated().any():
        raise ValueError('External biological group maps to multiple canonical models')
    if set(external.gdsc_drug_id.astype(int)) != set(drugs):
        raise ValueError('Frozen external drug panel differs from protocol')
    panel = main_drugs.sort_values('primary_drug_rank', kind='stable')
    if panel.gdsc_drug_id.astype(int).tolist() != list(drugs):
        raise ValueError('Main-drug order differs from frozen protocol')
    compounds = dict(zip(panel.gdsc_drug_id.astype(int), panel.ctrp_compound_id.astype(int)))
    for drug, block in external.groupby(external.gdsc_drug_id.astype(int)):
        if not block.ctrp_compound_id.astype(int).eq(compounds[drug]).all():
            raise ValueError('CTRP compound ID differs from frozen eight-drug map')
        if block.biological_group_id.duplicated().any():
            raise ValueError('Duplicate biological group within an external drug')
    if not (external.role.eq('external_primary').all() and external.cohort_version.eq('v2').all() and
            external.identity_qc_statuses.eq('SNP-matched-reference').all() and
            external.future_response_aggregation.eq('median_of_included_experiments').all()):
        raise ValueError('External role/QC/aggregation policy differs from frozen v2')
    external_ids = set(external.canonical_model_id)
    external_groups = set(external.biological_group_id)
    if external_ids & set(development.canonical_model_id) or external_groups & set(development.biological_group_id):
        raise ValueError('Development and external model/group overlap')
    candidate = global_split[global_split.global_role.eq('GDSC_development_pool_never_external')]
    if external_ids & set(candidate.canonical_model_id) or external_groups & set(candidate.biological_group_id):
        raise ValueError('External model/group overlaps known GDSC candidate group')
    registered = global_split.set_index('canonical_model_id')
    if registered.index.duplicated().any():
        raise ValueError('Global split has duplicate canonical identity')
    for row in external[['canonical_model_id', 'biological_group_id', 'expression_row_index', 'tissue', 'ccle_id']].drop_duplicates().itertuples(index=False):
        if row.canonical_model_id not in registered.index:
            raise ValueError('External model absent from frozen global split')
        known = registered.loc[row.canonical_model_id]
        if (known.biological_group_id != row.biological_group_id or
                int(known.expression_row_index) != int(row.expression_row_index) or
                known.tissue != row.tissue or known.ccle_id != row.ccle_id or
                str(known.in_selected_external_primary).lower() != 'true'):
            raise ValueError('External identity differs from frozen global split')
    samples = expression_samples.copy()
    samples['expression_row_index'] = pd.to_numeric(samples.expression_row_index, errors='raise').astype(int)
    samples = samples.set_index('expression_row_index')
    if samples.index.duplicated().any():
        raise ValueError('Duplicate expression row index')
    for row in external[['canonical_model_id', 'ccle_id', 'expression_row_index']].drop_duplicates().itertuples(index=False):
        idx = int(row.expression_row_index)
        if idx not in samples.index or samples.loc[idx].depmap_id != row.canonical_model_id or samples.loc[idx].ccle_id != row.ccle_id:
            raise ValueError('External identity/expression row mapping mismatch')
    pointers: dict[int, tuple[int, str]] = {}
    count = 0
    for row in external.itertuples(index=False):
        lines = parse_pipe_numbers(row.ctrp_source_row_numbers, field='source row')
        experiments = parse_pipe_numbers(row.ctrp_experiment_ids, field='experiment')
        if len(lines) != len(experiments) or len(lines) != int(row.n_source_experiments):
            raise ValueError('Source pointer/experiment count mismatch')
        if any(number < 2 for number in lines):
            raise ValueError('Source pointer must include the header offset')
        for line in lines:
            if line in pointers:
                raise ValueError('Duplicate source row pointer across external roster')
            pointers[line] = (int(row.gdsc_drug_id), row.canonical_model_id)
        count += len(lines)
    if count != 1201:
        raise ValueError('Frozen source-row pointer count changed')
    return {'external_unique_models': len(external_ids), 'external_unique_groups': len(external_groups),
            'eligible_drug_model_rows': len(external), 'selected_source_rows': count,
            'development_model_overlap': 0, 'development_group_overlap': 0,
            'known_candidate_group_overlap': 0,
            'per_drug_eligible': {str(drug): len(external[external.gdsc_drug_id.astype(int).eq(drug)]) for drug in drugs}}


def frozen_gene_schema(genes_path: Path, matrix_path: Path) -> tuple[tuple[str, ...], np.ndarray, np.ndarray]:
    genes = pd.read_csv(genes_path, sep='\t', dtype=str, keep_default_na=False)
    protein_idx = genes.index[genes.gene_type.eq('protein_coding')].to_numpy(dtype=int)
    stable = tuple(genes.iloc[protein_idx].ensembl_id)
    if len(stable) == 0 or len(set(stable)) != len(stable):
        raise ValueError('Frozen protein-coding stable-gene schema is invalid')
    if not np.array_equal(pd.to_numeric(genes.expression_column_index).to_numpy(dtype=int), np.arange(len(genes))):
        raise ValueError('Expression gene column index/order mismatch')
    matrix = np.load(matrix_path, mmap_mode='r')
    if matrix.dtype != np.dtype('float64') or matrix.ndim != 2 or matrix.shape[1] != len(genes):
        raise ValueError('Frozen expression TPM matrix dtype/shape mismatch')
    return stable, protein_idx, matrix


def _load_manifest(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={'gdsc_drug_id': str, 'ctrp_compound_id': str,
                                    'canonical_model_id': str, 'biological_group_id': str,
                                    'ccle_id': str, 'tissue': str}, keep_default_na=False)


def _source_hashes(first: dict) -> dict[str, str]:
    names = first['execution_files']
    required = {'src/run_external_validation_v1.py', 'src/external_validation_v1.py'}
    if not required.issubset(names):
        raise ValueError('First P6 lock does not bind the external scoring/statistics source')
    return {name: first['files'][name] for name in names}


def _check_final_run(first_sha: str, config: dict, ordinal: dict) -> tuple[dict, list[dict], dict]:
    meta_path = FINAL_RUN/'final_model_metadata.json'
    manifest_path = FINAL_RUN/'model_manifest.jsonl'
    meta = _json(meta_path)
    summary = _json(FINAL_RUN/'final_model_summary.json')
    audit = _json(FINAL_RUN/'final_model_independent_audit.json')
    if (meta.get('status'), meta.get('completed'), meta.get('completed_tasks'), meta.get('expected_fits')) != ('completed', True, 64, 1744):
        raise ValueError('Full-development final-model run is incomplete')
    if meta.get('freeze_sha256') != first_sha or meta.get('external_response_accessed') is not False:
        raise ValueError('Final-model run lock/response contract mismatch')
    if (audit.get('status'), audit.get('artifacts'), audit.get('fit_events'),
            audit.get('pre_fit_lock_sha256'), audit.get('model_manifest_sha256'),
            audit.get('external_response_accessed')) != (
            'completed_contract_validated', 64, 1744, first_sha, sha(manifest_path), False):
        raise ValueError('Independent final-model audit has not passed')
    if (summary.get('artifacts'), summary.get('fits'), summary.get('model_manifest_sha256')) != (64, 1744, sha(manifest_path)):
        raise ValueError('Final-model summary differs from frozen manifest')
    records = [json.loads(line) for line in manifest_path.read_text(encoding='utf-8').splitlines() if line]
    tasks = task_list(config, ordinal)
    by_task = {(int(r['task']['drug_id']), r['task']['model']): r for r in records}
    if len(records) != 64 or len(by_task) != 64 or set(by_task) != {(task['drug_id'], task['model']) for task in tasks}:
        raise ValueError('Final-model manifest omits or duplicates a frozen drug/arm')
    model_hashes = {}
    for task in tasks:
        record = by_task[(task['drug_id'], task['model'])]
        if record['task'] != task or not record.get('fresh_process_reload_verified'):
            raise ValueError('Final-model task/reload manifest mismatch')
        artifact_path = FINAL_RUN/record['artifact_path']
        artifact_path.resolve().relative_to((FINAL_RUN/'models').resolve())
        if record['artifact_path'] != 'models/'+task_stem(task)+'.joblib' or sha(artifact_path) != record['artifact_sha256']:
            raise ValueError('Final-model artifact SHA/task mismatch')
        key = artifact_path.relative_to(ROOT).as_posix()
        model_hashes[key] = record['artifact_sha256']
        checkpoint = read_checkpoint(FINAL_RUN/'checkpoints'/(task_stem(task)+'.json'), task, first_sha)
        if checkpoint['status'] != 'completed' or checkpoint['artifact_sha256'] != record['artifact_sha256']:
            raise ValueError('Final-model checkpoint/manifest mismatch')
        if meta.get('checkpoint_hashes', {}).get(task_stem(task)+'.json') != sha(FINAL_RUN/'checkpoints'/(task_stem(task)+'.json')):
            raise ValueError('Final-model checkpoint registered SHA mismatch')
    return meta, records, model_hashes


def preflight() -> dict:
    started = time.perf_counter()
    _, blocked = install_pre_external_guard()
    first, _, _ = verify_freeze(FIRST_LOCK, EXTERNAL_CONFIG, ORDINAL_CONFIG)
    first_sha = sha(FIRST_LOCK)
    config = _json(EXTERNAL_CONFIG)
    ordinal = _json(ORDINAL_CONFIG)
    if config.get('protocol_status') != 'frozen_pre_final_fit' or config.get('pre_fit_gates_passed') is not True:
        raise ValueError('P6 scientific protocol is not frozen before external preflight')
    source_hashes = _source_hashes(first)
    _same_file_hashes(first['files'], tuple(source_hashes))
    final_meta, records, model_hashes = _check_final_run(first_sha, config, ordinal)
    if EXTERNAL_RUN.exists() and any(EXTERNAL_RUN.iterdir()):
        raise FileExistsError('External preflight output already exists; preserve it for audit')
    external = _load_manifest(EXTERNAL_MANIFEST)
    development = _load_manifest(DEVELOPMENT_MANIFEST)
    global_split = pd.read_csv(GLOBAL_SPLIT, dtype=str, keep_default_na=False)
    main_drugs = pd.read_csv(MAIN_DRUGS, keep_default_na=False)
    samples = pd.read_csv(EXPRESSION_SAMPLES, sep='\t', dtype=str, keep_default_na=False)
    identity_audit = validate_frozen_roster(external, development, global_split,
                                            main_drugs, samples, config)
    stable, gene_idx, matrix = frozen_gene_schema(EXPRESSION_GENES, EXPRESSION_TPM)
    expression_genes_sha = sha(EXPRESSION_GENES)
    expression_tpm_sha = sha(EXPRESSION_TPM)
    global_tissue_code = _json(FINAL_RUN/'worker_context.json')['tissue_code']
    by_task = {(int(r['task']['drug_id']), r['task']['model']): r for r in records}
    development_tissues = {drug: set(development.loc[development.gdsc_drug_id.eq(str(drug)), 'tissue'])
                           for drug in config['drug_ids']}
    forbidden_ids = set(development.canonical_model_id)
    forbidden_groups = set(development.biological_group_id)
    candidate = global_split[global_split.global_role.eq('GDSC_development_pool_never_external')]
    forbidden_ids.update(candidate.canonical_model_id)
    forbidden_groups.update(candidate.biological_group_id)
    predictions = []
    coverage = []
    for drug in config['drug_ids']:
        block = external.loc[external.gdsc_drug_id.eq(str(drug))].sort_values('canonical_model_id', kind='stable')
        rows = block.expression_row_index.astype(int).to_numpy()
        if rows.min() < 0 or rows.max() >= matrix.shape[0] or len(set(rows)) != len(rows):
            raise ValueError('External expression-row index is out of range or duplicated')
        identity = block[['canonical_model_id', 'biological_group_id', 'tissue']].reset_index(drop=True)
        tpm = pd.DataFrame(matrix[np.ix_(rows, gene_idx)], index=identity.canonical_model_id,
                           columns=stable, dtype=np.float64)
        local_tissue_map = {name: global_tissue_code[name] for name in development_tissues[drug]}
        trained_codes = tuple(sorted(local_tissue_map.values()))
        batch = prepare_gateway(identity, tpm, stable_gene_ids=stable,
                                tissue_to_code=local_tissue_map,
                                training_tissue_codes=trained_codes,
                                forbidden_training_canonical_ids=forbidden_ids,
                                forbidden_training_group_ids=forbidden_groups)
        for arm in config['models']:
            record = by_task[(drug, arm)]
            artifact_path = FINAL_RUN/record['artifact_path']
            artifact = joblib.load(artifact_path)  # local SHA checked above
            task = record['task']
            if (artifact.get('task') != task or artifact.get('freeze_sha256') != first_sha or
                    tuple(artifact.get('stable_gene_ids', ())) != stable or
                    artifact.get('tissue_code') != global_tissue_code or
                    artifact.get('expression_genes_sha256') != expression_genes_sha or
                    artifact.get('expression_tpm_sha256') != expression_tpm_sha):
                raise ValueError('Reloaded final-model artifact source/identity mismatch')
            model = artifact['model']
            checkpoint = read_checkpoint(FINAL_RUN/'checkpoints'/(task_stem(task)+'.json'), task, first_sha)
            feature_audit = checkpoint['result']['feature_audit']
            selected = feature_audit['gene_ids'] if feature_audit is not None else []
            features = feature_audit['feature_names'] if feature_audit is not None else []
            out = predict_external_arm(model, batch, drug_id=drug, arm=arm,
                representation=model_spec(arm, ordinal)['representation'],
                expected_selected_gene_ids=selected, expected_feature_names=features,
                artifact_sha256=record['artifact_sha256'])
            predictions.append(out)
            coverage.append({'drug_id': drug, 'model': arm, 'eligible': len(block),
                             'predicted': int(out.status.eq('predicted').sum()),
                             'abstained': int(out.status.eq('abstained').sum()),
                             'unknown_tissue': int(out.unknown_training_tissue.sum()),
                             'zero_molecular_solution': int(out.zero_molecular_solution.astype(str).str.lower().eq('true').sum()),
                             'noncoverage_reasons': dict(Counter(
                                 reason for value in out.loc[out.status.eq('abstained'), 'abstain_reason']
                                 for reason in str(value).split(';') if reason))})
            del artifact, model, checkpoint, out
            gc.collect()
    result = pd.concat(predictions, ignore_index=True)
    if len(result) != 9480 or result.duplicated(['drug_id', 'model', 'canonical_model_id']).any():
        raise ValueError('Response-free 1,185 x 8 prediction coverage contract failed')
    if blocked:
        raise ValueError('Preflight guard blocked an attempted external-response access')
    EXTERNAL_RUN.mkdir(parents=True, exist_ok=True)
    prediction_path = EXTERNAL_RUN/'prediction_manifest.csv'
    result.loc[:, list(PREDICTION_FIELDS)].to_csv(prediction_path, index=False,
                                                  float_format='%.17g', lineterminator='\n')
    metadata = {'schema_version': '1.0', 'status': 'ready_no_response',
        'created_utc': now_utc(), 'first_lock_sha256': first_sha,
        'final_model_metadata_sha256': sha(FINAL_RUN/'final_model_metadata.json'),
        'final_model_manifest_sha256': sha(FINAL_RUN/'model_manifest.jsonl'),
        'final_model_independent_audit_sha256': sha(FINAL_RUN/'final_model_independent_audit.json'),
        'external_manifest_sha256': sha(EXTERNAL_MANIFEST),
        'global_split_manifest_sha256': sha(GLOBAL_SPLIT),
        'expression_genes_sha256': expression_genes_sha,
        'expression_tpm_sha256': expression_tpm_sha,
        'model_artifact_hashes': model_hashes,
        'source_hashes': source_hashes,
        'prediction_manifest_sha256': sha(prediction_path),
        'prediction_rows': len(result), 'identity_audit': identity_audit,
        'prediction_coverage': coverage,
        'resource_usage': resource_snapshot(started),
        'external_response_accessed': False, 'external_response_index_accessed': False}
    (EXTERNAL_RUN/'preflight_metadata.json').write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return metadata


def verify_second_lock_before_response() -> dict:
    """Validate every nonresponse dependency before touching the raw curve file."""
    if not SECOND_LOCK.is_file():
        raise PermissionError('Second P6 execution lock is absent; CTRP response remains sealed')
    saved = _json(SECOND_LOCK)
    if (saved.get('status'), saved.get('external_response_accessed'),
            saved.get('response_source_path'), saved.get('response_column'),
            saved.get('model_count'), saved.get('fit_events'), saved.get('prediction_rows')) != (
            'frozen_pre_external_response', False, RAW_RESPONSE_NAME,
            'area_under_curve', 64, 1744, 9480):
        raise PermissionError('Second P6 execution lock has not passed the hard gate')
    # The builder is intentionally response-free. It independently rechecks
    # first-lock Git/source hashes, all 64 artifacts, the audited final run and
    # the SHA-bound no-response prediction manifest before returning.
    scripts = ROOT/'scripts'
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    from freeze_external_execution_v1 import build_execution_lock
    current = build_execution_lock()
    for key, value in current.items():
        if key != 'created_utc' and saved.get(key) != value:
            raise PermissionError('Second P6 execution lock changed or a bound artifact differs: ' + key)
    preflight_meta = _json(EXTERNAL_RUN/'preflight_metadata.json')
    if (preflight_meta.get('status'), preflight_meta.get('external_response_accessed'),
            preflight_meta.get('prediction_manifest_sha256')) != (
            'ready_no_response', False, sha(EXTERNAL_RUN/'prediction_manifest.csv')):
        raise PermissionError('Response-free preflight state/hash differs')
    if saved['preflight_metadata_sha256'] != sha(EXTERNAL_RUN/'preflight_metadata.json'):
        raise PermissionError('Second lock does not bind the exact preflight metadata')
    for name, expected in saved['source_hashes'].items():
        if sha(ROOT/name) != expected:
            raise PermissionError('Second lock scientific source hash mismatch: ' + name)
    for name, expected in saved['model_artifact_hashes'].items():
        path = (ROOT/name).resolve()
        path.relative_to((FINAL_RUN/'models').resolve())
        if sha(path) != expected:
            raise PermissionError('Second lock final-model artifact hash mismatch: ' + name)
    return saved


def read_frozen_curve_responses(curve_path: Path, external: pd.DataFrame,
                                *, endpoint_field: str = 'area_under_curve') -> tuple[pd.DataFrame, dict]:
    """After the second gate only, read AUC at the 1,201 frozen line pointers.

    Nonselected rows are streamed past without converting a response. One bad
    selected experiment invalidates the *whole* drug/model label; we do not
    quietly take the median of just the remaining favorable replicate.
    """
    if endpoint_field != 'area_under_curve':
        raise ValueError('Only the frozen CTRPv2 post-QC area_under_curve endpoint is allowed')
    expected: dict[tuple[int, str], dict] = {}
    by_line: dict[int, tuple[int, str]] = {}
    for record in external.itertuples(index=False):
        key = (int(record.gdsc_drug_id), record.canonical_model_id)
        if key in expected:
            raise ValueError('Duplicate frozen external drug/model response unit')
        lines = parse_pipe_numbers(record.ctrp_source_row_numbers, field='source row')
        exps = parse_pipe_numbers(record.ctrp_experiment_ids, field='experiment')
        if len(lines) != len(exps) or len(lines) != int(record.n_source_experiments):
            raise ValueError('Source row/experiment count changed')
        expected[key] = {'compound': int(record.ctrp_compound_id),
                         'experiment_ids': set(exps), 'seen_experiments': set(),
                         'values': [], 'invalid': False,
                         'biological_group_id': record.biological_group_id}
        for line in lines:
            if line in by_line:
                raise ValueError('Duplicate frozen source row pointer')
            by_line[line] = key
    seen_lines: set[int] = set()
    with Path(curve_path).open('r', encoding='utf-8-sig', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, None)
        if header is None or len(header) != len(set(header)) or not {'experiment_id', 'master_cpd_id', endpoint_field}.issubset(header):
            raise ValueError('Unexpected or duplicated CTRPv2 post-QC endpoint header')
        positions = {name: header.index(name) for name in ('experiment_id', 'master_cpd_id', endpoint_field)}
        for line, fields in enumerate(reader, 2):
            if line not in by_line:
                continue
            if len(fields) != len(header):
                raise ValueError('Selected CTRPv2 source row has malformed column count')
            seen_lines.add(line)
            key = by_line[line]
            record = expected[key]
            try:
                experiment = int(fields[positions['experiment_id']])
                compound = int(fields[positions['master_cpd_id']])
            except ValueError as exc:
                raise ValueError('Selected source row has invalid experiment/compound ID') from exc
            if (experiment not in record['experiment_ids'] or
                    experiment in record['seen_experiments'] or compound != record['compound']):
                raise ValueError('Selected source row experiment/compound mismatches frozen pointer')
            record['seen_experiments'].add(experiment)
            try:
                response = float(fields[positions[endpoint_field]])
            except ValueError:
                response = float('nan')
            if not np.isfinite(response):
                record['invalid'] = True
            else:
                record['values'].append(response)
            if len(seen_lines) == len(by_line):
                break
    if seen_lines != set(by_line):
        raise ValueError('Frozen CTRPv2 source row pointer is absent from source file')
    rows = []
    invalid = []
    for (drug, model), item in expected.items():
        if item['seen_experiments'] != item['experiment_ids']:
            raise ValueError('Frozen source experiment set differs from selected rows')
        if item['invalid'] or len(item['values']) != len(item['experiment_ids']):
            invalid.append({'drug_id': drug, 'canonical_model_id': model,
                            'reason': 'nonfinite_or_unparseable_selected_source_auc'})
            continue
        rows.append({'drug_id': drug, 'canonical_model_id': model,
                     'biological_group_id': item['biological_group_id'],
                     'observed_auc': float(np.median(item['values']))})
    return pd.DataFrame(rows, columns=['drug_id', 'canonical_model_id',
                                       'biological_group_id', 'observed_auc']), {
        'selected_source_rows': len(by_line), 'response_valid_drug_model_rows': len(rows),
        'invalid_response_drug_model_rows': len(invalid),
        'invalid_response_reason_counts': dict(Counter(item['reason'] for item in invalid)),
        'response_aggregation': 'median_all_eligible_experiments_per_compound_and_canonical_model',
        'source_endpoint_field': endpoint_field,
    }


def _arm_descriptive(responses: pd.DataFrame, predictions: pd.DataFrame,
                     *, drug: int, arm: str, min_groups: int,
                     frozen_ids: tuple[str, ...]) -> dict:
    base = responses.loc[responses.drug_id.eq(drug)].set_index('canonical_model_id')
    rows = predictions.loc[predictions.drug_id.eq(drug) & predictions.model.eq(arm)].set_index('canonical_model_id')
    if rows.index.duplicated().any() or set(rows.index) != set(frozen_ids):
        raise ValueError('Response-free arm lacks exactly the frozen eligible IDs')
    if base.index.duplicated().any() or base.biological_group_id.duplicated().any() or not set(base.index) <= set(frozen_ids):
        raise ValueError('External labels duplicate groups or exceed the eligible roster')
    if any(base.loc[model, 'biological_group_id'] != rows.loc[model, 'biological_group_id'] for model in base.index):
        raise ValueError('External prediction/response biological group mismatch')
    common = base.join(rows[['predicted_auc', 'status']], how='inner')
    predicted = pd.to_numeric(common.predicted_auc, errors='coerce').to_numpy(dtype=float)
    valid = common.status.eq('predicted').to_numpy() & np.isfinite(predicted)
    subset = common.iloc[np.flatnonzero(valid)]
    rho, reason = spearman_checked(subset.observed_auc, subset.predicted_auc,
                                   subset.biological_group_id, min_groups=min_groups)
    return {'drug_id': drug, 'model': arm, 'frozen_eligible_ids': len(frozen_ids),
            'response_valid_ids': len(base), 'predicted_on_response_valid_ids': len(subset),
            'rho': rho, 'undefined_reason': reason,
            'unknown_training_tissue_eligible': int(rows.unknown_training_tissue.sum())
                if 'unknown_training_tissue' in rows else None,
            'zero_molecular_solution_eligible': int(rows.zero_molecular_solution.astype(str).str.lower().eq('true').sum())
                if 'zero_molecular_solution' in rows else None}


def _paired_even_if_no_valid_response(responses: pd.DataFrame, prediction: pd.DataFrame,
                                      *, drug: int, candidate: str, reference: str,
                                      frozen_ids: tuple[str, ...], design_min: int) -> tuple[pd.DataFrame, dict]:
    """Retain a zero-label drug as an explicit undefined fixed-panel member."""
    if responses.drug_id.eq(drug).any():
        return paired_common_table(responses, prediction, drug_id=drug,
            arm_a=candidate, arm_b=reference, design_min_groups=design_min,
            frozen_eligible_ids=frozen_ids)
    by_arm = {}
    for arm in (candidate, reference):
        sub = prediction.loc[prediction.drug_id.eq(drug) & prediction.model.eq(arm)]
        if (sub.canonical_model_id.duplicated().any() or
                set(sub.canonical_model_id) != set(frozen_ids) or
                not sub.status.isin(('predicted', 'abstained')).all()):
            raise ValueError('Response-free zero-label arm is incomplete or malformed')
        by_arm[arm] = sub
    empty = pd.DataFrame(columns=['canonical_model_id', 'biological_group_id',
                                  'observed_auc', 'prediction_a', 'prediction_b'])
    report = {'drug_id': drug, 'arm_a': candidate, 'arm_b': reference,
              'n_frozen_eligible_ids': len(frozen_ids),
              'frozen_eligible_roster_supplied': True,
              'n_response_valid_ids': 0, 'n_missing_or_invalid_response': len(frozen_ids),
              'n_response_valid_groups': 0,
              'n_arm_a_predicted_all_frozen': int(by_arm[candidate].status.eq('predicted').sum()),
              'n_arm_b_predicted_all_frozen': int(by_arm[reference].status.eq('predicted').sum()),
              'n_arm_a_predicted': 0, 'n_arm_b_predicted': 0,
              'n_common_ids': 0, 'n_common_groups': 0,
              'below_design_coverage': True,
              'n_arm_a_abstain_or_missing': 0, 'n_arm_b_abstain_or_missing': 0,
              'arm_a_noncoverage_reasons': {}, 'arm_b_noncoverage_reasons': {}}
    return empty, report


def score_locked_panel(config: dict, external: pd.DataFrame,
                       prediction: pd.DataFrame, responses: pd.DataFrame,
                       *, expected_global_groups: int = 175) -> dict:
    """Pure-data scoring API; tests use synthetic inputs and reduced B only."""
    drugs = tuple(map(int, config['drug_ids']))
    arms = tuple(config['models'])
    if len(drugs) != 8 or len(arms) != 8 or prediction.duplicated(['drug_id', 'model', 'canonical_model_id']).any():
        raise ValueError('External scoring panel/row uniqueness changed')
    global_groups = tuple(sorted(set(external.biological_group_id)))
    if len(global_groups) != expected_global_groups:
        raise ValueError('Frozen external global biological-group universe changed')
    stats = config['statistics']
    boot_config = stats['bootstrap']
    shared = make_shared_bootstrap(global_groups, draws=int(boot_config['draws']),
                                   seed=int(boot_config['seed']),
                                   expected_groups=expected_global_groups)
    min_groups = int(stats['minimum_independent_common_groups_for_spearman'])
    design_min = 80
    arm_rows = []
    frozen_ids = {}
    for drug in drugs:
        frozen_ids[drug] = tuple(external.loc[external.gdsc_drug_id.astype(int).eq(drug),
                                              'canonical_model_id'])
        for arm in arms:
            arm_rows.append(_arm_descriptive(responses, prediction, drug=drug,
                                             arm=arm, min_groups=min_groups,
                                             frozen_ids=frozen_ids[drug]))
    contrasts = [('primary', stats['primary_comparison']['candidate'],
                  stats['primary_comparison']['reference'],
                  int(stats['primary_comparison']['family_size']))]
    contrasts += [('representation_' + item['candidate'].removesuffix('_pair_full'),
                   item['candidate'], item['reference'],
                   int(stats['same_candidate_representation_family']['family_size']))
                  for item in stats['same_candidate_representation_family']['contrasts']]
    comparison_results = {}
    for name, candidate, reference, family_size in contrasts:
        pairs = {}
        coverage = {}
        for drug in drugs:
            pair, report = _paired_even_if_no_valid_response(responses, prediction,
                drug=drug, candidate=candidate, reference=reference,
                design_min=design_min, frozen_ids=frozen_ids[drug])
            pairs[drug], coverage[str(drug)] = pair, report
        inference = bootstrap_paired_rho(pairs, fixed_drug_ids=drugs,
            shared=shared, min_groups=min_groups,
            min_valid_fraction=float(boot_config['minimum_valid_fraction']),
            family_size=family_size)
        comparison_results[name] = {'candidate': candidate, 'reference': reference,
                                    'family_size': family_size, 'coverage': coverage,
                                    'inference': inference}
    descriptive_definitions = [(arm, 'tissue_null') for arm in arms if arm != 'tissue_null']
    descriptive_definitions += [
        ('mad_rank', 'mad_expression'), ('mad_pair_full', 'mad_rank'),
        ('response_association_rank', 'response_association_expression'),
        ('response_association_pair_full', 'response_association_rank')]
    descriptive_pairs = []
    for candidate, reference in descriptive_definitions:
        for drug in drugs:
            pair, coverage = _paired_even_if_no_valid_response(responses, prediction,
                drug=drug, candidate=candidate, reference=reference,
                design_min=design_min, frozen_ids=frozen_ids[drug])
            descriptive_pairs.append({
                'drug_id': drug, 'candidate': candidate, 'reference': reference,
                'descriptive_only': True, 'coverage': coverage,
                **paired_rho_difference(pair, min_groups=min_groups)})
    return {'schema_version': '1.0', 'scope': 'fixed_model_retrospective_cross_screen_rank_transfer',
            'drug_ids': drugs, 'model_arms': arms,
            'per_arm_descriptive': arm_rows, 'contrasts': comparison_results,
            'paired_descriptive': descriptive_pairs,
            'no_cross_screen_rmse_mae_or_r2': True, 'no_external_score_based_selection': True,
            'no_external_direction_flip': True}


def _clean_for_json(value):
    if isinstance(value, dict):
        return {str(k): _clean_for_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clean_for_json(v) for v in value]
    if isinstance(value, np.generic):
        return _clean_for_json(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _same_process_alive(started: dict) -> bool:
    """Reject concurrent execution; a recycled PID is not the original process."""
    if started.get('host') != platform.node():
        raise PermissionError('Cannot establish liveness of an attempt on another host')
    pid = started.get('pid')
    created = started.get('process_create_time')
    if not isinstance(pid, int) or not isinstance(created, (int, float)):
        raise PermissionError('Started marker lacks process identity for safe recovery')
    try:
        observed = psutil.Process(pid).create_time()
    except psutil.NoSuchProcess:
        return False
    except psutil.AccessDenied as exc:
        raise PermissionError('Cannot establish liveness of the previous process') from exc
    return abs(float(observed) - float(created)) <= 1e-6


def _next_execution_attempt(second_sha: str, *, resume: bool = False,
                            resume_unclean: bool = False) -> tuple[Path, Path]:
    """Append-only same-lock ledger; unclean replay needs dead-PID evidence."""
    if resume and resume_unclean:
        raise ValueError('Choose one explicit recovery mode')
    folder = EXTERNAL_RUN/'execution_attempts'
    summary = EXTERNAL_RUN/'evaluation_summary.json'
    if summary.exists():
        raise FileExistsError('External evaluation summary already exists; never evaluate twice')
    started = sorted(folder.glob('attempt_*_started.json')) if folder.exists() else []
    if started and not (resume or resume_unclean):
        raise FileExistsError('External response attempt exists; explicit same-lock recovery is required')
    if (resume or resume_unclean) and not started:
        raise ValueError('No failed external response attempt exists to resume')
    recognized = set()
    for index, marker in enumerate(started, 1):
        expected = folder/f'attempt_{index:04d}_started.json'
        original = _json(marker)
        if marker != expected or original.get('second_lock_sha256') != second_sha:
            raise PermissionError('Execution history is incomplete or belongs to another lock')
        failed = folder/f'attempt_{index:04d}_failed.json'
        completed = folder/f'attempt_{index:04d}_completed.json'
        unclean = folder/f'attempt_{index:04d}_unclean_abandoned.json'
        if completed.exists() or (failed.exists() and unclean.exists()):
            raise PermissionError('A previous attempt is successful or has contradictory terminal records')
        if failed.exists() or unclean.exists():
            terminal = failed if failed.exists() else unclean
            if _json(terminal).get('second_lock_sha256') != second_sha:
                raise PermissionError('Prior attempt terminal belongs to another lock')
            if index == len(started) and resume_unclean:
                raise PermissionError('Last attempt has a terminal record; use --resume')
            recognized.update((marker, terminal))
            continue
        if index != len(started) or not resume_unclean:
            raise PermissionError('A previous attempt is running or unaudited; cannot resume')
        if _same_process_alive(original):
            raise PermissionError('Previous external evaluation process is still alive')
        _write_json_exclusive(unclean, {
            'status': 'unclean_abandoned', 'second_lock_sha256': second_sha,
            'original_started_sha256': sha(marker),
            'original_pid': original['pid'],
            'original_process_create_time': original['process_create_time'],
            'reason': 'previous_process_absent_without_terminal_record',
            'recorded_utc': now_utc()})
        recognized.update((marker, unclean))
    if folder.exists() and set(folder.iterdir()) != recognized:
        raise PermissionError('Execution attempt directory contains an unexpected or orphaned file')
    folder.mkdir(parents=True, exist_ok=True)
    number = len(started) + 1
    return (folder/f'attempt_{number:04d}_started.json',
            folder/f'attempt_{number:04d}')


def _write_json_exclusive(path: Path, payload: dict) -> None:
    with path.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write('\n')


def execute(*, resume: bool = False, resume_unclean: bool = False) -> dict:
    started = time.perf_counter()
    second = verify_second_lock_before_response()
    second_sha = sha(SECOND_LOCK)
    start_path, terminal_stem = _next_execution_attempt(
        second_sha, resume=resume, resume_unclean=resume_unclean)
    process = psutil.Process(os.getpid())
    _write_json_exclusive(start_path, {
        'status': 'started_under_second_lock', 'started_utc': now_utc(),
        'second_lock_sha256': second_sha,
        'resumed_after_failed_same_lock_attempt': resume or resume_unclean,
        'pid': process.pid, 'process_create_time': process.create_time(),
        'host': platform.node()})
    try:
        raw_path = (ROOT/second['response_source_path']).resolve()
        raw_path.relative_to((ROOT/'data_raw/ctrp').resolve())
        if raw_path != (ROOT/RAW_RESPONSE_NAME).resolve():
            raise ValueError('Second lock permits only the frozen CTRPv2 post-QC curve file')
        actual_source_sha = sha(raw_path)  # first response-file access, after both locks
        if actual_source_sha != second['response_source_sha256']:
            raise ValueError('Raw CTRPv2 response source SHA mismatch')
        external = _load_manifest(EXTERNAL_MANIFEST)
        predictions = pd.read_csv(EXTERNAL_RUN/'prediction_manifest.csv',
            dtype={'canonical_model_id': str, 'biological_group_id': str,
                   'model': str, 'status': str}, keep_default_na=False)
        if sha(EXTERNAL_RUN/'prediction_manifest.csv') != second['prediction_manifest_sha256']:
            raise ValueError('Response-free prediction manifest changed after second lock')
        responses, response_audit = read_frozen_curve_responses(raw_path, external,
            endpoint_field=second['response_column'])
        config = _json(EXTERNAL_CONFIG)
        results = score_locked_panel(config, external, predictions, responses)
        preflight_meta = _json(EXTERNAL_RUN/'preflight_metadata.json')
        output = {'schema_version': '1.0', 'status': 'completed_one_time_external_evaluation',
                  'second_lock_sha256': second_sha,
                  'preflight_metadata_sha256': second['preflight_metadata_sha256'],
                  'prediction_manifest_sha256': second['prediction_manifest_sha256'],
                  'source_curve_sha256': actual_source_sha,
                  'response_audit': response_audit,
                  'identity_audit': preflight_meta['identity_audit'],
                  'prediction_coverage': preflight_meta['prediction_coverage'],
                  'evaluation': results,
                  'resource_usage': resource_snapshot(started),
                  'response_values_written_to_repository': False,
                  'individual_response_values_in_output': False}
        summary_path = EXTERNAL_RUN/'evaluation_summary.json'
        with summary_path.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(_clean_for_json(output), handle, ensure_ascii=False, indent=2,
                      allow_nan=False)
            handle.write('\n')
        _write_json_exclusive(Path(str(terminal_stem) + '_completed.json'), {
            'status': 'completed', 'started_under_second_lock': second['status'],
            'second_lock_sha256': second_sha,
            'evaluation_summary_sha256': sha(summary_path),
            'ended_utc': now_utc()})
        return output
    except BaseException:
        _write_json_exclusive(Path(str(terminal_stem) + '_failed.json'), {
            'status': 'failed_preserved', 'second_lock_sha256': second_sha,
            'error': traceback.format_exc(), 'ended_utc': now_utc()})
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--preflight', action='store_true', help='No CTRP response/index access')
    modes.add_argument('--execute', action='store_true', help='Only after the second execution lock')
    recovery = parser.add_mutually_exclusive_group()
    recovery.add_argument('--resume', action='store_true',
                          help='Same-lock replay after a preserved failed execution')
    recovery.add_argument('--resume-unclean', action='store_true',
                          help='Same-lock replay only after a dead prior process and no terminal record')
    args = parser.parse_args()
    if (args.resume or args.resume_unclean) and not args.execute:
        parser.error('Recovery flags require --execute')
    if args.preflight:
        outcome = preflight()
        print(json.dumps({'status': outcome['status'], 'prediction_rows': outcome['prediction_rows'],
                          'prediction_manifest_sha256': outcome['prediction_manifest_sha256'],
                          'external_response_accessed': False}, ensure_ascii=False))
    else:
        outcome = execute(resume=args.resume, resume_unclean=args.resume_unclean)
        print(json.dumps({'status': outcome['status'],
                          'second_lock_sha256': outcome['second_lock_sha256'],
                          'source_curve_sha256': outcome['source_curve_sha256']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
