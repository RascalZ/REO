"""Prediction-only replay of finite numeric portable REO model state.

This module never trains, scores, loads pickle/joblib, or resolves external
cohort eligibility. Identity and tissue labels must already be canonicalized.
Only NumPy/SciPy numeric operations and explicit JSON inputs are supported.
"""
from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import rankdata

MODEL_SCHEMA = 'reo-portable-model-1.0'
INPUT_SCHEMA = 'reo-inference-input-1.0'
OUTPUT_SCHEMA = 'reo-inference-output-1.0'
INPUT_KEYS = {'schema_version', 'sample_ids', 'biological_group_ids',
              'tissues', 'gene_ids', 'raw_tpm'}
REQUIRED_BUNDLE_KEYS = {
    'schema_version', 'task', 'source', 'input_gene_ids', 'tissue_code',
    'model_kind', 'fitted_kind', 'selection_rule', 'candidate_gene_indices',
    'candidate_gene_ids', 'feature_names', 'pair_indices', 'feature_mean',
    'feature_scale', 'feature_tissue_means', 'target_tissue_means',
    'training_tissue_codes', 'target_mean', 'target_sd', 'coefficients',
    'n_nonzero_molecular_coeffs', 'declared_molecular_arm', 'fallback_reason',
    'solver_reason', 'chosen_params',
}


def checked_response_free_path(path: str | os.PathLike[str]) -> Path:
    """Reject known sealed locations before opening any user-specified file.

    This is a local convenience guard, not an OS sandbox or a content scanner.
    It also rejects external run outputs and all old response-table basenames.
    Resolving first catches symlink/junction aliases on the current OS.
    """
    result = Path(path).expanduser().resolve()
    normalized = str(result).replace('\\', '/').casefold()
    parts = normalized.split('/')
    for index, part in enumerate(parts):
        remaining = parts[index+1:]
        if part == 'data_raw' and remaining and remaining[0] == 'ctrp':
            raise PermissionError('Sealed CTRP raw paths are forbidden in response-free inference')
        if part == 'data_processed' and remaining and remaining[0].startswith('ctrp'):
            raise PermissionError('Sealed CTRP processed paths are forbidden in response-free inference')
        if part == 'runs' and any('external' in segment for segment in remaining):
            raise PermissionError('External run paths are forbidden in response-free inference')
    basename = parts[-1]
    if basename == 'ctrp--drug--cellline--auc.txt' or (
        'ctrp' in basename and any(token in basename for token in ('response', 'curve', 'auc', 'index'))
    ):
        raise PermissionError('Legacy CTRP response/index table paths are forbidden')
    return result


def _reject_constant(value: str) -> None:
    raise ValueError('Nonstandard JSON numeric constant is forbidden: '+value)


def _unique_json_pairs(items: list[tuple[str, Any]]) -> dict:
    result: dict = {}
    for key, value in items:
        if key in result:
            raise ValueError('Duplicate JSON object key: '+key)
        result[key] = value
    return result


def load_json_response_free(path: str | os.PathLike[str]) -> dict:
    """Read only the explicit JSON/JSON.gz path after the sealed-path check."""
    target = checked_response_free_path(path)
    if not (target.name.casefold().endswith('.json') or target.name.casefold().endswith('.json.gz')):
        raise ValueError('Only JSON or JSON.gz inputs are accepted')
    opener = gzip.open if target.name.casefold().endswith('.gz') else open
    with opener(target, 'rt', encoding='utf-8-sig') as handle:
        result = json.load(handle, parse_constant=_reject_constant, object_pairs_hook=_unique_json_pairs)
    if not isinstance(result, dict):
        raise ValueError('JSON root must be an object')
    return result


def _clean_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value) and value.strip() == value


def _list(value: Any, name: str) -> list:
    if not isinstance(value, list):
        raise ValueError(name+' must be a JSON list')
    return value


def _string_list(value: Any, name: str, *, allow_empty: bool = False) -> list[str]:
    values = _list(value, name)
    if (not allow_empty and not values) or any(not _clean_string(item) for item in values):
        raise ValueError(name+' must contain nonempty canonical strings')
    if len(set(values)) != len(values):
        raise ValueError(name+' contains duplicates')
    return values


def _finite_array(value: Any, name: str, shape: tuple[int, ...]) -> np.ndarray:
    if value is None and 0 in shape:
        return np.empty(shape, dtype=np.float64)
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(name+' must be a finite numeric array') from exc
    if array.size == 0 and 0 in shape:
        array = array.reshape(shape)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(name+' has invalid shape or nonfinite values')
    return array


def _finite_scalar(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name+' must be a finite numeric scalar')
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(name+' must be finite')
    return result


def validate_bundle(bundle: dict) -> dict:
    """Validate portable dimensions and return ephemeral numeric state."""
    if not isinstance(bundle, dict) or not REQUIRED_BUNDLE_KEYS.issubset(bundle):
        raise ValueError('Portable model is missing required fields')
    if bundle['schema_version'] != MODEL_SCHEMA:
        raise ValueError('Unsupported portable model schema')
    if not isinstance(bundle['task'], dict) or not isinstance(bundle['source'], dict):
        raise ValueError('Model task and source must be objects')
    genes = _string_list(bundle['input_gene_ids'], 'input_gene_ids')
    mapping = bundle['tissue_code']
    if not isinstance(mapping, dict) or any(not _clean_string(k) or isinstance(v, bool) or not isinstance(v, int) for k, v in mapping.items()):
        raise ValueError('tissue_code must map canonical labels to integer codes')
    if len(set(mapping.values())) != len(mapping):
        raise ValueError('tissue_code has ambiguous duplicate codes')
    training = _list(bundle['training_tissue_codes'], 'training_tissue_codes')
    if not training or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(v) for v in training):
        raise ValueError('Training tissue codes must be finite numeric values')
    if len(set(training)) != len(training) or any(v not in mapping.values() for v in training):
        raise ValueError('Training tissue codes are duplicated or absent from vocabulary')
    if training != sorted(training):
        raise ValueError('Training tissue codes must retain sorted encoder category order')
    kind = bundle['fitted_kind']
    if kind not in {'expression', 'rank', 'pair_full', 'tissue_null', 'mean'}:
        raise ValueError('Unsupported fitted prediction branch')
    declared = bundle['declared_molecular_arm']
    if not isinstance(declared, bool):
        raise ValueError('declared_molecular_arm must be boolean')
    if not declared and (kind != 'tissue_null' or bundle['model_kind'] != 'tissue_null'):
        raise ValueError('Only the explicit tissue-null arm can bypass molecular validity')
    if kind in {'expression', 'rank', 'pair_full'} and bundle['model_kind'] != kind:
        raise ValueError('Molecular model_kind/fitted_kind mismatch')
    if not isinstance(bundle['chosen_params'], dict):
        raise ValueError('chosen_params must be an object')
    for name in ('selection_rule', 'fallback_reason', 'solver_reason'):
        if not isinstance(bundle[name], str):
            raise ValueError(name+' must be a string')
    indices = _list(bundle['candidate_gene_indices'], 'candidate_gene_indices') if bundle['candidate_gene_indices'] is not None else []
    if any(isinstance(i, bool) or not isinstance(i, int) or i < 0 or i >= len(genes) for i in indices) or len(set(indices)) != len(indices):
        raise ValueError('Candidate indices must be unique integer input positions')
    candidate_ids = [] if bundle['candidate_gene_ids'] is None else _string_list(bundle['candidate_gene_ids'], 'candidate_gene_ids', allow_empty=True)
    if candidate_ids != [genes[i] for i in indices] or candidate_ids != sorted(candidate_ids):
        raise ValueError('Candidate IDs must match indices in stable-ID output order')
    names = [] if bundle['feature_names'] is None else _string_list(bundle['feature_names'], 'feature_names', allow_empty=True)
    p = len(names)
    pairs = [] if bundle['pair_indices'] is None else _list(bundle['pair_indices'], 'pair_indices')
    if kind in {'tissue_null', 'mean'}:
        if indices or names or pairs:
            raise ValueError('Nonmolecular fitted branches must not invent encoded features')
    elif kind in {'expression', 'rank'}:
        prefix = 'log2TPM:' if kind == 'expression' else 'sample_rank:'
        if not indices or names != [prefix+gene for gene in candidate_ids] or pairs:
            raise ValueError('Expression/rank feature order is inconsistent')
    else:
        expected_pairs = [[indices[a], indices[b]] for a in range(len(indices)) for b in range(a+1, len(indices))]
        if not expected_pairs or pairs != expected_pairs:
            raise ValueError('Pair features must be the complete, oriented candidate upper triangle')
        if names != ['REO:'+genes[a]+'>'+genes[b] for a, b in pairs]:
            raise ValueError('Pair names do not match full-input pair positions')
    mean = _finite_array(bundle['feature_mean'], 'feature_mean', (p,))
    scale = _finite_array(bundle['feature_scale'], 'feature_scale', (p,))
    if np.any(scale <= 0):
        raise ValueError('All fitted feature scales must be positive')
    projection = _finite_array(bundle['feature_tissue_means'], 'feature_tissue_means', (len(training), p))
    tissue_mean = _finite_array(bundle['target_tissue_means'], 'target_tissue_means', (len(training),))
    coef = _finite_array(bundle['coefficients'], 'coefficients', (p,))
    nonzero = bundle['n_nonzero_molecular_coeffs']
    if isinstance(nonzero, bool) or not isinstance(nonzero, int) or nonzero != int(np.count_nonzero(coef)):
        raise ValueError('Exact nonzero coefficient count mismatch')
    target_mean = _finite_scalar(bundle['target_mean'], 'target_mean')
    target_sd = _finite_scalar(bundle['target_sd'], 'target_sd')
    if target_sd <= 0:
        raise ValueError('Fitted target SD must be positive')
    return dict(genes=genes, mapping=mapping, training=training, kind=kind, declared=declared,
                indices=indices, pairs=pairs, mean=mean, scale=scale, projection=projection,
                tissue_mean=tissue_mean, coef=coef, target_mean=target_mean, target_sd=target_sd)


def load_bundle(path: str | os.PathLike[str]) -> dict:
    """Load and validate a portable model, never a serialized estimator."""
    bundle = load_json_response_free(path)
    validate_bundle(bundle)
    return bundle


def predict_bundle(bundle: dict, sample: dict) -> dict:
    """Replay a fixed numeric model without fitting or reading response labels.

    Hard failures are structural/schema/duplicate-identity faults. Row-level
    identity, tissue and TPM failures are retained as explicit abstentions.
    This does not certify external eligibility, compound mapping or genealogy.
    """
    state = validate_bundle(bundle)
    if not isinstance(sample, dict) or set(sample) != INPUT_KEYS:
        raise ValueError('Inference input must have exactly the response-free input keys; labels are forbidden')
    if sample['schema_version'] != INPUT_SCHEMA:
        raise ValueError('Unsupported inference input schema')
    if sample['gene_ids'] != state['genes']:
        raise ValueError('Input gene IDs must match the complete ordered model schema')
    row_fields = ('sample_ids', 'biological_group_ids', 'tissues', 'raw_tpm')
    values = {name: _list(sample[name], name) for name in row_fields}
    n = len(values['sample_ids'])
    if n == 0 or any(len(v) != n for v in values.values()):
        raise ValueError('Aligned nonempty sample rows are required')
    for field in ('sample_ids', 'biological_group_ids'):
        valid_ids = [value for value in values[field] if _clean_string(value)]
        if len(set(valid_ids)) != len(valid_ids):
            raise ValueError('Duplicate '+field+' is a hard identity failure')
    gene_count = len(state['genes'])
    for row in values['raw_tpm']:
        if row is not None and (not isinstance(row, list) or len(row) != gene_count):
            raise ValueError('Every TPM row must have the complete gene count or be null')
        if row is not None and any(isinstance(value, (list, dict)) for value in row):
            raise ValueError('TPM entries must be scalar values, not nested dimensions')
    result = dict(schema_version=OUTPUT_SCHEMA, task=bundle['task'],
                  sample_ids=list(values['sample_ids']), predictions=[None]*n,
                  statuses=['abstained']*n, reasons=['']*n, unknown_tissue=[None]*n,
                  zero_molecular_solution=(bundle['n_nonzero_molecular_coeffs'] == 0) if state['declared'] else None,
                  no_response_labels_used=True)
    valid_rows, tpm_rows, tissue_rows = [], [], []
    for index in range(n):
        if not _clean_string(values['sample_ids'][index]) or not _clean_string(values['biological_group_ids'][index]):
            result['reasons'][index] = 'missing_or_invalid_identity'
            continue
        tissue = values['tissues'][index]
        if not _clean_string(tissue):
            result['reasons'][index] = 'missing_or_invalid_tissue'
            continue
        code = state['mapping'].get(tissue)
        unknown = code not in state['training']
        result['unknown_tissue'][index] = unknown
        d = np.asarray([float(code == item) for item in state['training']], dtype=np.float64)
        if state['declared']:
            raw = values['raw_tpm'][index]
            if raw is None or any(value is None for value in raw):
                result['reasons'][index] = 'missing_molecular_input'
                continue
            if any(isinstance(value, bool) for value in raw):
                result['reasons'][index] = 'nonnumeric_molecular_input'
                continue
            try:
                row = np.asarray(raw, dtype=np.float64)
            except (ValueError, TypeError, OverflowError):
                result['reasons'][index] = 'nonnumeric_molecular_input'
                continue
            if not np.isfinite(row).all():
                result['reasons'][index] = 'nonfinite_molecular_input'
                continue
            if np.any(row < 0):
                result['reasons'][index] = 'negative_molecular_input'
                continue
            tpm_rows.append(row)
        valid_rows.append(index)
        tissue_rows.append(d)
    if valid_rows:
        d = np.asarray(tissue_rows, dtype=np.float64)
        if state['kind'] == 'mean':
            predictions = np.full(len(valid_rows), state['target_mean'], dtype=np.float64)
        elif state['kind'] == 'tissue_null':
            predictions = state['target_mean']+state['target_sd']*(d@state['tissue_mean'])
        else:
            raw = np.asarray(tpm_rows, dtype=np.float64)
            if state['kind'] == 'expression':
                encoded = np.log2(raw[:, state['indices']]+1.0)
            elif state['kind'] == 'rank':
                encoded = rankdata(raw[:, state['indices']], axis=1, method='average').astype(np.float64)
            else:
                pairs = np.asarray(state['pairs'], dtype=np.int64)
                left, right = raw[:, pairs[:, 0]], raw[:, pairs[:, 1]]
                encoded = np.greater(left, right).astype(np.int8)-np.less(left, right).astype(np.int8)
            z = (encoded-state['mean'])/state['scale']
            standardized = (z-d@state['projection'])@state['coef']+d@state['tissue_mean']
            predictions = state['target_mean']+state['target_sd']*standardized
        if not np.isfinite(predictions).all():
            raise ValueError('Nonfinite fixed-model prediction is a hard failure')
        for index, prediction in zip(valid_rows, predictions):
            result['predictions'][index] = float(prediction)
            result['statuses'][index] = 'predicted'
    return result
