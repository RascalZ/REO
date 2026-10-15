"""One frozen full-development fit per drug and model arm.

This module has no file access.  It reuses the diagnostic estimator and the
repeated experiment's arm definitions; the only new operation is selecting a
single grid point on three grouped, tissue-stratified development folds before
refitting on every available development observation for that drug.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error

from modeling.nested_cv import make_grouped_splits, candidate_grid, _grid_boundary_details
from ordinal_diagnostic_models import DiagnosticRegressor
from ordinal_repeated_models import model_spec


def make_global_inner_manifest(identities, *, n_folds, seed):
    """Split the full 626-group identity frame once, before drug intersection."""
    required = {'canonical_model_id', 'biological_group_id', 'tissue'}
    if not required.issubset(identities.columns):
        raise ValueError('Global identity columns missing')
    ids = identities.loc[:, list(('canonical_model_id', 'biological_group_id', 'tissue'))].copy()
    if ids.isna().any().any() or ids.astype(str).eq('').any().any():
        raise ValueError('Empty global identity field')
    if ids.canonical_model_id.duplicated().any() or ids.biological_group_id.duplicated().any():
        raise ValueError('The frozen development cohort requires one identity per biological group')
    ids = ids.sort_values('canonical_model_id', kind='stable').reset_index(drop=True)
    if not isinstance(n_folds, int) or n_folds < 2 or not isinstance(seed, int):
        raise ValueError('Invalid final-fit split configuration')
    rows = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        splits = make_grouped_splits(ids.biological_group_id, ids.tissue, n_folds, seed)
    for fold, (train, validation) in enumerate(splits):
        for role, positions in [('train', train), ('validation', validation)]:
            for index in positions:
                rows.append({**ids.iloc[int(index)].to_dict(), 'inner_fold': fold, 'role': role})
    manifest = pd.DataFrame(rows, columns=[*ids.columns, 'inner_fold', 'role'])
    validate_global_inner_manifest(ids, manifest, n_folds)
    return manifest, [str(item.message) for item in caught]


def validate_global_inner_manifest(identities, manifest, n_folds):
    expected = identities.set_index('canonical_model_id')
    if len(expected) != expected.index.nunique() or expected.biological_group_id.duplicated().any():
        raise ValueError('Nonunique global identities')
    if set(manifest.inner_fold) != set(range(n_folds)) or set(manifest.role) != {'train', 'validation'}:
        raise ValueError('Incomplete global inner-fold manifest')
    for fold in range(n_folds):
        block = manifest.loc[manifest.inner_fold.eq(fold)]
        train = block.loc[block.role.eq('train')]
        valid = block.loc[block.role.eq('validation')]
        if train.empty or valid.empty or len(block) != len(expected):
            raise ValueError('Empty or incomplete global inner fold')
        if set(train.canonical_model_id) & set(valid.canonical_model_id):
            raise ValueError('Global inner-fold identity overlap')
        if set(train.biological_group_id) & set(valid.biological_group_id):
            raise ValueError('Global inner-fold biological-group overlap')
        if set(block.canonical_model_id) != set(expected.index):
            raise ValueError('Global inner fold does not cover all identities')
        for row in block.itertuples():
            reference = expected.loc[row.canonical_model_id]
            if row.biological_group_id != reference.biological_group_id or row.tissue != reference.tissue:
                raise ValueError('Global inner-fold identity metadata changed')
    counts = manifest.loc[manifest.role.eq('validation')].groupby('canonical_model_id').size()
    if counts.to_dict() != {key: 1 for key in expected.index}:
        raise ValueError('Global inner validation coverage is not exactly once')


def local_inner_splits(drug_data, manifest, n_folds):
    """Intersect the same outcome-free global assignment with one drug."""
    data = drug_data.reset_index(drop=True)
    if data.canonical_model_id.duplicated().any() or data.biological_group_id.duplicated().any():
        raise ValueError('A drug has duplicate model or biological-group identity')
    result = []
    for fold in range(n_folds):
        block = manifest.loc[manifest.inner_fold.eq(fold)]
        train_ids = set(block.loc[block.role.eq('train'), 'canonical_model_id'])
        valid_ids = set(block.loc[block.role.eq('validation'), 'canonical_model_id'])
        it = np.flatnonzero(data.canonical_model_id.isin(train_ids).to_numpy())
        iv = np.flatnonzero(data.canonical_model_id.isin(valid_ids).to_numpy())
        if not len(it) or not len(iv) or len(it) + len(iv) != len(data):
            raise ValueError('Drug inner fold has empty or incomplete intersection')
        if set(data.iloc[it].biological_group_id) & set(data.iloc[iv].biological_group_id):
            raise ValueError('Drug inner-fold biological-group overlap')
        result.append((it, iv))
    coverage = np.zeros(len(data), dtype=np.int8)
    for _, iv in result:
        coverage[iv] += 1
    if not np.all(coverage == 1):
        raise ValueError('Drug inner OOF validation coverage is not exactly once')
    return result


def prediction_sha256(predictions):
    values = np.ascontiguousarray(predictions, dtype='<f8')
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError('Predictions must be finite float64 vector')
    return hashlib.sha256(values.tobytes()).hexdigest()


def fit_full_development_model(X, y, groups, inner_splits, model_name,
                               ordinal_config, stable_gene_ids, audit_sink=None):
    """Select the frozen grid on pooled inner OOF RMSE, then refit all rows.

    The returned model is complete (builder, gene IDs, pair orientation,
    standardizer, FWL tissue state, coefficient vector).  Null solutions and
    fallbacks are returned as ordinary scientific outcomes, never filtered.
    """
    X = np.asarray(X)
    y = np.asarray(y, dtype=np.float64)
    groups = np.asarray(groups)
    if X.dtype != np.dtype(np.float64) or X.ndim != 2 or X.shape != (len(y), len(stable_gene_ids)+2):
        raise ValueError('Full-development X must be aligned raw float64 TPM plus tissue and audit ID')
    if not len(y) or y.ndim != 1 or groups.shape != (len(y),) or not np.isfinite(y).all():
        raise ValueError('Invalid aligned development response or group')
    if len(np.unique(X[:, -1])) != len(y) or len(np.unique(groups)) != len(y):
        raise ValueError('Duplicate development audit ID or biological group')
    if not np.isfinite(X).all() or (X[:, :-2] < 0).any():
        raise ValueError('Invalid development TPM')
    spec = model_spec(model_name, ordinal_config)
    adapted = deepcopy(ordinal_config)
    if model_name != 'tissue_null':
        adapted['feature']['n_genes'] = spec['n_genes']
        adapted['feature']['diagnostic_genes'] = spec['n_genes']
    grid = [{}] if model_name == 'tissue_null' else candidate_grid('expression', ordinal_config)
    if model_name != 'tissue_null' and len(grid) != 10:
        raise ValueError('Frozen molecular grid must have exactly ten entries')
    expected_folds = int(ordinal_config['inner_folds'])
    if len(inner_splits) != expected_folds:
        raise ValueError('Final inner-fold count differs from frozen ordinal model')
    coverage = np.zeros(len(y), dtype=np.int8)
    for it, iv in inner_splits:
        it, iv = np.asarray(it), np.asarray(iv)
        if not len(it) or not len(iv) or len(set(it)) != len(it) or len(set(iv)) != len(iv):
            raise ValueError('Empty or duplicate inner positions')
        if min(it.min(), iv.min()) < 0 or max(it.max(), iv.max()) >= len(y):
            raise ValueError('Inner position outside development data')
        if set(it) & set(iv) or set(it) | set(iv) != set(range(len(y))):
            raise ValueError('Every inner fold must partition full development data')
        if set(groups[it]) & set(groups[iv]):
            raise ValueError('Inner biological-group overlap')
        coverage[iv] += 1
    if not np.all(coverage == 1):
        raise ValueError('Inner validation coverage must be exactly once')

    events, scores = [], []
    best_score, best_params = float('inf'), None
    if model_name != 'tissue_null':
        for candidate, params in enumerate(grid):
            oof = np.full(len(y), np.nan, dtype=np.float64)
            fallbacks = []
            for fold, (it, iv) in enumerate(inner_splits):
                model = DiagnosticRegressor(spec['representation'], spec['selection_rule'], adapted,
                                            stable_gene_ids, **params).fit(X[it], y[it])
                event = {**model.audit(), 'declared_model': model_name,
                         'candidate_gene_budget': spec['n_genes'], 'stage': 'inner_training',
                         'candidate': candidate, 'inner_fold': fold, 'params': params}
                if set(event['fit_sample_ids']) != set(X[it, -1].astype(int)):
                    raise ValueError('Inner fit membership audit mismatch')
                events.append(event)
                if audit_sink is not None:
                    audit_sink(event)
                pred = model.predict(X[iv])
                if not np.isfinite(pred).all():
                    raise ValueError('Nonfinite inner prediction')
                oof[iv] = pred
                if model.fallback_reason_:
                    fallbacks.append(model.fallback_reason_)
            if not np.isfinite(oof).all():
                raise ValueError('Incomplete inner OOF predictions')
            score = float(np.sqrt(mean_squared_error(y, oof)))
            scores.append({'candidate': candidate, 'params': params,
                           'pooled_inner_rmse': score, 'fallbacks': fallbacks,
                           'inner_oof_prediction_sha256': prediction_sha256(oof)})
            if score < best_score:  # exact ties retain the first frozen grid point
                best_score, best_params = score, params.copy()
    else:
        best_score, best_params = None, {}
    model = DiagnosticRegressor(spec['representation'], spec['selection_rule'], adapted,
                                stable_gene_ids, **best_params).fit(X, y)
    event = {**model.audit(), 'declared_model': model_name,
             'candidate_gene_budget': spec['n_genes'], 'stage': 'full_development_refit',
             'candidate': -1, 'inner_fold': -1, 'params': best_params}
    if set(event['fit_sample_ids']) != set(X[:, -1].astype(int)):
        raise ValueError('Full-development fit membership audit mismatch')
    events.append(event)
    if audit_sink is not None:
        audit_sink(event)
    predictions = model.predict(X)
    if not np.isfinite(predictions).all():
        raise ValueError('Nonfinite full-development prediction')
    boundary = _grid_boundary_details(best_params, grid)
    return model, {'declared_model': model_name, 'selection_rule': spec['selection_rule'],
                   'representation': spec['representation'], 'candidate_gene_budget': spec['n_genes'],
                   'chosen_params': best_params, 'best_inner_rmse': best_score,
                   'inner_scores': scores, 'chosen_grid_boundary_details': boundary,
                   'chosen_grid_boundary': bool(any(v['at_boundary'] and v['search_extent_diagnostic']
                                                    for v in boundary.values())),
                   'full_fit_audit': model.audit(), 'feature_audit': model.feature_audit(),
                   'fit_events': events, 'fit_count': len(events),
                   'development_prediction_sha256': prediction_sha256(predictions),
                   'development_response_sha256': hashlib.sha256(np.ascontiguousarray(y, dtype='<f8').tobytes()).hexdigest(),
                   'development_sample_ids': X[:, -1].astype(int).tolist()}
