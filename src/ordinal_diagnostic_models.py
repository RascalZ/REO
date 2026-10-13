"""Versioned development diagnostic: two gene rules and three representations.

This module does not alter the frozen v1 pipeline. All learned quantities use
only the rows passed to the current fit. No file or external-data access occurs.
The complete pair representation retains constant and high-tie columns so its
raw representation is information-equivalent to sample-wise average ranks.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
from numbers import Integral, Real
import time

import numpy as np
from scipy.stats import rankdata
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import OneHotEncoder
from sklearn.utils.validation import check_is_fitted

from modeling.features import FoldFeatureBuilder
from modeling.nested_cv import TrainOnlyRegressor, candidate_grid, metrics, _grid_boundary_details


def _names_hash(names):
    return hashlib.sha256('\0'.join(names).encode('utf-8')).hexdigest()


class DiagnosticFeatureBuilder(TransformerMixin, BaseEstimator):
    """Train-only gene selection followed by one of three frozen encodings.

    ``response_association`` ranks absolute Pearson association after centering
    log2(TPM+1) and continuous response within each training tissue. Detection
    and nonzero within-tissue expression variance define eligibility. If the
    response is exactly constant within every tissue, scores are zero and
    stable IDs break ties; this condition is explicitly audited. No MAD filter
    or second MAD ranking is applied to the response-association rule.

    All selected columns are ordered by stable gene ID, in every arm. This is
    a deliberately versioned column convention, not bitwise replication of the
    old expression arm's decreasing-MAD column order.
    """
    _BLOCK_SIZE = 512

    def __init__(self, kind='expression', selection_rule='mad', n_genes=200,
                 min_tpm=1.0, min_detection_fraction=0.1, stable_gene_ids=None):
        self.kind = kind
        self.selection_rule = selection_rule
        self.n_genes = n_genes
        self.min_tpm = min_tpm
        self.min_detection_fraction = min_detection_fraction
        self.stable_gene_ids = stable_gene_ids

    def _validate(self, X, fitting):
        if self.kind not in {'expression', 'rank', 'pair_full'}:
            raise ValueError('kind must be expression, rank, or pair_full')
        if self.selection_rule not in {'mad', 'response_association'}:
            raise ValueError('Unsupported selection_rule')
        if isinstance(self.n_genes, bool) or not isinstance(self.n_genes, Integral) or self.n_genes < 1:
            raise ValueError('n_genes must be a positive integer')
        for name, value, upper in [('min_tpm', self.min_tpm, None),
                                   ('min_detection_fraction', self.min_detection_fraction, 1.0)]:
            if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or value < 0 or (upper is not None and value > upper):
                raise ValueError('Invalid '+name)
        if not isinstance(self.stable_gene_ids, tuple) or any(not isinstance(g, str) or not g.strip() for g in self.stable_gene_ids) or len(set(self.stable_gene_ids)) != len(self.stable_gene_ids):
            raise ValueError('stable_gene_ids must be a tuple of unique nonempty strings')
        X = np.asarray(X)
        if X.ndim != 2 or X.dtype != np.dtype(np.float64) or min(X.shape) < 1:
            raise ValueError('X must be a nonempty two-dimensional float64 raw TPM matrix')
        expected = len(self.stable_gene_ids) if fitting else self.n_features_in_
        if X.shape[1] != expected:
            raise ValueError('Input columns must match aligned stable gene IDs')
        for start in range(0, X.shape[1], self._BLOCK_SIZE):
            block = X[:, start:start+self._BLOCK_SIZE]
            if not np.isfinite(block).all() or np.any(block < 0):
                raise ValueError('Raw TPM must be finite and nonnegative')
        return X

    def fit(self, X, y=None, tissue=None):
        for name in list(vars(self)):
            if name.endswith('_'):
                delattr(self, name)
        X = self._validate(X, fitting=True)
        n, p = X.shape
        ids = np.asarray(self.stable_gene_ids, dtype=str)
        detection = np.empty(p, dtype=np.float64)
        scores = np.zeros(p, dtype=np.float64)
        eligibility = np.zeros(p, dtype=bool)
        signed = np.zeros(p, dtype=np.float64)
        exact_null_response = None
        if self.selection_rule == 'response_association':
            y = np.asarray(y, dtype=np.float64)
            tissue = np.asarray(tissue)
            if y.shape != (n,) or tissue.shape != (n,) or not np.isfinite(y).all():
                raise ValueError('Aligned finite training response and tissue are required')
            encoder = OneHotEncoder(handle_unknown='ignore', sparse_output=False, dtype=np.float64)
            d = encoder.fit_transform(tissue.reshape(-1, 1))
            tissue_rows = [d[:, j] == 1.0 for j in range(d.shape[1])]
            counts = d.sum(axis=0)
            yr = y-d@((d.T@y)/counts)
            exact_null_response = bool(all(np.ptp(y[rows]) == 0.0 for rows in tissue_rows))
            if exact_null_response:
                yr = np.zeros_like(y)
            y_norm = float(np.linalg.norm(yr))
        for start in range(0, p, self._BLOCK_SIZE):
            stop = min(start+self._BLOCK_SIZE, p)
            block = X[:, start:stop]
            detection[start:stop] = np.mean(block >= self.min_tpm, axis=0)
            logblock = np.log2(block+1.0)
            if self.selection_rule == 'mad':
                mad = np.median(np.abs(logblock-np.median(logblock, axis=0)), axis=0)
                scores[start:stop] = mad
                eligibility[start:stop] = mad > 0.0
            else:
                residual = logblock-d@((d.T@logblock)/counts[:, None])
                # Exact tissue constants have no association; avoid introducing
                # tiny artificial variance through repeated sum/count roundoff.
                exact_constant = np.logical_and.reduce([np.ptp(logblock[rows], axis=0) == 0.0 for rows in tissue_rows])
                residual[:, exact_constant] = 0.0
                norms = np.linalg.norm(residual, axis=0)
                valid = norms > 0.0
                eligibility[start:stop] = valid
                local = np.zeros(stop-start, dtype=np.float64)
                if y_norm > 0.0:
                    local[valid] = np.clip((residual[:, valid].T@yr)/(norms[valid]*y_norm), -1.0, 1.0)
                signed[start:stop] = local
                scores[start:stop] = np.abs(local)
        detected = detection >= self.min_detection_fraction
        eligible = np.flatnonzero(detected & eligibility)
        order = np.lexsort((ids[eligible], -scores[eligible]))
        ranked = eligible[order[:self.n_genes]]
        selected = ranked[np.argsort(ids[ranked], kind='stable')]
        genes = tuple(self.stable_gene_ids[i] for i in selected)
        pairs = np.empty((0, 2), dtype=np.int64)
        if self.kind == 'pair_full':
            left, right = np.triu_indices(len(selected), k=1)
            pairs = np.column_stack((selected[left], selected[right]))
            names = tuple(f'REO:{ids[a]}>{ids[b]}' for a, b in pairs)
        elif self.kind == 'rank':
            names = tuple('sample_rank:'+g for g in genes)
        else:
            names = tuple('log2TPM:'+g for g in genes)
        self.n_features_in_ = p
        self.gene_indices_ = selected.astype(np.int64, copy=True)
        self.gene_ids_ = genes
        self.selection_ranked_gene_ids_ = tuple(self.stable_gene_ids[i] for i in ranked)
        self.selected_scores_ = scores[selected].copy()
        self.selected_signed_associations_ = signed[selected].copy()
        self.pair_indices_ = pairs
        self.feature_names_ = np.asarray(names, dtype=str)
        self.n_output_features_ = len(names)
        self.audit_metadata_ = {
            'kind': self.kind, 'selection_rule': self.selection_rule,
            'training_samples': int(n), 'input_genes': int(p),
            'requested_genes': int(self.n_genes), 'selected_genes': len(genes),
            'gene_shortfall': max(0, int(self.n_genes)-len(genes)),
            'detection_eligible_genes': int(np.count_nonzero(detected)),
            'eligible_genes': int(len(eligible)), 'output_features': len(names),
            'min_tpm': float(self.min_tpm), 'min_detection_fraction': float(self.min_detection_fraction),
            'response_used_for_feature_construction': self.selection_rule == 'response_association',
            'selection_score': 'training log2(TPM+1) MAD' if self.selection_rule == 'mad' else 'absolute training within-tissue Pearson(log2(TPM+1), response)',
            'selection_order': 'score descending; stable gene ID ascending for exact score ties',
            'output_gene_order': 'stable gene ID ascending in all three arms',
            'selected_gene_ids_sha256': _names_hash(genes),
            'built_feature_names_sha256': _names_hash(names),
            'exact_tissue_constant_selection_response': exact_null_response,
            'zero_response_selection_policy': 'eligible genes retained with zero scores and stable-ID ties',
            'rank_rule': 'per sample, only selected genes, average midrank, exact raw float64 TPM ties',
            'pair_rule': 'all unordered selected-gene pairs; ID ascending orientation; exact raw TPM sign -1/0/+1; no pair filters',
            'constant_pairs_removed': 0, 'high_tie_pairs_removed': 0,
            'feature_construction_precedes_standardization_and_tissue_projection': True,
        }
        if self.kind == 'pair_full':
            n_constant = 0
            high_tie = 0
            for start in range(0, len(pairs), self._BLOCK_SIZE):
                z = FoldFeatureBuilder._encode_pair_block(X, pairs[start:start+self._BLOCK_SIZE])
                n_constant += int(np.count_nonzero(np.ptp(z, axis=0) == 0))
                high_tie += int(np.count_nonzero(np.mean(z == 0, axis=0) > 0.9))
            self.audit_metadata_.update(candidate_pairs=len(pairs), retained_pairs=len(pairs),
                                        training_constant_pairs_retained=n_constant,
                                        training_high_tie_pairs_retained=high_tie)
        return self

    def transform(self, X):
        check_is_fitted(self, 'feature_names_')
        X = self._validate(X, fitting=False)
        if self.kind == 'expression':
            return np.log2(X[:, self.gene_indices_]+1.0)
        if self.kind == 'rank':
            return rankdata(X[:, self.gene_indices_], axis=1, method='average').astype(np.float64)
        z = np.empty((len(X), self.n_output_features_), dtype=np.int8)
        for start in range(0, self.n_output_features_, self._BLOCK_SIZE):
            z[:, start:start+self._BLOCK_SIZE] = FoldFeatureBuilder._encode_pair_block(X, self.pair_indices_[start:start+self._BLOCK_SIZE])
        return z

    def get_feature_names_out(self, input_features=None):
        check_is_fitted(self, 'feature_names_')
        return self.feature_names_.copy()

    def get_audit_metadata(self, *, include_feature_names=False):
        check_is_fitted(self, 'feature_names_')
        result = deepcopy(self.audit_metadata_)
        if include_feature_names:
            result.update(gene_ids=list(self.gene_ids_), gene_indices=self.gene_indices_.tolist(),
                          feature_names=self.feature_names_.tolist(),
                          selection_ranked_gene_ids=list(self.selection_ranked_gene_ids_),
                          selected_gene_scores=self.selected_scores_.tolist(),
                          selected_gene_signed_associations=self.selected_signed_associations_.tolist(),
                          pair_indices_original_input_columns=self.pair_indices_.tolist())
        return result


class DiagnosticRegressor(TrainOnlyRegressor):
    """New builder with frozen v1 EN arithmetic and same-policy tissue null.

    Input is raw TPM, followed by tissue code and audit-only unique sample ID.
    The inherited EN numerical method and prediction method are not modified.
    All molecular arms use alpha_max_relative and training-only FWL. The null
    arm uses exactly their standardized-target tissue means, including the same
    global-mean prediction for previously unseen tissue labels.
    """
    def __init__(self, model_kind='expression', selection_rule='mad', config=None,
                 stable_gene_ids=None, alpha_ratio=1.0, l1_ratio=0.5):
        super().__init__(model_kind=model_kind, config=config, stable_gene_ids=stable_gene_ids,
                         alpha_ratio=alpha_ratio, l1_ratio=l1_ratio)
        self.selection_rule = selection_rule

    def fit(self, X, y):
        for name in list(vars(self)):
            if name.endswith('_'):
                delattr(self, name)
        start = time.perf_counter()
        X = np.asarray(X)
        y = np.asarray(y, dtype=np.float64)
        if X.ndim != 2 or X.dtype != np.dtype(np.float64) or X.shape != (len(y), len(self.stable_gene_ids)+2):
            raise ValueError('Expected float64 TPM + tissue + audit ID columns')
        if y.ndim != 1 or not len(y) or not np.isfinite(y).all() or not np.isfinite(X[:, -2:]).all() or len(np.unique(X[:, -1])) != len(y):
            raise ValueError('Finite aligned outcomes and unique audit IDs are required')
        if self.model_kind not in {'expression', 'rank', 'pair_full', 'tissue_null'}:
            raise ValueError('Unsupported diagnostic model_kind')
        if self.config['elastic_net'].get('regularization_mode') != 'alpha_max_relative':
            raise ValueError('Diagnostic protocol requires alpha_max_relative')
        self.fit_sample_ids_ = X[:, -1].astype(int).tolist()
        self.fit_response_sha256_ = hashlib.sha256(np.ascontiguousarray(y, dtype='<f8').tobytes()).hexdigest()
        self.prediction_sample_ids_ = []
        self.target_mean_ = float(y.mean())
        self.target_sd_ = float(y.std())
        self.fallback_reason_ = ''
        self.warning_messages_ = []
        self.builder_ = None
        self.pair_selector_ = None
        self.n_iter_ = 0
        self.predict_seconds_ = 0.0
        self.last_predict_seconds_ = 0.0
        self.prediction_calls_ = 0
        self.actual_alpha_ = self.alpha_max_ = self.path_ratio_ = self.dual_gap_ = None
        self.null_endpoint_ = False
        self.solver_reason_ = 'not_applicable'
        self.exact_tissue_constant_target_ = None
        self.n_nonzero_molecular_coeffs_ = 0
        self.n_molecular_features_used_ = 0
        self.rf_molecular_feature_names_ = ()
        self.rf_tissue_policy_ = None
        # An exactly repeated float (for example 0.7) can acquire nonzero
        # np.std through mean-reduction roundoff. Detect observed constancy.
        if np.ptp(y) == 0.0:
            self.target_mean_ = float(y[0])
            self.target_sd_ = 1.0
            self.fallback_reason_ = 'constant_training_response'
        ys = (y-self.target_mean_)/self.target_sd_
        self.encoder_ = OneHotEncoder(handle_unknown='ignore', sparse_output=False, dtype=np.float64)
        d = self.encoder_.fit_transform(X[:, -2].reshape(-1, 1))
        self.training_tissue_codes_ = self.encoder_.categories_[0].tolist()
        self.n_features_in_ = X.shape[1]
        self.target_tissue_means_ = (d.T@ys)/d.sum(axis=0)
        if self.model_kind == 'tissue_null':
            self.fitted_kind_ = 'tissue_null'
            self.solver_reason_ = 'same_policy_zero_molecular_tissue_means'
        elif self.fallback_reason_:
            self.fitted_kind_ = 'mean'
        else:
            f = self.config['feature']
            self.builder_ = DiagnosticFeatureBuilder(kind=self.model_kind, selection_rule=self.selection_rule,
                n_genes=f.get('diagnostic_genes', f.get('n_genes', 200)), min_tpm=f['min_tpm'],
                min_detection_fraction=f['min_detection_fraction'], stable_gene_ids=self.stable_gene_ids)
            self.builder_.fit(X[:, :-2], y, tissue=X[:, -2])
            z = self.builder_.transform(X[:, :-2])
            if z.shape[1] == 0:
                self.fitted_kind_ = 'tissue_null'
                self.fallback_reason_ = 'empty_features_same_policy_tissue_null'
                self.solver_reason_ = 'empty_features_same_policy_tissue_null'
            else:
                self.fitted_kind_ = self.model_kind
                self._fit_elastic_net(z, d, ys)
        self.fit_seconds_ = time.perf_counter()-start
        return self

    def predict(self, X):
        check_is_fitted(self, 'fitted_kind_')
        X = np.asarray(X)
        if X.ndim != 2 or X.shape[1] != self.n_features_in_:
            raise ValueError('Input columns do not match fitted design')
        self.prediction_sample_ids_ = X[:, -1].astype(int).tolist()
        if self.fitted_kind_ != 'tissue_null':
            return super().predict(X)
        start = time.perf_counter()
        X = np.asarray(X)
        if X.ndim != 2 or X.shape[1] != self.n_features_in_:
            raise ValueError('Input columns do not match fitted design')
        d = self.encoder_.transform(X[:, -2].reshape(-1, 1))
        pred = self.target_mean_+self.target_sd_*(d@self.target_tissue_means_)
        self.last_predict_seconds_ = time.perf_counter()-start
        self.predict_seconds_ += self.last_predict_seconds_
        self.prediction_calls_ += 1
        return pred

    def feature_audit(self):
        if self.builder_ is None:
            return None
        result = self.builder_.get_audit_metadata(include_feature_names=True)
        active = []
        if hasattr(self, 'regressor_') and hasattr(self.regressor_, 'coef_'):
            for j in np.flatnonzero(self.regressor_.coef_ != 0.0):
                coefficient = float(self.regressor_.coef_[j])
                active.append({'feature_name': str(self.builder_.feature_names_[j]),
                               'coefficient_standardized_feature_standardized_response': coefficient,
                               'coefficient_encoded_feature_original_response': float(self.target_sd_*coefficient/self.scaler_.scale_[j]),
                               'training_feature_scale': float(self.scaler_.scale_[j]),
                               'direction': 'positive' if coefficient > 0 else 'negative'})
        result['active_molecular_coefficients'] = active
        result['coefficient_interpretation'] = 'conditional on tissue; encoded feature means log2(TPM+1), within-sample average rank, or full three-valued pair sign; standardized target uses training standard deviation'
        result['coefficient_zero_rule'] = 'exact coefficient != 0, no threshold'
        return result

    def audit(self):
        result = super().audit()
        result['fit_response_sha256'] = self.fit_response_sha256_
        result['prediction_sample_ids'] = list(self.prediction_sample_ids_)
        result['selected_gene_ids'] = list(self.builder_.gene_ids_) if self.builder_ is not None else []
        result['training_encoded_matrix_shape'] = [len(self.fit_sample_ids_), int(result['n_candidate_molecular_features'] or 0)]
        if hasattr(self, 'regressor_') and hasattr(self.regressor_, 'coef_'):
            active = np.flatnonzero(self.regressor_.coef_ != 0.0)
            result['active_molecular_features'] = [str(self.builder_.feature_names_[j]) for j in active]
            result['active_molecular_coefficients_standardized'] = self.regressor_.coef_[active].tolist()
        else:
            result['active_molecular_features'] = []
            result['active_molecular_coefficients_standardized'] = []
        result['selection_rule'] = self.selection_rule
        result['tissue_policy'] = 'unpenalized training tissue means after FWL; unseen tissue uses training global mean plus any molecular residual prediction'
        result['full_pair_no_filter'] = self.model_kind == 'pair_full'
        return result


def evaluate_diagnostic_outer_fold(X, y, groups, train_idx, test_idx, inner_splits,
                                   model_kind, selection_rule, config, stable_gene_ids,
                                   audit_sink=None):
    """Choose parameters by pooled inner RMSE; evaluate outer outcomes last."""
    X = np.asarray(X)
    y = np.asarray(y)
    groups = np.asarray(groups)
    train_idx = np.asarray(train_idx)
    test_idx = np.asarray(test_idx)
    if set(groups[train_idx]) & set(groups[test_idx]) or set(X[train_idx, -1]) & set(X[test_idx, -1]):
        raise ValueError('Outer biological group or sample-ID overlap')
    xt, yt, gt = X[train_idx], y[train_idx], groups[train_idx]
    events, scores = [], []
    best_score, best_params = float('inf'), None
    def fit_one(indices, params, stage, candidate, inner_fold):
        model = DiagnosticRegressor(model_kind=model_kind, selection_rule=selection_rule,
                                    config=config, stable_gene_ids=stable_gene_ids, **params)
        model.fit(xt[indices], yt[indices])
        event = {**model.audit(), 'stage': stage, 'candidate': candidate, 'inner_fold': inner_fold, 'params': params}
        if set(event['fit_sample_ids']) != set(xt[indices, -1].astype(int)) or set(event['fit_sample_ids']) & set(X[test_idx, -1].astype(int)):
            raise ValueError('Fit membership audit mismatch')
        events.append(event)
        return model
    def finish_prediction(model):
        events[-1].update(model.audit())
        if audit_sink is not None:
            audit_sink(events[-1])
    grid = [{}] if model_kind == 'tissue_null' else candidate_grid('expression', config)
    if model_kind == 'tissue_null':
        best_score, best_params = None, {}
    else:
        for ci, params in enumerate(grid):
            oof = np.full(len(yt), np.nan)
            covered = np.zeros(len(yt), dtype=int)
            fallbacks = []
            for ii, (it, iv) in enumerate(inner_splits):
                it, iv = np.asarray(it), np.asarray(iv)
                if set(gt[it]) & set(gt[iv]) or set(it) & set(iv):
                    raise ValueError('Inner biological group or sample overlap')
                if set(it) | set(iv) != set(range(len(yt))):
                    raise ValueError('Each inner split must partition outer training rows')
                model = fit_one(it, params, 'inner_training', ci, ii)
                oof[iv] = model.predict(xt[iv])
                covered[iv] += 1
                finish_prediction(model)
                if model.fallback_reason_:
                    fallbacks.append(model.fallback_reason_)
            if not np.all(covered == 1) or not np.isfinite(oof).all():
                raise ValueError('Incomplete or nonfinite inner OOF predictions')
            score = float(np.sqrt(mean_squared_error(yt, oof)))
            scores.append({'candidate': ci, 'params': params, 'pooled_inner_rmse': score, 'fallbacks': fallbacks})
            if score < best_score:
                best_score, best_params = score, params.copy()
    model = fit_one(np.arange(len(yt)), best_params, 'outer_training_refit', -1, -1)
    pred = model.predict(X[test_idx])
    finish_prediction(model)
    boundary_details = _grid_boundary_details(best_params, grid)
    return {'chosen_params': best_params, 'best_inner_rmse': best_score, 'predictions': pred,
            'metrics': metrics(y[test_idx], pred), 'inner_scores': scores,
            'feature_audit': model.feature_audit(), 'fit_events': events,
            'outer_unknown_tissue_samples': int((~np.isin(X[test_idx, -2], model.training_tissue_codes_)).sum()),
            'chosen_grid_boundary': bool(any(d['at_boundary'] and d['search_extent_diagnostic'] for d in boundary_details.values())),
            'chosen_grid_boundary_details': boundary_details, 'outer_fit_audit': model.audit()}
