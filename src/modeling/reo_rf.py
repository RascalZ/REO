"""Training-only pair association selection for the REO random forest.

The association calculation controls training tissue means. The selected
columns are then returned in their original three-value representation: the
random forest never receives tissue-residualized TPM or REO observations.
This module performs no file access and is fitted afresh in every inner fit.
"""
from __future__ import annotations

from copy import deepcopy
from numbers import Integral

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted


class TrainOnlyReoPairSelector(TransformerMixin, BaseEstimator):
    """Rank REO pairs by training within-tissue Pearson association.

    ``fit(Z, y, tissue_design=D, feature_names=names)`` accepts the original REO
    matrix and a full training-fitted one-hot tissue design. Pearson correlation
    is computed after projecting both variables off the column space of D.
    Ranking is absolute correlation descending, then stable pair name ascending.
    Undefined associations (zero within-tissue pair variance) are excluded.
    A constant tissue-adjusted response retains zero pairs and records that
    reason; the enclosing random forest may then use the tissue columns alone.

    ``selected_indices_`` indexes the input Z columns, ``selected_feature_names_``
    gives the corresponding stable pair names. ``transform`` selects original
    columns without residualizing their values. Audit contents are JSON types.
    """
    _BLOCK_SIZE = 512

    def __init__(self, max_pairs=500, association="within_tissue_pearson"):
        self.max_pairs = max_pairs
        self.association = association

    def fit(self, Z, y, *, tissue_design, feature_names):
        for name in list(vars(self)):
            if name.endswith("_"):
                delattr(self, name)
        if isinstance(self.max_pairs, bool) or not isinstance(self.max_pairs, Integral) or self.max_pairs < 1:
            raise ValueError("max_pairs must be a positive integer.")
        if self.association != "within_tissue_pearson":
            raise ValueError("Only training within_tissue_pearson association is supported.")
        Z = np.asarray(Z)
        y = np.asarray(y, dtype=np.float64)
        D = np.asarray(tissue_design, dtype=np.float64)
        names = tuple(feature_names)
        if Z.ndim != 2 or y.shape != (Z.shape[0],) or D.ndim != 2 or D.shape[0] != len(y):
            raise ValueError("REO, target and tissue design must have aligned training rows.")
        if len(y) == 0 or D.shape[1] == 0 or len(names) != Z.shape[1] or len(set(names)) != len(names):
            raise ValueError("Non-empty aligned rows, tissue columns and unique feature names are required.")
        if not np.isfinite(Z).all() or not np.isfinite(y).all() or not np.isfinite(D).all():
            raise ValueError("Association inputs must be finite.")
        if not np.isin(D, (0.0, 1.0)).all() or not np.all(D.sum(axis=1) == 1.0):
            raise ValueError("tissue_design must be a full one-hot training design.")
        counts = D.sum(axis=0)
        if np.any(counts == 0):
            raise ValueError("Training tissue design cannot contain an empty tissue column.")

        yr = y - D @ ((D.T @ y) / counts)
        y_norm = float(np.linalg.norm(yr))
        correlation = np.zeros(Z.shape[1], dtype=np.float64)
        eligible = np.zeros(Z.shape[1], dtype=bool)
        reason = ""
        if y_norm == 0.0:
            reason = "constant_tissue_adjusted_training_response"
        else:
            for start in range(0, Z.shape[1], self._BLOCK_SIZE):
                stop = min(start + self._BLOCK_SIZE, Z.shape[1])
                block = np.asarray(Z[:, start:stop], dtype=np.float64)
                residual = block - D @ ((D.T @ block) / counts[:, None])
                norms = np.linalg.norm(residual, axis=0)
                valid = norms > 0.0
                eligible[start:stop] = valid
                local = np.zeros(stop - start, dtype=np.float64)
                local[valid] = (residual[:, valid].T @ yr) / (norms[valid] * y_norm)
                correlation[start:stop] = np.clip(local, -1.0, 1.0)
        candidates = np.flatnonzero(eligible)
        names_array = np.asarray(names, dtype=str)
        order = np.lexsort((names_array[candidates], -np.abs(correlation[candidates])))
        selected = candidates[order[: self.max_pairs]]
        if selected.size == 0 and not reason:
            reason = "no_pairs_with_training_within_tissue_variance"
        self.n_features_in_ = int(Z.shape[1])
        self.selected_indices_ = selected.astype(np.int64, copy=True)
        self.selected_feature_names_ = tuple(names[i] for i in selected)
        self.signed_associations_ = correlation[selected].copy()
        self.fallback_reason_ = reason
        self.audit_metadata_ = {
            "association": self.association,
            "association_description": "absolute Pearson correlation after training-only projection of REO and target off full tissue one-hot design",
            "training_samples": int(len(y)),
            "candidate_pairs": int(Z.shape[1]),
            "association_eligible_pairs": int(candidates.size),
            "requested_max_pairs": int(self.max_pairs),
            "selected_pairs": int(selected.size),
            "rank_tie_rule": "stable gene-pair feature name ascending",
            "selected_feature_indices": selected.tolist(),
            "selected_feature_names": list(self.selected_feature_names_),
            "selected_signed_associations": self.signed_associations_.tolist(),
            "fallback_reason": reason,
            "tree_feature_values": "original REO codes; no residualization",
        }
        return self

    def transform(self, Z):
        check_is_fitted(self, "selected_indices_")
        Z = np.asarray(Z)
        if Z.ndim != 2 or Z.shape[1] != self.n_features_in_:
            raise ValueError("REO transform columns differ from fitted candidate pairs.")
        if not np.isfinite(Z).all():
            raise ValueError("REO transform input must be finite.")
        return Z[:, self.selected_indices_]

    def get_feature_names_out(self, input_features=None):
        check_is_fitted(self, "selected_indices_")
        if input_features is not None and len(input_features) != self.n_features_in_:
            raise ValueError("input_features must match fitted candidate pair columns.")
        return np.asarray(self.selected_feature_names_, dtype=str)

    def get_audit_metadata(self):
        check_is_fitted(self, "selected_indices_")
        return deepcopy(self.audit_metadata_)
