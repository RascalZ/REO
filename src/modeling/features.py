"""Leakage-safe expression and relative-expression-order features.

Every distribution-dependent decision is learned by ``fit`` from the rows it
receives. Put this transformer inside the pipeline passed to the inner-CV
search; fitting it before that search would still leak inner-validation data.
The input columns must have been restricted by annotation, not by expression
distributions, and must retain the same stable gene IDs and order in each fold.
"""

from __future__ import annotations

from copy import deepcopy
from numbers import Integral, Real

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted


class FoldFeatureBuilder(TransformerMixin, BaseEstimator):
    """Select genes in a training fold and create frozen expression/REO features.

    Parameters
    ----------
    kind : {"expression", "reo"}, default="expression"
        Expression output is log2(TPM + 1); REO output is the exact ordering of
        raw TPM for each retained pair, using -1, 0, +1 for lower, tie, higher.
    n_genes : int, default=2000
        Maximum number of training-selected genes. If fewer genes meet the
        filter, all eligible genes are retained and the shortfall is recorded.
    min_tpm : float, default=1.0
        A training observation is detected when its raw TPM is >= this value.
    min_detection_fraction : float, default=0.1
        Minimum training-only detected fraction for a gene to be eligible.
    max_pair_tie_fraction : float, default=0.9
        REO pairs with a larger training tie fraction are removed. Constant
        pairs are always removed, including pairs with all-zero observations.
    stable_gene_ids : tuple of str, required
        Unique, non-empty identifiers aligned exactly with the input columns.
        They resolve equal-MAD gene ranks deterministically and define pair
        orientation. Do not select this tuple using expression distributions.

    Notes
    -----
    ``gene_indices_`` indexes the original input columns. For REO these indices
    are ordered by stable gene ID after selection; expression columns remain
    in decreasing training MAD order. ``pair_indices_`` also indexes original
    input columns (not positions within ``gene_indices_``).

    No response-dependent feature selection is performed by this transformer;
    ``y`` is deliberately ignored. Both the unsupervised distribution filters
    and all downstream scaling/model selection nevertheless belong inside
    each inner-training fit. The public audit dictionaries contain only JSON
    primitive types, and contain no individual observations or response values.
    """

    _GENE_BLOCK_SIZE = 512
    _PAIR_BLOCK_SIZE = 512

    def __init__(
        self,
        kind="expression",
        n_genes=2000,
        min_tpm=1.0,
        min_detection_fraction=0.1,
        max_pair_tie_fraction=0.9,
        stable_gene_ids=None,
    ):
        # Leave constructor parameters unchanged for sklearn.clone/GridSearchCV.
        self.kind = kind
        self.n_genes = n_genes
        self.min_tpm = min_tpm
        self.min_detection_fraction = min_detection_fraction
        self.max_pair_tie_fraction = max_pair_tie_fraction
        self.stable_gene_ids = stable_gene_ids

    def _validate_parameters(self):
        if self.kind not in {"expression", "reo"}:
            raise ValueError("kind must be 'expression' or 'reo'.")
        if isinstance(self.n_genes, bool) or not isinstance(self.n_genes, Integral):
            raise ValueError("n_genes must be a positive integer.")
        if self.n_genes <= 0:
            raise ValueError("n_genes must be a positive integer.")
        for name, value, upper in (
            ("min_tpm", self.min_tpm, None),
            ("min_detection_fraction", self.min_detection_fraction, 1.0),
            ("max_pair_tie_fraction", self.max_pair_tie_fraction, 1.0),
        ):
            if (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not np.isfinite(value)
                or value < 0
                or (upper is not None and value > upper)
            ):
                limit = "a finite non-negative number" if upper is None else "in [0, 1]"
                raise ValueError(f"{name} must be {limit}.")
        if not isinstance(self.stable_gene_ids, tuple):
            raise ValueError("stable_gene_ids must be a tuple aligned with the input columns.")
        if any(not isinstance(g, str) or not g.strip() for g in self.stable_gene_ids):
            raise ValueError("stable_gene_ids must contain non-empty strings.")
        if len(set(self.stable_gene_ids)) != len(self.stable_gene_ids):
            raise ValueError("stable_gene_ids must be unique.")

    def _validate_X(self, X, *, fitting):
        X = np.asarray(X)
        if X.ndim != 2:
            raise ValueError("X must be a two-dimensional raw-TPM matrix.")
        if X.dtype != np.dtype(np.float64):
            raise ValueError("X must have dtype float64; preserve raw-TPM ordering precision.")
        if X.shape[0] == 0 or X.shape[1] == 0:
            raise ValueError("X must contain at least one sample and one gene.")
        expected = len(self.stable_gene_ids) if fitting else self.n_features_in_
        if X.shape[1] != expected:
            raise ValueError(f"X has {X.shape[1]} columns; expected {expected} aligned gene columns.")
        # Check all input columns, including genes not selected in a previous fit.
        for start in range(0, X.shape[1], self._GENE_BLOCK_SIZE):
            block = X[:, start : start + self._GENE_BLOCK_SIZE]
            if not np.isfinite(block).all():
                raise ValueError("X contains non-finite raw-TPM values.")
            if np.any(block < 0):
                raise ValueError("X contains negative raw-TPM values.")
        return X

    @staticmethod
    def _encode_pair_block(X, pair_indices):
        left = X[:, pair_indices[:, 0]]
        right = X[:, pair_indices[:, 1]]
        # Comparison avoids both subtraction overflow and tolerance-based ties.
        encoded = np.greater(left, right).astype(np.int8)
        encoded -= np.less(left, right).astype(np.int8)
        return encoded

    @staticmethod
    def _summary(values):
        return {
            "minimum": float(np.min(values)),
            "median": float(np.median(values)),
            "maximum": float(np.max(values)),
        }

    def fit(self, X, y=None):
        """Learn filters, gene ranks and pair filters using these training rows."""
        # A failed refit must not leave the estimator looking fitted to old data.
        for name in list(vars(self)):
            if name.endswith("_"):
                delattr(self, name)
        self._validate_parameters()
        X = self._validate_X(X, fitting=True)

        n_samples, n_features = X.shape
        detection = np.empty(n_features, dtype=np.float64)
        mad = np.empty(n_features, dtype=np.float64)
        for start in range(0, n_features, self._GENE_BLOCK_SIZE):
            stop = min(start + self._GENE_BLOCK_SIZE, n_features)
            block = X[:, start:stop]
            detection[start:stop] = np.mean(block >= self.min_tpm, axis=0)
            log_block = np.log2(block + 1.0)
            medians = np.median(log_block, axis=0)
            mad[start:stop] = np.median(np.abs(log_block - medians), axis=0)

        detected = detection >= self.min_detection_fraction
        eligible_indices = np.flatnonzero(detected & (mad > 0.0))
        if eligible_indices.size == 0:
            raise ValueError("No genes pass training detection and positive log2TPM MAD filters.")
        ids = np.asarray(self.stable_gene_ids, dtype=str)
        rank = np.lexsort((ids[eligible_indices], -mad[eligible_indices]))
        selected = eligible_indices[rank[: self.n_genes]]
        if self.kind == "reo":
            selected = selected[np.argsort(ids[selected], kind="stable")]

        # Build the local result before publishing fitted attributes.
        actual_n_genes = int(selected.size)
        stats = {
            "kind": self.kind,
            "training_samples": int(n_samples),
            "input_genes": int(n_features),
            "requested_genes": int(self.n_genes),
            "selected_genes": actual_n_genes,
            "gene_shortfall": int(max(0, self.n_genes - actual_n_genes)),
            "detection_eligible_genes": int(np.count_nonzero(detected)),
            "detection_and_positive_mad_eligible_genes": int(eligible_indices.size),
            "min_tpm": float(self.min_tpm),
            "min_detection_fraction": float(self.min_detection_fraction),
            "max_pair_tie_fraction": float(self.max_pair_tie_fraction),
            "gene_rank_rule": "training log2(TPM+1) MAD descending; stable gene ID ascending for ties",
            "response_used_for_feature_construction": False,
            "selected_gene_detection_fraction": self._summary(detection[selected]),
            "selected_gene_log2tpm_mad": self._summary(mad[selected]),
        }

        if self.kind == "reo":
            if actual_n_genes < 2:
                raise ValueError("REO requires at least two training-eligible genes.")
            left_pos, right_pos = np.triu_indices(actual_n_genes, k=1)
            candidate_pairs = np.column_stack((selected[left_pos], selected[right_pos]))
            n_candidates = int(candidate_pairs.shape[0])
            tie_fraction = np.empty(n_candidates, dtype=np.float64)
            nonconstant = np.empty(n_candidates, dtype=bool)
            for start in range(0, n_candidates, self._PAIR_BLOCK_SIZE):
                stop = min(start + self._PAIR_BLOCK_SIZE, n_candidates)
                codes = self._encode_pair_block(X, candidate_pairs[start:stop])
                tie_fraction[start:stop] = np.mean(codes == 0, axis=0)
                nonconstant[start:stop] = np.min(codes, axis=0) != np.max(codes, axis=0)
            passes_tie = tie_fraction <= self.max_pair_tie_fraction
            retained = nonconstant & passes_tie
            pairs = candidate_pairs[retained]
            if pairs.shape[0] == 0:
                raise ValueError("No REO pairs pass training nonconstant and tie-fraction filters.")
            names = np.asarray(
                [f"REO:{ids[a]}>{ids[b]}" for a, b in pairs], dtype=str
            )
            stats.update(
                {
                    "candidate_pairs": n_candidates,
                    "constant_pairs_removed": int(np.count_nonzero(~nonconstant)),
                    "high_tie_nonconstant_pairs_removed": int(np.count_nonzero(nonconstant & ~passes_tie)),
                    "retained_pairs": int(pairs.shape[0]),
                    "retained_pair_training_tie_fraction": self._summary(tie_fraction[retained]),
                    "pair_orientation": "stable gene ID ascending; raw TPM(first) > TPM(second) gives +1",
                    "tie_rule": "exact raw-TPM equality gives 0; no rounding or tolerance",
                    "output_dtype": "int8",
                }
            )
        else:
            pairs = np.empty((0, 2), dtype=np.int64)
            names = np.asarray([f"log2TPM:{g}" for g in ids[selected]], dtype=str)
            tie_fraction = np.empty(0, dtype=np.float64)
            retained = np.empty(0, dtype=bool)
            stats["output_dtype"] = "float64"

        self.n_features_in_ = int(n_features)
        self.n_samples_seen_ = int(n_samples)
        self.gene_indices_ = selected.astype(np.int64, copy=True)
        self.gene_ids_ = tuple(self.stable_gene_ids[i] for i in selected)
        self.selected_gene_mad_ = mad[selected].copy()
        self.selected_gene_detection_fraction_ = detection[selected].copy()
        self.pair_indices_ = pairs.astype(np.int64, copy=True)
        self.pair_tie_fraction_ = tie_fraction[retained].copy()
        self.feature_names_ = names
        self.n_output_features_ = int(names.size)
        stats["output_features"] = self.n_output_features_
        self.training_statistics_ = stats
        self.audit_metadata_ = deepcopy(stats)
        return self

    def transform(self, X):
        """Apply the fitted features without learning from these rows."""
        check_is_fitted(self, "feature_names_")
        X = self._validate_X(X, fitting=False)
        if self.kind == "expression":
            return np.log2(X[:, self.gene_indices_] + 1.0)
        output = np.empty((X.shape[0], self.n_output_features_), dtype=np.int8)
        for start in range(0, self.n_output_features_, self._PAIR_BLOCK_SIZE):
            stop = min(start + self._PAIR_BLOCK_SIZE, self.n_output_features_)
            output[:, start:stop] = self._encode_pair_block(X, self.pair_indices_[start:stop])
        return output

    def get_feature_names_out(self, input_features=None):
        """Return the ordered stable feature names learned in the training fold."""
        check_is_fitted(self, "feature_names_")
        if input_features is not None and len(input_features) != self.n_features_in_:
            raise ValueError("input_features must have the same length as the fitted input columns.")
        return self.feature_names_.copy()

    def get_audit_metadata(self, *, include_feature_names=False):
        """Return a JSON-serializable audit copy, optionally with the feature list."""
        check_is_fitted(self, "feature_names_")
        result = deepcopy(self.audit_metadata_)
        if include_feature_names:
            result["gene_ids"] = list(self.gene_ids_)
            result["gene_indices"] = self.gene_indices_.tolist()
            result["feature_names"] = self.feature_names_.tolist()
            result["pair_indices_original_input_columns"] = self.pair_indices_.tolist()
        return result
