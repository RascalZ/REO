"""Response-blind P6 input gateway and fixed-model external statistics.

This module accepts in-memory, caller-supplied tables only. It neither locates
nor opens an external response file. In particular, ``prepare_gateway`` cannot
fit, select genes, impute TPM, or inspect outcomes. A separate, locked caller
must prove the source/identity contracts before using these helpers.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
import hashlib
from itertools import combinations
from math import ceil
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analyze_formal_results import generate_global_bootstrap_weights, weighted_metrics


PREDICTED = "predicted"
ABSTAINED = "abstained"


def names_sha256(names: Sequence[str]) -> str:
    """Match the frozen feature builder's NUL-delimited name digest."""
    return hashlib.sha256("\0".join(names).encode("utf-8")).hexdigest()


def _unique_names(names: Sequence[str], label: str) -> tuple[str, ...]:
    names = tuple(names)
    if not names or any(not isinstance(x, str) or not x.strip() for x in names):
        raise ValueError(f"{label} must contain nonempty string IDs")
    if len(set(names)) != len(names):
        raise ValueError(f"{label} contains duplicate IDs")
    return names


def _feature_names(selected: tuple[str, ...], representation: str) -> tuple[str, ...]:
    if representation == "expression":
        return tuple("log2TPM:" + gene for gene in selected)
    if representation == "rank":
        return tuple("sample_rank:" + gene for gene in selected)
    if representation == "pair_full":
        return tuple(f"REO:{left}>{right}" for left, right in combinations(selected, 2))
    if representation == "tissue_null":
        return ()
    raise ValueError("Unknown representation")


def validate_model_contract(model, *, representation: str,
                            stable_gene_ids: Sequence[str],
                            training_tissue_codes: Sequence[int],
                            expected_gene_sha256: str | None = None,
                            expected_selected_gene_ids: Sequence[str] | None = None,
                            expected_feature_names: Sequence[str] | None = None,
                            expected_feature_sha256: str | None = None) -> None:
    """Reject reordered genes/features and a model built for another design.

    A fitted constant-target or empty-feature model may have no builder. Its
    zero molecular signal remains a *valid* model, provided its input schema and
    fallback state agree with the frozen contract.
    """
    genes = _unique_names(stable_gene_ids, "Stable gene order")
    if expected_gene_sha256 is not None and names_sha256(genes) != expected_gene_sha256:
        raise ValueError("Stable gene-order hash mismatch")
    if tuple(getattr(model, "stable_gene_ids", ())) != genes:
        raise ValueError("Fitted model stable gene order mismatch")
    if getattr(model, "n_features_in_", None) != len(genes) + 2:
        raise ValueError("Fitted model input width mismatch")
    if getattr(model, "model_kind", None) != representation:
        raise ValueError("Fitted model representation mismatch")
    trained = tuple(getattr(model, "training_tissue_codes_", ()))
    if trained != tuple(training_tissue_codes):
        raise ValueError("Fitted model training tissue code mismatch")
    builder = getattr(model, "builder_", None)
    selected = tuple(getattr(builder, "gene_ids_", ())) if builder is not None else ()
    if len(set(selected)) != len(selected) or selected != tuple(sorted(selected)):
        raise ValueError("Selected gene order is not unique ascending stable ID")
    if any(gene not in genes for gene in selected):
        raise ValueError("Selected gene absent from the frozen input order")
    if builder is not None:
        positions = {gene: j for j, gene in enumerate(genes)}
        expected_indices = tuple(positions[gene] for gene in selected)
        if tuple(map(int, getattr(builder, "gene_indices_", ()))) != expected_indices:
            raise ValueError("Selected gene indices disagree with frozen raw gene order")
        if representation == "pair_full":
            expected_pairs = tuple(combinations(expected_indices, 2))
            actual_pairs = tuple(map(tuple, np.asarray(getattr(builder, "pair_indices_", ())).tolist()))
            if actual_pairs != expected_pairs:
                raise ValueError("Pair input indices/order disagree with frozen genes")
    if expected_selected_gene_ids is not None and selected != tuple(expected_selected_gene_ids):
        raise ValueError("Selected gene order differs from saved model manifest")
    actual_features = tuple(map(str, getattr(builder, "feature_names_", ()))) if builder is not None else ()
    if actual_features != _feature_names(selected, representation):
        raise ValueError("Encoded feature order/direction mismatch")
    if expected_feature_names is not None and actual_features != tuple(expected_feature_names):
        raise ValueError("Feature order differs from saved model manifest")
    if expected_feature_sha256 is not None and names_sha256(actual_features) != expected_feature_sha256:
        raise ValueError("Feature-order hash differs from saved model manifest")
    if builder is not None and tuple(getattr(builder, "stable_gene_ids", ())) != genes:
        raise ValueError("Feature builder raw gene order mismatch")
    if representation == "tissue_null" and builder is not None:
        raise ValueError("Tissue-only model unexpectedly carries molecular features")
    if representation != "tissue_null" and builder is None:
        if getattr(model, "fitted_kind_", None) not in ("mean", "tissue_null"):
            raise ValueError("Absent molecular builder without documented fallback")


@dataclass(frozen=True)
class GatewayBatch:
    identity: pd.DataFrame
    numeric_tpm: np.ndarray
    tissue_codes: np.ndarray
    audit_ids: np.ndarray
    identity_valid: np.ndarray
    molecular_valid: np.ndarray
    abstain_reason: tuple[str, ...]
    unknown_tissue: np.ndarray
    stable_gene_ids: tuple[str, ...]
    training_tissue_codes: tuple[int, ...]
    unknown_tissue_code: int


def prepare_gateway(identity: pd.DataFrame, tpm: pd.DataFrame, *,
                    stable_gene_ids: Sequence[str], tissue_to_code: Mapping[str, int],
                    training_tissue_codes: Sequence[int],
                    unknown_tissue_code: int | None = None,
                    expected_gene_sha256: str | None = None,
                    forbidden_training_canonical_ids: Sequence[str] = (),
                    forbidden_training_group_ids: Sequence[str] = ()) -> GatewayBatch:
    """Validate full TPM vector and produce explicit molecular abstentions.

    ``identity`` is the predeclared eligible sample roster (canonical model ID,
    biological group ID, tissue). ``tpm`` is a raw TPM table indexed by those
    canonical IDs. Missing *rows* are recorded; unexpected rows or any missing,
    reordered, duplicated, or extra *gene columns* are structural errors.
    ``tissue_to_code`` must contain exactly the tissues observed in this
    *drug's fitted model*, using their original global-development integer
    codes. It must not include global cohort tissues absent from this drug's
    training rows. Unknown but nonblank tissues, including those seen globally
    but absent for this drug, get a fresh code and remain evaluable.
    """
    required = {"canonical_model_id", "biological_group_id", "tissue"}
    if not required.issubset(identity.columns):
        raise ValueError("Identity roster lacks required columns")
    ident = identity.loc[:, ["canonical_model_id", "biological_group_id", "tissue"]].copy().reset_index(drop=True)
    ids = _unique_names(ident.canonical_model_id.tolist(), "Canonical model IDs")
    if set(ids) & set(forbidden_training_canonical_ids):
        raise ValueError("External canonical model overlaps frozen training identity")
    if not isinstance(tpm, pd.DataFrame):
        raise TypeError("TPM input must be a DataFrame with named ordered gene columns")
    genes = _unique_names(stable_gene_ids, "Stable gene order")
    if expected_gene_sha256 is not None and names_sha256(genes) != expected_gene_sha256:
        raise ValueError("Stable gene-order hash mismatch")
    if tuple(tpm.columns) != genes:
        raise ValueError("Raw TPM columns must exactly match frozen stable gene order")
    if not tpm.index.is_unique or any(not isinstance(x, str) or not x.strip() for x in tpm.index):
        raise ValueError("TPM row IDs must be unique nonempty canonical strings")
    if set(tpm.index) - set(ids):
        raise ValueError("TPM contains an unexpected canonical model ID")
    if not tissue_to_code or any(not isinstance(k, str) or not k.strip() for k in tissue_to_code):
        raise ValueError("Invalid training tissue name map")
    mapped_codes = tuple(tissue_to_code.values())
    if any(isinstance(v, bool) or not isinstance(v, (int, np.integer)) for v in mapped_codes):
        raise ValueError("Tissue codes must be distinct integers")
    if len(set(mapped_codes)) != len(mapped_codes):
        raise ValueError("Training tissue code collision")
    if any(abs(int(v)) >= 2**53 for v in mapped_codes):
        raise ValueError("Tissue codes must remain exact as float64")
    trained = tuple(training_tissue_codes)
    if set(mapped_codes) != set(trained) or len(trained) != len(set(trained)):
        raise ValueError("Tissue map and fitted training category codes disagree")
    if unknown_tissue_code is None:
        unknown_tissue_code = max(mapped_codes) + 1
    if isinstance(unknown_tissue_code, bool) or not isinstance(unknown_tissue_code, (int, np.integer)):
        raise ValueError("Unknown tissue code must be an integer")
    if unknown_tissue_code in set(trained):
        raise ValueError("Unknown tissue code collides with a trained category")
    if abs(int(unknown_tissue_code)) >= 2**53:
        raise ValueError("Unknown tissue code cannot be represented exactly as float64")
    tissue = ident.tissue.tolist()
    valid_group = [isinstance(x, str) and bool(x.strip()) for x in ident.biological_group_id]
    valid_groups = [ident.biological_group_id.iloc[i] for i, valid in enumerate(valid_group) if valid]
    if len(set(valid_groups)) != len(valid_groups):
        raise ValueError("Duplicate biological group in external identity roster")
    if set(valid_groups) & set(forbidden_training_group_ids):
        raise ValueError("External biological group overlaps frozen training/candidate groups")
    valid_tissue = [isinstance(x, str) and bool(x.strip()) for x in tissue]
    identity_valid = np.asarray([a and b for a, b in zip(valid_group, valid_tissue)], dtype=bool)
    unknown = np.asarray([valid_tissue[i] and name not in tissue_to_code for i, name in enumerate(tissue)], dtype=bool)
    codes = np.asarray([tissue_to_code.get(name, unknown_tissue_code) if valid_tissue[i]
                        else unknown_tissue_code for i, name in enumerate(tissue)], dtype=np.float64)
    # A unique, exactly representable local audit ID, never a response-derived ID.
    if len(ids) >= 2**53:
        raise ValueError("Audit IDs cannot be represented exactly as float64")
    audit_ids = np.arange(len(ids), dtype=np.float64)
    present = np.asarray([name in tpm.index for name in ids], dtype=bool)
    raw = tpm.reindex(ids).to_numpy()
    missing = pd.isna(raw)
    if np.issubdtype(raw.dtype, np.number):
        numeric = raw.astype(np.float64, copy=False)
    else:
        numeric = pd.to_numeric(pd.Series(raw.ravel()), errors="coerce").to_numpy(dtype=np.float64).reshape(raw.shape)
    nonnumeric = (~missing) & np.isnan(numeric)
    nonfinite = (~missing) & (~nonnumeric) & (~np.isfinite(numeric))
    negative = np.isfinite(numeric) & (numeric < 0)
    reasons: list[str] = []
    for i in range(len(ids)):
        problems = []
        if not valid_group[i]:
            problems.append("invalid_biological_group")
        if not valid_tissue[i]:
            problems.append("invalid_tissue")
        if not present[i]:
            problems.append("missing_expression_row")
        else:
            if missing[i].any():
                problems.append("missing_tpm")
            if nonnumeric[i].any():
                problems.append("nonnumeric_tpm")
            if nonfinite[i].any():
                problems.append("nonfinite_tpm")
            if negative[i].any():
                problems.append("negative_tpm")
        reasons.append(";".join(problems))
    valid = np.asarray([not reason for reason in reasons], dtype=bool)
    return GatewayBatch(ident, numeric, codes, audit_ids, identity_valid, valid,
                        tuple(reasons), unknown, genes, trained, int(unknown_tissue_code))


def predict_external_arm(model, batch: GatewayBatch, *, drug_id: int, arm: str,
                         representation: str,
                         expected_selected_gene_ids: Sequence[str] | None = None,
                         expected_feature_names: Sequence[str] | None = None,
                         expected_feature_sha256: str | None = None,
                         artifact_sha256: str | None = None) -> pd.DataFrame:
    """Predict only eligible rows; never pass bad TPM to a molecular model."""
    if representation == "tissue_null" and arm != "tissue_null":
        raise ValueError("Tissue-null representation assigned to a molecular arm")
    if representation != "tissue_null" and arm == "tissue_null":
        raise ValueError("Molecular representation assigned to tissue_null")
    if artifact_sha256 is not None and (len(artifact_sha256) != 64 or
            any(character not in "0123456789abcdef" for character in artifact_sha256)):
        raise ValueError("Artifact SHA-256 must be lowercase hexadecimal")
    validate_model_contract(model, representation=representation,
                            stable_gene_ids=batch.stable_gene_ids,
                            training_tissue_codes=batch.training_tissue_codes,
                            expected_selected_gene_ids=expected_selected_gene_ids,
                            expected_feature_names=expected_feature_names,
                            expected_feature_sha256=expected_feature_sha256)
    valid = batch.identity_valid if representation == "tissue_null" else batch.molecular_valid
    pred = np.full(len(batch.identity), np.nan, dtype=np.float64)
    if valid.any():
        if representation == "tissue_null":
            # Null fit ignores molecular columns; placeholders are *only* sent
            # to this audited null model, never interpreted as measured zeros.
            molecular = np.zeros((int(valid.sum()), len(batch.stable_gene_ids)), dtype=np.float64)
        else:
            molecular = batch.numeric_tpm[valid, :]
            if not np.isfinite(molecular).all() or (molecular < 0).any():
                raise AssertionError("Molecular abstention gateway passed invalid TPM")
        design = np.column_stack((molecular, batch.tissue_codes[valid], batch.audit_ids[valid]))
        obtained = np.asarray(model.predict(design), dtype=np.float64)
        if obtained.shape != (int(valid.sum()),):
            raise ValueError("Fitted model returned wrong prediction shape")
        if not np.isfinite(obtained).all():
            raise ValueError("Nonfinite fitted-model prediction is a hard external gate failure")
        pred[valid] = obtained
    status = np.where(valid & np.isfinite(pred), PREDICTED, ABSTAINED)
    reason = []
    for i in range(len(pred)):
        if status[i] == PREDICTED:
            reason.append("")
        elif not batch.identity_valid[i]:
            reason.append(";".join(x for x in batch.abstain_reason[i].split(";")
                                   if x in {"invalid_biological_group", "invalid_tissue"}))
        elif representation != "tissue_null" and not batch.molecular_valid[i]:
            reason.append(batch.abstain_reason[i])
        else:
            reason.append("nonfinite_model_prediction")
    n_nonzero = getattr(model, "n_nonzero_molecular_coeffs_", None)
    zero_molecular_solution = (None if representation == "tissue_null" or n_nonzero is None
                               else int(n_nonzero) == 0)
    return pd.DataFrame({"drug_id": int(drug_id), "model": arm,
                         "canonical_model_id": batch.identity.canonical_model_id,
                         "biological_group_id": batch.identity.biological_group_id,
                         "tissue": batch.identity.tissue,
                         "eligible_denominator_member": True,
                         "predicted_auc": pred, "status": status,
                         "abstain_reason": reason,
                         "unknown_training_tissue": batch.unknown_tissue,
                         "zero_molecular_solution": zero_molecular_solution,
                         "artifact_sha256": artifact_sha256})


def spearman_checked(y: Sequence[float], prediction: Sequence[float],
                     biological_group_ids: Sequence[str], *, min_groups: int = 3) -> tuple[float | None, str]:
    """Return a signed, unreversed Spearman coefficient or an explicit NA."""
    y = np.asarray(y, dtype=np.float64)
    p = np.asarray(prediction, dtype=np.float64)
    g = tuple(biological_group_ids)
    if y.ndim != 1 or p.shape != y.shape or len(g) != len(y):
        raise ValueError("Outcome, prediction, and group arrays must align")
    if len(set(g)) != len(g):
        raise ValueError("Spearman input repeats an independent biological group")
    if not np.isfinite(y).all() or not np.isfinite(p).all():
        return None, "nonfinite_response_or_prediction"
    if len(set(g)) < min_groups:
        return None, "insufficient_independent_groups"
    if np.ptp(y) == 0:
        return None, "constant_response"
    if np.ptp(p) == 0:
        return None, "constant_prediction"
    value = float(spearmanr(y, p).statistic)
    return (value, "") if np.isfinite(value) else (None, "undefined_spearman")


def paired_common_table(responses: pd.DataFrame, predictions: pd.DataFrame, *,
                        drug_id: int, arm_a: str, arm_b: str,
                        design_min_groups: int = 80,
                        frozen_eligible_ids: Sequence[str] | None = None) -> tuple[pd.DataFrame, dict]:
    """Join two arms on exact canonical ID and the same single response label."""
    response_columns = {"drug_id", "canonical_model_id", "biological_group_id", "observed_auc"}
    prediction_columns = {"drug_id", "model", "canonical_model_id", "biological_group_id", "predicted_auc", "status"}
    if not response_columns.issubset(responses) or not prediction_columns.issubset(predictions):
        raise ValueError("External response/prediction table lacks required columns")
    base = responses[responses.drug_id.eq(drug_id)].copy()
    if base.canonical_model_id.duplicated().any() or base.biological_group_id.duplicated().any():
        raise ValueError("External response roster must contain one ID per independent group")
    _unique_names(base.canonical_model_id.tolist(), "External response canonical IDs")
    _unique_names(base.biological_group_id.tolist(), "External response biological groups")
    if not np.isfinite(pd.to_numeric(base.observed_auc, errors="coerce").to_numpy(dtype=float)).all():
        raise ValueError("Eligible external response roster contains nonfinite label")
    response_ids = set(base.canonical_model_id)
    frozen_ids = (_unique_names(frozen_eligible_ids, "Frozen eligible external IDs")
                  if frozen_eligible_ids is not None else tuple(base.canonical_model_id))
    if not response_ids.issubset(set(frozen_ids)):
        raise ValueError("External response ID absent from frozen eligible roster")
    ids = set(frozen_ids)
    indexed = base.set_index("canonical_model_id")
    by_arm = {}
    reason_counts = {}
    all_frozen_predicted = {}
    for arm in (arm_a, arm_b):
        sub = predictions[predictions.drug_id.eq(drug_id) & predictions.model.eq(arm)].copy()
        if sub.canonical_model_id.duplicated().any() or set(sub.canonical_model_id) - ids:
            raise ValueError("Duplicate/unexpected prediction ID in external arm")
        if not sub.status.isin((PREDICTED, ABSTAINED)).all():
            raise ValueError("Unknown external prediction status")
        numeric_prediction = pd.to_numeric(sub.predicted_auc, errors="coerce").to_numpy(dtype=float)
        is_predicted = sub.status.eq(PREDICTED).to_numpy()
        all_frozen_predicted[arm] = int(is_predicted.sum())
        if not np.isfinite(numeric_prediction[is_predicted]).all():
            raise ValueError("Predicted external row has nonfinite model output")
        if np.isfinite(numeric_prediction[~is_predicted]).any():
            raise ValueError("Abstained external row carries a numeric prediction")
        scored_sub = sub[sub.canonical_model_id.isin(response_ids)]
        joined_group = indexed.biological_group_id.reindex(scored_sub.canonical_model_id).to_numpy()
        if not np.array_equal(joined_group, scored_sub.biological_group_id.to_numpy()):
            raise ValueError("Prediction biological group differs from eligible roster")
        by_arm[arm] = sub.set_index("canonical_model_id").reindex(indexed.index)
        absent = by_arm[arm].status.isna()
        reasons = (by_arm[arm].abstain_reason if "abstain_reason" in by_arm[arm]
                   else pd.Series("abstained_unspecified", index=indexed.index))
        reasons = reasons.where(~absent, "missing_prediction_row").fillna("abstained_unspecified")
        reason_counts[arm] = dict(Counter(
            part for combined in reasons[~by_arm[arm].status.eq(PREDICTED)].tolist()
            for part in str(combined).split(";") if part))
    left, right = by_arm[arm_a], by_arm[arm_b]
    a_pred = pd.to_numeric(left.predicted_auc, errors="coerce").to_numpy(dtype=float)
    b_pred = pd.to_numeric(right.predicted_auc, errors="coerce").to_numpy(dtype=float)
    a_ok = left.status.eq(PREDICTED).fillna(False).to_numpy() & np.isfinite(a_pred)
    b_ok = right.status.eq(PREDICTED).fillna(False).to_numpy() & np.isfinite(b_pred)
    common = a_ok & b_ok
    paired = indexed.iloc[np.flatnonzero(common)].reset_index().copy()
    paired["prediction_a"] = a_pred[common]
    paired["prediction_b"] = b_pred[common]
    n_common = int(common.sum())
    return paired, {"drug_id": int(drug_id), "arm_a": arm_a, "arm_b": arm_b,
                    "n_frozen_eligible_ids": len(frozen_ids),
                    "frozen_eligible_roster_supplied": frozen_eligible_ids is not None,
                    "n_response_valid_ids": len(base),
                    "n_missing_or_invalid_response": len(frozen_ids) - len(base),
                    "n_response_valid_groups": base.biological_group_id.nunique(),
                    "n_arm_a_predicted_all_frozen": all_frozen_predicted[arm_a],
                    "n_arm_b_predicted_all_frozen": all_frozen_predicted[arm_b],
                    "n_arm_a_predicted": int(a_ok.sum()), "n_arm_b_predicted": int(b_ok.sum()),
                    "n_common_ids": n_common, "n_common_groups": paired.biological_group_id.nunique(),
                    "below_design_coverage": paired.biological_group_id.nunique() < design_min_groups,
                    "n_arm_a_abstain_or_missing": int((~a_ok).sum()),
                    "n_arm_b_abstain_or_missing": int((~b_ok).sum()),
                    "arm_a_noncoverage_reasons": reason_counts[arm_a],
                    "arm_b_noncoverage_reasons": reason_counts[arm_b]}


def paired_rho_difference(paired: pd.DataFrame, *, min_groups: int = 3) -> dict:
    needed = {"biological_group_id", "observed_auc", "prediction_a", "prediction_b"}
    if not needed.issubset(paired):
        raise ValueError("Paired table lacks required columns")
    g = paired.biological_group_id.tolist()
    rho_a, reason_a = spearman_checked(paired.observed_auc, paired.prediction_a, g, min_groups=min_groups)
    rho_b, reason_b = spearman_checked(paired.observed_auc, paired.prediction_b, g, min_groups=min_groups)
    return {"rho_a": rho_a, "rho_b": rho_b,
            "delta_rho": None if rho_a is None or rho_b is None else rho_a - rho_b,
            "undefined_reason": ";".join(x for x in (reason_a, reason_b) if x),
            "n_common_groups": len(set(g))}


def strict_fixed_eight_macro(per_drug: Mapping[int, float | None],
                             fixed_drug_ids: Sequence[int]) -> tuple[float | None, str]:
    ids = tuple(fixed_drug_ids)
    if len(ids) != 8 or len(set(ids)) != 8:
        raise ValueError("Strict macro requires exactly eight distinct frozen drug IDs")
    if set(per_drug) != set(ids):
        return None, "incomplete_fixed_eight_panel"
    values = [per_drug[drug] for drug in ids]
    if any(value is None or not np.isfinite(value) for value in values):
        return None, "at_least_one_drug_undefined"
    return float(np.mean(values)), ""


@dataclass(frozen=True)
class SharedBootstrap:
    groups: tuple[str, ...]
    counts: np.ndarray
    seed: int


def shared_bootstrap_sha256(shared: SharedBootstrap) -> str:
    """Digest group order, matrix dimensions and counts with fixed byte order."""
    h = hashlib.sha256()
    h.update(names_sha256(shared.groups).encode("ascii"))
    h.update(np.asarray(shared.counts.shape, dtype="<u8").tobytes())
    h.update(np.ascontiguousarray(shared.counts, dtype="<u8").tobytes())
    return h.hexdigest()


def make_shared_bootstrap(global_group_ids: Sequence[str], *, draws: int,
                          seed: int, expected_groups: int | None = None) -> SharedBootstrap:
    groups = _unique_names(global_group_ids, "Global biological group universe")
    if expected_groups is not None and len(groups) != expected_groups:
        raise ValueError("Frozen external biological-group count changed")
    if draws < 1:
        raise ValueError("Bootstrap draws must be positive")
    names, counts = generate_global_bootstrap_weights(groups, draws, seed)
    return SharedBootstrap(tuple(names.tolist()), counts, int(seed))


def _interval(point: float | None, draws: np.ndarray, *, min_valid_fraction: float,
              family_size: int | None = None) -> dict:
    valid = draws[np.isfinite(draws)]
    requested = len(draws)
    if not 0 < min_valid_fraction <= 1 or (family_size is not None and family_size < 1):
        raise ValueError("Invalid interval validity/family configuration")
    enough = len(valid) >= ceil(requested * min_valid_fraction)
    result = {"point": point, "bootstrap_requested": requested,
              "bootstrap_valid": len(valid), "bootstrap_invalid": requested - len(valid),
              "minimum_valid_required": ceil(requested * min_valid_fraction),
              "nominal_95_lower": None, "nominal_95_upper": None,
              "family_adjusted_lower": None, "family_adjusted_upper": None,
              "family_size": family_size,
              "ci_status": "reported" if point is not None and enough else
                  "withheld_undefined_point" if point is None else "withheld_insufficient_valid_draws",
              "interval_scope": "conditional_on_fixed_trained_models_and_external_biological_groups"}
    if point is not None and enough:
        result["nominal_95_lower"], result["nominal_95_upper"] = map(float, np.quantile(valid, [.025, .975], method="linear"))
        if family_size is not None:
            tail = .05 / (2 * family_size)
            result["family_adjusted_lower"], result["family_adjusted_upper"] = map(
                float, np.quantile(valid, [tail, 1 - tail], method="linear"))
    return result


def bootstrap_paired_rho(pairs_by_drug: Mapping[int, pd.DataFrame], *,
                         fixed_drug_ids: Sequence[int], shared: SharedBootstrap,
                         min_groups: int = 3, min_valid_fraction: float = .99,
                         family_size: int = 9, chunk_size: int = 256) -> dict:
    """Use one global group resampling schedule for all drugs and both arms.

    Exact duplicate-observation midranks are recomputed within each bootstrap
    draw. A strict macro draw is valid only if *every* fixed drug's paired
    difference is defined in that same draw; no available-case nanmean is used.
    """
    ids = tuple(fixed_drug_ids)
    if len(ids) != 8 or len(set(ids)) != 8 or set(pairs_by_drug) != set(ids):
        raise ValueError("Paired bootstrap requires the exact fixed-eight panel")
    if shared.counts.ndim != 2 or shared.counts.shape[1] != len(shared.groups):
        raise ValueError("Shared bootstrap count matrix has wrong shape")
    if (not shared.groups or len(set(shared.groups)) != len(shared.groups) or
            len(shared.counts) < 1 or not np.issubdtype(shared.counts.dtype, np.integer) or
            np.any(shared.counts < 0) or
            not np.all(shared.counts.sum(axis=1) == len(shared.groups))):
        raise ValueError("Invalid fixed-universe shared bootstrap multiplicities")
    if chunk_size < 1:
        raise ValueError("Bootstrap chunk size must be positive")
    group_index = {group: index for index, group in enumerate(shared.groups)}
    b = len(shared.counts)
    point = {}
    boot = {}
    for drug in ids:
        paired = pairs_by_drug[drug]
        if paired.biological_group_id.duplicated().any():
            raise ValueError("External paired set has repeated biological group")
        groups = paired.biological_group_id.tolist()
        if any(group not in group_index for group in groups):
            raise ValueError("Paired biological group absent from shared global universe")
        point[drug] = paired_rho_difference(paired, min_groups=min_groups)
        a = np.full(b, np.nan)
        r = np.full(b, np.nan)
        invalid_reason = np.full(b, "undefined_numeric", dtype="<U32")
        invalid_a = np.full(b, "undefined_numeric", dtype="<U32")
        invalid_r = np.full(b, "undefined_numeric", dtype="<U32")
        if groups:
            columns = [group_index[g] for g in groups]
            y = paired.observed_auc.to_numpy(dtype=np.float64)
            pa = paired.prediction_a.to_numpy(dtype=np.float64)
            pr = paired.prediction_b.to_numpy(dtype=np.float64)
            for start in range(0, b, chunk_size):
                stop = min(start + chunk_size, b)
                weights = shared.counts[start:stop, columns]
                group_n = np.count_nonzero(weights, axis=1)
                ma = weighted_metrics(y, pa, weights)["spearman"]
                mr = weighted_metrics(y, pr, weights)["spearman"]
                ma[group_n < min_groups] = np.nan
                mr[group_n < min_groups] = np.nan
                # Exact-value ranges among groups present in this draw identify
                # why a rank correlation is undefined without resampling again.
                present = weights > 0
                def constant(values):
                    low = np.min(np.where(present, values[None, :], np.inf), axis=1)
                    high = np.max(np.where(present, values[None, :], -np.inf), axis=1)
                    return low == high
                cy, ca, cr = constant(y), constant(pa), constant(pr)
                local_reason = invalid_reason[start:stop]
                local_reason[group_n < min_groups] = "insufficient_groups"
                enough = group_n >= min_groups
                local_reason[enough & cy] = "constant_response"
                local_reason[enough & ~cy & ca & cr] = "both_predictions_constant"
                local_reason[enough & ~cy & ca & ~cr] = "candidate_constant"
                local_reason[enough & ~cy & ~ca & cr] = "reference_constant"
                local_reason[np.isfinite(ma) & np.isfinite(mr)] = ""
                for local_arm, is_constant, metric in (
                        (invalid_a[start:stop], ca, ma),
                        (invalid_r[start:stop], cr, mr)):
                    local_arm[group_n < min_groups] = "insufficient_groups"
                    local_arm[enough & cy] = "constant_response"
                    local_arm[enough & ~cy & is_constant] = "constant_prediction"
                    local_arm[np.isfinite(metric)] = ""
                a[start:stop], r[start:stop] = ma, mr
        boot[drug] = (a, r, invalid_reason, invalid_a, invalid_r)
    per_drug = {}
    for drug in ids:
        a, r, reasons, reasons_a, reasons_r = boot[drug]
        current = point[drug]
        difference = a - r
        per_drug[drug] = {**current,
                          "invalid_paired_draw_reasons": dict(Counter(reasons[~np.isfinite(difference)].tolist())),
                          "invalid_candidate_rho_draw_reasons": dict(Counter(reasons_a[~np.isfinite(a)].tolist())),
                          "invalid_reference_rho_draw_reasons": dict(Counter(reasons_r[~np.isfinite(r)].tolist())),
                          "rho_a_interval": _interval(current["rho_a"], a, min_valid_fraction=min_valid_fraction),
                          "rho_b_interval": _interval(current["rho_b"], r, min_valid_fraction=min_valid_fraction),
                          "delta_interval": _interval(current["delta_rho"], difference,
                              min_valid_fraction=min_valid_fraction, family_size=family_size)}
    macro_point, macro_reason = strict_fixed_eight_macro({drug: point[drug]["delta_rho"] for drug in ids}, ids)
    all_differences = np.stack([boot[drug][0] - boot[drug][1] for drug in ids])
    valid_all = np.isfinite(all_differences).all(axis=0)
    macro_draws = np.full(b, np.nan)
    macro_draws[valid_all] = all_differences[:, valid_all].mean(axis=0)
    macro = {"delta_rho": macro_point, "undefined_reason": macro_reason,
             "invalid_draw_reasons": {"at_least_one_drug_undefined": int((~valid_all).sum())},
             "delta_interval": _interval(macro_point, macro_draws,
                 min_valid_fraction=min_valid_fraction, family_size=family_size)}
    return {"fixed_drug_ids": ids, "bootstrap_seed": shared.seed,
            "shared_global_groups": len(shared.groups),
            "shared_group_order": shared.groups,
            "shared_multiplicity_sha256": shared_bootstrap_sha256(shared),
            "per_drug": per_drug, "strict_macro": macro}
