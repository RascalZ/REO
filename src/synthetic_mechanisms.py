"""Deterministic, explicitly bounded synthetic mechanism cohorts for P5e."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import numpy as np


def truncated_standard_normal(rng: np.random.Generator, size, bound: float = 3.0) -> np.ndarray:
    """Draw standard normals truncated by rejection, with frozen draw order."""
    values = np.asarray(rng.standard_normal(size), dtype=np.float64)
    flat = values.reshape(-1)
    invalid = np.flatnonzero((flat < -bound) | (flat > bound))
    while invalid.size:
        flat[invalid] = rng.standard_normal(invalid.size)
        invalid = np.flatnonzero((flat < -bound) | (flat > bound))
    return values


def generate_cohort(config: dict, scenario: str, replicate: int) -> dict:
    if scenario not in config['scenario_indices']:
        raise ValueError(f'Unknown synthetic scenario: {scenario}')
    if not 0 <= replicate < int(config['replicates_per_scenario']):
        raise ValueError('Replicate index is outside the frozen range')
    scenario_index = int(config['scenario_indices'][scenario])
    seed = 2026101400 + 100 * scenario_index + int(replicate)
    rng = np.random.Generator(np.random.PCG64(seed))
    n = int(config['n_samples']); p = int(config['n_genes'])
    tissue = np.repeat(np.arange(int(config['n_tissues']), dtype=np.int16),
                       int(config['samples_per_tissue']))
    h = (tissue.astype(np.float64) - 3.5) / 3.5
    j = np.arange(p, dtype=np.float64)
    gene_offset = 3.0 + 0.008 * j
    load = np.linspace(-1.0, 1.0, p, dtype=np.float64)
    eta = None

    if scenario == 'tissue_only':
        epsilon = truncated_standard_normal(rng, (n, p))
        log2_tpm = gene_offset[None, :] + 0.20 * h[:, None] * load[None, :] + 0.35 * epsilon
        response = 0.50 + 0.28 * h
    elif scenario == 'weak_molecular':
        epsilon = truncated_standard_normal(rng, (n, p))
        eta = truncated_standard_normal(rng, n)
        z = np.mean(epsilon[:, :5], axis=1) * np.sqrt(5.0)
        log2_tpm = gene_offset[None, :] + 0.15 * h[:, None] * load[None, :] + 0.35 * epsilon
        response = 0.50 + 0.05 * h + 0.04 * z + 0.14 * eta
    elif scenario == 'common_sample_amplitude':
        amplitude = truncated_standard_normal(rng, n)
        eta = truncated_standard_normal(rng, n)
        log2_tpm = gene_offset[None, :] + 0.20 * amplitude[:, None]
        response = 0.50 + 0.04 * h + 0.12 * amplitude + 0.02 * eta
    elif scenario == 'predefined_pair_order':
        epsilon = truncated_standard_normal(rng, (n, p))
        eta = truncated_standard_normal(rng, n)
        log2_tpm = gene_offset[None, :] + epsilon
        log2_tpm[:, 0] = 3.5 + 0.35 * epsilon[:, 0]
        log2_tpm[:, 1] = 3.5 + 0.35 * epsilon[:, 1]
        response = (0.50 + 0.04 * h
                    + 0.15 * np.sign(log2_tpm[:, 0] - log2_tpm[:, 1])
                    + 0.03 * eta)
    else:  # pragma: no cover - guarded above
        raise AssertionError(scenario)

    response = np.clip(np.asarray(response, dtype=np.float64), 0.01, 0.99)
    if log2_tpm.shape != (n, p) or not np.isfinite(log2_tpm).all() or np.any(log2_tpm < 0.0):
        raise ValueError('Synthetic log2(TPM) must be finite, nonnegative and correctly shaped')
    tpm = np.exp2(log2_tpm) - 1.0
    if (tpm.shape != (n, p) or not np.isfinite(tpm).all() or np.any(tpm < 0.0)
            or response.shape != (n,) or not np.isfinite(response).all()
            or np.any(response < 0.01) or np.any(response > 0.99)):
        raise ValueError('Synthetic TPM or response violates the frozen input contract')

    checks = {
        'scenario': scenario,
        'scenario_index': scenario_index,
        'replicate': int(replicate),
        'seed': seed,
        'n_samples': int(n),
        'n_genes': int(p),
        'n_tissues': int(np.unique(tissue).size),
        'samples_per_tissue_min': int(np.min(np.bincount(tissue))),
        'samples_per_tissue_max': int(np.max(np.bincount(tissue))),
        'tpm_min': float(np.min(tpm)),
        'tpm_max': float(np.max(tpm)),
        'response_min': float(np.min(response)),
        'response_max': float(np.max(response)),
        'log2_tpm_min': float(np.min(log2_tpm)),
        'within_tissue_response_max_ptp': max(
            float(np.ptp(response[tissue == value])) for value in np.unique(tissue)),
        'common_amplitude_rank_invariant': None,
        'common_amplitude_pair_invariant': None,
    }
    if scenario == 'common_sample_amplitude':
        order = np.argsort(tpm, axis=1, kind='stable')
        checks['common_amplitude_rank_invariant'] = bool(np.all(order == order[0]))
        left, right = np.triu_indices(p, k=1)
        pair_sign = np.greater(tpm[:, left], tpm[:, right]).astype(np.int8)
        pair_sign -= np.less(tpm[:, left], tpm[:, right]).astype(np.int8)
        checks['common_amplitude_pair_invariant'] = bool(np.all(pair_sign == pair_sign[0]))
        if not checks['common_amplitude_rank_invariant'] or not checks['common_amplitude_pair_invariant']:
            raise ValueError('Common-amplitude DGP did not preserve its exact rank/pair invariant')
    if scenario == 'tissue_only' and checks['within_tissue_response_max_ptp'] != 0.0:
        raise ValueError('Tissue-only DGP must have exactly constant responses within each tissue')

    return {
        'tpm': np.ascontiguousarray(tpm, dtype=np.float64),
        'response': np.ascontiguousarray(response, dtype=np.float64),
        'tissue': np.ascontiguousarray(tissue, dtype=np.int16),
        'sample_id': np.arange(1, n + 1, dtype=np.int32),
        'checks': checks,
    }


def array_sha256(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(value.dtype.str.encode('ascii'))
    digest.update(str(value.shape).encode('ascii'))
    digest.update(value.tobytes(order='C'))
    return digest.hexdigest()


def atomic_savez(path: Path, **arrays) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('wb') as stream:
        np.savez_compressed(stream, **arrays)
    os.replace(temp, path)
