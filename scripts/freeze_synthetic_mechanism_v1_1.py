"""Freeze the bounded synthetic mechanism execution amendment before model fitting."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK_PATH = ROOT / 'config/synthetic_mechanism_execution_v1_1.lock.json'
INPUTS = [
    'docs/SYNTHETIC_MECHANISM_EXECUTION_V1_1.md',
    'config/synthetic_mechanisms_v1_1.json',
    'config/synthetic_mechanisms_v1.json',
    'config/p5e_protocol_v1.lock.json',
    'config/ordinal_repeated_v1.json',
    'config/ordinal_repeated_v1.lock.json',
    'src/run_synthetic_mechanisms.py',
    'src/analyze_synthetic_mechanisms.py',
    'src/synthetic_mechanisms.py',
    'src/ordinal_repeated_models.py',
    'src/ordinal_diagnostic_models.py',
    'src/modeling/nested_cv.py',
    'src/modeling/features.py',
    'src/modeling/reo_rf.py',
]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding='utf-8'))


def build() -> dict:
    cfg = read(ROOT / 'config/synthetic_mechanisms_v1_1.json')
    old_cfg = read(ROOT / 'config/synthetic_mechanisms_v1.json')
    p5e_lock = read(ROOT / 'config/p5e_protocol_v1.lock.json')
    model = read(ROOT / 'config/ordinal_repeated_v1.json')
    model_lock = read(ROOT / 'config/ordinal_repeated_v1.lock.json')
    if p5e_lock['inputs']['config/synthetic_mechanisms_v1.json'] != sha(ROOT / 'config/synthetic_mechanisms_v1.json'):
        raise ValueError('Original frozen synthetic config was modified')
    if model_lock['files']['config/ordinal_repeated_v1.json'] != sha(ROOT / 'config/ordinal_repeated_v1.json'):
        raise ValueError('Fixed-eight source model configuration was modified')
    if cfg['expected_contract'] != {
        'cohorts': 40, 'tasks': 960, 'fits': 26160,
        'per_sample_model_oof_predictions': 76800,
        'new_real_development_model_fits': 0, 'synthetic_model_fits': 26160,
    }:
        raise ValueError('Synthetic fit/task/OOF contract differs from frozen P5e v1')
    if old_cfg['expected_contract'] != cfg['expected_contract']:
        raise ValueError('Execution amendment must preserve the original P5e budget')
    if (cfg['scenario_order'] != old_cfg['scenario_order']
            or cfg['replicates_per_scenario'] != old_cfg['replicates_per_scenario']
            or cfg['n_samples'] != old_cfg['n_samples']
            or cfg['n_tissues'] != old_cfg['n_tissues']
            or cfg['samples_per_tissue'] != old_cfg['samples_per_tissue']
            or cfg['n_genes'] != old_cfg['n_genes']):
        raise ValueError('Execution amendment may not change the original synthetic cohort contract')
    if (model['elastic_net']['alpha_ratios'] != cfg['elastic_net']['alpha_ratios']
            or model['elastic_net']['l1_ratios'] != cfg['elastic_net']['l1_ratios']
            or model['feature']['n_genes'] != cfg['elastic_net']['feature_budget']):
        raise ValueError('Execution amendment differs from the pre-existing model/feature grid')
    inputs = {name: sha(ROOT / name) for name in INPUTS}
    return {
        'schema_version': '1.0',
        'status': 'locked_before_synthetic_fits',
        'protocol_version': 'P5e-synthetic-v1.1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'upstream_p5e_lock_sha256': sha(ROOT / 'config/p5e_protocol_v1.lock.json'),
        'upstream_model_lock_sha256': sha(ROOT / 'config/ordinal_repeated_v1.lock.json'),
        'original_synthetic_config_sha256': sha(ROOT / 'config/synthetic_mechanisms_v1.json'),
        'amendment_reason': 'The original DGP allowed negative raw TPM under rare Gaussian tails and left fold-seed scope ambiguous; v1.1 freezes rejection-truncated normal draws and outer/inner seed formulas without changing the task budget or substantive scenarios.',
        'inputs': inputs,
        'execution_contract': cfg['expected_contract'],
        'scenario_order': cfg['scenario_order'],
        'rng': cfg['rng'],
        'noise_support': cfg['scenario_noise_support'],
        'outer_split': cfg['outer_split'],
        'inner_split': cfg['inner_split'],
        'external_response_access': 'forbidden',
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lock', action='store_true')
    args = parser.parse_args()
    current = build()
    if LOCK_PATH.exists():
        previous = read(LOCK_PATH)
        for key, value in current.items():
            if key != 'created_utc' and previous.get(key) != value:
                raise ValueError('Frozen synthetic execution lock mismatch: ' + key)
        print(json.dumps({'status': 'verified', 'lock_sha256': sha(LOCK_PATH)}, ensure_ascii=False))
        return
    if not args.lock:
        print(json.dumps({'status': 'ready_for_lock', 'execution_contract': current['execution_contract']}, ensure_ascii=False))
        return
    LOCK_PATH.write_text(json.dumps(current, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'locked', 'lock_sha256': sha(LOCK_PATH)}, ensure_ascii=False))


if __name__ == '__main__':
    main()
