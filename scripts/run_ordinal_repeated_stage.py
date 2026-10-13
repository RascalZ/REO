"""Run frozen repeated evaluation, then post-exit QA, statistics and reporting.

No live run metadata is opened. Exit 75 means a preserved resource pause and
never authorizes downstream analysis; explicitly resume the same frozen run.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
THREAD_VARIABLES = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                    'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
RESOURCE_PAUSE_EXIT = 75


def child_environment(execution_commit):
    env = dict(os.environ, REO_EXECUTION_COMMIT=execution_commit,
               PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8')
    env.update({key: '1' for key in THREAD_VARIABLES})
    return env


def stage_commands(config, lock, run_dir, resume=False):
    common = ['--config', str(config), '--run-dir', str(run_dir)]
    return [
        ('fit', ['src/run_ordinal_repeated.py', *common, '--lock', str(lock)]
         + (['--resume'] if resume else [])),
        ('validate', ['scripts/validate_ordinal_repeated.py', *common, '--lock', str(lock)]),
        ('statistics', ['src/analyze_ordinal_repeated.py', *common]),
        ('report', ['src/report_ordinal_repeated.py', *common]),
    ]


def execute_stages(steps, execution_commit, logs=None):
    """Stop on any nonzero child exit; retain an exclusive log per stage."""
    logs = ROOT/'logs' if logs is None else Path(logs)
    logs.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    env = child_environment(execution_commit)
    for name, arguments in steps:
        started = time.perf_counter()
        log = logs/f'ordinal_repeated_{stamp}_{name}.log'
        print(json.dumps(dict(event='stage_start', stage=name, log=str(log))), flush=True)
        with log.open('xb') as handle:
            result = subprocess.run([sys.executable, '-u', *arguments], cwd=ROOT,
                                    env=env, stdout=handle, stderr=subprocess.STDOUT)
        event = dict(event='stage_finished', stage=name, exit_code=result.returncode,
                     wall_seconds=time.perf_counter()-started)
        if name == 'fit' and result.returncode == RESOURCE_PAUSE_EXIT:
            event.update(status='resource_paused', downstream_started=False,
                         action='Resume explicitly using the same scientific commit/config/lock; no task is omitted.')
        print(json.dumps(event), flush=True)
        if result.returncode:
            return result.returncode
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execution-commit', required=True)
    parser.add_argument('--config', type=Path, default=ROOT/'config/ordinal_repeated_v1.json')
    parser.add_argument('--lock', type=Path, default=ROOT/'config/ordinal_repeated_v1.lock.json')
    parser.add_argument('--run-dir', type=Path)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    config = args.config.resolve()
    # Configuration is immutable; never derive progress from live metadata.
    cfg = json.loads(config.read_text(encoding='utf-8'))
    run = args.run_dir.resolve() if args.run_dir else ROOT/'runs'/cfg['run_name']
    return execute_stages(stage_commands(config, args.lock.resolve(), run, args.resume),
                          args.execution_commit)


if __name__ == '__main__':
    sys.exit(main())
