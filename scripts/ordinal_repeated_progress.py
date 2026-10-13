"""Read append-only repeated-evaluation logs/resources; never live metadata."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def progress_from_log(path, expected_tasks):
    completed, failed, last_task, heartbeat, pressure = 0, 0, None, None, None
    status = 'running_or_interrupted'
    with Path(path).open(encoding='utf-8', errors='replace') as handle:
        for line in handle:
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            event = obj.get('event')
            if event == 'ordinal_repeated_start':
                completed = obj['total_tasks'] - obj['pending_tasks']
            elif event == 'ordinal_repeated_task_finished':
                completed += obj['status'] == 'completed'
                failed += obj['status'] == 'failed'
                last_task = {key: obj.get(key) for key in ('task', 'status', 'wall_seconds')}
                if obj['status'] == 'failed':
                    status = 'failed'
            elif event == 'heartbeat':
                heartbeat = obj
            elif event in ('resource_wait', 'resource_paused'):
                pressure = obj
                if event == 'resource_paused':
                    status = 'resource_paused'
            elif event == 'ordinal_repeated_completed':
                status = 'training_completed_pending_independent_audit'
    return dict(status=status, completed_tasks=completed, expected_tasks=expected_tasks,
                failed_tasks_this_session=failed, last_finished_task=last_task,
                last_heartbeat=heartbeat, last_resource_event=pressure,
                notice='Progress only. Safe pause is not model failure. Scientific completion requires post-exit independent validation.')


def last_complete_resource_row(path):
    last = None
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            if None not in row and None not in row.values():
                last = row
    return last


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'config/ordinal_repeated_v1.json')
    parser.add_argument('--run-dir', type=Path)
    parser.add_argument('--log', type=Path)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text(encoding='utf-8'))
    expected = len(cfg['drug_ids'])*cfg['outer_repeats']*cfg['outer_folds']*len(cfg['models'])
    logs = [args.log] if args.log else sorted((ROOT/'logs').glob('ordinal_repeated_*_fit.log'))
    if not logs:
        print(json.dumps(dict(status='no_fit_log', expected_tasks=expected)))
        return
    latest = logs[-1]
    output = dict(log=str(latest), **progress_from_log(latest, expected))
    run = args.run_dir if args.run_dir else ROOT/'runs'/cfg['run_name']
    resources = sorted(run.glob('resource_samples_session*.csv'),
                       key=lambda p: int(p.stem.split('session')[-1]))
    if resources:
        output['last_resource_sample'] = last_complete_resource_row(resources[-1])
    print(json.dumps(output, ensure_ascii=False))


if __name__ == '__main__':
    main()
