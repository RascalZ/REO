"""Read append-only diagnostic logs; never open live run metadata."""
from __future__ import annotations
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    cfg = json.loads((ROOT / 'config/ordinal_diagnostic_v1.json').read_text(encoding='utf-8'))
    logs = sorted((ROOT / 'logs').glob('ordinal_*_fit.log'))
    if not logs:
        print(json.dumps({'status': 'no_fit_log'}))
        return
    latest = logs[-1]
    completed = 0
    last_task = None
    failed = 0
    heartbeat = None
    for line in latest.read_text(encoding='utf-8', errors='replace').splitlines():
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            continue
        if obj.get('event') == 'ordinal_start':
            completed = obj['total_tasks'] - obj['pending_tasks']
        elif obj.get('event') == 'ordinal_task_finished':
            completed += obj['status'] == 'completed'
            failed += obj['status'] == 'failed'
            last_task = {'task': obj['task'], 'status': obj['status'],
                         'wall_seconds': obj['wall_seconds']}
        elif obj.get('event') == 'heartbeat':
            heartbeat = obj
    output = {'log': latest.relative_to(ROOT).as_posix(), 'completed_tasks': completed,
              'expected_tasks': 840, 'failed_tasks_this_session': failed,
              'last_finished_task': last_task, 'last_heartbeat': heartbeat,
              'notice': 'Progress only; completeness requires post-exit independent validation.'}
    resources = sorted((ROOT / 'runs' / cfg['run_name']).glob('resource_samples_session*.csv'),
                       key=lambda p: int(p.stem.split('session')[-1]))
    if resources:
        with resources[-1].open(encoding='utf-8', newline='') as handle:
            rows = csv.DictReader(handle)
            last_row = None
            for row in rows:
                if None not in row.values():
                    last_row = row
        output['last_resource_sample'] = last_row
    print(json.dumps(output, ensure_ascii=False))


if __name__ == '__main__':
    main()
