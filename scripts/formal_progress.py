"""Compact read-only progress; no individual responses or performance values."""
import argparse
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--run-dir',default='runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1')
args=ap.parse_args()
run=Path(args.run_dir);run=run if run.is_absolute() else ROOT/run
meta=json.loads((run/'formal_metadata.json').read_text(encoding='utf-8'))
resource=sorted(run.glob('resource_samples_session*.csv'))[-1]
with resource.open('rb') as stream:
    stream.seek(max(0,resource.stat().st_size-1000))
    last=stream.read().decode('utf-8').strip().splitlines()[-1].split(',')
print(json.dumps({'status':meta['status'],'completed_tasks':meta['completed_tasks'],
    'expected_tasks':meta['expected_tasks'],'elapsed_session_seconds':float(last[0]),
    'process_tree_rss_GiB':float(last[1])/2**30,'available_ram_GiB':float(last[2])/2**30,
    'inflight':int(last[3]),'execution_commit':meta['execution_commit'][:7]},ensure_ascii=False))
