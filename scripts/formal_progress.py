"""Compact read-only progress; no individual responses or performance values."""
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--run-dir',default='runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1')
args=ap.parse_args()
run=Path(args.run_dir);run=run if run.is_absolute() else ROOT/run
# Do not open live formal_metadata.json: Windows can deny the writer's atomic
# replace while any reader holds that path, including FILE_SHARE_DELETE readers.
config=json.loads((ROOT/'config/formal_model_v1.json').read_text(encoding='utf-8'))
logs=[ROOT/'logs/formal_v1.console.log',ROOT/'logs/formal_v1_resume.console.log']
finished=set();base=0;status='running';commit='';session_starts=[]
for log in logs:
    if not log.exists():continue
    for line in log.read_text(encoding='utf-8',errors='replace').splitlines():
        if line.startswith('Traceback (most recent call last)'):status='failed'
        if not line.startswith('{'):continue
        try:event=json.loads(line)
        except json.JSONDecodeError:continue
        if event.get('event')=='formal_start':
            status='running';base=840-int(event['pending_tasks']);finished.clear()
            commit=event['execution_commit'][:7]
        elif event.get('event')=='formal_task_finished':
            if event.get('status')=='completed':
                task=event['task'];finished.add((task['drug_id'],task['repeat'],task['outer_fold'],task['model']))
            else:status='failed'
        elif event.get('event')=='formal_completed':status='completed'
resource=max(run.glob('resource_samples_session*.csv'),key=lambda p:int(p.stem.removeprefix('resource_samples_session')))
with resource.open('rb') as stream:
    stream.seek(max(0,resource.stat().st_size-1000))
    last=stream.read().decode('utf-8').strip().splitlines()[-1].split(',')
print(json.dumps({'status':status,'completed_tasks':base+len(finished),
    'expected_tasks':840,'elapsed_session_seconds':float(last[0]),
    'process_tree_rss_GiB':float(last[1])/2**30,'available_ram_GiB':float(last[2])/2**30,
    'inflight':int(last[3]),'execution_commit':commit,'live_metadata_opened':False},ensure_ascii=False))
