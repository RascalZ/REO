"""Finish the authorized formal development QA/statistics/report once fits end.

This local companion does not launch fits, change a scientific lock, access CTRP,
retry failures or publish anything. Every downstream command fails closed.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--run-dir',default='runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1')
args=ap.parse_args()
run=Path(args.run_dir);run=run.resolve() if run.is_absolute() else (ROOT/run).resolve()
run.relative_to(ROOT)
state_path=run/'postprocess_status.json'
def save(state):
    state['updated_utc']=datetime.now(timezone.utc).isoformat()
    temp=state_path.with_suffix('.json.tmp')
    temp.write_text(json.dumps(state,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    temp.replace(state_path)
state={'status':'checking_completed_fit_process','completed_steps':[],'external_response_accessed':False}
save(state)
try:
    # The sequential launcher calls this only after the training child exits.
    # No live metadata polling is permitted on this Windows filesystem.
    meta=json.loads((run/'formal_metadata.json').read_text(encoding='utf-8'))
    if meta['completed'] is not True:raise RuntimeError('All fits and a successful training-process exit are required')
    directory=ROOT/'results/formal_development'
    commands=[('full_run_validation',['scripts/validate_formal_run.py','--run-dir',str(run)]),
              ('frozen_paired_statistics',['src/analyze_formal_results.py','--run-dir',str(run)]),
              ('statistics_source_validation',['scripts/validate_formal_run.py','--run-dir',str(run),'--statistics-dir',str(directory)]),
              ('aggregate_report_and_figures',['src/report_formal_development.py'])]
    with (ROOT/'logs/formal_postprocess.log').open('ab') as log:
        for label,command in commands:
            state.update(status='running',current_step=label);save(state)
            subprocess.run([sys.executable,*command],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
            state['completed_steps'].append(label);save(state)
    state.update(status='completed',current_step=None);save(state)
except BaseException as error:
    state.update(status='failed',error_type=type(error).__name__,error=str(error));save(state)
    raise
