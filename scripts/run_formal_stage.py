"""Run the frozen numerical source, then QA/inference/report sequentially.

Live progress reads append-only logs, never the atomically replaced metadata.
This orchestration repair changes no scientific source, input or algorithm.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--resume',action='store_true')
ap.add_argument('--execution-commit',required=True)
args=ap.parse_args()
environment=os.environ.copy();environment['REO_EXECUTION_COMMIT']=args.execution_commit
log=ROOT/('logs/formal_v1_resume.console.log' if args.resume else 'logs/formal_v1.console.log')
mode='ab' if args.resume else 'xb'
command=[sys.executable,'src/run_formal_evaluation.py']+(['--resume'] if args.resume else [])
with log.open(mode) as output:
    subprocess.run(command,cwd=ROOT,env=environment,stdout=output,stderr=subprocess.STDOUT,check=True)
subprocess.run([sys.executable,'scripts/finish_formal_stage.py'],cwd=ROOT,env=environment,check=True)
