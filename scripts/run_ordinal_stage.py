"""Sequential diagnostic fit, independent audit, statistics and report launcher."""
from __future__ import annotations
import argparse,json,os,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--resume',action='store_true')
    ap.add_argument('--execution-commit',default='e9b576d');args=ap.parse_args()
    env=dict(os.environ,REO_EXECUTION_COMMIT=args.execution_commit,PYTHONUNBUFFERED='1',PYTHONIOENCODING='utf-8')
    # Apply the already frozen single-thread policy before numerical imports.
    # Runtime threadpool_limits alone can leave large Windows commit reservations.
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
                'BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[key]='1'
    logs=ROOT/'logs';logs.mkdir(exist_ok=True)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    steps=[('fit',['src/run_ordinal_diagnostic.py']+(['--resume'] if args.resume else [])),
           ('validate',['scripts/validate_ordinal_diagnostic.py']),
           ('statistics',['src/analyze_ordinal_diagnostic.py']),
           ('report',['src/report_ordinal_diagnostic.py'])]
    for name,arguments in steps:
        start=time.perf_counter();log=logs/f'ordinal_{stamp}_{name}.log'
        print(json.dumps(dict(event='stage_start',stage=name,log=log.relative_to(ROOT).as_posix())),flush=True)
        with log.open('wb') as handle:
            result=subprocess.run([sys.executable,'-u',*arguments],cwd=ROOT,env=env,stdout=handle,stderr=subprocess.STDOUT)
        print(json.dumps(dict(event='stage_finished',stage=name,exit_code=result.returncode,
                              wall_seconds=time.perf_counter()-start)),flush=True)
        if result.returncode:raise SystemExit(result.returncode)

if __name__=='__main__':main()
