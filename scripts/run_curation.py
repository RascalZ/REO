"""Run phase-1 preparation sequentially, or verify an existing lock without rewriting.

Complete reconstruction needs the original inputs described in DATA_AVAILABILITY.
This command never runs prediction models or examines external response magnitudes.
"""
from pathlib import Path
import argparse,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from project_paths import prepare_output_directories

def run(script,*args):
    subprocess.run([sys.executable,str(ROOT/'src'/script),*args],check=True,cwd=ROOT)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--download',action='store_true',help='Fetch documented reference data and missing CTRP members')
    ap.add_argument('--lock',action='store_true',help='Lock the cohort only after complete validation')
    args=ap.parse_args()
    if (ROOT/'config/cohort_freeze_v1.lock.json').exists():
        run('validate_cohorts.py');return
    prepare_output_directories()
    run('inventory_and_gdsc.py')
    run('audit_expression.py',*(['--download-references'] if args.download else []))
    if args.download:run('fetch_sanger_metadata.py')
    run('curate_ctrp.py',*(['--download'] if args.download else []))
    run('curate_drugs.py','--inventory-only')
    if args.download:run('curate_drugs.py','--pubchem-candidates')
    run('curate_drugs.py','--build-map')
    run('curate_drugs.py','--self-check')
    run('curate_cells.py');run('map_ctrp_cells.py')
    # Environment/register must exist before dependency snapshot; workbook uses finalized counts.
    import build_review_artifacts
    import csv,json,platform,importlib.metadata
    with (ROOT/'data_manifest/download_register.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f);w.writerow(['dataset','status','purpose','action_or_limit','path','source_url']);w.writerows(build_review_artifacts.DOWNLOADS)
    env={'python':platform.python_version(),'platform':platform.platform(),'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','openpyxl','rdkit']}}
    (ROOT/'config/curation_environment.json').write_text(json.dumps(env,indent=2),encoding='utf-8')
    run('freeze_cohorts.py')
    run('build_review_artifacts.py')
    run('validate_cohorts.py',*(['--lock'] if args.lock else []))

if __name__=='__main__':main()
