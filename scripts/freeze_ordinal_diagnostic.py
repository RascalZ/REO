"""Create or verify an additive diagnostic freeze; preserve every v1 source byte."""
from __future__ import annotations
import argparse, hashlib, json, sys
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from run_development_pilot import sha

def build_contract():
    cfg=json.loads((ROOT/'config/ordinal_diagnostic_v1.json').read_text(encoding='utf-8'))
    assert cfg['drug_ids']==[1091,1930,1010,1199]
    assert (cfg['outer_folds'],cfg['outer_repeats'],cfg['inner_folds'])==(3,1,3)
    assert cfg['permutation']['count']==9 and cfg['permutation']['seed']==2026101301
    assert cfg['selection_rules']==['mad','response_association']
    assert cfg['representations']==['expression','rank','pair_full']
    assert cfg['models']==['tissue_null']+[rule+'_'+kind for rule in cfg['selection_rules'] for kind in cfg['representations']]
    assert cfg['feature']['n_genes']==200 and cfg['repeat_seeds']==[20261013]
    assert cfg['elastic_net']['alpha_ratios']==[1.,.3,.1,.03,.01] and cfg['elastic_net']['l1_ratios']==[.5,1.]
    assert cfg['resources']['parallel_jobs']==cfg['resources']['blas_threads']==1
    assert cfg['statistics']['bootstrap_draws']==10000 and cfg['external_response_access']=='forbidden'
    manifest=pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv',dtype=str,keep_default_na=False)
    assert set(manifest.role)=={'development'}
    selected=manifest[manifest.gdsc_drug_id.astype(int).isin(cfg['drug_ids'])]
    assert selected.groupby('gdsc_drug_id').size().tolist()==[620]*4
    groups=manifest.drop_duplicates('canonical_model_id')
    assert len(groups)==groups.biological_group_id.nunique()==626
    oldlock=json.loads((ROOT/'config/formal_analysis_v1.lock.json').read_text(encoding='utf-8'))
    for name,expected in oldlock['files'].items():
        if sha(ROOT/name)!=expected:raise ValueError('Historical v1 frozen input changed: '+name)
    prefix='runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1/'
    source=json.loads((ROOT/(prefix+'formal_metadata.json')).read_text(encoding='utf-8'))
    assert source['completed'] and source['status']=='completed'
    assert sha(ROOT/(prefix+'development_responses.csv'))==source['prepared_input_hashes']['development_responses.csv']
    execution=['config/ordinal_diagnostic_v1.json','docs/ORDINAL_DIAGNOSTIC_V1_PROTOCOL.md',
        'src/run_ordinal_diagnostic.py','src/ordinal_diagnostic_models.py','src/analyze_ordinal_diagnostic.py',
        'src/run_formal_evaluation.py','src/run_development_pilot.py','src/project_paths.py','src/analyze_formal_results.py',
        'scripts/freeze_ordinal_diagnostic.py',*[str(p.relative_to(ROOT)).replace('\\','/') for p in sorted((ROOT/'src/modeling').glob('*.py'))]]
    files={name:sha(ROOT/name) for name in execution}
    for name in ['config/formal_analysis_v1.lock.json','config/cohort_freeze_v2.json','config/cohort_freeze_v2.lock.json',
                 'results/cohort_v2/development_sample_manifest.csv','results/cohort_v2/main_drugs.csv',
                 prefix+'formal_metadata.json',prefix+'development_responses.csv']:
        files[name]=sha(ROOT/name)
    return cfg,execution,files,oldlock['files'],dict(tasks=840,fits=22440,oof_predictions=173600,
        drug_model_observations=2480,global_development_groups=626,observed_tasks=84,permutation_tasks=756,
        observed_fits=2244,permutation_fits=20196)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--lock',action='store_true');ap.add_argument('--verify-only',action='store_true');args=ap.parse_args()
    cfg,execution,files,old,contract=build_contract();path=ROOT/'config/ordinal_diagnostic_v1.lock.json'
    if args.verify_only or path.exists():
        previous=json.loads(path.read_text(encoding='utf-8'))
        assert previous['files']==files and previous['v1_files']==old and previous['contract']==contract
        print(json.dumps(dict(status='verified',lock_sha256=sha(path),contract=contract)));return
    if not args.lock:
        print(json.dumps(dict(status='ready_for_lock',contract=contract)));return
    obj=dict(schema_version='1.0',status='frozen_pre_diagnostic_fits',created_utc=datetime.now(timezone.utc).isoformat(),
        model_config='config/ordinal_diagnostic_v1.json',protocol='docs/ORDINAL_DIAGNOSTIC_V1_PROTOCOL.md',
        execution_files=execution,files=files,v1_files=old,contract=contract,
        historical_outcomes_seen=True,external_response_accessed=False,mutable_plan_and_reports_excluded=True)
    with path.open('x',encoding='utf-8',newline='\n') as handle:json.dump(obj,handle,ensure_ascii=False,indent=2);handle.write('\n')
    print(json.dumps(dict(status='locked',lock_sha256=sha(path),contract=contract)))

if __name__=='__main__':main()
