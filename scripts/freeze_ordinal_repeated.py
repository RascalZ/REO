"""Freeze additive repeated evaluation; independently preserve both historical locks."""
import argparse,json,sys
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from ordinal_repeated_runtime import sha


def build_contract():
    import pandas as pd
    cfg=json.loads((ROOT/'config/ordinal_repeated_v1.json').read_text(encoding='utf-8'))
    expected={'tasks':960,'fits':26160,'oof_predictions':119208,'drug_model_observations':4967,
        'global_development_groups':626,'molecular_outer_fits':840,'tissue_outer_fits':120,'observed_tasks':960,'permutation_tasks':0}
    assert cfg['expected_contract']==expected
    assert cfg['drug_ids']==[1931,1036,1191,1930,1199,1091,1010,1922]
    assert (cfg['outer_folds'],cfg['outer_repeats'],cfg['inner_folds'])==(5,3,3)
    assert cfg['repeat_seeds']==[20261013,20261014,20261015]
    assert cfg['models']==['tissue_null','expression_strong']+[r+'_'+k for r in ['mad','response_association'] for k in ['expression','rank','pair_full']]
    assert cfg['feature']['n_genes']==200 and cfg['feature']['strong_expression_genes']==2000
    assert cfg['permutation']['count']==0 and cfg['external_response_access']=='forbidden'
    assert cfg['statistics']['bootstrap_draws']==10000 and cfg['statistics']['bootstrap_seed']==2026101602
    assert cfg['statistics']['noninferiority']['margin'] is None and cfg['statistics']['noninferiority']['tests_enabled'] is False
    assert cfg['elastic_net']['alpha_ratios']==[1.,.3,.1,.03,.01] and cfg['elastic_net']['l1_ratios']==[.5,1.]
    assert cfg['resources']['parallel_jobs']==cfg['resources']['blas_threads']==1 and cfg['resources']['task_process_isolation'] is True
    manifest=pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv',dtype=str,keep_default_na=False)
    assert len(manifest)==4967 and set(manifest.role)=={'development'} and manifest.biological_group_id.nunique()==626
    assert manifest.canonical_model_id.nunique()==626
    old=json.loads((ROOT/'config/formal_analysis_v1.lock.json').read_text(encoding='utf-8'))
    diag=json.loads((ROOT/'config/ordinal_diagnostic_v1.lock.json').read_text(encoding='utf-8'))
    for preserved in (old['files'],diag['files'],diag['v1_files']):
        for name,expected_sha in preserved.items():
            if sha(ROOT/name)!=expected_sha:raise ValueError('Historical frozen file changed: '+name)
    execution=['config/ordinal_repeated_v1.json','docs/ORDINAL_REPEATED_V1_PROTOCOL.md',
        'src/run_ordinal_repeated.py','src/ordinal_repeated_models.py','src/ordinal_repeated_runtime.py',
        'src/analyze_ordinal_repeated.py','scripts/freeze_ordinal_repeated.py',
        'src/ordinal_diagnostic_models.py','src/run_formal_evaluation.py','src/run_development_pilot.py',
        'src/project_paths.py','src/analyze_formal_results.py']
    execution += [str(p.relative_to(ROOT)).replace('\\','/') for p in sorted((ROOT/'src/modeling').glob('*.py'))]
    files={name:sha(ROOT/name) for name in execution}
    for name in ['config/formal_analysis_v1.lock.json','config/ordinal_diagnostic_v1.lock.json',
        'config/cohort_freeze_v2.json','config/cohort_freeze_v2.lock.json',
        'results/cohort_v2/development_sample_manifest.csv','results/cohort_v2/main_drugs.csv',
        'runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1/formal_metadata.json',
        'runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1/development_responses.csv']:
        files[name]=sha(ROOT/name)
    return cfg,execution,files,old['files'],diag['files'],expected


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--lock',action='store_true');parser.add_argument('--verify-only',action='store_true');args=parser.parse_args()
    cfg,execution,files,old,diag,contract=build_contract();path=ROOT/'config/ordinal_repeated_v1.lock.json'
    if args.verify_only or path.exists():
        previous=json.loads(path.read_text(encoding='utf-8'))
        assert previous['files']==files and previous['v1_files']==old and previous['diagnostic_files']==diag and previous['contract']==contract
        print(json.dumps(dict(status='verified',lock_sha256=sha(path),contract=contract)));return
    if not args.lock:
        print(json.dumps(dict(status='ready_for_lock',contract=contract)));return
    obj=dict(schema_version='1.0',status='frozen_pre_repeated_fits',created_utc=datetime.now(timezone.utc).isoformat(),
        model_config='config/ordinal_repeated_v1.json',protocol='docs/ORDINAL_REPEATED_V1_PROTOCOL.md',execution_files=execution,
        files=files,v1_files=old,diagnostic_files=diag,contract=contract,historical_outcomes_seen=True,
        external_response_accessed=False,mutable_plan_reports_validators_excluded=True)
    with path.open('x',encoding='utf-8',newline='\n') as handle:json.dump(obj,handle,ensure_ascii=False,indent=2);handle.write('\n')
    print(json.dumps(dict(status='locked',lock_sha256=sha(path),contract=contract)))


if __name__=='__main__':main()
