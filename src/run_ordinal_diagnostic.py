"""Frozen four-case ordinal diagnostic with nine full-pipeline null controls.

New version only. Reuses immutable v1 I/O and split helpers, never v1 outputs.
Monitor append-only stdout/resources; never read live metadata on Windows.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.metadata, json, os, platform, time, traceback, threading
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_info
import run_formal_evaluation as base
from run_development_pilot import ROOT, sha, execution_source_snapshot, verify_development_inputs, install_data_guard
from ordinal_diagnostic_models import evaluate_diagnostic_outer_fold

SOURCE_RUN = ROOT/'runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1'
OOF_FIELDS = base.OOF_FIELDS + ['permutation_id','selection_rule','representation']

def atomic_json(path, obj):
    """Bounded retry for transient Windows handles; same payload, no fit retry."""
    for attempt in range(21):
        try:
            base.atomic_json(path,obj)
            return
        except PermissionError:
            if attempt==20: raise
            time.sleep(.1)

def checkpoint_name(task):
    return f"p{task['permutation_id']}_d{task['drug_id']}_r{task['repeat']}_f{task['outer_fold']}_{task['model']}.json"

def task_list(config):
    tasks=[]
    for permutation in range(config['permutation']['count']+1):
        for drug in config['drug_ids']:
            for repeat in range(config['outer_repeats']):
                for fold in range(config['outer_folds']):
                    for model in config['models']:
                        if model=='tissue_null': rule,kind='none','tissue_null'
                        else:
                            rule=next(r for r in config['selection_rules'] if model.startswith(r+'_'))
                            kind=model[len(rule)+1:]
                        tasks.append(dict(drug_id=drug,repeat=repeat,outer_fold=fold,model=model,
                                          permutation_id=permutation,selection_rule=rule,representation=kind))
    return tasks

def make_permutations(responses, config):
    records=[];coverage=[]
    for drug in config['drug_ids']:
        data=responses[responses.gdsc_drug_id.astype(int).eq(drug)].sort_values('canonical_model_id').reset_index(drop=True)
        if len(data)!=data.biological_group_id.nunique() or data.canonical_model_id.duplicated().any():
            raise ValueError('Permutation requires one unique model per independent group')
        values=data.response_auc.to_numpy(dtype=np.float64)
        strata=[np.flatnonzero(data.tissue.eq(t).to_numpy()) for t in sorted(set(data.tissue))]
        for permutation in range(config['permutation']['count']+1):
            donor=np.arange(len(data))
            rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([config['permutation']['seed'],drug,permutation])))
            if permutation:
                for positions in strata: donor[positions]=rng.permutation(positions)
            for target,source in enumerate(donor):
                records.append(dict(drug_id=drug,permutation_id=permutation,
                    target_canonical_model_id=data.iloc[target].canonical_model_id,
                    donor_canonical_model_id=data.iloc[source].canonical_model_id,
                    tissue=data.iloc[target].tissue,assigned_auc_hex=float(values[source]).hex()))
            coverage.append(dict(drug_id=drug,drug_name=data.iloc[0].drug_name,permutation_id=permutation,
                n_groups=len(data),n_tissues=len(strata),singleton_groups=sum(len(s)==1 for s in strata),
                identity_donors=int(np.sum(donor==np.arange(len(data)))),
                unchanged_response_count=int(np.sum(values[donor]==values)),
                changed_response_fraction=float(np.mean(values[donor]!=values))))
    return pd.DataFrame(records),pd.DataFrame(coverage)

def worker_init(run_dir, config, freeze_sha):
    base.worker_init(run_dir,config,freeze_sha)
    manifest=pd.read_csv(Path(run_dir)/'permutation_manifest.csv',keep_default_na=False)
    base.WORKER['permutation_lookup']={(int(d),int(p)):dict(zip(b.target_canonical_model_id,b.assigned_auc_hex))
                                     for (d,p),b in manifest.groupby(['drug_id','permutation_id'])}

def assigned_y(data, task, lookup):
    mapping=lookup[(task['drug_id'],task['permutation_id'])]
    return np.asarray([float.fromhex(mapping[name]) for name in data.canonical_model_id],dtype=np.float64)

def execute_task(task):
    w=base.WORKER;started=time.perf_counter();process=psutil.Process();cpu0=process.cpu_times()
    stop=threading.Event();peaks=[]
    def sample():
        while not stop.is_set():
            peaks.append(process.memory_info().rss);stop.wait(w['config']['resources']['rss_sampling_interval_seconds'])
    sampler=threading.Thread(target=sample,daemon=True);sampler.start()
    path=w['run']/'checkpoints'/checkpoint_name(task);events=[]
    obj=dict(schema_version='1.0',task=task,freeze_sha256=w['freeze_sha'],status='failed',oof_rows=[],
             result=None,fit_events=events,error=None)
    try:
        data,x=base.get_drug(task['drug_id']);y=assigned_y(data,task,w['permutation_lookup'])
        train,test,inner=base.local_split_indices(data,w['folds'],task,w['config']['inner_folds'])
        obj.update(expected_test_ids=sorted(data.iloc[test].canonical_model_id.tolist()),
                   expected_training_expression_ids=sorted(x[train,-1].astype(int).tolist()),
                   inner_training_expression_ids=[sorted(x[train[it],-1].astype(int).tolist()) for it,iv in inner],
                   outer_test_expression_ids=sorted(x[test,-1].astype(int).tolist()),
                   assigned_response_sha256=hashlib.sha256(y.tobytes()).hexdigest())
        result=evaluate_diagnostic_outer_fold(x,y,data.biological_group_id.to_numpy(),train,test,inner,
            task['representation'],task['selection_rule'],w['config'],w['stable'],audit_sink=events.append)
        predictions=result.pop('predictions')
        if not np.isfinite(predictions).all(): raise ValueError('Nonfinite heldout prediction')
        for index,prediction in zip(test,predictions):
            row=data.iloc[index]
            obj['oof_rows'].append(dict(**task,drug_name=row.drug_name,canonical_model_id=row.canonical_model_id,
                biological_group_id=row.biological_group_id,tissue=row.tissue,
                observed_auc=float(y[index]),predicted_auc=float(prediction)))
        result.pop('fit_events',None);obj.update(result=result,status='completed')
    except Exception: obj['error']=traceback.format_exc()
    finally:
        stop.set();sampler.join();cpu=process.cpu_times();mem=process.memory_info()
        obj['resources']=dict(wall_seconds=time.perf_counter()-started,
            cpu_seconds=cpu.user+cpu.system-cpu0.user-cpu0.system,sampled_peak_rss_bytes=max(peaks+[mem.rss]),
            windows_lifetime_peak_working_set_bytes=getattr(mem,'peak_wset',None),pid=process.pid,
            memory_scope='worker including libraries, mapped TPM and current drug; lifetime peak may precede task',
            threadpools=threadpool_info())
        obj['data_guard']=dict(external_response_accessed=False,blocked_attempts=len(w['blocked']),accessed_path_count=len(w['access']))
        atomic_json(path,obj)
    return dict(task=task,status=obj['status'],path=str(path),wall_seconds=obj['resources']['wall_seconds'],
                fit_count=len(events),error=obj['error'])

def validate_checkpoint(path,task,freeze_sha):
    return base.validate_resume_checkpoint(path,task,freeze_sha)

def validate_prepared_checkpoint(path,task,response,folds,lookup,config):
    """Rebuild registered or crash-landed checkpoint identities from bound inputs."""
    obj=json.loads(Path(path).read_text(encoding='utf-8'))
    data=response[response.gdsc_drug_id.astype(int).eq(task['drug_id'])].sort_values('canonical_model_id').reset_index(drop=True)
    y=assigned_y(data,task,lookup)
    train,test,inner=base.local_split_indices(data,folds,task,config['inner_folds'])
    expression_ids=data.expression_row_index.astype(int).to_numpy()
    expected=dict(expected_test_ids=sorted(data.iloc[test].canonical_model_id.tolist()),
        expected_training_expression_ids=sorted(expression_ids[train].tolist()),
        outer_test_expression_ids=sorted(expression_ids[test].tolist()),
        inner_training_expression_ids=[sorted(expression_ids[train[it]].tolist()) for it,iv in inner],
        assigned_response_sha256=hashlib.sha256(y.tobytes()).hexdigest())
    if any(obj.get(key)!=value for key,value in expected.items()):raise ValueError('Checkpoint differs from prepared identity/label inputs')
    by_id={row['canonical_model_id']:row for row in obj['oof_rows']}
    for index in test:
        row=by_id[data.iloc[index].canonical_model_id]
        if row['observed_auc']!=float(y[index]) or row['biological_group_id']!=data.iloc[index].biological_group_id or row['tissue']!=data.iloc[index].tissue:
            raise ValueError('Checkpoint scoring response/identity mismatch')

def aggregate(out,tasks,freeze_sha):
    counts=dict(tasks=0,fit_events=0,oof_predictions=0,failed=0)
    with (out/'oof_predictions.csv').open('w',newline='',encoding='utf-8-sig') as fp, \
         (out/'fit_events.jsonl').open('w',encoding='utf-8') as fe, \
         (out/'outer_fit_audits.jsonl').open('w',encoding='utf-8') as fa:
        writer=csv.DictWriter(fp,fieldnames=OOF_FIELDS);writer.writeheader()
        for task in tasks:
            path=out/'checkpoints'/checkpoint_name(task)
            if not validate_checkpoint(path,task,freeze_sha):raise ValueError('Failed task cannot be omitted')
            obj=json.loads(path.read_text(encoding='utf-8'));counts['tasks']+=1
            writer.writerows(obj['oof_rows']);counts['oof_predictions']+=len(obj['oof_rows'])
            for event in obj['fit_events']:
                fe.write(json.dumps({**task,**event},ensure_ascii=False,allow_nan=False)+'\n');counts['fit_events']+=1
            fa.write(json.dumps(dict(task=task,**obj['result'],resources=obj['resources'],data_guard=obj['data_guard']),ensure_ascii=False,allow_nan=False)+'\n')
    return counts

def verify_freeze(path):
    lock=json.loads(path.read_text(encoding='utf-8'))
    for name,expected in lock['files'].items():
        if sha(ROOT/name)!=expected:raise ValueError('Frozen diagnostic source/input changed: '+name)
    for name,expected in lock['v1_files'].items():
        if sha(ROOT/name)!=expected:raise ValueError('Preserved formal v1 changed: '+name)
    revision,hashes=execution_source_snapshot([ROOT/name for name in lock['execution_files']])
    return lock,revision,hashes

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=ROOT/'config/ordinal_diagnostic_v1.json')
    ap.add_argument('--lock',type=Path,default=ROOT/'config/ordinal_diagnostic_v1.lock.json')
    ap.add_argument('--run-dir',type=Path);ap.add_argument('--resume',action='store_true');args=ap.parse_args()
    cfg_path=args.config.resolve();lock_path=args.lock.resolve();cfg=json.loads(cfg_path.read_text(encoding='utf-8'))
    access,blocked=install_data_guard();lock,revision,hashes=verify_freeze(lock_path);freeze_sha=sha(lock_path)
    if sha(cfg_path)!=lock['files'][cfg_path.relative_to(ROOT).as_posix()]:raise ValueError('Wrong config')
    verification=verify_development_inputs('v2');out=(args.run_dir or ROOT/'runs'/cfg['run_name']).resolve()
    if not out.is_relative_to(ROOT/'runs'):raise ValueError('Run outputs must remain in repository runs directory')
    if out.exists()!=args.resume:raise ValueError('Existing run requires explicit --resume; new run must not resume')
    out.mkdir(parents=True,exist_ok=True);(out/'checkpoints').mkdir(exist_ok=True)
    metadata_path=out/'ordinal_metadata.json';tasks=task_list(cfg)
    if args.resume:
        meta=json.loads(metadata_path.read_text(encoding='utf-8'))
        if meta['freeze_sha256']!=freeze_sha or meta['execution_commit']!=revision:raise ValueError('Resume version mismatch')
        for name,expected in meta['prepared_input_hashes'].items():
            if sha(out/name)!=expected:raise ValueError('Prepared input changed: '+name)
        for name,expected in meta['checkpoint_hashes'].items():
            if sha(out/'checkpoints'/name)!=expected:raise ValueError('Registered checkpoint changed')
    else:
        manifest=pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv',dtype=str,keep_default_na=False)
        if set(manifest.role)!={'development'}:raise ValueError('Non-development identity')
        identities=manifest.drop_duplicates('canonical_model_id').sort_values('canonical_model_id').reset_index(drop=True)
        if len(identities)!=identities.biological_group_id.nunique():raise ValueError('New group structure requires new protocol')
        folds,split_warnings=base.make_global_split_records(identities,cfg)
        folds.to_csv(out/'global_fold_manifest.csv',index=False,encoding='utf-8-sig')
        previous=json.loads((SOURCE_RUN/'formal_metadata.json').read_text(encoding='utf-8'))
        if not previous['completed'] or sha(SOURCE_RUN/'development_responses.csv')!=previous['prepared_input_hashes']['development_responses.csv']:
            raise ValueError('Formal-v1 prepared response provenance failed')
        response=pd.read_csv(SOURCE_RUN/'development_responses.csv',dtype={'gdsc_drug_id':str},keep_default_na=False,float_precision='round_trip')
        response=response[response.gdsc_drug_id.astype(int).isin(cfg['drug_ids'])].copy()
        if len(response)!=lock['contract']['drug_model_observations']:raise ValueError('Diagnostic response coverage changed')
        response.to_csv(out/'development_responses.csv',index=False,encoding='utf-8-sig')
        permutations,coverage=make_permutations(response,cfg)
        permutations.to_csv(out/'permutation_manifest.csv',index=False,encoding='utf-8-sig')
        coverage.to_csv(out/'permutation_coverage.csv',index=False)
        packages=['numpy','pandas','scikit-learn','scipy','psutil','threadpoolctl']
        prepared=['global_fold_manifest.csv','development_responses.csv','permutation_manifest.csv','permutation_coverage.csv']
        meta=dict(schema_version='1.0',run_name=cfg['run_name'],created_utc=datetime.now(timezone.utc).isoformat(),
            config_sha256=sha(cfg_path),freeze_sha256=freeze_sha,execution_commit=revision,source_hashes=hashes,
            configuration=cfg,development_input_verification=verification,cohort_sha256=verification['freeze_sha256'],
            python=platform.python_version(),packages={p:importlib.metadata.version(p) for p in packages},
            expected_tasks=lock['contract']['tasks'],expected_fit_events=lock['contract']['fits'],
            expected_oof_predictions=lock['contract']['oof_predictions'],completed_tasks=0,completed=False,status='prepared',
            external_response_accessed=False,global_split_warnings=split_warnings,
            prepared_input_hashes={name:sha(out/name) for name in prepared},checkpoint_hashes={},sessions=[],
            source_prepared_sha256=sha(SOURCE_RUN/'development_responses.csv'),source_formal_metadata_sha256=sha(SOURCE_RUN/'formal_metadata.json'))
        atomic_json(metadata_path,meta)
    response=pd.read_csv(out/'development_responses.csv',keep_default_na=False,float_precision='round_trip')
    folds=pd.read_csv(out/'global_fold_manifest.csv',keep_default_na=False)
    permutation_manifest=pd.read_csv(out/'permutation_manifest.csv',keep_default_na=False)
    lookup={(int(d),int(p)):dict(zip(b.target_canonical_model_id,b.assigned_auc_hex))
            for (d,p),b in permutation_manifest.groupby(['drug_id','permutation_id'])}
    pending=[]
    for task in tasks:
        path=out/'checkpoints'/checkpoint_name(task)
        if path.exists():
            if validate_checkpoint(path,task,freeze_sha):
                validate_prepared_checkpoint(path,task,response,folds,lookup,cfg)
                meta['checkpoint_hashes'][path.name]=sha(path);continue
            archive=out/'failed_attempts';archive.mkdir(exist_ok=True)
            destination=archive/(path.stem+'_'+sha(path)+'.json')
            if destination.exists():raise ValueError('Repeated failed checkpoint; preserve and inspect manually')
            os.replace(path,destination);meta['checkpoint_hashes'].pop(path.name,None)
        pending.append(task)
    session=dict(started_utc=datetime.now(timezone.utc).isoformat(),pending_at_start=len(pending),
                 completed_tasks_this_session=0,status='running',max_inflight_observed=0)
    meta['sessions'].append(session);meta.update(completed=False,status='running',completed_tasks=len(tasks)-len(pending))
    atomic_json(metadata_path,meta);started=time.perf_counter();peak=0;minimum=psutil.virtual_memory().available
    print(json.dumps(dict(event='ordinal_start',pending_tasks=len(pending),total_tasks=len(tasks),execution_commit=revision,
                          available_GiB=minimum/2**30)),flush=True)
    first_error=None
    try:
        with (out/f'resource_samples_session{len(meta["sessions"])}.csv').open('w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=['elapsed_seconds','process_tree_rss_bytes','available_ram_bytes','inflight']);writer.writeheader()
            with ProcessPoolExecutor(max_workers=cfg['resources']['parallel_jobs'],initializer=worker_init,initargs=(str(out),cfg,freeze_sha)) as pool:
                inflight={};cursor=0;last_sample=-10.;last_heartbeat=-60.
                while cursor<len(pending) or inflight:
                    while cursor<len(pending) and len(inflight)<cfg['resources']['parallel_jobs'] and first_error is None:
                        task=pending[cursor];cursor+=1;inflight[pool.submit(execute_task,task)]=task
                    session['max_inflight_observed']=max(session['max_inflight_observed'],len(inflight))
                    elapsed=time.perf_counter()-started
                    if elapsed-last_sample>=1:
                        rss=0
                        for process in [psutil.Process()]+psutil.Process().children(recursive=True):
                            try:rss+=process.memory_info().rss
                            except psutil.Error:pass
                        available=psutil.virtual_memory().available;peak=max(peak,rss);minimum=min(minimum,available)
                        writer.writerow(dict(elapsed_seconds=elapsed,process_tree_rss_bytes=rss,available_ram_bytes=available,inflight=len(inflight)))
                        handle.flush();last_sample=elapsed
                    done,_=wait(list(inflight),timeout=1,return_when=FIRST_COMPLETED)
                    for future in done:
                        task=inflight.pop(future)
                        try:result=future.result()
                        except Exception:
                            first_error=traceback.format_exc();result=dict(task=task,status='failed',error=first_error)
                        if result['status']=='completed':session['completed_tasks_this_session']+=1
                        else:first_error=result['error']
                        if 'path' in result:meta['checkpoint_hashes'][checkpoint_name(task)]=sha(Path(result['path']))
                        meta['completed_tasks']=len(tasks)-len(pending)+session['completed_tasks_this_session']
                        print(json.dumps(dict(event='ordinal_task_finished',**result),ensure_ascii=False),flush=True)
                        atomic_json(metadata_path,meta)
                    if first_error and not inflight:raise RuntimeError(first_error)
                    if elapsed-last_heartbeat>=60:
                        print(json.dumps(dict(event='heartbeat',completed_tasks=meta['completed_tasks'],elapsed_seconds=round(elapsed,1),inflight=len(inflight))),flush=True)
                        last_heartbeat=elapsed
        counts=aggregate(out,tasks,freeze_sha)
        if (counts['tasks'],counts['fit_events'],counts['oof_predictions'])!=(meta['expected_tasks'],meta['expected_fit_events'],meta['expected_oof_predictions']):
            raise ValueError('Full diagnostic counts differ: '+str(counts))
        meta.update(completed=True,status='completed',output_counts=counts,
            output_hashes={name:sha(out/name) for name in ['oof_predictions.csv','fit_events.jsonl','outer_fit_audits.jsonl','permutation_coverage.csv']},
            preserved_failed_attempts=len(list((out/'failed_attempts').glob('*.json'))))
        session['status']='completed'
    except BaseException:
        meta.update(completed=False,status='failed',failure=traceback.format_exc());session['status']='failed';raise
    finally:
        session.update(wall_seconds=time.perf_counter()-started,sampled_process_tree_peak_rss_bytes=peak,
            minimum_available_ram_bytes=minimum,finished_utc=datetime.now(timezone.utc).isoformat(),
            rss_scope='sum parent and child RSS; shared pages may be double-counted')
        meta['data_guard']=dict(external_response_accessed=False,blocked_attempts=len(blocked),accessed_path_count=len(access))
        atomic_json(metadata_path,meta)
    print(json.dumps(dict(event='ordinal_completed',counts=meta['output_counts'],wall_seconds=session['wall_seconds'])),flush=True)

if __name__=='__main__':main()
