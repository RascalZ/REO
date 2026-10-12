"""Frozen, resumable eight-drug development evaluation. CTRP is never read."""
from __future__ import annotations
import argparse, csv, hashlib, importlib.metadata, json, os, platform, subprocess
import threading, time, traceback, warnings
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_limits, threadpool_info
from run_development_pilot import (ROOT, sha, execution_source_snapshot,
    verify_development_inputs, install_data_guard, read_selected_responses)
from project_paths import original_input
from modeling.nested_cv import make_grouped_splits, evaluate_outer_fold

WORKER = {}
OOF_FIELDS = ['drug_id','drug_name','model','repeat','outer_fold','canonical_model_id',
              'biological_group_id','tissue','observed_auc','predicted_auc']

def atomic_json(path, obj):
    path = Path(path); temp = path.with_suffix(path.suffix+'.tmp')
    if 'task' in obj and 'oof_rows' in obj and 'freeze_sha256' in obj:
        obj=dict(obj);obj.pop('payload_sha256',None)
        obj['payload_sha256']=hashlib.sha256(json.dumps(obj,sort_keys=True,ensure_ascii=False,
                          allow_nan=False,separators=(',',':')).encode('utf-8')).hexdigest()
    temp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    os.replace(temp,path)

def checkpoint_name(task):
    return f"d{task['drug_id']}_r{task['repeat']}_f{task['outer_fold']}_{task['model']}.json"

def validate_resume_checkpoint(path, expected_task, freeze_sha):
    obj = json.loads(Path(path).read_text(encoding='utf-8'))
    digest=obj.pop('payload_sha256',None)
    actual=hashlib.sha256(json.dumps(obj,sort_keys=True,ensure_ascii=False,
                        allow_nan=False,separators=(',',':')).encode('utf-8')).hexdigest()
    if digest!=actual:raise ValueError('Checkpoint payload checksum mismatch: '+str(path))
    if obj.get('task') != expected_task or obj.get('freeze_sha256') != freeze_sha:
        raise ValueError('Checkpoint task/freeze mismatch: '+str(path))
    if obj.get('status') != 'completed':
        return False
    rows = obj.get('oof_rows',[])
    ids = [r['canonical_model_id'] for r in rows]
    if not rows or len(ids)!=len(set(ids)) or sorted(ids)!=obj['expected_test_ids']:
        raise ValueError('Checkpoint prediction coverage invalid: '+str(path))
    for r in rows:
        if any(r[k] != expected_task[k] for k in expected_task):
            raise ValueError('Checkpoint row context invalid')
        if not np.isfinite([r['observed_auc'],r['predicted_auc']]).all():
            raise ValueError('Checkpoint contains nonfinite predictions')
    return True

def find_pending_tasks(task_list, checkpoint_dir, freeze_sha):
    pending = []
    for task in task_list:
        path=Path(checkpoint_dir)/checkpoint_name(task)
        if not path.exists() or not validate_resume_checkpoint(path,task,freeze_sha):
            pending.append(task)
    return pending

def task_list(config):
    return [{'drug_id':int(d),'repeat':r,'outer_fold':f,'model':m}
            for d in config['drug_ids'] for r in range(config['outer_repeats'])
            for f in range(config['outer_folds']) for m in config['models']]

def make_global_split_records(ids, config):
    rows=[]; split_warnings=[]
    for repeat,seed in enumerate(config['repeat_seeds']):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            outer=make_grouped_splits(ids.biological_group_id,ids.tissue,config['outer_folds'],seed)
            for f,(train,test) in enumerate(outer):
                for role,index in [('train',train),('test',test)]:
                    for i in index:
                        rows.append({**ids.iloc[i][['canonical_model_id','biological_group_id','tissue']].to_dict(),
                                     'repeat':repeat,'outer_fold':f,'inner_fold':-1,'role':role})
                inner=make_grouped_splits(ids.iloc[train].biological_group_id,ids.iloc[train].tissue,
                                         config['inner_folds'],seed+100+f)
                for j,(it,iv) in enumerate(inner):
                    for role,index in [('train',it),('validation',iv)]:
                        for i in train[index]:
                            rows.append({**ids.iloc[i][['canonical_model_id','biological_group_id','tissue']].to_dict(),
                                         'repeat':repeat,'outer_fold':f,'inner_fold':j,'role':role})
        split_warnings.extend({'repeat':repeat,'warning':str(w.message)} for w in caught)
    return pd.DataFrame(rows),split_warnings

def local_split_indices(data, folds, task, inner_folds):
    block=folds[(folds.repeat==task['repeat'])&(folds.outer_fold==task['outer_fold'])]
    def locate(inner,role,frame):
        names=set(block[(block.inner_fold==inner)&(block.role==role)].canonical_model_id)
        return np.flatnonzero(frame.canonical_model_id.isin(names).to_numpy())
    train=locate(-1,'train',data); test=locate(-1,'test',data)
    if len(train)+len(test)!=len(data):raise ValueError('Incomplete drug outer partition')
    outer=data.iloc[train]
    inner=[(locate(j,'train',outer),locate(j,'validation',outer)) for j in range(inner_folds)]
    return train,test,inner

def worker_init(run_dir,config,freeze_sha):
    WORKER.clear()
    access,blocked=install_data_guard()
    limits=threadpool_limits(limits=config['resources']['blas_threads'])
    genes=pd.read_csv(ROOT/'data_processed/expression_genes.tsv',sep='\t',dtype=str,keep_default_na=False)
    gene_idx=genes.index[genes.gene_type.eq('protein_coding')].to_numpy()
    stable=tuple(genes.iloc[gene_idx].ensembl_id)
    if len(stable)!=len(set(stable)):raise ValueError('Duplicate stable gene ID')
    responses=pd.read_csv(Path(run_dir)/'development_responses.csv',dtype={'gdsc_drug_id':str,
                         'canonical_model_id':str,'biological_group_id':str,'tissue':str},keep_default_na=False)
    tissue_names=sorted(set(responses.tissue))
    WORKER.update(run=Path(run_dir),config=config,freeze_sha=freeze_sha,access=access,blocked=blocked,
                  thread_limit=limits,gene_idx=gene_idx,stable=stable,responses=responses,
                  tissue_code={t:i for i,t in enumerate(tissue_names)},
                  folds=pd.read_csv(Path(run_dir)/'global_fold_manifest.csv',keep_default_na=False),
                  matrix=np.load(ROOT/'data_processed/expression_tpm.npy',mmap_mode='r'),drug=None)
    if WORKER['matrix'].dtype!=np.float64:raise ValueError('Raw TPM precision changed')

def get_drug(drug):
    if WORKER['drug']!=drug:
        WORKER.pop('x',None); WORKER.pop('data',None)
        data=WORKER['responses'][WORKER['responses'].gdsc_drug_id.eq(str(drug))].sort_values('canonical_model_id').reset_index(drop=True)
        indices=data.expression_row_index.astype(int).to_numpy(); genes=WORKER['gene_idx']
        x=np.empty((len(data),len(genes)+2),dtype=np.float64)
        x[:,:-2]=WORKER['matrix'][np.ix_(indices,genes)]
        x[:,-2]=data.tissue.map(WORKER['tissue_code']); x[:,-1]=indices
        if not np.isfinite(x).all() or (x[:,:-2]<0).any():raise ValueError('Unsupported TPM input')
        WORKER.update(drug=drug,data=data,x=x)
    return WORKER['data'],WORKER['x']

def execute_task(task):
    start=time.perf_counter(); process=psutil.Process(); cpu0=process.cpu_times()
    stop=threading.Event(); peaks=[]
    def sample():
        while not stop.is_set():
            peaks.append(process.memory_info().rss);stop.wait(.1)
    thread=threading.Thread(target=sample,daemon=True);thread.start()
    path=WORKER['run']/'checkpoints'/checkpoint_name(task)
    events=[]; obj={'schema_version':'1.0','task':task,'freeze_sha256':WORKER['freeze_sha'],
                  'status':'failed','oof_rows':[],'result':None,'fit_events':events,'error':None}
    try:
        data,x=get_drug(task['drug_id']); y=data.response_auc.to_numpy(dtype=np.float64)
        train,test,inner=local_split_indices(data,WORKER['folds'],task,WORKER['config']['inner_folds'])
        obj['expected_test_ids']=sorted(data.iloc[test].canonical_model_id.tolist())
        obj['expected_training_expression_ids']=sorted(x[train,-1].astype(int).tolist())
        obj['inner_training_expression_ids']=[sorted(x[train[it],-1].astype(int).tolist()) for it,iv in inner]
        obj['outer_test_expression_ids']=sorted(x[test,-1].astype(int).tolist())
        result=evaluate_outer_fold(x,y,data.biological_group_id.to_numpy(),train,test,inner,
                                   task['model'],WORKER['config'],WORKER['stable'],audit_sink=events.append)
        if not np.isfinite(result['predictions']).all():raise ValueError('Nonfinite outer prediction')
        for index,pred in zip(test,result.pop('predictions')):
            row=data.iloc[index]
            obj['oof_rows'].append({**task,'drug_name':row.drug_name,
                 'canonical_model_id':row.canonical_model_id,'biological_group_id':row.biological_group_id,
                 'tissue':row.tissue,'observed_auc':float(y[index]),'predicted_auc':float(pred)})
        result.pop('fit_events');obj['result']=result;obj['status']='completed'
    except Exception:
        obj['error']=traceback.format_exc()
    finally:
        stop.set();thread.join();cpu=process.cpu_times();mem=process.memory_info()
        obj['resources']={'wall_seconds':time.perf_counter()-start,
            'cpu_seconds':cpu.user+cpu.system-cpu0.user-cpu0.system,
            'sampled_peak_rss_bytes':max(peaks+[mem.rss]),
            'windows_lifetime_peak_working_set_bytes':getattr(mem,'peak_wset',None),
            'memory_scope':'worker RSS including libraries, mapped TPM pages and current drug cache; lifetime OS peak can precede this task',
            'threadpools':threadpool_info(),'pid':process.pid}
        obj['data_guard']={'external_response_accessed':False,'blocked_attempts':len(WORKER['blocked']),
                           'accessed_path_count':len(WORKER['access'])}
        atomic_json(path,obj)
    return {'task':task,'status':obj['status'],'path':str(path),'wall_seconds':obj['resources']['wall_seconds'],
            'fit_count':len(events),'error':obj['error']}

def verify_freeze(lock_path):
    lock=json.loads(lock_path.read_text(encoding='utf-8'))
    for name,expected in lock['files'].items():
        if sha(ROOT/name)!=expected:raise ValueError('Formal freeze file changed: '+name)
    revision,hashes=execution_source_snapshot([ROOT/name for name in lock['execution_files']])
    return lock,revision,hashes

def aggregate_checkpoints(out,tasks,freeze_sha):
    counts={'tasks':0,'fit_events':0,'oof_predictions':0,'failed':0}
    with (out/'oof_predictions.csv').open('w',newline='',encoding='utf-8-sig') as predfile, \
         (out/'fit_events.jsonl').open('w',encoding='utf-8') as eventfile, \
         (out/'outer_fit_audits.jsonl').open('w',encoding='utf-8') as auditfile:
        writer=csv.DictWriter(predfile,fieldnames=OOF_FIELDS);writer.writeheader()
        for task in tasks:
            path=out/'checkpoints'/checkpoint_name(task)
            if not path.exists():raise ValueError('Missing checkpoint '+str(task))
            obj=json.loads(path.read_text(encoding='utf-8'))
            if not validate_resume_checkpoint(path,task,freeze_sha):
                counts['failed']+=1;continue
            counts['tasks']+=1;writer.writerows(obj['oof_rows']);counts['oof_predictions']+=len(obj['oof_rows'])
            for e in obj['fit_events']:
                eventfile.write(json.dumps({**task,**e},ensure_ascii=False,allow_nan=False)+'\n')
                counts['fit_events']+=1
            auditfile.write(json.dumps({'task':task,**obj['result'],'resources':obj['resources']},ensure_ascii=False,allow_nan=False)+'\n')
    return counts

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=ROOT/'config/formal_model_v1.json')
    ap.add_argument('--lock',type=Path,default=ROOT/'config/formal_analysis_v1.lock.json')
    ap.add_argument('--run-dir',type=Path)
    ap.add_argument('--resume',action='store_true')
    args=ap.parse_args();args.config=args.config.resolve();args.lock=args.lock.resolve()
    config=json.loads(args.config.read_text(encoding='utf-8'));out=(args.run_dir or ROOT/'runs'/config['run_name']).resolve()
    access,blocked=install_data_guard();lock,revision,source_hashes=verify_freeze(args.lock)
    freeze_sha=sha(args.lock)
    if sha(args.config)!=lock['files'][args.config.relative_to(ROOT).as_posix()]:raise ValueError('Wrong formal config')
    verification=verify_development_inputs('v2');tasks=task_list(config)
    if out.exists() and not args.resume:raise ValueError('Run already exists; use explicit --resume')
    if args.resume and not out.exists():raise ValueError('Resume run missing')
    out.mkdir(parents=True,exist_ok=True);(out/'checkpoints').mkdir(exist_ok=True)
    meta_path=out/'formal_metadata.json'
    if args.resume:
        meta=json.loads(meta_path.read_text(encoding='utf-8'))
        if meta['freeze_sha256']!=freeze_sha or meta['execution_commit']!=revision:
            raise ValueError('Resume scientific version mismatch; pin REO_EXECUTION_COMMIT to recorded version')
        for name,expected in meta['prepared_input_hashes'].items():
            if sha(out/name)!=expected:raise ValueError('Prepared development run input changed: '+name)
        for name,expected in meta['checkpoint_hashes'].items():
            if sha(out/'checkpoints'/name)!=expected:raise ValueError('Registered checkpoint changed: '+name)
    else:
        manifest=pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv',dtype=str,keep_default_na=False)
        if set(manifest.role)!={'development'}:raise ValueError('Nondevelopment role in manifest')
        panel=pd.read_csv(ROOT/'results/cohort_v2/main_drugs.csv',dtype=str)
        if panel.gdsc_drug_id.astype(int).tolist()!=config['drug_ids']:raise ValueError('Fixed drug panel mismatch')
        ids=manifest.drop_duplicates('canonical_model_id').sort_values('canonical_model_id').reset_index(drop=True)
        if len(ids)!=ids.biological_group_id.nunique():raise ValueError('Multiple models per biological group require new protocol')
        folds,split_warnings=make_global_split_records(ids,config)
        folds.to_csv(out/'global_fold_manifest.csv',index=False,encoding='utf-8-sig')
        source=original_input('gdsc2');expected=pd.read_csv(ROOT/'data_manifest/original_source_sha256.csv',dtype=str)
        expected=expected[expected.path.eq(str(source))]
        if len(expected)!=1 or sha(source)!=expected.iloc[0].sha256:raise ValueError('Original GDSC input changed')
        read_selected_responses(manifest,source,out/'development_responses.csv')
        packages=['numpy','pandas','scikit-learn','scipy','psutil','threadpoolctl']
        meta={'schema_version':'1.0','run_name':config['run_name'],'created_utc':datetime.now(timezone.utc).isoformat(),
              'config_sha256':sha(args.config),'statistics_sha256':sha(ROOT/'config/formal_statistics_v1.json'),
              'protocol_sha256':sha(ROOT/'docs/FORMAL_ANALYSIS_PROTOCOL.md'),'freeze_sha256':freeze_sha,
              'execution_commit':revision,'source_hashes':source_hashes,'configuration':config,
              'development_input_verification':verification,'cohort_sha256':verification['freeze_sha256'],
              'python':platform.python_version(),'packages':{p:importlib.metadata.version(p) for p in packages},
              'expected_tasks':len(tasks),'expected_fit_events':18120,'expected_oof_predictions':104307,
              'expected_drug_observations':4967,'completed_tasks':0,'completed':False,'status':'prepared',
              'external_response_accessed':False,'global_split_warnings':split_warnings,
              'prepared_input_hashes':{n:sha(out/n) for n in ['global_fold_manifest.csv','development_responses.csv']},
              'checkpoint_hashes':{},'sessions':[]}
        atomic_json(meta_path,meta)
    pending=find_pending_tasks(tasks,out/'checkpoints',freeze_sha)
    # Recover only completed checkpoints that landed just before a parent crash.
    # Their identities/labels are independently rebuilt from bound prepared data.
    response=None;folds=None
    for task in tasks:
        name=checkpoint_name(task);path=out/'checkpoints'/name
        if path.exists() and name not in meta['checkpoint_hashes']:
            obj=json.loads(path.read_text(encoding='utf-8'))
            if obj['status']=='completed':
                if response is None:
                    response=pd.read_csv(out/'development_responses.csv',dtype={'gdsc_drug_id':str},keep_default_na=False)
                    folds=pd.read_csv(out/'global_fold_manifest.csv',keep_default_na=False)
                data=response[response.gdsc_drug_id.eq(str(task['drug_id']))].sort_values('canonical_model_id').reset_index(drop=True)
                train,test,inner=local_split_indices(data,folds,task,config['inner_folds'])
                if sorted(data.iloc[test].canonical_model_id.tolist())!=obj['expected_test_ids']:
                    raise ValueError('Unregistered checkpoint frozen test identities mismatch')
                if sorted(data.iloc[train].expression_row_index.astype(int).tolist())!=obj['expected_training_expression_ids']:
                    raise ValueError('Unregistered checkpoint frozen train identities mismatch')
                by_id={r['canonical_model_id']:r for r in obj['oof_rows']}
                for row in data.iloc[test].itertuples():
                    if by_id[row.canonical_model_id]['observed_auc']!=float(row.response_auc):
                        raise ValueError('Unregistered checkpoint development outcome mismatch')
            meta['checkpoint_hashes'][name]=sha(path)
    # Preserve a previous failed attempt before an explicitly authorized resume.
    # Completed tasks are never overwritten; failed attempts remain reviewable.
    for task in pending:
        previous=out/'checkpoints'/checkpoint_name(task)
        if previous.exists():
            archive=out/'failed_attempts';archive.mkdir(exist_ok=True)
            archived=archive/(previous.stem+'_'+sha(previous)+'.json')
            if archived.exists():
                if sha(archived)!=sha(previous):raise ValueError('Failed-attempt archive mismatch')
                previous.unlink()
            else:os.replace(previous,archived)
            meta['checkpoint_hashes'].pop(previous.name,None)
    session={'started_utc':datetime.now(timezone.utc).isoformat(),'pending_at_start':len(pending),
             'max_inflight_observed':0,'completed_tasks_this_session':0,'status':'running'}
    meta['sessions'].append(session);meta.update(status='running',completed=False)
    atomic_json(meta_path,meta);start=time.perf_counter();first_error=None
    samples=out/f"resource_samples_session{len(meta['sessions'])}.csv"
    print(json.dumps({'event':'formal_start','pending_tasks':len(pending),'max_workers':config['resources']['parallel_jobs'],
                      'available_GiB':psutil.virtual_memory().available/2**30,'execution_commit':revision}),flush=True)
    try:
        with samples.open('w',newline='',encoding='utf-8') as resource_file:
            writer=csv.DictWriter(resource_file,fieldnames=['elapsed_seconds','process_tree_rss_bytes','available_ram_bytes','inflight'])
            writer.writeheader();peak=0;min_available=psutil.virtual_memory().available
            with ProcessPoolExecutor(max_workers=config['resources']['parallel_jobs'],initializer=worker_init,
                                     initargs=(str(out),config,freeze_sha)) as pool:
                inflight={};next_task=0;last_sample=-10.;last_heartbeat=-60.
                while next_task<len(pending) or inflight:
                    available=psutil.virtual_memory().available
                    # The second worker is admitted only with a fixed RAM reserve;
                    # this changes dispatch, never seeds, data or model choices.
                    allowed=1 if available<config['resources']['second_worker_min_available_bytes'] else config['resources']['parallel_jobs']
                    while next_task<len(pending) and len(inflight)<allowed and first_error is None:
                        task=pending[next_task];next_task+=1;inflight[pool.submit(execute_task,task)]=task
                    session['max_inflight_observed']=max(session['max_inflight_observed'],len(inflight))
                    now=time.perf_counter()-start
                    if now-last_sample>=1:
                        processes=[psutil.Process()]+psutil.Process().children(recursive=True)
                        rss=0
                        for process in processes:
                            try:rss+=process.memory_info().rss
                            except psutil.Error:pass
                        peak=max(peak,rss);min_available=min(min_available,available)
                        writer.writerow({'elapsed_seconds':now,'process_tree_rss_bytes':rss,'available_ram_bytes':available,'inflight':len(inflight)})
                        resource_file.flush();last_sample=now
                    done,_=wait(list(inflight),timeout=1,return_when=FIRST_COMPLETED)
                    for future in done:
                        task=inflight.pop(future)
                        try:result=future.result()
                        except Exception:
                            first_error=traceback.format_exc();result={'task':task,'status':'failed','error':first_error}
                        if result['status']!='completed':first_error=result['error']
                        else:session['completed_tasks_this_session']+=1
                        if 'path' in result:meta['checkpoint_hashes'][checkpoint_name(task)]=sha(Path(result['path']))
                        print(json.dumps({'event':'formal_task_finished',**result,'error':None if not first_error else 'see preserved failure checkpoint'},ensure_ascii=False),flush=True)
                        meta['completed_tasks']=len(tasks)-len(pending)+session['completed_tasks_this_session']
                        atomic_json(meta_path,meta)
                    if first_error and not inflight:raise RuntimeError(first_error)
                    if now-last_heartbeat>=60:
                        print(json.dumps({'event':'heartbeat','completed_tasks':meta['completed_tasks'],'elapsed_seconds':round(now,1),
                                          'inflight':len(inflight),'available_GiB':round(available/2**30,2)}),flush=True)
                        last_heartbeat=now
        counts=aggregate_checkpoints(out,tasks,freeze_sha)
        if counts['tasks']!=len(tasks) or counts['fit_events']!=meta['expected_fit_events'] or counts['oof_predictions']!=meta['expected_oof_predictions']:
            raise ValueError('Formal output count mismatch: '+str(counts))
        meta.update(completed=True,status='completed',completed_tasks=counts['tasks'],output_counts=counts,
                    output_hashes={n:sha(out/n) for n in ['oof_predictions.csv','fit_events.jsonl','outer_fit_audits.jsonl']})
        meta['preserved_failed_attempts']=len(list((out/'failed_attempts').glob('*.json')))
        session['status']='completed'
    except BaseException:
        meta.update(completed=False,status='failed',failure=traceback.format_exc());session['status']='failed'
        raise
    finally:
        session.update(wall_seconds=time.perf_counter()-start,
                       sampled_process_tree_peak_rss_bytes=locals().get('peak'),
                       minimum_available_ram_bytes=locals().get('min_available'),
                       rss_scope='sum of parent and child RSS; shared mapped pages can be double-counted',
                       finished_utc=datetime.now(timezone.utc).isoformat())
        meta['data_guard']={'external_response_accessed':False,'blocked_attempts':len(blocked),'accessed_path_count':len(access)}
        atomic_json(meta_path,meta)
    print(json.dumps({'event':'formal_completed','counts':meta['output_counts'],'wall_seconds':session['wall_seconds']},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
