"""Frozen fixed-eight evaluation with one fresh child process per outer task."""
from __future__ import annotations
from ordinal_repeated_runtime import set_preimport_threads
EARLY_THREADS = set_preimport_threads()
import argparse, csv, hashlib, importlib.metadata, json, os, platform, subprocess, sys
import threading, time, traceback
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_info
import run_formal_evaluation as base
from run_development_pilot import ROOT, execution_source_snapshot, verify_development_inputs, install_data_guard
from ordinal_repeated_models import model_spec, evaluate_repeated_outer_fold
from ordinal_repeated_runtime import atomic_json, sha, system_memory, process_tree_memory, pressure_reasons

SOURCE_RUN = ROOT/'runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1'
OOF_FIELDS = base.OOF_FIELDS + ['permutation_id', 'selection_rule', 'representation', 'n_genes']


def task_list(config):
    return [dict(drug_id=int(drug), repeat=repeat, outer_fold=fold, model=model,
                 permutation_id=0, **model_spec(model, config))
            for drug in config['drug_ids'] for repeat in range(config['outer_repeats'])
            for fold in range(config['outer_folds']) for model in config['models']]


def checkpoint_name(task):
    return f"d{task['drug_id']}_r{task['repeat']}_f{task['outer_fold']}_{task['model']}.json"


def validate_checkpoint(path, task, freeze_sha):
    return base.validate_resume_checkpoint(path, task, freeze_sha)


def execute_task(task, run, config, freeze_sha, source_context=None):
    started = time.perf_counter()
    process = psutil.Process(); cpu0 = process.cpu_times()
    events, samples, access, blocked = [], [], [], []
    stop = threading.Event()
    obj = dict(schema_version='1.0', task=task, freeze_sha256=freeze_sha,
               status='failed', oof_rows=[], result=None, fit_events=events, error=None)
    if source_context is not None:
        obj.update(execution_commit=source_context['execution_commit'], source_hashes=source_context['source_hashes'])
    initialization_seconds = fit_work_seconds = None
    def sample():
        while not stop.is_set():
            try:
                m = process.memory_info()
                samples.append(dict(rss=int(m.rss), private=int(getattr(m, 'private', m.vms)), **system_memory()))
            except Exception:
                pass
            stop.wait(config['resources']['rss_sampling_interval_seconds'])
    sampler = threading.Thread(target=sample, daemon=True); sampler.start()
    try:
        base.worker_init(str(run), config, freeze_sha)
        w = base.WORKER; access, blocked = w['access'], w['blocked']
        # The immutable helper's default CSV parser can round a last bit.
        # Match prepared-source/controller/independent exact-float verification.
        w['responses'] = pd.read_csv(run/'development_responses.csv',
            dtype={'gdsc_drug_id':str,'canonical_model_id':str,'biological_group_id':str,'tissue':str},
            keep_default_na=False, float_precision='round_trip')
        data, x = base.get_drug(task['drug_id']); y = data.response_auc.to_numpy(dtype=np.float64)
        train, test, inner = base.local_split_indices(data, w['folds'], task, config['inner_folds'])
        initialization_seconds = time.perf_counter()-started
        obj.update(expected_test_ids=sorted(data.iloc[test].canonical_model_id.tolist()),
                   expected_training_expression_ids=sorted(x[train, -1].astype(int).tolist()),
                   inner_training_expression_ids=[sorted(x[train[it], -1].astype(int).tolist()) for it, iv in inner],
                   outer_test_expression_ids=sorted(x[test, -1].astype(int).tolist()),
                   assigned_response_sha256=hashlib.sha256(y.tobytes()).hexdigest())
        fit_start = time.perf_counter()
        result = evaluate_repeated_outer_fold(x, y, data.biological_group_id.to_numpy(), train, test,
                    inner, task['model'], config, w['stable'], audit_sink=events.append)
        fit_work_seconds = time.perf_counter()-fit_start
        prediction = result.pop('predictions')
        if not np.isfinite(prediction).all():
            raise ValueError('Nonfinite heldout prediction')
        for index, value in zip(test, prediction):
            row = data.iloc[index]
            obj['oof_rows'].append(dict(**task, drug_name=row.drug_name,
                canonical_model_id=row.canonical_model_id, biological_group_id=row.biological_group_id,
                tissue=row.tissue, observed_auc=float(y[index]), predicted_auc=float(value)))
        result.pop('fit_events', None); obj.update(status='completed', result=result)
    except BaseException:
        obj['error'] = traceback.format_exc()
    finally:
        stop.set(); sampler.join(); cpu = process.cpu_times(); mem = process.memory_info()
        if not samples:
            samples.append(dict(rss=int(mem.rss), private=int(getattr(mem, 'private', mem.vms)), **system_memory()))
        obj['resources'] = dict(wall_seconds=time.perf_counter()-started,
            initialization_wall_seconds=initialization_seconds, fit_work_wall_seconds=fit_work_seconds,
            cpu_seconds=cpu.user+cpu.system-cpu0.user-cpu0.system,
            sampled_peak_rss_bytes=max(s['rss'] for s in samples),
            sampled_peak_private_bytes=max(s['private'] for s in samples),
            minimum_available_ram_bytes=min(s['available_ram_bytes'] for s in samples),
            minimum_commit_free_bytes=min(s['system_commit_free_bytes'] for s in samples),
            windows_lifetime_peak_working_set_bytes=getattr(mem, 'peak_wset', None),
            pid=process.pid, process_create_time=process.create_time(), controller_pid=process.ppid(),
            task_process_isolated=True, preimport_threads_applied=EARLY_THREADS,
            memory_scope='Fresh task worker including imports, TPM mapping pages and preparation; sampled peaks not true upper bounds',
            threadpools=threadpool_info())
        obj['data_guard'] = dict(external_response_accessed=False, blocked_attempts=len(blocked), accessed_path_count=len(access))
        atomic_json(run/'checkpoints'/checkpoint_name(task), obj)
    return obj['status'] == 'completed'


def aggregate(out, tasks, freeze_sha):
    counts = dict(tasks=0, fit_events=0, oof_predictions=0, failed=0)
    with (out/'oof_predictions.csv').open('w', newline='', encoding='utf-8-sig') as fp, \
         (out/'fit_events.jsonl').open('w', encoding='utf-8', newline='\n') as fe, \
         (out/'outer_fit_audits.jsonl').open('w', encoding='utf-8', newline='\n') as fa:
        writer = csv.DictWriter(fp, fieldnames=OOF_FIELDS); writer.writeheader()
        for task in tasks:
            path = out/'checkpoints'/checkpoint_name(task)
            if not validate_checkpoint(path, task, freeze_sha):
                raise ValueError('Failed task cannot be omitted')
            obj = json.loads(path.read_text(encoding='utf-8')); counts['tasks'] += 1
            writer.writerows(obj['oof_rows']); counts['oof_predictions'] += len(obj['oof_rows'])
            for event in obj['fit_events']:
                fe.write(json.dumps({**task, **event}, ensure_ascii=False, allow_nan=False)+'\n')
                counts['fit_events'] += 1
            fa.write(json.dumps(dict(task=task, **obj['result'], resources=obj['resources'],
                data_guard=obj['data_guard']), ensure_ascii=False, allow_nan=False)+'\n')
    return counts


def verify_freeze(path):
    lock = json.loads(path.read_text(encoding='utf-8'))
    for category in ('files', 'v1_files', 'diagnostic_files'):
        for name, expected in lock[category].items():
            if sha(ROOT/name) != expected:
                raise ValueError('Frozen or historical scientific file changed: '+name)
    revision, hashes = execution_source_snapshot([ROOT/name for name in lock['execution_files']])
    return lock, revision, hashes


def verify_prepared_checkpoint(path, task, response, folds, config):
    obj = json.loads(path.read_text(encoding='utf-8'))
    data = response[response.gdsc_drug_id.astype(int).eq(task['drug_id'])].sort_values('canonical_model_id').reset_index(drop=True)
    y = data.response_auc.to_numpy(dtype=np.float64)
    train, test, inner = base.local_split_indices(data, folds, task, config['inner_folds'])
    indices = data.expression_row_index.astype(int).to_numpy()
    expected = dict(expected_test_ids=sorted(data.iloc[test].canonical_model_id.tolist()),
        expected_training_expression_ids=sorted(indices[train].tolist()), outer_test_expression_ids=sorted(indices[test].tolist()),
        inner_training_expression_ids=[sorted(indices[train[it]].tolist()) for it, iv in inner],
        assigned_response_sha256=hashlib.sha256(y.tobytes()).hexdigest())
    if any(obj.get(k) != v for k, v in expected.items()):
        raise ValueError('Checkpoint differs from bound prepared identity/label inputs')
    records = {r['canonical_model_id']: r for r in obj['oof_rows']}
    for index in test:
        source = data.iloc[index]; row = records[source.canonical_model_id]
        if row['observed_auc'] != float(y[index]) or row['tissue'] != source.tissue or row['biological_group_id'] != source.biological_group_id:
            raise ValueError('Checkpoint outcome/identity altered')


def recovery_snapshot(out, meta):
    records = {p.name: dict(status=json.loads(p.read_text(encoding='utf-8'))['status'], sha256=sha(p))
               for p in sorted((out/'checkpoints').glob('*.json'))}
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    atomic_json(out/f'recovery_before_{stamp}.json', dict(metadata=meta, checkpoints=records))


def child_command(out, task_index, context_sha):
    return [sys.executable, '-u', str(Path(__file__).resolve()), '--run-dir', str(out),
            '--worker-task', str(task_index), '--context-sha', context_sha]


def load_worker_context(out, expected_sha):
    """Bind each fresh interpreter to small immutable scientific/prepared files."""
    context_path = Path(out)/'worker_context.json'
    if not expected_sha or sha(context_path) != expected_sha:
        raise ValueError('Immutable worker context binding failed')
    context = json.loads(context_path.read_text(encoding='utf-8'))
    for name, expected in context['prepared_input_hashes'].items():
        if sha(Path(out)/name) != expected:
            raise ValueError('Worker prepared input changed: '+name)
    for name, expected in context['source_hashes'].items():
        if sha(ROOT/name) != expected:
            raise ValueError('Worker scientific source changed: '+name)
    return context


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'config/ordinal_repeated_v1.json')
    parser.add_argument('--lock', type=Path, default=ROOT/'config/ordinal_repeated_v1.lock.json')
    parser.add_argument('--run-dir', type=Path); parser.add_argument('--resume', action='store_true')
    parser.add_argument('--worker-task', type=int); parser.add_argument('--context-sha')
    args = parser.parse_args()
    if args.worker_task is not None:
        out = args.run_dir.resolve(); out.relative_to(ROOT/'runs')
        context = load_worker_context(out, args.context_sha)
        tasks = task_list(context['configuration'])
        return 0 if execute_task(tasks[args.worker_task], out, context['configuration'], context['freeze_sha256'], context) else 1
    config_path = args.config.resolve(); lock_path = args.lock.resolve()
    config = json.loads(config_path.read_text(encoding='utf-8'))
    access, blocked = install_data_guard(); lock, revision, hashes = verify_freeze(lock_path); freeze_sha = sha(lock_path)
    if sha(config_path) != lock['files'][config_path.relative_to(ROOT).as_posix()]:
        raise ValueError('Wrong frozen configuration')
    verification = verify_development_inputs('v2')
    out = (args.run_dir or ROOT/'runs'/config['run_name']).resolve(); out.relative_to(ROOT/'runs')
    if out.exists() != args.resume:
        raise ValueError('Existing run requires explicit --resume; new run must not resume')
    out.mkdir(parents=True, exist_ok=True); (out/'checkpoints').mkdir(exist_ok=True)
    metadata_path = out/'ordinal_repeated_metadata.json'; tasks = task_list(config)
    if args.resume:
        meta = json.loads(metadata_path.read_text(encoding='utf-8'))
        if meta['status'] == 'completed':
            raise ValueError('Completed research output is preserved; do not resume it')
        if meta['freeze_sha256'] != freeze_sha or meta['execution_commit'] != revision:
            raise ValueError('Resume scientific version mismatch')
        for name, expected in meta['prepared_input_hashes'].items():
            if sha(out/name) != expected:
                raise ValueError('Prepared input changed: '+name)
        for name, expected in meta['checkpoint_hashes'].items():
            if sha(out/'checkpoints'/name) != expected:
                raise ValueError('Registered checkpoint changed: '+name)
        recovery_snapshot(out, meta)
    else:
        manifest = pd.read_csv(ROOT/'results/cohort_v2/development_sample_manifest.csv', dtype=str, keep_default_na=False)
        if set(manifest.role) != {'development'}:
            raise ValueError('Nondevelopment input')
        identities = manifest.drop_duplicates('canonical_model_id').sort_values('canonical_model_id').reset_index(drop=True)
        if not (len(identities) == identities.biological_group_id.nunique() == 626):
            raise ValueError('Independent group structure changed')
        folds, split_warnings = base.make_global_split_records(identities, config)
        folds.to_csv(out/'global_fold_manifest.csv', index=False, encoding='utf-8-sig', lineterminator='\n')
        previous = json.loads((SOURCE_RUN/'formal_metadata.json').read_text(encoding='utf-8'))
        if not previous['completed'] or sha(SOURCE_RUN/'development_responses.csv') != previous['prepared_input_hashes']['development_responses.csv']:
            raise ValueError('Formal-v1 prepared response provenance failed')
        response = pd.read_csv(SOURCE_RUN/'development_responses.csv', dtype={'gdsc_drug_id': str}, keep_default_na=False, float_precision='round_trip')
        response = response[response.gdsc_drug_id.astype(int).isin(config['drug_ids'])].copy()
        if len(response) != lock['contract']['drug_model_observations']:
            raise ValueError('Fixed-eight response coverage changed')
        response.to_csv(out/'development_responses.csv', index=False, encoding='utf-8-sig', lineterminator='\n')
        prepared = {name: sha(out/name) for name in ('global_fold_manifest.csv', 'development_responses.csv')}
        atomic_json(out/'worker_context.json', dict(configuration=config, freeze_sha256=freeze_sha,
            execution_commit=revision, source_hashes=hashes, prepared_input_hashes=prepared))
        prepared['worker_context.json'] = sha(out/'worker_context.json')
        meta = dict(schema_version='1.0', run_name=config['run_name'], created_utc=datetime.now(timezone.utc).isoformat(),
            config_sha256=sha(config_path), freeze_sha256=freeze_sha, execution_commit=revision, source_hashes=hashes,
            configuration=config, development_input_verification=verification, cohort_sha256=verification['freeze_sha256'],
            python=platform.python_version(), packages={p: importlib.metadata.version(p) for p in ('numpy','pandas','scikit-learn','scipy','psutil','threadpoolctl')},
            expected_tasks=lock['contract']['tasks'], expected_fit_events=lock['contract']['fits'],
            expected_oof_predictions=lock['contract']['oof_predictions'], completed_tasks=0, completed=False, status='prepared',
            external_response_accessed=False, global_split_warnings=split_warnings, prepared_input_hashes=prepared,
            checkpoint_hashes={}, sessions=[], source_prepared_sha256=sha(SOURCE_RUN/'development_responses.csv'),
            source_formal_metadata_sha256=sha(SOURCE_RUN/'formal_metadata.json'))
        atomic_json(metadata_path, meta)
    response = pd.read_csv(out/'development_responses.csv', keep_default_na=False, float_precision='round_trip')
    folds = pd.read_csv(out/'global_fold_manifest.csv', keep_default_na=False)
    pending = []
    for index, task in enumerate(tasks):
        path = out/'checkpoints'/checkpoint_name(task)
        if path.exists():
            if validate_checkpoint(path, task, freeze_sha):
                verify_prepared_checkpoint(path, task, response, folds, config)
                meta['checkpoint_hashes'][path.name] = sha(path); continue
            archive = out/'failed_attempts'; archive.mkdir(exist_ok=True)
            destination = archive/(path.stem+'_'+sha(path)+'.json')
            if destination.exists():
                raise ValueError('Repeated failure payload already archived; inspect explicitly')
            os.replace(path, destination); meta['checkpoint_hashes'].pop(path.name, None)
        pending.append((index, task))
    session = dict(started_utc=datetime.now(timezone.utc).isoformat(), pending_at_start=len(pending),
                   completed_tasks_this_session=0, status='running', max_inflight_observed=1,
                   task_process_isolation=True, resource_wait_seconds=0., child_wall_seconds_sum=0.)
    meta['sessions'].append(session); meta.update(completed=False, status='running', completed_tasks=len(tasks)-len(pending))
    atomic_json(metadata_path, meta); started = time.perf_counter(); peak = private_peak = 0
    minimum = system_memory()['available_ram_bytes']; commit_minimum = system_memory()['system_commit_free_bytes']
    print(json.dumps(dict(event='ordinal_repeated_start', pending_tasks=len(pending), total_tasks=len(tasks), execution_commit=revision)), flush=True)
    exit_code = 0; child = None
    fields = ['elapsed_seconds','process_tree_rss_bytes','process_tree_private_bytes','process_tree_count',
              'available_ram_bytes','system_commit_free_bytes','system_commit_bytes','system_commit_limit_bytes','commit_measurement','inflight']
    try:
        with (out/f'resource_samples_session{len(meta["sessions"])}.csv').open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
            def record(inflight):
                nonlocal peak, private_peak, minimum, commit_minimum
                memory = process_tree_memory(); peak=max(peak,memory['process_tree_rss_bytes']); private_peak=max(private_peak,memory['process_tree_private_bytes'])
                minimum=min(minimum,memory['available_ram_bytes']); commit_minimum=min(commit_minimum,memory['system_commit_free_bytes'])
                writer.writerow(dict(elapsed_seconds=time.perf_counter()-started, inflight=inflight, **memory)); handle.flush()
                return memory
            for index, task in pending:
                wait_start = time.perf_counter()
                while True:
                    memory = record(0); reasons = pressure_reasons(memory, config['resources'])
                    if not reasons:
                        break
                    waited = time.perf_counter()-wait_start
                    print(json.dumps(dict(event='resource_wait', reasons=reasons, wait_seconds=waited, **memory)), flush=True)
                    if waited >= config['resources']['pressure_max_wait_seconds']:
                        session['resource_wait_seconds'] += waited
                        session['status']='resource_paused'; meta['status']='resource_paused'; exit_code=75
                        print(json.dumps(dict(event='resource_paused', completed_tasks=meta['completed_tasks'], pending_task=task)), flush=True)
                        break
                    time.sleep(config['resources']['pressure_poll_seconds'])
                if exit_code:
                    break
                session['resource_wait_seconds'] += time.perf_counter()-wait_start
                child_start = time.perf_counter()
                env=dict(os.environ, PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8')
                child=subprocess.Popen(child_command(out,index,meta['prepared_input_hashes']['worker_context.json']),cwd=ROOT,env=env)
                last_heartbeat=-60.
                while child.poll() is None:
                    record(1); elapsed=time.perf_counter()-started
                    if elapsed-last_heartbeat >= 60:
                        print(json.dumps(dict(event='heartbeat', completed_tasks=meta['completed_tasks'], elapsed_seconds=round(elapsed,1), inflight=1)),flush=True)
                        last_heartbeat=elapsed
                    time.sleep(1)
                child_wall=time.perf_counter()-child_start; session['child_wall_seconds_sum']+=child_wall
                path=out/'checkpoints'/checkpoint_name(task)
                if not path.exists():
                    atomic_json(path,dict(schema_version='1.0',task=task,freeze_sha256=freeze_sha,status='failed',oof_rows=[],result=None,fit_events=[],
                        error=f'Worker exited {child.returncode} without a complete checkpoint; partial-fit count unavailable',
                        resources=dict(wall_seconds=child_wall,cpu_seconds=0.,sampled_peak_rss_bytes=0),
                        data_guard=dict(external_response_accessed=False,blocked_attempts=0)))
                obj=json.loads(path.read_text(encoding='utf-8')); meta['checkpoint_hashes'][path.name]=sha(path)
                if child.returncode or not validate_checkpoint(path,task,freeze_sha):
                    print(json.dumps(dict(event='ordinal_repeated_task_finished',task=task,status='failed',
                        wall_seconds=obj['resources']['wall_seconds'],fit_count=len(obj['fit_events']),error=obj['error'])),flush=True)
                    raise RuntimeError(obj['error'] or 'Worker exited nonzero')
                session['completed_tasks_this_session']+=1
                meta['completed_tasks']=len(tasks)-len(pending)+session['completed_tasks_this_session']
                print(json.dumps(dict(event='ordinal_repeated_task_finished',task=task,status='completed',wall_seconds=obj['resources']['wall_seconds'],
                    child_wall_seconds=child_wall,fit_count=len(obj['fit_events']),sampled_peak_rss_bytes=obj['resources']['sampled_peak_rss_bytes'])),flush=True)
                del obj; atomic_json(metadata_path,meta); child=None
        if not exit_code:
            counts=aggregate(out,tasks,freeze_sha)
            if (counts['tasks'],counts['fit_events'],counts['oof_predictions']) != (meta['expected_tasks'],meta['expected_fit_events'],meta['expected_oof_predictions']):
                raise ValueError('Complete repeated task/fit/OOF contract mismatch')
            meta.update(completed=True,status='completed',output_counts=counts,
                output_hashes={name:sha(out/name) for name in ('oof_predictions.csv','fit_events.jsonl','outer_fit_audits.jsonl')},
                preserved_failed_attempts=len(list((out/'failed_attempts').glob('*.json'))))
            session['status']='completed'
    except BaseException:
        # Allow an already started child to finish before metadata is declared stopped.
        if child is not None and child.poll() is None:
            child.wait()
        meta.update(completed=False,status='failed',failure=traceback.format_exc());session['status']='failed';raise
    finally:
        session.update(wall_seconds=time.perf_counter()-started,finished_utc=datetime.now(timezone.utc).isoformat(),
            sampled_process_tree_peak_rss_bytes=peak,sampled_process_tree_peak_private_bytes=private_peak,
            minimum_available_ram_bytes=minimum,minimum_commit_free_bytes=commit_minimum,
            rss_scope='Sum parent/child RSS may double-count shared pages; sampled private/system commit recorded separately')
        meta['data_guard']=dict(external_response_accessed=False,blocked_attempts=len(blocked),accessed_path_count=len(access))
        atomic_json(metadata_path,meta)
    if not exit_code:
        print(json.dumps(dict(event='ordinal_repeated_completed',counts=meta['output_counts'],wall_seconds=session['wall_seconds'])),flush=True)
    return exit_code


if __name__ == '__main__':
    import sys
    sys.exit(main())
