"""Aggregate completed synthetic mechanism controls without inferential tests."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from ordinal_repeated_runtime import atomic_json, sha

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN = ROOT / 'runs/fixed8_synthetic_mechanism_controls_v1_1'
RESULTS = ROOT / 'results/synthetic_mechanisms'
REPORT = ROOT / 'reports/阶段9_合成机制控制.md'


def load_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding='utf-8'))


def atomic_csv(path: Path, rows: list[dict], fields: list[str]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='raise')
        writer.writeheader(); writer.writerows(rows)
    os.replace(temp, path)
    return sha(path)


def rho(x, y):
    x = np.asarray(x, dtype=np.float64); y = np.asarray(y, dtype=np.float64)
    if len(x) < 2 or np.ptp(x) == 0.0 or np.ptp(y) == 0.0:
        return None
    return float(spearmanr(x, y).statistic)


def point_metrics(y, pred):
    y = np.asarray(y, dtype=np.float64); pred = np.asarray(pred, dtype=np.float64)
    error = y - pred
    return dict(n=len(y), spearman=rho(y, pred), rmse=float(np.sqrt(np.mean(error ** 2))),
                mae=float(np.mean(np.abs(error))),
                r2=None if np.ptp(y) == 0.0 else float(1.0 - np.sum(error ** 2) / np.sum((y - np.mean(y)) ** 2)))


def summarize(values: list[float | None]) -> dict:
    valid = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    if not valid:
        return dict(n_defined=0, mean=None, median=None, minimum=None, maximum=None)
    return dict(n_defined=len(valid), mean=statistics.fmean(valid), median=statistics.median(valid),
                minimum=min(valid), maximum=max(valid))


def within_tissue_increment(y, prediction, tissue_null, tissue):
    y = np.asarray(y, dtype=np.float64); prediction = np.asarray(prediction, dtype=np.float64)
    tissue_null = np.asarray(tissue_null, dtype=np.float64); tissue = np.asarray(tissue)
    e0 = y - tissue_null; q = prediction - tissue_null
    per_rho = []; per_rmse = []
    for label in np.unique(tissue):
        rows = tissue == label
        per_rho.append(rho(e0[rows], q[rows]))
        model_error = y[rows] - prediction[rows]
        null_error = y[rows] - tissue_null[rows]
        per_rmse.append(float(np.sqrt(np.mean(null_error ** 2)) - np.sqrt(np.mean(model_error ** 2))))
    rho_available = [value for value in per_rho if value is not None]
    return dict(
        residual_rho_pooled=rho(e0, q),
        residual_rho_tissue_defined_count=len(rho_available),
        residual_rho_tissue_mean_available=(statistics.fmean(rho_available) if rho_available else None),
        residual_rho_tissue_macro_strict=(statistics.fmean(per_rho) if len(per_rho) == 8 and all(v is not None for v in per_rho) else None),
        residual_rho_tissue_macro_coverage=f'{len(rho_available)}/8',
        rmse_gain_vs_tissue_null=float(np.sqrt(np.mean((y - tissue_null) ** 2))
                                       - np.sqrt(np.mean((y - prediction) ** 2))),
        rmse_gain_tissue_macro=statistics.fmean(per_rmse) if len(per_rmse) == 8 else None,
    )


def analyze(run_dir: Path) -> dict:
    meta_path = run_dir / 'synthetic_mechanism_metadata.json'
    if not meta_path.is_file():
        raise FileNotFoundError('Synthetic run metadata is missing')
    metadata = load_json(meta_path)
    config = metadata['synthetic_config']
    expected = config['expected_contract']
    if metadata.get('status') != 'completed' or not metadata.get('completed'):
        raise ValueError('Do not analyze an incomplete synthetic run')
    if metadata.get('external_response_accessed') or metadata.get('real_development_data_accessed') or metadata.get('ctrp_accessed'):
        raise ValueError('Unexpected real/external data access recorded in synthetic run')
    tasks_by_name = metadata['tasks_by_name']
    tasks = list(tasks_by_name.values())
    if len(tasks) != expected['tasks']:
        raise ValueError('Task list differs from the frozen run contract')

    task_rows = []
    cohort_predictions = {}
    fit_groups = {}
    selected_outer = {}
    checkpoint_hashes = {}
    total_fit_count = total_oof_count = 0
    for task in tasks:
        name = (f"{task['scenario']}_r{task['replicate']:02d}_f{task['outer_fold']}_"
                f"{task['model']}.json")
        path = run_dir / 'checkpoints' / name
        expected_sha = metadata['registered_checkpoint_hashes'].get(name)
        if not expected_sha or sha(path) != expected_sha:
            raise ValueError('Checkpoint missing or hash differs from completed run metadata: ' + name)
        obj = load_json(path)
        if obj.get('status') != 'completed' or obj.get('task') != task:
            raise ValueError('Failed/incompatible task checkpoint: ' + name)
        if len(obj['fit_events']) != int(obj['fit_count']):
            raise ValueError('Fit event audit count differs from task checkpoint: ' + name)
        total_fit_count += int(obj['fit_count']); total_oof_count += len(obj['oof_rows'])
        checkpoint_hashes[name] = expected_sha
        result = obj['result']; metrics = result['metrics']; chosen = result.get('chosen_params') or {}
        feature = result.get('outer_feature_audit_compact') or {}
        outer = result.get('outer_fit_audit_compact') or {}
        task_rows.append(dict(scenario=task['scenario'], scenario_index=task['scenario_index'],
            replicate=task['replicate'], outer_fold=task['outer_fold'], model=task['model'],
            n=metrics.get('n'), spearman=metrics.get('spearman'), rmse=metrics.get('rmse'),
            mae=metrics.get('mae'), r2=metrics.get('r2'),
            chosen_alpha_ratio=chosen.get('alpha_ratio'), chosen_l1_ratio=chosen.get('l1_ratio'),
            selected_gene_count=feature.get('selected_genes'), selected_gene_ids_sha256=feature.get('selected_gene_ids_sha256'),
            pair_feature_count=feature.get('output_features'),
            outer_nonzero_molecular_coefficients=outer.get('n_nonzero_molecular_coeffs'),
            outer_zero_molecular_solution=(None if task['model'] == 'tissue_null' else outer.get('n_nonzero_molecular_coeffs') == 0),
            fit_count=obj['fit_count'], wall_seconds=obj['resources']['wall_seconds'],
            sampled_peak_rss_bytes=obj['resources']['sampled_peak_rss_bytes'],
            sampled_peak_private_bytes=obj['resources']['sampled_peak_private_bytes']))
        cohort_key = (task['scenario'], int(task['replicate']), task['model'])
        cohort_predictions.setdefault(cohort_key, []).extend(obj['oof_rows'])
        fit_groups.setdefault((task['scenario'], task['model']), []).extend(obj['fit_events'])
        genes = feature.get('gene_ids', [])
        selected_outer.setdefault((task['scenario'], int(task['replicate']), task['model']), {})[task['outer_fold']] = set(genes)

    if (total_fit_count, total_oof_count) != (expected['fits'], expected['per_sample_model_oof_predictions']):
        raise ValueError('Aggregated fit/OOF count differs from frozen synthetic contract')

    # Verify each cohort/model has complete one-time outer coverage and common labels.
    cohort_rows = {}
    for scenario in config['scenario_order']:
        for replicate in range(int(config['replicates_per_scenario'])):
            cohort_key = (scenario, replicate)
            by_model = {}
            canonical_y = canonical_tissue = canonical_null = None
            for model in config['models']:
                rows = cohort_predictions[(scenario, replicate, model)]
                rows = sorted(rows, key=lambda row: int(row['synthetic_sample_id']))
                ids = [int(row['synthetic_sample_id']) for row in rows]
                if ids != list(range(1, int(config['n_samples']) + 1)):
                    raise ValueError(f'Outer OOF predictions do not cover each synthetic sample once: {scenario} replicate {replicate} {model}')
                y = np.asarray([row['observed_response'] for row in rows], dtype=np.float64)
                pred = np.asarray([row['predicted_response'] for row in rows], dtype=np.float64)
                tissue = np.asarray([row['tissue'] for row in rows], dtype=np.int16)
                if canonical_y is None:
                    canonical_y = y; canonical_tissue = tissue
                elif not np.array_equal(y, canonical_y) or not np.array_equal(tissue, canonical_tissue):
                    raise ValueError('Model arms do not share exact synthetic labels and tissues')
                by_model[model] = dict(rows=rows, y=y, prediction=pred, metrics=point_metrics(y, pred))
            null_model = by_model['tissue_null']
            null_pred = null_model['prediction']
            for model in config['models']:
                item = by_model[model]
                increments = within_tissue_increment(canonical_y, item['prediction'], null_pred, canonical_tissue)
                metric = item['metrics']
                cohort_rows[(scenario, replicate, model)] = dict(
                    scenario=scenario, replicate=replicate, model=model, n=metric['n'],
                    spearman=metric['spearman'], rmse=metric['rmse'], mae=metric['mae'], r2=metric['r2'],
                    residual_rho_pooled=increments['residual_rho_pooled'],
                    residual_rho_tissue_defined_count=increments['residual_rho_tissue_defined_count'],
                    residual_rho_tissue_mean_available=increments['residual_rho_tissue_mean_available'],
                    residual_rho_tissue_macro_strict=increments['residual_rho_tissue_macro_strict'],
                    residual_rho_tissue_macro_coverage=increments['residual_rho_tissue_macro_coverage'],
                    rmse_gain_vs_tissue_null=increments['rmse_gain_vs_tissue_null'],
                    rmse_gain_tissue_macro=increments['rmse_gain_tissue_macro'],
                )

    task_fields = list(task_rows[0])
    result_hashes = {}
    result_hashes['performance_by_outer_task.csv'] = atomic_csv(
        RESULTS / 'performance_by_outer_task.csv', task_rows, task_fields)
    cohort_table = list(cohort_rows.values())
    cohort_fields = list(cohort_table[0])
    result_hashes['performance_by_cohort.csv'] = atomic_csv(
        RESULTS / 'performance_by_cohort.csv', cohort_table, cohort_fields)

    summary_rows = []
    metrics_to_summarize = ('spearman', 'rmse', 'mae', 'r2', 'residual_rho_pooled',
        'residual_rho_tissue_mean_available', 'residual_rho_tissue_macro_strict',
        'rmse_gain_vs_tissue_null', 'rmse_gain_tissue_macro')
    for scenario in config['scenario_order']:
        for model in config['models']:
            rows = [cohort_rows[(scenario, replicate, model)] for replicate in range(10)]
            summary = dict(scenario=scenario, model=model, cohorts=len(rows))
            for metric in metrics_to_summarize:
                values = summarize([row[metric] for row in rows])
                prefix = metric + '_'
                for stat_name, stat_value in values.items():
                    summary[prefix + stat_name] = stat_value
            summary_rows.append(summary)
    summary_fields = list(summary_rows[0])
    result_hashes['scenario_summary.csv'] = atomic_csv(
        RESULTS / 'scenario_summary.csv', summary_rows, summary_fields)

    fit_rows = []
    for scenario in config['scenario_order']:
        for model in config['models']:
            events = fit_groups[(scenario, model)]
            molecular = model != 'tissue_null'
            nonzero = [event.get('n_nonzero_molecular_coeffs') for event in events
                       if event.get('n_nonzero_molecular_coeffs') is not None]
            path_ratios = [event.get('path_ratio') for event in events if event.get('path_ratio') is not None]
            fit_rows.append(dict(scenario=scenario, model=model, total_fits=len(events),
                molecular_fits=len(nonzero) if molecular else 0,
                zero_molecular_fits=sum(int(value == 0) for value in nonzero) if molecular else None,
                zero_molecular_fit_fraction=(sum(int(value == 0) for value in nonzero) / len(nonzero) if molecular and nonzero else None),
                nonzero_molecular_fits=sum(int(value > 0) for value in nonzero) if molecular else None,
                alpha_max_null_fits=sum(int(event.get('solver_reason') == 'zero_training_alpha_max') for event in events),
                theoretical_null_endpoint_fits=sum(int(bool(event.get('null_endpoint'))) for event in events),
                fit_warning_count=sum(len(event.get('warnings') or []) for event in events),
                fit_fallback_count=sum(int(event.get('fitted_kind') == 'mean'
                    or str(event.get('fallback_reason') or '').startswith('empty_features')) for event in events),
                fit_nonzero_alpha_ratio_count=sum(int(value < 1.0) for value in path_ratios),
                fit_alpha_ratio_1_count=sum(int(value == 1.0) for value in path_ratios),
                median_fit_seconds=statistics.median([float(event.get('fit_seconds') or 0.0) for event in events]),
                median_prediction_seconds=statistics.median([float(event.get('predict_seconds') or 0.0) for event in events])))
    result_hashes['fit_diagnostics.csv'] = atomic_csv(RESULTS / 'fit_diagnostics.csv', fit_rows, list(fit_rows[0]))

    selection_rows = []
    for scenario in config['scenario_order']:
        for model in config['models']:
            model_outer = [(replicate, selected_outer[(scenario, replicate, model)]) for replicate in range(10)]
            all_sets = [genes for _, folds in model_outer for genes in folds.values()]
            driver1 = 'SIMGENE_0001'; driver2 = 'SIMGENE_0002'
            per_cohort_jaccard = []
            for replicate, folds in model_outer:
                fold_sets = [folds[index] for index in sorted(folds)]
                pairs = []
                for i in range(len(fold_sets)):
                    for j in range(i + 1, len(fold_sets)):
                        union = fold_sets[i] | fold_sets[j]
                        pairs.append(len(fold_sets[i] & fold_sets[j]) / len(union) if union else None)
                valid = [value for value in pairs if value is not None]
                if valid:
                    per_cohort_jaccard.append(statistics.fmean(valid))
            denom = len(all_sets)
            selection_rows.append(dict(scenario=scenario, model=model, outer_fits=denom,
                mean_selected_gene_count=statistics.fmean([len(values) for values in all_sets]),
                driver_1_selected_count=sum(driver1 in values for values in all_sets),
                driver_1_selected_frequency=(sum(driver1 in values for values in all_sets) / denom if denom else None),
                driver_2_selected_count=sum(driver2 in values for values in all_sets),
                driver_2_selected_frequency=(sum(driver2 in values for values in all_sets) / denom if denom else None),
                within_cohort_outer_fold_jaccard_mean=(statistics.fmean(per_cohort_jaccard) if per_cohort_jaccard else None),
                within_cohort_outer_fold_jaccard_min=(min(per_cohort_jaccard) if per_cohort_jaccard else None),
                within_cohort_outer_fold_jaccard_max=(max(per_cohort_jaccard) if per_cohort_jaccard else None)))
    result_hashes['gene_selection_summary.csv'] = atomic_csv(
        RESULTS / 'gene_selection_summary.csv', selection_rows, list(selection_rows[0]))

    data_checks = []
    for cohort_key, record in metadata['cohort_records'].items():
        check = record['checks']
        data_checks.append(dict(cohort=cohort_key, scenario=check['scenario'], replicate=check['replicate'],
            seed=check['seed'], n_samples=check['n_samples'], n_genes=check['n_genes'], n_tissues=check['n_tissues'],
            samples_per_tissue_min=check['samples_per_tissue_min'], samples_per_tissue_max=check['samples_per_tissue_max'],
            tpm_min=check['tpm_min'], tpm_max=check['tpm_max'], log2_tpm_min=check['log2_tpm_min'],
            response_min=check['response_min'], response_max=check['response_max'],
            within_tissue_response_max_ptp=check['within_tissue_response_max_ptp'],
            common_amplitude_rank_invariant=check['common_amplitude_rank_invariant'],
            common_amplitude_pair_invariant=check['common_amplitude_pair_invariant'],
            cohort_file_sha256=record['cohort_file_sha256'], split_file_sha256=record['split_file_sha256']))
    result_hashes['data_generation_checks.csv'] = atomic_csv(
        RESULTS / 'data_generation_checks.csv', data_checks, list(data_checks[0]))

    output_counts = dict(tasks=len(task_rows), fits=total_fit_count, oof_predictions=total_oof_count,
        cohort_summaries=len(cohort_table), scenario_model_summaries=len(summary_rows),
        all_checkpoints_hash_verified=True)
    if output_counts['tasks'] != expected['tasks']:
        raise ValueError('Aggregate task count mismatch')
    write_report(run_dir, metadata, summary_rows, fit_rows, selection_rows, output_counts, result_hashes)
    result_hashes['reports/阶段9_合成机制控制.md'] = sha(REPORT)
    manifest = dict(schema_version='1.0', status='completed', analysis_time_utc=datetime.now(timezone.utc).isoformat(),
        run_name=metadata['run_name'], execution_lock_sha256=metadata['freeze_sha256'],
        execution_commit=metadata['execution_commit'], metadata_sha256=sha(meta_path),
        checkpoint_hashes=checkpoint_hashes, outputs=result_hashes, counts=output_counts,
        inferential_tests_performed=False, p_values_reported=False, confidence_intervals_reported=False,
        external_response_accessed=False, real_development_data_accessed=False, ctrp_accessed=False,
        interpretation_limit='Finite-sample behavior under prespecified synthetic DGPs; no biological inference or evidence of general REO superiority')
    atomic_json(RESULTS / 'synthetic_mechanism_manifest.json', manifest)
    result_hashes['synthetic_mechanism_manifest.json'] = sha(RESULTS / 'synthetic_mechanism_manifest.json')
    return dict(counts=output_counts, outputs=result_hashes, summary_rows=summary_rows,
                fit_rows=fit_rows, selection_rows=selection_rows, runtime=metadata.get('sessions', []))


def write_report(run_dir: Path, metadata: dict, summaries: list[dict], fit_rows: list[dict], selections: list[dict],
                  counts: dict, outputs: dict) -> None:
    session_seconds = sum(float(session.get('wall_seconds', 0.0)) for session in metadata.get('sessions', []))
    max_rss = max((int(row.get('sampled_peak_rss_bytes') or 0)
                   for row in metadata.get('tasks_by_name', {}).values()), default=0)
    # Per-task resource maxima are read from checkpoints, avoiding assumptions about shared memory.
    task_rss = []
    for name in metadata.get('registered_checkpoint_hashes', {}):
        try:
            obj = load_json(run_dir / 'checkpoints' / name)
            task_rss.append(int(obj.get('resources', {}).get('sampled_peak_rss_bytes', 0)))
        except (OSError, ValueError):
            continue
    if task_rss:
        max_rss = max(task_rss)
    max_tree_rss = max((int(session.get('sampled_process_tree_peak_rss_bytes', 0))
                        for session in metadata.get('sessions', [])), default=0)
    max_tree_private = max((int(session.get('sampled_process_tree_peak_private_bytes', 0))
                            for session in metadata.get('sessions', [])), default=0)
    minimum_available = min((int(session.get('minimum_available_ram_bytes', 0))
                            for session in metadata.get('sessions', [])), default=0)
    minimum_commit = min((int(session.get('minimum_commit_free_bytes', 0))
                         for session in metadata.get('sessions', [])), default=0)
    model_map = {(row['scenario'], row['model']): row for row in summaries}
    lines = [
        '# 阶段9：有限合成机制控制', '',
        f"执行锁：`{metadata['freeze_sha256']}`；执行提交：`{metadata['execution_commit']}`。",
        f"合成数据队列{len(metadata['cohort_records'])}个，模型task {counts['tasks']}个，fit {counts['fits']:,}次，OOF {counts['oof_predictions']:,}条；真实开发/外部数据拟合0次，CTRP访问0次。",
        '该结果只描述预设有限数据生成机制下的实现行为，不构成对真实药物机制或REO一般优势的推断；未计算p值和推断置信区间。', '',
        '## 模拟行为摘要', '',
        '|情形|模型|OOF Spearman均值（10队列）|组织内残差ρ可用组织均值|RMSE均值|组织参照RMSE改善均值|零分子fit比例|',
        '|---|---|---:|---:|---:|---:|---:|'
    ]
    fit_map = {(row['scenario'], row['model']): row for row in fit_rows}
    for scenario in metadata['synthetic_config']['scenario_order']:
        for model in metadata['synthetic_config']['models']:
            row = model_map[(scenario, model)]; fits = fit_map[(scenario, model)]
            values = [row[f'{key}_mean'] for key in ('spearman', 'residual_rho_tissue_mean_available', 'rmse', 'rmse_gain_vs_tissue_null')]
            rendered = [f'{value:.4f}' if value is not None else 'NA' for value in values]
            zero = fits['zero_molecular_fit_fraction']
            lines.append(f"|{scenario}|{model}|{rendered[0]}|{rendered[1]}|{rendered[2]}|{rendered[3]}|{f'{zero:.3f}' if zero is not None else '—'}|")
    lines.extend(['', '## DGP完整性与特征选择', '',
        f"- 纯组织响应在每个组织内的最大range严格为0；每个队列均为8组织×30样本，TPM有限非负。",
        f"- 共同幅度情形中，10/10队列的220基因样本内秩和全部pair符号逐样本完全不变（由生成阶段直接验收）。",
        '- 预设顺序情形的两个driver在外层训练候选基因集中的入选频率见`gene_selection_summary.csv`；这是该固定DGP下的选择诊断，不作为成功门槛。',
        '- 零解、非零解、解析零端点、求解警告、回退和全部内/外层候选拟合均保留在本地任务检查点；公开聚合计数见`fit_diagnostics.csv`。', '',
        '## 执行与资源', '',
        f"控制器累计运行约{session_seconds / 3600:.2f}小时；单worker、每task独立进程、BLAS单线程。worker采样RSS最大{max_rss / 1024**3:.3f} GiB；控制器/worker进程树RSS采样最大{max_tree_rss / 1024**3:.3f} GiB，private最大{max_tree_private / 1024**3:.3f} GiB。系统观测的最低可用物理内存{minimum_available / 1024**3:.3f} GiB，最低commit余量{minimum_commit / 1024**3:.3f} GiB。采样值是操作记录，不是进程隔离性能基准。", '',
        '## 解释和下一步', '',
        '情形间差异仅可用于解释当前编码与训练流程在明确机制下可能如何退化：共同样本幅度特意令秩/pair无法观察该幅度；预设pair顺序特意包含比较驱动；纯组织情形检查模型是否回到组织参照；弱分子信号展示有限样本筛选和拟合波动。不得把任一情形的预期表现外推为真实药物机制。',
        '按总计划，下一阶段是固定八药、23组织逐一留出的LOTO敏感性分析。不得据合成排名改候选规则、网格或药物；完成LOTO后另行冻结最终模型与外验缺失/弃判方案，通过gate前CTRP继续封存。', '',
        '### 聚合文件', '',
        '- `../results/synthetic_mechanisms/performance_by_outer_task.csv`: 每个外层task指标/参数/资源。',
        '- `../results/synthetic_mechanisms/performance_by_cohort.csv`: 每个合成队列×模型完整240样本OOF指标、组织残差增量和误差改善。',
        '- `../results/synthetic_mechanisms/scenario_summary.csv`: 每情形×模型跨10队列的描述性均值、中位数、范围和定义数。',
        '- `../results/synthetic_mechanisms/fit_diagnostics.csv`: 分子零解、警告、回退及拟合耗时。',
        '- `../results/synthetic_mechanisms/gene_selection_summary.csv`: 预设driver选择频率与外折间Jaccard稳定性。',
        '- `../results/synthetic_mechanisms/data_generation_checks.csv`: 40个队列的seed、输入范围及DGP不变量。', '',
        '结果、代码和报告的SHA256记录在`../results/synthetic_mechanisms/synthetic_mechanism_manifest.json`。本地逐样本OOF与checkpoint保留在忽略的`runs/`目录，不上传Git。',
    ])
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()
    result = analyze(args.run_dir.resolve())
    print(json.dumps(dict(status='completed', counts=result['counts'], outputs=result['outputs']), ensure_ascii=False))


if __name__ == '__main__':
    main()
