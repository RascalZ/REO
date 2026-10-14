#!/usr/bin/env python3
"""Summarize zero-endpoint selection and fixed-panel heterogeneity from v1 OOF."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import statistics
import subprocess
import sys
import threading
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

try:
    import psutil
except ImportError as exc:  # pragma: no cover - depends on managed analysis env
    raise SystemExit("psutil is required; run with the repository .venv Python") from exc


REPO = Path(__file__).resolve().parents[1]
RUN_NAME = (
    "gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-"
    "nested5x3-inner3-repeated-v1"
)
RESULT_DIR = REPO / "results" / "ordinal_repeated"
RUN_DIR = REPO / "runs" / RUN_NAME
OUT_DIR = REPO / "results" / "development_diagnostics"
REPORT_PATH = REPO / "reports" / "阶段7_零解与药物异质性诊断.md"
DRUG_ORDER = [1931, 1036, 1191, 1930, 1199, 1091, 1010, 1922]
MODEL_ORDER = [
    "expression_strong",
    "mad_expression",
    "mad_rank",
    "mad_pair_full",
    "response_association_expression",
    "response_association_rank",
    "response_association_pair_full",
]
MOLECULAR_ORDER = [m for m in MODEL_ORDER if m != "expression_strong"]
MODEL_LABELS = {
    "tissue_null": "组织参照",
    "expression_strong": "强表达2000",
    "mad_expression": "MAD表达200",
    "mad_rank": "MAD秩200",
    "mad_pair_full": "MAD pair",
    "response_association_expression": "关联表达200",
    "response_association_rank": "关联秩200",
    "response_association_pair_full": "关联pair",
}
DRUG_LABELS = {
    1931: "MIRA-1",
    1036: "PLX-4720",
    1191: "Bortezomib",
    1930: "Telomerase Inhibitor IX",
    1199: "Tamoxifen",
    1091: "BMS-536924",
    1010: "Gefitinib",
    1922: "Cediranib",
}


class PeakRSS:
    def __init__(self) -> None:
        self.process = psutil.Process()
        self.stop_event = threading.Event()
        self.peak = 0
        self.thread = threading.Thread(target=self._sample, daemon=True)

    def _sample(self) -> None:
        while not self.stop_event.is_set():
            try:
                self.peak = max(self.peak, self.process.memory_info().rss)
            except psutil.Error:
                pass
            self.stop_event.wait(0.1)

    def start(self) -> None:
        self.thread.start()

    def stop(self) -> int:
        self.stop_event.set()
        self.thread.join()
        return max(self.peak, self.process.memory_info().rss)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=columns, extrasaction="raise", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fnum(value: str | float | int | None) -> float | None:
    if value in (None, ""):
        return None
    return float(value)


def mean_metric(rows: list[dict[str, str]], drug: int, model: str, metric: str) -> float:
    values = [
        float(row[metric])
        for row in rows
        if int(row["gdsc_drug_id"]) == drug
        and row["model"] == model
        and row[metric] != ""
    ]
    if len(values) != 3:
        raise ValueError(f"Expected all three repeats for {drug}/{model}/{metric}, got {len(values)}")
    return statistics.fmean(values)


def render_report(
    profiles: list[dict],
    zero_summary: list[dict],
    repeat_summary: list[dict],
    decomposition: list[dict],
    path_margin_rows: list[dict],
    manifest: dict,
) -> str:
    totals = {
        "outer_fits": sum(int(row["outer_fits"]) for row in zero_summary),
        "zero": sum(int(row["zero_molecular_fits"]) for row in zero_summary),
        "nonzero": sum(int(row["nonzero_molecular_fits"]) for row in zero_summary),
        "alpha1": sum(int(row["alpha_ratio_1_selected_fits"]) for row in zero_summary),
        "warn": sum(int(row["warning_outer_fits"]) for row in zero_summary),
        "fallback": sum(int(row["fallback_outer_fits"]) for row in zero_summary),
    }
    tissues = sorted({int(row["n_tissues"]) for row in decomposition})
    lines = [
        "# 阶段7：零解与药物异质性诊断",
        "",
        "日期：2026-10-14。依据开发集固定八药重复评价已完成的fit审计、公开聚合结果和开发响应矩阵；未新拟合模型，未读取CTRP响应、外部表达或外验索引。计算规则见[开发集诊断协议v1](../docs/DEVELOPMENT_DIAGNOSTICS_PROTOCOL_V1.md)，唯一总计划为[RESEARCH_PLAN.md](../RESEARCH_PLAN.md)。",
        "",
        "## 零解主要来自预设内层零端点",
        "",
        f"七个分子臂共有{totals['outer_fits']}个外层fit，其中{totals['zero']}个零分子解、{totals['nonzero']}个非零解。审计台账显示{totals['alpha1']}个fit选择`alpha_ratio=1`解析零端点；solver warning={totals['warn']}、fallback={totals['fallback']}。因此这批零解记录的是当前数据/候选/内层RMSE路径所作的模型选择，不是数值失败或真实药物缺少分子效应的证明。",
        "",
        "每个药物×分子臂的表同时报告零解频率和内层“最佳非零候选RMSE−最佳零端点RMSE”。该差值为负表示本组内层OOF中某个非零候选较低；它同样参与了超参数选择，因此只用于解释选择竞争，不能作为独立验证或P值。外层上下文彼此重叠，不能把15个context当独立样本。",
        "",
        "|药物|模型|零解/外层fit|选alpha_ratio=1|非零路径低于零端点的内层context|最佳非零−零RMSE中位数|",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in zero_summary:
        lines.append(
            f"|{row['drug_name']}|{MODEL_LABELS[row['model']]}|"
            f"{row['zero_molecular_fits']}/{row['outer_fits']}|"
            f"{row['alpha_ratio_1_selected_fits']}/{row['outer_fits']}|"
            f"{row['inner_nonzero_better_contexts']}/{row['outer_fits']}|"
            f"{float(row['median_best_nonzero_minus_zero_rmse']):+.6g}|"
        )
    lines += [
        "",
        "按repeat拆分的结果与模型失败/回退字段保存在[`zero_solution_by_repeat.csv`](../results/development_diagnostics/zero_solution_by_repeat.csv)；内层竞争摘要保存在[`inner_path_competition.csv`](../results/development_diagnostics/inner_path_competition.csv)。",
        "",
        "## 固定药物面板存在不同的组织和表现构成",
        "",
        "组织间/组织内平方和按每药开发响应直接拆分，是当前开发样本的描述性构成。小组织保留；该比例不是药物生物机制归因。主比较点值沿既有三个完整repeat分别计算后取均值；不新算区间或P值。RMSE gain定义为强表达RMSE减关联pair RMSE，正值有利于关联pair。",
        "",
        "|药物|N|组织数|组织间SS/总SS|主Δρ关联pair−强表达|三repeat主Δρ|关联pair−组织参照Δρ|RMSE改善（强表达−关联pair）|关联pair零解/15|",
        "|---|---:|---:|---:|---:|---|---:|---:|---:|",
    ]
    for row in profiles:
        reps = ", ".join(f"{float(row[f'main_delta_rho_repeat_{i}']):+.3f}" for i in range(1, 4))
        lines.append(
            f"|{row['drug_name']}|{row['n_observations']}|{row['n_tissues']}|"
            f"{float(row['between_tissue_ss_fraction']):.3f}|"
            f"{float(row['main_delta_rho']):+.3f}|{reps}|"
            f"{float(row['association_pair_minus_tissue_rho']):+.3f}|"
            f"{float(row['strong_expression_minus_association_pair_rmse']):+.5f}|"
            f"{row['association_pair_zero_fits']}/15|"
        )
    lines += [
        "",
        "固定面板药物之间的排序、RMSE、组织贡献和零解并不一致；不得以少数较高点差概括所有药物。完整表示×药物长表包括强表达/组织null/两候选规则六臂的Spearman、RMSE、零解、RMSE与ρ配对差、稳定性中位数及预测偏离；见[`drug_heterogeneity.csv`](../results/development_diagnostics/drug_heterogeneity.csv)和[`response_decomposition.csv`](../results/development_diagnostics/response_decomposition.csv)。原始条件区间继续以阶段6冻结结果为准。",
        "",
        "## 解释边界和下一步",
        "",
        "此结果支持对“模型何时选择组织参照”和药物间差异作更清楚的描述；不支持宣称模型确诊了某个药物的真实生物信号。候选基因、模型臂和数据来自已经看过的开发结果，这一诊断属于假设生成。15个外层拟合上下文不是独立观测；表中效应点值不增加原统计家族，也不替代冻结区间。",
        "",
        "下一步按总计划进入有限合成统计机制检验与组织敏感性/LOTO。响应关联候选规则的八药真实拟合已经完成；只有出现明确可区分的新假设才另行提议第三规则。留出单个组织时组织null为常数，其组织内Spearman为NA；该NA不编码为0。外验继续封存。",
        "",
        "## 运行与来源",
        "",
        f"分析程序SHA256：`{manifest['analysis_script_sha256']}`；执行Git提交：`{manifest['execution_commit']}`；处理器峰值RSS（100ms进程采样）：{manifest['peak_rss_bytes']:,} bytes；运行{manifest['elapsed_seconds']:.2f}秒。SHA/行数清单见[`diagnostics_manifest.json`](../results/development_diagnostics/diagnostics_manifest.json)。输入校验与总计必须闭合：零分子外层fit {totals['zero']}，选择零端点 {totals['alpha1']}，响应拆分覆盖8药。",
        "",
    ]
    if len(tissues) != 1 or tissues[0] != 23:
        raise ValueError(f"Expected 23 development tissues for every drug; got {tissues}")
    return "\n".join(lines)


def run() -> dict:
    started = time.perf_counter()
    peak = PeakRSS()
    peak.start()

    input_paths = {
        "config_lock": REPO / "config" / "ordinal_repeated_v1.lock.json",
        "analysis_config": REPO / "config" / "ordinal_repeated_v1.json",
        "analysis_provenance": RESULT_DIR / "analysis_provenance.json",
        "repeat_metrics": RESULT_DIR / "repeat_metrics.csv",
        "paired_effects": RESULT_DIR / "paired_effect_intervals.csv",
        "zero_summary": RESULT_DIR / "molecular_and_failure_summary.csv",
        "chosen_parameters": RESULT_DIR / "chosen_parameters.csv",
        "feature_stability": RESULT_DIR / "feature_stability.csv",
        "prediction_deviation": RESULT_DIR / "prediction_deviation_from_tissue.csv",
        "outer_fit_audits": RUN_DIR / "outer_fit_audits.jsonl",
        "development_responses": RUN_DIR / "development_responses.csv",
    }
    for path in input_paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    provenance = json.loads(input_paths["analysis_provenance"].read_text(encoding="utf-8"))
    lock = json.loads(input_paths["config_lock"].read_text(encoding="utf-8"))
    if provenance["status"] != "completed" or provenance["external_response_accessed"]:
        raise ValueError("The fixed-eight development run must be complete and CTRP access must remain false")
    if provenance["config_sha256"] != lock["files"]["config/ordinal_repeated_v1.json"]:
        raise ValueError("Analysis provenance is not bound to the fixed-eight lock")
    if provenance["lock_sha256"] != sha256(input_paths["config_lock"]):
        raise ValueError("The fixed-eight config lock hash differs from the completed provenance")
    config = json.loads(input_paths["analysis_config"].read_text(encoding="utf-8"))
    if config["run_name"] != RUN_NAME or config["external_response_access"] != "forbidden":
        raise ValueError("Analysis config does not match the sealed fixed-eight development run")
    if int(provenance["new_model_fits"]) != 0:
        raise ValueError("This diagnostic stage must not fit new models")

    metrics = read_csv(input_paths["repeat_metrics"])
    zero_source = read_csv(input_paths["zero_summary"])
    params_source = read_csv(input_paths["chosen_parameters"])
    stability_source = read_csv(input_paths["feature_stability"])
    deviation_source = read_csv(input_paths["prediction_deviation"])

    zero_by_key = {(int(r["gdsc_drug_id"]), r["model"], int(r["repeat"])): r for r in zero_source}
    selected_by_key = {(int(r["gdsc_drug_id"]), r["model"], int(r["repeat"]), int(r["outer_fold"])): r for r in params_source}
    if len(zero_by_key) != 8 * 8 * 3 or len(selected_by_key) != 960:
        raise ValueError("The pre-existing zero/parameter summaries do not cover the frozen design")

    context_rows: list[dict] = []
    audit_task_count = 0
    tissue_null_task_count = 0
    with input_paths["outer_fit_audits"].open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            item = json.loads(line)
            audit_task_count += 1
            task = item["task"]
            guard = item.get("data_guard") or {}
            if guard.get("external_response_accessed") or guard.get("blocked_attempts"):
                raise ValueError(f"Unexpected external access guard state in audit line {line_number}")
            if int(task.get("permutation_id", 0)) != 0:
                continue
            model = task["model"]
            if model == "tissue_null":
                tissue_null_task_count += 1
                continue
            drug = int(task["drug_id"])
            repeat = int(task["repeat"])
            outer_fold = int(task["outer_fold"])
            key = (drug, model, repeat, outer_fold)
            chosen = selected_by_key.get(key)
            if chosen is None:
                raise ValueError(f"Missing selected-parameter record for outer audit line {line_number}")
            scores = item["inner_scores"]
            if len(scores) != 10:
                raise ValueError(f"Expected 10 inner candidates in audit line {line_number}")
            null_scores = [float(r["pooled_inner_rmse"]) for r in scores if float(r["params"]["alpha_ratio"]) == 1.0]
            nonzero_scores = [float(r["pooled_inner_rmse"]) for r in scores if float(r["params"]["alpha_ratio"]) < 1.0]
            if len(null_scores) != 2 or len(nonzero_scores) != 8:
                raise ValueError(f"Unexpected zero/nonzero grid composition in audit line {line_number}")
            best_zero = min(null_scores)
            best_nonzero = min(nonzero_scores)
            chosen_params = json.loads(chosen["chosen_params"])
            audit = item["outer_fit_audit"]
            zero_selected = float(chosen_params["alpha_ratio"]) == 1.0
            exact_zero = int(audit["n_nonzero_molecular_coeffs"]) == 0
            null_endpoint = bool(audit["null_endpoint"])
            if zero_selected != exact_zero or zero_selected != null_endpoint:
                raise ValueError(f"Selected alpha, coefficients and null-endpoint state disagree: {key}")
            context_rows.append({
                "gdsc_drug_id": drug,
                "drug_name": DRUG_LABELS[drug],
                "model": model,
                "repeat": repeat,
                "outer_fold": outer_fold,
                "selected_alpha_ratio": float(chosen_params["alpha_ratio"]),
                "selected_l1_ratio": float(chosen_params["l1_ratio"]),
                "zero_molecular_solution": zero_selected,
                "analytical_null_endpoint": null_endpoint,
                "warning_count": len(audit.get("warnings") or []),
                "fallback_reason": audit.get("fallback_reason") or "",
                "best_zero_inner_rmse": best_zero,
                "best_nonzero_inner_rmse": best_nonzero,
                "best_nonzero_minus_zero_rmse": best_nonzero - best_zero,
                "nonzero_candidate_strictly_better": best_nonzero < best_zero,
                "inner_candidate_count": len(scores),
            })
    expected_molecular_outer = 8 * 7 * 15
    if (len(context_rows), tissue_null_task_count, audit_task_count) != (expected_molecular_outer, 120, 960):
        raise ValueError("Outer fit audit rows do not cover the complete frozen 960-task design")

    context_groups: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for row in context_rows:
        context_groups[(row["gdsc_drug_id"], row["model"])].append(row)
    inner_rows = []
    zero_rows = []
    repeat_rows = []
    for drug in DRUG_ORDER:
        for model in MOLECULAR_ORDER:
            contexts = context_groups[(drug, model)]
            if len(contexts) != 15:
                raise ValueError(f"Expected 15 audit contexts for {drug}/{model}; got {len(contexts)}")
            diffs = [row["best_nonzero_minus_zero_rmse"] for row in contexts]
            inner_rows.append({
                "gdsc_drug_id": drug,
                "drug_name": DRUG_LABELS[drug],
                "model": model,
                "outer_contexts": len(contexts),
                "nonzero_candidate_better_contexts": sum(row["nonzero_candidate_strictly_better"] for row in contexts),
                "best_nonzero_minus_zero_rmse_mean": statistics.fmean(diffs),
                "best_nonzero_minus_zero_rmse_median": statistics.median(diffs),
                "best_nonzero_minus_zero_rmse_min": min(diffs),
                "best_nonzero_minus_zero_rmse_max": max(diffs),
                "interpretation": "negative_favors_nonzero_within_same_inner_selection_data_only",
            })
            summary = [zero_by_key[(drug, model, repeat)] for repeat in range(3)]
            alpha_one = sum(row["zero_molecular_solution"] for row in contexts)
            nzero = sum(int(row["zero_molecular_fits"]) for row in summary)
            nouter = sum(int(row["outer_fits"]) for row in summary)
            if nouter != 15 or alpha_one != nzero:
                raise ValueError(f"Zero summaries disagree with fit-level selected endpoint for {drug}/{model}")
            zero_rows.append({
                "gdsc_drug_id": drug,
                "drug_name": DRUG_LABELS[drug],
                "model": model,
                "outer_fits": nouter,
                "zero_molecular_fits": nzero,
                "nonzero_molecular_fits": sum(int(row["nonzero_molecular_fits"]) for row in summary),
                "alpha_ratio_1_selected_fits": alpha_one,
                "warning_outer_fits": sum(int(row["warning_outer_fits"]) for row in summary),
                "fallback_outer_fits": sum(int(row["fallback_outer_fits"]) for row in summary),
                "inner_nonzero_better_contexts": sum(row["nonzero_candidate_strictly_better"] for row in contexts),
                "median_best_nonzero_minus_zero_rmse": statistics.median(diffs),
                "scope": "15_overlapping_outer_training_contexts_descriptive",
            })
            for repeat in range(3):
                subset = [row for row in contexts if row["repeat"] == repeat]
                if len(subset) != 5:
                    raise ValueError(f"Expected 5 folds in repeat {repeat} for {drug}/{model}")
                sr = zero_by_key[(drug, model, repeat)]
                repeat_rows.append({
                    "gdsc_drug_id": drug,
                    "drug_name": DRUG_LABELS[drug],
                    "model": model,
                    "repeat": repeat,
                    "outer_fits": int(sr["outer_fits"]),
                    "zero_molecular_fits": int(sr["zero_molecular_fits"]),
                    "nonzero_molecular_fits": int(sr["nonzero_molecular_fits"]),
                    "alpha_ratio_1_selected_fits": sum(row["zero_molecular_solution"] for row in subset),
                    "nonzero_better_inner_contexts": sum(row["nonzero_candidate_strictly_better"] for row in subset),
                    "warning_outer_fits": int(sr["warning_outer_fits"]),
                    "fallback_outer_fits": int(sr["fallback_outer_fits"]),
                })

    # A sample-level response scan is restricted to the fixed GDSC development file.
    groups: dict[int, list[float]] = defaultdict(list)
    tissues_by_drug: dict[int, dict[str, list[float]]] = {drug: defaultdict(list) for drug in DRUG_ORDER}
    names: dict[int, str] = {}
    with input_paths["development_responses"].open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["role"] != "development" or row["eligibility_status"] == "excluded":
                continue
            drug = int(row["gdsc_drug_id"])
            if drug not in DRUG_LABELS:
                continue
            value = float(row["response_auc"])
            tissue = row["tissue"]
            groups[drug].append(value)
            tissues_by_drug[drug][tissue].append(value)
            names[drug] = row["drug_name"]
    decomp_rows = []
    for drug in DRUG_ORDER:
        values = groups[drug]
        by_tissue = tissues_by_drug[drug]
        if not values or len(by_tissue) != 23:
            raise ValueError(f"Unexpected development response coverage for {drug}")
        grand = statistics.fmean(values)
        total_ss = sum((value - grand) ** 2 for value in values)
        between_ss = sum(len(v) * (statistics.fmean(v) - grand) ** 2 for v in by_tissue.values())
        within_ss = sum(sum((value - statistics.fmean(v)) ** 2 for value in v) for v in by_tissue.values())
        if not math.isclose(total_ss, between_ss + within_ss, rel_tol=1e-11, abs_tol=1e-10):
            raise ValueError(f"Tissue sums of squares do not close for drug {drug}")
        tissue_ns = [len(v) for v in by_tissue.values()]
        decomp_rows.append({
            "gdsc_drug_id": drug,
            "drug_name": DRUG_LABELS[drug],
            "n_observations": len(values),
            "n_biological_groups": len(values),
            "n_tissues": len(by_tissue),
            "response_mean": grand,
            "response_population_sd": statistics.pstdev(values),
            "response_min": min(values),
            "response_max": max(values),
            "total_sum_squares": total_ss,
            "between_tissue_sum_squares": between_ss,
            "within_tissue_sum_squares": within_ss,
            "between_tissue_ss_fraction": between_ss / total_ss if total_ss else None,
            "within_tissue_ss_fraction": within_ss / total_ss if total_ss else None,
            "tissues_n_ge_20": sum(n >= 20 for n in tissue_ns),
            "tissues_n_ge_30": sum(n >= 30 for n in tissue_ns),
            "min_tissue_n": min(tissue_ns),
            "median_tissue_n": statistics.median(tissue_ns),
            "max_tissue_n": max(tissue_ns),
            "scope": "descriptive_fixed_development_sample_sum_of_squares",
        })

    stability_map: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for row in stability_source:
        stability_map[(int(row["gdsc_drug_id"]), row["model"])].append(row)
    deviation_map = {(int(row["gdsc_drug_id"]), row["model"]): row for row in deviation_source}
    zero_count = {(int(row["gdsc_drug_id"]), row["model"]): int(row["zero_molecular_fits"]) for row in zero_source if int(row["repeat"]) == 0}
    profiles = []
    comparison_specs = [
        ("main_delta_rho", "response_association_pair_full", "expression_strong"),
        ("association_pair_minus_tissue_rho", "response_association_pair_full", "tissue_null"),
        ("association_pair_minus_expression200_rho", "response_association_pair_full", "response_association_expression"),
        ("association_pair_minus_rank200_rho", "response_association_pair_full", "response_association_rank"),
        ("association_minus_mad_expression_rho", "response_association_expression", "mad_expression"),
        ("association_minus_mad_rank_rho", "response_association_rank", "mad_rank"),
        ("association_minus_mad_pair_rho", "response_association_pair_full", "mad_pair_full"),
    ]
    for drug in DRUG_ORDER:
        item = {
            "gdsc_drug_id": drug,
            "drug_name": DRUG_LABELS[drug],
            "n_observations": len(groups[drug]),
            "n_tissues": len(tissues_by_drug[drug]),
            "between_tissue_ss_fraction": decomp_rows[DRUG_ORDER.index(drug)]["between_tissue_ss_fraction"],
            "within_tissue_ss_fraction": decomp_rows[DRUG_ORDER.index(drug)]["within_tissue_ss_fraction"],
        }
        for model in ["tissue_null", *MODEL_ORDER]:
            item[f"{model}_spearman"] = mean_metric(metrics, drug, model, "spearman")
            item[f"{model}_rmse"] = mean_metric(metrics, drug, model, "rmse")
        for name, candidate, reference in comparison_specs:
            for metric in ["spearman"]:
                item[name] = item[f"{candidate}_{metric}"] - item[f"{reference}_{metric}"]
        item["strong_expression_minus_association_pair_rmse"] = item["expression_strong_rmse"] - item["response_association_pair_full_rmse"]
        for repeat in range(3):
            item[f"main_delta_rho_repeat_{repeat + 1}"] = (
                mean_metric([r for r in metrics if int(r["repeat"]) == repeat], drug, "response_association_pair_full", "spearman")
                - mean_metric([r for r in metrics if int(r["repeat"]) == repeat], drug, "expression_strong", "spearman")
            )
        for model in MOLECULAR_ORDER:
            item[f"{model}_zero_fits"] = zero_count[(drug, model)]
            st = stability_map[(drug, model)]
            vals = [float(row["active_feature_jaccard"]) for row in st if row["active_feature_jaccard"] != ""]
            item[f"{model}_active_feature_jaccard_median"] = statistics.median(vals) if vals else None
        dev = deviation_map[(drug, "response_association_pair_full")]
        item["association_pair_prediction_rms_deviation_from_tissue"] = fnum(dev["repeat_mean"])
        profiles.append(item)

    # Write immutable-source summaries first; hashes are attached after all products exist.
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    zero_columns = list(zero_rows[0])
    inner_columns = list(inner_rows[0])
    repeat_columns = list(repeat_rows[0])
    decomp_columns = list(decomp_rows[0])
    profile_columns = list(profiles[0])
    write_csv(OUT_DIR / "zero_solution_diagnostics.csv", zero_rows, zero_columns)
    write_csv(OUT_DIR / "inner_path_competition.csv", inner_rows, inner_columns)
    write_csv(OUT_DIR / "zero_solution_by_repeat.csv", repeat_rows, repeat_columns)
    write_csv(OUT_DIR / "response_decomposition.csv", decomp_rows, decomp_columns)
    write_csv(OUT_DIR / "drug_heterogeneity.csv", profiles, profile_columns)

    molecular_total = sum(int(row["outer_fits"]) for row in zero_rows)
    zero_total = sum(int(row["zero_molecular_fits"]) for row in zero_rows)
    alpha1_total = sum(int(row["alpha_ratio_1_selected_fits"]) for row in zero_rows)
    warnings_total = sum(int(row["warning_outer_fits"]) for row in zero_rows)
    fallback_total = sum(int(row["fallback_outer_fits"]) for row in zero_rows)
    if (molecular_total, zero_total, alpha1_total, warnings_total, fallback_total) != (840, 430, 430, 0, 0):
        raise ValueError("Fit-level outcomes disagree with the completed fixed-eight validation totals")
    if [row["n_observations"] for row in decomp_rows] != [625, 622, 620, 620, 620, 620, 620, 620]:
        raise ValueError("Development observation counts disagree with frozen panel coverage")
    expected_dev_observations = int(provenance["oof_rows"]) // (int(provenance["outer_repeats"]) * len(config["models"]))
    if sum(int(row["n_observations"]) for row in decomp_rows) != expected_dev_observations:
        raise ValueError("Development response panel total disagrees with the completed run provenance")

    elapsed = time.perf_counter() - started
    peak_bytes = peak.stop()
    script_hash = sha256(Path(__file__))
    input_hashes = {key: sha256(path) for key, path in input_paths.items()}
    output_paths = [
        OUT_DIR / "zero_solution_diagnostics.csv",
        OUT_DIR / "inner_path_competition.csv",
        OUT_DIR / "zero_solution_by_repeat.csv",
        OUT_DIR / "response_decomposition.csv",
        OUT_DIR / "drug_heterogeneity.csv",
    ]
    manifest = {
        "status": "completed",
        "analysis_name": "fixed_eight_development_zero_solution_and_drug_heterogeneity",
        "analysis_protocol": "docs/DEVELOPMENT_DIAGNOSTICS_PROTOCOL_V1.md",
        "analysis_protocol_version": "v1",
        "analysis_script": "src/analyze_development_diagnostics.py",
        "analysis_script_sha256": script_hash,
        "execution_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, check=True, capture_output=True, text=True
        ).stdout.strip(),
        "source_run": RUN_NAME,
        "source_execution_commit": provenance["execution_commit"],
        "source_lock_sha256": provenance["lock_sha256"],
        "source_metadata_sha256": provenance["metadata_sha256"],
        "source_validation_sha256": provenance["validation_sha256"],
        "source_oof_sha256": provenance["OOF_source_sha256"],
        "development_groups": int(provenance["n_global_biological_groups"]),
        "molecular_outer_fits": molecular_total,
        "zero_molecular_outer_fits": zero_total,
        "alpha_ratio_1_selected_outer_fits": alpha1_total,
        "warning_outer_fits": warnings_total,
        "fallback_outer_fits": fallback_total,
        "external_response_accessed": False,
        "new_model_fits": 0,
        "input_sha256": input_hashes,
        "output_sha256": {},
        "output_rows": {
            "zero_solution_diagnostics.csv": len(zero_rows),
            "inner_path_competition.csv": len(inner_rows),
            "zero_solution_by_repeat.csv": len(repeat_rows),
            "response_decomposition.csv": len(decomp_rows),
            "drug_heterogeneity.csv": len(profiles),
        },
        "elapsed_seconds": elapsed,
        "peak_rss_bytes": peak_bytes,
        "python": sys.version,
        "platform": platform.platform(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "post_hoc_development_descriptive_diagnostic; no new inference family",
    }
    manifest_path = OUT_DIR / "diagnostics_manifest.json"
    manifest["output_sha256"] = {path.name: sha256(path) for path in output_paths}
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report = render_report(profiles, zero_rows, repeat_rows, decomp_rows, inner_rows, manifest)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8", newline="\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=REPO, help="repository root (default: source checkout)")
    args = parser.parse_args()
    if args.root.resolve() != REPO:
        raise SystemExit("For provenance safety, run this script from the repository source location")
    manifest = run()
    print(json.dumps({k: manifest[k] for k in ["status", "execution_commit", "molecular_outer_fits", "zero_molecular_outer_fits", "external_response_accessed", "elapsed_seconds", "peak_rss_bytes", "output_sha256"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
