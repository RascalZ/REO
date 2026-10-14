"""Summarize the frozen 23-tissue LOTO sensitivity run after validation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import run_tissue_loto_v1_3 as loto
from validate_tissue_loto_v1_3 import validate_run


ROOT = loto.ROOT
RESULTS = ROOT / "results/tissue_loto"
REPORT = ROOT / "reports/阶段10_LOTO组织外推敏感性.md"


def spearman(y, prediction):
    y = np.asarray(y, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    if len(y) < 2 or np.ptp(y) == 0:
        return None, "constant_response"
    if np.ptp(prediction) == 0:
        return None, "constant_prediction"
    value = spearmanr(y, prediction).statistic
    if not np.isfinite(value):
        return None, "nonfinite_correlation"
    return float(value), ""


def rmse(y, prediction):
    delta = np.asarray(y, dtype=np.float64) - np.asarray(prediction, dtype=np.float64)
    return float(np.sqrt(np.mean(delta * delta)))


def atomic_csv(path: Path, rows: list[dict], fields: list[str]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(path)
    return loto.sha(path)


def fmt(value, digits=4):
    if value is None or not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):+.{digits}f}"


def build_outputs(run: Path) -> dict:
    lock, lock_sha, _, _ = loto.load_execution_lock()
    config = loto.json_load(loto.EXEC_CONFIG_PATH)
    validation = validate_run(run, config, loto.json_load(loto.MODEL_CONFIG_PATH), lock, lock_sha)
    run_meta = loto.json_load(run / "loto_metadata.json")
    oof = pd.read_csv(
        run / "oof_predictions.csv",
        dtype={"canonical_model_id": str, "biological_group_id": str, "tissue": str,
               "heldout_tissue": str, "drug_name": str, "model": str},
        keep_default_na=False, float_precision="round_trip",
    )
    fit_events = loto.read_jsonl(run / "fit_events.jsonl")
    outer = loto.read_jsonl(run / "outer_fit_audits.jsonl")
    outer_lookup = {
        (int(row["task"]["drug_id"]), row["task"]["heldout_tissue"], row["task"]["model"]): row
        for row in outer
    }
    if len(outer_lookup) != len(outer):
        raise ValueError("Duplicate outer-fit audit mapping")

    metric_rows = []
    raw_by_key = {}
    for (drug_id, drug_name, heldout_tissue, model), subset in oof.groupby(
        ["drug_id", "drug_name", "heldout_tissue", "model"], sort=True
    ):
        y = subset.observed_auc.to_numpy(dtype=np.float64)
        prediction = subset.predicted_auc.to_numpy(dtype=np.float64)
        rho, reason = spearman(y, prediction)
        error = y - prediction
        score_rmse = float(np.sqrt(np.mean(error * error)))
        score_mae = float(np.mean(np.abs(error)))
        score_r2 = None if np.ptp(y) == 0 else float(1.0 - np.sum(error * error) / np.sum((y - y.mean()) ** 2))
        audit = outer_lookup[(int(drug_id), heldout_tissue, model)]
        fit_audit = audit.get("outer_fit_audit", {})
        n_nonzero = fit_audit.get("n_nonzero_molecular_coeffs")
        zero_solution = None if model == "tissue_null" else (None if n_nonzero is None else int(n_nonzero) == 0)
        row = dict(
            drug_id=int(drug_id), drug_name=drug_name, heldout_tissue=heldout_tissue,
            model=model, n_samples=len(subset), spearman=rho, spearman_unavailable_reason=reason,
            rmse=score_rmse, mae=score_mae, r2=score_r2,
            unknown_tissue_samples=int(audit["outer_unknown_tissue_samples"]),
            outer_zero_molecular_solution=zero_solution,
            outer_nonzero_molecular_coefficients=n_nonzero,
            task_wall_seconds=float(audit["resources"]["wall_seconds"]),
            task_peak_rss_bytes=int(audit["resources"]["sampled_peak_rss_bytes"]),
        )
        metric_rows.append(row)
        raw_by_key[(int(drug_id), heldout_tissue, model)] = (y, prediction, row)

    contrast_rows = []
    per_drug_macro = []
    primary_names = loto.json_load(ROOT / "config/tissue_sensitivity_loto_v1.json")["expected_primary_display_tissues"]
    for drug_id in config["drug_ids"]:
        for tissue in config["tissue_order"]:
            strong = raw_by_key[(int(drug_id), tissue, "expression_strong")]
            pair = raw_by_key[(int(drug_id), tissue, "response_association_pair_full")]
            rank = raw_by_key[(int(drug_id), tissue, "response_association_rank")]
            null = raw_by_key[(int(drug_id), tissue, "tissue_null")]
            strong_rho, pair_rho, rank_rho = strong[2]["spearman"], pair[2]["spearman"], rank[2]["spearman"]
            pair_delta = None if strong_rho is None or pair_rho is None else pair_rho - strong_rho
            rank_delta = None if strong_rho is None or rank_rho is None else rank_rho - strong_rho
            row = dict(
                drug_id=int(drug_id), drug_name=strong[2]["drug_name"], heldout_tissue=tissue,
                n_samples=pair[2]["n_samples"],
                pair_minus_expression_strong_delta_rho=pair_delta,
                pair_minus_expression_strong_defined=pair_delta is not None,
                rank_minus_expression_strong_delta_rho=rank_delta,
                rmse_gain_pair_vs_tissue_null=null[2]["rmse"] - pair[2]["rmse"],
                rmse_gain_pair_vs_expression_strong=strong[2]["rmse"] - pair[2]["rmse"],
                pair_rho=pair_rho, expression_strong_rho=strong_rho,
                tissue_null_rho=null[2]["spearman"],
                pair_zero_molecular_solution=pair[2]["outer_zero_molecular_solution"],
            )
            contrast_rows.append(row)

        primary = [row for row in contrast_rows if row["drug_id"] == int(drug_id)
                   and row["heldout_tissue"] in primary_names]
        delta_values = [row["pair_minus_expression_strong_delta_rho"] for row in primary
                        if row["pair_minus_expression_strong_delta_rho"] is not None]
        gain_values = [row["rmse_gain_pair_vs_expression_strong"] for row in primary]
        gain_null_values = [row["rmse_gain_pair_vs_tissue_null"] for row in primary]
        per_drug_macro.append(
            dict(
                unit=str(int(drug_id)), drug_id=int(drug_id),
                drug_name=raw_by_key[(int(drug_id), config["tissue_order"][0], "expression_strong")][2]["drug_name"],
                primary_tissues_expected=10, primary_tissues_observed=len(primary),
                primary_delta_rho_defined=len(delta_values),
                primary_delta_rho_strict_mean=(statistics.fmean(delta_values) if len(primary) == 10 and len(delta_values) == 10 else None),
                primary_delta_rho_defined_cell_mean=(statistics.fmean(delta_values) if delta_values else None),
                primary_rmse_gain_pair_vs_expression_strong_mean=statistics.fmean(gain_values),
                primary_rmse_gain_pair_vs_tissue_null_mean=statistics.fmean(gain_null_values),
                all23_delta_rho_defined=sum(
                    row["pair_minus_expression_strong_delta_rho"] is not None for row in contrast_rows if row["drug_id"] == int(drug_id)
                ),
                all23_delta_rho_defined_cell_mean=(statistics.fmean(
                    row["pair_minus_expression_strong_delta_rho"] for row in contrast_rows
                    if row["drug_id"] == int(drug_id) and row["pair_minus_expression_strong_delta_rho"] is not None
                ) if any(row["drug_id"] == int(drug_id) and row["pair_minus_expression_strong_delta_rho"] is not None for row in contrast_rows) else None),
            )
        )

    macro_values = [row["primary_delta_rho_strict_mean"] for row in per_drug_macro]
    per_drug_macro.append(
        dict(
            unit="equal_weight_macro_8_drugs", drug_id="", drug_name="Equal-weight macro",
            primary_tissues_expected=80, primary_tissues_observed=sum(row["primary_tissues_observed"] for row in per_drug_macro),
            primary_delta_rho_defined=sum(row["primary_delta_rho_defined"] for row in per_drug_macro),
            primary_delta_rho_strict_mean=(statistics.fmean(macro_values) if len(macro_values) == 8 and all(value is not None for value in macro_values) else None),
            primary_delta_rho_defined_cell_mean=None,
            primary_rmse_gain_pair_vs_expression_strong_mean=statistics.fmean(row["primary_rmse_gain_pair_vs_expression_strong_mean"] for row in per_drug_macro),
            primary_rmse_gain_pair_vs_tissue_null_mean=statistics.fmean(row["primary_rmse_gain_pair_vs_tissue_null_mean"] for row in per_drug_macro),
            all23_delta_rho_defined=sum(row["all23_delta_rho_defined"] for row in per_drug_macro),
            all23_delta_rho_defined_cell_mean=None,
        )
    )

    global_rows = []
    for (drug_id, drug_name, model), subset in oof.groupby(["drug_id", "drug_name", "model"], sort=True):
        rho, reason = spearman(subset.observed_auc, subset.predicted_auc)
        global_rows.append(dict(drug_id=int(drug_id), drug_name=drug_name, model=model, n_samples=len(subset),
                                concatenated_cross_tissue_spearman=rho, unavailable_reason=reason,
                                interpretation="cross-tissue calibration-sensitive; not a within-tissue effect"))

    result_hashes = {}
    result_hashes["metrics_by_drug_tissue_model.csv"] = atomic_csv(
        RESULTS / "metrics_by_drug_tissue_model.csv", metric_rows,
        list(metric_rows[0].keys()),
    )
    result_hashes["contrasts_by_drug_tissue.csv"] = atomic_csv(
        RESULTS / "contrasts_by_drug_tissue.csv", contrast_rows,
        list(contrast_rows[0].keys()),
    )
    result_hashes["drug_macro.csv"] = atomic_csv(
        RESULTS / "drug_macro.csv", per_drug_macro,
        list(per_drug_macro[0].keys()),
    )
    result_hashes["concatenated_oof_diagnostics.csv"] = atomic_csv(
        RESULTS / "concatenated_oof_diagnostics.csv", global_rows,
        list(global_rows[0].keys()),
    )
    fit_warning_count = sum(bool(row.get("fit_warning") or row.get("warning")) for row in fit_events)
    fallback_count = sum(bool(row.get("fallback_reason")) for row in fit_events)
    zero_summary = {
        model: sum(row["model"] == model and row["outer_zero_molecular_solution"] is True for row in metric_rows)
        for model in config["models"] if model != "tissue_null"
    }
    split_warning_counts = run_meta.get("preflight", {}).get("inner_splitter_warning_message_counts", {})
    split_warning_text = "; ".join(f"{count}×“{message}”" for message, count in split_warning_counts.items()) or "none"
    main_row = per_drug_macro[-1]
    lines = [
        "# 阶段10：固定八药23组织留一组织外层敏感性分析",
        "",
        f"执行锁：`{lock_sha}`；LOTO执行commit：`{run_meta['execution_commit']}`。",
        f"有效任务/拟合/外层预测：{validation['tasks']:,}/{validation['fits']:,}/{validation['oof_predictions']:,}；新模型运行耗时{sum(session.get('wall_seconds', 0) for session in run_meta['sessions'])/3600:.2f}小时。",
        "此分析属于开发集敏感性检查，不新增确认性检验、P值或置信区间；不访问CTRP。所有预测按药物×完整留出组织评估，不能把组织fold当作来自组织总体的独立重复。",
        "",
        "## 主要对照",
        "",
        "主对照为每个留出组织内的 `response_association_pair_full − expression_strong` Spearman差值；主展示固定为预先锁定的10个每药均≥20独立组的组织。药物内10组织等权，八药再等权。任何主组织的相关未定义时，该药的strict均值记NA；任何一药strict值NA时八药macro记NA。可用格均值只作描述，不能替代strict值。",
        "",
        "|药物|Δρ strict（10组织）|定义组织格|pair RMSE改善 vs 强表达|pair RMSE改善 vs 组织null|",
        "|---|---:|---:|---:|---:|",
    ]
    for row in per_drug_macro[:-1]:
        lines.append(
            f"|{row['drug_name']}|{fmt(row['primary_delta_rho_strict_mean'])}|"
            f"{row['primary_delta_rho_defined']}/10|"
            f"{fmt(row['primary_rmse_gain_pair_vs_expression_strong_mean'])}|"
            f"{fmt(row['primary_rmse_gain_pair_vs_tissue_null_mean'])}|"
        )
    lines.append(
        f"|八药等权宏值|{fmt(main_row['primary_delta_rho_strict_mean'])}|"
        f"{main_row['primary_delta_rho_defined']}/80|"
        f"{fmt(main_row['primary_rmse_gain_pair_vs_expression_strong_mean'])}|"
        f"{fmt(main_row['primary_rmse_gain_pair_vs_tissue_null_mean'])}|"
    )
    lines += [
        "",
        "组织null对所有未见组织均使用外层训练响应的全局均值；在单个留出组织内是常数，因此其Spearman按协议记NA，不能将NA填成0。RMSE改善按同一外层测试样本配对计算，正数代表pair误差更小。",
        "",
        "## 模型零解与资源",
        "",
        f"无拟合预检记录到{run_meta.get('preflight', {}).get('inner_splitter_warning_contexts', 0)}/184个内层分层上下文有稀有组织警告（{split_warning_text}）。三折继续按锁定策略保留，所有外层训练样本均恰好作为一次内层验证样本；罕见组织的调参证据仍不稳定。",
        "外层零分子系数解数量：" + ", ".join(f"`{key}` {value}/184" for key, value in zero_summary.items()) + ".",
        f"内层/外层拟合事件中记录的警告数为{fit_warning_count}，回退数为{fallback_count}；零解、常数预测、NA及所有模型结果均保留。",
        f"控制器及worker进程树观测峰值RSS {max(session.get('sampled_process_tree_peak_rss_bytes', 0) for session in run_meta['sessions']) / 1024**3:.3f} GiB；worker单任务峰值RSS {max(row['task_peak_rss_bytes'] for row in metric_rows) / 1024**3:.3f} GiB。系统采样下限：可用内存 {min(session.get('minimum_available_ram_bytes', float('inf')) for session in run_meta['sessions']) / 1024**3:.3f} GiB，commit余量 {min(session.get('minimum_commit_free_bytes', float('inf')) for session in run_meta['sessions']) / 1024**3:.3f} GiB。采样是操作记录，不是隔离基准或严格内存上界。",
        "",
        "## 解释边界",
        "",
        "LOTO评估模型对训练中未出现组织的迁移压力，不等价于独立实验室/平台验证，也不能据此证明临床可转化。留出组织的标签只用于最后计分；内层调参在其余组织内按固定的训练内分层grouped三折完成。不同组织的细胞系组成、组织标签和药敏测量仍可能构成域偏移。",
        "",
        "全23组织的药物×组织×模型指标、主/次对照、定义格、RMSE改善和跨组织汇总见下列聚合文件。样本级ID、表达/响应、逐样本预测和模型检查点仅保留在本地忽略的`runs/`目录。",
        "",
        "- `../results/tissue_loto/metrics_by_drug_tissue_model.csv`：23组织×八模型逐药汇总指标、覆盖和零解。",
        "- `../results/tissue_loto/contrasts_by_drug_tissue.csv`：每药×组织主要Spearman差值与配对RMSE改善。",
        "- `../results/tissue_loto/drug_macro.csv`：10组织严格药物均值及八药等权宏值，保留NA/定义格。",
        "- `../results/tissue_loto/concatenated_oof_diagnostics.csv`：拼接预测相关，仅作跨组织校准敏感诊断。",
        "- `../results/tissue_loto/loto_analysis_manifest.json`：输入/输出哈希、验证状态及访问声明。",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    result_hashes["reports/阶段10_LOTO组织外推敏感性.md"] = loto.sha(REPORT)
    manifest = {
        "schema_version": "1.0", "status": "completed_descriptive_loto_sensitivity",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "execution_lock_sha256": lock_sha, "execution_commit": run_meta["execution_commit"],
        "run_metadata_sha256": loto.sha(run / "loto_metadata.json"),
        "preflight_manifest_sha256": loto.sha(run / "preflight_manifest.json"),
        "validation": validation,
        "counts": run_meta["output_counts"],
        "output_hashes": result_hashes,
        "external_response_accessed": False, "ctrp_accessed": False,
        "inferential_tests_performed": False, "confidence_intervals_reported": False,
        "interpretation_limit": "Development-set leave-one-tissue-out sensitivity analysis; no independent platform validation or clinical inference",
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "loto_analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    print(json.dumps({"status": manifest["status"], "counts": manifest["counts"],
                      "primary_macro_strict_delta_rho": main_row["primary_delta_rho_strict_mean"],
                      "primary_macro_defined_cells": main_row["primary_delta_rho_defined"],
                      "output_hashes": result_hashes}, ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path)
    args = parser.parse_args()
    config = loto.json_load(loto.EXEC_CONFIG_PATH)
    run = (args.run_dir or ROOT / "runs" / config["run_name"]).resolve()
    run.relative_to(ROOT / "runs")
    build_outputs(run)


if __name__ == "__main__":
    main()
