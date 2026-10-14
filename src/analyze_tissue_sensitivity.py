"""Describe tissue-stratified molecular increments in the frozen eight-drug OOF."""
from __future__ import annotations

import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import hashlib
import json
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"
OUTPUT = ROOT / "results/tissue_sensitivity"
REPORT = ROOT / "reports/阶段8_组织分层OOF诊断.md"
PROTOCOL_LOCK = ROOT / "config/p5e_protocol_v1.lock.json"
ANALYSIS_LOCK = ROOT / "config/tissue_oof_analysis_v1.lock.json"
DRUGS = [1931, 1036, 1191, 1930, 1199, 1091, 1010, 1922]
MODELS = [
    "tissue_null", "expression_strong", "mad_expression", "mad_rank", "mad_pair_full",
    "response_association_expression", "response_association_rank", "response_association_pair_full",
]
PRIMARY_TISSUES = [
    "Bowel", "Breast", "CNS/Brain", "Esophagus/Stomach", "Lung", "Lymphoid", "Myeloid",
    "Ovary/Fallopian Tube", "Pancreas", "Skin",
]
TISSUE_SETS = {"primary10": PRIMARY_TISSUES}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def spearman(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) < 2 or len(y) != len(x):
        return None
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Nonfinite values reached tissue OOF scoring")
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    value = float(spearmanr(x, y).statistic)
    return value if np.isfinite(value) else None


def rmse(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(y - p))))


def verify_source() -> tuple[dict, dict]:
    lock = json.loads(PROTOCOL_LOCK.read_text(encoding="utf-8"))
    if lock["status"] != "frozen_before_p5e_derived_analysis_and_model_fits":
        raise ValueError("P5e protocol is not in its pre-analysis frozen state")
    for relative, expected in lock["inputs"].items():
        if sha(ROOT / relative) != expected:
            raise ValueError("P5e frozen input has changed: " + relative)
    if sha(PROTOCOL_LOCK) != "deb15cd16258b85c86f1605a3099b85319b33a69f086c8ffe912ba12543e0e1b":
        raise ValueError("Unexpected P5e protocol lock")
    metadata = json.loads((RUN / "ordinal_repeated_metadata.json").read_text(encoding="utf-8"))
    if metadata["external_response_accessed"] is not False or metadata["status"] != "completed":
        raise ValueError("Source OOF run is incomplete or external data was accessed")
    if sha(RUN / "oof_predictions.csv") != lock["source_oof_sha256"]:
        raise ValueError("Frozen repeated OOF source hash mismatch")
    return lock, metadata


def aggregate_metrics(metrics: pd.DataFrame, tissue_set: str, tissues: list[str]) -> list[dict]:
    rows: list[dict] = []
    for drug in DRUGS:
        for repeat in range(3):
            d = metrics.loc[
                metrics.drug_id.eq(drug) & metrics.repeat.eq(repeat) & metrics.tissue.isin(tissues)
            ]
            for model in MODELS:
                part = d.loc[d.model.eq(model)].copy()
                expected = len(tissues)
                if len(part) != expected:
                    raise ValueError(f"Missing tissue cells for {drug}/{repeat}/{model}/{tissue_set}")
                row = {
                    "drug_id": drug, "drug_name": part.drug_name.iloc[0], "repeat": repeat,
                    "tissue_set": tissue_set, "model": model,
                    "expected_tissues": expected, "rho_y_pred_defined_tissues": int(part.rho_y_pred.notna().sum()),
                    "rho_y_pred_available_mean": part.rho_y_pred.mean(),
                    "rho_y_pred_strict_mean": part.rho_y_pred.mean() if part.rho_y_pred.notna().all() else np.nan,
                    "rho_residual_defined_tissues": int(part.rho_residual.notna().sum()),
                    "rho_residual_available_mean": part.rho_residual.mean(),
                    "rho_residual_strict_mean": part.rho_residual.mean() if part.rho_residual.notna().all() else np.nan,
                    "rmse_gain_vs_tissue_null": part.rmse_gain_vs_tissue_null.mean(),
                    "rmse_gain_vs_expression_strong": part.rmse_gain_vs_expression_strong.mean(),
                    "delta_rho_vs_expression_strong_defined_tissues": int(part.delta_rho_vs_expression_strong.notna().sum()),
                    "delta_rho_vs_expression_strong_available_mean": part.delta_rho_vs_expression_strong.mean(),
                    "delta_rho_vs_expression_strong_strict_mean": (
                        part.delta_rho_vs_expression_strong.mean()
                        if part.delta_rho_vs_expression_strong.notna().all() else np.nan
                    ),
                }
                rows.append(row)
    return rows


def build_report(metrics: pd.DataFrame, macros: pd.DataFrame, source_sha: str, script_sha: str,
                 elapsed: float, peak_rss: int, coverage: dict) -> str:
    lines = [
        "# 阶段8：固定八药OOF组织分层诊断",
        "",
        "日期：2026-10-14。此分析复用固定八药5折×3重复OOF，不新增拟合，也不访问CTRP。它是开发集组织分层描述，不构成新的确认性检验，不增加置信区间，也不能证明泛化或无泄漏。",
        "",
        f"输入OOF SHA-256：`{source_sha}`；分析程序SHA-256：`{script_sha}`。固定八药共{coverage['oof_rows']:,}条预测、{coverage['n_groups']}个跨药物生物组、{coverage['n_tissues']}个组织；覆盖完整。资源记录：{elapsed:.2f}秒，采样峰值RSS {peak_rss / 1024**3:.3f} GiB。",
        "",
        "## 指标解释",
        "",
        "每药×组织×repeat×模型均报告绝对OOF Spearman、RMSE、`rho_residual = Spearman(y - tissue_null, prediction - tissue_null)`及RMSE相对同政策组织参照的改善量。组织null预测在外折之间可因训练组织均值不同而形成分块；其绝对Spearman只作诊断，不能解释为组织内个体排序能力。组织内残差Spearman若任一输入为常数则记NA；不赋0、不删除该组织。主展示固定为每药均≥20个独立组的10个组织，其余13个组织仍在完整明细表中。",
        "",
        "模型对照Δρ为同一药物/组织/repeat内的`Spearman(y, candidate)-Spearman(y, expression_strong)`，RMSE改善定义为`RMSE(reference)-RMSE(candidate)`，正值有利候选。组织等权宏值同时给完整覆盖的strict值和定义格可用均值；不完整strict值为NA。repeat独立列示，不把3次重复合并为更大的样本量。",
        "",
        "## 主展示组织上的八药汇总",
        "",
        "表内Δρ、ΔRMSE为10个组织等权均值；Δρ strict要求10/10组织格定义。组织null模型的组织内残差Spearman按常数预测记NA。以下为描述性汇总，无P值或新增CI。",
        "",
        "|药物|关联候选完整pair−强表达 Δρ|定义组织格/10|pair组织残差ρ（可用组织均值）|pair RMSE改善 vs null|pair RMSE改善 vs 强表达|",
        "|---|---:|---:|---:|---:|---:|",
    ]
    primary = macros.loc[macros.tissue_set.eq("primary10")]
    for drug in DRUGS:
        part = primary.loc[
            primary.drug_id.eq(drug) & primary.model.eq("response_association_pair_full")
        ]
        # Average repeat-level tissue summaries only when all 3 repeat values exist.
        delta = part.delta_rho_vs_expression_strong_strict_mean
        strict_delta = delta.mean() if delta.notna().all() else np.nan
        n_delta = int(part.delta_rho_vs_expression_strong_defined_tissues.min()) if len(part) else 0
        residual = part.rho_residual_available_mean.mean()
        gain_null = part.rmse_gain_vs_tissue_null.mean()
        gain_strong = part.rmse_gain_vs_expression_strong.mean()
        name = part.drug_name.iloc[0] if len(part) else str(drug)
        fmt = lambda x: "NA" if pd.isna(x) else f"{x:+.4f}"
        lines.append(f"|{name}|{fmt(strict_delta)}|{n_delta}/10|{fmt(residual)}|{fmt(gain_null)}|{fmt(gain_strong)}|")
    lines += [
        "",
        "## 全量产物",
        "",
        "- `results/tissue_sensitivity/oof_tissue_metrics.csv`：药物×组织×repeat×模型聚合指标，含全部23组织。",
        "- `results/tissue_sensitivity/tissue_macro_by_repeat.csv`：每药×repeat×模型的10组织等权摘要，未定义格与严格NA规则均保留。",
        "- `results/tissue_sensitivity/tissue_sensitivity_manifest.json`：输入/输出哈希、覆盖、运行时间、内存及数据访问声明。",
        "",
        "局限：重复外折之间高度重叠；组织样本量不均衡且小组织相关极不稳定；残差相关与RMSE是不同估计目标；该分层汇总未纳入训练折/特征筛选不确定性。LOTO仍是下一步的真正整组织压力分析，不能由此处结果替代。",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    start = time.perf_counter()
    process = psutil.Process()
    cpu0 = process.cpu_times()
    rss_values = [process.memory_info().rss]
    stop = threading.Event()

    def sample_rss() -> None:
        while not stop.wait(0.05):
            try:
                rss_values.append(process.memory_info().rss)
            except psutil.Error:
                pass

    sampler = threading.Thread(target=sample_rss, daemon=True)
    sampler.start()
    lock, source_metadata = verify_source()
    if not ANALYSIS_LOCK.exists():
        raise FileNotFoundError("Frozen tissue OOF analysis lock is required before scoring")
    analysis_lock = json.loads(ANALYSIS_LOCK.read_text(encoding="utf-8"))
    if analysis_lock["analysis_script_sha256"] != sha(Path(__file__)):
        raise ValueError("Frozen tissue OOF analysis script hash mismatch")
    if analysis_lock["protocol_lock_sha256"] != sha(PROTOCOL_LOCK):
        raise ValueError("Tissue OOF analysis points at a different protocol lock")

    source_path = RUN / "oof_predictions.csv"
    frame = pd.read_csv(source_path, dtype={
        "drug_name": str, "model": str, "canonical_model_id": str,
        "biological_group_id": str, "tissue": str, "selection_rule": str,
        "representation": str,
    }, keep_default_na=False, float_precision="round_trip")
    expected_rows = 4967 * len(MODELS) * 3
    if len(frame) != expected_rows:
        raise ValueError(f"Expected {expected_rows} OOF rows, observed {len(frame)}")
    if sorted(frame.drug_id.unique().tolist()) != sorted(DRUGS):
        raise ValueError("Fixed drug panel mismatch")
    if sorted(frame.model.unique().tolist()) != sorted(MODELS):
        raise ValueError("Fixed model panel mismatch")
    if sorted(frame["repeat"].unique().tolist()) != [0, 1, 2] or set(frame.permutation_id) != {0}:
        raise ValueError("Unexpected repeat or permutation coverage")
    if not np.isfinite(frame[["observed_auc", "predicted_auc"]].to_numpy(dtype=float)).all():
        raise ValueError("Nonfinite source response or prediction")

    sample_check = frame.groupby(["drug_id", "repeat", "canonical_model_id"], sort=False).agg(
        model_count=("model", "nunique"), y_count=("observed_auc", "nunique"),
        tissue_count=("tissue", "nunique"), group_count=("biological_group_id", "nunique"),
    )
    if len(sample_check) != 4967 * 3 or not sample_check.eq({"model_count": 8, "y_count": 1, "tissue_count": 1, "group_count": 1}).all().all():
        raise ValueError("Model arms do not share complete sample/response/tissue identities")

    tissue_counts = frame.loc[frame.model.eq("tissue_null")].groupby(
        ["drug_id", "repeat", "tissue"], sort=False
    ).biological_group_id.nunique().reset_index(name="n_groups")
    observed_tissues = sorted(tissue_counts.tissue.unique().tolist())
    if len(observed_tissues) != 23 or not set(PRIMARY_TISSUES).issubset(observed_tissues):
        raise ValueError("Frozen tissue inventory changed")
    for drug in DRUGS:
        per_drug = tissue_counts.loc[tissue_counts.drug_id.eq(drug)]
        # Every tissue belongs in all repeats; the locked display set has at least 20 groups.
        if set(per_drug.tissue.unique()) != set(observed_tissues):
            raise ValueError(f"Tissue coverage differs across the frozen panel for drug {drug}")
        for tissue in PRIMARY_TISSUES:
            n_by_repeat = per_drug.loc[per_drug.tissue.eq(tissue), "n_groups"].unique()
            if len(n_by_repeat) != 1 or int(n_by_repeat[0]) < 20:
                raise ValueError(f"Locked display tissue threshold no longer holds: {drug}/{tissue}")

    metric_rows: list[dict] = []
    grouped = frame.groupby(["drug_id", "repeat", "tissue"], sort=True)
    for (drug, repeat, tissue), subset in grouped:
        wide = subset.pivot(index="canonical_model_id", columns="model", values="predicted_auc")
        observed = subset.drop_duplicates("canonical_model_id").set_index("canonical_model_id")
        wide = wide.sort_index()
        observed = observed.loc[wide.index]
        if wide.columns.tolist() != sorted(MODELS) and set(wide.columns) != set(MODELS):
            raise ValueError("An OOF tissue cell is missing a model arm")
        y = observed.observed_auc.to_numpy(dtype=np.float64)
        null_pred = wide["tissue_null"].to_numpy(dtype=np.float64)
        strong_pred = wide["expression_strong"].to_numpy(dtype=np.float64)
        e0 = y - null_pred
        rmse_null = rmse(y, null_pred)
        rmse_strong = rmse(y, strong_pred)
        rho_null = spearman(y, null_pred)
        for model in MODELS:
            pred = wide[model].to_numpy(dtype=np.float64)
            q = pred - null_pred
            model_rmse = rmse(y, pred)
            rho_model = spearman(y, pred)
            rho_resid = spearman(e0, q)
            rho_strong = spearman(y, strong_pred)
            delta_rho = None if rho_model is None or rho_strong is None else rho_model - rho_strong
            metric_rows.append({
                "drug_id": int(drug), "drug_name": observed.drug_name.iloc[0],
                "repeat": int(repeat), "tissue": tissue, "model": model, "n_samples": len(y),
                "rho_y_pred": rho_model, "rho_tissue_null_y_pred_diagnostic": rho_null,
                "rho_residual": rho_resid, "rmse": model_rmse,
                "rmse_tissue_null": rmse_null, "rmse_gain_vs_tissue_null": rmse_null - model_rmse,
                "rmse_expression_strong": rmse_strong,
                "rmse_gain_vs_expression_strong": rmse_strong - model_rmse,
                "delta_rho_vs_expression_strong": delta_rho,
                "prediction_sd": float(np.std(pred, ddof=0)),
                "molecular_component_sd_vs_null": float(np.std(q, ddof=0)),
                "display_tissue_ge20_groups": tissue in PRIMARY_TISSUES,
            })

    metrics = pd.DataFrame(metric_rows)
    expected_metric_rows = len(DRUGS) * 3 * len(observed_tissues) * len(MODELS)
    if len(metrics) != expected_metric_rows:
        raise ValueError("Aggregated tissue metric grid is incomplete")
    macro_rows = []
    for tissue_set, tissues in TISSUE_SETS.items():
        macro_rows.extend(aggregate_metrics(metrics, tissue_set, tissues))
    macros = pd.DataFrame(macro_rows)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    metrics_path = OUTPUT / "oof_tissue_metrics.csv"
    macros_path = OUTPUT / "tissue_macro_by_repeat.csv"
    metrics.to_csv(metrics_path, index=False, encoding="utf-8-sig", float_format="%.12g")
    macros.to_csv(macros_path, index=False, encoding="utf-8-sig", float_format="%.12g")
    stop.set()
    sampler.join()
    elapsed = time.perf_counter() - start
    cpu = process.cpu_times()
    peak_rss = max(rss_values)
    source_sha = sha(source_path)
    script_sha = sha(Path(__file__))
    coverage = {
        "oof_rows": int(len(frame)), "n_groups": int(frame.biological_group_id.nunique()),
        "n_tissues": int(len(observed_tissues)), "metric_rows": int(len(metrics)),
        "macro_rows": int(len(macros)),
        "undefined_rho_y_pred_cells": int(metrics.rho_y_pred.isna().sum()),
        "undefined_rho_residual_cells": int(metrics.rho_residual.isna().sum()),
        "undefined_delta_rho_vs_expression_strong_cells": int(metrics.delta_rho_vs_expression_strong.isna().sum()),
        "primary_display_tissues": PRIMARY_TISSUES,
        "all_other_tissues_retained": sorted(set(observed_tissues) - set(PRIMARY_TISSUES)),
    }
    report_text = build_report(metrics, macros, source_sha, script_sha, elapsed, peak_rss, coverage)
    REPORT.write_text(report_text, encoding="utf-8", newline="\n")
    outputs = {
        str(metrics_path.relative_to(ROOT)).replace("\\", "/"): sha(metrics_path),
        str(macros_path.relative_to(ROOT)).replace("\\", "/"): sha(macros_path),
        str(REPORT.relative_to(ROOT)).replace("\\", "/"): sha(REPORT),
    }
    manifest = {
        "schema_version": "1.0", "status": "completed_descriptive_oof_tissue_stratification",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_lock_sha256": sha(PROTOCOL_LOCK),
        "analysis_lock_sha256": sha(ANALYSIS_LOCK),
        "source_run_freeze_sha256": lock["source_run_freeze_sha256"],
        "source_oof_sha256": source_sha, "analysis_script_sha256": script_sha,
        "external_response_accessed": False, "new_model_fits": 0,
        "wall_seconds": elapsed, "cpu_seconds": cpu.user + cpu.system - cpu0.user - cpu0.system,
        "sampled_peak_rss_bytes": peak_rss,
        "memory_scope": "Single analysis process; RSS sampled every 50ms and is an observed rather than guaranteed upper bound",
        "coverage": coverage, "output_hashes": outputs,
    }
    manifest_path = OUTPUT / "tissue_sensitivity_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": manifest["status"], "metric_rows": len(metrics),
                      "undefined_residual_rho": coverage["undefined_rho_residual_cells"],
                      "wall_seconds": round(elapsed, 3), "sampled_peak_rss_bytes": peak_rss,
                      "source_oof_sha256": source_sha}, ensure_ascii=False))


if __name__ == "__main__":
    main()
