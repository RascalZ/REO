"""Recompute the prespecified tissue OOF contrasts on null-residual scale."""
from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/tissue_sensitivity"
SOURCE = RESULTS / "oof_tissue_metrics.csv"
PER_DRUG_SOURCE = RESULTS / "tissue_macro_by_repeat.csv"
MANIFEST_SOURCE = RESULTS / "tissue_sensitivity_manifest.json"
CONFIG = ROOT / "config/tissue_oof_residual_increment_v1_1.json"
LOCK = ROOT / "config/tissue_residual_increment_v1_2.lock.json"
OUT_TISSUE = RESULTS / "tissue_residual_contrasts_by_repeat.csv"
OUT_DRUG = RESULTS / "tissue_residual_macro_by_drug.csv"
OUT_FIXED8 = RESULTS / "tissue_residual_fixed8_macro_by_repeat.csv"
OUT_MANIFEST = RESULTS / "tissue_residual_increment_manifest.json"
REPORT = ROOT / "reports/阶段8b_组织残差增量澄清.md"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def verify_inputs() -> tuple[dict, dict, pd.DataFrame]:
    lock = load_json(LOCK)
    config = load_json(CONFIG)
    if lock["status"] != "frozen_before_residual_increment_scoring":
        raise ValueError("Residual-increment scoring is not frozen")
    for relative, expected in lock["inputs"].items():
        if sha(ROOT / relative) != expected:
            raise ValueError("Frozen residual-increment input changed: " + relative)
    if config["external_response_access"] != "forbidden":
        raise ValueError("External response access must remain forbidden")
    data = pd.read_csv(SOURCE)
    if len(data) != 4416 or data.drug_id.nunique() != 8 or data["repeat"].nunique() != 3:
        raise ValueError("Expected frozen 4,416-row stage8 tissue metric grid")
    if data.model.nunique() != 8 or data.tissue.nunique() != 23:
        raise ValueError("Stage8 model/tissue panel changed")
    return lock, config, data


def value(cell: pd.DataFrame, model: str, column: str) -> float:
    row = cell.loc[cell.model.eq(model), column]
    if len(row) != 1:
        raise ValueError(f"Expected exactly one row for {model}/{column}")
    v = row.iloc[0]
    return np.nan if pd.isna(v) else float(v)


def make_report(tissue: pd.DataFrame, drug_macro: pd.DataFrame, fixed8: pd.DataFrame,
                source_sha: str, script_sha: str, elapsed: float) -> str:
    lines = [
        "# 阶段8b：组织参照残差尺度的分子增量",
        "",
        "阶段8的绝对预测ρ仍作为诊断产物保存；其折间组织null波动可能影响同组织Spearman。此处按预先写入的P5e澄清v1.1，重算组织null残差与模型分子分量的相关及预设配对差。它不改任何OOF预测，不重拟合模型，不触及CTRP，也不新增CI/P值。",
        "",
        f"聚合输入SHA-256：`{source_sha}`；脚本SHA-256：`{script_sha}`；运行时间{elapsed:.2f}秒。",
        "",
        "定义：`e0 = y - tissue_null_prediction`；`q = model_prediction - tissue_null_prediction`；组织残差ρ为`Spearman(e0,q)`。对比Δρ为候选残差ρ减参考残差ρ。一个组织×repeat任一臂残差向量为常数/无定义时，该配对Δρ记NA。RMSE改善为同一批OOF行上的`RMSE(reference)-RMSE(candidate)`。主展示为10个预设组织；全部23组织的逐格结果均保留。",
        "",
        "## 八药×十组织严格覆盖结果",
        "",
        "下表按repeat分别列出固定八药等权均值。Δρ strict要求每药10个组织、再要求全部8药可定义；任何缺格均为NA。RMSE无相关系数的常数限制，正常成对计算。可用格均值和分母见CSV，不取代strict值。",
        "",
        "|候选−参考|repeat|strict Δρ|定义药物/8|可定义组织格/80|可用格Δρ（诊断）|RMSE改善|",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    selected = fixed8.loc[fixed8.contrast.eq("response_association_pair_full_minus_expression_strong")].sort_values("repeat")
    for _, row in selected.iterrows():
        strict = "NA" if pd.isna(row.delta_rho_strict_equal_drug_macro) else f"{row.delta_rho_strict_equal_drug_macro:+.4f}"
        available = "NA" if pd.isna(row.delta_rho_available_tissue_cell_mean) else f"{row.delta_rho_available_tissue_cell_mean:+.4f}"
        rmseg = "NA" if pd.isna(row.rmse_gain_strict_equal_drug_macro) else f"{row.rmse_gain_strict_equal_drug_macro:+.4f}"
        lines.append(f"|assoc pair−strong|{int(row.repeat)}|{strict}|{int(row.delta_rho_defined_drugs)}/8|{int(row.delta_rho_defined_tissue_cells)}/80|{available}|{rmseg}|")
    lines += [
        "",
        "## 逐药摘要",
        "",
        "|药物|assoc pair−strong Δρ strict/available|覆盖|pair−同候选表达 Δρ strict/available|覆盖|pair−同候选rank Δρ strict/available|覆盖|pair RMSE改善 vs strong|",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    primary = drug_macro.loc[drug_macro.tissue_set.eq("primary10")]
    contrast_names = [
        "response_association_pair_full_minus_expression_strong",
        "response_association_pair_full_minus_response_association_expression",
        "response_association_pair_full_minus_response_association_rank",
    ]
    for drug in sorted(primary.drug_id.unique()):
        d = primary.loc[primary.drug_id.eq(drug)]
        parts = []
        for contrast in contrast_names:
            c = d.loc[d.contrast.eq(contrast)].sort_values("repeat")
            vals = c.delta_rho_strict_tissue_mean
            strict_value = vals.mean() if vals.notna().all() else np.nan
            available_value = c.delta_rho_available_case_mean.mean()
            denom = int(c.delta_rho_defined_tissues.min()) if len(c) else 0
            strict_text = "NA" if pd.isna(strict_value) else f"{strict_value:+.4f}"
            available_text = "NA" if pd.isna(available_value) else f"{available_value:+.4f}"
            parts.extend([f"{strict_text} / {available_text}", f">={denom}/10 per repeat"])
        c_strong = d.loc[d.contrast.eq(contrast_names[0])]
        gain = c_strong.rmse_gain.mean() if len(c_strong) else np.nan
        gain_text = "NA" if pd.isna(gain) else f"{gain:+.4f}"
        name = d.drug_name.iloc[0]
        lines.append(f"|{name}|{parts[0]}|{parts[1]}|{parts[2]}|{parts[3]}|{parts[4]}|{parts[5]}|{gain_text}|")
    lines += [
        "",
        "## 解读边界",
        "",
        "此处相关是模型相对组织参照分量与组织参照残差的关联，区别于绝对预测rho；它仍条件于已有开发OOF与响应关联候选筛选，不能解释为独立验证、因果机制或无偏显著增量。重复训练集相互重叠，三repeat均值/范围不作为置信区间。stage8/8a的绝对相关对比原样保留，现明确降为包含外折组织基线波动的诊断量。",
        "",
        "实现记录：v1.1首次汇总在写出三个聚合CSV后于报告字段名称处停止，未写manifest/报告，也未访问外部数据或拟合模型；v1.2只修正报告字段绑定和可用格分母汇总，所有模型/相关/误差公式不变。",
        "",
        "机器可读结果：",
        "",
        "- `results/tissue_sensitivity/tissue_residual_contrasts_by_repeat.csv`：药物×组织×repeat四个预设配对对比。",
        "- `results/tissue_sensitivity/tissue_residual_macro_by_drug.csv`：药物内10组织strict与available-case汇总。",
        "- `results/tissue_sensitivity/tissue_residual_fixed8_macro_by_repeat.csv`：固定八药等权macro及可用分母。",
        "- `results/tissue_sensitivity/tissue_residual_increment_manifest.json`：来源/输出哈希、覆盖及禁读记录。",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    start = time.perf_counter()
    lock, config, data = verify_inputs()
    contrasts = config["primary_contrasts"]
    display_tissues = config["tissue_sets"]["primary10"]
    tissue_rows = []
    for (drug, repeat, tissue), cell in data.groupby(["drug_id", "repeat", "tissue"], sort=True):
        for contrast in contrasts:
            candidate, reference = contrast["candidate"], contrast["reference"]
            r_c = value(cell, candidate, "rho_residual")
            r_r = value(cell, reference, "rho_residual")
            delta = r_c - r_r if np.isfinite(r_c) and np.isfinite(r_r) else np.nan
            rmse_c = value(cell, candidate, "rmse")
            rmse_r = value(cell, reference, "rmse")
            rho_gain = value(cell, candidate, "rmse_gain_vs_tissue_null")
            rho_gain_ref = value(cell, reference, "rmse_gain_vs_tissue_null")
            # Pairwise RMSE on identical OOF rows; for the strong reference the source stores it directly.
            rmse_gain = rmse_r - rmse_c
            tissue_rows.append({
                "drug_id": int(drug), "drug_name": cell.drug_name.iloc[0], "repeat": int(repeat),
                "tissue": tissue, "primary_display_tissue": tissue in display_tissues,
                "contrast": contrast["name"], "candidate": candidate, "reference": reference,
                "candidate_residual_rho": r_c, "reference_residual_rho": r_r,
                "delta_rho_residual": delta, "rmse_reference": rmse_r, "rmse_candidate": rmse_c,
                "rmse_gain_reference_minus_candidate": rmse_gain,
                "candidate_gain_vs_tissue_null": rho_gain,
                "reference_gain_vs_tissue_null": rho_gain_ref,
                "n_samples": int(cell.n_samples.iloc[0]),
            })
    tissue = pd.DataFrame(tissue_rows)
    expected_tissue_rows = 8 * 3 * 23 * len(contrasts)
    if len(tissue) != expected_tissue_rows:
        raise ValueError(f"Expected {expected_tissue_rows} tissue contrast rows")

    drug_rows = []
    primary = tissue.loc[tissue.primary_display_tissue]
    for (drug, repeat, contrast), part in primary.groupby(["drug_id", "repeat", "contrast"], sort=True):
        if len(part) != len(display_tissues):
            raise ValueError("Missing a prespecified primary tissue cell")
        defined = part.delta_rho_residual.notna()
        drug_rows.append({
            "drug_id": int(drug), "drug_name": part.drug_name.iloc[0], "repeat": int(repeat),
            "tissue_set": "primary10", "contrast": contrast,
            "expected_tissues": len(display_tissues), "delta_rho_defined_tissues": int(defined.sum()),
            "delta_rho_available_case_mean": part.delta_rho_residual.mean(),
            "delta_rho_strict_tissue_mean": part.delta_rho_residual.mean() if defined.all() else np.nan,
            "rmse_gain": part.rmse_gain_reference_minus_candidate.mean(),
            "n_samples_min": int(part.n_samples.min()), "n_samples_max": int(part.n_samples.max()),
        })
    drug_macro = pd.DataFrame(drug_rows)
    fixed8_rows = []
    for (repeat, contrast), part in drug_macro.groupby(["repeat", "contrast"], sort=True):
        if len(part) != 8:
            raise ValueError("Fixed eight-drug summary is incomplete")
        n_drugs = int(part.delta_rho_strict_tissue_mean.notna().sum())
        defined_cells = int(part.delta_rho_defined_tissues.sum())
        tissue_cells = tissue.loc[
            tissue.repeat.eq(repeat) & tissue.contrast.eq(contrast) & tissue.primary_display_tissue
        ]
        if len(tissue_cells) != 80:
            raise ValueError("Fixed-eight tissue-cell comparison denominator is not 80")
        available_drugs = int(part.delta_rho_available_case_mean.notna().sum())
        fixed8_rows.append({
            "repeat": int(repeat), "tissue_set": "primary10", "contrast": contrast,
            "expected_drugs": 8, "expected_tissue_cells": 80,
            "delta_rho_defined_drugs": n_drugs, "delta_rho_defined_tissue_cells": defined_cells,
            "delta_rho_available_case_equal_drug_mean": part.delta_rho_available_case_mean.mean(),
            "delta_rho_available_drugs": available_drugs,
            "delta_rho_available_tissue_cell_mean": tissue_cells.delta_rho_residual.mean(),
            "delta_rho_strict_equal_drug_macro": (
                part.delta_rho_strict_tissue_mean.mean() if n_drugs == 8 else np.nan
            ),
            "rmse_gain_defined_drugs": int(part.rmse_gain.notna().sum()),
            "rmse_gain_strict_equal_drug_macro": part.rmse_gain.mean() if part.rmse_gain.notna().all() else np.nan,
        })
    fixed8 = pd.DataFrame(fixed8_rows)
    if len(drug_macro) != 8 * 3 * len(contrasts) or len(fixed8) != 3 * len(contrasts):
        raise ValueError("Per-drug/fixed-eight contrast summary has wrong coverage")

    tissue.to_csv(OUT_TISSUE, index=False, encoding="utf-8-sig", float_format="%.12g")
    drug_macro.to_csv(OUT_DRUG, index=False, encoding="utf-8-sig", float_format="%.12g")
    fixed8.to_csv(OUT_FIXED8, index=False, encoding="utf-8-sig", float_format="%.12g")
    elapsed = time.perf_counter() - start
    script_sha = sha(Path(__file__))
    source_sha = sha(SOURCE)
    report = make_report(tissue, drug_macro, fixed8, source_sha, script_sha, elapsed)
    REPORT.write_text(report, encoding="utf-8", newline="\n")
    output_paths = [OUT_TISSUE, OUT_DRUG, OUT_FIXED8, REPORT]
    output_hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in output_paths}
    manifest = {
        "schema_version": "1.0", "status": "completed_residual_increment_summary",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_lock_sha256": lock["protocol_lock_sha256"], "analysis_lock_sha256": sha(LOCK),
        "source_tissue_metrics_sha256": source_sha,
        "source_per_drug_macro_sha256": sha(PER_DRUG_SOURCE),
        "source_stage8_manifest_sha256": sha(MANIFEST_SOURCE),
        "analysis_script_sha256": script_sha, "external_response_accessed": False,
        "new_model_fits": 0,
        "coverage": {"tissue_rows": len(tissue), "drug_macro_rows": len(drug_macro),
                     "fixed8_macro_rows": len(fixed8), "n_drugs": 8, "n_repeats": 3,
                     "n_tissues": 23, "primary_tissues": display_tissues,
                     "undefined_residual_contrasts_primary10": int(primary.delta_rho_residual.isna().sum()),
                     "undefined_fixed8_strict_repeat_contrasts": int(fixed8.delta_rho_strict_equal_drug_macro.isna().sum())},
        "wall_seconds": elapsed, "output_hashes": output_hashes,
        "implementation_attempts": [
            {"version": "v1.1", "outcome": "report generation failed after aggregate CSV writes",
             "reason": "report referenced a nonexistent available-case field",
             "scientific_results_changed": False, "model_fits": 0, "external_access": False},
            {"version": "v1.2", "outcome": "report field corrected; explicit available tissue-cell mean added",
             "scientific_results_changed": False, "model_fits": 0, "external_access": False},
        ],
        "interpretation": "Prespecified residual-scale descriptive contrasts; no inference or model refitting",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": manifest["status"], "tissue_rows": len(tissue),
                      "undefined_primary_contrasts": manifest["coverage"]["undefined_residual_contrasts_primary10"],
                      "fixed8_strict_NA": manifest["coverage"]["undefined_fixed8_strict_repeat_contrasts"],
                      "wall_seconds": round(elapsed, 3), "source_sha256": source_sha}, ensure_ascii=False))


if __name__ == "__main__":
    main()
