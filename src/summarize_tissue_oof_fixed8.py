"""Build equal-drug tissue-stratified macro summaries from the frozen OOF tables."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/tissue_sensitivity"
INPUT = RESULTS / "tissue_macro_by_repeat.csv"
OUTPUT = RESULTS / "fixed8_macro_by_repeat.csv"
REPORT = ROOT / "reports/阶段8a_固定八药组织宏汇总.md"
MANIFEST = RESULTS / "fixed8_macro_manifest.json"
ANALYSIS_LOCK = ROOT / "config/tissue_oof_fixed8_macro_v1.lock.json"
MODELS = [
    "tissue_null", "expression_strong", "mad_expression", "mad_rank", "mad_pair_full",
    "response_association_expression", "response_association_rank", "response_association_pair_full",
]
DRUG_IDS = [1931, 1036, 1191, 1930, 1199, 1091, 1010, 1922]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify() -> tuple[dict, pd.DataFrame]:
    lock = json.loads(ANALYSIS_LOCK.read_text(encoding="utf-8"))
    if lock["status"] != "frozen_before_fixed8_macro_aggregation":
        raise ValueError("Fixed-eight macro aggregation is not frozen")
    for relative, expected in lock["inputs"].items():
        if sha(ROOT / relative) != expected:
            raise ValueError("Frozen fixed-eight macro input changed: " + relative)
    data = pd.read_csv(INPUT)
    if sorted(data.drug_id.unique().tolist()) != sorted(DRUG_IDS):
        raise ValueError("Fixed eight drug identities changed")
    if sorted(data.model.unique().tolist()) != sorted(MODELS) or sorted(data.repeat.unique().tolist()) != [0, 1, 2]:
        raise ValueError("Expected tissue macro rows are incomplete")
    if set(data.tissue_set) != {"primary10"}:
        raise ValueError("Unexpected tissue macro set")
    if len(data) != 8 * 3 * 8:
        raise ValueError("Expected 192 fixed drug/repeat/model summaries")
    return lock, data


def main() -> None:
    start = time.perf_counter()
    lock, data = verify()
    rows = []
    metrics = [
        "rho_y_pred", "rho_residual", "rmse_gain_vs_tissue_null",
        "rmse_gain_vs_expression_strong", "delta_rho_vs_expression_strong",
    ]
    sources = {
        "rho_y_pred": ("rho_y_pred_strict_mean", "rho_y_pred_defined_tissues"),
        "rho_residual": ("rho_residual_strict_mean", "rho_residual_defined_tissues"),
        "rmse_gain_vs_tissue_null": ("rmse_gain_vs_tissue_null", None),
        "rmse_gain_vs_expression_strong": ("rmse_gain_vs_expression_strong", None),
        "delta_rho_vs_expression_strong": (
            "delta_rho_vs_expression_strong_strict_mean",
            "delta_rho_vs_expression_strong_defined_tissues",
        ),
    }
    for repeat in range(3):
        for model in MODELS:
            cell = data.loc[data.repeat.eq(repeat) & data.model.eq(model)]
            if len(cell) != 8:
                raise ValueError(f"Expected all eight drugs for {repeat}/{model}")
            row = {"repeat": repeat, "tissue_set": "primary10", "model": model,
                   "expected_drugs": 8, "expected_tissue_cells": 80}
            for metric in metrics:
                value_col, count_col = sources[metric]
                values = cell[value_col]
                defined_drugs = int(values.notna().sum())
                if count_col:
                    defined_tissue_cells = int(cell[count_col].sum())
                else:
                    defined_tissue_cells = 80
                strict = float(values.mean()) if defined_drugs == 8 else np.nan
                row[metric + "_defined_drugs"] = defined_drugs
                row[metric + "_defined_tissue_cells"] = defined_tissue_cells
                row[metric + "_available_case_equal_drug_mean"] = float(values.mean()) if defined_drugs else np.nan
                row[metric + "_strict_equal_drug_mean"] = strict
            rows.append(row)
    result = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, index=False, encoding="utf-8-sig", float_format="%.12g")

    # Descriptive summaries retain the three repeat-specific values and show their range.
    lines = [
        "# 阶段8a：固定八药组织分层宏汇总",
        "",
        "本补充汇总按已冻结协议要求，在每个repeat先对10个主展示组织等权汇总，再对固定8药等权平均。它只重排阶段8聚合表，不访问样本级新预测、不新增模型拟合或统计推断。三个repeat单独列示；其均值/范围仅描述运行间变化，不视作独立抽样区间。",
        "",
        f"源阶段8宏表SHA-256：`{sha(INPUT)}`；代码SHA-256：`{sha(Path(__file__))}`；本次用时{time.perf_counter()-start:.2f}秒。CTRP访问：否。",
        "",
        "|模型|repeat 0 Δρ|repeat 1 Δρ|repeat 2 Δρ|三repeat均值[范围]|定义药物/8|组织残差ρ可用均值|ΔRMSE vs null|ΔRMSE vs strong|",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for model in MODELS:
        part = result.loc[result.model.eq(model)].sort_values("repeat")
        d = part.delta_rho_vs_expression_strong_strict_equal_drug_mean
        if d.notna().all():
            d_summary = f"{d.mean():+.4f} [{d.min():+.4f},{d.max():+.4f}]"
        else:
            d_summary = "NA [" + ",".join("NA" if pd.isna(x) else f"{x:+.4f}" for x in d) + "]"
        d_values = ["NA" if pd.isna(x) else f"{x:+.4f}" for x in d]
        defined = int(part.delta_rho_vs_expression_strong_defined_drugs.min())
        rho = part.rho_residual_available_case_equal_drug_mean.mean()
        n_rho = int(part.rho_residual_defined_drugs.min())
        gain_null = part.rmse_gain_vs_tissue_null_strict_equal_drug_mean.mean()
        gain_strong = part.rmse_gain_vs_expression_strong_strict_equal_drug_mean.mean()
        fmt = lambda x: "NA" if pd.isna(x) else f"{x:+.4f}"
        lines.append(
            f"|{model}|{d_values[0]}|{d_values[1]}|{d_values[2]}|{d_summary}|{defined}/8|"
            f"{fmt(rho)} ({n_rho}/8)|{fmt(gain_null)}|{fmt(gain_strong)}|"
        )
    lines += [
        "",
        "Δρ定义为各组织内`Spearman(y, candidate)-Spearman(y, expression_strong)`；strict值要求固定8药及10个组织格完整。正的RMSE改善代表相对参照误差较低。组织残差ρ栏为每药可定义组织的均值再作等药权平均，分母另列；它不等同于strict宏值。",
        "",
        "完整repeat×模型×药物数据见`results/tissue_sensitivity/fixed8_macro_by_repeat.csv`。无P值或CI；不能据此宣称REO优势、非劣或等效。",
    ]
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    manifest = {
        "schema_version": "1.0", "status": "completed_fixed8_tissue_macro_summary",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_lock_sha256": sha(ANALYSIS_LOCK), "protocol_lock_sha256": lock["protocol_lock_sha256"],
        "source_macro_sha256": sha(INPUT), "script_sha256": sha(Path(__file__)),
        "external_response_accessed": False, "new_model_fits": 0,
        "wall_seconds": time.perf_counter() - start,
        "output_hashes": {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path)
                          for path in [OUTPUT, REPORT]},
        "interpretation": "Descriptive equal-drug summary over frozen ten-tissue macros; repeats are not independent inferential replicates",
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": manifest["status"], "rows": len(result),
                      "macro_lock_sha256": sha(ANALYSIS_LOCK), "report": str(REPORT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
