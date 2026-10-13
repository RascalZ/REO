"""Post-hoc paired path decomposition of frozen v1 OOF, without any fitting.

All three effects use the identical saved biological-group bootstrap draws.
The path is expression2000 -> expression200 -> REO200. It is an algebraic
performance decomposition, not an identified causal allocation of information.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import psutil
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from analyze_formal_results import (clean, dump, generate_global_bootstrap_weights,
    install_response_source_guard, interval_summary, sha, weighted_metrics)

MODELS = ("expression", "expression_matched", "reo")
METRICS = ("spearman", "rmse")
COMPONENTS = {
    "candidate_range": ("expression_matched", "expression"),
    "representation_and_fitting": ("reo", "expression_matched"),
    "total": ("reo", "expression"),
}


def decompose_arrays(expression2000, expression200, reo200, direction=1):
    """Preserve NA, and check the identity for every finite coordinate."""
    e, m, p = [np.asarray(x, dtype=float) for x in (expression2000, expression200, reo200)]
    if not e.shape == m.shape == p.shape:
        raise ValueError("Decomposition arrays must have identical paired shapes")
    values = {
        "candidate_range": direction * (m - e),
        "representation_and_fitting": direction * (p - m),
        "total": direction * (p - e),
    }
    finite = np.isfinite(e) & np.isfinite(m) & np.isfinite(p)
    error = values["total"] - values["candidate_range"] - values["representation_and_fitting"]
    maximum = float(np.max(np.abs(error[finite]))) if finite.any() else 0.0
    if maximum > 1e-12:
        raise AssertionError("Paired decomposition identity failed")
    return values, {"n_coordinates": int(e.size), "n_finite_coordinates": int(finite.sum()),
                    "n_invalid_coordinates_retained": int((~finite).sum()), "max_abs_identity_error": maximum}


def analyze_decomposition(oof, manifest, config, groups, counts):
    data = oof.rename(columns={"drug_id": "gdsc_drug_id"}).copy()
    groups = np.asarray(groups, dtype=str)
    counts = np.asarray(counts)
    expected_groups = np.asarray(sorted(set(manifest.biological_group_id.astype(str))), dtype=str)
    if not np.array_equal(groups, expected_groups):
        raise ValueError("Bootstrap group order must equal frozen global manifest order")
    if counts.shape != (config["bootstrap"]["draws"], len(groups)):
        raise ValueError("Wrong shared bootstrap shape")
    if not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all() or not (counts.sum(1) == len(groups)).all():
        raise ValueError("Invalid shared bootstrap multiplicities")
    repeats = list(range(config["outer_repeats"]))
    drugs = config["drug_ids_in_frozen_rank"]
    if data.duplicated(["gdsc_drug_id", "model", "repeat", "canonical_model_id"]).any():
        raise ValueError("Duplicate OOF identity")
    position = {g: i for i, g in enumerate(groups)}
    points, boot = {}, {}
    for drug in drugs:
        source = manifest[manifest.gdsc_drug_id.eq(drug)].sort_values("canonical_model_id")
        if source.empty or source.canonical_model_id.duplicated().any():
            raise ValueError("Expected unique per-drug manifest identities")
        ids = source.canonical_model_id.tolist()
        group_map = source.set_index("canonical_model_id").biological_group_id
        group_positions = np.array([position[g] for g in source.biological_group_id])
        subset = data[data.gdsc_drug_id.eq(drug) & data.model.isin(MODELS)]
        if subset.empty or not subset.groupby("canonical_model_id").observed_auc.nunique(dropna=False).eq(1).all():
            raise ValueError("Outcome must be constant across model/repeat records")
        reference_y = subset.groupby("canonical_model_id").observed_auc.first().reindex(ids).to_numpy(float)
        for model in MODELS:
            point_array = np.full((len(repeats), len(METRICS)), np.nan)
            boot_array = np.full((len(repeats), len(counts), len(METRICS)), np.nan)
            for repeat in repeats:
                record = subset[subset.model.eq(model) & subset["repeat"].eq(repeat)].set_index("canonical_model_id")
                if set(record.index) != set(ids):
                    raise ValueError("Incomplete or extra paired OOF identity; do not drop it")
                record = record.reindex(ids)
                if not record.biological_group_id.eq(group_map.reindex(ids)).all():
                    raise ValueError("OOF group identity mismatch")
                if not np.array_equal(record.observed_auc.to_numpy(float), reference_y):
                    raise ValueError("OOF outcomes changed")
                p = record.predicted_auc.to_numpy(float)
                values = weighted_metrics(reference_y, p, np.ones((1, len(ids)), dtype=int))
                point_array[repeat] = [values[metric][0] for metric in METRICS]
                block = config["bootstrap"]["computational_block_draws"]
                for start in range(0, len(counts), block):
                    stop = min(start + block, len(counts))
                    values = weighted_metrics(reference_y, p, counts[start:stop, group_positions])
                    for index, metric in enumerate(METRICS):
                        boot_array[repeat, start:stop, index] = values[metric]
            points[(drug, model)] = point_array
            boot[(drug, model)] = boot_array
    effects, stability, audit = [], [], []
    effect_points, effect_boot = {}, {}
    for drug in drugs:
        for index, metric in enumerate(METRICS):
            direction = -1 if metric == "rmse" else 1
            p, pa = decompose_arrays(*[points[(drug, model)][:, index] for model in MODELS], direction)
            b, ba = decompose_arrays(*[boot[(drug, model)][:, :, index] for model in MODELS], direction)
            audit += [{"gdsc_drug_id": drug, "metric": metric, "level": level, **check}
                      for level, check in [("repeat_points", pa), ("every_repeat_and_bootstrap_draw", ba)]]
            for component in COMPONENTS:
                effect_points[(drug, metric, component)] = p[component]
                effect_boot[(drug, metric, component)] = np.mean(b[component], axis=0)
    for metric in METRICS:
        for component in COMPONENTS:
            effect_points[("macro_equal_fixed_drugs", metric, component)] = np.mean(
                np.stack([effect_points[(drug, metric, component)] for drug in drugs]), axis=0)
            effect_boot[("macro_equal_fixed_drugs", metric, component)] = np.mean(
                np.stack([effect_boot[(drug, metric, component)] for drug in drugs]), axis=0)
        for scope in [*drugs, "macro_equal_fixed_drugs"]:
            # Validate identities again after fixed-repeat / fixed-drug averaging.
            for level, table in [("averaged_bootstrap_draws", effect_boot), ("averaged_repeat_points", effect_points)]:
                a, b, total = [table[(scope, metric, c)] for c in COMPONENTS]
                finite = np.isfinite(a) & np.isfinite(b) & np.isfinite(total)
                errors = np.abs(total - a - b)
                maximum = float(errors[finite].max()) if finite.any() else 0.0
                if maximum > 1e-12:
                    raise AssertionError("Averaged decomposition identity failed")
                audit.append({"gdsc_drug_id": scope, "metric": metric, "level": level,
                    "n_coordinates": int(total.size), "n_finite_coordinates": int(finite.sum()),
                    "n_invalid_coordinates_retained": int((~finite).sum()), "max_abs_identity_error": maximum})
            for component, (candidate, reference) in COMPONENTS.items():
                p = effect_points[(scope, metric, component)]
                b = effect_boot[(scope, metric, component)]
                common = {"gdsc_drug_id": scope, "component": component, "metric": metric,
                    "candidate": candidate, "reference": reference,
                    "effect_direction": "candidate_minus_reference" if metric == "spearman" else "reference_minus_candidate",
                    "positive_favors_candidate": True,
                    "analysis_status": "post_hoc_descriptive_not_confirmatory_or_causal"}
                effects.append({**common, **interval_summary(np.mean(p), b, config, reason="undefined_fixed_repeat_metric"),
                    "n_unique_biological_groups": len(groups) if isinstance(scope, str) else manifest[manifest.gdsc_drug_id.eq(scope)].biological_group_id.nunique()})
                stability.append({**common, "repeat_values": json.dumps(clean(p.tolist()), allow_nan=False),
                    "repeat_mean": np.mean(p), "repeat_min": np.min(p), "repeat_max": np.max(p),
                    "repeat_range": np.ptp(p), "repeat_sample_sd": np.std(p, ddof=1),
                    "positive_repeats": int(np.sum(p > 0)), "negative_repeats": int(np.sum(p < 0)),
                    "zero_repeats": int(np.sum(p == 0)), "undefined_repeats": int(np.sum(~np.isfinite(p)))})
    return {"paired_decomposition_intervals": pd.DataFrame(effects),
            "repeat_stability": pd.DataFrame(stability), "identity_audit": pd.DataFrame(audit)}


def write_report(result, manifest, provenance, path):
    rows = result["paired_decomposition_intervals"]
    names = dict(zip(manifest.gdsc_drug_id, manifest.drug_name))
    names["macro_equal_fixed_drugs"] = "固定八药等权宏平均"
    lines = ["# v1 OOF差值分解与配对置信区间", "",
        "状态：事后描述性分析；未重新拟合、未重新调参；正式v1结果保持不变，CTRP继续封存。", "",
        "比较路径为表达2000 → 表达200 → REO200。总差 = 候选范围步差 + 表示与拟合步差，逐次重复、逐个bootstrap抽样以及最终宏平均均核验这一恒等式。", "",
        "Spearman采用后一步减前一步；RMSE采用前一步减后一步，均为正值有利于后一步。RMSE绝对差可加；相对RMSE改善率的分母不同，不纳入加法分解。", "",
        "每个药先计算三个完整OOF重复的指标再求均值。重用正式v1的626个全局生物学组、10,000次共同抽样及原PCG64种子2026101201。各药保留固定样本覆盖，八药等权；重复次数不计入独立样本量。区间从配对差的抽样分布计算，绝不做区间端点相减。", "",
        "以下均为未校正的95%条件bootstrap区间，限定于现有OOF记录；不涵盖重新训练、重新选特征和调参产生的不确定性，不作新增确认性显著性、因果占比、模型等效或普遍REO优势声明。", ""]
    for metric, title in [("spearman", "ΔSpearman"), ("rmse", "RMSE改善量（原药敏AUC单位）")]:
        lines += ["## " + title, "", "| 药物 | 表达200相对表达2000 | REO200相对表达200 | 总差：REO200相对表达2000 |", "|---|---:|---:|---:|"]
        for drug in [*manifest.gdsc_drug_id.drop_duplicates().tolist(), "macro_equal_fixed_drugs"]:
            selected = rows[(rows.gdsc_drug_id.eq(drug)) & rows.metric.eq(metric)].set_index("component")
            cells = []
            for component in COMPONENTS:
                row = selected.loc[component]
                cells.append(f"{row.point:+.6f} [{row.conditional_95_lower:+.6f}, {row.conditional_95_upper:+.6f}]")
            lines.append("| " + names[drug] + " | " + " | ".join(cells) + " |")
        lines.append("")
    lines += ["## 解读与边界", "",
        "三个诊断案例的点估计支持不同的工作假设：BMS-536924在两个步骤均有负差；Telomerase Inhibitor IX的负总差主要出现在缩小候选范围一步；Gefitinib在缩小范围时没有负点差，负差主要出现在表示与拟合一步。这是沿既定比较路径的描述，不是对总体机制的因果识别。", "",
        "第二步同时改变连续幅度信息、特征展开、可表达的函数和EN正则化几何，不能称为纯编码信息损失。两步骤还可能相互作用；本分解不估计‘候选范围或编码各造成百分之多少’，也不能替代新增秩臂的受控比较。", "",
        "点差为零且区间为[0,0]时，仅表示这批OOF相同预测的配对差为零，不表示未来训练或新域中模型等效。零分子预测及失败结果不作筛除。", "",
        f"代数审计：{len(result['identity_audit'])}条检查，最大绝对浮点误差 {result['identity_audit'].max_abs_identity_error.max():.3g}。所有10,000个draw在所有估计量上均有效；具体重复值和范围见repeat_stability.csv。", "",
        f"运行耗时 {provenance['wall_seconds']:.2f}秒；进程生命周期峰值工作集 {provenance['process_lifetime_peak_rss_bytes'] / 1024**2:.1f}MiB。这里是只读统计过程，不含v1训练开销。", "",
        "复现：`.venv/Scripts/python.exe src/analyze_oof_decomposition.py`；输出目录已存在时默认停止，避免覆盖先前证据。", "",
        "输入哈希、共同权重哈希、代码哈希和禁止数据访问记录见 `results/oof_decomposition/analysis_provenance.json`。", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default="runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1")
    parser.add_argument("--output-dir", default="results/oof_decomposition")
    args = parser.parse_args()
    run, output = [(ROOT / value).resolve() for value in (args.run_dir, args.output_dir)]
    run.relative_to(ROOT)
    output.relative_to(ROOT)
    if output.exists():
        raise FileExistsError("Use a new output directory; existing analysis is retained")
    blocked = install_response_source_guard()
    started = time.perf_counter()
    config_path = ROOT / "config/formal_statistics_v1.json"
    lock_path = ROOT / "config/formal_analysis_v1.lock.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    for filename, checksum in lock["files"].items():
        if sha(ROOT / filename) != checksum:
            raise RuntimeError("Frozen v1 input changed: " + filename)
    meta = json.loads((run / "formal_metadata.json").read_text(encoding="utf-8"))
    if meta.get("completed") is not True or meta.get("completed_tasks") != 840 or meta.get("external_response_accessed") is not False or meta.get("data_guard", {}).get("blocked_attempts") != 0:
        raise RuntimeError("Complete, sealed formal v1 run required")
    if meta["freeze_sha256"] != sha(lock_path) or meta["statistics_sha256"] != sha(config_path):
        raise RuntimeError("Formal run lock/config changed")
    previous_path = ROOT / "results/formal_development/analysis_provenance.json"
    previous = json.loads(previous_path.read_text(encoding="utf-8"))
    oof_path = run / "oof_predictions.csv"
    manifest_path = ROOT / "results/cohort_v2/development_sample_manifest.csv"
    for expected, path in [(previous["OOF_source_sha256"], oof_path),
                           (meta["output_hashes"]["oof_predictions.csv"], oof_path),
                           (previous["development_manifest_sha256"], manifest_path),
                           (previous["analyzer_sha256"], ROOT / "src/analyze_formal_results.py"),
                           (previous["formal_analysis_lock_sha256"], lock_path)]:
        if sha(path) != expected:
            raise RuntimeError("Bound source checksum mismatch: " + path.name)
    group_path = ROOT / "results/formal_development/bootstrap_group_manifest.csv"
    weight_path = ROOT / "results/formal_development/bootstrap_group_multiplicities.npz"
    group_frame = pd.read_csv(group_path, dtype={"biological_group_id": str})
    if not np.array_equal(group_frame.bootstrap_group_position, np.arange(len(group_frame))):
        raise RuntimeError("Saved group positions changed")
    groups = group_frame.biological_group_id.to_numpy(str)
    with np.load(weight_path, allow_pickle=False) as saved:
        counts = saved["multiplicities"]
    group_hash = hashlib.sha256("\0".join(groups).encode()).hexdigest()
    count_hash = hashlib.sha256(counts.tobytes(order="C")).hexdigest()
    if group_hash != previous["bootstrap_group_order_sha256"] or count_hash != previous["bootstrap_count_matrix_sha256"]:
        raise RuntimeError("Saved global bootstrap weights/order changed")
    regenerated_groups, regenerated = generate_global_bootstrap_weights(groups, config["bootstrap"]["draws"], config["bootstrap"]["seed"])
    if not np.array_equal(regenerated_groups, groups) or not np.array_equal(regenerated, counts):
        raise RuntimeError("Saved weights do not reproduce the frozen generator")
    manifest = pd.read_csv(manifest_path, dtype={"canonical_model_id": str, "biological_group_id": str})
    oof = pd.read_csv(oof_path, dtype={"canonical_model_id": str, "biological_group_id": str})
    with threadpool_limits(limits=1):
        result = analyze_decomposition(oof, manifest, config, groups, counts)
    # Previously reported paired quantities must match; only the candidate-range
    # interval is new. This is an additional audit, not the source of new CIs.
    previous_effects = pd.read_csv(ROOT / "results/formal_development/paired_effect_intervals.csv")
    max_difference = 0.0
    for _, row in result["paired_decomposition_intervals"].iterrows():
        comparison = {"total": "primary_EN", "representation_and_fitting": "matched_gene_EN_ablation"}.get(row.component)
        if comparison:
            previous_row = previous_effects[(previous_effects.gdsc_drug_id.astype(str).eq(str(row.gdsc_drug_id))) & previous_effects.comparison.eq(comparison) & previous_effects.metric.eq(row.metric)].iloc[0]
            for field in ("point", "conditional_95_lower", "conditional_95_upper"):
                difference = abs(float(row[field]) - float(previous_row[field]))
                max_difference = max(max_difference, difference)
                if difference > 1e-12:
                    raise AssertionError("Recomputed paired quantity differs from v1")
    output.mkdir(parents=True, exist_ok=False)
    for name, frame in result.items():
        frame.to_csv(output / (name + ".csv"), index=False, encoding="utf-8-sig", lineterminator="\n")
    memory = psutil.Process().memory_info()
    provenance = {
        "status": "completed", "analysis": "post_hoc_oof_path_decomposition_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(), "wall_seconds": time.perf_counter() - started,
        "process_rss_at_completion_bytes": memory.rss,
        "process_lifetime_peak_rss_bytes": getattr(memory, "peak_wset", memory.rss),
        "analyzer_sha256": sha(Path(__file__)), "frozen_metric_helper_sha256": sha(ROOT / "src/analyze_formal_results.py"),
        "source_run": run.relative_to(ROOT).as_posix(), "formal_analysis_lock_sha256": sha(lock_path),
        "v1_statistics_config_sha256": sha(config_path), "v1_statistics_provenance_sha256": sha(previous_path),
        "OOF_source_sha256": sha(oof_path), "development_manifest_sha256": sha(manifest_path),
        "bootstrap_count_matrix_sha256": count_hash, "bootstrap_group_order_sha256": group_hash,
        "saved_weights_reproduced_exactly": True, "bootstrap_draws": len(counts), "n_global_biological_groups": len(groups),
        "fixed_drug_count": 8, "outer_repeats": 3, "no_repeat_N_inflation": True,
        "maximum_previous_paired_statistic_difference": max_difference,
        "identity_checks": len(result["identity_audit"]),
        "maximum_identity_error": result["identity_audit"].max_abs_identity_error.max(),
        "minimum_valid_bootstrap_draw_fraction": result["paired_decomposition_intervals"].bootstrap_valid_fraction.min(),
        "intervals": "Unadjusted descriptive 95% conditional percentile intervals; not causal allocation or full-train-refit uncertainty",
        "new_model_fits": 0, "external_response_accessed": False, "blocked_source_access_attempts": blocked,
        "numpy_version": np.__version__, "pandas_version": pd.__version__,
        "output_sha256": {name + ".csv": sha(output / (name + ".csv")) for name in result},
    }
    dump(output / "analysis_provenance.json", provenance)
    write_report(result, manifest, provenance, ROOT / "docs/OOF_DECOMPOSITION.md")
    print(json.dumps(clean(provenance), ensure_ascii=False))


if __name__ == "__main__":
    main()
