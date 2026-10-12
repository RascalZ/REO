"""Publish aggregate-only REO200/RF pilot artifacts after a completed run.

Reads explicitly permitted run summaries, outer-fit/feature audits, choices and
cohort-v2 counts. Never opens individual outcomes, OOF prediction records,
GDSC/CTRP response sources or expression matrices. A running/failed run cannot
produce a final report. Mutable report artifacts are outside scientific locks.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN = "runs/gdsc2-cancer-v2-mira1-plx4720-nested-reo200-en-rf"
TITLE = "阶段3：癌模型队列v2与REO200/RF开发集试跑"


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def dump(path, value):
    text = json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    # Public JSON contains no machine paths or individual outcome/fit-ID records.
    assert not re.search(r"[A-Za-z]:[\\/]", text)
    assert not any(f'"{key}"' in text for key in ["fit_sample_ids", "observed_auc", "response_auc", "input_files_accessed"])
    Path(path).write_text(text, encoding="utf-8")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def json_read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def repository_path(path):
    path = Path(path)
    resolved = path.resolve() if path.is_absolute() else (ROOT / path).resolve()
    resolved.relative_to(ROOT)
    return resolved


def git_state():
    def run(*args):
        result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
        return result.stdout.strip() if result.returncode == 0 else None
    status = run("status", "--porcelain")
    return {"report_generation_commit": run("rev-parse", "HEAD"),
            "working_tree_status_at_report_generation": "unavailable" if status is None else "dirty" if status else "clean"}


def finite(value):
    return None if pd.isna(value) else float(value)


def fmt(value):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return "是" if value else "否"
    if isinstance(value, (int, np.integer)):
        return str(value)
    if isinstance(value, (float, np.floating)):
        return f"{value:.5g}"
    return str(value).replace("\n", " ")


def table_md(frame):
    cells = [[fmt(v).replace("|", "\\|") for v in row] for row in frame.itertuples(index=False, name=None)]
    header = "| " + " | ".join(map(str, frame.columns)) + " |"
    return "\n".join([header, "| " + " | ".join(["---"] * len(frame.columns)) + " |"] + ["| " + " | ".join(row) + " |" for row in cells])


def inline_html(text):
    escaped = html.escape(text)
    escaped = re.sub(r"\[([^\]]+)\]\((https://[^\s)]+)\)", r'<a href="\2">\1</a>', escaped)
    return escaped


def render_blocks(blocks):
    markdown = []
    rendered = []
    for kind, value in blocks:
        if kind == "h1":
            markdown.append("# " + value); rendered.append("<h1>" + inline_html(value) + "</h1>")
        elif kind == "h2":
            markdown.append("## " + value); rendered.append("<h2>" + inline_html(value) + "</h2>")
        elif kind == "p":
            markdown.append(value); rendered.append("<p>" + inline_html(value) + "</p>")
        elif kind == "table":
            markdown.append(table_md(value))
            display = value.map(fmt)
            rendered.append('<div class="table-scroll">' + display.to_html(index=False, border=0, escape=True) + "</div>")
        elif kind == "list":
            markdown.append("\n".join("- " + v for v in value))
            rendered.append("<ul>" + "".join("<li>" + inline_html(v) + "</li>" for v in value) + "</ul>")
    document = "\n\n".join(markdown) + "\n"
    page = '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>''' + html.escape(TITLE) + '''</title><style>
    body{margin:0;background:#f3f5f7;color:#1b263b;font-family:"Microsoft YaHei","Noto Sans SC",Arial,sans-serif;line-height:1.7}main{max-width:1150px;margin:30px auto;padding:34px;background:#fff;border-radius:12px}h1{font-size:28px;line-height:1.35}h2{margin-top:32px;border-bottom:2px solid #dce6ef;padding-bottom:6px;font-size:21px}p,li{font-size:15px}.table-scroll{overflow:auto;margin:16px 0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px 11px;text-align:left;border-bottom:1px solid #dfe6ed;white-space:nowrap}th{background:#edf3f8;color:#173b57}tr:nth-child(even){background:#f8fafc}a{color:#165c94}ul{padding-left:24px}@media(max-width:700px){main{margin:0;padding:20px;border-radius:0}h1{font-size:24px}}@media print{body{background:#fff}main{margin:0;padding:0}.table-scroll{overflow:visible}table{font-size:10px}th,td{white-space:normal}}
    </style></head><body><main>''' + "\n".join(rendered) + "</main></body></html>\n"
    return document, page


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default=DEFAULT_RUN)
    args = parser.parse_args()
    run = repository_path(args.run_dir)
    meta = json_read(run / "run_manifest.json")
    if meta.get("status") != "completed":
        raise RuntimeError("Final aggregation requires run status=completed; no report was written.")
    assert meta.get("external_response_accessed") is False
    assert meta.get("blocked_external_access_attempts") == [], "External access attempts require investigation before publication."
    config = meta["configuration"]
    assert config["cohort_version"] == "v2"
    summary = pd.read_csv(run / "summary_metrics.csv")
    folds = pd.read_csv(run / "fold_metrics.csv")
    within = pd.read_csv(run / "within_tissue_metrics.csv")
    audits_source = json_read(run / "outer_fit_audits.json")
    choices_source = json_read(run / "chosen_parameters.json")
    panel = pd.read_csv(ROOT / "results/cohort_v2/main_drugs.csv")
    cohort = json_read(ROOT / "results/cohort_v2/cohort_summary.json")
    validation = json_read(ROOT / "reports/cancer_cohort_v2_validation.json")
    run_validation = json_read(run / "pilot_validation.json")
    assert run_validation['status']=='passed'
    lock = json_read(ROOT / "config/cohort_freeze_v2.lock.json")
    assert validation["status"] == "passed"
    assert sha(ROOT / "config/cohort_freeze_v2.json") == lock["freeze_sha256"]
    assert sha(ROOT / lock["validation_report"]) == lock["validation_report_sha256"]
    expected_keys = {(d, m) for d in config["drug_ids"] for m in config["models"]}
    assert set(zip(summary.gdsc_drug_id, summary.model)) == expected_keys
    assert len(summary) == len(expected_keys) and not summary.duplicated(["gdsc_drug_id", "model"]).any()
    assert len(folds) == len(expected_keys) * config["outer_folds"]
    assert not folds.duplicated(["gdsc_drug_id", "model", "outer_fold"]).any()
    expected_n = panel.set_index("gdsc_drug_id").n_development_models.to_dict()
    assert all(int(row.n) == int(expected_n[row.gdsc_drug_id]) for row in summary.itertuples())
    assert len(audits_source) == len(folds) == len(choices_source)
    labels = {"mean": "训练均值", "tissue": "组织Ridge", "expression": f"表达{config['feature']['expression_genes']} EN",
              "expression_matched": f"表达{config['feature']['matched_expression_genes']} EN", "reo": f"REO{config['feature']['reo_genes']} EN",
              "expression_rf": f"表达{config['feature']['expression_genes']} RF", "reo_rf": f"REO{config['feature']['reo_genes']} RF"}
    drug_names = panel.set_index("gdsc_drug_id").gdsc_drug_name.to_dict()
    summary["drug_name"] = summary.gdsc_drug_id.map(drug_names)
    summary["model_label"] = summary.model.map(labels)
    selected_audit_fields = ["gdsc_drug_id", "model", "outer_fold", "fitted_kind", "fallback_reason", "fit_seconds", "predict_seconds", "n_iter", "actual_alpha", "alpha_max", "path_ratio", "null_endpoint", "solver_reason", "n_nonzero_molecular_coeffs", "dual_gap", "n_selected_genes", "n_features", "n_candidate_molecular_features", "n_molecular_features_used", "n_tissue_features", "rf_pair_selection_reason", "selected_gene_ids_sha256", "built_feature_names_sha256", "rf_selected_pair_names_sha256"]
    outer_rows = []
    for event in audits_source:
        row = {key: event.get(key) for key in selected_audit_fields}
        row["warning_count"] = len(event.get("warnings", []))
        row["molecular_component"] = "not_applicable_baseline" if row["model"] in ["mean", "tissue"] else "zero_coefficients_tissue_only" if row["n_nonzero_molecular_coeffs"] == 0 else "no_molecular_features_or_mean_fallback" if row["n_molecular_features_used"] == 0 or row["fitted_kind"] in ["mean", "tissue"] else "molecular_features_retained"
        outer_rows.append(row)
    outer = pd.DataFrame(outer_rows)
    choices = []
    for event in choices_source:
        choices.append({key: event.get(key) for key in ["gdsc_drug_id", "model", "outer_fold", "chosen_params", "best_inner_rmse", "grid_boundary", "grid_boundary_details", "inner_scores"]})
    comparisons = [("primary_EN", "reo", "expression"), ("secondary_RF", "reo_rf", "expression_rf"),
                   ("matched_gene_EN_ablation", "reo", "expression_matched"),
                   ("REO_EN_vs_tissue", "reo", "tissue"), ("expression_EN_vs_tissue", "expression", "tissue")]
    differences = []
    fold_differences = []
    indexed = summary.set_index(["gdsc_drug_id", "model"])
    indexed_folds = folds.set_index(["gdsc_drug_id", "model", "outer_fold"])
    for drug in config["drug_ids"]:
        for comparison, candidate, reference in comparisons:
            a, b = indexed.loc[(drug, candidate)], indexed.loc[(drug, reference)]
            assert int(a.n) == int(b.n)
            rho_a, rho_b = finite(a.spearman), finite(b.spearman)
            differences.append({"gdsc_drug_id": drug, "drug_name": drug_names[drug], "comparison": comparison, "candidate": candidate, "reference": reference, "n": int(a.n),
                "delta_spearman_candidate_minus_reference": None if rho_a is None or rho_b is None else rho_a - rho_b,
                "rmse_improvement_reference_minus_candidate": float(b.rmse - a.rmse),
                "relative_rmse_improvement": None if b.rmse == 0 else float((b.rmse - a.rmse) / b.rmse),
                "interpretation": "Paired-protocol OOF point difference only; no confidence interval, p-value, equivalence or superiority inference.",
                "pairing_check_scope": "equal expected cohort N and shared runner split contract; this aggregate generator does not reopen individual OOF IDs/outcomes"})
            for fold in range(config["outer_folds"]):
                a, b = indexed_folds.loc[(drug, candidate, fold)], indexed_folds.loc[(drug, reference, fold)]
                assert int(a.n) == int(b.n)
                rho_a, rho_b = finite(a.spearman), finite(b.spearman)
                fold_differences.append({"gdsc_drug_id": drug, "comparison": comparison, "outer_fold": fold, "n": int(a.n),
                    "delta_spearman_candidate_minus_reference": None if rho_a is None or rho_b is None else rho_a - rho_b,
                    "rmse_improvement_reference_minus_candidate": float(b.rmse - a.rmse)})
    feature_rows = []
    for path in sorted((run / "features").glob("*.json")):
        feature = json_read(path)
        # Feature names/indices and any selected-ID lists remain in private run audits.
        row = {key: feature.get(key) for key in ["kind", "training_samples", "input_genes", "requested_genes", "selected_genes", "gene_shortfall", "candidate_pairs", "constant_pairs_removed", "high_tie_nonconstant_pairs_removed", "retained_pairs", "output_features"]}
        row["source_audit_name"] = path.name
        row["retained_pair_training_tie_fraction"] = feature.get("retained_pair_training_tie_fraction")
        rf = feature.get("random_forest", {})
        row["rf_molecular_features_used"] = rf.get("molecular_features_used")
        row["rf_n_estimators"] = rf.get("n_estimators")
        pair = rf.get("pair_selection", {})
        row["rf_pair_selection_summary"] = {key: value for key, value in pair.items() if isinstance(value, (str, int, float, bool, type(None))) and key not in ["selected_feature_names", "selected_indices"]}
        feature_rows.append(row)
    resources = meta["resources"]
    public = {
        "schema_version": "1.0", "report_generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_run": run.relative_to(ROOT).as_posix(), "source_run_status": meta["status"],
        "source_run_created_utc": meta.get("created_utc"), "declared_analysis_date": meta.get("declared_analysis_date"),
        "code_execution_commit": meta.get("code_execution_commit"), **git_state(),
        "source_hashes": meta.get("source_hashes", {}), "config_sha256": meta["config_sha256"],
        "cohort_version": config["cohort_version"], "cohort_freeze_sha256": lock["freeze_sha256"],
        "cohort_validation_report_sha256": lock["validation_report_sha256"],
        "actual_run_validation": {key:run_validation[key] for key in ['status','source_revision','counts','within_tissue_response_variation','checks','limitations']},
        "python": meta["python"], "packages": meta["packages"], "logical_cpus": meta["logical_cpus"], "total_ram_bytes": meta["total_ram_bytes"],
        "configuration": config, "resources": resources, "n_fits": meta["n_fits"], "n_oof_predictions": meta["n_oof_predictions"],
        "external_response_accessed": False, "blocked_external_access_attempt_count": len(meta["blocked_external_access_attempts"]),
        "within_tissue_minimum_n": 20, "within_tissue_undefined_correlation_rows": int(within.spearman.isna().sum()),
        "interpretation": "Expanded engineering/development pilot only. Formal eight-drug repeated nested CV, uncertainty/permutation procedures and frozen external validation remain pending.",
        "report_inputs_sha256": [{"path": path.relative_to(ROOT).as_posix(), "sha256": sha(path)} for path in [run / "run_manifest.json", run / "summary_metrics.csv", run / "fold_metrics.csv", run / "outer_fit_audits.json", run / "chosen_parameters.json", run / "within_tissue_metrics.csv", ROOT / "results/cohort_v2/main_drugs.csv", ROOT / "reports/cancer_cohort_v2_validation.json"]],
    }
    output = ROOT / "results/pilot_reo200_rf"
    output.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output / "summary_metrics.csv", index=False, encoding="utf-8-sig")
    folds.to_csv(output / "fold_metrics.csv", index=False, encoding="utf-8-sig")
    outer.to_csv(output / "outer_fit_audits.csv", index=False, encoding="utf-8-sig")
    within.to_csv(output / "within_tissue_metrics.csv", index=False, encoding="utf-8-sig")
    difference_frame = pd.DataFrame(differences)
    difference_frame.to_csv(output / "paired_metric_differences.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(fold_differences).to_csv(output / "paired_fold_metric_differences.csv", index=False, encoding="utf-8-sig")
    dump(output / "chosen_parameters.json", choices)
    dump(output / "feature_audit_summary.json", feature_rows)
    dump(output / "resource_and_provenance_summary.json", public)
    null_rows = outer[outer.model.isin(["reo", "expression", "expression_matched"]) & outer.n_nonzero_molecular_coeffs.eq(0)]
    main_reo_outer = outer[outer.model.eq("reo")]
    main_reo_zero = int(main_reo_outer.n_nonzero_molecular_coeffs.eq(0).sum())
    if main_reo_zero == len(main_reo_outer):
        leading_finding = f"本轮主要REO EN在全部{len(main_reo_outer)}个外层重拟合中选择零分子系数解，预测实质上仅来自训练组织效应。因此，即使其误差比含表达特征的EN略低，也不能将差异称为REO分子信息的增量。辅助RF的点比较应结合下表判断，但不能挽救为主EN已经有效的结论。"
    else:
        leading_finding = f"本轮主要REO EN有{main_reo_zero}/{len(main_reo_outer)}个外层重拟合为零分子系数解。其余外折保留分子项不等于已证实可靠增量；须结合同ID点差、组织内描述和后续正式重复/外验完整判断。"
    by_model = outer[outer.model.isin(["reo", "expression", "expression_matched"])].groupby(["gdsc_drug_id", "model"]).agg(outer_folds=("outer_fold", "size"), zero_molecular_folds=("n_nonzero_molecular_coeffs", lambda values: int(values.eq(0).sum())), max_nonzero=("n_nonzero_molecular_coeffs", "max"), max_dual_gap=("dual_gap", "max"), warnings=("warning_count", "sum")).reset_index()
    choice_counts = pd.DataFrame(choices).groupby(["gdsc_drug_id", "model"]).grid_boundary.agg(["size", "sum"]).reset_index().rename(columns={"size": "外折数", "sum": "候选边界外折数", "model": "模型", "gdsc_drug_id": "药物ID"})
    main_table = summary[["drug_name", "model_label", "n", "spearman", "rmse", "mae", "r2", "wall_seconds", "convergence_warning_fits", "outer_folds_without_molecular_component"]].rename(columns={"drug_name": "药物", "model_label": "模型", "n": "OOF N", "spearman": "Spearman", "rmse": "RMSE", "mae": "MAE", "r2": "R²", "wall_seconds": "流程秒", "convergence_warning_fits": "警告fits", "outer_folds_without_molecular_component": "无分子项外折"})
    difference_table = difference_frame[difference_frame.comparison.isin(["primary_EN", "secondary_RF", "matched_gene_EN_ablation"])][["drug_name", "comparison", "n", "delta_spearman_candidate_minus_reference", "relative_rmse_improvement"]].copy()
    difference_table["relative_rmse_improvement"] *= 100
    difference_table.columns = ["药物", "比较", "N", "ΔSpearman", "相对RMSE改善%"]
    audit_table = outer[outer.model.isin(["reo", "expression", "expression_matched", "reo_rf", "expression_rf"])][["gdsc_drug_id", "model", "outer_fold", "path_ratio", "actual_alpha", "n_nonzero_molecular_coeffs", "n_candidate_molecular_features", "n_molecular_features_used", "dual_gap", "warning_count"]].rename(columns={"gdsc_drug_id": "药物ID", "model": "模型", "outer_fold": "外折", "path_ratio": "alpha比", "actual_alpha": "实际alpha", "n_nonzero_molecular_coeffs": "非零分子系数", "n_candidate_molecular_features": "候选分子特征", "n_molecular_features_used": "模型分子特征", "dual_gap": "dual gap", "warning_count": "收敛警告"})
    cohort_table = panel[["gdsc_drug_name", "n_development_models", "n_external_primary_models"]].rename(columns={"gdsc_drug_name": "固定主分析化合物", "n_development_models": "v2开发N", "n_external_primary_models": "v2主外验N"})
    resource_table = pd.DataFrame([{"总墙钟秒": resources["wall_seconds"], "CPU秒": resources["cpu_seconds"], "采样峰值RSS GiB": resources["sampled_peak_rss_bytes"] / 1024**3, "Windows峰值working set GiB": None if resources.get("windows_peak_working_set_bytes") is None else resources["windows_peak_working_set_bytes"] / 1024**3, "RSS采样间隔秒": resources["sample_interval_seconds"], "实际fits": meta["n_fits"], "BLAS线程": config["resources"]["blas_threads"], "候选并行数": config["resources"]["parallel_jobs"]}])
    within_table = within[within.model.isin(["tissue", "expression", "reo", "expression_rf", "reo_rf"])][["gdsc_drug_id", "model", "tissue", "n", "spearman", "rmse", "spearman_unavailable_reason"]].rename(columns={"gdsc_drug_id": "药物ID", "model": "模型", "tissue": "组织", "n": "N", "spearman": "Spearman", "rmse": "RMSE", "spearman_unavailable_reason": "相关不可用原因"})
    blocks = [("h1", TITLE), ("p", f"报告生成UTC：{public['report_generated_utc']}；运行声明日期：{public['declared_analysis_date']}；运行状态：{meta['status']}。执行代码提交：{public['code_execution_commit']}；报告生成提交：{public['report_generation_commit']}（工作区{public['working_tree_status_at_report_generation']}）。数值均来自完成运行的聚合输出，不读取个体响应。"),
        ("h2", "本阶段完成与解释边界"),
        ("p", f"本轮在癌模型队列v2的两种预设化合物上运行{len(config['models'])}种模型；外层{config['outer_folds']}折×{config['outer_repeats']}次、内层{config['inner_folds']}折，共{meta['n_fits']}次fit。REO候选基因预算为{config['feature']['reo_genes']}，普通表达主基线预算为{config['feature']['expression_genes']}；匹配表达预算为{config['feature']['matched_expression_genes']}。这是扩展工程pilot及内部点估计，不是正式8药5折×3次结果，也不构成可靠增量或等效证明。"),
        ("p", leading_finding),
        ("p", f"CTRP响应访问标志为False，访问防护拦截尝试数为{len(meta['blocked_external_access_attempts'])}；外验标签尚未用于拟合、选参或评价。所有模型使用相同全局生物组划分；聚合报告检查N与预设名单一致，精确逐ID预测配对由独立检查验收，报告脚本不重新打开个体OOF表。"),
        ("h2", "癌模型队列v2"),
        ("p", f"资格修订仅依据metadata：排除{cohort['removed_development_unique_models']}个开发及{cohort['removed_external_primary_unique_models']}个主外验明确非癌模型，得到{cohort['development_unique_models']}个开发、{cohort['external_unique_models']}个主外验模型。HCC1588/ACH-001078由[Cellosaurus CVCL_A351](https://www.cellosaurus.org/CVCL_A351)及[KCLB71588](https://cellbank.snu.ac.kr/cellline/search?q=71588&sc=y)核实为肺鳞癌模型，精确Sanger标识闭合后纳入。资格修订发生在旧开发pilot之后、CTRP标签开放之前，不能称为旧pilot之前的预注册。v1与旧运行保留，固定8药排序、ID和剂量窗保持。"),
        ("table", cohort_table),
        ("p", f"v2独立资格/身份验证状态为{validation['status']}；主外验与全部已知GDSC2候选生物组交集为{cohort['external_overlap_with_all_known_GDSC_candidate_groups']}。只按权威精确casefold词表修正{cohort['tissue_capitalization_normalized_global_rows']}处组织大小写，保留v1_tissue。队列freeze SHA256：{lock['freeze_sha256']}。"),
        ("h2", "内部OOF模型比较"), ("table", main_table),
        ("p", "主要比较为REO EN相对表达EN2000；RF比较为辅助、同模型家族内的整条pipeline比较。每药单独报告OOF Spearman及原始AUC尺度误差；不合并药物计算一个总相关。训练均值的Spearman按基线定义记NA；不同折均值拼接后产生的排序不能作为个体预测能力。"),
        ("table", difference_table),
        ("p", "ΔSpearman=候选REO−参考表达，正值表示名义方向有利；相对RMSE改善=(参考RMSE−候选RMSE)/参考RMSE。匹配基因数表达模型是消融，不替换2000基因强表达主基线。此处只有点差，无置信区间或确证性P值；3个相关外折不能按3个独立重复作普通t检验，非显著也不证明等效。"),
        ("h2", "正则路径、零分子解与组织贡献"),
        ("p", f"每个EN fit独立从自身训练标准化及组织投影后的X/y计算alpha_max，再使用固定相对路径{config['elastic_net']['alpha_ratios']}和L1比例{config['elastic_net']['l1_ratios']}。训练时组织项采用不受惩罚的FWL投影；REO先由原TPM构造，不能用gene-wise组织残差生成顺序。本轮三种EN的{len(outer[outer.model.isin(['reo', 'expression', 'expression_matched'])])}个外层重拟合中有{len(null_rows)}个零分子系数解。"),
        ("table", by_model), ("table", audit_table),
        ("p", "若某个EN外折的非零分子系数为0，其预测仅来自组织效应，不能把全局Spearman解释为转录组或基因对信号。ratio=1是理论alpha_max零解端点：实现显式采用解析零分子系数，防止浮点epsilon污染组织内排序；该修复只适用于理论端点，没有把较弱正则路径的真实非零系数任意截断。旧pilot与历史文件保持原样。"),
        ("table", choice_counts),
        ("p", "相对路径最大端点被选中可以表示训练内CV偏好零分子项，具有明确的null含义；最弱端点被反复选中则提示搜索下界可能尚不足。候选边界和收敛证据应完整保留，不能据外验结果再拓展路径。同样候选次数并不保证不同表示达到相同求解难度或调参充分度。"),
        ("h2", "RF比较所包含的额外步骤"),
        ("p", f"RF固定{config['random_forest']['n_estimators']}棵树，候选max_features为{config['random_forest']['max_features']}、min_samples_leaf为{config['random_forest']['min_samples_leaf']}。REO RF在每个训练折内按控制组织后的连续响应关联筛至最多{config['random_forest']['reo_max_pairs']}对，普通表达RF没有相应pair监督筛选。因此两种RF的差异包含表示、候选预算和监督筛选等pipeline差异，不能纯归因于编码。两者追加相同训练拟合的组织one-hot；它们的组织政策与EN的FWL不同，跨家族比较也是pipeline比较。"),
        ("h2", "组织内描述与相关不可用情况"),
        ("p", f"仅展示每药组织N≥{public['within_tissue_minimum_n']}的描述性OOF记录，共{len(within)}个模型×组织记录，Spearman为NA的记录数为{public['within_tissue_undefined_correlation_rows']}。组织/均值基线以及所有外折均为零分子项的EN没有组织内个体排序能力，相关记NA并保留原因。稀有组织没有该表的记录不等于被从主分析删除；这些数值也不是独立组织外验或机制证据。"),
        ("table", within_table),
        ("h2", "运行资源、实现修复及可复现记录"), ("table", resource_table),
        ("p", "内存为整个单进程口径，包含输入准备、库、映射页面和拟合；采样峰值RSS及Windows历史峰值working set分列，均不能称为仅模型矩阵内存。流程时间包括内层搜索、外层重拟合和预测。不同模型重复做训练内特征构造是该资源测量的一部分。记录线程限制与并行预算；若线程池信息在限制上下文退出后采集，它反映恢复状态。"),
        ("list", ["理论null端点采用显式解析零系数，已保留solver_reason、actual_alpha、非零数及dual gap证据；不修改旧pilot数字。", f"新维护边界：若合成target仅含组织均值，浮点舍入可引出伪非零EN系数或RF pair筛选。运行结束后已用逐组织精确常量检测修复，38项回归/防泄漏验收通过。实际906fits没有触发该边界；最小全训练集组织内响应残差RMS={run_validation['within_tissue_response_variation']['minimum_positive_fit_residual_rms']:.6f}，远高于舍入残差。独立审计只读已记录OOF/成员，维护修复不改写执行5d05349或数值。", "命令行相对配置路径修复只影响路径解析和可移植执行，不改变科学参数；首个探针在任何拟合前失败并保留失败/资源证据，路径修复后完成8次探针。运行config和源码SHA保存在来源JSON中。", "原异常句柄关闭、失败状态和资源收尾、逐外层检查点修复应继续用于正式长运行；成功pilot尚不能代替完整续跑验收。", "报告只发布汇总指标、outer审计和配置，不公开个体响应或fit成员ID；run目录中的详细内部审计由既定数据管理政策保留。"]),
        ("h2", "下一gate与尚未完成的证据"),
        ("list", ["根据本轮求解、路径边界与资源证据冻结正式模型/统计配置；若调整配置，写明工程或统计理由并保留本运行，不能追求REO胜出。", "按固定8药执行外层5折×3次、内层3折，完整报告主EN、辅助RF、组织/均值基线及匹配表达消融。重复CV改善划分稳定性证据，不增加独立样本N。", "实施预设生物组重采样、配对效应及全流程置换；固定检验族与多重性规则，不把非显著差异称为等效。", "正式开发结束后在全开发拟合并冻结模型/特征/代码，才开放CTRP标签进行一次预设排序迁移评价。不同药敏剂量和共享表达来源的限制仍成立。", "新旧pilot同时改变资格、全局fold、基因预算及alpha路径，不能把分数变化归因于单一因素。REO可能有增量、接近或更差，三种结果均须按预设任务完整报告。"]),
        ("p", "可机读产物：results/pilot_reo200_rf/summary_metrics.csv、fold_metrics.csv、paired_metric_differences.csv、paired_fold_metric_differences.csv、outer_fit_audits.csv、chosen_parameters.json、feature_audit_summary.json、within_tissue_metrics.csv及resource_and_provenance_summary.json。")]
    markdown, page = render_blocks(blocks)
    reports = ROOT / "reports"
    reports.mkdir(exist_ok=True)
    (reports / "阶段3_癌队列v2与REO200_RF试跑.md").write_text(markdown, encoding="utf-8")
    (reports / "阶段3_癌队列v2与REO200_RF试跑.html").write_text(page, encoding="utf-8")
    print(json.dumps({"status": "completed", "source_run": run.relative_to(ROOT).as_posix(), "aggregate_models": len(summary), "outer_audits": len(outer), "reports": ["reports/阶段3_癌队列v2与REO200_RF试跑.md", "reports/阶段3_癌队列v2与REO200_RF试跑.html"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
