"""Aggregate-only report/figures after formal development and QA complete.

Never reads OOF rows, raw response values, CTRP, expression matrices, bootstrap
group identities or multiplicities. This is mutable reporting code, separate
from the locked numerical engine and statistical protocol.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import html
import itertools
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from report_reo200_rf_pilot import fmt, table_md, inline_html

ROOT = Path(__file__).resolve().parents[1]
MACRO = "macro_equal_fixed_drugs"
DRUG_NAMES = {1931: "MIRA-1", 1036: "PLX-4720", 1191: "Bortezomib", 1930: "Telomerase Inhibitor IX",
              1199: "Tamoxifen", 1091: "BMS-536924", 1010: "Gefitinib", 1922: "Cediranib"}
MODEL_NAMES = {"mean": "训练均值", "tissue": "组织Ridge", "expression": "表达2000 EN", "expression_matched": "表达200 EN",
               "reo": "REO200 EN", "expression_rf": "表达2000 RF", "reo_rf": "REO200 RF"}
MODEL_PLOT = {"mean": "Training mean", "tissue": "Tissue Ridge", "expression": "Expression2000 EN", "expression_matched": "Expression200 EN",
              "reo": "REO200 EN", "expression_rf": "Expression2000 RF", "reo_rf": "REO200 RF"}
COLORS = {"tissue": "#555555", "expression": "#0072B2", "expression_matched": "#56B4E9", "reo": "#D55E00", "expression_rf": "#009E73", "reo_rf": "#CC79A7"}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def local(path):
    value = Path(path)
    result = value.resolve() if value.is_absolute() else (ROOT / value).resolve()
    result.relative_to(ROOT)
    return result


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def name(drug):
    return "固定8药等权宏平均" if str(drug) == MACRO else DRUG_NAMES.get(int(drug), str(drug))


def plot_name(drug):
    return "Equal-drug macro (8)" if str(drug) == MACRO else DRUG_NAMES.get(int(drug), str(drug))


def jaccard(a, b):
    union = a | b
    return np.nan if not union else len(a & b) / len(union)


def descriptive_audit_tables(audits):
    """Stream task audits; cache small gene and RF-selected-pair sets only."""
    resources, parameters, summaries = [], [], []
    feature_sets = defaultdict(list)
    for event in audits:
        task = event["task"]
        context = {"drug_id": int(task["drug_id"]), "model": task["model"], "repeat": int(task["repeat"]), "outer_fold": int(task["outer_fold"])}
        if context["model"] not in ["mean", "tissue"]:
            feature_sets.setdefault((context["drug_id"], context["model"], "training_MAD_candidate_genes"), [])
        if context["model"] == "reo_rf":
            feature_sets.setdefault((context["drug_id"], context["model"], "RF_supervised_selected_REO_pairs"), [])
        audit = event["outer_fit_audit"]
        res = event.get("resources", {})
        resources.append({**context, "task_wall_seconds": res.get("wall_seconds"), "task_cpu_seconds": res.get("cpu_seconds"),
            "worker_task_peak_rss_bytes": res.get("sampled_peak_rss_bytes"), "worker_lifetime_peak_working_set_bytes": res.get("windows_lifetime_peak_working_set_bytes")})
        parameters.append({**context, "chosen_parameters": json.dumps(event.get("chosen_params", {}), sort_keys=True),
            "best_inner_rmse": event.get("best_inner_rmse"), "grid_boundary": bool(event.get("chosen_grid_boundary", False)),
            "grid_boundary_details": json.dumps(event.get("chosen_grid_boundary_details", {}), sort_keys=True),
            "actual_alpha": audit.get("actual_alpha"), "alpha_max": audit.get("alpha_max"), "alpha_ratio": audit.get("path_ratio"),
            "null_endpoint": audit.get("null_endpoint"), "solver_reason": audit.get("solver_reason"), "n_iter": audit.get("n_iter"), "dual_gap": audit.get("dual_gap")})
        nzero = audit.get("n_nonzero_molecular_coeffs")
        summaries.append({**context, "zero_molecular_coefficients": int(nzero == 0),
            "no_molecular_features": int(audit.get("n_molecular_features_used") == 0),
            "nonzero_molecular_coefficients": nzero, "outer_warning_count": len(audit.get("warnings", [])),
            "outer_warning_fit": int(bool(audit.get("warnings"))), "outer_fallback_fit": int(bool(audit.get("fallback_reason"))),
            "fallback_reason": audit.get("fallback_reason", ""), "exact_tissue_constant_target": bool(audit.get("exact_tissue_constant_target", False)),
            "selected_genes": audit.get("n_selected_genes"), "candidate_molecular_features": audit.get("n_candidate_molecular_features"),
            "used_molecular_features": audit.get("n_molecular_features_used")})
        feature = event.get("feature_audit")
        if feature:
            genes = frozenset(feature.get("gene_ids", []))
            if genes:
                feature_sets[(context["drug_id"], context["model"], "training_MAD_candidate_genes")].append(genes)
            if context["model"] == "reo_rf":
                selected = frozenset(feature.get("random_forest", {}).get("molecular_feature_names", []))
                if selected:
                    feature_sets[(context["drug_id"], context["model"], "RF_supervised_selected_REO_pairs")].append(selected)
    task_resources = pd.DataFrame(resources)
    grouped = task_resources.groupby(["drug_id", "model"], sort=False).agg(outer_tasks=("outer_fold", "size"),
        task_wall_seconds_sum=("task_wall_seconds", "sum"), task_cpu_seconds_sum=("task_cpu_seconds", "sum"),
        worker_task_peak_rss_bytes_max=("worker_task_peak_rss_bytes", "max"),
        worker_lifetime_peak_working_set_bytes_max=("worker_lifetime_peak_working_set_bytes", "max")).reset_index()
    grouped["worker_task_peak_rss_GiB_max"] = grouped.worker_task_peak_rss_bytes_max / 1024**3
    grouped["worker_lifetime_peak_working_set_GiB_max"] = grouped.worker_lifetime_peak_working_set_bytes_max / 1024**3
    detailed = pd.DataFrame(summaries)
    zeros = detailed.groupby(["drug_id", "model"], sort=False).agg(outer_fits=("outer_fold", "size"),
        zero_molecular_coefficient_outer_fits=("zero_molecular_coefficients", "sum"),
        no_molecular_feature_outer_fits=("no_molecular_features", "sum"), outer_warning_fits=("outer_warning_fit", "sum"),
        outer_warning_messages=("outer_warning_count", "sum"), outer_fallback_fits=("outer_fallback_fit", "sum"),
        exact_tissue_constant_outer_fits=("exact_tissue_constant_target", "sum"),
        nonzero_coefficients_min=("nonzero_molecular_coefficients", "min"), nonzero_coefficients_max=("nonzero_molecular_coefficients", "max"),
        candidate_features_min=("candidate_molecular_features", "min"), candidate_features_max=("candidate_molecular_features", "max"),
        used_features_min=("used_molecular_features", "min"), used_features_max=("used_molecular_features", "max")).reset_index()
    stability = []
    expected_set_counts = detailed.groupby(["drug_id", "model"]).size().to_dict()
    for (drug, model, kind), sets in feature_sets.items():
        values = np.asarray([jaccard(a, b) for a, b in itertools.combinations(sets, 2)])
        valid = values[np.isfinite(values)]
        stability.append({"drug_id": drug, "model": model, "feature_scope": kind, "outer_feature_sets": len(sets),
            "expected_outer_fits": expected_set_counts[(drug, model)], "missing_feature_audit_sets": expected_set_counts[(drug, model)] - len(sets),
            "pairwise_comparisons": len(values), "undefined_pairwise_comparisons": len(values) - len(valid),
            "jaccard_mean": np.nan if not len(valid) else valid.mean(), "jaccard_median": np.nan if not len(valid) else np.median(valid),
            "jaccard_min": np.nan if not len(valid) else valid.min(), "jaccard_max": np.nan if not len(valid) else valid.max(),
            "interpretation": "Descriptive overlapping-training-set stability; not independent replicates, not a CI, and EN candidate genes are not nonzero biomarkers."})
    feature_columns = ["drug_id", "model", "feature_scope", "outer_feature_sets", "expected_outer_fits", "missing_feature_audit_sets", "pairwise_comparisons", "undefined_pairwise_comparisons", "jaccard_mean", "jaccard_median", "jaccard_min", "jaccard_max", "interpretation"]
    return grouped, pd.DataFrame(stability, columns=feature_columns), pd.DataFrame(parameters), zeros


def save_figure(fig, folder, stem):
    fig.savefig(folder / (stem + ".svg"), bbox_inches="tight")
    fig.savefig(folder / (stem + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def draw_summary_figures(models, effects, repeats, stability, drug_ids, folder):
    """Scientific static figures from aggregate data only, including NA rows."""
    folder.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    keys = [*map(str, drug_ids), MACRO]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
    for ax, comparison, title in zip(axes, ["primary_EN", "secondary_RF"], ["Primary: REO200 EN - Expression2000 EN", "Exploratory: REO200 RF - Expression2000 RF"]):
        part = effects[(effects.comparison == comparison) & (effects.metric == "spearman")].set_index("gdsc_drug_id").reindex(keys)
        for i, (_, row) in enumerate(part.iterrows()):
            if np.isfinite(row.point):
                if np.isfinite(row.conditional_95_lower) and np.isfinite(row.conditional_95_upper):
                    ax.plot([row.conditional_95_lower, row.conditional_95_upper], [i, i], color="#b5bdc6", linewidth=6, solid_capstyle="butt")
                if np.isfinite(row.family_adjusted_conditional_lower) and np.isfinite(row.family_adjusted_conditional_upper):
                    ax.plot([row.family_adjusted_conditional_lower, row.family_adjusted_conditional_upper], [i, i], color="#0072B2" if comparison == "primary_EN" else "#009E73", linewidth=2)
                ax.scatter(row.point, i, c="#D55E00", s=28, zorder=3)
                if row.ci_status != "reported":
                    ax.annotate("interval withheld", (row.point, i), xytext=(6, 5), textcoords="offset points", fontsize=8)
            else:
                ax.text(0.02, i, "NA (task retained)", transform=ax.get_yaxis_transform(), fontsize=8, color="#555555")
        ax.axvline(0, linestyle="--", color="#777777", linewidth=1)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Mean repeat Spearman difference")
        ax.grid(axis="x", alpha=.2)
    axes[0].set_yticks(np.arange(len(keys)), [plot_name(x) for x in keys]);axes[0].invert_yaxis()
    fig.suptitle("Conditional OOF bootstrap intervals: marginal 95% (thick grey), family-9 adjusted (thin color)", fontsize=11)
    fig.tight_layout();save_figure(fig, folder, "paired_spearman_forest")
    order = ["mean", "tissue", "expression", "expression_matched", "reo", "expression_rf", "reo_rf"]
    rho = models[models.metric.eq("spearman")].pivot(index="model", columns="gdsc_drug_id", values="point").reindex(index=order, columns=keys)
    fig, ax = plt.subplots(figsize=(14, 5))
    matrix = np.ma.masked_invalid(rho.to_numpy(dtype=float))
    cmap = plt.colormaps["RdBu_r"].copy();cmap.set_bad("#eeeeee")
    im = ax.imshow(matrix, cmap=cmap, vmin=-1, vmax=1, aspect="auto")
    ax.set_yticks(range(len(order)), [MODEL_PLOT[x] for x in order]);ax.set_xticks(range(len(keys)), [plot_name(x) for x in keys], rotation=35, ha="right")
    for y in range(len(order)):
        for x in range(len(keys)):
            value = rho.iloc[y, x]
            ax.text(x, y, "NA" if not np.isfinite(value) else f"{value:.3f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="Mean of three complete-OOF Spearman correlations")
    ax.set_title("Descriptive model correlations; macro weights all eight drugs equally")
    fig.tight_layout();save_figure(fig, folder, "model_spearman_heatmap")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
    for ax, comparison, title in zip(axes, ["primary_EN", "secondary_RF"], ["Primary EN repeat effects", "Exploratory RF repeat effects"]):
        part = repeats[repeats.comparison.eq(comparison) & repeats.metric.eq("spearman")].set_index("gdsc_drug_id").reindex(keys)
        for i, (_, row) in enumerate(part.iterrows()):
            values = np.asarray([np.nan if value is None else value for value in json.loads(row.repeat_values)], dtype=float)
            finite = np.isfinite(values)
            for repeat, value in enumerate(values):
                if finite[repeat]:ax.scatter(value, i + (repeat - 1) * .15, color=["#0072B2", "#D55E00", "#009E73"][repeat], s=24)
            if finite.all():ax.plot([values.min(), values.max()], [i, i], color="#bbbbbb", linewidth=1)
            elif not finite.any():ax.text(.02, i, "NA retained", transform=ax.get_yaxis_transform(), fontsize=8)
        ax.axvline(0, color="#777777", linestyle="--", linewidth=1);ax.set_title(title);ax.set_xlabel("Repeat-level paired Spearman difference");ax.grid(axis="x", alpha=.2)
    axes[0].set_yticks(range(len(keys)), [plot_name(x) for x in keys]);axes[0].invert_yaxis()
    fig.suptitle("Three complete repeats: colored points and descriptive range, not independent-study CI", fontsize=11)
    fig.tight_layout();save_figure(fig, folder, "paired_repeat_stability")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [2, 1]})
    genes = stability[stability.feature_scope.eq("training_MAD_candidate_genes")]
    for model in ["expression", "expression_matched", "reo", "expression_rf", "reo_rf"]:
        part = genes[genes.model.eq(model)].set_index("drug_id").reindex(drug_ids)
        axes[0].plot(range(len(drug_ids)), part.jaccard_mean, "o-", label=MODEL_PLOT[model], color=COLORS[model], linewidth=1)
    pairs = stability[stability.feature_scope.eq("RF_supervised_selected_REO_pairs")].set_index("drug_id").reindex(drug_ids)
    axes[1].plot(range(len(drug_ids)), pairs.jaccard_mean, "o-", color=COLORS["reo_rf"])
    for ax in axes:
        ax.set_xticks(range(len(drug_ids)), [plot_name(d) for d in drug_ids], rotation=45, ha="right", fontsize=8)
        ax.set_ylim(0, 1.02);ax.set_ylabel("Mean pairwise Jaccard");ax.grid(axis="y", alpha=.2)
    axes[0].set_title("Training-MAD candidate genes (not EN nonzero features)");axes[0].legend(fontsize=8, ncol=2)
    axes[1].set_title("REO RF supervised-selected pairs")
    fig.suptitle("Descriptive feature membership stability across overlapping outer training sets", fontsize=11)
    fig.tight_layout();save_figure(fig, folder, "descriptive_feature_jaccard")


def report_document(blocks):
    md, content = [], []
    for kind, value in blocks:
        if kind in ["h1", "h2"]:
            level = 1 if kind == "h1" else 2
            md.append("#" * level + " " + value);content.append(f"<{kind}>" + inline_html(value) + f"</{kind}>")
        elif kind == "p":
            md.append(value);content.append("<p>" + inline_html(value) + "</p>")
        elif kind == "table":
            md.append(table_md(value));content.append('<div class="scroll">' + value.map(fmt).to_html(index=False, border=0, escape=True) + "</div>")
        elif kind == "figure":
            caption, stem = value
            md.append(f"![{caption}](figures_formal/{stem}.png)")
            content.append(f'<figure><img src="figures_formal/{stem}.png" alt="{html.escape(caption)}"><figcaption>{inline_html(caption)}</figcaption></figure>')
        elif kind == "list":
            md.append("\n".join("- " + item for item in value));content.append("<ul>" + "".join("<li>" + inline_html(item) + "</li>" for item in value) + "</ul>")
    page = '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>固定八药重复嵌套开发评价</title><style>body{background:#f3f5f7;color:#1b263b;font-family:"Microsoft YaHei",Arial,sans-serif;line-height:1.7;margin:0}main{max-width:1200px;margin:28px auto;background:white;border-radius:12px;padding:35px}h1{font-size:28px;line-height:1.4}h2{font-size:21px;border-bottom:2px solid #dfe9f0;margin-top:34px;padding-bottom:6px}p,li{font-size:15px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px;margin:15px 0}th,td{padding:8px 10px;text-align:left;white-space:nowrap;border-bottom:1px solid #e0e7ed}th{background:#edf3f8}tr:nth-child(even){background:#f8fafc}figure{margin:22px 0}img{width:100%;height:auto}figcaption{color:#566573;font-size:13px}a{color:#165c94}@media(max-width:700px){main{margin:0;padding:18px;border-radius:0}h1{font-size:23px}}</style></head><body><main>''' + "\n".join(content) + "</main></body></html>\n"
    return "\n\n".join(md) + "\n", page


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", default="results/formal_development")
    parser.add_argument("--validation-report", default=None)
    args = parser.parse_args()
    directory = local(args.analysis_dir)
    provenance = read_json(directory / "analysis_provenance.json")
    run = local(provenance["source_run"])
    meta = read_json(run / "formal_metadata.json")
    if meta.get("completed") is not True or meta.get("completed_tasks") != meta.get("expected_tasks"):
        raise RuntimeError("Formal completion required; partial metrics are not consumed")
    validation_path = local(args.validation_report) if args.validation_report else run / "formal_validation.json"
    validation = read_json(validation_path)
    if validation.get("status") != "passed":
        raise RuntimeError("Independent full formal validation must pass before reporting")
    assert meta.get("external_response_accessed") is False and meta.get("data_guard", {}).get("blocked_attempts") == 0
    assert provenance["formal_analysis_lock_sha256"] == meta["freeze_sha256"]
    assert provenance["statistics_config_sha256"] == meta["statistics_sha256"]
    assert provenance["statistics_protocol_sha256"] == meta["protocol_sha256"]
    assert validation["formal_lock_sha256"] == meta["freeze_sha256"]
    assert validation["cohort_freeze_sha256"] == meta["cohort_sha256"]
    assert validation["execution_commit"] == meta["execution_commit"]
    assert validation["counts"]["completed_tasks_verified"] == meta["expected_tasks"]
    audits_path = run / "outer_fit_audits.jsonl"
    assert sha(audits_path) == meta["output_hashes"]["outer_fit_audits.jsonl"] == provenance["outer_fit_audits_source_sha256"]
    models = pd.read_csv(directory / "model_intervals.csv", dtype={"gdsc_drug_id": str})
    effects = pd.read_csv(directory / "paired_effect_intervals.csv", dtype={"gdsc_drug_id": str})
    repeats = pd.read_csv(directory / "repeat_stability.csv", dtype={"gdsc_drug_id": str})
    molecular = pd.read_csv(directory / "molecular_and_failure_audit_summary.csv")
    config = meta["configuration"]
    drug_ids = config["drug_ids"]
    assert len(drug_ids) == 8
    keys = [*map(str, drug_ids), MACRO]
    for comparison in ["primary_EN", "secondary_RF"]:
        expected = effects[effects.comparison.eq(comparison) & effects.metric.eq("spearman")]
        assert len(expected) == 9 and set(expected.gdsc_drug_id) == set(keys)
        assert expected.family_size.eq(9).all()
    def audit_stream():
        with audits_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.strip():yield json.loads(line)
    resources, feature_stability, choices, zeros = descriptive_audit_tables(audit_stream())
    assert int(zeros.outer_fits.sum()) == int(meta["expected_tasks"])
    recorded_molecular = molecular.groupby(["gdsc_drug_id", "model"], sort=False).outer_fits.sum().to_dict()
    assert all(recorded_molecular[(row.drug_id, row.model)] == row.outer_fits for row in zeros.itertuples())
    resources.to_csv(directory / "resource_summary.csv", index=False, encoding="utf-8-sig")
    feature_stability.to_csv(directory / "feature_stability.csv", index=False, encoding="utf-8-sig")
    choices.to_csv(directory / "chosen_parameter_summary.csv", index=False, encoding="utf-8-sig")
    zeros.to_csv(directory / "zero_and_warning_by_drug.csv", index=False, encoding="utf-8-sig")
    figure_folder = ROOT / "reports/figures_formal"
    draw_summary_figures(models, effects, repeats, feature_stability, drug_ids, figure_folder)
    def effect_table(comparison, metric="spearman"):
        data = effects[effects.comparison.eq(comparison) & effects.metric.eq(metric)].set_index("gdsc_drug_id").reindex(keys).reset_index()
        data["gdsc_drug_id"] = data.gdsc_drug_id.map(name)
        return data[["gdsc_drug_id", "point", "conditional_95_lower", "conditional_95_upper", "family_adjusted_conditional_lower", "family_adjusted_conditional_upper", "bootstrap_valid", "bootstrap_invalid", "ci_status", "undefined_reason"]].rename(columns={"gdsc_drug_id": "固定任务", "point": "平均配对差", "conditional_95_lower": "条件95%下限", "conditional_95_upper": "条件95%上限", "family_adjusted_conditional_lower": "家族9条件下限", "family_adjusted_conditional_upper": "家族9条件上限", "bootstrap_valid": "有效draws", "bootstrap_invalid": "无定义draws", "ci_status": "区间状态", "undefined_reason": "NA原因"})
    primary = effects[effects.comparison.eq("primary_EN") & effects.metric.eq("spearman") & effects.gdsc_drug_id.ne(MACRO)]
    positive, negative, equal, undefined = int(primary.point.gt(0).sum()), int(primary.point.lt(0).sum()), int(primary.point.eq(0).sum()), int(primary.point.isna().sum())
    main_zero = zeros[zeros.model.eq("reo")]
    zero_count, outer_count = int(main_zero.zero_molecular_coefficient_outer_fits.sum()), int(main_zero.outer_fits.sum())
    stability_table = repeats[repeats.comparison.isin(["primary_EN", "secondary_RF"]) & repeats.metric.eq("spearman")].copy()
    stability_table["gdsc_drug_id"] = stability_table.gdsc_drug_id.map(name)
    stability_table = stability_table[["gdsc_drug_id", "comparison", "repeat_values", "repeat_mean", "repeat_min", "repeat_max", "repeat_range", "repeat_sample_sd"]]
    zero_table = zeros[zeros.model.isin(["reo", "expression", "expression_matched", "reo_rf", "expression_rf"])].copy()
    zero_table["drug_id"] = zero_table.drug_id.map(name);zero_table["model"] = zero_table.model.map(MODEL_NAMES)
    zero_table = zero_table[["drug_id", "model", "outer_fits", "zero_molecular_coefficient_outer_fits", "outer_warning_fits", "outer_fallback_fits", "nonzero_coefficients_min", "nonzero_coefficients_max", "used_features_min", "used_features_max"]]
    sessions = meta.get("sessions", [])
    session_table = pd.DataFrame([{ "会话": i + 1, "状态": row.get("status"), "墙钟秒": row.get("wall_seconds"),
        "进程树采样峰值RSS GiB": None if row.get("sampled_process_tree_peak_rss_bytes") is None else row["sampled_process_tree_peak_rss_bytes"] / 1024**3,
        "会话完成任务": row.get("completed_tasks_this_session"), "最大inflight": row.get("max_inflight_observed")} for i, row in enumerate(sessions)])
    resource_table = resources.copy();resource_table["drug_id"] = resource_table.drug_id.map(name);resource_table["model"] = resource_table.model.map(MODEL_NAMES)
    resource_table = resource_table[["drug_id", "model", "outer_tasks", "task_wall_seconds_sum", "task_cpu_seconds_sum", "worker_task_peak_rss_GiB_max", "worker_lifetime_peak_working_set_GiB_max"]]
    aux_rows = effects[effects.comparison.isin(["REO_EN_vs_tissue", "expression_EN_vs_tissue", "matched_gene_EN_ablation"]) & effects.metric.eq("spearman")].copy()
    aux_rows["gdsc_drug_id"] = aux_rows.gdsc_drug_id.map(name)
    aux_table = aux_rows[["gdsc_drug_id", "comparison", "point", "conditional_95_lower", "conditional_95_upper", "ci_status"]]
    generated = datetime.now(timezone.utc).isoformat()
    primary_macro=effects[effects.comparison.eq('primary_EN') & effects.metric.eq('spearman') & effects.gdsc_drug_id.astype(str).eq(MACRO)].iloc[0]
    secondary_macro=effects[effects.comparison.eq('secondary_RF') & effects.metric.eq('spearman') & effects.gdsc_drug_id.astype(str).eq(MACRO)].iloc[0]
    run_wall=sum(session['wall_seconds'] for session in meta['sessions'])
    failed_sessions=sum(session['status']=='failed' for session in meta['sessions'])
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False).stdout.strip()
    blocks = [("h1", "阶段4：固定八药重复嵌套开发评价"),
        ("p", f"报告生成UTC：{generated}；正式执行提交：{meta['execution_commit']}；报告生成提交：{commit}。完整正式run和独立validator均通过后生成；此脚本仅读聚合统计、任务审计和metadata，未读取个体OOF、原响应、CTRP或表达矩阵。"),
        ("h2", "完整性与当前证据"),
        ("p", f"固定8药宏平均配对ΔSpearman：主EN {primary_macro.point:.6f}，家族9条件区间[{primary_macro.family_adjusted_conditional_lower:.6f}, {primary_macro.family_adjusted_conditional_upper:.6f}]；辅助RF {secondary_macro.point:.6f}，家族9条件区间[{secondary_macro.family_adjusted_conditional_lower:.6f}, {secondary_macro.family_adjusted_conditional_upper:.6f}]。本次固定流程的方向不支持REO优于强表达；结论限当前面板与算法，区间不是完整训练算法泛化保证。研究价值、创新与下一诊断建议详见docs/FORMAL_SCIENTIFIC_ASSESSMENT.md及唯一总计划v2.5。"),
        ("p", f"固定8药、7模型、外层{config['outer_folds']}折×{config['outer_repeats']}次、内层{config['inner_folds']}折已完成{meta['completed_tasks']}/{meta['expected_tasks']}个外层任务。主REO EN相对表达2000 EN的逐药平均ΔSpearman点方向：正{positive}项、负{negative}项、恰为0有{equal}项、无定义{undefined}项；这些方向计数本身不是确证性检验。主REO EN有{zero_count}/{outer_count}个外层重拟合为零分子系数解，须结合组织贡献解释，不能将其全局相关归为分子REO信息。"),
        ("p", f"全局独立生物组N为{provenance['n_global_biological_groups']}；每个repeat为完整的一次OOF，不把3个repeat叠成三倍N。bootstrap为{provenance['bootstrap_draws']}次同一global BIOGROUP抽样，跨药物、模型和repeat共用权重。CTRP响应访问False，访问防护尝试数0；独立药敏来源迁移和临床作用仍未验证。"),
        ("h2", "主要配对结果：REO200 EN与表达2000 EN"), ("table", effect_table("primary_EN")),
        ("p", "ΔSpearman为每药三个完整OOF相关差的算术均值；宏摘要固定等权8药。家族9包括8药与宏摘要，以Bonferroni双尾端点0.05/(2×9)提供名义家族校正条件区间，同时展示边际95%条件区间。NA任务、无定义draw及区间不足状态保留，不删除药物重算有利宏平均。"),
        ("h2", "辅助架构结果：RF"), ("table", effect_table("secondary_RF")),
        ("p", "RF结果属于独立secondary family9的探索性证据，不能在主EN弱时升为主结论；两个家族没有合并的全研究覆盖保证。REO RF包含训练内控制组织后关联筛至最多500对，表达RF没有对应pair监督筛选；差异是两条pipeline的比较，不能纯归因于排序编码。EN与RF对组织的处理也不同，跨家族比较不能证明某单一因素的作用。"),
        ("figure", ("图1：主EN与辅助RF的条件OOF重采样区间；薄彩线为家族9校正，宽灰线为边际95%，均不覆盖完整重训练算法变异。", "paired_spearman_forest")),
        ("h2", "误差指标、组织增量及模型分解"),
        ("table", effect_table("primary_EN", "rmse")),
        ("p", "RMSE配对差为参考表达减候选REO，正值表示名义误差改善；相关差则是REO减表达，正值表示名义排序改善。误差、MAE、R²及相对RMSE比例的详细区间在paired_effect_intervals.csv，属于描述性条件区间，不和主Spearman家族混为同一确证检验。回归损失用于内层选参，Spearman用于药物内排序评价，两者不保证同方向。"),
        ("table", aux_table), ("table", zero_table),
        ("p", "零分子系数EN仍可由组织均值产生全局相关及误差改善；这不是基因对标志物或组织内个体信号。没有零系数也不自动证明独立增量。上表警告与fallback仅统计外层重拟合；不将0外层warning误读为全部内层候选fit均无warning。任何有限预测的求解警告、空特征回退和真实失败都需按完整validator及运行记录继续审阅。"),
        ("figure", ("图2：各模型三个完整重复的平均OOF Spearman；训练均值NA保留，宏摘要等权固定8药，不混合药物响应计算总相关。", "model_spearman_heatmap")),
        ("h2", "重复稳定性与特征成员稳定性"), ("table", stability_table),
        ("figure", ("图3：三个重复的配对效应与描述性range；颜色表示重复，不将3次重复或15个相关折当独立研究。", "paired_repeat_stability")),
        ("figure", ("图4：训练MAD候选基因及REO RF监督筛出的pair的Jaccard；共享训练数据会影响稳定性，候选基因稳定不等于非零EN特征或生物标志物稳定。", "descriptive_feature_jaccard")),
        ("p", "feature_stability.csv报告不同外层训练集合的pairwise Jaccard及比较数。EN审计仅保存候选gene/pair名单和非零数，未保存非零系数对应的feature列表，故本报告不伪造EN active-feature Jaccard。REO RF的pair成员稳定性包含监督筛选步骤。高成员稳定性不能替代预测增量、独立验证或机制证据。"),
        ("h2", "运行资源、选参与来源绑定"), ("table", session_table), ("table", resource_table),
        ("p", f"两会话计算墙钟合计{run_wall:.1f}秒（{run_wall/3600:.2f}小时），不含中断空档、首轮输入准备及后处理。保留{failed_sessions}次控制器失败会话；模型failed为{meta['output_counts']['failed']}，全拟合solver warnings为{validation['diagnostics']['warning_fits_retained']}，记录的回退为{validation['diagnostics']['documented_fallback_fits_retained']}。673任务后Windows元数据替换WinError5已保存，非冻结监控改用追加日志、串行后处理，原数值源码未改。原673检查点逐个SHA保持的独立恢复证据见results/formal_development/recovery_and_execution_summary.json。"),
        ("p", "会话墙钟及进程树RSS采样从任务调度阶段开始，未覆盖此前输入准备，亦不包含后处理。逐模型task时间含内层搜索、外层重拟合、输入缓存及预测，sum不可当并行时的总墙钟。worker RSS含库、映射页面和当前/历史缓存；Windows lifetime peak可早于当前任务，不能把它当某模型独占或新增内存。不同进程RSS求和还可能重复计入共享页面。"),
        ("p", "chosen_parameter_summary.csv保留每个外层任务的参数、alpha/alpha_max、相对路径、理论null状态、dual gap与边界。ratio=1可以是合理null选择；最弱端点常被选则提示搜索边界。预设网格及有限迭代算法保持，不能据本结果或外验事后换网格追求REO获胜。不同候选次数及维数不能直接解释为同等调参充分度。"),
        ("p", f"正式统计协议SHA256：{provenance['statistics_protocol_sha256']}；统计配置SHA256：{provenance['statistics_config_sha256']}；formal lock SHA256：{provenance['formal_analysis_lock_sha256']}；bootstrap count matrix SHA256：{provenance['bootstrap_count_matrix_sha256']}。validator来源：{validation_path.relative_to(ROOT).as_posix()}，SHA256：{sha(validation_path)}。"),
        ("h2", "证据的实际边界与下一决定"),
        ("list", ["全部区间条件于已生成OOF记录，没有重新特征选择、调参和训练，不覆盖完整算法不确定性或保证CV泛化性能覆盖；Bonferroni调整不能修复这一基础限制。", "重复CV改善划分稳定性描述，不增加独立样本N；非显著或区间跨0不证明等效，当前没有事先建立的实际意义/等效界值。", "不存在来自bootstrap正负尾部比例的确证P值；全部8药、宏摘要、失败及NA应一并判断，不挑最佳任务。", "先审查主要效应幅度、重复方向、组织/零分子贡献和辅助架构一致性，再科学评估方法合理性、创新与发表价值。此报告不自动宣称普遍优势、临床成功或可投稿定稿。", "冻结CTRP评价、全开发最终模型锁及必要的后续研究方向判断仍是独立gate；不能因开发结果弱而恢复旧泄漏流程，或临时添加机制/深度模型包装优势。"]),
        ("p", "聚合统计与新增描述表均在results/formal_development；图形SVG/PNG在reports/figures_formal。公开输出不包含个体响应、fit成员ID或bootstrap组名单。")]
    markdown, page = report_document(blocks)
    reports = ROOT / "reports";reports.mkdir(exist_ok=True)
    (reports / "阶段4_固定八药重复嵌套评价.md").write_text(markdown, encoding="utf-8")
    (reports / "阶段4_固定八药重复嵌套评价.html").write_text(page, encoding="utf-8")
    print(json.dumps({"status": "completed", "source_run": provenance["source_run"], "report": "reports/阶段4_固定八药重复嵌套评价.html", "outer_tasks": int(zeros.outer_fits.sum())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
