"""Render aggregate-only four-case ordinal diagnostic report and scientific plots."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import html
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from analyze_formal_results import sha
from analyze_ordinal_diagnostic import MACRO

LABELS = {"tissue_null": "Tissue means", "mad_expression": "MAD / Expression",
          "mad_rank": "MAD / Rank", "mad_pair_full": "MAD / Full pair",
          "response_association_expression": "Association / Expression",
          "response_association_rank": "Association / Rank", "response_association_pair_full": "Association / Full pair"}
COLORS = {"mad_expression": "#315f9a", "mad_rank": "#258a81", "mad_pair_full": "#c16c30",
          "response_association_expression": "#315f9a", "response_association_rank": "#258a81", "response_association_pair_full": "#c16c30"}


def ci_text(row, digits=4):
    if row.ci_status != "reported":
        return "NA（" + str(row.ci_status) + "）"
    return f"{row.point:+.{digits}f} [{row.conditional_95_lower:+.{digits}f}, {row.conditional_95_upper:+.{digits}f}]"


def markdown_table(frame):
    def cell(value):
        return str(value).replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(map(cell, frame.columns)) + " |",
                      "| " + " | ".join(["---"] * len(frame.columns)) + " |"] +
                     ["| " + " | ".join(map(cell, row)) + " |" for row in frame.itertuples(index=False, name=None)])


def figures(frames, config, names, directory):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "none", "savefig.facecolor": "white"})
    directory.mkdir(parents=True, exist_ok=True)
    effects, audits, resource = frames["paired_effect_intervals"], frames["molecular_and_failure_summary"], frames["resource_summary"]
    drugs = [str(drug) for drug in config["drug_ids"]]
    models = [model for model in config["models"] if model != "tissue_null"]
    saved = []
    def save(fig, name):
        for extension in ("png", "svg"):
            fig.savefig(directory / (name + "." + extension), dpi=190, bbox_inches="tight")
        plt.close(fig)
        saved.append(name)
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.5), sharex=True)
    for ax, drug in zip(axes.flat, drugs):
        records = effects[effects.gdsc_drug_id.eq(drug) & effects.reference.eq("tissue_null") & effects.metric.eq("spearman")].set_index("candidate")
        for y, model in enumerate(models):
            row = records.loc[model]
            if np.isfinite(row.point):
                if row.ci_status == "reported":
                    ax.plot([row.conditional_95_lower, row.conditional_95_upper], [y, y], color=COLORS[model], linewidth=1.5)
                ax.scatter(row.point, y, color=COLORS[model], marker="s" if model.startswith("response") else "o", zorder=3, s=32)
        ax.axvline(0, color="#777777", linestyle="--", linewidth=.8)
        ax.set_yticks(range(len(models)), [LABELS[m] for m in models])
        ax.invert_yaxis()
        ax.set_title(names[drug])
        ax.set_xlabel("Spearman increment over same-policy tissue means")
        ax.grid(axis="x", alpha=.15)
    fig.suptitle("Observed-label OOF: paired increments and conditional 95% intervals", fontsize=13)
    fig.text(.5, .005, "Exploratory 3-fold pilot; no multiplicity adjustment or full training-refit uncertainty.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .03, 1, .97))
    save(fig, "tissue_increment_forest")

    matrix = np.full((len(drugs), len(models)), np.nan)
    for i, drug in enumerate(drugs):
        for j, model in enumerate(models):
            row = audits[audits.gdsc_drug_id.eq(drug) & audits.permutation_id.eq(0) & audits.model.eq(model)].iloc[0]
            matrix[i, j] = row.zero_molecular_fits
    fig, ax = plt.subplots(figsize=(11.5, 4.1))
    im = ax.imshow(matrix, cmap="YlOrBr", vmin=0, vmax=config["outer_folds"], aspect="auto")
    ax.set_xticks(range(len(models)), [LABELS[m].replace(" / ", "\n") for m in models])
    ax.set_yticks(range(len(drugs)), [names[d] for d in drugs])
    for i in range(len(drugs)):
        for j in range(len(models)):
            ax.text(j, i, f"{int(matrix[i, j])}/{config['outer_folds']}", ha="center", va="center",
                    color="white" if matrix[i, j] >= config["outer_folds"] * .7 else "#292929", fontsize=11)
    ax.set_title("Zero molecular coefficients in observed-label outer fits")
    fig.colorbar(im, ax=ax, label="Zero-molecular outer fits", ticks=range(config["outer_folds"] + 1))
    fig.tight_layout()
    save(fig, "zero_molecular_outer_fits")

    observed = resource[resource.permutation_id.eq(0)].groupby("model", sort=False).agg(
        seconds=("task_wall_seconds_sum", "sum"), peak=("sampled_peak_rss_bytes_max", "max")).reindex(config["models"])
    null = resource[resource.permutation_id.gt(0)].groupby("model", sort=False).task_wall_seconds_sum.sum().reindex(config["models"])
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
    x = np.arange(len(config["models"]))
    axes[0].bar(x, observed.seconds / 60, label="Observed labels", color="#315f9a")
    axes[0].bar(x, null / 60, bottom=observed.seconds / 60, label="All nine null runs", color="#a6bcd8")
    axes[0].set_ylabel("Sum of full outer-task wall time (minutes)")
    axes[0].set_title("Cost includes inner selection and outer refits")
    axes[0].legend(frameon=False)
    axes[1].bar(x, observed.peak / 2**30, color="#258a81")
    axes[1].set_ylabel("Maximum sampled worker RSS (GiB)")
    axes[1].set_title("Observed-label tasks; includes libraries and TPM pages")
    for ax in axes:
        ax.set_xticks(x, [LABELS[m].replace(" / ", "\n") for m in config["models"]], rotation=45, ha="right")
        ax.grid(axis="y", alpha=.15)
        ax.set_axisbelow(True)
    fig.tight_layout()
    save(fig, "runtime_and_memory")
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/ordinal_diagnostic_v1.json")
    parser.add_argument("--results", default="results/ordinal_diagnostic")
    args = parser.parse_args()
    config = json.loads((ROOT / args.config).read_text(encoding="utf-8"))
    directory = ROOT / args.results
    provenance = json.loads((directory / "analysis_provenance.json").read_text(encoding="utf-8"))
    if provenance["status"] != "completed" or provenance["external_response_accessed"] is not False:
        raise RuntimeError("Completed sealed statistical analysis required")
    metadata_path = ROOT / provenance["source_run"] / "ordinal_metadata.json"
    if sha(metadata_path) != provenance["metadata_sha256"]:
        raise RuntimeError("Completed run metadata changed since statistical analysis")
    run_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    frames = {}
    for filename, checksum in provenance["output_sha256"].items():
        path = directory / filename
        if sha(path) != checksum:
            raise RuntimeError("Statistical output changed: " + filename)
        frames[Path(filename).stem] = pd.read_csv(path, dtype={"gdsc_drug_id": str}, float_precision="round_trip")
    panel = pd.read_csv(ROOT / "results/cohort_v2/main_drugs.csv", dtype={"gdsc_drug_id": str})
    names = dict(zip(panel.gdsc_drug_id, panel.drug_name))
    names[MACRO] = "四诊断案例等权均值"
    report_directory = ROOT / "reports"
    stems = figures(frames, config, names, report_directory / "figures_ordinal")
    blocks, markup = [], []
    def paragraph(text):
        blocks.extend([text, ""])
        markup.append("<p>" + html.escape(text) + "</p>")
    def heading(text, level=2):
        blocks.extend(["#" * level + " " + text, ""])
        markup.append(f"<h{level}>" + html.escape(text) + f"</h{level}>")
    def table(frame):
        blocks.extend([markdown_table(frame), ""])
        markup.append('<div class="table-wrap">' + frame.to_html(index=False, escape=True, border=0) + "</div>")
    def figure(stem, caption):
        blocks.extend([f"![{caption}](figures_ordinal/{stem}.png)", "", caption, ""])
        markup.append(f'<figure><img src="figures_ordinal/{stem}.png" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>')

    heading("阶段5：候选规则、秩向量与完整基因对的诊断试跑", 1)
    paragraph("两种候选规则×三种表示，同政策组织均值参照，三种诊断药加Tamoxifen，以及九套全流程组织内置换标签。协议先于本轮拟合冻结，属于见过正式v1结果后的探索性研究。原正式v1结果保留，CTRP继续封存。")
    paragraph(f"完整执行{provenance['expected_outer_tasks']}个外层任务、{provenance['actual_fit_events']:,}次真实拟合，保存{provenance['oof_rows']:,}条OOF记录。真实独立观察为四药各620个模型，不是OOF行数。仅3折外层×1次、3折内层；没有重复CV稳定性证据。")
    heading("预设问题与解释边界")
    paragraph("在同一个训练划分内，MAD或组织调整响应关联规则各选最多200个基因；该规则的表达、样本内平均秩及完整三值pair臂共享基因。完整pair与平均秩有相同排序信息，但线性rank只对应受约束的pair权重，EN惩罚不同。因此两者的性能差不是纯维度损失；表达与rank差也不是可识别的纯幅度信息损失。")
    paragraph("同政策组织参照使用训练组织响应均值，训练未见组织用训练总体均值；这是分子系数为零时的同政策比较。Tamoxifen是历史弱信号案例，不能当作已知真阴性，出现增益本身不证明泄漏。三个主诊断案例由旧结果启发，四药均值不替代固定八药正式结论。")
    heading("真实标签下的性能与分子使用")
    point = frames["point_metrics"]
    audit = frames["molecular_and_failure_summary"]
    effects = frames["paired_effect_intervals"]
    rows = []
    for drug in map(str, config["drug_ids"]):
        for model in config["models"]:
            p = point[point.gdsc_drug_id.eq(drug) & point.permutation_id.eq(0) & point.model.eq(model)].iloc[0]
            a = audit[audit.gdsc_drug_id.eq(drug) & audit.permutation_id.eq(0) & audit.model.eq(model)].iloc[0]
            rows.append({"药物": names[drug], "模型": LABELS[model], "Spearman": f"{p.spearman:.4f}",
                "RMSE": f"{p.rmse:.5f}", "MAE": f"{p.mae:.5f}", "R²": f"{p.r2:.4f}",
                "零分子拟合": "参照" if model == "tissue_null" else f"{int(a.zero_molecular_fits)}/3",
                "预测偏离组织参照RMS": f"{p.prediction_rms_difference_from_tissue_null:.6f}"})
    table(pd.DataFrame(rows))
    paragraph("非零系数不等于有效分子信号；同时看相对同政策组织参照的外层增量、区间与预测偏离程度。零解未筛除。空活跃集合的Jaccard记为NA，不用1制造虚假的稳定性；三个重叠训练折的特征一致性仅作描述。")
    figure("zero_molecular_outer_fits", "图1：真实标签外层拟合的零分子解。")
    heading("配对增量与条件置信区间")
    paragraph("全部15种预设比较及四指标均保留。Spearman/R²正向为候选减参照，RMSE/MAE为参照减候选。用全部626个开发生物学组的同一套10,000次抽样重算配对指标；区间是未校正95%条件OOF区间，不涵盖训练、特征选择或调参重做的算法不确定性。不能据此作等效、因果比例或四药之外的普遍优势主张。")
    figure("tissue_increment_forest", "图2：六个分子臂相对同政策组织参照的配对Spearman增量。圆形MAD，方形响应关联。")
    comparison_rows = []
    for drug in [*map(str, config["drug_ids"]), MACRO]:
        selected = effects[effects.gdsc_drug_id.eq(drug)]
        for contrast in selected.contrast.drop_duplicates():
            rows_for_contrast = selected[selected.contrast.eq(contrast)].set_index("metric")
            comparison_rows.append({"药物": names[drug], "比较": contrast,
                "ΔSpearman [95%CI]": ci_text(rows_for_contrast.loc["spearman"]),
                "RMSE改善量 [95%CI]": ci_text(rows_for_contrast.loc["rmse"], 5)})
    table(pd.DataFrame(comparison_rows))
    heading("九套全流程置换控制")
    paragraph("每药在组织内置换独立生物学组的响应donor，训练与评分都使用同一套assigned标签；重新筛基因、选择参数并拟合，未沿用真实标签特征。组织关联可被保留，因此null评价的是分子臂相对同政策组织参照的增量。九次只作工程和机制诊断，不报告置换P值，也不因为某次正差或Tamoxifen改善就宣布泄漏。")
    nulls = frames["null_contrast_summary"]
    null_rows = []
    for drug in map(str, config["drug_ids"]):
        for model in config["models"]:
            if model == "tissue_null":
                continue
            selected = nulls[nulls.gdsc_drug_id.eq(drug) & nulls.candidate.eq(model) & nulls.reference.eq("tissue_null") & nulls.metric.eq("spearman")].iloc[0]
            observed = effects[effects.gdsc_drug_id.eq(drug) & effects.candidate.eq(model) & effects.reference.eq("tissue_null") & effects.metric.eq("spearman")].iloc[0]
            null_rows.append({"药物": names[drug], "模型": LABELS[model], "真实Δρ": f"{observed.point:+.4f}",
                "九次null最小": f"{selected.null_min:+.4f}", "null中位": f"{selected.null_median:+.4f}", "null最大": f"{selected.null_max:+.4f}",
                "有效次数": f"{int(selected.null_runs_valid)}/{int(selected.null_runs_requested)}"})
    table(pd.DataFrame(null_rows))
    paragraph("每次置换的全部点差见null_contrasts.csv；完整9次内有未定义指标时，摘要保持NA，不删掉该次再计算范围。singleton组织、donor未改变、响应数值未改变等实际置换覆盖如下；组织内交换性是控制成立的前提，不能假设所有组织或培养条件均可交换。")
    table(frames["permutation_coverage"])
    heading("资源、完整性和失败记录")
    resources = frames["resource_summary"]
    events = frames["all_fit_events_summary"]
    paragraph(f"所有任务累计墙钟时间{resources.task_wall_seconds_sum.sum()/3600:.2f}小时；真实标签任务{resources[resources.permutation_id.eq(0)].task_wall_seconds_sum.sum()/60:.1f}分钟；九套控制任务{resources[resources.permutation_id.gt(0)].task_wall_seconds_sum.sum()/3600:.2f}小时。每任务包含内层调参与外层重拟合；完整控制运行绝不是对已拟合模型改标签评分。")
    paragraph(f"单worker任务采样RSS最大{resources.sampled_peak_rss_bytes_max.max()/2**30:.3f}GiB，含库、表达映射页与缓存；生命周期峰值可能早于当前任务，不等于某个模型独占内存。全部拟合中的警告事件{int(events.all_fit_warning_events.sum())}、回退事件{int(events.all_fit_fallback_events.sum())}，均保留在审计文件。只读统计耗时{provenance['wall_seconds']:.1f}秒。")
    sessions = run_metadata.get("sessions", [])
    if sessions:
        session_rows = [{"运行会话": i + 1, "状态": session.get("status", "unknown"),
            "会话墙钟秒": f"{session.get('wall_seconds', float('nan')):.1f}",
            "本会话完成任务": session.get("completed_tasks_this_session", "见原始审计")}
            for i, session in enumerate(sessions)]
        table(pd.DataFrame(session_rows))
    paragraph(f"运行元数据记录的历史失败尝试归档数为{run_metadata.get('preserved_failed_attempts', '未单列，请参阅运行元数据')}；未将历史失败默认为成功或从审计中移除。CTRP响应访问为False，阻止访问尝试为{run_metadata.get('data_guard', {}).get('blocked_attempts', '未记录')}。")
    figure("runtime_and_memory", "图3：各模型完整任务的耗时与实际采样RSS。资源比较属于整条管线，不是表示本身的固有成本。")
    heading("本轮可以支持什么")
    for model in [m for m in config["models"] if m != "tissue_null"]:
        row = effects[effects.gdsc_drug_id.eq(MACRO) & effects.candidate.eq(model) & effects.reference.eq("tissue_null") & effects.metric.eq("spearman")].iloc[0]
        paragraph(f"{LABELS[model]}：四案例平均组织以上Δρ {ci_text(row)}。这只是本轮案例摘要，须与逐药效应、RMSE、零解和全流程null一起解释。")
    paragraph("既有v1差值分解见docs/OOF_DECOMPOSITION.md。新版采用不同外层划分和统一stable-ID特征列序，并且fullpair保留旧版删除的列，不能把新旧分数相减当作纯候选或编码效应。协议规定本轮不按成绩改变候选规则、基因数、调参网格或药物面板；下一步若做固定八药评价，应另行冻结完整版本，仍不得称为独立确认。CTRP继续封存。")
    heading("可复现文件")
    paragraph("唯一总计划：RESEARCH_PLAN.md；本轮科学协议：docs/ORDINAL_DIAGNOSTIC_V1_PROTOCOL.md；配置及锁：config/ordinal_diagnostic_v1.json与.lock.json。输出聚合表位于results/ordinal_diagnostic/。个体OOF、标签donor映射、完整活跃系数及bootstrap组权重留在本地runs，未混入公共报告。")
    paragraph("统计表哈希、代码哈希、共同抽样哈希、来源运行和访问审计均见analysis_provenance.json。本报告生成于UTC " + datetime.now(timezone.utc).isoformat() + "。")
    stem = "阶段5_秩与基因对诊断试跑"
    (report_directory / (stem + ".md")).write_text("\n".join(blocks), encoding="utf-8")
    style = "body{font-family:system-ui,'Microsoft YaHei',sans-serif;max-width:1480px;margin:38px auto;padding:0 26px;color:#22313f;line-height:1.7}h1{font-size:30px}h2{margin-top:38px;border-bottom:1px solid #ccd6df;padding-bottom:8px}p{max-width:1200px}table{border-collapse:collapse;min-width:750px;font-size:13px}td,th{padding:8px 11px;border-bottom:1px solid #dce4ec;text-align:left;white-space:nowrap}th{background:#eaf0f6;position:sticky;top:0}tbody tr:nth-child(even){background:#f6f8fa}.table-wrap{overflow:auto;margin:16px 0 24px;max-height:720px;border:1px solid #dce4ec}img{width:100%;max-width:1380px}figure{margin:26px 0}figcaption{color:#526171;font-size:14px}"
    page = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>' + stem + '</title><style>' + style + '</style><body>' + "\n".join(markup) + '</body></html>'
    (report_directory / (stem + ".html")).write_text(page, encoding="utf-8")
    print(json.dumps({"status": "completed", "report": "reports/" + stem + ".html", "figure_groups": stems}, ensure_ascii=False))


if __name__ == "__main__":
    main()
