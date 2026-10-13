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
            path = directory / (name + "." + extension)
            fig.savefig(path, dpi=190, bbox_inches="tight")
            if extension == "svg":
                path.write_text("\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8", newline="\n")
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
    axes[0].set_ylabel("Successful outer-task wall-time sum (minutes)")
    axes[0].set_title("Successful checkpoints: inner selection + outer refit")
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
    validation_path = ROOT / provenance["source_run"] / "ordinal_validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if validation.get("status") != "passed" or validation.get("diagnostic_lock_sha256") != provenance["lock_sha256"]:
        raise RuntimeError("Independent final validation of the same frozen run is required")
    overhead = validation.get("operational_overhead")
    if not isinstance(overhead, dict):
        raise RuntimeError("Operational overhead audit required before reporting total resource use")
    if overhead["successful_pipeline_fits"] != provenance["actual_fit_events"]:
        raise RuntimeError("Successful scientific fit budget differs from operational audit")
    frames = {}
    for filename, checksum in provenance["output_sha256"].items():
        path = directory / filename
        if sha(path) != checksum:
            raise RuntimeError("Statistical output changed: " + filename)
        frames[Path(filename).stem] = pd.read_csv(path, dtype={"gdsc_drug_id": str}, float_precision="round_trip")
    panel = pd.read_csv(ROOT / "results/cohort_v2/main_drugs.csv", dtype={"gdsc_drug_id": str})
    names = dict(zip(panel.gdsc_drug_id, panel.gdsc_drug_name))
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
    paragraph(f"最终成功的{provenance['expected_outer_tasks']}个外层任务包含{provenance['actual_fit_events']:,}次完整拟合，保存{provenance['oof_rows']:,}条OOF记录。该拟合数是最终纳入既定评价的预算，不含历史失败尝试中的已完成拟合及独立资源排障probe；额外开销在后文单列。真实独立观察为四药各620个模型，不是OOF行数。仅3折外层×1次、3折内层；没有重复CV稳定性证据。")
    effects = frames["paired_effect_intervals"]
    macro_pair_expression = effects[effects.gdsc_drug_id.eq(MACRO) & effects.contrast.eq("response_association:pair_full-expression") & effects.metric.eq("spearman")].iloc[0]
    paragraph("当前判读：训练内响应关联候选改善了三种表示的四案例平均预测；但该规则下pair相对表达200的宏ΔSpearman为" + ci_text(macro_pair_expression, 6) + "，仍不能证明REO可靠优势或等价。值得继续固定八药评价，研究候选规则与表示约束的适用边界；不能以本轮案例替代旧正式v1结论。")
    paragraph("新规则的活跃基因对仍不稳定，预测增量不能直接写成机制性耐药签名。完整科学解读见docs/ORDINAL_DIAGNOSTIC_ASSESSMENT.md，下一工作包和唯一维护计划见RESEARCH_PLAN.md。以下保留所有药物、表示、指标、零解与失败。")
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
    if not np.isclose(overhead["successful_task_wall_seconds_sum"], resources.task_wall_seconds_sum.sum(), atol=1e-5, rtol=1e-12):
        raise RuntimeError("Successful task-time total differs between statistical and operational audits")
    paragraph(f"最终成功checkpoint内的任务墙钟时间之和为{resources.task_wall_seconds_sum.sum()/3600:.2f}小时；其中真实标签任务{resources[resources.permutation_id.eq(0)].task_wall_seconds_sum.sum()/60:.1f}分钟，九套控制任务{resources[resources.permutation_id.gt(0)].task_wall_seconds_sum.sum()/3600:.2f}小时。每任务包含内层调参与外层重拟合；这些和数不含失败任务、排障probe、控制器准备时间或停机间隔，不能当作从首次启动到全部结束的总历时。")
    paragraph(f"成功任务内的单worker采样RSS最大{resources.sampled_peak_rss_bytes_max.max()/2**30:.3f}GiB，含库、表达映射页与缓存；生命周期峰值可能早于当前任务，不等于某个模型独占内存。最终纳入评价的拟合中，警告事件{int(events.all_fit_warning_events.sum())}、回退事件{int(events.all_fit_fallback_events.sum())}，均保留在审计文件。只读统计耗时{provenance['wall_seconds']:.1f}秒；历史失败事件另行核计。")
    sessions = run_metadata.get("sessions", [])
    if sessions:
        session_rows = [{"运行会话": i + 1, "状态": session.get("status", "unknown"),
            "会话墙钟秒": f"{session.get('wall_seconds', float('nan')):.1f}",
            "本会话完成任务": session.get("completed_tasks_this_session", "见原始审计")}
            for i, session in enumerate(sessions)]
        table(pd.DataFrame(session_rows))
    paragraph(f"运行元数据记录的历史失败尝试归档数为{run_metadata.get('preserved_failed_attempts', '未单列，请参阅运行元数据')}；未将历史失败默认为成功或从审计中移除。CTRP响应访问为False，阻止访问尝试为{run_metadata.get('data_guard', {}).get('blocked_attempts', '未记录')}。")
    table(pd.DataFrame([
        {"拟合计数类别": "最终成功评价任务中的完整fits", "数量": overhead["successful_pipeline_fits"], "统计用途": "构成冻结嵌套评价与OOF"},
        {"拟合计数类别": "失败任务中在报错前已完成的额外fits", "数量": overhead["failed_attempt_completed_fit_events"], "统计用途": "保留操作审计，不重复纳入OOF"},
        {"拟合计数类别": "独立资源排障probe的完整fits", "数量": overhead["resource_probe_fit_events"], "统计用途": "验证原失败调用的资源恢复，不是nested评价或额外外层折"},
        {"拟合计数类别": "已完成fit事件总数", "数量": overhead["completed_fit_events_total"], "统计用途": "前三项之和，包括额外运行开销"},
        {"拟合计数类别": "确认进入但未完成的失败fit调用", "数量": overhead["failed_incomplete_fit_attempts"], "统计用途": "单列，不伪装成完成事件"},
        {"拟合计数类别": "尝试fit调用数下限", "数量": overhead["attempted_fit_total_lower_bound"], "统计用途": "已完成事件加确认进入的未完成调用；是下限"},
    ]))
    def duration(value):
        return "未可靠计量" if value is None else f"{float(value):.1f}秒（{float(value)/3600:.3f}小时）"
    table(pd.DataFrame([
        {"耗时口径": "成功评价任务墙钟之和", "记录": duration(overhead["successful_task_wall_seconds_sum"]), "范围及重叠": "模型比较用，含内层调参与外层refit"},
        {"耗时口径": "失败任务墙钟之和", "记录": duration(overhead["failed_attempt_wall_seconds_sum"]), "范围及重叠": "属于训练会话内开销，未重复计入成功任务"},
        {"耗时口径": "资源排障probe", "记录": duration(overhead["resource_probe_wall_seconds"]), "范围及重叠": "独立排障开销，可能落在训练会话间隔内"},
        {"耗时口径": "训练控制器各会话墙钟之和", "记录": duration(overhead["controller_session_wall_seconds_sum"]), "范围及重叠": "包含成功/失败任务及调度、监控等；不与前两行再相加"},
        {"耗时口径": "首次会话启动至最后会话结束总历时", "记录": duration(overhead["elapsed_first_start_to_final_finish_seconds"]), "范围及重叠": "包括会话及会话之间的全部间隔，不等于持续计算时间"},
        {"耗时口径": "会话之间的间隔", "记录": duration(overhead["inter_session_gap_seconds"]), "范围及重叠": "含故障调查、probe、环境准备及其他等待，不能称为纯暂停或休眠"},
    ]))
    paragraph("额外fits与时间依据完成后的独立operational_overhead审计报告填写，不使用预测预算或进度快照替代实际计数。排障probe不产生新的评价折、不增加独立样本量，也不作为效果确认。资源恢复只应改变运行条件；模型、候选规则、参数网格和冻结源哈希是否保持不变由独立验证记录证明。")
    paragraph("耗时和内存是本机真实运行条件下的观察记录，不是隔离计算环境中的效率基准。任务顺序、表达页和库的缓存、其他应用的系统内存占用及故障后进程重启均可能影响数值；未随机化任务顺序或重复进行资源基准测量。因此这些记录用于复现预算与解释故障，不能独立证明某种表示具有普遍计算效率优势。")
    figure("runtime_and_memory", "图3：各模型最终成功任务的耗时与实际采样RSS，额外失败/probe不在图内。资源比较属于整条管线，不是表示本身的固有成本。")
    heading("本轮可以支持什么")
    for model in [m for m in config["models"] if m != "tissue_null"]:
        row = effects[effects.gdsc_drug_id.eq(MACRO) & effects.candidate.eq(model) & effects.reference.eq("tissue_null") & effects.metric.eq("spearman")].iloc[0]
        paragraph(f"{LABELS[model]}：四案例平均组织以上Δρ {ci_text(row)}。这只是本轮案例摘要，须与逐药效应、RMSE、零解和全流程null一起解释。")
    paragraph("既有v1差值分解见docs/OOF_DECOMPOSITION.md。新版采用不同外层划分和统一stable-ID特征列序，并且fullpair保留旧版删除的列，不能把新旧分数相减当作纯候选或编码效应。协议规定本轮不按成绩改变候选规则、基因数、调参网格或药物面板；下一步若做固定八药评价，应另行冻结完整版本，仍不得称为独立确认。CTRP继续封存。")
    heading("旧OOF与新诊断如何共同约束论文主张")
    paragraph("旧正式八药OOF回答已冻结v1管线的平均预测比较。路径差值分解只描述从表达2000到表达200再到REO200的两步性能变化，不能分配因果损失百分比。新四案例2×3设计回答预设表示/候选规则变化的诊断问题，不能用新三折结果覆盖旧五折×三次结论。")
    table(pd.DataFrame([
        {"可能观察到的模式": "同一候选集内rank较好、fullpair较弱", "可支持的线索": "相同序信息的模型约束与正则化方式可能影响泛化", "不能推出": "已定量证明维度或共线性造成多少损失"},
        {"可能观察到的模式": "连续表达较好，rank和fullpair均较弱", "可支持的线索": "当前序表示与EN组合可能未保留可利用的预测信息", "不能推出": "所有排序算法的上限均低于连续表达"},
        {"可能观察到的模式": "响应关联候选在多个表示中提高组织以上增量", "可支持的线索": "原MAD候选范围可能限制了当前管线", "不能推出": "训练内筛选相关性是独立机制证据"},
        {"可能观察到的模式": "非零系数增加，但外层增量不稳定或接近null范围", "可支持的线索": "更多分子项并未自动带来有用泛化", "不能推出": "仅凭非零系数增加就值得宣称预测创新"},
        {"可能观察到的模式": "某案例改善，其他案例弱或反向", "可支持的线索": "表示适用边界存在异质性，须保留全部案例", "不能推出": "删除不利药物后形成普遍优势结论"},
    ]))
    paragraph("响应关联规则使用连续log表达与组织调整药敏的关联进行选基因，本身不是对所有表示完全中立的筛选目标。共享候选集解决了三个表示输入基因不一致的问题，却不能使选择规则与表示的相互作用消失。若观察到表达臂获益更多，需要连同此设计特点解释，不能只归功于连续幅度。")
    paragraph("工程完整性和科学价值分别评价。信息流、身份、共同基因、零分子组织参照一致性和全流程null重新拟合必须通过审计；成绩是否理想不决定这些检查是否通过。真实效应高于九个null点最多只是诊断线索，不转成置换显著性。当前最有可能承载论文价值的问题是‘候选规则、表示约束和组织来源如何共同影响药敏预测’，实际创新性仍需完成固定八药受控评价及锁定外验后判断。")
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
