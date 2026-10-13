"""Aggregate-only report for completed fixed-eight repeated ordinal analysis."""
from __future__ import annotations
import argparse
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
from analyze_ordinal_repeated import MACRO

PRIMARY = "primary:association_pair-strong_expression"
LABELS = {"tissue_null": "Tissue means", "expression_strong": "Expression 2000 / MAD",
    "mad_expression": "Expression 200 / MAD", "mad_rank": "Rank 200 / MAD", "mad_pair_full": "Full pair / MAD",
    "response_association_expression": "Expression 200 / Association", "response_association_rank": "Rank 200 / Association",
    "response_association_pair_full": "Full pair / Association"}


def interval(row, adjusted=False):
    prefix = "family_adjusted_conditional_" if adjusted else "conditional_95_"
    lower, upper = getattr(row, prefix + "lower"), getattr(row, prefix + "upper")
    if not all(np.isfinite(value) for value in (row.point, lower, upper)):
        return "NA（" + row.ci_status + "）"
    return f"{row.point:+.5f} [{lower:+.5f}, {upper:+.5f}]"


def make_figures(frames, config, names, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "svg.fonttype": "none"})
    output.mkdir(parents=True, exist_ok=True)
    effects = frames["paired_effect_intervals"]
    scopes = [*map(str, config["drug_ids"]), MACRO]
    selected = effects[effects.contrast.eq(PRIMARY) & effects.metric.eq("spearman")].set_index("gdsc_drug_id")
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    for index, drug in enumerate(scopes):
        row = selected.loc[drug]
        if not np.isfinite(row.point):
            ax.text(0, index, "undefined / retained", va="center", color="#777777")
            continue
        for lower, upper, color, width in ((row.family_adjusted_conditional_lower, row.family_adjusted_conditional_upper, "#9fb9d7", 4.),
            (row.conditional_95_lower, row.conditional_95_upper, "#315f9a", 1.5)):
            if np.isfinite(lower) and np.isfinite(upper):
                ax.plot([lower, upper], [index, index], color=color, linewidth=width)
        ax.scatter(row.point, index, color="#315f9a", marker="D" if drug == MACRO else "o", zorder=3)
    ax.axvline(0, color="#777777", linestyle="--", linewidth=.8)
    ax.set_yticks(range(len(scopes)), [names[drug] for drug in scopes])
    ax.invert_yaxis()
    ax.set_xlabel("Mean of three complete-OOF Spearman differences")
    ax.set_title("Primary: association full pair - strong expression (MAD 2000)")
    fig.text(.5, .01, "Dark: 95% conditional interval. Pale: separate family-9 nominal interval. No retraining coverage guarantee.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, 1))
    for extension in ("png", "svg"):
        fig.savefig(output / ("primary_paired_effects." + extension), dpi=190, bbox_inches="tight")
    plt.close(fig)
    audit = frames["molecular_and_failure_summary"]
    models = [model for model in config["models"] if model != "tissue_null"]
    matrix = np.empty((len(config["drug_ids"]), len(models)))
    for i, drug in enumerate(map(str, config["drug_ids"])):
        for j, model in enumerate(models):
            matrix[i, j] = audit[audit.gdsc_drug_id.eq(drug) & audit.model.eq(model)].zero_molecular_fits.sum()
    fig, ax = plt.subplots(figsize=(12, 6.5))
    im = ax.imshow(matrix, cmap="YlOrBr", vmin=0, vmax=config["outer_folds"] * config["outer_repeats"], aspect="auto")
    ax.set_xticks(range(len(models)), [LABELS[model].replace(" / ", "\n") for model in models], rotation=30, ha="right")
    ax.set_yticks(range(len(config["drug_ids"])), [names[str(drug)] for drug in config["drug_ids"]])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{int(matrix[i,j])}/15", ha="center", va="center", color="white" if matrix[i,j] >= 11 else "#292929")
    ax.set_title("Zero molecular coefficients across five outer folds x three repeats")
    fig.colorbar(im, ax=ax, label="Zero-molecular outer fits")
    fig.tight_layout()
    for extension in ("png", "svg"):
        fig.savefig(output / ("zero_molecular_fits." + extension), dpi=190, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/ordinal_repeated_v1.json")
    parser.add_argument("--run-dir")
    parser.add_argument("--results", default="results/ordinal_repeated")
    args = parser.parse_args()
    config_path = ROOT / args.config
    config = json.loads(config_path.read_text(encoding="utf-8"))
    directory = ROOT / args.results
    provenance = json.loads((directory / "analysis_provenance.json").read_text(encoding="utf-8"))
    if provenance.get("status") != "completed" or provenance.get("config_sha256") != sha(config_path):
        raise RuntimeError("Completed statistics from the current frozen config required")
    run = ROOT / (args.run_dir or provenance["source_run"])
    metadata_path, validation_path = run / "ordinal_repeated_metadata.json", run / "ordinal_repeated_validation.json"
    if sha(metadata_path) != provenance["metadata_sha256"] or sha(validation_path) != provenance["validation_sha256"]:
        raise RuntimeError("Completed metadata or independent validation changed")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if metadata.get("completed") is not True or validation.get("status") != "passed":
        raise RuntimeError("Paused/partial/failed runs cannot produce a completed research report")
    frames, overhead = {}, None
    for filename, checksum in provenance["output_sha256"].items():
        path = directory / filename
        if sha(path) != checksum:
            raise RuntimeError("Changed statistical aggregate: " + filename)
        if filename.endswith(".csv"):
            frames[Path(filename).stem] = pd.read_csv(path, dtype={"gdsc_drug_id": str}, float_precision="round_trip")
        elif filename == "operational_overhead.json":
            overhead = json.loads(path.read_text(encoding="utf-8"))
    if overhead is None or overhead["successful_pipeline_fits"] != provenance["actual_fit_events"]:
        raise RuntimeError("Operational fit budget missing or inconsistent")
    panel = pd.read_csv(ROOT / "results/cohort_v2/main_drugs.csv", dtype={"gdsc_drug_id": str})
    names = dict(zip(panel.gdsc_drug_id, panel.gdsc_drug_name))
    names[MACRO] = "Equal eight-drug macro"
    output = ROOT / "reports"
    make_figures(frames, config, names, output / "figures_ordinal_repeated")
    md, markup = [], []
    def heading(text, level=2):
        md.extend(["#" * level + " " + text, ""])
        markup.append(f"<h{level}>" + html.escape(text) + f"</h{level}>")
    def paragraph(text):
        md.extend([text, ""])
        markup.append("<p>" + html.escape(text) + "</p>")
    def table(frame):
        def cell(value):
            return str(value).replace("|", "\\|").replace("\n", " ")
        md.extend(["| " + " | ".join(map(cell, frame.columns)) + " |", "| " + " | ".join(["---"] * len(frame.columns)) + " |"])
        md.extend("| " + " | ".join(map(cell, row)) + " |" for row in frame.itertuples(index=False, name=None))
        md.append("")
        markup.append('<div class="table-wrap">' + frame.to_html(index=False, escape=True, border=0) + "</div>")
    def figure(stem, caption):
        path = "figures_ordinal_repeated/" + stem + ".png"
        md.extend(["![" + caption + "](" + path + ")", "", caption, ""])
        markup.append('<figure><img src="' + path + '" alt="' + html.escape(caption) + '"><figcaption>' + html.escape(caption) + '</figcaption></figure>')
    heading("阶段6：固定八药候选规则与表示的重复嵌套评价", 1)
    effects = frames["paired_effect_intervals"]
    primary = effects[effects.contrast.eq(PRIMARY) & effects.metric.eq("spearman")]
    macro = primary[primary.gdsc_drug_id.eq(MACRO)].iloc[0]
    lower, upper = macro.family_adjusted_conditional_lower, macro.family_adjusted_conditional_upper
    if not np.isfinite(lower) or not np.isfinite(upper):
        conclusion = "主比较存在无定义组成或有效抽样不足，区间扣留，不能从成功子集宣称优势。"
    elif lower > 0:
        conclusion = "主宏比较的名义条件区间位于正向，支持当前固定流程的开发集平均排序增量；效果是否足够有用及能否迁移仍待独立判断。"
    elif upper < 0:
        conclusion = "主宏比较的名义条件区间位于负向，当前REO流程平均排序弱于强表达基线。"
    else:
        conclusion = "主宏比较的名义条件区间跨零，现有证据不足以认定REO优于强表达基线；也不能认定非劣或等价。"
    paragraph(conclusion + " 主宏ΔSpearman及family9区间：" + interval(macro, adjusted=True))
    paragraph(f"固定八药×8臂×5外层×3重复、3内层已完成：{provenance['expected_outer_tasks']}个成功任务、{provenance['actual_fit_events']:,}次成功评价拟合、{provenance['oof_rows']:,}条OOF。独立生物组{provenance['n_observed_biological_groups']}个；重复与模型记录不能扩大独立N。既有v1和四案例结果保留，CTRP继续封存，本轮没有新增置换。")
    heading("指标、配对与推断范围")
    paragraph("每种药物、每个模型先分别计算三份完整OOF的指标，再平均三个指标；不先平均预测再求相关，不将15个相关外折当独立试验。所有药物、模型和重复共享626组的10,000次bootstrap，Spearman按每次抽样的重复数重新计算精确并列秩。任何组成无定义则保留NA，不能删除药物或重复后重新定义宏平均。")
    paragraph("主比较是响应关联候选pair相对MAD2000强表达的Spearman，8药+宏平均构成family9。两种规则各自pair−表达200的Spearman构成单独family18。其它比较及RMSE/MAE/R²为未校正描述；两个家族不提供全研究共同覆盖保证。区间均条件于现有OOF，不包含重新选特征、调参及训练的不确定性；family18双尾每端约14个bootstrap次序量，名义尾部区间不能解读为精确算法覆盖或P值。")
    paragraph("新版受历史开发结果启发，不是独立确证。未论证非劣界值、未执行非劣或等效检验；区间跨零不证明信息完全保留。三个候选×表示差分之差在同次抽样和重复内计算、存在代数关联，不是幅度、维数、共线性或惩罚的因果分摊。")
    figure("primary_paired_effects", "图1：主比较的三次完整OOF指标平均差；两种区间均为条件区间。")
    rows = []
    stability = frames["repeat_stability"]
    for drug in [*map(str, config["drug_ids"]), MACRO]:
        row = primary[primary.gdsc_drug_id.eq(drug)].iloc[0]
        s = stability[stability.gdsc_drug_id.eq(drug) & stability.comparison.eq(PRIMARY) & stability.metric.eq("spearman")].iloc[0]
        error = effects[effects.gdsc_drug_id.eq(drug) & effects.contrast.eq(PRIMARY) & effects.metric.eq("rmse")].iloc[0]
        rows.append({"药物": names[drug], "Δρ及family9区间": interval(row, True), "边际95%区间": interval(row),
            "三个重复Δρ": s.repeat_values, "RMSE改善量及95%CI": interval(error)})
    table(pd.DataFrame(rows))
    heading("匹配200基因的表示比较与候选规则作用")
    rows = []
    for drug in [*map(str, config["drug_ids"]), MACRO]:
        for rule in ("mad", "response_association"):
            row = effects[effects.gdsc_drug_id.eq(drug) & effects.contrast.eq(rule + ":pair_full-expression") & effects.metric.eq("spearman")].iloc[0]
            rows.append({"药物": names[drug], "候选规则": rule, "pair−表达200及family18区间": interval(row, True), "边际95%区间": interval(row)})
    table(pd.DataFrame(rows))
    macro_rows = effects[effects.gdsc_drug_id.eq(MACRO) & effects.metric.eq("spearman")]
    table(pd.DataFrame([{"宏比较": row.contrast, "Δρ及边际95%CI": interval(row), "家族": row.family_name} for row in macro_rows.itertuples()]))
    paragraph("完整25种比较×9个单位×4指标均保存于paired_effect_intervals.csv。共享候选基因使三种表示输入范围相同；响应关联规则以连续log表达为筛选基础，规则与表示仍可能交互。完整pair与平均秩承载相同序信息，但线性函数范围和EN惩罚不等价。")
    heading("零分子解、重复波动与特征稳定性")
    figure("zero_molecular_fits", "图2：每药每臂15次外层拟合中的零分子解，包括强表达基线。零解未剔除。")
    audit = frames["molecular_and_failure_summary"]
    molecular_audit = audit[audit.model.ne("tissue_null")]
    paragraph(f"840个分子外层模型中零分子解{int(molecular_audit.zero_molecular_fits.sum())}个。各重复预测相对同政策组织均值的RMS偏离单列；非零系数和预测变化不自动等于有效分子增量。恰好相同预测造成的[0,0]区间只是本批记录相同，不是算法等价证据。")
    paragraph("repeat_stability.csv给出每个模型/比较的三个完整重复指标、范围及SD，不把三个重复当三项独立研究。feature_stability.csv比较15个外层模型的105组组合，区分同重复不同折与不同重复，记录训练样本Jaccard以显示重叠。双空活跃集合记NA；不同维数空间的活跃Jaccard不能直接当算法优劣排名，预测增量也不能代替稳定机制签名。")
    heading("计算、暂停和失败的实际记录")
    resources = frames["resource_summary"]
    events = frames["all_fit_events_summary"]
    paragraph(f"最终成功任务墙钟和{resources.task_wall_seconds_sum.sum()/3600:.2f}小时，最大任务采样RSS {resources.sampled_peak_rss_bytes_max.max()/2**30:.3f}GiB；外层任务包括初始化/内层调参/重拟合。全成功预算内求解警告{int(events.all_fit_warning_events.sum())}、回退{int(events.all_fit_fallback_events.sum())}。资源是本机实际观察，受缓存、任务顺序、系统负载和重启影响，不能声称隔离计算效率优势。")
    keys = ("successful_pipeline_fits", "failed_attempt_completed_fit_events", "resource_probe_fit_events", "completed_fit_events_total",
        "failed_incomplete_fit_attempts", "attempted_fit_total_lower_bound", "successful_task_wall_seconds_sum", "failed_attempt_wall_seconds_sum",
        "resource_probe_wall_seconds", "controller_session_wall_seconds_sum", "elapsed_first_start_to_final_finish_seconds", "inter_session_gap_seconds",
        "resource_paused_session_count", "failed_session_count")
    table(pd.DataFrame([{"操作审计项目": key, "实际值": overhead.get(key, "未单列")} for key in keys]))
    paragraph("成功任务、失败尝试、资源probe和未完成调用分别计数。任务wall从execute_task入口计时，包括worker初始化、数据准备及嵌套评价，但不包括此前子进程启动和数值库导入；不能将它称为fresh-child完整端到端耗时。若会话审计另列child_wall_seconds_sum，则包含子进程启动与等待。控制器会话时间包含任务时间，不能重复累加；首末历时还包含会话间隔，间隔可能含资源等待、调查和人工暂停，不等于纯休眠。资源不足的安全暂停不算模型失败，未完成运行不得生成最终成绩报告。")
    heading("下一步边界与可复现材料")
    paragraph("是否继续外验应同时考虑逐药与宏效应、误差指标、组织以上增量、重复波动和特征稳定性，不能只挑排序点估计或表现好的药物。当前评价不证明临床转化，也未检验RNA-seq至芯片迁移。最终全开发模型与外验统计方案另行冻结后，才能执行既定CTRP跨药敏筛选评价。")
    paragraph("统计来源与代码哈希见results/ordinal_repeated/analysis_provenance.json；候选/活跃特征、参数边界、零解、资源及覆盖有独立聚合表。个体OOF、训练成员和bootstrap组权重仍留在本地runs；本报告只读取完成后的元数据、独立审计和聚合结果。唯一总计划为RESEARCH_PLAN.md。")
    stem = "阶段6_固定八药表示与候选规则重复评价"
    (output / (stem + ".md")).write_bytes(("\n".join(md) + "\n").encode("utf-8"))
    style = "body{font-family:system-ui,'Microsoft YaHei',sans-serif;max-width:1450px;margin:36px auto;padding:0 24px;color:#22313f;line-height:1.7}h2{margin-top:34px;border-bottom:1px solid #ccd6df}table{border-collapse:collapse;font-size:13px;min-width:750px}th,td{padding:8px 10px;border-bottom:1px solid #dce4ec;text-align:left;white-space:nowrap}th{background:#eaf0f6}tbody tr:nth-child(even){background:#f6f8fa}.table-wrap{overflow:auto;max-height:740px;margin:20px 0}img{width:100%}figure{margin:24px 0}figcaption{font-size:14px;color:#526171}"
    page = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>' + stem + '</title><style>' + style + '</style><body>' + "\n".join(markup) + '</body></html>'
    (output / (stem + ".html")).write_bytes(page.encode("utf-8"))
    print(json.dumps({"status": "completed", "report": "reports/" + stem + ".html"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
