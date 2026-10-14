# 阶段8：固定八药OOF组织分层诊断

日期：2026-10-14。此分析复用固定八药5折×3重复OOF，不新增拟合，也不访问CTRP。它是开发集组织分层描述，不构成新的确认性检验，不增加置信区间，也不能证明泛化或无泄漏。

输入OOF SHA-256：`b2ffeeb13083340fdc2602a1d40f228615d175bf173f20a42a20c6b611e73a9e`；分析程序SHA-256：`4d641b6013657bfebdc0972dfa8915e43dc990568defd7e32a751301c5784d1b`。固定八药共119,208条预测、626个跨药物生物组、23个组织；覆盖完整。资源记录：5.07秒，采样峰值RSS 0.159 GiB。

## 指标解释

每药×组织×repeat×模型均报告绝对OOF Spearman、RMSE、`rho_residual = Spearman(y - tissue_null, prediction - tissue_null)`及RMSE相对同政策组织参照的改善量。组织null预测在外折之间可因训练组织均值不同而形成分块；其绝对Spearman只作诊断，不能解释为组织内个体排序能力。组织内残差Spearman若任一输入为常数则记NA；不赋0、不删除该组织。主展示固定为每药均≥20个独立组的10个组织，其余13个组织仍在完整明细表中。

模型对照Δρ为同一药物/组织/repeat内的`Spearman(y, candidate)-Spearman(y, expression_strong)`，RMSE改善定义为`RMSE(reference)-RMSE(candidate)`，正值有利候选。组织等权宏值同时给完整覆盖的strict值和定义格可用均值；不完整strict值为NA。repeat独立列示，不把3次重复合并为更大的样本量。

## 主展示组织上的八药汇总

表内Δρ、ΔRMSE为10个组织等权均值；Δρ strict要求10/10组织格定义。组织null模型的组织内残差Spearman按常数预测记NA。以下为描述性汇总，无P值或新增CI。

|药物|关联候选完整pair−强表达 Δρ|定义组织格/10|pair组织残差ρ（可用组织均值）|pair RMSE改善 vs null|pair RMSE改善 vs 强表达|
|---|---:|---:|---:|---:|---:|
|MIRA-1|+0.0219|10/10|+0.0404|-0.0001|+0.0002|
|PLX-4720|-0.0124|10/10|-0.1448|-0.0007|-0.0005|
|Bortezomib|-0.0536|10/10|+0.0069|-0.0011|-0.0008|
|Telomerase Inhibitor IX|+0.0204|10/10|+0.4061|+0.0110|+0.0009|
|Tamoxifen|+0.0375|10/10|+0.0232|-0.0001|-0.0001|
|BMS-536924|-0.0629|10/10|+0.1364|+0.0001|-0.0028|
|Gefitinib|-0.0023|10/10|+0.1743|-0.0000|-0.0007|
|Cediranib|+0.0454|10/10|+0.0439|-0.0011|-0.0010|

## 全量产物

- `results/tissue_sensitivity/oof_tissue_metrics.csv`：药物×组织×repeat×模型聚合指标，含全部23组织。
- `results/tissue_sensitivity/tissue_macro_by_repeat.csv`：每药×repeat×模型的10组织等权摘要，未定义格与严格NA规则均保留。
- `results/tissue_sensitivity/tissue_sensitivity_manifest.json`：输入/输出哈希、覆盖、运行时间、内存及数据访问声明。

局限：重复外折之间高度重叠；组织样本量不均衡且小组织相关极不稳定；残差相关与RMSE是不同估计目标；该分层汇总未纳入训练折/特征筛选不确定性。LOTO仍是下一步的真正整组织压力分析，不能由此处结果替代。
