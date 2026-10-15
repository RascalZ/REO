# 阶段10：固定八药23组织留一组织外层敏感性分析

执行锁：`709a5499f87fdd75d4fd202041c9200a08d068938c7190e3a915dbf65edb954e`；LOTO执行commit：`fedc745940374b6870cb6c4bcb66ec8f731ded27`。
有效任务/拟合/外层预测：1,472/40,112/39,736；新模型运行耗时8.22小时。
此分析属于开发集敏感性检查，不新增确认性检验、P值或置信区间；不访问CTRP。所有预测按药物×完整留出组织评估，不能把组织fold当作来自组织总体的独立重复。

## 主要对照

主对照为每个留出组织内的 `response_association_pair_full − expression_strong` Spearman差值；主展示固定为预先锁定的10个每药均≥20独立组的组织。药物内10组织等权，八药再等权。任何主组织的相关未定义时，该药的strict均值记NA；任何一药strict值NA时八药macro记NA。可用格均值只作描述，不能替代strict值。

|药物|Δρ strict（10组织）|定义组织格|pair RMSE改善 vs 强表达|pair RMSE改善 vs 组织null|
|---|---:|---:|---:|---:|
|MIRA-1|NA|3/10|-0.0001|-0.0006|
|PLX-4720|NA|0/10|-0.0002|-0.0002|
|Bortezomib|NA|1/10|-0.0013|-0.0004|
|Telomerase Inhibitor IX|-0.0124|10/10|+0.0003|+0.0132|
|Tamoxifen|NA|0/10|+0.0000|+0.0000|
|BMS-536924|-0.0526|10/10|-0.0033|-0.0060|
|Gefitinib|NA|8/10|-0.0007|+0.0033|
|Cediranib|NA|2/10|+0.0052|+0.0021|
|八药等权宏值|NA|34/80|-0.0000|+0.0014|

组织null对所有未见组织均使用外层训练响应的全局均值；在单个留出组织内是常数，因此其Spearman按协议记NA，不能将NA填成0。RMSE改善按同一外层测试样本配对计算，正数代表pair误差更小。

## 模型零解与资源

无拟合预检记录到184/184个内层分层上下文有稀有组织警告（8×“The least populated class in y has only 2 members, which is less than n_splits=3.”; 176×“The least populated class in y has only 1 members, which is less than n_splits=3.”）。三折继续按锁定策略保留，所有外层训练样本均恰好作为一次内层验证样本；罕见组织的调参证据仍不稳定。
外层零分子系数解数量：`expression_strong` 93/184, `mad_expression` 80/184, `mad_rank` 97/184, `mad_pair_full` 120/184, `response_association_expression` 95/184, `response_association_rank` 76/184, `response_association_pair_full` 94/184.
内层/外层拟合事件中记录的警告数为0，回退数为0；零解、常数预测、NA及所有模型结果均保留。
控制器及worker进程树观测峰值RSS 1.350 GiB；worker单任务峰值RSS 1.101 GiB。系统采样下限：可用内存 0.349 GiB，commit余量 0.440 GiB。采样是操作记录，不是隔离基准或严格内存上界。

## 解释边界

LOTO评估模型对训练中未出现组织的迁移压力，不等价于独立实验室/平台验证，也不能据此证明临床可转化。留出组织的标签只用于最后计分；内层调参在其余组织内按固定的训练内分层grouped三折完成。不同组织的细胞系组成、组织标签和药敏测量仍可能构成域偏移。

全23组织的药物×组织×模型指标、主/次对照、定义格、RMSE改善和跨组织汇总见下列聚合文件。样本级ID、表达/响应、逐样本预测和模型检查点仅保留在本地忽略的`runs/`目录。

- `../results/tissue_loto/metrics_by_drug_tissue_model.csv`：23组织×八模型逐药汇总指标、覆盖和零解。
- `../results/tissue_loto/contrasts_by_drug_tissue.csv`：每药×组织主要Spearman差值与配对RMSE改善。
- `../results/tissue_loto/drug_macro.csv`：10组织严格药物均值及八药等权宏值，保留NA/定义格。
- `../results/tissue_loto/concatenated_oof_diagnostics.csv`：拼接预测相关，仅作跨组织校准敏感诊断。
- `../results/tissue_loto/loto_analysis_manifest.json`：输入/输出哈希、验证状态及访问声明。
