# R1：既有合成DGP的可报告性核查

核查日期：2026-10-16。结论：**可以作为Supplement中的有限实现行为控制报告**，无需删除整个模块，也无需重新拟合。必须同时说明DGP的刻意构造、全部队列与缺失统计、警告以及源码和环境证据的边界。

本次只读冻结配置、源码、全部公开合成聚合表、阶段9/9b报告和该合成run的metadata。未执行生成器、拟合、原分析器、置换或重采样；未打开个体响应、个体预测、NPZ或pickle/checkpoint内容。只核对checkpoint文件是否存在及metadata/原manifest中已登记hash字典是否一致，不能把这一点写成当前重新核验了960个checkpoint payload。原分析器的完成manifest记录当时960个checkpoint hash均已验证。

## 完备性和来源链

- 冻结版本：P5e-synthetic-v1.1；锁创建于2026-10-14 08:32:46.998719 UTC，session于08:34:20.332103 UTC开始。
- 执行提交：`d151b668f916fda643e90ba1a26035cd103b9036`。
- 执行锁SHA256：`c149bb2d41d3d1a13c874e99ac73c087b25f12d4d5684c32582ff54d461043fd`。
- 14/14冻结输入hash、7/7原manifest输出hash、run metadata hash均与现存文件一致；无冻结输入被修改。
- 四场景各10队列，共40队列；每队列240个独立synthetic groups、8组织各30、220基因；3外折×1次、3内折、8臂。
- 全部960 task、26,160 fit events、76,800 OOF均在完成metadata中保留。320个队列×模型单元、960个外层task、32个场景×模型汇总均完整且无重复；所有队列行n=240、外层测试n=80。
- 对既有cohort聚合与既有summary的288个“场景×模型×指标”对应进行了标准库描述值核对，定义数/均值/中位数/最小/最大值一致；这是原汇总的数值核验，不是新增科研统计推断。
- 既有metadata记录一个已完成session，960/960任务成功；`failed_attempts`下无文件、无failure字段，fit-diagnostic fallback为0。该判断限于所保存执行记录，不声称完整机器历史从未发生其他操作。
- 原manifest和metadata均记载真实开发数据访问=false、外验响应访问=false、CTRP访问=false。本次也未访问这些材料。

逐项检查、准确值、原始文件hash和限制见 [synthetic_reportability_audit.json](synthetic_reportability_audit.json)。

## 可直接加入的完整表

| 建议编号 | 来源 | 规模 | 处理 |
|---|---|---:|---|
| S11 | `results/synthetic_mechanisms/performance_by_cohort.csv` | 320行×15列 | 完整逐字节复制，无筛选、无精度变化 |
| S12 | `results/synthetic_mechanisms/data_generation_checks.csv` | 40行×19列 | 完整逐字节复制，保留seed、全部检查与原文件hash |

复制件位于 `tables/scientific_revision_v1/`，S11的SHA256为`0c149e110e72e354705ba9513ae5ad3f52f36a54e80ff54322e64e5963f93e77`，S12为`49f6d3c609a9ae93895c12e8f8a2c0b7926d612168a18208ef0b5491f9b7bc4f`，均与原表相同。表编号为本次工作建议，最终整合时可调整显示编号而不改数据。

S11的pooled residual rho、available-tissue residual mean和strict tissue-macro residual rho各有235个空值，已全部保留；它们表示统计未定义，不是0或被删除的失败。Spearman、RMSE、MAE、R²四个point指标全部320行有定义。S12共同幅度以外场景的invariant空值表示不适用。

其余原始聚合全部保留在 `results/synthetic_mechanisms/`：外层任务960行、场景汇总32行、fit诊断32行、基因筛选32行，以及事后候选集合一致性审计240行。报告不得只展示pair表现有利的第4情形。

## 纳入稿件时必须避免的七个错误

1. **数学尺度名称。** 源变量和S12列名`log2_tpm`实际是L=log2(TPM+1)，因为输入转换是TPM=2^L−1；新Methods应按实际转换定义。旧冻结源码/表头不改。
2. **共同幅度的变换。** 方程为TPM+1乘样本因子，而非严格的TPM乘法。所有样本的gene order按构造相同，10/10队列rank/pair invariant为True。这是刻意不可见的幅度信号，不是对真实表达平台鲁棒性的检验。
3. **强表达实际基因数。** 配置预算2000，但此处只有220基因；全部强表达outer fits实际使用220。六个同预算臂为200；pair为19,900特征，表达/秩为200。
4. **零解分母。** 选定outer refit的零解为615/840，分场景为210/210、198/210、120/210、87/210。`fit_diagnostics`的13,575/26,040混合全部内层候选和外层refit，含人为必评估的alpha-ratio=1零端点，不能误称最终选择的模型零解率。
5. **警告必须保留。** 共267次fit warning，全部在共同幅度场景，三种连续表达臂各89；fallback=0。聚合文件未给警告类别，不能未经证据断言全为某一特定warning类型或把warning等同失败。
6. **第4DGP按设计有偏好。** 两个低方差driver在MAD的30个outer候选集中均未入选，在响应关联和强表达的30个中均入选。这可展示候选筛选影响信息可及性，不构成真实数据的因果机制分解或候选规则一般优劣证据。
7. **浮点与环境记录。** 纯组织情形within-tissue y range严格为0，但原OOF RMSE最大1.962615573354719e−17，Spearman不全为精确1；原值保留。合成run metadata/lock未单独保存实际package版本，不能把当前/外验环境冒充其运行实录。可准确报告源码执行版本、hash、固定requirements和实际已有资源记录；环境细节可列为复现记录局限。

## 文本整合建议

英文完整方法见 [synthetic_methods_draft.md](synthetic_methods_draft.md)，结果见 [synthetic_results_draft.md](synthetic_results_draft.md)。主文若仅保留一句，应写“在指定合成DGP中检查候选规则和表示的实现行为，全部队列结果见补充表”，而不写“验证了真实药物零解由信号稀释/有损表示造成”。

本控制没有p值、推断CI或真实标签全流程置换，不补足八药置换证据，不验证AIVC/基础模型，也不能把各效应解释为候选、幅度、维数、共线性和惩罚方式各自的因果份额。其可发表价值是可审计的有限对照及失效边界展示。

## 本次新增产物和原结果保护

仅新增本目录4份文档/审计JSON及两张精确复制表和其hash manifest；未改旧表、旧报告、冻结配置、源代码、主稿或总计划。原数据和原分析结论全部保留。
