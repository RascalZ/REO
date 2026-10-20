# 摘要与论文叙述：作者通读反馈讨论稿

日期：2026-10-20。状态：供讨论的写作方案，尚未替换当前作者通读Word/PDF；不是新统计方案或新科学结果。当前通读包、原R1/M08及全部锁定分析保留。不读取CTRP个体响应或预测，不训练、调参或评分。

## 1. CTRP v2名称澄清

CTRP v2是发布方使用的正式资源名：Cancer Therapeutics Response Portal, version 2。NCI CTD²将其明确标为“CTRP v2, 2015”，Broad官方门户也使用该名称，不是本项目自定义。`cohort-v2`与`rebuild_v2`才是本项目的队列/工作版本号。

本研究中CTRP v2提供外部药敏响应，表达输入沿用历史CCLE来源。因此本稿评价GDSC2到CTRP v2的跨药敏筛选迁移，不是第二个独立转录组平台验证。首次出现应展开全称，后续统一为“CTRP v2”；当前Word尚未作这一全文编辑。

官方证据：[NCI CTD²](https://ccgdata.cancer.gov/ctd2-data-project/broad-institute-screening-dependencies-cancer-cell-lines)、[Broad官方门户](https://portals.broadinstitute.org/ctrp.v2/)、[Seashore-Ludlow et al., Cancer Discovery 2015](https://doi.org/10.1158/2159-8290.CD-15-0235)。不把不同官方分析子集的规模混写为本研究源文件总规模；本稿175组/8药是自身资格筛选后的子集。

## 2. 当前稿的问题

当前Abstract约268个空白分词；问题主要在篇幅分配与论证顺序。Background首句就缩到本研究的REO增量，后续被模型臂、候选预算、分折、队列和合成控制占据。Results平列主比较、RMSE、外验绝对表现、同预算、零解、LOTO和合成七类发现，小数精度强化台账感。Conclusions主要解释“没有证明什么”，没有充分表达读者可带走的科学判断。

正文Background包含必要背景，但过早进入特征维数和信息流审计。Discussion各段的边界多数成立，不过反复收紧范围，遮住了主要解释。需要保留所有诚实限定，同时用“研究问题—关键发现—解释—使用含义”组织内容；不能以删除不利结果来改善叙述。

合适定位：指定药敏任务上的表示比较与预测来源诊断研究。工程规范支撑证据可信度，不能单独充当科学创新。主优越性假说未获支持、可直接核验的零系数状态、描述性候选政策差异和有限外验是不同证据类型，不合并成“整个项目阴性”或“所有信号消失”。

## 3. 三类写作风格

| 风格 | 叙述入口 | 最适合学习的原始研究 | 本稿采用范围 |
|---|---|---|---|
| 研究问题与评价导向 | 比较的得分是否回答了预定生物预测问题？ | DrEval；The specification game | 以组织参照、分子增量和表示差异构成主线；不宣称首次提出整个通用基准框架 |
| 模型行为与适用条件导向 | 在什么设置下模型保留分子项，在哪些设置下只输出参照预测？ | Xia跨研究比较；Herbert drug-blind研究 | 解释实际零解和跨筛选边界；不声称已识别真实药物机制或证明泛化失败的唯一原因 |
| 聚焦假说检验导向 | 匹配对照后，某种表示是否建立稳定优势？ | Ahlmann-Eltze线性扰动基线研究 | 清楚陈述本REO配置未建立平均优势；不写成等效、非劣或所有REO方法无效 |

建议以第一类组织摘要，用第二类解释结果，第三类保持结论直接。外验服务于主问题；本稿不足以写成跨平台鲁棒性研究。

三个可比较的英文入口（均为本稿新写的句子，不是论文原文引语）：

1. **问题导向：** Comparing transcriptomic representations requires separating the contribution of feature encoding from tissue structure and candidate-gene selection.
2. **行为诊断：** Molecular predictors can retain an overall prediction score when regularization selects a model with no molecular coefficients.
3. **迁移边界：** The value of a transcriptomic representation must be assessed separately from the transport of its predictions between drug-response screens.

## 4. 推荐英文摘要讨论稿

**Background:** Within-sample gene ordering is preserved under common strictly increasing transformations, but discards expression magnitude. Whether this trade-off improves drug-response prediction must be distinguished from the contributions of tissue origin and candidate-gene selection. We evaluated relative expression ordering (REO) by comparing expression, rank and gene-pair representations of the same 200 genes under two training-fold gene-selection rules, alongside tissue-only and 2,000-gene expression baselines in an elastic-net pipeline. Eight fixed compounds underwent repeated nested cross-validation in 626 cell-line groups and retrospective transfer to 175 non-overlapping groups from the Cancer Therapeutics Response Portal, version 2 (CTRP v2), with shared Cancer Cell Line Encyclopedia expression inputs.

**Results:** Gene-pair encoding did not establish an average ranking advantage over the strong-expression baseline. Mean drug-wise Spearman differences were −0.0059 in development and −0.0004 across screens, with family-adjusted conditional bootstrap intervals of −0.0187 to 0.0066 and −0.0422 to 0.0435, respectively. These intervals condition on saved predictions or fixed fitted models. Observed gains over the tissue reference in development were shared by the REO and expression models. Within matched gene sets, representation contrasts varied with the gene-selection rule without establishing an average pair advantage. In 430 of 840 molecular outer fits, the selected model had zero molecular coefficients and reproduced the tissue reference.

**Conclusions:** The results distinguish predictive signal from a representation-specific benefit. In this eight-compound setting, development increments above tissue were not unique to REO, and pair encoding did not establish an average gain over expression models. Accounting for gene-selection policy and zero molecular solutions helps interpret the value of order-based drug-response prediction.

该稿按BMC Genomics的Background/Results/Conclusions组织，低于350词，不在摘要内插文献。条件区间不是完整重训不确定性；此处没有以不显著差异证明等效。固定面板受历史开发启发、零解选择依赖调参目标、八药未新增全流程置换等限制仍必须保留在全文。文献不提供本研究数值，英文结果只来自既有M08/作者通读稿的锁定统计。

源数值到写作精度：开发主Δρ −0.005882及family-9区间[−0.018714,+0.006614]；外验主Δρ −0.00037164及family-9区间[−0.042166,+0.043505]。摘要只是展示精度变化，不修改原结果。零解430/840是拟合次数而非840个独立样本。

组织参照的开发比较来自正文既有描述性结果：主pair − tissue = +0.05654，强表达 − tissue = +0.06242，所附区间未经多重性调整；因此新摘要仅说观察到正向差值，不把这两项提升为额外确认性结论或外验已证实分子增量。数学措辞须写“共同严格递增变换”；当前正文泛写strictly monotone需在正式文字修订时改清，递减变换会翻转次序。此项不变性不等于已经证明跨平台鲁棒。

## 5. 正文重写顺序

Background四段：药敏预测中表示选择的科学意义；REO的顺序不变性动机及信息/学习权衡；候选政策、组织结构与表示难以分离的比较问题；本研究具体问题及匹配设计。将维数、折数、运行事件和审计实现放在Methods。

Results由发现领句，再给效应量及区间：先主比较，其次同候选预算的表示对照，再零分子状态，最后跨筛选的相对与绝对表现。LOTO不可估、负面结果和合成控制仍有明确位置，不因写作更简练而删除或转义。

Discussion建议依次解释：三层证据分别回答什么；候选政策如何限制表示比较的解释；零端点如何揭示预测来源（不称solver失败或真实机制证明）；跨筛选结果支持/不支持的范围；与DrEval/Specification game等已有研究的关系及REO特定贡献；集中陈述限制并给出可执行的后续检验。AIVC仅保留任务边界明确的短讨论，平衡反向文献。

可替换Discussion首段的英文示例（未写入正式稿）：

> Observed drug-response ranking and evidence for a representation-specific gain answer different questions. The evaluated models retained predictive association, yet matched comparisons did not establish an average advantage for REO over a strong expression baseline. Frequent selection of the zero molecular endpoint further showed why a nonzero prediction score cannot by itself establish a contribution from the molecular features. Together, these findings favor interpreting REO through its incremental value within a specified selection-and-learning procedure.

## 6. 原始文章与阅读范围

下列分析是结构和论证层次的借鉴，不把其他任务的实证结论当作本研究结果。

1. Bernett et al. **Critical evaluation of drug response prediction models with DrEval.** Nature Communications. 2026;17:4238. [原文/全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13168506/)，[DOI](https://doi.org/10.1038/s41467-026-72903-w)。已读摘要、背景、关键结果、Discussion及限制。可学“领域问题—评价障碍—框架—关键比较”，不照搬其通用平台贡献范围。
2. Codicè et al. **The specification game: rethinking the evaluation of drug response prediction for precision oncology.** Journal of Cheminformatics. 2025;17:33. [原文](https://link.springer.com/article/10.1186/s13321-025-00972-y)。已读摘要、主要结果与结论。可学“指标究竟测到了什么”；本稿没有做同样的全局跨药汇总/IC50浓度偏差试验。
3. Xia et al. **A cross-study analysis of drug response prediction in cancer cell lines.** Briefings in Bioinformatics. 2022;23(1):bbab356；在线2021。 [DOI](https://doi.org/10.1093/bib/bbab356)。已读可及摘要、Results和Discussion结论段，未逐项审查全部补充材料。可学将迁移差异归纳为条件性认识；本稿没有其重复测量变异和性能上界研究。
4. Herbert et al. **Monotherapy cancer drug-blind response prediction is limited to intraclass generalization.** PLOS Computational Biology. 2026;22(4):e1013232. [原始全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13128123/)。已读摘要、Introduction及Discussion，浏览关键结果。可学“现象—待检验解释—对照—适用条件”；其未见药物任务不同于本研究，机制结论不能移用。
5. Ahlmann-Eltze et al. **Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines.** Nature Methods. 2025;22:1657–1661. [原文](https://www.nature.com/articles/s41592-025-02772-6)。已核对摘要、主要结果和末尾讨论；后续直接期刊访问受限制，不声称逐项复核全部补充材料。可学主张直接且有任务边界；其遗传扰动任务不同于药敏回归。
6. 平衡来源：Miller et al. **Deep learning perturbation models can outperform baselines on calibrated metrics.** Nature Biotechnology. 2026，在线10月1日。[原文](https://www.nature.com/articles/s41587-026-03307-w)。已核对正式摘要和主要分析说明，未逐项审查补充材料。提醒基线胜负依赖任务、指标及模型；不宣称简单模型普遍胜出，零系数率本身也依赖调参目标与正则化政策。

期刊格式依据：[BMC Genomics Research Article instructions](https://link.springer.com/journal/12864/submission-guidelines/research-article)，2026-10-20核查：Abstract最多350词、三部分Background/Results/Conclusions，减少缩写且不加引用。
