# 阶段21：AIVC相邻领域讨论、证据边界与引文整合

日期：2026-10-15至16；文献检索/元数据访问截至2026-10-15。

用户批准将严格收缩后的AIVC联系加入讨论。本次仅修订两个英文工作稿的Discussion及参考文献，生成对应Word稿并更新引文台账、总计划；研究问题仍为固定八药的REO、表达量、秩及组织参照比较。未重新拟合、开展新的统计推断或读取/重跑CTRPv2响应。阶段20为26条参考文献时的历史格式审计；本记录接续维护当前稿件变化。

## Inserted Discussion paragraph

Recent single-cell perturbation benchmarks highlight the need for task-appropriate baselines and explicit evaluation targets. Simple linear baselines can match or exceed evaluated deep-learning models in some settings [27]. PerturBench examines collapsed predictions [28], while Systema addresses systematic variation in perturbation-prediction evaluation [29]. Gains over uninformative baselines have also been reported using calibrated metrics [30]. Our drug-response analysis provides a complementary, task-specific case study: descriptive representation contrasts varied with the candidate-gene policy, and approximately half of the outer-fold molecular fits selected the zero-coefficient endpoint. These observations motivate reporting covariate-reference performance, molecular increments, representation contrasts, and whether fitted predictions retain a molecular contribution. Selection of an exact zero solution in the present elastic-net pipeline depends on the tuning objective and regularization path and should not be equated with mode collapse in deep perturbation models. The synthetic controls demonstrate behavior under specified mechanisms without identifying the causes of zero solutions in real drug data. Extending these checks to single-cell foundation models requires direct, task-specific evaluation.

## 文献与论断对应

| 新编号 | 正式来源 | 本段用途 | 不支持的外推 |
|---|---|---|---|
|27|Ahlmann-Eltze et al. Nat Methods. 2025;22:1657–1661. [DOI](https://doi.org/10.1038/s41592-025-02772-6)|在被评估的扰动任务中，简单线性基线可以有竞争力|不等于全部AIVC无效；未检验药敏AUC或本研究EN零端点|
|28|Wu et al. PerturBench. NeurIPS 2025, Advances in Neural Information Processing Systems 38. [会议原文](https://proceedings.neurips.cc/paper_files/paper/2025/hash/8aee537279a66ced96319dfca3c00002-Abstract-Datasets_and_Benchmarks_Track.html)|已有预测塌缩和评价指标研究|不能声称本研究首次提供定量塌缩诊断；EN系数零与深度预测塌缩的共同机制未被证明|
|29|Viñas Torné et al. Systema. Nat Biotechnol. 2026;44:1050–1059；2025-08-25在线发表。[DOI](https://doi.org/10.1038/s41587-025-02777-8)|背景/系统变化与扰动特异评价的区分|本研究组织参照问题与单细胞共享对照伪相关不是同一数学偏倚|
|30|Miller et al. Nat Biotechnol. 2026-10-01正式在线发表。[DOI](https://doi.org/10.1038/s41587-026-03307-w)|校准指标下深度模型可超过无信息基线，呈现争论另一面|不表示所有指标、任务及未见细胞环境均已解决；尚未分配卷页时保留ahead-of-print说明|

新增元数据及完整来源存于`manuscript_v1/aivc_reference_metadata.json`与`manuscript_v1/aivc_references_additions.bib`。PerturBench按会议`inproceedings`条目处理，不冒充期刊论文；Miller的卷页以截至本日正式分配状态为准，不补造。未将scPertEval、VCBench、Nicol等2026预印本扩充为本段核心证据。

## 保留的科学边界

1. 430/840=51.19%的外层分子拟合及29/56=51.79%的最终分子臂选择合法解析零端点。两个分母是模型臂/拟合上下文，不是独立实验。将两个英文稿的“Nearly half of the outer-fold molecular fits”改为“Approximately half”；当前结果不能证明真实数据不存在分子信息。
2. association减MAD的pair政策差+0.02101不等于交互量。既有pair相对expression的差分之差+0.01539、条件区间[+0.00350,+0.02731]为未多重校正的描述性结果；本段仅写描述性表示比较随候选政策变化，不新增确证性政策或预训练语料结论。
3. 开发集pair相对组织参照宏差+0.05654、强表达+0.06242，不能转换为组织解释方差的比例。组织信息可以是当前预测任务的合法信息；是否构成混杂须结合泛化目标判断。LOTO严格宏NA表示预定统计量不能完整估计，不得当成零泛化。
4. 外验主差−0.000372是pair减强表达的相对差，不能称“全部预测信号消失”。pair减组织参照的外验描述性宏点差为+0.031466，强表达为+0.031838，未新增相应推断家族。不得把这些点差升级为已确证外部分子增量。
5. 合成DGP仅展示指定机制下固定实现的行为；尚未识别真实药物零解的原因。两道锁作为可审计执行措施保留，不主张治理方法优先权；本研究未直接评价AIVC或单细胞基础模型。标题、摘要、主终点及全部锁定结果不因本次讨论扩展而变更。

## 验收记录

- 两份正文的新增段落一致，共162个空格分隔词；独立只读审查通过，四处引用均有对应的正式原文支持。原文“Nearly half”在两稿中均改为“Approximately half”。
- 数字引用1–30连续且全部被引用；BibTeX为30个唯一条目（26篇期刊论文、1篇会议论文、3项资源），28个DOI与正文匹配，无重复。期刊Crossref原始/规范化档案均为26条，会议元数据单独归档；Systema规范化记录使用已分配的2026年6月卷期，原始提供方导出保留不改。26行PubMed/NLM缩写映射与正文、期刊DOI集合逐一匹配。
- 数字引文检查器不再硬编码23篇论文/24个DOI，分别校验期刊/会议/资源、连续编号及DOI唯一性；原26条、新30条和故障情形均经过检查。两份正式稿的17项结构检查全部通过。BibTeX元数据工具为30有效、0错误、0重复；两项建议字段警告均为Miller正式在线稿尚无volume/pages，已据Crossref/PubMed核实并明确留空，不能为消除警告而编造。
- `npm run build`通过；DOCX包含两个可编辑表格、新增段落及四个新增DOI。使用本机Microsoft Word实际渲染为28页A4；人工检查第20页讨论、第21页局限/Conclusions、第26页新参考文献、第28页末尾图注，未见裁切或遗漏。该主稿DOCX保持原有正文引用加文末图注结构，不宣称已嵌入图片。工作区另有主图嵌入审阅版及构建脚本，属于另一批修改，未纳入本次AIVC提交。
- 总计划更新为v3.33，README与证据映射同步。本次未改变模型、数据、配置、统计输出、主表、图件或外验锁；Git差异空白检查通过。机器可读验收为`manuscript_v1/aivc_integration_audit.json`。

正文仍含作者身份、声明、许可、归档DOI与AI使用说明的待确认占位；本次修改不将稿件标为投稿就绪。
