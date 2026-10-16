# Comments to authors

Author-requested AI-assisted working audit dated 2026-10-16; not an editorial decision. Eight open revision groups follow. No new model fitting or external scoring was performed.

### Major comment M01
- Location: Manuscript_BMC_Genomics_draft.md:29,35,179,181
- Observation: 创新性与最接近的DrEval、Specification Game及稀疏kTSP研究区别尚不具体。
- Evidence or criterion: REF05;REF06;REF14;METHODS-REPRESENTATIONS
- Why it matters: 通用基线、嵌套CV、组织参照及跨筛选已经有先例，若作为新颖性核心容易被认为重复；这不等于当前REO特定问题没有价值。
- Requested action: Background/Discussion明确继承既有benchmark原则，将贡献限定为同候选集合、两政策×三表示、同一组织调整EN与解析零端点的实证比较；注明未评价稀疏kTSP或全部REO方法。

### Major comment M02
- Location: Manuscript_BMC_Genomics_draft.md:61,83,89,95,99-107
- Observation: 癌资格/身份QC、完整分折种子与稀有组织处理、正式运行版本、相关/有效bootstrap判据仍未完整报告。
- Evidence or criterion: COHORT-RULES;DEVELOPMENT-CONFIG;EXTERNAL-LOCK;LOTO-CONFIG;NESTED-CODE;STATISTICS-CODE
- Why it matters: 当前数值能追溯，但读者不能仅由投稿包重建主要准入与估计步骤；无solver警告不等于无splitter警告。
- Requested action: 补充方法/运行版本表：Non-Cancerous/unknown/确证例外、全GDSC亲缘排除及主外验SNP-QC；StratifiedGroupKFold、内折种子repeat_seed+100+outer_fold、最终2026101501、LOTO2026101401+组织序号；稀有组织警告及未知组织政策；实际Python/包版本与执行commit；相关最少3组、恒定NA、至少99%有效draw、linear分位插值。

### Major comment M03
- Location: Manuscript_BMC_Genomics_draft.md:187,189
- Observation: Discussion引用了合成控制的结果性解释，但本文Methods、Results及8个补充文件没有提供其设计/结果出处。
- Evidence or criterion: MANUSCRIPT;SUPPLEMENT-LIST;STAGE09
- Why it matters: 审稿人无法在本投稿包核查其支持作用；写得谨慎不能替代证据可定位。
- Requested action: 优先删除依赖未报告合成实验的结果性句子，保留当前八药未做全流程置换的局限；若确需保留，则补简短Methods与已存在聚合结果/固定版本链接。

### Major comment M04
- Location: docs/DATA_AVAILABILITY.md:Original project compatibility inputs; manuscript Methods:107; README public-artifact scope
- Observation: 公众路径仍依赖legacy桥表/四列表；实际候选/active特征、模型预处理与折定义缺少公开或申请获取说明。
- Evidence or criterion: DATA-AVAILABILITY;INPUT-PATHS;PUBLIC-REPO;LICENSE
- Why it matters: 源码和聚合hash可审计不等于可从公开输入复现模型；编辑/审稿人需要实际测试路径。
- Requested action: 区分历史兼容审计与当前论文所需输入；列准确release/成员文件/URL/hash；整理允许分享的最小映射、fold定义、特征/系数/零解与模型schema，给公开或受限获取路径及最小推理示例。不得冒称64训练模型已公开。

### Major comment M05
- Location: Manuscript_BMC_Genomics_draft.md:9-13,109,199-223
- Observation: 作者/机构/通讯、伦理/同意、基金/利益/贡献/致谢、数据代码声明和AI说明仍是占位。
- Evidence or criterion: MANUSCRIPT;POLICY-RESEARCH;POLICY-AI
- Why it matters: 稿件明确不是投稿就绪版本，不能从公开数据或空白自动推定声明内容。
- Requested action: 按author_inputs.md一次收集作者事实；AI声明覆盖实际设计、代码、分析核查与写作辅助，最终由作者审核承担责任；另列Code Availability并核对数据版本及引用。

### Major comment M06
- Location: CITATION.cff:5-10; manuscript Methods:107 and Availability:207
- Observation: 代码为MIT，但CITATION仍v0.1.0/2026-10-10、通用贡献者，无当前投稿档案版本/持久标识；上游数据许可不能由MIT替代。
- Evidence or criterion: LICENSE;CITATION;POLICY-CODE;DATA-AVAILABILITY
- Why it matters: 浮动GitHub主分支不能替代可识别的归档版本；必须准确陈述发布范围及访问条件。
- Requested action: 科学整改完成后选定投稿快照tag/commit、作者与文件清单，并取得真实归档标识。DOI是推荐实现而非唯一允许选项；核对2018 CCLE与GDSC具体原文件再分发范围、CMP条款和CTRP/GENCODE归属。文章CC BY或CC BY-NC-ND由作者依资助要求选择。

### Major comment M07
- Location: Manuscript_BMC_Genomics_with_figures.md; main_figure_embedding_audit.json
- Observation: 本地含图衍生稿仅26条References却包含27–30引用；其绑定的主稿MD/DOCX SHA均已过期。
- Evidence or criterion: SNAPSHOT-AUDIT;AIVC-AUDIT
- Why it matters: 会把旧参考文献及当前讨论混装，不能作为新版投稿或发给合作者的最终稿。
- Requested action: 以30条参考文献的主BMC稿为权威；所有内容修订完后重建主DOCX，再从该版本生成含图审阅稿并重验图数、表格、正文、引文及source SHA。原旧副本保留并标明历史。

### Minor comment M08
- Location: BMC manuscript Background:29; Figure1-3 graphic titles; final References
- Observation: 数字引文覆盖通过，但首次出现顺序为1,2,3,22,23,4…；部分图件内仍有重复总图题。
- Evidence or criterion: SNAPSHOT-AUDIT;POLICY-FORMAT;NLM-CITATION-SEQUENCE
- Why it matters: 覆盖/元数据正确不等于最终BMC/Vancouver样式合规；图题应在正文图注，分面标题/键可留图中。
- Requested action: 用最终BMC/Vancouver输出统一作者缩写、et al规则及按首次引用的编号，同时更新NLM映射/引文检查；移除图内重复总图题，保留必要面板标识和图例。

# Confidential comments to editor

Not applicable: this is author-requested manuscript preparation. There is no journal assignment, no editor-only submission, and no independent editorial recommendation. The human authors must verify the audit and final manuscript.
