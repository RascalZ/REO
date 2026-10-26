contract_role: eic
## Dimension Scores

### D1: methodology_rigor
score: not_assessed

### D2: domain_accuracy
score: not_assessed

### D3: argumentative_coherence
score: not_assessed

### D4: cross_disciplinary_relevance
score: not_assessed

### D5: writing_and_structure
score: warn
trigger: "局部定义、段落顺序、图表说明、交叉引用或冗余表述增加阅读负担"

### D6: venue_fit_and_contribution
score: pass

## Review Body

本席采用计算基因组学与转录组方法研究的编辑视角，身份并非真实期刊编辑。已审阅冻结英文主文、Additional file 14 补充文本、附表说明和四幅原始 PNG；未访问其他席位报告、项目历史、个体 CTRP 数据或模型输出，也未重算结果。D1–D4 留给相应席位，本报告的分数仅覆盖 D5、D6。没有修改 Phase 1 判据。

criteria_binding_unavailable

所给 BMC 官方政策背景仅供理解作者目标，没有形成正式目标标准绑定。本席不作期刊合规认证、文献首创认证或基于预计影响力的门槛判断。D6 的通过依据是稿件存在可辨认的有限经验贡献：在明确的八药、选基因和学习流程内，区分总体预测、相对组织基线的分子增量与编码方式的额外收益；阴性比较结果本身不构成贡献不足。D5 的警示针对局部阅读负担和图注说明，不表示主要比较或核心结论无法辨认。

需要明确区分证据状态。本席未发现可由所给材料直接推翻的核心贡献主张。条件区间不覆盖重新选特征、调参和拟合的不确定性，外验共用 CCLE 表达源且属于回顾性跨屏评估，LOTO 严格宏值未定义，以及合成控制无法识别真实数据机制，均已在稿件中披露；这些限制需要保留在结论范围内，不能改写成稿件声称了等效、跨表达平台鲁棒性或临床有效性。下列 W2 是需要核对既有状态记录的图示查询，不是已证实的统计错误；W3、W4 是作者阅读稿明确延后的投稿准备事项，不是已发生的数据造假、来源虚报或科学结论被推翻。

### S1: 主问题、实际流程比较与匹配候选集比较可以清楚区分

**Evidence Anchor**: text: 主文 Background 最后一段（以 Here, we asked 开头），“The primary comparison evaluated a response-associated gene-pair model against a strong 2,000-gene expression baseline.”

引言明确提出实际 REO 流程与强表达基线的主比较，随后说明同一 200 基因集合内的三种表示比较。Methods 的 estimands 和 Results 的章节顺序延续了这个结构，避免把候选集变化与表示变化当作同一个问题。这个设计说明足以让读者识别稿件评估的对象，无需把文章包装成新算法论文。

### S2: 近零主比较被表述为有边界的经验结果

**Evidence Anchor**: table: 主文 Table 2，Equal-drug macro 行及表题对 development/external conditional intervals 的限定

主文保留开发与跨屏的全部八药结果、正负点估计、区间和三个精确零差值，宏差值分别约为 −0.0059、−0.0004。正文没有把区间包含零解释成等效或非劣效；同时报告绝对外验排序关联，避免将表示差值接近零误写成所有可迁移预测关联消失。对相关研究读者而言，新增价值在于这种比较解释与诊断，而非必须出现正向性能提升。这个贡献范围能从当前稿件直接辨认；其优先权与领域事实仍需领域席位作独立核验。

### S3: 零分子解与保留组织预测的区别呈现清楚

**Evidence Anchor**: figure: Figure 3B，七个分子模型臂的 full-development molecular coefficients 状态矩阵及灰色/青色图例

图 3B 明确展示的是最终拟合模型的分子系数状态，图注也说明它不是外层折零解频率。正文将 430/840 外层分子拟合与 29/56 最终分子模型分别给出，并说明零分子解仍含组织分量、非零系数本身不证明预测用途或机制。补充 SM4 给出预测恒等式，使这个诊断不仅是笼统的模型失败描述，也没有把它夸大成真实数据中不存在分子信号。

### S4: 主要图示与外验边界保持一致

**Evidence Anchor**: figure: Figure 1，底部 Locked CTRPv2 cross-screen evaluation 框及右上 Expression input 框

图 1 将表达源共用、175 个不重叠外部组和一次性回顾排序写入流程图；图 2 与图 4B 的横轴给出差值方向，并标注 family-adjusted conditional intervals。图 4A 明示 34/80 与宏值 NA。这些可直接看到的限定有助于使图示与正文的有限结论一致，而不是只在讨论末尾加入免责说明。

### W1: 主文执行日志细节过密，削弱方法层级

**Severity**: Minor
**Evidence Anchor**: text: 主文 Methods—Leave-one-tissue-out sensitivity analysis 第 2 段，“Before fitting, execution amendment v1.1 resolved a prose/config seed conflict in favor of the frozen machine value 2026101401.”
**Confidence**: 5 — 主文与补充 SM3 对读的编辑结构判断。

这是一项表达结构意见。主文在介绍 LOTO 训练和评估边界后，立即进入版本号、种子冲突、日志和行标识修复；Reproducibility and data governance 又连续列出三个科学源提交哈希。相同的可追溯细节已在 SM3、SM6 与 S13 说明。它们有记录价值，但在方法主文中占据的层级与核心设计近似，使读者较难快速识别哪些信息改变了科学流程、哪些仅描述执行修复。现有稿件仍能追踪问题和结论，因此不达到 D5 阻断程度。

最低修改是保留主文中有关训练边界、预先冻结、实际科学版本及修复未改变分析选择的一两句概述，把具体修复序列、逐项种子和长哈希集中在现有补充材料，并明确指向对应小节。不要删去影响科学选择的修订披露，也不需要为此重跑分析。

### W2: LOTO 灰色单元尚不能让读者区分未定义的统计原因

**Severity**: Minor
**Evidence Anchor**: figure: Figure 4A，标为 Undefined Δρ 的灰色单元、各药的 defined-count 行标与 34/80 注记
**Confidence**: 4 — 对现有图注可读性的直接判断，实际原因分布尚需作者核对。

这是一项需验证的编辑查询。图 4A 正确地报告 definedness，并未把 NA 填成零；但灰色将所有未定义情形合并。补充 SM5 列出样本不足、常数响应、常数预测或非有限相关等可产生 NA 的规则，主文没有说明在该十组织网格中实际出现哪些原因。由于读者已被告知这十个组织满足样本覆盖门槛，若只看到灰色，容易把未定义误认成数据缺失或笼统的模型失败。这里不应未经记录核对就声称所有灰格都由零分子解导致，更不能据此推断真实数据没有信号。

最低修改是核对已存在的 S1/S7 状态字段，在图注或结果段补一句实际未定义原因及其范围；若有多个原因，可用小表或图注简要分列，并继续保留严格宏值 NA。补充说明应使用已有记录；若原因尚未整理，应坦陈未分类，不能补造原因或以可用单元均值替代预定宏值。这个请求改善解释，不要求改变敏感性分析或补做完整跨组织验证。

### W3: 代码与完整派生模型资产的访问安排仍待落实

**Severity**: Minor
**Evidence Anchor**: text: 主文 Code availability，“Public or reviewer access and an archival release remain to be arranged before submission.”
**Confidence**: 5 — 稿件明确披露的当前访问状态。

这是延后的投稿准备工作。稿件诚实说明仓库仍为私有，身份映射、原始折分、64 个数值模型、系数和预处理状态虽已在本地准备，却未包含于版本控制仓库；因此不能把当前材料视为已公开、可直接复用的完整基准资产。这个访问限制并不反证主文的近零比较结论，本席也没有把本地准备写成经过独立复现。

最低完成条件是为代码和支撑结果的最低派生材料逐项说明可执行的获取方式、访问限制、责任联系人与版本，并在提交时落实审稿可访问安排。可按上游许可约束提供受控访问，不要求重新分发全部原始表达和药敏矩阵；持久归档是可选的加强方案，不能把某一特定平台或 DOI 当作当前冻结判据。访问安排落实后再更新声明，现阶段不应宣称已经公开完成。

### W4: 作者贡献与致谢仍为待填声明

**Severity**: Minor
**Evidence Anchor**: text: 主文 Authors' contributions，“Author-specific contribution roles will be completed by the listed authors before submission.”；Acknowledgements，“The acknowledgement statement will be finalized by the authors before submission.”
**Confidence**: 5 — 声明中的明确占位表述。

这是另一项投稿准备工作，与结果真伪或算法贡献强弱分开处理。稿件已经披露 AI 协助范围，并将判断、核验与批准责任归于人类作者，但作者贡献和致谢尚未完成。最低修改是由所列作者确认具体贡献、相应责任及致谢内容，再替换占位句；不能由本次审稿推定个人职责或代替最终作者批准。本项属于当前阅读稿的明确待办，不据此声称存在署名不当。
