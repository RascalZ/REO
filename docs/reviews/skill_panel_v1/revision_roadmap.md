# 修订路线图：不可变来源顺序 core

本路线图只综合五份已完成的审查卡，不改写原稿。总体机械决定为 **Major Revision**，源卡六条弱点仍全部为 **Minor**；EIC-W3 与 R1-W1 合并为同一个访问项后共五项：must_fix 2 项，should_fix 3 项。总体标签由冻结合同 F3 决定，不按问题数量、confidence 或本路线图的编辑门槛重新推导。

`schema_version: revision-roadmap/1.0`；`revision_round: 1`；`calibration_status: NOT_CALIBRATED`；`criteria_binding_unavailable`。

固定排序键为 EIC、R1、R2、R3、DA，再按 finding ordinal、subclaim ordinal 和 item id。以下顺序只表示来源可追溯性。R/S 是门槛分类内的运输引用，不是工作排名。当前没有作者选择 sidecar；拟议 block 和操作仅供作者随后选择，尚不构成执行授权。

## 机器绑定与来源

机器 core 为 `revision_roadmap.json`。绑定的原文标记副本为 `local_private/input/paper.anchored.md`，SHA256=`5761134e6c5e56e434665f4e26d1e3955438c35f88651de97b76ab713eba2c3e`；block manifest 为 `local_private/input/paper.block-manifest.json`，SHA256=`9c145bd72e6d32c9d05f794123028f29b707dc0feb08673ec0628e61434bb815`。标记副本仅用于精确定位；输入方已核对去标记后与冻结原文逐字节一致，原稿未改动。

实际源卡为 `local_private/phase2/{EIC,R1,R2,R3,DA}.md`，同步副本为 `review_reports/`。五席预承诺和实际 Phase 1→Phase 2 检查分别见 `precommitments/` 与 `phase_conformance.json`；机器结构检查不认证科学结论。每项证据锚、等级和 confidence 均由源 finding 搬运，不另设严重程度。confidence 是未校准的自报范围信息，不改变门槛、排序或机械决定。

## 完整来源顺序

| 顺序／ID | 来源位置／子问题 | 内容 | 原等级／confidence | 编辑门槛／运输引用 | 拟议 block 操作 |
| --- | --- | --- | --- | --- | --- |
| REV-001 | EIC finding 1／SC-1 | 方法主文的执行修复细节层级 | Minor／5 | should_fix／S1 | B0053、B0063：replace_block |
| REV-002 | EIC finding 2／SC-2 | LOTO 灰格实际未定义原因的编辑查询 | Minor／4 | should_fix／S2 | B0084、B0147：replace_block |
| REV-003 | EIC finding 3 + R1 finding 1／SC-3 | 代码及关键派生资产的实际访问安排 | 两源均 Minor／5 | must_fix／R1 | B0064、B0115、B0117：replace_block |
| REV-004 | EIC finding 4／SC-4 | 作者贡献与致谢占位声明 | Minor／5 | must_fix／R2 | B0123、B0125：replace_block |
| REV-005 | R2 finding 1／SC-5 | 两屏 AUC 实验、积分、归一化及方向定义 | Minor／5 | should_fix／S3 | B0032：insert_after |

访问项是两席相互印证（2/4），其余为单席意见；其他非 DA 席位的沉默不计为赞同。没有 CONSENSUS-3/4 弱点或 SPLIT。闭合 schema 不含“2/4 corroborated”枚举，因此 REV-003 由 `SINGLE-VERIFIER` 保守承载，不表示实际只有一席提出；实际支持数量和两套证据保存在 `source_refs`、`corroborating_sources` 及 `consensus_summary`。DA 的 Critical/Major 表为空，未另造 DA 修订项。

## REV-001／S1：压缩主文执行日志细节

来源：EIC-W1。原等级 **Minor**；confidence **5**，依据为“主文与补充 SM3 对读的编辑结构判断”。

**搬运证据锚 — text**：主文 Methods—Leave-one-tissue-out sensitivity analysis 第2段：“Before fitting, execution amendment v1.1 resolved a prose/config seed conflict in favor of the frozen machine value 2026101401.” 源意见同时指出 Reproducibility and data governance 中的三个长科学提交哈希，以及现有 SM3、SM6、S13 的对应追溯信息。

最低建议：主文保留训练边界、预先冻结、实际科学版本及修复未改变分析选择的一两句概述；将具体修复序列、逐项种子和长哈希集中指向已有补充小节。影响科学选择的修订披露须保留，无需重新执行分析。

- **Acceptance criteria**: 主文仍可直接辨认训练边界、冻结及实际科学版本和修复未改分析选择；具体序列／种子／哈希通过准确补充交叉引用可查；不改变科学结果或执行新拟合。

`obligation_class=should_fix`；`cost_scope=section`：LOTO 与 Reproducibility 方法段及现有补充交叉引用；`consequence_if_unaddressed=interpretive_ambiguity_remains`，目标为 Methods 核心流程与执行修复层级。拟议范围仅为 B0053、B0063 的 `replace_block`。

## REV-002／S2：核对并解释 LOTO 灰格含义

来源：EIC-W2。原等级 **Minor**；confidence **4**，依据为“对现有图注可读性的直接判断，实际原因分布尚需作者核对”。这是一项待验证的编辑查询，未确认为统计错误。

**搬运证据锚 — figure**：Figure 4A，标为 Undefined Δρ 的灰色单元、各药的 defined-count 行标与 34/80 注记。

最低建议：核对已存在 S1/S7 状态字段，在图注或结果段简述十组织网格实际发生的未定义原因及范围；必要时用小表分列。若尚未分类，明确说明未分类。继续保留严格宏值 NA，不推定全部灰格由零分子解造成，不把 NA 填零或以可用格均值代替预定宏值；无需改变敏感性分析或补完整跨组织验证。

- **Acceptance criteria**: 补充解释逐项可追溯至现有 S1/S7 状态记录并标明范围；无法分类者如实标注；正文与图注仍保留 34/80 和严格宏值 NA；未造原因、未换统计量。

`obligation_class=should_fix`；`cost_scope=section`：Figure 4 图注或 LOTO Results 的说明；`consequence_if_unaddressed=interpretive_ambiguity_remains`，目标为 Figure 4A 灰格含义。拟议范围为 B0084、B0147 的 `replace_block`，具体采用哪处须由后续作者选择。

## REV-003／R1：提交前落实访问／重建路径

来源：EIC-W3 与 R1-W1，合并同一访问安排问题。两份源 finding 均为 **Minor**、confidence **5**。EIC 依据为“稿件明确披露的当前访问状态”；R1 依据为“数据与代码声明直接确认当前访问状态”。

**搬运证据锚 — EIC text**：主文 Code availability：“Public or reviewer access and an archival release remain to be arranged before submission.”

**搬运证据锚 — R1 text**：Declarations—Availability of data and materials 第2段：“They are not currently included in the version-controlled repository”；Code availability 第1段：“Public or reviewer access and an archival release remain to be arranged before submission”。

最低建议：提交前为精确代码版本、身份／分组衍生规则及原始折定义提供审稿人可用入口或确定重建步骤；逐项交代固定模型、系数和预处理状态的许可访问或重建方式、获取条件、限制、版本及责任联系人。核对引用版本与实际内容一致，再同步相关声明。允许受上游许可约束的受控访问；不要求重新分发全部原始矩阵、取得某一平台 DOI 或重拟合模型。

- **Acceptance criteria**: 各类最低必需材料均有可执行入口或可追踪重建步骤、精确版本、联系人与限制；按所声明权限可取得相应内容或完成所列获取步骤；声明及引用版本与实际一致，受限内容有明确替代路径；落实前不声称公开完成或独立复现。

`obligation_class=must_fix`；`cost_scope=other/access_arrangements`：访问和声明安排；`consequence_if_unaddressed=method_reproducibility_unresolved`，目标为关键代码／派生资产的实际可获得性。拟议范围为 B0064、B0115、B0117 的 `replace_block`。这些是已延后的提交前收尾，不反证主比较，也不构成数据泄漏或数值失真的指控。

可选增强保留源卡的范围：带清单、校验和、输入来源及许可边界的归档包，或第三方按说明重建主要图表。若未来选择从源数据重训练，需要另立新拟合提案；它不在本项最低修复中，本轮未执行。

## REV-004／R2：由作者完成贡献与致谢声明

来源：EIC-W4。原等级 **Minor**；confidence **5**，依据为“声明中的明确占位表述”。

**搬运证据锚 — text**：主文 Authors' contributions：“Author-specific contribution roles will be completed by the listed authors before submission.”；Acknowledgements：“The acknowledgement statement will be finalized by the authors before submission.”

最低建议：所列作者确认具体贡献、相应责任和致谢内容后，再替换两处占位句。审查者不推断个人分工、填写作者事实或代替作者最终批准。本项是已披露的投稿准备工作，不据此指控署名不当。

- **Acceptance criteria**: 两处未来完成的占位句均被作者确认的实际声明替换，作者确认有可追溯记录；具体内容由作者提供及核实；当前 core 不记录或假定确认已经发生。

`obligation_class=must_fix`；`cost_scope=other/author_declarations`：作者确认和两处声明；`consequence_if_unaddressed=reporting_requirement_unmet`，目标为贡献与致谢小节。拟议范围为 B0123、B0125 的 `replace_block`。

## REV-005／S3：补两屏归档终点对照

来源：R2-W1。原等级 **Minor**；confidence **5**，依据为“已逐项核查所列终点说明与支持表字段”。当前是说明不足，未证实终点方向或数据使用错误。

**搬运证据锚 — absence**：定位为 Methods 的 Data sources、Response outcome 和 Full-development models 小节及 Table S8；预期内容为“两屏所用归档终点的读出与暴露时长、AUC 积分自变量和归一化定义，以及 CTRP 方向说明”；源卡检查范围为“完整主文、Additional file 14 的 SM1–SM8、Table S8 全部字段及 Figure 1/4 图注”。

最低建议：依据实际使用的归档与原始资源，增加简短两屏终点对照，列出归档文件及原字段名、读出、暴露时长、积分轴、积分区间／归一化约定、数值方向及依据。CTRP 剂量下限或积分细节无法核验时明确标为未重建，不补未经核对的数字。源卡引有 GDSC 官方说明以支持可核验信息；本综合没有开展新的外部核验。

- **Acceptance criteria**: 两屏各项信息及来源齐列，无法从实际归档核验的字段明确标为未重建；新增说明不冒称方向错误已确认，也不把药内排序改成绝对 AUC 校准或实验条件的因果解释。

`obligation_class=should_fix`；`cost_scope=section`：Methods 简短终点对照；`consequence_if_unaddressed=interpretive_ambiguity_remains`，目标为 GDSC2／CTRPv2 AUC 测量及方向定义。拟议范围为 B0032 的 `insert_after`。无需训练、评分或重抽样。

源卡的可选增强为统一剂量窗或统一曲线处理比较；这会改变冻结评价问题，涉及原始曲线及另行范围设计，不是本项最低条件。

## 来源完整性与作者回应

六条原始弱点已全部进入上述五项；R3 未列 D4 弱点，DA 未列 Critical/Major，未以数量配额新增条目。已披露且结论未越界的条件区间、LOTO、共同表达源、有限合成控制和原环境记录限制继续保留，没有转化为新增研究重做任务。未完成的全部题录核验、个体标签／预测和原始训练复现仍是本轮核验范围之外。

后续作者应对每个 REV 项明确选择及说明，作者选择须另建并与本 core 哈希绑定的 sidecar；本路线图不预填接受、拒绝、显示／执行次序或旁及修改授权。回应可采用 skill 的 `templates/revision_response_template.md` 逐项记录，并引用对应可检验最低修复条件。
