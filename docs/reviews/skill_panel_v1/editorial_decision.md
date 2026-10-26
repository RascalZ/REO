# 编辑综合决定

稿件：*Evaluating the incremental value of relative expression ordering for drug-response prediction across eight compounds*。本轮依据冻结的 `reviewer/reviewer_full/v2` 合同、五份已完成的审查卡和执行来源记录进行综合，日期为 2026-10-26。原稿未改动。

`calibration_status: NOT_CALIBRATED`

**总体机械决定：Major Revision。** 这是冻结合同中 F3 的执行结果；六条来源弱点的单项等级均为 **Minor**，本综合没有将其中任何一条升级为 Major，也没有提出重训练或新增队列的最低要求。

五席均披露 `criteria_binding_unavailable`。本轮为明确 unbound 的审查；作者请求的 BMC Genomics Research article 及所给官方政策仅作背景。D6 的 pass 表示 EIC 在本合同范围内辨认到有限经验贡献，不构成正式期刊标准认证、真实接收概率或文献首创认证。

## 一、六维角色范围与机械综合

只纳入合同规定的 eligible role 且实际评估的分数。角色外 `not_assessed` 和弃权不进入分子或分母。六维均有合格席位作出判断，无未评估维度。

| 维度 | 优先级 | 合格角色／实际席位 | 纳入的分数 | 审计判定 |
| --- | --- | --- | --- | --- |
| D1 methodology_rigor | mandatory | methodology／R1 | R1=warn | warn |
| D2 domain_accuracy | mandatory | domain／R2 | R2=warn | warn |
| D3 argumentative_coherence | mandatory | da、methodology／DA、R1 | DA=pass；R1=pass | pass |
| D4 cross_disciplinary_relevance | high | perspective／R3 | R3=pass | pass |
| D5 writing_and_structure | normal | eic／EIC | EIC=warn | warn |
| D6 venue_fit_and_contribution | mandatory | eic／EIC | EIC=pass | pass |

各维审计判定为合格且已评估分数中的最差值；本轮无 block、fatal block。条件求值仍逐维应用各条件的 `cross_reviewer_quantifier`，随后应用条件的维度量词，未以维度审计表代替第二阶段求值。

| 条件 | 合同条件及跨席量词 | fired | 条件 severity | 合同 action |
| --- | --- | --- | --- | --- |
| F1 | 任一 mandatory 维度有 fatal block；any | false | 95 | reject |
| F2 | 任一 mandatory 维度为 block；any | false | 90 | major_revision |
| F3 | 至少两个 mandatory 维度为 warn 或更差；majority | true | 70 | major_revision |
| F4 | 任一 high 维度为 block；any | false | 60 | major_revision |
| F5 | 任一维度为 warn 或更差；any | true | 40 | minor_revision |
| F0 | 每个维度均为 pass；all | false | 10 | accept |

F3 中，D1、D2 均只有一个合格已评估席位，majority 由各自 owner 的 warn 满足；D3 有两个合格席位，须两者均满足才算该维命中，而本轮两者均为 pass；D6 为 pass。因此恰有两个 mandatory 维度命中。F3 的 severity=70 高于 F5 的 40，按冻结优先规则选择 F3。

dimension_verdicts: [D1=warn, D2=warn, D3=pass, D4=pass, D5=warn, D6=pass]
fired_conditions: [F3, F5]
da_critical_adjudications: []
editorial_decision=major_revision

这里有四种不同语义：维度分数用于条件求值；条件的数值 severity 仅决定触发条件之间的优先级；单项 Critical/Major/Minor 从源卡原样搬运；路线图的 must_fix/should_fix 是编辑完成门槛。它们不相互换算。自报 confidence 只表示审查范围与不确定性，未用于加权、投票或覆盖异议。总体 Major Revision 不能反推存在单项 Major 科学错误，单项 Minor 也不能覆盖 F3 所得总体决定。

## 二、已有意见的综合与保留

EIC、R1、R2、R3 及 DA 从各自范围认可：主问题限于固定八药、指定候选规则及学习程序；实际强表达流程比较与同候选集的表示比较得到区分；近零差值未被写成等效；零分子解仍可保留组织预测；共同 CCLE 表达来源下的跨屏排序没有被提升为跨表达平台或临床验证。对应证据分别见 EIC S1–S4、R1 S1–S5、R2 S1–S4、R3 S1–S6 和 DA 的反论证记录。这些判断仅限各卡所述阅读和核对范围。

以下六条源弱点合并为五项。EIC-W3 与 R1-W1 是相同访问安排问题，保留两份来源和各自证据。正式路线图按 EIC、R1、R2、R3、DA，再按来源 finding ordinal 排列；顺序表示溯源，不表示工作优先级。

| 子问题 | 来源、原等级、自报 confidence | 综合处置与路线图 |
| --- | --- | --- |
| SC-1：执行日志细节与长哈希占用主文的方法层级 | EIC-W1；Minor；5 | 单席意见，should_fix；REV-001。保留科学边界和实际版本概述，以现有补充小节承载具体修复细节。 |
| SC-2：LOTO 灰格的实际未定义原因未集中说明 | EIC-W2；Minor；4 | 单席编辑查询，should_fix；REV-002。先核对现有 S1/S7 状态，再简述实际原因；未确认原因前保留未分类状态。 |
| SC-3：代码及关键派生资产缺实际访问安排 | EIC-W3、R1-W1；两者均 Minor、5 | 两席相互印证（2/4），不达 CONSENSUS-3/4；must_fix；REV-003。提交前落实可执行访问或重建路径。 |
| SC-4：作者贡献与致谢仍为占位声明 | EIC-W4；Minor；5 | 单席意见，must_fix；REV-004。由作者确认后替换占位内容。 |
| SC-5：两屏 AUC 终点的实验、积分与方向定义仍不集中完整 | R2-W1；Minor；5 | 单席意见，should_fix；REV-005。按实际归档和原始资源补终点对照，无法核验部分明确标注。 |

这些弱点未出现需要仲裁的存在、严重程度或修复方向冲突。未提到某条弱点的席位保持沉默，未计为赞同。R1 对 34/80 定义域及严格宏值 NA 的一致性肯定，与 EIC 要求解释灰格原因的局部查询相容。R3 的 D4 pass 和 DA 的 D3 pass 也未否定角色外的编辑及访问事项。本轮没有 CONSENSUS-3/4 弱点或 SPLIT；2/4 访问项由闭合 schema 的 `SINGLE-VERIFIER` 枚举保守承载，不表示实际只有一席提出。其两席来源完整保存于 `source_refs` 与 `corroborating_sources`，不借此改写支持数量。

### 科学报告不足及编辑查询

R2-W1 是已确认的局部终点说明不足：现稿缺少可集中读取的两屏读出、暴露时长、AUC 积分轴及区间／归一化约定，CTRP 的具体方向也须与实际所用字段核对。来源没有确认方向错误或错误使用数据；最低修复是来源核对和简短报告补充。EIC-W2 需要依据已有状态记录回答实际 NA 成因，是待验证的图示查询，不能预先写成统计错误，不能把全部灰格归因于零分子解。EIC-W1 属于主文结构与可读性调整。三项均保留 Minor。

### 明确延后的投稿准备事项

REV-003、REV-004 是提交前完成门槛，与核心结论真伪分开。访问条款应逐项区分代码、身份／分组衍生规则、折定义、固定模型和预处理状态，并写明版本、入口或重建步骤、限制及责任联系人；允许受上游条款约束的受控访问。特定平台、DOI、重新分发所有上游矩阵或第三方重训练不是最低修复条件。作者贡献、相应责任与致谢内容须由所列作者核实，本综合不填写个人分工、推断许可事实或代替作者批准。

### 持续保留的范围边界

五卡没有把已披露且未被结论越过的局限重列为缺陷：开发比较受先前探索影响；固定预测的条件区间未包含重新选特征、调参及拟合的不确定性；LOTO 仅 34/80 单元有定义且严格宏值为 NA；外验共用 CCLE 表达源；有限合成控制及其告警、环境记录只支持限定行为。当前意见不要求为这些已限定结论增加大规模训练、真实标签全流程置换、患者队列或独立表达平台。R1 的完整归档／第三方重建和 R2 的统一剂量窗比较保留为可选增强，若未来选择，应另行明确研究范围与执行授权。

## 三、DA Critical 裁决

DA 的 CRITICAL 和 MAJOR 终末表均为空，未报告 C1…Cn。本轮无需构造 Critical 裁决或作者回应条目，故机器审计行保留空数组。DA 已检查近零与等效、候选筛选和模型几何、调参目标、探索性开发、条件区间、组织零解、LOTO 与外推边界；其 pass 限于现稿实际命题，不能转化为个体数据、原始训练或所有参考文献已被复验的声明。

## 四、提交前 must_fix 明细（仅按来源顺序）

| 运输引用 | 事项 | 来源与原等级 | 搬运的证据锚 | 路线图 |
| --- | --- | --- | --- | --- |
| R1 | 落实代码及关键派生资产的审稿访问／重建路径 | EIC-W3、R1-W1；Minor | text：Code availability “Public or reviewer access and an archival release remain to be arranged before submission.”；R1 另引 Availability 第2段及 Code 第1段 | REV-003 |
| R2 | 作者确认贡献与致谢并替换占位句 | EIC-W4；Minor | text：Authors' contributions “Author-specific contribution roles will be completed by the listed authors before submission.”；Acknowledgements “The acknowledgement statement will be finalized by the authors before submission.” | REV-004 |

R1/R2 为运输引用，不是工作排名。完整五项、所有原证据、自报 confidence、bounded cost/consequence、精确拟议 block 范围及可检验最低修复条件见 `revision_roadmap.md` 和闭合 core `revision_roadmap.json`。当前没有作者选择 sidecar，任何 `proposed_targets` 都尚非改稿授权。

## Review Panel Provenance (#540/#740)

来源为 `panel-provenance.json`，carrier 为 `panel-carrier.json`（status=valid）。执行方已实际运行 provenance build、validate、build-carrier、validate-carrier 并获 exit 0；本文逐轴搬运记录，不从角色名推断独立性。

| 轴 | 值 | 解释范围 |
| --- | --- | --- |
| role_separated | true | 五席分别承担 eic、methodology、domain、perspective、da。 |
| fresh_context | true | 固定范围为 `within_panel_attempt_only`；不据此证明重试或后续轮次相对全部历史均使用新上下文。 |
| blind_to_peer_outputs | true | 席位记录均为 `peer_outputs_visible=false`。 |
| human_distinct | false | 五席 actor_type 均为 model，human_reviewer_id 为空。 |
| model_family_distinct | unknown | model_family 均未记录，无法确认跨模型家族分离。 |
| provider_distinct | unknown | provider 均未记录，无法确认跨提供商分离。 |

Model-family provenance is incomplete; cross-family separation is unknown, so correlated-error risk cannot be ruled out.

据此，模型家族来源不完整，相关错误风险仍在；角色分离、新上下文和互盲不构成独立错误过程保证。本综合阶段在五卡完成后读取全部意见，不属于五席互盲声明的范围。没有计算跨家族汇总或任何数值独立性指标，也未执行未获同意的外部跨模型传输。

绑定信息：contract SHA256=`e9712090d2469fea15a37b8e22d4e137afbcb2bf38d5789939c5df56738ef7af`；normalized manifest SHA256=`0a55bed6f6d43a70e7a7750b1092d8923e134bcdab61937263ed7f9314136f69`；execution topology SHA256=`0eff4322b1af13719a8a18769f247a0023e9d9c1d60aefa34ae778012c8886ee`；provenance artifact raw SHA256=`ca79984498fa233b640506b94f71429e710cbc51ff91cf2a509f193fae15ebd6`。`precommitment_audit.json` 记录了安装包 full.json 的 CRLF 与 LF pin 差异及隔离副本仅还原 LF 后精确匹配的处理，未改合同字段或校验器。`phase_conformance.json` 记录五席真实 Phase 1→Phase 2 检查通过及中文 trigger 词数提示；这些是程序／结构核对，不能认证科学正确。

## 来源及验证限度

源卡固定保存在 `review_reports/EIC.md`、`review_reports/R1.md`、`review_reports/R2.md`、`review_reports/R3.md`、`review_reports/DA.md`，预承诺位于 `precommitments/`，本次综合输入为对应 `local_private/phase2/` 卡。输入版本与资产散列见 `intake_manifest.json`。R1 的数值检查限于所列汇总表的计数、相减及描述统计；R2 的公开来源阅读有明确的全文取回边界。本轮没有完成全部 31 项参考文献的逐条原文验证，没有取得并审计个体标签／预测、全部身份关系、原始日志或固定模型，也没有独立重算区间端点或复现原始训练。`NOT_CALIBRATED` 持续有效，检查器通过不改变上述科学核验边界。
