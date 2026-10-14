# 阶段10i：LOTO v1.3执行进度快照

日期：2026-10-14  
状态：正式冻结运行进行中；本报告仅记录已完成任务、拟合/OOF计数、检查点状态和资源摘要，不查看或解释预测表现。

总计划：[`RESEARCH_PLAN.md`](../RESEARCH_PLAN.md) v3.18

## 冻结执行身份

- 运行目录：`runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_3`
- 执行源码提交：`fedc745940374b6870cb6c4bcb66ec8f731ded27`
- 锁SHA-256：`709a5499f87fdd75d4fd202041c9200a08d068938c7190e3a915dbf65edb954e`
- 固定全预算：1,472任务、40,112 fit events、39,736条OOF预测

## 执行快照

对前12个完整组织区块（组织索引0–11）的最终JSON检查点进行只读盘点：768/1,472项完成，占52.2%；累计20,928个fit events和23,856条OOF。Biliary Tract、Bladder/Urinary Tract、Bone、Bowel、Breast、CNS/Brain、Cervix、Esophagus/Stomach、Head and Neck、Kidney、Liver和Lung各64项全部完成。768个最终检查点状态均为`completed`，且资源字段完整。控制台随后确认已进入Lymphoid区块，达到至少780/1,472项；正式全run验证尚未开始。

已完成检查点中的最大采样任务RSS/私有内存为1.101/0.877 GiB；任务内采样的最低可用RAM/commit余量为0.494/0.306 GiB。冻结配置的**下一任务启动门槛**分别为0.500/1.000 GiB；该两项是运行中任务资源采样的低点，不表示改变门槛或算法。控制器仍以单worker串行调度，截至该快照未报告资源暂停或任务失败；后续继续监测门槛等待与错误事件。checkpoint data guard均未访问外部响应（`external_response_accessed=false`），阻止访问次数为0。该盘点只读取任务身份、状态、fit/OOF条数、data guard和资源字段，没有查看预测值、标签、指标、模型参数、零解率或显著性。CTRP继续封存。

## 后续验收顺序

保持锁定执行版本、输入、任务顺序、统计方案和单worker资源策略，继续完成剩余任务，不按中途表现改变任何规则。全部任务完成后，先独立验证执行锁和输入/来源哈希、任务/拟合/OOF覆盖、种子、checkpoint payload、失败状态及资源事件；验证通过后再运行冻结分析器，并严格保留NA、覆盖分母、零解、失败和警告。内部LOTO结果验收后，才冻结CTRP的missing/abstain、coverage及一次性统计协议；此前不访问CTRP响应。
