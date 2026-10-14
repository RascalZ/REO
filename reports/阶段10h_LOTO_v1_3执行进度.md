# 阶段10h：LOTO v1.3执行进度快照

日期：2026-10-14

状态：正式冻结运行进行中；本报告仅记录执行状态和资源摘要，不解释预测表现。

总计划：[`RESEARCH_PLAN.md`](../RESEARCH_PLAN.md) v3.17

## 冻结执行身份

- 运行目录：`runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_3`
- 执行源码提交：`fedc745940374b6870cb6c4bcb66ec8f731ded27`
- 锁SHA-256：`709a5499f87fdd75d4fd202041c9200a08d068938c7190e3a915dbf65edb954e`
- 固定全预算：1,472任务、40,112 fit events、39,736条OOF预测

## 执行快照

对已完成的最终JSON检查点进行只读盘点：644/1,472项完成，占43.8%；累计17,534个fit events和14,684条OOF。Biliary Tract、Bladder/Urinary Tract、Bone、Bowel、Breast、CNS/Brain、Cervix、Esophagus/Stomach、Head and Neck、Kidney十个组织均完成64项；Liver已完成4/64项。644个已发现检查点状态均为`completed`，控制器仍按锁定顺序串行运行，尚未执行全run独立验证。

已完成检查点中的最大采样任务RSS/私有内存为1.101/0.877 GiB，最小可用RAM/commit余量为0.758/0.637 GiB。这些数字只覆盖已完成的任务采样。任务和数据守卫记录未访问CTRP响应；本快照只读状态、行数、任务标识和资源字段，没有查看或汇总预测值、标签、药物性能、零解率、区间或显著性。

## 后续验收顺序

继续完成剩余828项，并保持既定任务顺序与资源政策。全部任务完成后，独立验证锁与输入哈希、来源、任务/拟合/OOF覆盖、种子、checkpoint payload和错误状态；通过后再运行冻结分析器，严格保留NA、覆盖分母、零解、失败/警告。内部LOTO结果验收后，才冻结CTRP的missing/abstain、coverage及一次性统计协议；在此之前不访问外验响应。
