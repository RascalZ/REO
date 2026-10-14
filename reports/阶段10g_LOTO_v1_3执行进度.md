# 阶段10g：LOTO v1.3执行进度快照

日期：2026-10-14

状态：正式冻结运行进行中；仅汇报执行进度与资源，不解释性能。

总计划：[`RESEARCH_PLAN.md`](../RESEARCH_PLAN.md) v3.16

## 冻结执行身份

- 运行目录：`runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_3`
- 执行源码提交：`fedc745940374b6870cb6c4bcb66ec8f731ded27`
- 锁SHA-256：`709a5499f87fdd75d4fd202041c9200a08d068938c7190e3a915dbf65edb954e`
- 全预算：1,472任务、40,112 fit events、39,736条OOF预测

## 当前覆盖快照

只读检查已完成任务的最终检查点后，计得513/1,472项完成（34.9%），累计13,953个fit events和12,472条OOF。Biliary Tract、Bladder/Urinary Tract、Bone、Bowel、Breast、CNS/Brain、Cervix和Esophagus/Stomach八个组织各64项完成；Head and Neck已完成1/64项。513个已存在检查点的状态均为`completed`；控制器仍串行执行，尚未做全运行独立完整性验收。

检查点资源字段显示最大采样任务RSS/私有内存为1.101/0.877 GiB；已观测最低可用RAM/commit余量为0.758/0.637 GiB。统计仅从资源、状态、任务标识、fit计数和OOF行数得到，不读取模型分数或响应标签。控制器记录无CTRP响应访问。

## 继续执行与验收

1. 继续原v1.3任务计划和检查点，不缩减区块或模型臂。
2. 全部任务结束后，先核对锁、输入/源码哈希、任务集与顺序、种子、fit事件/OOF总数和checkpoint payload，再运行独立结果验证器。
3. 验证通过后按冻结规则汇总LOTO，保留无定义指标、覆盖分母、零解、模型失败与资源事件；内部结果验收后才冻结CTRP评价所需的missing/abstain与coverage政策。

LOTO执行期间CTRP保持封存，外部结果不用于调参。
