# 阶段10f：LOTO v1.3执行进度快照

日期：2026-10-14

状态：正式冻结运行进行中；本报告不包含性能解释。

总计划：[`RESEARCH_PLAN.md`](../RESEARCH_PLAN.md) v3.15

## 执行标识

- 运行目录：`runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_3`
- 执行源码提交：`fedc745940374b6870cb6c4bcb66ec8f731ded27`
- v1.3锁SHA-256：`709a5499f87fdd75d4fd202041c9200a08d068938c7190e3a915dbf65edb954e`
- 原定全预算：1,472任务、40,112 fit events、39,736条OOF预测

## 当前任务覆盖

只读取已完成最终JSON检查点中的任务标识、状态、行数和资源字段后，快照为450/1,472项已完成，占30.6%；累计12,240个fit events和9,652条OOF。Biliary Tract、Bladder/Urinary Tract、Bone、Bowel、Breast、CNS/Brain和Cervix七个组织均完成64项；Esophagus/Stomach已完成2/64项。450个现有检查点状态均为`completed`。控制器仍单进程顺序运行，正式全运行独立验证尚未启动。

## 资源和数据边界

截至已完成检查点，任务采样最大RSS为1.101 GiB，最大私有内存0.877 GiB；最小可用RAM/commit余量分别为0.758/0.637 GiB。资源摘要只描述已经发生的采样，不是保证性上限。已完成检查点的数据守卫中没有外部响应文件访问记录；CTRP仍封存。

本快照不打开或汇总模型预测值、响应标签、rho/RMSE、模型系数、零解比例或任何药物效能结果，不产生效应区间或显著性判断。所有已有任务计数都属于原1,472任务预算，没有新增任务或改变任务选择。

## 后续顺序

继续完成剩余1,022项。完成后先核对全run的锁、输入与源码哈希，任务/fit/OOF总数、任务顺序、内外层种子、checkpoint payload和错误状态；之后再运行冻结分析器，保留strict NA、组织覆盖分母、零解及失败/警告记录。内部验收前不准备/读取CTRP响应文件；内部结果验收后另行冻结外验missing/abstain政策。
