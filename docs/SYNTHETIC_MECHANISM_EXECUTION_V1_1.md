# 有限合成机制验证执行协议 v1.1

日期：2026-10-14。该执行修订在任何合成模型拟合前冻结，承接[P5e v1总协议](TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md)，保留原始合成方案和锁[`synthetic_mechanisms_v1.json`](../config/synthetic_mechanisms_v1.json)、[`p5e_protocol_v1.lock.json`](../config/p5e_protocol_v1.lock.json)。仅解决随机噪声支持及分折种子定义不足，不读取真实响应、真实表达或CTRP，不改真实模型、候选规则、表示、网格、种子总数或拟合预算。

## 冻结的输入与随机数

配置为[`synthetic_mechanisms_v1_1.json`](../config/synthetic_mechanisms_v1_1.json)。共有4种数据生成情形，每种10个独立生成队列；每队列240个样本、8个组织各30个样本、220个基因。每个样本对应唯一synthetic biological group；stable gene IDs为`SIMGENE_0001`至`SIMGENE_0220`。响应连续且截断至[0.01, 0.99]。TPM由`2**log2TPM - 1`变换；值必须有限且非负。

每队列用`numpy.random.Generator(PCG64(seed))`，种子固定为`2026101400 + 100*scenario_index + replicate_index`，情形顺序及编号为`0=tissue_only, 1=weak_molecular, 2=common_sample_amplitude, 3=predefined_pair_order`，replicate编号0–9。所有标准正态随机量均通过拒绝采样限制在[-3,3]：先按C行优先顺序生成标准正态数组，按扁平索引升序重抽越界元素，直至数组都在界内。该截断是预先确定、与标签无关的极小尾部支持修订，使第三情形中`log2TPM`保持正值；不可根据运行结果更改。

## 数据生成方程

令组织编号`t=0,...,7`，`h=(t-3.5)/3.5`；基因编号`j=0,...,219`，`load_j`为[-1,1]的220个等距值；样本噪声`epsilon`及响应噪声`eta`均为上述截断标准正态。

1. **纯组织信号**：`y=0.50+0.28h`；`log2TPM_ij=3+0.008j+0.20*h_i*load_j+0.35*epsilon_ij`。样本噪声不进入响应。
2. **弱分子信号**：`log2TPM_ij=3+0.008j+0.15*h_i*load_j+0.35*epsilon_ij`；`z_i=mean(epsilon_i0,...,epsilon_i4)*sqrt(5)`；`y_i=0.50+0.05h_i+0.04z_i+0.14eta_i`。
3. **共同样本幅度**：`a_i`为截断标准正态；`log2TPM_ij=3+0.008j+0.20a_i`；`y_i=0.50+0.04h_i+0.12a_i+0.02eta_i`。任何样本内的全部220基因乘同一正比例因子；精确验证rank与完整pair编码不因该转换改变。
4. **预设pair顺序**：先生成`epsilon`矩阵；`log2TPM_i0=3.5+0.35epsilon_i0`、`log2TPM_i1=3.5+0.35epsilon_i1`，其余基因为`3+0.008j+epsilon_ij`；`y_i=0.50+0.04h_i+0.15sign(log2TPM_i0-log2TPM_i1)+0.03eta_i`。

前3种情形的响应方程按[原冻结协议](TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md#b-有限合成统计机制检验)执行；上方完整重述避免DGP解释依赖外部记忆。只对`y`按[0.01,0.99]截断。生成器在所有表达转换后断言shape、finite、TPM非负、响应范围及8×30组织构成。

## 嵌套评价与任务预算

使用原冻结的7个分子模型及`expression_strong`、同政策`tissue_null`共8臂；每臂使用相同训练内特征选择和200/2000基因预算、三种表达/秩/完整三值pair表示、相对alpha-max EN网格、内层RMSE选择、零端点和组织FWL政策。仅对表达输入规则作必要的配置适配；不拟合随机森林。

外层为`StratifiedGroupKFold(3, shuffle=True)`，seed=`cohort_seed+1000`，在240个唯一synthetic group上按8组织分层；内层同为3折grouped/tissue-stratified，限制于各外层训练样本，seed=`cohort_seed+1100+outer_fold`。所有八臂共享同队列分折。每个分子外层task做3×10个内层候选fit及1个外层refit，共31 fits；组织null每task 1 fit。因此40队列×3外折×(7×31+1)=26,160 fits、960个独立task、76,800条样本×模型外层预测。按原始协议，完备性/信息流是执行验收项，预测胜负不是。

## 冻结检查与解释界限

运行必须保留所有失败、零解、回退、警告、参数、OOF、训练成员审计、耗时及采样RSS/系统内存；每个task独立进程、单worker、BLAS单线程、可原子恢复。公开Git只同步代码、锁、聚合表与报告；合成值、样本级预测和完整检查点留在忽略的本地`runs/`。

分析报告按DGP×模型给出每个队列的完整OOF统计，并汇总10队列的均值/中位数/范围、同政策组织参照的组织以上排序相关`Spearman(y-null, prediction-null)`、RMSE改善、分子零解与候选基因稳定性。报告driver基因是否进入训练选择，但不设性能或选择频率成功线；不给生物学p值或推断CI。共同幅度情形要求逐样本秩及完整pair严格不变；纯组织情形记录组织null及分子解；预设顺序情形仅作为已知DGP的表征压力检查。

以上均为假设明确的有限合成控制，只能核对实现对特定数据生成机制的行为；不证明真实药物反应的生成机制，不补足固定八药证据，也不确证REO更优。完成后按总计划进入23组织LOTO；CTRP继续封存。
