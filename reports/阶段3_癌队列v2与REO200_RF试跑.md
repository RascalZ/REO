# 阶段3：癌模型队列v2与REO200/RF开发集试跑

报告生成UTC：2026-10-12T03:57:48.406732+00:00；运行声明日期：2026-10-12；运行状态：completed。执行代码提交：5d053490f4a4b025245eee3c880fb42827154e2f；报告生成提交：5d053490f4a4b025245eee3c880fb42827154e2f（工作区dirty）。数值均来自完成运行的聚合输出，不读取个体响应。

## 本阶段完成与解释边界

本轮在癌模型队列v2的两种预设化合物上运行7种模型；外层3折×1次、内层3折，共906次fit。REO候选基因预算为200，普通表达主基线预算为2000；匹配表达预算为200。这是扩展工程pilot及内部点估计，不是正式8药5折×3次结果，也不构成可靠增量或等效证明。

本轮主要REO EN在全部6个外层重拟合中选择零分子系数解，预测实质上仅来自训练组织效应。因此，即使其误差比含表达特征的EN略低，也不能将差异称为REO分子信息的增量。辅助RF的点比较应结合下表判断，但不能挽救为主EN已经有效的结论。

CTRP响应访问标志为False，访问防护拦截尝试数为0；外验标签尚未用于拟合、选参或评价。所有模型使用相同全局生物组划分；聚合报告检查N与预设名单一致，精确逐ID预测配对由独立检查验收，报告脚本不重新打开个体OOF表。

## 癌模型队列v2

资格修订仅依据metadata：排除1个开发及5个主外验明确非癌模型，得到626个开发、175个主外验模型。HCC1588/ACH-001078由[Cellosaurus CVCL_A351](https://www.cellosaurus.org/CVCL_A351)及[KCLB71588](https://cellbank.snu.ac.kr/cellline/search?q=71588&sc=y)核实为肺鳞癌模型，精确Sanger标识闭合后纳入。资格修订发生在旧开发pilot之后、CTRP标签开放之前，不能称为旧pilot之前的预注册。v1与旧运行保留，固定8药排序、ID和剂量窗保持。

| 固定主分析化合物 | v2开发N | v2主外验N |
| --- | --- | --- |
| MIRA-1 | 625 | 108 |
| PLX-4720 | 622 | 167 |
| Bortezomib | 620 | 169 |
| Telomerase Inhibitor IX | 620 | 168 |
| Tamoxifen | 620 | 165 |
| BMS-536924 | 620 | 165 |
| Gefitinib | 620 | 161 |
| Cediranib | 620 | 82 |

v2独立资格/身份验证状态为passed；主外验与全部已知GDSC2候选生物组交集为0。只按权威精确casefold词表修正2处组织大小写，保留v1_tissue。队列freeze SHA256：fe1b725aa6acced46e1e5574930bb4c816c4414fc21c08d46cdb553d1fdc2d98。

## 内部OOF模型比较

| 药物 | 模型 | OOF N | Spearman | RMSE | MAE | R² | 流程秒 | 警告fits | 无分子项外折 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MIRA-1 | 训练均值 | 625 | NA | 0.024157 | 0.013647 | -0.0088705 | 0.15185 | 0 | 3 |
| MIRA-1 | 组织Ridge | 625 | 0.076464 | 0.022872 | 0.013407 | 0.095558 | 1.3295 | 0 | 3 |
| MIRA-1 | 表达2000 EN | 625 | 0.099947 | 0.023533 | 0.013991 | 0.042511 | 36.029 | 0 | 2 |
| MIRA-1 | 表达200 EN | 625 | 0.083997 | 0.023038 | 0.013696 | 0.08242 | 30.329 | 0 | 2 |
| MIRA-1 | REO200 EN | 625 | 0.099853 | 0.022837 | 0.013401 | 0.09831 | 97.04 | 0 | 3 |
| MIRA-1 | 表达2000 RF | 625 | 0.14392 | 0.023079 | 0.013229 | 0.079128 | 290.48 | 0 | 0 |
| MIRA-1 | REO200 RF | 625 | 0.045039 | 0.023341 | 0.013593 | 0.058105 | 56.755 | 0 | 0 |
| PLX-4720 | 训练均值 | 622 | NA | 0.066232 | 0.038444 | -0.0068068 | 0.14178 | 0 | 3 |
| PLX-4720 | 组织Ridge | 622 | 0.21906 | 0.053024 | 0.032612 | 0.3547 | 1.2466 | 0 | 3 |
| PLX-4720 | 表达2000 EN | 622 | 0.23363 | 0.053862 | 0.033119 | 0.33415 | 36.507 | 0 | 2 |
| PLX-4720 | 表达200 EN | 622 | 0.22496 | 0.053057 | 0.032618 | 0.35391 | 29.183 | 0 | 3 |
| PLX-4720 | REO200 EN | 622 | 0.22496 | 0.053057 | 0.032618 | 0.35391 | 126.32 | 0 | 3 |
| PLX-4720 | 表达2000 RF | 622 | 0.32233 | 0.052258 | 0.032797 | 0.37322 | 332.83 | 0 | 0 |
| PLX-4720 | REO200 RF | 622 | 0.30544 | 0.053059 | 0.033502 | 0.35386 | 61.051 | 0 | 0 |

主要比较为REO EN相对表达EN2000；RF比较为辅助、同模型家族内的整条pipeline比较。每药单独报告OOF Spearman及原始AUC尺度误差；不合并药物计算一个总相关。训练均值的Spearman按基线定义记NA；不同折均值拼接后产生的排序不能作为个体预测能力。

| 药物 | 比较 | N | ΔSpearman | 相对RMSE改善% |
| --- | --- | --- | --- | --- |
| MIRA-1 | primary_EN | 625 | -9.3834e-05 | 2.9576 |
| MIRA-1 | secondary_RF | 625 | -0.098881 | -1.135 |
| MIRA-1 | matched_gene_EN_ablation | 625 | 0.015856 | 0.86966 |
| PLX-4720 | primary_EN | 622 | -0.0086763 | 1.4953 |
| PLX-4720 | secondary_RF | 622 | -0.016892 | -1.5325 |
| PLX-4720 | matched_gene_EN_ablation | 622 | 0 | 0 |

ΔSpearman=候选REO−参考表达，正值表示名义方向有利；相对RMSE改善=(参考RMSE−候选RMSE)/参考RMSE。匹配基因数表达模型是消融，不替换2000基因强表达主基线。此处只有点差，无置信区间或确证性P值；3个相关外折不能按3个独立重复作普通t检验，非显著也不证明等效。

## 正则路径、零分子解与组织贡献

每个EN fit独立从自身训练标准化及组织投影后的X/y计算alpha_max，再使用固定相对路径[1.0, 0.3, 0.1, 0.03, 0.01]和L1比例[0.5, 1.0]。训练时组织项采用不受惩罚的FWL投影；REO先由原TPM构造，不能用gene-wise组织残差生成顺序。本轮三种EN的18个外层重拟合中有15个零分子系数解。

| gdsc_drug_id | model | outer_folds | zero_molecular_folds | max_nonzero | max_dual_gap | warnings |
| --- | --- | --- | --- | --- | --- | --- |
| 1036 | expression | 3 | 2 | 124 | 0.00021804 | 0 |
| 1036 | expression_matched | 3 | 3 | 0 | 0 | 0 |
| 1036 | reo | 3 | 3 | 0 | 0 | 0 |
| 1931 | expression | 3 | 2 | 70 | 4.2977e-05 | 0 |
| 1931 | expression_matched | 3 | 2 | 37 | 6.8951e-05 | 0 |
| 1931 | reo | 3 | 3 | 0 | 0 | 0 |

| 药物ID | 模型 | 外折 | alpha比 | 实际alpha | 非零分子系数 | 候选分子特征 | 模型分子特征 | dual gap | 收敛警告 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1931 | expression | 0 | 0.3 | 0.11428 | 70 | 2000 | 2000 | 4.2977e-05 | 0 |
| 1931 | expression_matched | 0 | 1 | 0.2262 | 0 | 200 | 200 | 0 | 0 |
| 1931 | reo | 0 | 1 | 0.36076 | 0 | 19900 | 19900 | 0 | 0 |
| 1931 | expression_rf | 0 | NA | NA | NA | 2000 | 2000 | NA | 0 |
| 1931 | reo_rf | 0 | NA | NA | NA | 19900 | 500 | NA | 0 |
| 1931 | expression | 1 | 1 | 0.30112 | 0 | 2000 | 2000 | 0 | 0 |
| 1931 | expression_matched | 1 | 1 | 0.2188 | 0 | 200 | 200 | 0 | 0 |
| 1931 | reo | 1 | 1 | 0.35604 | 0 | 19900 | 19900 | 0 | 0 |
| 1931 | expression_rf | 1 | NA | NA | NA | 2000 | 2000 | NA | 0 |
| 1931 | reo_rf | 1 | NA | NA | NA | 19900 | 500 | NA | 0 |
| 1931 | expression | 2 | 1 | 0.40495 | 0 | 2000 | 2000 | 0 | 0 |
| 1931 | expression_matched | 2 | 0.3 | 0.067733 | 37 | 200 | 200 | 6.8951e-05 | 0 |
| 1931 | reo | 2 | 1 | 0.38805 | 0 | 19899 | 19899 | 0 | 0 |
| 1931 | expression_rf | 2 | NA | NA | NA | 2000 | 2000 | NA | 0 |
| 1931 | reo_rf | 2 | NA | NA | NA | 19899 | 500 | NA | 0 |
| 1036 | expression | 0 | 0.3 | 0.062559 | 124 | 2000 | 2000 | 0.00021804 | 0 |
| 1036 | expression_matched | 0 | 1 | 0.18883 | 0 | 200 | 200 | 0 | 0 |
| 1036 | reo | 0 | 1 | 0.27736 | 0 | 19900 | 19900 | 0 | 0 |
| 1036 | expression_rf | 0 | NA | NA | NA | 2000 | 2000 | NA | 0 |
| 1036 | reo_rf | 0 | NA | NA | NA | 19900 | 500 | NA | 0 |
| 1036 | expression | 1 | 1 | 0.27778 | 0 | 2000 | 2000 | 0 | 0 |
| 1036 | expression_matched | 1 | 1 | 0.18462 | 0 | 200 | 200 | 0 | 0 |
| 1036 | reo | 1 | 1 | 0.26094 | 0 | 19900 | 19900 | 0 | 0 |
| 1036 | expression_rf | 1 | NA | NA | NA | 2000 | 2000 | NA | 0 |
| 1036 | reo_rf | 1 | NA | NA | NA | 19900 | 500 | NA | 0 |
| 1036 | expression | 2 | 1 | 0.27326 | 0 | 2000 | 2000 | 0 | 0 |
| 1036 | expression_matched | 2 | 1 | 0.20414 | 0 | 200 | 200 | 0 | 0 |
| 1036 | reo | 2 | 1 | 0.28437 | 0 | 19899 | 19899 | 0 | 0 |
| 1036 | expression_rf | 2 | NA | NA | NA | 2000 | 2000 | NA | 0 |
| 1036 | reo_rf | 2 | NA | NA | NA | 19899 | 500 | NA | 0 |

若某个EN外折的非零分子系数为0，其预测仅来自组织效应，不能把全局Spearman解释为转录组或基因对信号。ratio=1是理论alpha_max零解端点：实现显式采用解析零分子系数，防止浮点epsilon污染组织内排序；该修复只适用于理论端点，没有把较弱正则路径的真实非零系数任意截断。旧pilot与历史文件保持原样。

| 药物ID | 模型 | 外折数 | 候选边界外折数 |
| --- | --- | --- | --- |
| 1036 | expression | 3 | 2 |
| 1036 | expression_matched | 3 | 3 |
| 1036 | expression_rf | 3 | 2 |
| 1036 | mean | 3 | 0 |
| 1036 | reo | 3 | 3 |
| 1036 | reo_rf | 3 | 2 |
| 1036 | tissue | 3 | 0 |
| 1931 | expression | 3 | 2 |
| 1931 | expression_matched | 3 | 2 |
| 1931 | expression_rf | 3 | 1 |
| 1931 | mean | 3 | 0 |
| 1931 | reo | 3 | 3 |
| 1931 | reo_rf | 3 | 2 |
| 1931 | tissue | 3 | 0 |

相对路径最大端点被选中可以表示训练内CV偏好零分子项，具有明确的null含义；最弱端点被反复选中则提示搜索下界可能尚不足。候选边界和收敛证据应完整保留，不能据外验结果再拓展路径。同样候选次数并不保证不同表示达到相同求解难度或调参充分度。

## RF比较所包含的额外步骤

RF固定300棵树，候选max_features为['sqrt', 0.1]、min_samples_leaf为[1, 5, 10]。REO RF在每个训练折内按控制组织后的连续响应关联筛至最多500对，普通表达RF没有相应pair监督筛选。因此两种RF的差异包含表示、候选预算和监督筛选等pipeline差异，不能纯归因于编码。两者追加相同训练拟合的组织one-hot；它们的组织政策与EN的FWL不同，跨家族比较也是pipeline比较。

## 组织内描述与相关不可用情况

仅展示每药组织N≥20的描述性OOF记录，共140个模型×组织记录，Spearman为NA的记录数为70。组织/均值基线以及所有外折均为零分子项的EN没有组织内个体排序能力，相关记NA并保留原因。稀有组织没有该表的记录不等于被从主分析删除；这些数值也不是独立组织外验或机制证据。

| 药物ID | 模型 | 组织 | N | Spearman | RMSE | 相关不可用原因 |
| --- | --- | --- | --- | --- | --- | --- |
| 1931 | tissue | Ovary/Fallopian Tube | 29 | NA | 0.010765 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Myeloid | 31 | NA | 0.026808 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Skin | 33 | NA | 0.016049 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Breast | 44 | NA | 0.021994 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Lymphoid | 75 | NA | 0.054613 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Lung | 132 | NA | 0.0099881 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Pancreas | 28 | NA | 0.0091035 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Esophagus/Stomach | 46 | NA | 0.010021 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | CNS/Brain | 34 | NA | 0.011406 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | tissue | Bowel | 40 | NA | 0.0086331 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | expression | Ovary/Fallopian Tube | 29 | -0.12221 | 0.011145 | NA |
| 1931 | expression | Myeloid | 31 | -0.50196 | 0.027525 | NA |
| 1931 | expression | Skin | 33 | -0.10187 | 0.016684 | NA |
| 1931 | expression | Breast | 44 | -0.3809 | 0.023749 | NA |
| 1931 | expression | Lymphoid | 75 | -0.2776 | 0.055836 | NA |
| 1931 | expression | Lung | 132 | 0.043197 | 0.010536 | NA |
| 1931 | expression | Pancreas | 28 | -0.027016 | 0.0094965 | NA |
| 1931 | expression | Esophagus/Stomach | 46 | -0.26209 | 0.010203 | NA |
| 1931 | expression | CNS/Brain | 34 | 0.078928 | 0.011042 | NA |
| 1931 | expression | Bowel | 40 | 0.049978 | 0.0083325 | NA |
| 1931 | reo | Ovary/Fallopian Tube | 29 | NA | 0.011018 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Myeloid | 31 | NA | 0.026619 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Skin | 33 | NA | 0.016088 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Breast | 44 | NA | 0.0224 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Lymphoid | 75 | NA | 0.054408 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Lung | 132 | NA | 0.0099601 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Pancreas | 28 | NA | 0.008961 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Esophagus/Stomach | 46 | NA | 0.010249 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | CNS/Brain | 34 | NA | 0.011549 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | reo | Bowel | 40 | NA | 0.0085532 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1931 | expression_rf | Ovary/Fallopian Tube | 29 | 0.23596 | 0.0099805 | NA |
| 1931 | expression_rf | Myeloid | 31 | -0.23952 | 0.027333 | NA |
| 1931 | expression_rf | Skin | 33 | 0.19352 | 0.015725 | NA |
| 1931 | expression_rf | Breast | 44 | 0.17477 | 0.021644 | NA |
| 1931 | expression_rf | Lymphoid | 75 | -0.17055 | 0.055399 | NA |
| 1931 | expression_rf | Lung | 132 | -0.082576 | 0.010286 | NA |
| 1931 | expression_rf | Pancreas | 28 | 0.29009 | 0.0085543 | NA |
| 1931 | expression_rf | Esophagus/Stomach | 46 | 0.30768 | 0.0097336 | NA |
| 1931 | expression_rf | CNS/Brain | 34 | 0.062185 | 0.011561 | NA |
| 1931 | expression_rf | Bowel | 40 | -0.2349 | 0.0091992 | NA |
| 1931 | reo_rf | Ovary/Fallopian Tube | 29 | -0.36404 | 0.011601 | NA |
| 1931 | reo_rf | Myeloid | 31 | -0.43185 | 0.029335 | NA |
| 1931 | reo_rf | Skin | 33 | -0.039773 | 0.016568 | NA |
| 1931 | reo_rf | Breast | 44 | -0.076815 | 0.021721 | NA |
| 1931 | reo_rf | Lymphoid | 75 | -0.053883 | 0.055323 | NA |
| 1931 | reo_rf | Lung | 132 | -0.019921 | 0.010389 | NA |
| 1931 | reo_rf | Pancreas | 28 | 0.03503 | 0.0093766 | NA |
| 1931 | reo_rf | Esophagus/Stomach | 46 | 0.29534 | 0.0090494 | NA |
| 1931 | reo_rf | CNS/Brain | 34 | -0.33079 | 0.0119 | NA |
| 1931 | reo_rf | Bowel | 40 | -0.22533 | 0.0093664 | NA |
| 1036 | tissue | Ovary/Fallopian Tube | 27 | NA | 0.038987 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Myeloid | 31 | NA | 0.057364 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Skin | 33 | NA | 0.15743 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Breast | 44 | NA | 0.070164 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Lymphoid | 74 | NA | 0.040168 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Lung | 131 | NA | 0.028764 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Pancreas | 28 | NA | 0.030563 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Esophagus/Stomach | 46 | NA | 0.026892 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | CNS/Brain | 34 | NA | 0.037783 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | tissue | Bowel | 41 | NA | 0.044082 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | expression | Ovary/Fallopian Tube | 27 | -0.14581 | 0.039249 | NA |
| 1036 | expression | Myeloid | 31 | 0.15046 | 0.057607 | NA |
| 1036 | expression | Skin | 33 | -0.067452 | 0.15488 | NA |
| 1036 | expression | Breast | 44 | -0.21166 | 0.078198 | NA |
| 1036 | expression | Lymphoid | 74 | 0.069823 | 0.039921 | NA |
| 1036 | expression | Lung | 131 | 0.012456 | 0.030558 | NA |
| 1036 | expression | Pancreas | 28 | -0.035516 | 0.029982 | NA |
| 1036 | expression | Esophagus/Stomach | 46 | -0.21319 | 0.029212 | NA |
| 1036 | expression | CNS/Brain | 34 | -0.01649 | 0.036309 | NA |
| 1036 | expression | Bowel | 41 | 0.15227 | 0.043592 | NA |
| 1036 | reo | Ovary/Fallopian Tube | 27 | NA | 0.03919 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Myeloid | 31 | NA | 0.057381 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Skin | 33 | NA | 0.15729 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Breast | 44 | NA | 0.070224 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Lymphoid | 74 | NA | 0.04017 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Lung | 131 | NA | 0.028768 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Pancreas | 28 | NA | 0.030496 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Esophagus/Stomach | 46 | NA | 0.026954 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | CNS/Brain | 34 | NA | 0.037839 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | reo | Bowel | 41 | NA | 0.044078 | baseline_or_null_molecular_model_has_no_within_tissue_individual_ranking |
| 1036 | expression_rf | Ovary/Fallopian Tube | 27 | 0.35592 | 0.035373 | NA |
| 1036 | expression_rf | Myeloid | 31 | -0.26855 | 0.060424 | NA |
| 1036 | expression_rf | Skin | 33 | 0.23295 | 0.15323 | NA |
| 1036 | expression_rf | Breast | 44 | 0.10162 | 0.0659 | NA |
| 1036 | expression_rf | Lymphoid | 74 | 0.14447 | 0.039752 | NA |
| 1036 | expression_rf | Lung | 131 | 0.16899 | 0.029581 | NA |
| 1036 | expression_rf | Pancreas | 28 | 0.34154 | 0.028899 | NA |
| 1036 | expression_rf | Esophagus/Stomach | 46 | 0.1292 | 0.027327 | NA |
| 1036 | expression_rf | CNS/Brain | 34 | -0.11719 | 0.041683 | NA |
| 1036 | expression_rf | Bowel | 41 | 0.12178 | 0.043558 | NA |
| 1036 | reo_rf | Ovary/Fallopian Tube | 27 | 0.14042 | 0.037882 | NA |
| 1036 | reo_rf | Myeloid | 31 | -0.2754 | 0.058749 | NA |
| 1036 | reo_rf | Skin | 33 | 0.12667 | 0.15454 | NA |
| 1036 | reo_rf | Breast | 44 | 0.043693 | 0.070451 | NA |
| 1036 | reo_rf | Lymphoid | 74 | 0.0025028 | 0.040284 | NA |
| 1036 | reo_rf | Lung | 131 | 0.30489 | 0.028145 | NA |
| 1036 | reo_rf | Pancreas | 28 | 0.12917 | 0.032636 | NA |
| 1036 | reo_rf | Esophagus/Stomach | 46 | 0.032377 | 0.028213 | NA |
| 1036 | reo_rf | CNS/Brain | 34 | -0.21161 | 0.045333 | NA |
| 1036 | reo_rf | Bowel | 41 | 0.19077 | 0.042442 | NA |

## 运行资源、实现修复及可复现记录

| 总墙钟秒 | CPU秒 | 采样峰值RSS GiB | Windows峰值working set GiB | RSS采样间隔秒 | 实际fits | BLAS线程 | 候选并行数 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1125.8 | 1081.5 | 0.84004 | 0.87965 | 0.1 | 906 | 1 | 1 |

内存为整个单进程口径，包含输入准备、库、映射页面和拟合；采样峰值RSS及Windows历史峰值working set分列，均不能称为仅模型矩阵内存。流程时间包括内层搜索、外层重拟合和预测。不同模型重复做训练内特征构造是该资源测量的一部分。记录线程限制与并行预算；若线程池信息在限制上下文退出后采集，它反映恢复状态。

- 理论null端点采用显式解析零系数，已保留solver_reason、actual_alpha、非零数及dual gap证据；不修改旧pilot数字。
- 新维护边界：若合成target仅含组织均值，浮点舍入可引出伪非零EN系数或RF pair筛选。运行结束后已用逐组织精确常量检测修复，38项回归/防泄漏验收通过。实际906fits没有触发该边界；最小全训练集组织内响应残差RMS=0.020150，远高于舍入残差。独立审计只读已记录OOF/成员，维护修复不改写执行5d05349或数值。
- 命令行相对配置路径修复只影响路径解析和可移植执行，不改变科学参数；首个探针在任何拟合前失败并保留失败/资源证据，路径修复后完成8次探针。运行config和源码SHA保存在来源JSON中。
- 原异常句柄关闭、失败状态和资源收尾、逐外层检查点修复应继续用于正式长运行；成功pilot尚不能代替完整续跑验收。
- 报告只发布汇总指标、outer审计和配置，不公开个体响应或fit成员ID；run目录中的详细内部审计由既定数据管理政策保留。

## 下一gate与尚未完成的证据

- 根据本轮求解、路径边界与资源证据冻结正式模型/统计配置；若调整配置，写明工程或统计理由并保留本运行，不能追求REO胜出。
- 按固定8药执行外层5折×3次、内层3折，完整报告主EN、辅助RF、组织/均值基线及匹配表达消融。重复CV改善划分稳定性证据，不增加独立样本N。
- 实施预设生物组重采样、配对效应及全流程置换；固定检验族与多重性规则，不把非显著差异称为等效。
- 正式开发结束后在全开发拟合并冻结模型/特征/代码，才开放CTRP标签进行一次预设排序迁移评价。不同药敏剂量和共享表达来源的限制仍成立。
- 新旧pilot同时改变资格、全局fold、基因预算及alpha路径，不能把分数变化归因于单一因素。REO可能有增量、接近或更差，三种结果均须按预设任务完整报告。

可机读产物：results/pilot_reo200_rf/summary_metrics.csv、fold_metrics.csv、paired_metric_differences.csv、paired_fold_metric_differences.csv、outer_fit_audits.csv、chosen_parameters.json、feature_audit_summary.json、within_tissue_metrics.csv及resource_and_provenance_summary.json。
