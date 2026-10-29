# Stage 7 零解计数展示勘误

日期：2026-10-29。范围仅限既有聚合计数的汇总与展示；未重训、评分、运行原诊断分析器，也未读取开发集或CTRP逐样本响应、OOF预测或模型检查点。

## 原因与影响

原`src/analyze_development_diagnostics.py`生成药物异质性汇总时，仅取`repeat == 0`的零分子拟合计数，却在报告第二张表中以全部三个repeat的15次外层拟合作为分母。错误影响`results/development_diagnostics/drug_heterogeneity.csv`中的七个`*_zero_fits`列及报告第二张表最后一列；第一张零解诊断表、Table S6、840次分子外层拟合及430个零解总计正确。性能、区间、模型系数与预测均不受影响。

## 修正与验证

修正从`results/ordinal_repeated/molecular_and_failure_summary.csv`按药物和模型合计repeat 0、1、2；每个模型的分母仍为15次外层拟合。56个药物×分子臂计数全部与`tables/manuscript_v1/TableS6_zero_solution_audit.csv`一致。其中40个CSV单元格发生数值改变，其余16个本来即与全repeat总计相同。CSV所有非零解计数字段保持原字符串值。

| 药物 | 原第二表显示 | 正确关联pair零解/15 |
| --- | ---: | ---: |
| MIRA-1 | 3/15 | 7/15 |
| PLX-4720 | 5/15 | 13/15 |
| Bortezomib | 4/15 | 11/15 |
| Telomerase Inhibitor IX | 0/15 | 0/15 |
| Tamoxifen | 5/15 | 14/15 |
| BMS-536924 | 1/15 | 1/15 |
| Gefitinib | 0/15 | 1/15 |
| Cediranib | 2/15 | 7/15 |

修正后的分析器使用独立聚合函数，拒绝缺失repeat、重复repeat及超出外层拟合数的零解计数。四个纯聚合回归测试通过，分别覆盖原repeat-0遗漏、分母不完整、重复计数和计数越界。测试没有读科研数据。首次使用内置Python因缺少`psutil`未执行测试，随后使用项目已有`.venv`成功；未安装或更改环境。

## 来源与历史保留

- 原分析器、CSV、报告和manifest完整保存在`local_private/stage7_before/`相应相对路径下。
- `results/development_diagnostics/diagnostics_manifest.json`中的原执行提交、分析器SHA、运行日期、耗时、内存及输入哈希仍描述2026-10-14原执行；不能将它们理解为本次修正的执行信息。
- manifest的`original_execution_output_sha256`保留原输出绑定；`output_sha256`绑定当前文件，只有被更正的CSV哈希变化；`artifact_correction`指向本次独立勘误凭据。
- [机器可读勘误凭据](stage7_correction_audit.json)记录输入、修正前后SHA、40项数值修改及验证结果；[一次性修正程序](correct_stage7_aggregate_counts.py)仅读取现存聚合表及备份，拒绝覆盖已存在的凭据。

本次勘误没有改变430/840、任何逐药性能、置信区间、调参选择或最终科学结论。后续引用药物异质性中的零解数，应使用修正表或Table S6的三个repeat合计，不能引用原第二张表的repeat-0计数。
