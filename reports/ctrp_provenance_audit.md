# CTRPv2.0 原始来源与身份完整性审计

日期：2026-10-10。执行：`python -X utf8 rebuild_v2/src/curate_ctrp.py --download`。

本次已取得原始发布包中的完整药敏曲线表及其细胞、化合物、实验、培养板元数据，建立了可追溯到原文行号的索引。全过程没有展示 AUC 原值、计算响应分布、评估预测表现，亦未按响应大小选择药物或细胞。原项目只读。

## 1. 下载路径、版本与完整性

历史 NCI 路径 `https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/` 当前重定向/连接异常。DepMap 下载列表接口返回浏览器验证 HTML，未作为数据读取；Synapse 公共文件元数据可读，但匿名实体文件请求被拒绝，未绕过权限。

可用归档地址由 [BHKLAB-DataProcessing 的 CTRPv2 管道配置](https://github.com/BHKLAB-DataProcessing/CTRPv2-Pharmacoset_Snakemake/blob/main/config/pipeline.yaml) 明确列出：

- [CTRPv2.0_2015_ctd2_ExpandedDataset.zip](https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset.zip)
- 归档 HTTP Content-Length：342,737,645 字节。没有下载完整 ZIP，而是读取中央目录并通过 HTTP Range 提取必要成员。
- 本次固定 **v2.0 原发布包**；没有将 v2.1、其他重拟合版本或 PharmacoGx 衍生数值混用。CDN 的文件更新时间不等于原研究数据采集时间。
- 9 个提取成员均通过 ZIP CRC32 检查；除 `MANIFEST.txt` 自身外，8 个成员另外通过归档内发布的 MD5 校验；每个成员记录本地 SHA256。传输校验不等于独立证明数据科学有效性。
- 来源层级准确表述为：**由可追溯的学术整理者仓库指出的历史 CTD² 归档 CDN 内容**。不将 CloudFront 主机擅自描述为当前 NCI 官网目录。

记录保存在 `data_raw/ctrp/download_manifest.json`、`archive_members.json`、`archive_provenance.json`。发现路径所用配置和下载规则已保存为 `bhklab_pipeline.yaml`、`bhklab_download.smk`，下载失败/验证页也与真正数据分开记录。

|成员|解压后字节数|用途|
|---|---:|---|
|v20.data.curves_post_qc.txt|59,801,043|全部已发布 post-QC 曲线、原 AUC 和拟合字段|
|v20.meta.per_cell_line.txt|71,591|master_ccl_id、原始细胞名称、组织/组织学|
|v20.meta.per_experiment.txt|95,862|experiment→cell、日期、run、SNP 身份状态|
|v20.meta.per_compound.txt|116,632|master_cpd_id、化合物名、SMILES、剂量设计；药物身份由独立映射流程审计|
|v20.meta.per_assay_plate.txt|422,350|实验→培养板和 DMSO 统计字段|
|v20._COLUMNS.txt|9,480|发布方字段说明|
|v20._README.txt|1,706|发布方数据文件说明|
|CTRPv2.0._README.docx|13,313|同版 Word 说明|
|MANIFEST.txt|886|发布方成员 MD5|

曲线文件 SHA256：`b97bb9004f97ff7d5910c0f66c6d26163a2ce88529eaefc5df0e51624e5ac6a9`。

## 2. 覆盖与主键检查：仅身份/可用性

|检查|结果|
|---|---:|
|曲线记录|395,263|
|有曲线的 master_cpd_id|545|
|有曲线的 master_ccl_id|887|
|有曲线的 experiment_id|907|
|唯一 experiment_id × master_cpd_id|395,263；无重复|
|唯一 compound × cell|387,130|
|含重复实验的 compound × cell|7,708|
|细胞元数据记录|1,107|
|实验元数据记录|1,061|
|实验元数据唯一 experiment_id|907|
|有多条元数据的 experiment_id|153|
|同 experiment_id 映射多个细胞|0|
|曲线找不到细胞或化合物元数据|0|
|AUC 字段缺失/非有限数值|0；仅检查可读取性，未统计数值大小|

**实验元数据不能直接 join。** 153 个 experiment_id 对应多条元数据，日期不同，其中 1 个的 run_id 也有不同值。但各自 master_ccl_id、培养条件及 SNP 状态一致。脚本先验证同 ID 的细胞身份唯一，再输出一行身份记录，并以 JSON 保留全部源元数据及行号。直接 many-to-many 合并会扩增曲线记录。

原文件没有另一个名为 `screen_id` 的独立键。输出的 `CTRPv2.0:experiment:<id>` 是明确声明的命名空间包装，不是新发现的实验标识，也不代表一个新的生物学个体。

## 3. SNP 身份 QC 与曲线 QC

发布方 COLUMNS 将 `snp_fp_status` 定义为确认细胞系身份的 SNP 基因分型状态，更新日期为 2015-12-15。按唯一 experiment_id 计数：

|原始身份状态|实验数|
|---|---:|
|SNP-matched-reference|770|
|SNP-unconfirmed|36|
|SNP-not-tested|101|

同一 experiment_id 没有状态冲突；当前文件中没有明示 mismatch/fail 的标签。`unconfirmed` 不能解释为已经确认身份错误；`not-tested` 也不能解释为已匹配。README 没有给出相反定义。

重构主流程已提出：主外验只纳入 `SNP-matched-reference`，其余两类保留为预先冻结的敏感性分析。这属于新增的保守身份资格规则，应报告过滤后的细胞、药物和组织覆盖变化。curation 文件本身保留全部记录和状态，由冻结的下游资格规则执行筛选，防止不可追溯地删行。

`conc_pts_fit` 是参与曲线拟合的浓度-反应点数，`fit_num_param` 是拟合自由参数数。原包化合物列说明给出名义 16 点、2 倍连续稀释设计；实际 post-QC 曲线的拟合点数不是恒等于 16。

- 全部发布曲线的拟合点数至少为 8；319,236 条为 16 点；835 条超过 16 点，最大为 29。完整点数计数见 `ctrp_curation_summary.json`。
- 8, 9, 10, 11, 12, 13, 14, 15 点的记录分别为 1,431、907、1,398、2,255、5,485、7,652、16,910、39,154 条。
- README/COLUMNS 没有明确声明可直接移植的“本版统一最低拟合点数”准入规则。观察到最少 8 点，不等于证明作者采用了一个通用的 ≥8 规则。
- 不新增“必须恰好 16 点”筛选，也不根据曲线反应大小筛选。暂沿用已发布 post-QC 曲线及其拟合信息。
- 本次未下载庞大的 raw-well、pre-QC、post-QC 逐剂量表，因此没有重新运行原始实验 QC，也未复原逐点 censoring 的 `qc_type`。不得在 Methods 宣称已重做这些步骤。若后续科学问题确需统一剂量区间重积分，应作为明确的新 endpoint 分析另立版本，而不是静默替换当前终点。

## 4. AUC 定义与响应隔离

`area_under_curve` 的本版字段说明是 sigmoid 浓度-反应拟合曲线的积分面积，属于药物反应曲线面积，**不是分类 ROC-AUC**。所有原值及 curve-fit 置信区间/参数都在原始曲线表中保留，索引以物理行号可精确定位。

不要仅因字段名为 AUC 就假定它采用 GDSC 的数值尺度或必定标准化至 [0,1]。CTRP 作者的 [2018 更新计算实现](https://github.com/remontoire-pac/ctrp-reference/blob/master/auc/README.md) 明确另含统一积分区间和 [0,1] 归一化步骤，这不能反向证明 2015 原曲线表已经执行同一终点变换。当前选择保留原版，不混入重新拟合值。跨库主要外验使用药物内秩关联；不把不一致的积分定义直接用于跨库绝对 RMSE。

到本次审计结束，没有输出任何原始药敏 AUC 值、分位数、范围、相关性或模型成绩。机器处理仅判断字段是否为有限数值；没有将其作为筛选条件。下一步应先冻结 compound 映射、独立细胞/生物组划分、资格规则、模型和统计方案，再做外部预测及评估。

重复实验保留原记录；解盲阶段按预先固定的规则在 **同 compound × 同 canonical cell** 内取 median AUC。n 按独立细胞/生物组计，不能把重复 experiment 当成独立外验样本。median 规则尚未在本次 curation 中执行，也没有查看聚合值。

## 5. 旧本地表的可用范围

旧 `3.Project/Pan-Drugs-REOs/data/CTRP--Drug--Cellline--AUC.txt` 只有 `Drug / Cellline / CCLE_primary_site / AUC` 四列。

对原 archive 以 ID 关联官方名称后，其 **(药名、细胞名、组织标签) 含重复次数的完整集合** 与旧表逐项一致：395,263 行，没有任何多余/缺少的身份行。这说明旧表身份覆盖与本次 archive 相符；不能仅据此证明历史脚本来源、响应数值未变、行顺序相同或已正确处理重复实验。本次没有比较 AUC 值是否相等。

新流程应使用带 experiment/master ID 的原始 archive 索引。不要继续把旧四列表当作具有独立外验资格的终点表，也不要静默按药名/细胞名去重。化合物同名、混合药、立体结构问题由 `drug_*` 身份审计独立解决。

## 6. 下游接口

|输出|主键/用途|
|---|---|
|`data_processed/ctrp_curve_index.csv`|一行一条曲线；source_row_number 从原 TSV 第 2 行起；experiment_id、派生 screen_id、master_cpd_id、master_ccl_id、官方 ccl_name、组织、response_available、SNP 状态、fit 点数/参数数；不含响应值|
|`data_processed/ctrp_experiment_metadata.csv`|唯一 experiment_id，归并身份但不丢弃原元数据 JSON|
|`data_processed/ctrp_cell_metadata.csv`|唯一 master_ccl_id；全部 1,107 条，另标记是否有药敏曲线|
|`data_processed/ctrp_pair_availability.csv`|compound × cell 的实验 ID 列表和可用记录数；不聚合响应|
|`data_processed/ctrp_drug_availability.csv`|每 compound ID 的身份覆盖数；不提供反应统计|
|`data_processed/ctrp_curation_summary.json`|身份完整性、QC 计数、原文件 hash 和处理范围|
|`data_processed/ctrp_legacy_identity_summary.json`|旧表与 archive 身份集合比较，不比较响应|
|`data_raw/ctrp/v20.data.curves_post_qc.txt`|唯一的本版原始响应存储；由记录 ID 与源行号连接|

官方细胞元数据只有 `ccl_name`，没有 `CCLE_Name` 或 DepMap `ModelID`。输出不会伪造一个“官方 ccle_name”。统一 ACH、Sanger ID、同来源生物组映射由另一个独立身份流程完成，应记录映射来源与歧义，而不凭部分字符串猜测。

## 7. 验证记录与限度

完成：下载成员 CRC32 + 随包 MD5 + 本地 SHA256；所有源 ID 引用有效；experiment→cell 无歧义；experiment×compound 曲线唯一；没有因多行实验元数据导致扩增；原始 archive 不覆盖；程序可独立再次运行；数据索引与报告不暴露外验反应大小。脚本只使用 Python 标准库。

尚未完成且没有声称完成：canonical ACH/生物组最终资格、跨库药物结构最终裁决、样本拆分冻结、跨库外验、原始逐剂量拟合复算，以及最终投稿所需的完整独立复现。这些属于后续已分工任务。
