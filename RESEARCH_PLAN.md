# 基因对药物反应预测项目重构科研计划

**计划v2.0，更新日期：2026-10-10。唯一维护版本为本文件。工作周期：10月10日至11月8日，11月9日至10日保留投稿检查缓冲。**

本计划以最新投稿稿件、主要分析代码、保存的结果矩阵/模型台账、补充表、图2/图4和公开数据源核查为基础。目标是一个月内完成可审查、可复现、论断与证据相符的新稿并投稿；接收和正式上线还取决于编辑及外审，不能由一个月计划保证。用户已确认Naunyn-Schmiedeberg's Archives of Pharmacology亦已拒稿；按1名主要执行者、与助手协作、现有电脑和公共数据设计，版面费不是主要限制。

## 0 执行台账与当前版本（每次工作先读此节）

**首选期刊：BMC Genomics。当前：数据、身份、覆盖和队列构建完成；队列v1正式验证通过并已锁定。模型拟合为0，尚无新的预测性能。**

本文件是唯一持续维护的总计划。原计划是历史快照，阶段报告是证据；发现错误、版本/范围/方法变化，应先更新本节及变更表，再建立对应配置版本。不得因外验成绩不好更换主药、身份规则或测试样本。

|工作包|状态|依据/下一步|
|---|---|---|
|P0 数据/代码审计|完成|[阶段1报告](reports/阶段1_数据核验与队列锁定.md)、[代码处置](reports/现有代码处置清单.md)|
|P1 下载/表达|必要部分完成|本地TPM 1,019×57,820、float64；补CTRP原包成员、结构/身份/GENCODE；[下载登记](data_manifest/download_register.csv)|
|P2 细胞/药物/准入|完成|64个GDSC IDs对应63种结构支持CTRP单药；全部295 IDs有样本量台账；冲突排除|
|P3 队列v1/验证|完成并锁定|627开发、180外验唯一模型/生物组；药物×模型4,975/1,222不能累加当独立N|
|P4 开发专用试跑|进行中|配置[pilot_v0.1.json](config/pilot_v0.1.json)已在首次真实拟合前写定；MIRA-1/PLX-4720，开发专用嵌套流程、资源和失败审计|
|P5 模型/统计冻结和正式运行|未开始|外层5折×3重复/内层3折；所有监督筛选仅在内层训练；生物学单位CI/置换|
|P6 冻结外验|未开始|最终模型/统计配置冻结后评估CTRP响应，完整报告8药|
|P7 图表/正文/投稿|未开始|结果决定结论；不承诺旧高准确率或接收|

### 本轮输入与面板

表达为2018 CCLE RSEM TPM；2024 Broad Model和2026 Sanger快照只补身份/同源关系，未混用表达版本。开发来自GDSC2 27Oct23；外验为原始CTRPv2 expanded post-QC、SNP核验通过、在全部GDSC2及已知关联组之外的模型。自然外验池183，8药覆盖180；未启用20%共有模型保留。

|化合物|GDSC DRUG_ID|CTRP master_cpd_id|开发独立组|主外验独立组|
|---|---:|---:|---:|---:|
|MIRA-1|1931|580922|626|113|
|PLX-4720|1036|606256|623|172|
|Bortezomib|1191|411737|621|173|
|Telomerase Inhibitor IX|1930|411865|621|173|
|Tamoxifen|1199|26972|621|170|
|BMS-536924|1091|411859|621|169|
|Gefitinib|1010|52926|621|165|
|Cediranib|1922|606136|621|87|

逐药浓度/组织覆盖见[main_drugs.csv](results/main_drugs.csv)。面板含研究探针，不是8种均获批药物；药理类别只描述，未用于改面板。

### 解释边界与未完成项

- 旧全标签DEG/Fisher/RFE及下游结果须重算。CTRP历史用过，本次为回顾性重新设计的身份隔离验证，不能称前瞻盲测。
- 两库均用CCLE表达，当前检验跨药敏screen及未知细胞系排序迁移，不证明跨表达平台或患者效用。
- 每药GDSC2固定浓度区间，按独立组/模型覆盖及浓度数值排序；仅解决库内混合。跨库积分窗口不同，Spearman也不能消除其改变排序的可能。暂不做跨库绝对RMSE/校准主张。
- 共同浓度窗重积分先查拟合参数可行性；若增加，先设端点版本再计算，不能事后挑更好端点。
- 独立范围只指全部本地GDSC2及已知关联组；未知亲缘仍有限制，不能扩称所有GDSC1/2。
- 原项目兼容输入仍需存在；公开数据独立重放和干净环境全流程安装是发布前待验项。[复现状态](docs/REPRODUCIBILITY.md)写清完成/未完成。

### 计划变更记录

### 阶段2试跑协议（首次真实拟合前确定）

固定已冻结面板前两种药物MIRA-1/PLX-4720，不根据响应或性能挑药。使用各自完整开发样本626/623、同一627模型全集的全局分组；外层3折×1次、内层3折。正式分析的5折×3次尚未执行。模型为训练均值、组织Ridge、2,000高变基因表达EN、同80基因表达EN、80基因REO EN；正式200基因只在首外层训练做不使用标签的构造资源测时，RF未进入此次试跑。

低检测规则TPM≥1且至少10%训练样本检测、log2TPM MAD>0、稳定Ensembl ID并列排序；REO恒定pair及tie>90%剔除。每次内层fit重建全部筛选、编码、标准化、组织设计矩阵及目标均值/SD。组织项在EN中不惩罚：先从TPM构造REO，再在训练内对特征和目标做组织投影，预测时恢复组织项；未见组织的效应为0，组织基线回退训练全局均值。现有TPM无缺失；出现缺失即报告失败，尚未实现正式外部missing/abstain规则。

每个非均值模型固定6个候选，以原AUC尺度的合并内层OOF RMSE选择，并列选配置列表首项；记录网格边界和收敛。α=0.1/0.03/0.01、L1比例=0.5/1，目标SD及特征缩放仅在当前训练拟合。试跑分数只作探索描述，不进行P值/CI、不按分数换药或宣称REO胜出。均值基线每外折恒定，其合并OOF排序可能产生折均值假象，总体Spearman明确NA，误差照常报告。

只读取开发清单所列GDSC源行响应并核对药物、曲线ID、精确浓度窗和0–1范围；CTRP原包、旧响应及响应索引在运行时禁止打开。官方AUC方向与定义见[Sanger药敏文档](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/)、[GDSC拟合说明](https://github.com/CancerRxGene/gdscIC50/blob/master/vignettes/gdscIC50.Rmd)。发布AUC来自上游联合拟合，本轮严格训练内拟合的范围是下游预测流程，不声称上游端点估计也已按CV重拟合。

预算462次真实拟合、6,245条预期OOF预测（626+623，每个模型1次OOF）；单进程、BLAS单线程，100ms采样完整进程RSS。首个完整外层模型完成后据实测时间估算总耗时。配置、源代码和输入哈希、实际fit成员、每外层特征、OOF预测、收敛/失败、耗时和资源都写入本地runs；阶段报告和聚合结果再同步仓库。试跑前独立合成测试已验证持留标签/表达变化不能改变超参数或训练产物，原队列v1保持锁定。

|日期|变更|依据|解释影响|
|---|---|---|---|
|2026-10-10|首选BMC Genomics，排除全部12本已拒期刊|用户确认及范围核查|不再将已拒BMC Bioinformatics/Human Genomics列备选|
|2026-10-10|固定现有2018 TPM，无需重下载全表达|找到完整矩阵、数值/秩序核验|未宣称远端压缩包字节同一性|
|2026-10-10|纳入patient/parent图及SNP matched准入|权威元数据可得|按独立生物组计数，未知亲缘仍保留|
|2026-10-10|优先外验≥80，固定8药|58种达80；未看响应大小或效能|计算前覆盖/精度修订，不按成绩挑药|
|2026-10-10|自然CTRP-only外验，不预留GDSC2的20%|183自然组覆盖足够|开发可用N较多，全局角色隔离|
|2026-10-10|加入库内浓度窗筛选及跨库限制|76个GDSC IDs多窗，跨库上限不同|外验保留剂量设计差异，不以Spearman掩盖|
|2026-10-10|队列v1正式锁定|1,607,707条断言通过，输入/实现/样本哈希完整|后续变化必须先记录并建立新队列版本|
|2026-10-10|Git仓库与唯一总计划|用户指定RascalZ/REO|原稿/大数据留本地；代码/聚合队列版本化|
|2026-10-10|干净clone的最小公开检查通过|隔离Python标准库运行公开哈希/计数检查，全部源码语法检查通过|仅验公开入口；完整源数据重放及模型结果复现仍待完成|
|2026-10-10|阶段2工程pilot采用两药/80gene/3×3 CV|用户授权小规模试跑；先验证状态边界和算力|非正式200gene/5×3/RF研究结果；6candidate inner RMSE预先固定|
|2026-10-10|明确组织不惩罚与上游端点范围|FWL实现及官方联合拟合说明|无泄漏声明限下游嵌套管道；未知组织回退显式验收|

---

## 1 核心判断与发表可行性

**项目可以重构，但需要重新完成主要建模与验证，不能按“润色后换刊”处理。** 可以保留的是基因对表示思想、数据来源线索、部分ID映射与绘图函数；不能直接保留的是目前的高准确率/AUROC、已训练模型的泛化结论、核心基因对的确定性和临床转化宣称。

推荐把论文中心从“构建泛癌高准确率耐药图谱”调整为：

> 在未见过的癌细胞系上，样本内基因间相对表达排序能否提供超出癌种与普通表达量模型的药物反应预测信息？它在独立药敏筛选来源中何时可迁移、何时失效？

这是**待检验研究问题** ，不是预定REO一定优于基线的结论。第一周冻结模型与统计方案后，无论结果支持或不支持REO，都完整报告预设任务；不能再从大量任务中留下高分条件来证明普遍有效。

| 层面 | 当前判断 | 重构后的机会与边界 |
|---|---|---|
| 科学问题 | 基因对可作为表示方法，但“REO+随机森林”本身创新有限 | 有价值的贡献应是公平比较、可靠的跨来源结果和清楚的适用条件 |
| 现有性能可信度 | 主要流程存在监督筛选先于划分，旧泛化效能证据不足 | 全流程嵌套评估、外部冻结测试后可以修复；结果可能显著下降 |
| 现有样本量 | 逐药物–癌种模型过小，不适合原来的大规模强结论 | 主分析改为按药物汇集癌种，并控制组织来源；局部癌种结果为次要分析 |
| 外部证据 | CTRP有细胞系重叠；GEO是少数细胞系重复/处理状态，非患者验证 | 新细胞系、独立药敏来源的测试可建立；临床效益仍不能据此证明 |
| 一个月可行性 | 同时重建上千分类器、扩展多组学/临床/单细胞不现实 | 本轮固定8种主要化合物，2种表示×2类轻量模型、1个冻结外部来源，较可执行 |
| 发表前景 | 原稿现状不建议立即再投 | 满足本计划证据门槛后，具有进入范围匹配期刊评审的合理基础；无法给出可信的录用百分比 |

2026年已有专门讨论药敏模型评估的DrEval，也已有REO用于药物反应预测的工作。新稿应承认这些先例，避免“首次全面”“天然平台无关”等无范围限定的优先性宣称。DrEval可作为评价实现参考，不应直接照抄其默认数据版本、端点或全部配置。[Bernett等，2026](https://doi.org/10.1038/s41467-026-72903-w)

## 2 本地审计揭示的问题

以下数字来自**当前保存文件和脚本** ，不等于已经复现了历史投稿运行。历史证据在原计划归档；当前处置见[现有代码处置清单](reports/现有代码处置清单.md)。

| 优先级 | 已发现的问题 | 为什么影响发表 | 必须采取的动作 |
|---|---|---|---|
| P0 | Step2全样本DEG、Step3全样本Fisher、Step4全样本RFE，Step5才划分训练/测试 | 测试标签已参与特征形成，不能估计未知样本性能 | 从原始表达与药敏表重新划分；每个内层训练折重新执行所有数据驱动步骤 |
| P0 | 1,923个匹配后条件N中位16，IQR12–31，范围5–51；1,410个至少一类≤10 | “两组各>10”在表达匹配后并未维持；重复CV不能创造新样本 | 基于匹配完成后的独立细胞系数选任务；改变主分析单位 |
| P0 | 现存流程计数为1,923→353 Fisher输出→332稳定对→255最终特征矩阵→231模型 | 稿件混淆候选条件、可建模条件及成功模型，易产生幸存者偏差 | 统一纳入台账和排除代码；报告失败条件与预测覆盖率 |
| P0 | Step5计算概率但ROC传入硬分类标签；227个非缺失AUC均等于balanced accuracy | 不是通常意义的概率排序AUROC，且自动选择方向可能遮蔽反向结果 | 使用事先指定正类的概率/连续评分，固定方向；保留低于0.5的结果 |
| P0 | GEO专用脚本前段drugname为AZD7762，后段重训才为Lapatinib，输出名均写Lapatinib | 当前版本有错药模型风险，无法确认图4历史来源 | 由drug ID、模型manifest、pair schema三重检查；重新生成外验，不能只改字符串 |
| P0 | CTRP元数据26条件均有与原特征筛选样本重叠；574/922条条件×细胞系观测重叠 | 独立数据库不等于独立生物个体；该比例不是去重个体比例，也不是精确70%训练集重叠率 | 用canonical cell-line ID重新确定开发/测试集合；所有特征选择排除测试个体 |
| P1 | 稿件500/800/1000树、sqrt(p)/2sqrt(p)、5折/5×3/Bootstrap等混用 | 读者无法判断实际算法和重现效能 | 配置文件为唯一参数来源，Methods从最终执行配置生成 |
| P1 | 最终255个矩阵特征数中位10，稿件写312；最新S1写DepMap24Q1而代码来源指向2018资料 | 数字和版本来源无法追溯 | 新版锁定输入与输出哈希、版本和独立样本数 |
| P1 | 同名药物多ID：Docetaxel24个条件实际是1007的15个+1819的9个，只有15个独立癌种 | 把实验ID重复当作化合物/癌种覆盖会夸大范围 | 区分compound ID、assay ID、batch，不按名称盲目合并 |
| P1 | 单细胞0值/相等表达均被归为同一顺序，未证明独立重复；正文混用pseudo-bulk与细胞级 | 可能检测深度和处理状态，而非真实耐药；细胞数不能充当样本量 | 主稿删除验证地位，按生物重复聚合并报告可测率，细胞图只作探索 |
| P1 | S3B与S3C标题不同，之后的数据完全相同；S2文件dimension及drawing关系损坏 | 补充材料不能支持对应分析或可靠机器读取 | 全部表格由结果脚本生成，发布XLSX+CSV并回读核对 |
| P1 | 已发现GSE156246错引来源及若干基因/机制与引用对象不符 | DOI真实也可能不支持正文句子 | 对核心论断逐条做claim–citation核对；见稿件审计例证 |

**不要把尚未证实的问题写成事实。** 当前代码中确有count输入分支，抽查为整数，DESeq2的设计式为`~ condition`；所以不能直接断言“TPM输入DESeq2”。应核对全量来源、转换、重复基因处理和实际运行路径后再下结论。同样，某个保存脚本写错药不自动证明历史图件一定用了错模型，但足以要求重跑和追溯。

## 3 从目标期刊反推新稿

### 3.1 排除既往拒稿期刊后的主推荐

用户列出的12本已拒期刊全部排除本轮候选，包括已确认拒稿的Naunyn-Schmiedeberg's Archives of Pharmacology。BMC Bioinformatics、Human Genomics已拒；下列BMC Genomics、BMC Medical Genomics是不同期刊，并未出现在拒稿清单中。同属出版社或使用同一平台，不等于同一本期刊；本轮也不以改题后再投原刊作为默认路径。

推荐顺序为**BMC Genomics → Scientific Reports → BMC Medical Genomics** 。三刊都在Springer Nature的Snapp平台生态，两本BMC格式最接近。排序考虑新研究问题与证据类型的契合，不是按影响因子或未经核实的“容易中”排序。

| 顺位 | 期刊 | 本稿最合适的定位 | 核实的时间指标 | 当前官网标准APC |
|---|---|---|---|---|
| 1 | BMC Genomics | 转录组表示的系统比较、药物反应预测、可重复性与适用边界 | 首次决定中位4天；本次未取得可直接比较的投稿至接收/接收至上线中位数 | USD3,190，税及机构协议另计 |
| 2 | Scientific Reports | 严谨的跨来源实证比较；技术有效、外部证据完整 | 官方2025指标：首编辑决定12天；投稿至接收138天 | USD2,850，税及机构协议另计 |
| 3 | BMC Medical Genomics | 药物基因组学问题和药物特异结果更突出的版本 | 首次决定中位3天；接收/上线时间未核实 | USD3,290，税及机构协议另计 |

来源：[BMC Genomics](https://link.springer.com/journal/12864)、[BMC Medical Genomics](https://link.springer.com/journal/12920)、[Nature Portfolio官方2025指标](https://www.nature.com/nature-portfolio/about-journals/journal-metrics)、[Genomics费用](https://link.springer.com/journal/12864/how-to-publish-with-us)、[Scientific Reports费用](https://www.nature.com/srep/open-access)、[Medical Genomics费用](https://link.springer.com/journal/12920/how-to-publish-with-us)。价格以接收日及具体安排为准。

**3/4/12天是首决口径，可能包括快速拒稿，不能当作外审完成或录用速度。** Scientific Reports的138天是已接收稿件的历史中位数，也不是本稿承诺。此次没有找到三个期刊可比的当前接收至上线数据，因此不编造“最快见刊”排序。修复证据和减少因不匹配而被拒的轮次，比追求几天的首决数字更影响总见刊时间。

相对机会判断：新稿若完成本轮8种预设化合物的预设公平比较、冻结外部评估及开源复现，BMC Genomics最契合；若主要为清楚的实证结果且算法创新弱，Scientific Reports仍可考虑；若得到明确的药物基因组学新发现，可将BMC Medical Genomics上调。只有小样本内部CV或仍有来源不清的外验时，三者都不应视为高成功率目标。

本轮不优先BMC Cancer：其现行范围对仅公共数据库计算/预测以及标志物验证有明确要求，不能把同一小GEO内重训包装为满足要求。既往拒过的BMC Bioinformatics、Human Genomics等也不作为默认重投。完整官方依据、实际投稿入口和格式差异见`journal_evidence.md`。[BMC Cancer范围](https://link.springer.com/journal/12885/aims-and-scope)、[Snapp官方介绍](https://www.springernature.com/gp/snapp)

### 3.2 Frontiers三刊：主题匹配与收稿门槛分别判断

本次补查改变的是对Frontiers的优先级判断：**不能把“过去发表过类似文章”当作现行政策下纯公共数据模型仍容易接收的证据，也不能把Genetics称为保底。** 出版商2026年发布的2025 JIF汇总表列为：Pharmacology 5.4/Q1、Oncology 3.4/Q2、Genetics 3.0/Q2。这里的Q是其官网展示的JIF rank quartile，不是中科院分区；若涉及单位考核，应再核对指定年份与学科类别，不沿用旧约数。[官方指标及年份说明](https://www.frontiersin.org/about/impact)

| 期刊/建议专栏 | 主题关联 | 核查后的主要障碍 | 对本月纯公共数据方案的判断 |
|---|---|---|---|
| Frontiers in Pharmacology；Pharmacogenetics and Pharmacogenomics | 药物反应与表达基因组学相关，专题方向确实接近 | 全刊要求适当验证，并排除没有药理机制/药物开发新见解的纯统计或预测研究 | 有条件候选，当前不排在BMC Genomics之前；严谨预测比较本身未必满足其药理学门槛 |
| 同刊；Experimental Pharmacology and Drug Discovery | 名称包含药物发现，但不代表接收所有预测算法 | 专栏明确不收仅in silico研究 | 本月以公共数据和计算为主，不推荐此专栏 |
| Frontiers in Oncology；Cancer Genetics | 肿瘤与转录组主题相关 | 全刊列明非公开独立临床/患者队列或体内/体外生物验证；Cancer Genetics专栏还明确排除无实验或原位验证的纯生信 | 公开GDSC→CTRP不能自行认定符合这一要求；不列入本月Top3 |
| 同刊；Cancer Immunity and Immunotherapy | 当前项目没有免疫治疗结局或免疫机制主问题 | 专栏问题与现有证据不匹配，且仍受全刊要求约束 | 不为投稿强加免疫分析，不推荐 |
| Frontiers in Genetics；Computational Genomics | 方法、表达分析、药物反应题目最接近该专栏 | 全刊对公共数据纯生信的验证和创新有明确限制；常规nested CV/SMOTE不是算法创新 | 若未来考虑Frontiers，先核对具体研究类型适用性；不可作为自动保底 |

官方依据：[Pharmacology总范围](https://www.frontiersin.org/journals/pharmacology/about)、[Pharmacogenetics and Pharmacogenomics](https://www.frontiersin.org/journals/pharmacology/sections/pharmacogenetics-and-pharmacogenomics/about)、[Experimental Pharmacology and Drug Discovery](https://www.frontiersin.org/journals/pharmacology/sections/experimental-pharmacology-and-drug-discovery/about)、[Oncology总范围](https://www.frontiersin.org/journals/oncology/about)、[Genetics总范围](https://www.frontiersin.org/journals/genetics/about)、[Computational Genomics](https://www.frontiersin.org/journals/genetics/sections/computational-genomics/about)。

专栏依据另见：[Cancer Genetics](https://www.frontiersin.org/journals/oncology/sections/cancer-genetics/about)、[Cancer Immunity and Immunotherapy](https://www.frontiersin.org/journals/oncology/sections/cancer-immunity-and-immunotherapy/about)。Pharmacology与Genetics现行Original Research/Methods均列A型标准APC CHF3,150，税费及协议另计；费用不是本轮排序的主要矛盾。[Pharmacology费用](https://www.frontiersin.org/journals/pharmacology/for-authors/publishing-fees)、[Genetics费用](https://www.frontiersin.org/journals/genetics/for-authors/publishing-fees)

Frontiers通用方法政策列出独立临床/患者队列、体内/体外生物验证，另将多数据库验证限定于文献计量研究；Oncology本刊范围比通用页进一步写明临床/患者队列应为非公开。两者不能混成“所有Frontiers只收非公开临床数据”。同样，不能把已有公共细胞系药敏测量重新命名为本项目新完成的独立生物验证。**我的判断是存在显著范围风险，而非声称这三刊绝不可能接收任何纯计算方法稿。** 若实际获得符合要求的验证、或编辑针对完整研究摘要确认适合，再重排；本次没有代用户联系编辑。[Frontiers方法政策](https://www.frontiersin.org/guidelines/policies-and-publication-ethics)

### 3.3 速度、系统和最终顺序

同一Frontiers平台确实能减少格式切换，但相同系统不消除专栏和全刊的证据门槛。当前Frontiers作者指南给出的全平台平均投稿至接收为96天，生产指南称最终批准后通常2–3周发表；另一指南仍写平均审稿少于90天且未清楚定义终点，不能将这些数字拼接成三本目标刊的速度承诺。96天是平台平均、BMC的3/4天是本刊首决中位数、Scientific Reports的138天是本刊已接收稿投稿至接收中位数，覆盖对象和统计量不同。没有当前同口径、同稿型数据，就不写“Frontiers肯定更快见刊”，也不由全刊历史接收率推算本稿成功概率。[Frontiers 96天指南](https://www.frontiersin.org/for-authors/guide-to-platforms/using-peer-review-forum)、[生产时长](https://www.frontiersin.org/for-authors/submitting-research/what-is-the-production-process)、[另一审稿指南](https://www.frontiersin.org/for-authors/submitting-research/what-is-the-peer-review-process)

BMC Genomics与BMC Medical Genomics的研究稿指南均容纳新的计算方法；本次查阅的范围、研究稿指南与共享编辑政策未见上述Frontiers式公共数据纯生信的统一实验验证条款。这只说明证据类型的范围障碍较少，不代表简单拼接公共数据、缺乏新信息或验证不足也能接收。[BMC Genomics研究稿要求](https://link.springer.com/journal/12864/submission-guidelines/research-article)、[BMC Medical Genomics研究稿要求](https://link.springer.com/journal/12920/submission-guidelines/research-article)

因此，在**1人、1个月、公共数据、无已落实新增湿实验或临床队列** 的条件下，Top3保持BMC Genomics、Scientific Reports、BMC Medical Genomics。依据是对新稿证据类型的相对匹配；不是对录用速度或结果的保证。若坚持考察Frontiers，优先核对Genetics–Computational Genomics的研究类型适用性，其次Pharmacology–Pharmacogenetics and Pharmacogenomics；Oncology及Experimental Pharmacology专栏不作为这一个月方案的默认目标。不要同时在多刊正式投稿。

## 4 重构后的研究设计

### 4.1 主研究与三层证据

**主要研究：** 每个药物建立跨癌种的连续药物反应预测模型，输入未经药物处理的基线表达，主表示为REO，和普通表达量及癌种基线公平比较。分析单位为独立癌细胞系，不是药物–细胞系测量行、培养重复或单个细胞。

**主要结局：** 需在模型阶段完成单位/方向验收的连续药敏曲线AUC，统一约定数值越大表示体外敏感性越低；各来源原始定义和剂量范围单独保留。正文写`drug-response AUC`，分类指标写`AUROC`，避免混淆。不要称连续AUC某一高值自动证明“获得性耐药”。

**主要比较：** 事先指定200候选基因形成的REO elastic net，与2,000高变基因表达量elastic net的药物内预测关联之差；同200基因的表达量模型只作表示匹配消融。RF是必要的非线性基线和辅助REO模型。模型胜负不能靠外部结果挑选。

| 证据层 | 数据与用途 | 能支持的结论 |
|---|---|---|
| A 主证据 | 开发集重复嵌套CV；冻结的CTRP细胞系测试 | 对未见过细胞系的预测，以及独立药敏筛选来源的迁移 |
| B 次要证据 | 按癌种评估；组织留出；表征/特征预算/端点敏感性 | 哪些条件有信号、组织差异贡献及模型失效边界 |
| C 探索证据 | GSE38376，GSE156246 bulk与单细胞；有限富集 | 与处理或持留状态的关联；用于提出后续实验假设 |

不再把A层细胞系敏感性、C层获得性耐药/药物耐受持留细胞、患者疗效混为同一个结局。未经实验验证，不写“机制确认”“耐药驱动基因”“可用于伴随诊断”。

### 4.2 药物与样本纳入规则

第1–4天只根据药物身份、表达可用性、样本覆盖、已知检测质控和数据来源设计筛选；不查看候选模型在外部数据的相关性、P值或ROC。

1. 建立全部GDSC2与CTRP重叠化合物台账，优先相同化合物、可辨认盐型/结构、非混合处理和明确剂量反应端点。禁止仅按模糊名称匹配。
2. 本轮按覆盖规则确定**8种主分析化合物**；队列v1锁定后不按模型效果扩大或替换。按预定规则排序：可用独立细胞系数、冻结外验覆盖、组织覆盖；并列按固定化合物ID。生物学类型只作预定分层，不按预测高低决定纳入。
3. 按独立生物组计数：开发≥200、组织≥4、外验最低≥50；足够8种均达到≥80时优先在该池排序。本轮58种达到80目标，故主面板全部≥80。**这些是本项目的可行性目标，不是统计学保证或普适标准。** 第4天根据实际N与预期CI宽度确定最终任务集。
4. 外验N不足或表型变异很小的药物不能偷偷删除；若冻结前已不达预设纳入条件，则进入可机读排除表。冻结后失败保留结果、原因及覆盖率。
5. 原始逐药物–癌种分类器不作为主分析。只有确有足够独立样本、两类且评价可解释的条件才作次要分析；原数据最大N51，不能假装提高CV重复次数即可解决不足。
6. Lapatinib可作为预定生物学案例；如果不满足主任务覆盖，明确列为探索案例，不为保留旧图而放松主分析标准。

以80个独立测试样本、观察相关约0.30为例，简单Fisher-z近似给出的区间仍约0.09–0.49；N50时约0.02–0.53。这只是相关系数精度的粗略量级示例，实际Spearman分布和聚类结构需用适合的重采样/模拟检查，不能用这个例子代替样本量论证。

### 4.3 要重新下载什么、保留什么

**执行更新：** 必要补充已下载，完整本地TPM已确认。最新状态以§0及下载登记为准；下表扩展并非本月全部必做。

| 数据/产物 | 操作 | 科学理由与一月内优先级 |
|---|---|---|
| 本地GDSC2拟合药敏表 | 优先保留原始文件，核对版本、字段、测量单位、唯一ID、重复曲线、质量标签与SHA256；缺少元数据才补下载 | 旧文件不因年份旧而自动无效；首要是来源可追溯。主表无需从全部原始曲线重新拟合起步 |
| 完整基线CCLE/DepMap表达与Model注释 | 先核验2018 counts/TPM是否完整；推荐冻结一个来源说明清楚的公开TPM/log2(TPM+1)发布版本与配套Model表，补下载缺失的版本文件 | 不能从原“NR/R条件子集”还原连续任务；需要全体可用细胞系。版本选择依据覆盖/数据定义，不根据模型成绩 |
| DESeq2用counts | 主模型不再强制DEG门槛；若做探索DE，另外核实同版本原始/期望counts及恰当导入方式 | TPM、FPKM、log表达不能冒充原始counts。RNA-seq差异表达与跨基因排序使用的量纲要分开 |
| CTRPv2药敏/化合物/细胞系注释 | 重新整理完整可追溯表；本地文件完整可复用，缺药物结构/浓度/身份则补元数据 | 核心外部证据；清除别名误配及训练个体重叠 |
| GDSC1 | 可作额外筛选来源复测，优先级低于CTRP | GDSC1/GDSC2是实验体系差异，不能写成随机拆分的两个独立患者队列；需同样做ID隔离 |
| PRISM secondary dose-response | 第14天前若主流程已稳定、共有药物和N充足，再补下载；不列为最低交付的强制项 | 第三来源可提高说服力，但primary单浓度viability与dose-response AUC不能混用；先查质量和身份匹配 |
| Sanger表达测量 | 仅作为可选的真正跨表达测量来源测试，审查merged表的来源标记 | 现行Sanger整合表可能同时含Broad来源，不能把下载站点不同当独立表达平台；见官方说明 |
| GSE38376 | 保留矩阵、补齐GSM处理/重复注释和平台注释，不必重下全部原始array | 仅一对亲本/耐药衍生系；按处理状态描述，不充当临床外验 |
| GSE156246 | 优先补齐39个样本的metadata和既有bulk计数/表达数据，核查本地单细胞矩阵是否完整 | 除已用六个单细胞条件外，还含四个细胞系的bulk处理与亲本样本；作为补充探索更有针对性 |
| 全部FASTQ/SRA、多组学、TCGA生存 | 本月不作为必需下载/分析 | 大幅扩大任务但不能直接回答药物预测问题；TCGA OS不是药物响应金标准 |
| DEG、pair、模型、评估表、主图、机制富集 | 主要结果全部重算；旧结果放只读归档 | 不在泄漏筛选出的特征上补做CV后声称问题已修复 |

数据入口与核验：[Broad历史CCLE文件](https://data.broadinstitute.org/ccle/)、[DepMap数据页](https://depmap.org/portal/data_page/)、[Sanger药敏定义](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/)、[Sanger表达来源说明](https://depmap.sanger.ac.uk/documentation/datasets/expression/)、[GSE38376](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38376)、[GSE156246](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156246)。本次传统CancerRxGene下载页工具访问失败；执行时从官方现行入口定位版本，不以第三方镜像静默替换。

### 4.4 两种外部测试设计与选择顺序

**执行更新：** 队列v1采用自然CTRP-only已知生物组。下面共有细胞保留为未启用备用设计。

**优先检查自然外部独有细胞系是否足够。** 如果CTRP-only且表达可用的细胞系足以覆盖预定药物，可直接作为冻结测试。若数量不足，采用事先保留的共同细胞系设计：

- 在所有药物共同的canonical cell-line全集上，按组织分层、固定随机种子预留约20%的共同细胞系。比例由第4天可用N检查确定，不按结果调整；保存唯一`split_manifest`。
- 这些细胞系在GDSC中的所有药物观测、特征筛选、拟合和调参中均排除，不只是当前药物的一行。其CTRP结果只在流程冻结后用于测试。
- 同时封存这些共同细胞系的GDSC响应。冻结后在相同ID上分别评估GDSC与CTRP，以区分“新细胞系难度”与“药敏来源改变的难度”；CTRP-only个体单列，不能把不同细胞组成的差异都归因于screen变化。
- 剩余共同细胞系可进入开发，但它们的CTRP数据不进入“新细胞系外验”。若展示同一细胞系跨筛选的一致性，单独标为cross-assay reproducibility。
- 各模型使用同一外验名单。测试表达只做预先规定的逐样本变换和预测，不能联合标准化、拟合ComBat、选择可变基因、补齐训练缺失、寻找最佳阈值或重新调药物方向。
- 使用同一CCLE表达配不同药敏库，只能证明**跨药敏筛选来源** ，不能声称跨转录组平台。真正跨表达来源的测试需另外核验样本/测量来源。

**不能把整个CTRP细胞系全集都剔出GDSC** 而不检查剩余规模，两库高度重叠，可能导致开发集过小。若自然独有/预留设计都不能支撑足够N，降低主药物数量或调整研究范围，不能降格后仍写成强外部验证。

**历史使用必须披露：** CTRP已用于旧稿。本次冻结的是重新设计后不再参与选择的测试集，并非从未被研究者查看过的前瞻性盲法队列。历史使用不自动构成本次代码泄漏，但留下研究者选择的影响，需在局限中说明。若时间与覆盖允许，PRISM中按覆盖预先选定的少量药物作为此前未使用来源，其优先级高于重新包装GEO或PPI；不能看CTRP结果后再挑PRISM中表现好的药物。

### 4.5 REO、预处理与对照模型

主特征从可追溯TPM或log2(TPM+1)构造，同一严格单调变换保持样本内顺序；但gene-wise标准化、基因长度校正差异、基因特异测量偏倚和dropout会改变顺序。**REO不等于无需处理、无需质控或自动免疫批次效应。**

预定起始配置：在每个训练折内过滤低检测基因、按MAD选前200个基因，生成19,900个无向基因对；300基因/44,850对作为预算敏感性方案。选择共同基因时可依据平台注释，不得用外验表达分布筛选。每个内层训练折都重新计算这些筛选，不缓存全数据计算出的基因名单。

有向编码事先规定：A>B为1，A<B为−1，A=B为0；相等（包括双零）是tie，缺失基因是missing，两者不能混用。保存每对的tie率、缺失率及顺序方向。主RNA-seq数据应采用预定共同基因集减少结构性缺失；超过预定特征缺失比例（建议5%，冻结前决定）的外部样本标为不可评估并报告覆盖率，不凭测试集均值/众数补值以挽救成绩。

每个训练折内去除恒定pair，固定高tie过滤规则（起始建议tie比例>90%剔除，pilot后在结果评估前冻结），记录低变化率。基因对按稳定Ensembl ID规范排列并保存A/B方向，不能由输入文件行序决定。200候选基因是本次方案边界，不能把结果推广为所有REO候选生成方法均有效/无效。

低于弃判门槛的零星missing，使用该特征在相应训练折中拟合并保存的众数/中位数规则填补；表达量基线同样使用各自训练折逐基因中位数。缺失不改写为tie，不按新样本其他基因对的均值填补。主要结果另外报告完整可测样本敏感性，以及所有模型共同可评估集合；若填补对结论影响大，则不能宣称稳健迁移。

| 模型 | 输入与目的 | 调参/公平比较 |
|---|---|---|
| B0 训练均值 | 每药只预测开发训练部分平均反应 | 下限基线；常量预测的Spearman不可定义，应记NA，不写成0；用MAE/RMSE评价 |
| B1 癌种基线 | 组织one-hot的正则化模型/收缩组织均值 | 测组织差异能解释多少；未见组织回退训练全局均值 |
| B2 表达量elastic net | train-fit标准化的logTPM；2,000高变基因与组织协变量 | 唯一主要表达比较对象；同200基因版本为预设消融，不替代主基线 |
| B3 表达量RF | 同一预设来源的2,000高变基因与组织信息 | 防止只与较弱线性模型比较；同样调参预算 |
| M1 REO elastic net | REO特征+同样组织协变量 | 预先指定的主要REO模型；在训练折内正则化，选择规则固定 |
| M2 REO RF | 训练折内筛出的pair+同样组织信息 | 检验非线性作用；不得在外验后选为“主模型” |

为控制计算，RF前的pair排序在各内层训练折内按与连续反应的关联强度进行，最多选500对；该筛选规则对所有折、所有药物一致。它是预测筛选，不把筛选P值当确证性生物发现。elastic net的整个候选pair集合可由正则化处理。探索性旧分类分支的DEG/Fisher/SRP可保留为消融，但所有监督步骤必须嵌套；只认一个方向的SRP还须纠正为两种方向对称评价。

**为什么不照搬“DESeq2→DEG→Fisher/SRP→RFE”作为主流程：** 这些步骤回答二分类组间差异，不能原样套用连续主结局；多层筛选还会增加小样本不稳定性与计算量。DESeq2不是使用转录组机器学习的必经步骤，`|log2FC|≥2、FDR<0.05、SRP≥70%`也不是普适发表标准。若预定案例做旧流程消融，使用可追溯counts、训练折内二分阈值与可估的设计式（如`~ tissue + batch + sensitivity`，以元数据和非混杂性为前提）；DEG、Fisher多重校正、双方向SRP、RFE特征数及RF调参均在每个内层训练折重做。可选择锁定门槛或把有限门槛纳入内层网格，不能看外验后调宽；空特征折按预定回退规则保留。解释性DE结果仍不等于机制验证。

RF建议300棵树固定，`max_features`与`min_samples_leaf`共不超过8组候选；elastic net以固定的少量L1比例和正则化路径调节。具体参数及软件版本第7天前在`analysis_config`中冻结。不要为了复刻原稿数值去试500/800/1000树后选择外验最好的配置。额外复杂模型不属于本月必需项。

所有EN模型对组织协变量采用相同处理，推荐组织项不受惩罚、基因/基因对项正则化，并记录参考组与未见组织回退规则。RF增加训练内控制组织后关联排序的预设敏感性：只在筛选统计量中去除组织效应，REO编码仍从原始统一量纲表达生成，不能把gene-wise组织残差当REO输入。加入组织协变量本身并不能保证消除组织混杂。

### 4.6 嵌套评价、指标与不确定性

开发集采用**5折外层×3次重复，内层3折** 。所有划分按canonical cell-line/可识别亲缘模型组进行，同源重复和衍生模型不能分跨训练/测试。多个药物共用同一全局细胞系划分，以免模型选择间接接触外部个体。尽量平衡组织覆盖，不能保证的稀有组织写入记录。

**64%训练/16%验证/20%测试的一次划分是三分留出，不自动成为nested CV。** 数字比例本身不能说明是否嵌套；关键是外层反复估计泛化、内层选择配置，且每个内层训练折从原始输入重新拟合整条流程。仅在整个外层80%上先做DEG/RFE，再用它们的结果跑内层CV，仍会让内层验证信息参与选择。该错误不等于外层测试已泄漏，但会扭曲内部模型选择；应统一用完整pipeline纠正。[nestedcv原始论文](https://doi.org/10.1093/bioadv/vbad048)

每个外层的顺序严格如下：

1. 留出外层测试，禁止读取其结局参与任何选择。
2. 在外层训练中建立内层训练/验证；每个内层从低表达过滤、MAD筛选、pair构造/筛选、标准化/缺失处理开始拟合。
3. 仅用内层指标选超参数，使用同一预设目标与同一搜索预算；训练折无可用特征时回退到预先定义基线并记失败，不能丢弃该折。
4. 最佳配置在整个外层训练重拟合全部步骤，预测外层测试并保存每一条out-of-fold结果。
5. 完成开发后仅根据开发证据选定最终配置/预定模型族，在全部开发数据重拟合；冻结代码、输入、特征、方向及模型哈希，才打开外部标签并一次性评估。

**内部主报告：** 每药每次重复的OOF Spearman、MAE、RMSE、R²（允许为负）；报告分布及样本数。训练均值可作为误差基线；不把所有药物AUC直接堆叠算一个很高的相关系数。

**冻结外部主报告：** 每药Spearman及95%区间、预设REO与表达量模型的配对差、覆盖率、全部失败结果。由于GDSC和CTRP的剂量范围、AUC积分定义可能不同，未完成事先规定的同剂量区间标准化前，不直接用跨库RMSE/R²声称准确回归绝对药敏值；也不在测试库做分位数映射来提高一致性。

主要配对差必须在两模型均可评价的共同ID集合上计算，并同时报告全部外验ID、各模型覆盖率、交集N和排除原因。常量响应/常量预测导致相关不可定义时记NA和原因，不将其静默移出“成功率”分母；主报告保留逐药结果，不依赖删除失败任务后的单一宏平均。多个药物的主要配对推断须预先定义检验族与多重性处理。CV每次重复分别形成一份完整OOF预测；重复间不直接堆成更大的独立测试集。

**组织混杂：** 比较REO+组织与组织基线、表达+组织；展示样本充足组织内相关与误差，补充leave-one-tissue-out。组织留出评价的是对新组织的迁移，不等于新药迁移；单药模型不能宣称预测未见药物。

**区间与统计：** 外部固定模型的配对比较以独立细胞系/生物组为重采样单位，建议2,000次bootstrap；跨药物汇总时同一次全局ID重采样保留同一细胞系在不同药物中的相关结构。该区间主要反映固定模型在测试群体中的不确定性，不涵盖全部训练算法变异。重复CV各折相依，禁止把15折当15个独立研究直接算普通t检验CI。CV重复间波动与测试bootstrap区间分开命名。

**置换检验两部分：** (a)最低交付锁定2个按覆盖而非成绩选择的药物，只针对主要REO-EN管道，采用外层3折一次、内层3折及缩小但预先固定的EN网格，进行99次全管道标签置换；观察统计与置换使用这套完全相同的小协议，每次重新特征筛选/调参。是否扩为199次在看置换结果前按运行预算决定；RF置换为可选。(b)冻结外部预测与标签之间的关联可做999次或更多ID层面的置换；这是固定预测的外验检验，不冒充全流程置换。全局标签置换检验“没有可预测关联”，不能检验“REO与表达模型同样好”；后者主要用外部配对差及bootstrap。要检验超出组织的信号时使用恰当的组织内交换/残差置换并说明交换性假设；组织过小则不强行检验。所有P值用(1+b)/(1+B)，预定检验族做BH校正。

把冻结测试标签打乱1,000次，是可以定义的关联检验，但固定模型的预测不会因标签改变而改变，通常无需每次“重新预测”。它不能检查DEG/RFE/调参产生的全流程择优效应，更不能修复既有泄漏。99次全流程置换的最小P值为0.01，只提供有限分辨率；不能报告P=0或用低分辨率结果作过强显著性宣称。软件中的`permutation_test_score`需包装整个拟合pipeline和嵌套调参过程，单独传最终RF不等于已检验完整方法。[官方置换检验文档](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.permutation_test_score.html)

实现时分开“CV分组”“重采样单位”“置换允许的交换层”。不能把每个样本唯一的cell-line ID直接当作该函数限制置换的`groups`：组内只有一个样本就无法打乱标签。独立细胞系ID是重排单位；若采用组织内置换，组织才是交换层，同源模型的依赖另按预定块结构处理并验证。保存每次置换索引，断言实际发生重排。

### 4.7 阈值问题如何收口

连续主任务直接消除“人为二分决定主结论”的依赖。分类仅作为连接旧研究的预设次要分析：训练集药物内中位数为主阈值，另评估训练集上下三分位极端组与原几何/瀑布规则。阈值在各训练折重新估计、冻结后施加到相应验证/测试数据；不得用全数据或外验中位数重定义类别。

上下三分位会排除中间样本，必须报告保留比例，不能拿其较高AUROC直接证明优于全人群任务。跨库AUC定义未对齐时，绝对分类阈值不能直接移植，外部以连续排序评价为主。分类报告AUROC、AUPRC及其阳性率基线、balanced accuracy、MCC，必要时Brier score；固定Sensitive/Low-sensitivity含义和概率方向，不再用含混的R/NR。

分类支线加入简洁kTSP/基因对投票对照，预先固定少量k候选（例如1、3、5），pair方向、特征和k均在训练内选择；它回答复杂分类器是否优于简单规则。kTSP不替代连续主任务的EN/RF比较，优先在预定案例完成，避免扩大成第二套全药物项目。其必要性来自直接REO药敏先例（文献R17、R18）。

**类别不平衡与SMOTE：** 本计划主任务是连续回归，不需要常规分类SMOTE。次要分类先报告真实类别比例，以不加权模型为基础、训练内类别权重为预设敏感性；如使用随机过采样，仅在内层训练中复制真实样本，且在统计DE/Fisher之后实施，绝不能把复制样本当独立重复做差异检验。验证、外层测试和外部测试维持自然类别比例；概率校准若开展也只能在训练内部完成，并在自然分布测试集评价。

普通SMOTE对−1/0/1的pair做连续插值，可能生成0.3等不对应离散排序的值，因而不默认用于REO。直接取符号也不构成通用解决方案：没有tie的两个±1向量用同一个插值系数时，取符号通常退回其中一个原样本；而逐pair独立合成或类别众数合成又未必满足所有基因间排序的一致性。SMOTENC/SMOTEN不会自动证明该生物表示有效。此处是由算法定义及REO结构作出的判断，不是声称任何数据增强都必然无效。[SMOTE官方定义](https://imbalanced-learn.org/stable/over_sampling.html)、[重采样泄漏示例](https://imbalanced-learn.org/stable/common_pitfalls.html)

重抽样改变训练权重或评估方式，**不会把20个细胞系变成200个独立细胞系** 。本地1,923条件匹配后N中位16的根本应对是调整分析单位、减少模型复杂度、获得更多独立数据并诚实呈现区间；下采样还会损失已有信息，不能作为小样本的默认补救。基线的作用是公平检验增量价值，结果允许简单模型获胜，不以“突出复杂REO流程优越性”为预设目标。

### 4.8 生物学解释、GEO及单细胞的去留

**GSE38376：** 18个样本来自同一SKBR3亲本与耐药衍生系，包含三种短时处理状态；未处理3+3最接近基线状态比较。只展示冻结模型的分数、可测pair比例和逐样本预测。不能通过把同一细胞系的重复随机拆分重训后，声称外部效能从0.67“恢复”到0.83。删除本地重训作为主结果；如果保留学习实验，明确它不能提供独立泛化证据。[GEO设计](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38376)

**GSE156246：** 主稿最多一个探索面板，主模型结果不稳时整块移到补充。先用实际GSM/生物重复形成pseudo-bulk；若一个条件只有一个独立捕获样本，则没有支持条件间总体推断的生物学重复，细胞bootstrap或人工分块不能补出重复。双零/低检测率必须显式记录，不用未经验证的插补创造基因顺序。样本内技术下采样只评估测序深度敏感性，不等价于生物学CI。基于分数的“高低”不能当已知耐药金标准，更不能由处理时间的单调上升推断动态预测效能。原研究的bulk样本可补充状态方向观察，但仍不证明患者疗效。[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156246)、[Chang等关联论文记录](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA657182)

pseudobulk具体按“独立生物样本×预定义细胞类型/群体”求和原始UMI，再进行样本级QC和事先规定的归一化；GSM是记录标识，必须核查其是否代表真正独立生物重复。不能把同一捕获样本随机拆成10份，就得到n=10。聚合降低稀疏性，但不能消除bulk TPM与单细胞UMI的基因长度和测量差异，也不能自动证明原REO模型可移植。UMI→CPM后再log并不等于bulk TPM；需报告pair可测率、tie率、深度敏感性和表达尺度差异，结论保留为状态探索。[单细胞重复与伪重复](https://doi.org/10.1038/s41467-021-21038-1)

**scDEAL已纳入文献R19，但本月不设为必做。** 它是有价值的bulk–single-cell深度迁移框架，不是修复原REO模型的即插即用步骤。加入它会形成另一套模型、超参数和比较问题；它利用目标单细胞表达作表征学习/域适配，因此该目标表达已参与训练过程，不能再称“完全未接触的单细胞外验”。如后续实施，应说明迁移设置、分离适配与评估样本、按真实样本分组、设置相同源/目标数据预算的对照，且不用目标真实标签调参。缺少独立耐药金标准时，输出仍是预测分数，不能用它验证自身。当前没有增加独立标签和重复，仅增加DaNN复杂度不会解决Major 6。[scDEAL原始论文，2022](https://doi.org/10.1038/s41467-022-34277-7)

**富集/PPI：** 只对开发阶段确定的稳定特征探索，不按外验成绩挑特征；背景为实际进入筛选的可测候选基因，绝不能任意使用全基因组。记录重复选择稳定性、gene frequency和共享组织因素，必要时匹配表达水平/网络度数的随机基因集。做GSEA需要预先定义完整的排序统计，不能把已经筛过的几十个显著基因当全基因排序。PPI hub是数据库网络中心性，不是耐药机制验证；预测验证与机制验证分开。

## 5 一个月的执行安排

设备抽查为i5-13500H、16逻辑线程、约15.7GB内存，F盘约269.7GB可用；因此先以1–2个药物作实测，不假设所有线程并行就更快。建议单药块处理、最多2个进程，避免模型内部与外层同时并行。pair使用紧凑数组，日志保存峰值内存和用时。

| 时间 | 工作与交付 | 当天应能回答的问题/验收门槛 |
|---|---|---|
| D1–2 10/10–11 | 只读归档旧版本；canonical sample/drug字典；完整数据版本表；复核原错误 | 哪些独立细胞系、化合物和表达版本真实可用？哪些旧结果撤下？ |
| D3–4 10/12–13 | GDSC–表达–CTRP交集；固定8个药物与外部ID；QC及排除台账 | 主药物开发与外验N是否足够？不能满足则现在缩题，不等最后一周补数据 |
| D5–7 10/14–16 | 1–2药物端到端pilot；6类模型；全管道泄漏检查；配置/统计计划冻结 | 一次嵌套评估能否运行？ID、方向、tie、模型药物是否都正确？估算总计算预算 |
| D8–12 10/17–21 | 全主任务重复嵌套CV，OOF预测、模型参数和失败记录 | 所有预定药物是否都有完整结果？REO相对基线差异是否受癌种驱动？ |
| D13–16 10/22–25 | 全管道置换；训练集样本量/特征预算/标签阈值敏感性；组织留出 | 主要结论是否依赖单一选择？置换与正常结果有无异常？ |
| D17–19 10/26–28 | 冻结模型后运行CTRP一次性外验；bootstrap配对差；可选第三来源 | 能否支持外部关联/表示比较？不能按外验失败换种子、换药物或再调参 |
| D20–21 10/29–30 | 选择性GEO/pseudo-bulk探索；必要时全部降补充；稿件定位复核 | 论文应是性能改进、适用边界，还是尚不足以投稿？ |
| D22–25 10/31–11/3 | 5张主图、3张主表；Results→Methods→Discussion→Introduction→Abstract | 所有句子是否对应新结果？Methods是否与配置逐项一致？ |
| D26–28 11/4–6 | 干净环境复现；稿件内部审查；逐条回应原审稿问题；引用和文件校验 | 另一环境能否重建核心表？所有旧问题能否指出已修复证据？ |
| D29–30 11/7–8 | 作者确认、cover letter、期刊格式、资料可用性、提交包 | 可审核的新稿、代码、补充文件与作者信息一致；满足门槛后投稿 |
| 11/9–10 | 下载中断、格式/校样式检查的缓冲 | 不用缓冲期进行新的外部结果驱动模型筛选 |

pilot后计算总任务量：8药物×15外层×3内层×8个RF配置×2种RF表示，约5,760次内层RF拟合，另有最终拟合和线性模型。该数是预算例子，不是已经执行。若实测超预算，**在正式结果查看前** 将药物数固定为8、RF网格压缩、减少非主敏感性；保留外部独立性和全流程嵌套。不得先计算到显著为止，也不牺牲基线公平性。主流程与可重复性优先于单细胞、第三数据源和PPI。

**职责分工：** 主要执行者确认药物身份/数据语义、临床相关措辞、贡献与最终科学解释；助手协助代码、日志审计、数据字典、结果生成、文献核验及稿件一致性。作者对数据和结论负责，不能将未执行的分析写成已完成。

## 6 图表与manuscript的具体重建

工作题目：**Relative expression ordering for cancer drug-response prediction: a leakage-controlled cross-screen evaluation** 。如果最终使用真正跨表达来源的设计且证据成立，再考虑cross-platform；没有就不写。

| 主图 | 内容 | 必须标出的信息 |
|---|---|---|
| Figure 1 | 数据流、独立个体与化合物去重、嵌套训练边界和冻结外验 | 原始测量行与独立样本数分开；所有监督步骤画在训练框内部 |
| Figure 2 | 药物覆盖、组织组成、响应分布与数据重叠 | 开发/外验N，全部预设药物，不只展示“显著”任务 |
| Figure 3 | REO vs表达模型与组织基线的内部比较 | 药物内指标、配对差；重复波动不冒充独立CI；去掉数百条件折线 |
| Figure 4 | 外部每药森林图及全部结果；适度展示散点与组织异质性 | 95%区间、真实N、方向、覆盖、失败；散点案例事先定义 |
| Figure 5 | 特征预算/训练量/组织留出/置换/稳定性 | 以主科学结论最相关的3–4面板为限；单细胞如保留移补充 |

主表1：来源版本、表达尺度、反应定义、独立细胞系、重复、用途与许可。主表2：模型/特征/参数预算和内部主要结果。主表3：冻结外部完整结果与配对效应差。**5图+3表** 方便适配Scientific Reports的8个display items；这是一份设计方案，不是所有期刊的强制要求。

补充材料至少包括：S1全任务纳入排除台账；S2样本/药物身份及split清单；S3完整OOF和外验逐样本预测；S4模型配置、特征列表与符号方向；S5所有敏感性/置换与失败；S6如保留的GEO与富集结果。以CSV/TSV为源，XLSX为可读副本，图件同时输出PDF/SVG与投稿分辨率栅格图。表中不以图片替代数值，保留未四舍五入的计算源表；图统一用Sensitive/Lower-sensitivity或continuous response，色盲友好，标记不同不确定性定义。

**推荐写作顺序：** Results → Methods → Discussion → Background/Introduction → Abstract → Title/Cover letter。数字、样本量和配置从脚本导出；不要先写“显著优于”再找支撑结果。

| 章节 | 论证任务 | 建议内容规模（写作预算，非期刊硬限） |
|---|---|---|
| Background | 说明跨来源药敏预测难题、相对排序的潜力与已知局限、研究空白 | 600–800英文词，结尾明确RQ和比较，不宣称普适创新 |
| Methods | 数据字典、身份去重、分割、每个训练步骤、对照、主次终点、统计与复现 | 1,500–2,200词，细节网格/配置可放补充 |
| Results | 覆盖→内部比较→外部迁移→失效边界；探索最后 | 1,200–1,600词，突出效应与区间，避免只有P值 |
| Discussion | 新结论、与直接竞争工作比较、为何有效/无效、外推限制 | 800–1,100词；承认细胞系、组织差异、剂量范围、选择的药物面板 |
| Conclusions | 一句可检验结论与一项明确边界 | 80–120词，不写临床决策价值已经证实 |

BMC Genomics母版采用结构摘要Background/Results/Conclusions，≤350词；BMC Medical Genomics增加Methods小标题；Scientific Reports准备≤200词不分节摘要，并按现行指南调整题目和正文。均准备数据/代码可用性、资金、利益冲突、作者贡献和适用伦理说明；公开细胞系数据不等于可以虚构伦理批准号。作者列表应核对用户提到的Wang Zhen与最新title page的差异。完整期刊要求见单独核验报告。

## 7 原审稿意见对应的完成证据

| 原意见 | 新版必须交出的证据 |
|---|---|
| Major1 数据泄漏 | 固定split；完整train-only管道；每折特征记录；测试ID零交集断言；不复用旧监督特征 |
| Major2 小样本、高性能、不确定性 | 改药物级单位；独立N；重复嵌套CV；外部CI、置换、简单/强表达基线；报告负结果 |
| Major3 RF参数冲突 | 一份配置、一套执行日志、导出的最终参数表；分类与回归参数分开 |
| Major4 1923→255 | 新全任务台账及筛选流程图；旧档案解释353/332/255/231各层，不混淆文件数与科学任务数 |
| Major5 二分阈值 | 连续主结局；训练折内阈值敏感性；极端组覆盖率；跨库量纲限制 |
| Major6 GEO重训及单细胞 | 主效能移至新细胞系冻结CTRP；删除小GEO重训的泛化主张；GSM单位pseudo-bulk、dropout和可测率 |
| Major7 循环富集 | exploratory定位、正确背景、稳定性、完整统计表；不以PPI/富集验证因果 |
| Minor1 术语/语言 | 全文只用定义明确的一组术语，自动核对AUC与AUROC、样本数和各数据层名称 |
| Minor2 夸张论断 | 删除exceptional robustness、clinical translatability、伴随诊断及驱动机制等未经验证措辞 |
| Minor3 完整开放 | 输入来源/哈希、版本环境、代码、split、配置、特征、模型、OOF预测、所有结果与README |
| Minor4 DESeq2 | 若保留DE：真实counts来源、导入、低表达过滤、设计式与可估协变量；不把变换后的量用作counts |

这是本次新稿的质量检查表，不冒充给原期刊的正式返修回复。投稿cover letter说明新的问题、严格评价和具体发现；不必把多次拒稿史写成新稿贡献。

## 8 复现结构与投稿门槛

建议另建`rebuild_v2/`，原工程只读保存。新项目包含`data_manifest`、`sample_dictionary`、`drug_dictionary`、`analysis_config`、`split_manifest`、`src`、`results`、`figures`、`tables`、`manuscript`和锁定依赖环境；文件名及model manifest必须带确切drug ID、训练样本集合哈希、特征顺序与预处理版本。

必须自动检查：样本ID不相交；同源重复不跨折；模型drug ID与测试药物相同；特征A/B与编码方向匹配；分类概率取正确列；常量预测、只有单类、全部tie/缺失都有明确处理；每项结果可从逐样本预测重算；不因AUC低于0.5或相关为负而反转方向；外部结果不能成为新调参输入。软件测试应集中在这些高风险行为，不能只检查代码能运行。

第21天作三种判断：

- **可以进入投稿写作：** 主要流程无泄漏、所有预定药物结果齐全、外部评价可解释、对照公平、结论有明确内容，数据代码能重建。REO不必在所有药物领先，但必须如实报告何时有用/无用及其不确定性。
- **转为适用边界/比较研究：** REO与表达相当或更差，但实验完整，结果能排除某些常见强宣称，且覆盖与稳定性足够。需重新评估相对现有benchmark的新信息；“阴性”不自动等于可发表。
- **暂不投稿：** 外部N过小、标签/药物身份不可靠、测试被反复用于调参，或主体只剩旧小样本GEO。此时延长一至两周补最关键证据，比再投一轮未经修复的稿件更合理；一个月是管理目标，不能越过科学门槛。

**本次交付是审计与研究方案，并未完成新模型计算、证实新性能或替用户提交稿件。** 接下来第一项实际计算应是建立完整身份与覆盖台账、锁定外验，再运行1–2药物pilot；不应从批量重跑旧Step2开始。

## 9 2018年以后的文献与方案决策对照

以下文献目录在本文件附后，逐项列出题名、年份、期刊、DOI/官方链接和本计划中的用途。文献用于约束具体决策，不能代替本项目自己的验证；较早的DESeq2、GDSC等原始方法论文在正式manuscript中仍应保留必要引用，不因本计划要求2017年后文献而错误删除经典来源。


### 文献如何进入方案

| 本项目需要改变的决定 | 优先证据 | 如何用，不能如何用 |
|---|---|---|
| 从“高AUC泛癌分类器”转为“基因对表征的增量价值、稳定性和跨库可迁移性评估” | R04、R06、R07、R17、R18 | REO和RF本身已有先例；应检验对简单基线的增益，不预设优越性 |
| 重新构建表达、药敏、细胞系和化合物映射，冻结版本 | R01–R03 | 标准化ID和剂量反应定义；不能将两个共享细胞系的数据库称作完全独立人群 |
| AUC连续结局作主分析，二分类作阈值敏感性分析 | R05–R07、R18 | 药敏曲线AUC与ROC-AUC严格区分；跨库AUC也需剂量范围、测量时间、方向核验 |
| 所有监督筛选和调参在训练内部完成 | R08–R11 | 外层评估、内层调参；DEG、Fisher、特征数和模型选择均属于拟合过程 |
| 所有按样本、药物、癌种的排除和失败均报告 | R06、R08、R14、R15 | 无有效特征也是预定义流程的失败结果，不能删掉再只报告成功组合 |
| 展示重采样稳定性、独立样本不确定性和置换对照 | R03、R10–R12 | 外层折/重复间并非独立样本，不能把几十次重复当几十倍样本量 |
| 设置表达EN/RF、REO-EN/RF、kTSP、均值/癌种基线 | R06、R07、R16–R18 | 同一数据拆分、同一调参预算；单药任务不存在可用于预测新细胞系的已知“该细胞系药敏均值” |
| 单细胞降为外推/探索，按独立生物样本聚合推断 | R19–R23 | 分离供者/细胞系/重复单位；数千细胞不能代替独立重复；处理后DTP状态不是治疗前临床药效标签 |
| 富集与PPI解释降级，不作机制验证 | R24、R23 | 用实际受测可入选基因背景，报告数据库版本；独立扰动证据比再次富集更有意义 |
| 公开数据清单、划分、特征、模型、配置和可运行脚本 | R08、R13–R15 | TRIPOD+AI/PROBAST+AI供透明度和偏倚检查借鉴，不将细胞系研究包装成已临床验证研究 |

### 24篇核心文献

**R01. Ghandi M, et al. (2019). Next-generation characterization of the Cancer Cell Line Encyclopedia. Nature, 569:503–508.** DOI：[10.1038/s41586-019-1186-3](https://www.nature.com/articles/s41586-019-1186-3)。

- 落地：下载明确版本的表达矩阵和样本信息，建立DepMap ID/细胞系别名/组织类型映射，保存数据来源和校验值。

- 边界：这篇论文不定义你当前下载文件的单位；raw counts、TPM、log2(TPM+1)必须从实际文件说明确认。使用公开矩阵无需默认重下FASTQ。

**R02. Corsello SM, et al. (2020). Discovering the anticancer potential of non-oncology drugs by systematic viability profiling. Nature Cancer, 1:235–248.** DOI：[10.1038/s43018-019-0018-6](https://www.nature.com/articles/s43018-019-0018-6)。[PMC正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7328899/)。

- 落地：优先核查多剂量secondary screen是否覆盖预定药物；PRISM作为可选第二外部药敏来源。

- 边界：不要把primary单剂量活性当剂量反应AUC；单细胞系条形码混池测定与GDSC测定存在实验差异。重叠细胞系测试的是测定跨平台一致性，不能替代新细胞系泛化。

**R03. Feizi N, et al. (2022; online 2021-11-26). PharmacoDB 2.0: improving scalability and transparency of in vitro pharmacogenomics analysis. Nucleic Acids Research, 50(D1):D1348–D1357.** DOI：[10.1093/nar/gkab1084](https://academic.oup.com/nar/article/50/D1/D1348/6438034)。[PubMed](https://pubmed.ncbi.nlm.nih.gov/34850112/)。

- 落地：作为跨库ID与药物重叠审计入口；用版本化PharmacoSets/来源文件，并核查化合物结构ID，避免只靠药物名称拼接。

- 边界：整合数据库不是新增独立实验；同一测定通过两种门户下载不构成两次验证。文中规模是该版本规模，不能当当前规模。

**R04. Xia F, et al. (2022; online 2021). A cross-study analysis of drug response prediction in cancer cell lines. Briefings in Bioinformatics, 23(1):bbab356.** DOI：[10.1093/bib/bbab356](https://pubmed.ncbi.nlm.nih.gov/34524425/)。[PMC正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8769697/)。

- 落地：至少保留一个锁定跨库评估；先用样本/药物重叠与测定可比性决定GDSC2→CTRP或反向方向，再冻结选择。

- 边界：某库在该文更优不保证在本项目药物子集更优；不能为追求更高外部成绩反复选择训练方向。

**R05. Bayer FP, Gander M, Kuster B, The M. (2023). CurveCurator: a recalibrated F-statistic to assess, classify, and explore significance of dose–response curves. Nature Communications, 14:7902.** DOI：[10.1038/s41467-023-43696-z](https://www.nature.com/articles/s41467-023-43696-z)。

- 落地：若当前AUC来源或拟合口径无法追溯，下载原始dose–response或复用已核验统一处理结果；评估曲线QC与AUC定义。

- 边界：四参数log-logistic模型有形状假设，非单调/双相曲线可能不适用；本月无需为所有库重新发明曲线拟合。

**R06. Bernett J, Iversen P, et al. (2026). Critical evaluation of drug response prediction models with DrEval. Nature Communications, 17:4238.** DOI：[10.1038/s41467-026-72903-w](https://www.nature.com/articles/s41467-026-72903-w)。官网正式发表日期：2026-05-12。

- 落地：参考其LCO方案、基线和公开处理数据，考虑给REO模型实现简洁的基准接口；使用每药物指标与宏平均，并做关键消融。

- 边界：不要声称REO一定优于DrEval；该文主要是连续药敏预测，分类和不确定性功能仍有边界。其发表也意味着“增加无泄漏CV”本身不足以构成方法创新。

**R07. Codicè F, et al. (2025). The specification game: rethinking the evaluation of drug response prediction for precision oncology. Journal of Cheminformatics, 17:33.** DOI：[10.1186/s13321-025-00972-y](https://link.springer.com/article/10.1186/s13321-025-00972-y)。

- 落地：主分析保留连续dose-response AUC，按药物分别汇报Spearman/MAE/R²，再作宏平均；IC50作为同口径敏感性分析，而非默认唯一主标签。

- 边界：AUC并非自动跨测定可比，仍受剂量范围、时间、测定技术影响；论文不能为临床敏感/耐药切点背书。

**R08. Kapoor S, Narayanan A. (2023). Leakage and the reproducibility crisis in machine-learning-based science. Patterns, 4(9):100804.** DOI：[10.1016/j.patter.2023.100804](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499856/)。

- 落地：为每个步骤写fit/transform界限，记录数据划分、信息流及模型说明表；确认DEG/Fisher/阈值/特征数/调参不接触外层测试标签。

- 边界：文中案例不能估计本项目虚高多少；必须用重算结果量化，不预先宣称旧结果下降固定百分比。

**R09. Bernett J, et al. (2024). Guiding questions to avoid data leakage in biological machine learning applications. Nature Methods, 21:1444–1453.** DOI：[10.1038/s41592-024-02362-y](https://pubmed.ncbi.nlm.nih.gov/39122953/)。

- 落地：细胞系重复/衍生株/跨数据库同一细胞系整体划分；外部验证明确是在新测定、还是新生物实体。

- 边界：本文提供问题框架，不提供本研究的最小样本量或固定CV折数。

**R10. Vabalas A, Gowen E, Poliakoff E, Casson AJ. (2019). Machine learning algorithm validation with a limited sample size. PLOS ONE, 14(11):e0224365.** DOI：[10.1371/journal.pone.0224365](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224365)。

- 落地：重复嵌套CV并保存外层OOF预测，避免在小癌种上切出极小单次测试集后只报告点估计。

- 边界：不要将模拟结论理解成嵌套CV可以消除小样本不确定性；也不能由本文推导“n≥50/80/100保证可靠”。

**R11. Lewis MJ, et al. (2023). nestedcv: an R package for fast implementation of nested cross-validation with embedded feature selection designed for transcriptomics and high-dimensional data. Bioinformatics Advances, 3(1):vbad048.** DOI：[10.1093/bioadv/vbad048](https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad048/7117540)。[PMC正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10125905/)。

- 落地：可参考此包搭建R实现；但本项目DEG→REO→Fisher是一条定制流水线，必须逐层检查训练索引，不能因调用nestedcv就默认完全安全。

- 边界：筛选阈值若参与内层调参，其筛选也要在每个内层训练集重做；只在外层训练集一次选好特征再内层调参的配置要明确评估其偏倚。

**R12. Nogueira S, Sechidis K, Brown G. (2018). On the Stability of Feature Selection Algorithms. Journal of Machine Learning Research, 18(174):1–54.** [期刊官网](https://www.jmlr.org/beta/papers/v18/17-514.html)；[正式PDF](https://www.jmlr.org/papers/volume18/17-514/17-514.pdf)。官网未列DOI，勿编造。

- 落地：报告基因对选择频率、Jaccard及可选的校正稳定性指标；最终模型的特征表与OOF过程中各折特征表分别保存。

- 边界：基因相关/互换特征可降低集合稳定性而不降低预测；稳定性高也不等于因果机制成立。PDF明确Published 4/18，符合2017年后要求。

**R13. Heil BJ, et al. (2021). Reproducibility standards for machine learning in the life sciences. Nature Methods, 18:1132–1135.** DOI：[10.1038/s41592-021-01256-7](https://www.nature.com/articles/s41592-021-01256-7)。

- 落地：本月最低完成数据manifest、config、环境锁、固定随机种子、每样本预测、完整特征、训练模型和一键重建主要图表。

- 边界：公开并可运行不等于数据分析正确，仍需R08–R11的信息流审计。

**R14. Collins GS, et al. (2024). TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ, 385:e078378.** DOI：[10.1136/bmj-2023-078378](https://www.bmj.com/content/385/bmj-2023-078378)。

- 落地：借鉴适用条目制成补充报告清单，并标明“adapted for a preclinical cell-line prediction study”。

- 边界：本项目不是临床诊断/预后模型，不能以“遵循TRIPOD+AI”替代临床验证。官网另有作者单位更正，使用更新后正文。

**R15. Moons KGM, et al. (2025). PROBAST+AI: an updated quality, risk of bias, and applicability assessment tool for prediction models using regression or artificial intelligence methods. BMJ, 388:e082505.** DOI：[10.1136/bmj-2024-082505](https://www.bmj.com/content/388/bmj-2024-082505)。

- 落地：在投稿前进行一次适用条目自查，特别关注样本量、缺失、筛选、结果选择与过拟合。

- 边界：源自医疗预测模型语境，细胞系项目借鉴其原则，不宣称完成临床适用性认证。

**R16. Maeser D, Gruener RF, Huang RS. (2021). oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data. Briefings in Bioinformatics, 22(6):bbab260.** DOI：[10.1093/bib/bbab260](https://pubmed.ncbi.nlm.nih.gov/34260682/)。[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8574972/)。

- 落地：作为相关工作与可选表达预测对照；主要基线仍用自己严格训练内拟合的EN/RF，明确跨域处理是否用目标队列信息。

- 边界：oncoPredict的预测值不能用作真实药敏标签，也不能用另一个算法的预测去证明本算法正确。

**R17. He J, et al. (2020; online 2019). Qualitative transcriptional signature for predicting pathological response of colorectal cancer to FOLFOX therapy. Cancer Science, 111(1):253–265.** DOI：[10.1111/cas.14263](https://pubmed.ncbi.nlm.nih.gov/31785020/)。[正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6942442/)。

- 落地：引言承认既有REO方法，分类支线加入kTSP/简洁投票对照，报告签名稀疏性与可解释性。

- 边界：FOLFOX临床反应和单药体外AUC是不同结局；生存相关不能自动证明治疗预测效应。不能继承旧文高性能作为本项目预期。

**R18. Tercan B. (2026). Robust predictors for drug response of patients with acute myeloid leukemia. PLOS ONE, 21(2):e0343422.** DOI：[10.1371/journal.pone.0343422](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0343422)。官网正式发表日期：2026-02-23。

- 落地：必须说明REO-RF相对kTSP的增量问题；可以考察连续AUC、严格跨细胞系验证、特征稳定性是否形成更可信贡献。

- 边界：本研究并非临床疗效试验；其规则亚组一致性部分使用全数据训练规则，不能将该部分当额外独立验证。文中“无需标准化”的措辞不应直接照搬，方法实际进行了TPM转换。不能把本项目更换到AML作为一月内默认必做项。

**R19. Chen J, Wang X, et al. (2022). Deep transfer learning of cancer drug responses by integrating bulk and single-cell RNA-seq data. Nature Communications, 13:6494.** DOI：[10.1038/s41467-022-34277-7](https://www.nature.com/articles/s41467-022-34277-7)。

- 落地：若保留单细胞，只在锁定bulk模型后做探索；有明确标签和独立重复时才比较scDEAL，且分清有无目标域适配。

- 边界：数千单细胞的细胞级分类成绩不能证明独立患者效能；训练/验证处理状态和治疗前预测是不同科学问题。一个月内不强制上复杂迁移网络。

**R20. Luecken MD, Theis FJ. (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. Molecular Systems Biology, 15:e8746.** DOI：[10.15252/msb.20188746](https://pmc.ncbi.nlm.nih.gov/articles/PMC6582955/)。

- 落地：记录每样本UMI/基因数/线粒体等QC；REO统计零值、并列和可判定基因对比例；以真实生物重复的pseudobulk作稳健性检验。

- 边界：CPM和统一单调log变换可保留同一样本内排序；gene-wise z-score、不同基因校正、UMI与全长TPM的基因长度因素可能改变排序。REO并不天然免疫dropout/批次。

**R21. Squair JW, et al. (2021). Confronting false discoveries in single-cell differential expression. Nature Communications, 12:5692.** DOI：[10.1038/s41467-021-25960-2](https://pubmed.ncbi.nlm.nih.gov/34584091/)。[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8479118/)。

- 落地：细胞级可视化与样本级推断分开；按供者/重复/细胞系聚合，不能把同一样本随机切成多个“伪重复”增加n。

- 边界：该文评估DE，不直接证明某种drug-response分类器更好；向预测评估的推断单位原则转用需要明说。

**R22. Zimmerman KD, Espeland MA, Langefeld CD. (2021). A practical solution to pseudoreplication bias in single-cell studies. Nature Communications, 12:738.** DOI：[10.1038/s41467-021-21038-1](https://www.nature.com/articles/s41467-021-21038-1)。

- 落地：划分、bootstrap和置换按生物样本/细胞系分组；若每组只有一个独立样本，单细胞结果只能作描述性案例。

- 边界：统计模型不能凭空补足缺失的生物重复，也不能用随机细胞分组制造独立验证。

**R23. Chang CA, Jen J, et al. (2022). Ontogeny and Vulnerabilities of Drug-Tolerant Persisters in HER2+ Breast Cancer. Cancer Discovery, 12(4):1022–1045.** DOI：[10.1158/2159-8290.CD-20-1265](https://pubmed.ncbi.nlm.nih.gov/34911733/)。[PMC正文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8983469/)。

- 落地：更正原稿数据引文；核查GEO样本表后考虑其bulk和scRNA部分作lapatinib/DTP状态探索，提前定义处理前/处理后及真实独立重复。

- 边界：处理后DTP vs parental的可区分性不等于治疗前药敏预测；持续存活不等于细胞系永久耐药；该组实验不能直接用来验证患者临床疗效。将公开实验结果作为背景证据，不冒充本项目新机制实验。

**R24. Reimand J, et al. (2019). Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap. Nature Protocols, 14:482–517.** DOI：[10.1038/s41596-018-0103-9](https://pmc.ncbi.nlm.nih.gov/articles/PMC6607905/)。

- 落地：只对训练发现的稳定特征作探索富集，用表达QC后可入选基因作为背景，完整报告全部结果；优先少量非冗余通路，不堆砌GO气泡图和PPI网络。

- 边界：富集不能作为对筛选过程的独立验证，PPI边与通路显著性不能证明耐药因果机制；独立队列重现方向或外部扰动证据才能增加证据层级。

文献均按正式发表口径核验，年份为2018—2026；完整核验范围、支持段落类型和更多细节见同目录 literature_evidence.md。这是针对本项目的定向证据包，不是系统综述。
