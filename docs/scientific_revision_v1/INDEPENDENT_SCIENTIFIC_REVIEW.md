# R1科学边界与合成报告的独立审读

日期：2026-10-16。对象：`manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md`，并核对已生成的 `manuscript_v1/Supplementary_methods_scientific_revision_v1.md`（SM1–SM8）。本文件是独立于R1主稿/补充文件生成步骤的一轮源证据审读，不是期刊正式外审、独立人体专家评议或对全部原始数据的重新复现。审读者此前整理了合成源证据，因此不主张与此前DGP核查具有独立错误过程。

## 判断

**科学边界通过，建议修正两处方法措辞后继续整合；未发现需要改变合成数值、删掉不利情形或重新运行实验的高优先级问题。** 四个DGP及全部结果可留在补充材料，当前主文的概括与源证据一致。创新定位已收缩到指定候选政策×表示、组织参照、零解与跨筛选任务的可审计比较，没有把研究升级为AIVC机制证实、通用REO优劣定理或原创外验治理框架。

本次审读不替代原始真实响应/OOF审计，也未重验个体预测或model pickle。真实数据的全数值核对由其他来源审计负责；这里的完整数值复核范围是既有合成聚合、生成方程、执行metadata/hash和主稿的推断边界。原合成可报告性审计的24项检查通过；14个冻结输入、7个原manifest输出和run metadata的hash均匹配。未运行拟合、原分析器、置换或重采样。

## 两处可直接修正的文字

| 编号/优先级 | 位置与问题 | 最小修改建议 | 是否影响结果 |
|---|---|---|---|
| S-R1-01，中：术语歧义 | 主稿L59将AUC写为“without ... outcome transformation”；L81及源码实际在训练内标准化target并逆变换。前者本意应是未把source AUC改为二分类、秩或其他科学终点，但字面上排除了实际内部标准化。 | 改为“The continuous endpoint was the source-reported fitted dose-response AUC, without dichotomization or thresholding. Training-only target standardization was used internally during fitting, and predictions were transformed back to the source AUC scale for evaluation.”保留GDSC方向说明。 | 否，准确表达现有代码。 |
| S-R1-02，低：警告范围 | 补充SM3 L62的“zero solver warnings and fallbacks reported for the completed fitting runs”在该段上下文指真实数据拟合，但单独看可能覆盖SM7/8中有267 warning的合成run。主稿L149/L167明确指repeated development，本身没有这个冲突。 | 将“completed fitting runs”限定成“the corresponding completed real-data fitting stages described in this section”，或仅写“splitter warnings were recorded separately from fit-level solver warnings and fallbacks”。不要删除SM8的267警告。 | 否，仅限定指代范围。 |

原子代码证据：`src/ordinal_diagnostic_models.py:272`计算标准化target，`:313`还原原尺度组织预测；主稿L81和补充SM4已给出正确原尺度公式。以上行号对应本次读到的R1版本，后续整合可能移动。

## 合成数值、设计与全部材料核查

| 核查项 | 主稿/补充位置 | 源证据 | 结论 |
|---|---|---|---|
| 四场景各10队列；每队列240组、8×30组织、220基因 | 主稿L117；SM7.1/SM7.2 | `config/synthetic_mechanisms_v1_1.json`；完整40行generation checks | 一致，未把10队列说成同队列10次CV。 |
| 3外折×1次、3内折、8臂；960 task、26,160 fits、76,800 OOF | SM7.3、SM8.1 | execution lock；run metadata；960行outer-task表 | 一致。未与真实八药的5×3设计混淆。 |
| 强表达实际220基因；同预算臂200；pair 19,900 | 主稿L117/L203；SM7.3/SM8.6 | outer-task与240行candidate-set audit | 一致，未声称本合成数据有2000个实际基因。 |
| L=log2(TPM+1)，X=2^L−1；全部噪声截断[-3,3]；响应clip[0.01,0.99] | SM7.2 | `src/synthetic_mechanisms.py`和冻结v1.1协议 | 方程与代码一致。原`log2_tpm`列名数学含义已正确说明。 |
| 各cohort seed与outer/inner公式 | SM7.2/SM7.3；S12 | 冻结公式和40个已记录seed | 一致，四seed区间2026101400–1409、1500–1509、1600–1609、1700–1709。数字seed不代表执行日。 |
| amplitude只缩放TPM+1；各样本order按设计相同 | 主稿L117；SM7.2/SM8.4 | 10/10 rank与pair invariant为True | 准确，不再误称精确TPM乘法或跨平台验证。 |
| 615/840 outer selected molecular zero | 主稿L201；SM8.6 | 完整outer-task聚合 | 分场景210/210、198/210、120/210、87/210相加正确。 |
| 13,575/26,040为全部grid/refit零解 | SM8.6 | `fit_diagnostics.csv` | 与outer selected率分开，未混用分母；含必评解析零端点。 |
| 267 warnings、0 fallback | 主稿L203/L217；SM8.4 | 三连续表达臂各89，全部amplitude场景 | 一致。正确保留警告，未臆断其类别，也未将其等同267失败任务。 |
| 纯组织ρ均值0.9997、全部210 outer molecular为0 | 主稿L203；SM8.2 | `scenario_summary.csv`和outer-task聚合 | 一致。补充保留最大RMSE 1.962615573354719e−17，没有把所有记录强制改成ρ=1。 |
| 弱信号分子均值0.1367–0.1451、组织0.1373 | 主稿L203；SM8.3 | 完整32行scenario summaries | 一致。没有从该DGP推断真实药物一定为低信号。 |
| amplitude ordinal 0.1541、连续0.9830 | 主稿L203；SM8.4 | scenario summaries；所有ordinal outer为0 | 一致且说明按设计的信息不可见性。 |
| pair-order关联表达/秩/pair 0.7573/0.7605/0.9055；MAD 0.2147/0.2069/0.2112 | 主稿L203；SM8.5 | scenario summaries | 与原值四位小数一致，未只挑有利pair情形。 |
| 两driver：关联30/30入选、MAD0/30，strong30/30 | 主稿L203；SM8.5 | `gene_selection_summary.csv` | 一致。正确将其表述为刻意设计的selector stress case。 |
| S11完整320×15，S12完整40×19；三个residual rho列各235 NA | 主稿L310–312；SM8.1 | 两复制表hash与原表一致；全部cohort行 | NA/不适用保留；Spearman、RMSE、MAE、R²点值均320行定义。 |
| failure、资源、版本 | SM8.1/SM8.6；主稿L117/L121/L217 | 单completed session；metadata、manifest、源hash | 记录范围内无failed_attempts文件、0fallback。资源限操作观测；actual package snapshot未记录的缺口未被补造。 |

合成S11/S12的工作复制hash分别为`0c149e110e72e354705ba9513ae5ad3f52f36a54e80ff54322e64e5963f93e77`与`49f6d3c609a9ae93895c12e8f8a2c0b7926d612168a18208ef0b5491f9b7bc4f`。Supplement SM7/8覆盖了源方程、修订、随机数顺序、seed、模型、全部四场景、warning、NA和限制，没有发现科学内容缺项。最终Additional_file_11/12分发复制仍应由打包检查比对相同hash；本次审读时主稿/补充Markdown已生成，但未发现该最终命名的CSV分发副本，属于后续包装状态而非数值错误。

## 核心科学边界核查

| 边界 | 位置 | 审读判断 |
|---|---|---|
| 没有可靠平均优势，不等于无价值/等效/非劣 | 摘要L19–21；Methods L109；Results L159/L175；Conclusion L221 | 符合区间跨0、未设margin的证据范围，没有将“未发现差异”改成信息完整保留。 |
| 真实零解不是拟合失败，非零系数也不是生物机制 | Results L167；Discussion L211；Figure 3 legend | 合理；真实430/840、29/56与合成615/840保持不同分母/任务。 |
| LOTO strict macro NA不应换available平均 | Methods L99/L111；Results L171；Figure 4 | 已明确保留34/80定义性及strict NA。没有以严格宏不可估计直接证明全部组织迁移失败。 |
| 当前八药没有新全流程置换 | Methods L109/L117；Discussion L217 | 如实保留；没有把有限DGP或历史四药九次置换升级为八药经验null。 |
| CTRP历史接触，当前一次性锁不等于项目整体盲态 | Methods L45；Discussion L213 | 已如实说明。共同CCLE表达、不同筛选/剂量的边界明确。 |
| 近零pair−strong不意味着全部信号消失 | 摘要L19；Results L179 | 新增绝对macro组织0.2284、strong0.2602、pair0.2599，明确仅描述性且无macro CI；没有声称统计上证明外部正分子增量。 |
| 候选政策×表示差异不是因果份额 | Results L163；Discussion L209/L215 | 记录描述性差异，明确不能分配幅度、维数、共线性、惩罚的因果份额。 |
| 候选规则不序中立 | Methods L71；Discussion L209 | 对连续表达相关筛选的条件性已保留；没有用它证明一般order-only优势。 |
| AIVC只是相邻任务背景 | Discussion L215 | 正反文献并列，elastic-net零解不等同深度mode collapse，要求直接任务评价；未把本项目称为虚拟细胞检验。 |
| 创新和治理定位 | Background L33–35；Discussion | 当前窄贡献是同组织处理下的匹配候选×表示比较与零解审计。未声称first、首创冻结治理、解决全领域混杂或一般机制证明。合理，但仍是有限八药/单外部screen的比较论文。 |
| 条件CI和重复CV限制 | Methods L105–113；Discussion L217 | 无train-refit coverage主张，重复分折不增加独立生物样本，family9/family18分别控制而非全研究family。措辞科学上审慎。 |

## 尚不能从本轮审读推出的结论

- 不能宣告整项真实数据分析已由新团队从原始数据独立复现；本轮未读取个体响应/预测，也没有重算任何评价指标。
- 合成DGP是明确构造的支持性说明，尤其amplitude与pair-order场景并非自然界机制的鉴别实验。它们提升可解释性，不把观察关联变成机制因果证据。
- 语义/源数值通过不等于BMC Genomics保证送审或接收；八药选择的项目历史、缺失完整算法重拟合不确定性、公开复现路径及作者声明仍是文章的实际边界。
- 作者、伦理、基金、贡献、AI披露和发布许可/获取安排的占位仍未完成；其存在已在稿首标明，不能把本次科学审读当作提交许可。

## 建议处置

保留R1新增的绝对外部表现、秩臂对照及完整合成补充材料。先完成上面两处文字修正，再做表格分发hash和主稿/补充DOCX排版核查；无需为了本轮内容审读重新打开CTRP、重新拟合、补算区间或修改候选规则。
