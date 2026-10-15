# 阶段20：BMC Genomics结构稿、Word排版与补充文件包审计

日期：2026-10-15
状态：结构稿和文件包已生成并完成格式/渲染检查；仍不是可直接提交的终稿。
范围：文稿组织、Word文件、图表与补充文件；本阶段未重新拟合模型、重算统计量或访问/重跑CTRPv2外验。

## 1. 本阶段交付

| 文件 | 用途 | 状态 |
|---|---|---|
| `manuscript_v1/Manuscript_draft.md` | 证据约束版英文全文 | 保留为独立源稿，未覆盖 |
| `manuscript_v1/Manuscript_BMC_Genomics_draft.md` | BMC Genomics结构化工作稿 | 已按研究论文组织，含作者/声明占位 |
| `manuscript_v1/Manuscript_BMC_Genomics_draft.docx` | 可编辑投稿排版工作稿 | Word成功渲染，27页A4；文件内明显标有草稿/作者待补，不可误作终稿 |
| `manuscript_v1/build_bmc_docx.js` | 可重复构建DOCX的脚本 | `npm run build`通过；两张正文表被构建为Word表格对象 |
| `manuscript_v1/bmc_additional_files/Additional_file_1.csv`–`Additional_file_8.csv` | Tables S1–S8机器可读补充文件 | 顺序与稿件引用一致，文件尺寸约2.7–311 KB |

## 2. 期刊结构与版式核对

本稿按BMC Genomics Research Article页面列出的正文和声明结构整理：Background、Methods、Results、Discussion、Conclusions，以及Ethics approval and consent to participate、Consent for publication、Availability of data and materials、Competing interests、Funding、Authors' contributions、Acknowledgements等声明小节。摘要采用Background/Results/Conclusions三个标签，约270个空格分词，低于350词上限；关键词为6项。补充文件1–8先按顺序在Methods中引用，文件名、格式、标题及内容简介逐条列出。

Word稿采用A4、12磅Times New Roman正文、双倍行距、页码及逐页行号；导出PDF确认27页，首页题名和红色作者占位清楚，Table 1/2是可编辑Word表格。表1跨页时表头重复；表2没有拆坏行。图注集中在文末，Figure 1单独完整显示、Figure 2–4在下一页完整显示，没有拆开的图注或只含末尾几行图注的空白页。Conclusions段落整体分页；Methods的数据治理段含明确的生成式AI披露占位，待作者确认工具/版本及准确范围。版面抽查仅证明当前Word渲染可读，投稿系统转换后的最终版仍应再检查。

官方要求参见[BMC Genomics Research Article guidelines](https://link.springer.com/journal/12864/submission-guidelines/research-article)及[general submission guidelines](https://link.springer.com/journal/12864/submission-guidelines)。其中摘要上限、关键词数量、章节/声明、行距与行号、表格和补充文件规则应以投稿当日官网说明为准。23篇期刊论文的期刊名已按PubMed/NLM的ISO缩写核对并写入正文参考文献；该核查通过NCBI E-utilities文章记录的`ISOAbbreviation`字段完成，对照表保存在`manuscript_v1/nlm_journal_abbreviations.csv`。DOI元数据和引用号覆盖不能代替逐句文献支持审查或最后的参考文献格式审读。

## 3. 图件与表格文件

| 图 | 当前PNG像素 | 文件大小 | 检查 |
|---|---:|---:|---|
| Figure 1 | 6198 × 2522 | 656 KB | 600 dpi；以170 mm、300 dpi等效宽度缩放抽查，标签可读、无裁切 |
| Figure 2 | 5005 × 3153 | 288 KB | 600 dpi；误差区间及文字可读 |
| Figure 3 | 7485 × 3496 | 600 KB | 600 dpi；面板顺序与图注一致 |
| Figure 4 | 7921 × 3252 | 655 KB | 600 dpi；LOTO未定义格与外验部分可辨 |

四张图各自低于1 MB；八个补充CSV各低于0.32 MB。图表源SHA和输出哈希沿用先前图表manifest。缩放预览用于可读性筛查；终稿应在最终图片宽度/投稿系统生成的PDF中再看一次字号和嵌入效果。

## 4. 文献与自动校验

参考文献共26项：23项期刊论文、1项数据集DOI资源和2项官方数据/API说明；期刊论文中21项为2018年及以后。正文序号1–26均被引用，未发现悬空引用或未被引用条目；23个期刊论文DOI与Crossref核对一致；BibTeX结构检查26/26通过，0错误、0警告、0重复。用于复核的脚本与结果为：

- `python src/validate_numeric_references.py`
- `python src/validate_numeric_references.py --manuscript manuscript_v1/Manuscript_BMC_Genomics_draft.md --output manuscript_v1/citation_validation_bmc.json`
- `python C:/Users/Hp/.codex/skills/citation-management/scripts/validate_citations.py manuscript_v1/references.bib --report manuscript_v1/citation_validation_metadata.json --verbose`
- `manuscript_v1/nlm_journal_abbreviations.csv` records the full journal name, NLM ISO abbreviation, DOI, and metadata field used for the 23 journal-article references.

DOCX由`manuscript_v1/build_bmc_docx.js`生成；`npm run build`通过。DOCX结构复核发现2个Word表格对象、238个段落以及页码和行号设置。Microsoft Word实际PDF渲染为27页A4并对首页、表格分页、Conclusions标题和末尾文件/图注页进行了视觉抽查。构建依赖已锁定于`package-lock.json`；当前环境为Node 20时，依赖安装可能给出某个传递依赖要求Node 22的EBADENGINE警告，但本次DOCX构建与Word渲染均成功。

## 5. 尚未关闭的投稿门槛

1. 作者必须提供并确认作者顺序、姓名拼写、单位地址、通讯作者信息、每位作者贡献以及致谢内容。现有红色占位保留，未从旧投稿文件推断或拼接作者信息。
2. 作者需如实核实伦理/同意声明的适用性、基金、利益冲突、数据/代码可用性、软件与上游数据许可、生成式AI使用声明。对纯公开数据计算研究也应按事实写伦理说明，不能自动写“无需审批”。
3. 完成逐条引文是否直接支持所在论断的最后人工审核；期刊缩写已按PubMed/NLM记录核对，但References仍需按最终BMC输出样式整体审读。锁定公开代码快照、许可证、运行环境，并取得可引用的归档DOI。GitHub主分支链接不能代替带版本的归档。
4. 再核投稿日BMC Genomics作者指南、图表编号/图注上限、文件名/大小和系统导入后版面。当前结构稿仍显示“DRAFT — AUTHOR INPUT REQUIRED — NOT FOR SUBMISSION”。
5. 八药评估未新增全流程标签置换分析；这项限制已由阶段17/19记录。是否另行补充须先单独立项并披露其在既有开发结果后提出的探索属性，不能用于改写主分析或再次利用CTRPv2调参。

## 6. 科学结论边界

阶段20仅整理已有证据。开发集REO pair减强表达宏差仍为−0.00588，CTRPv2一次性回顾性跨筛选外验为−0.000372；这些结果不支持REO平均增益，也不构成等效/非劣性检验。保留零解、组织参照、负面或未定义结果。CTRPv2继续封存，任何后续投稿优化都不得依据其性能更换药物、特征规则或统计口径。
