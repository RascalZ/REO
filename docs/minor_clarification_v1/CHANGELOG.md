# Minor 1 与 Minor 2 的最小补充

日期：2026-10-26。作者明确同意前一轮讨论的最小修改，不接受扩展重写。本轮形成英文通读v4及对应中文修订通读v3；原v3中英文副本保留。

## 实际修改

每种语言仅改两个内容块，其余195个非空块原样保留。

1. **Response outcome 段，中文第5页／英文第8页。** 增加两句：CTRP v2使用发布的post-QC sigmoid拟合曲线字段area_under_curve，保持发布尺度而不作额外归一化；本研究将较大细胞活力曲线AUC解释为较低敏感性，跨筛选比较药物内排序而非绝对AUC。方向表述保留为本研究的解释，没有写成发布方明确给定的尺度等同。没有宣称重新QC、积分或剂量统一，也未新建端点对照表。
2. **图4A图注，中文第15页／英文第25页。** 增加一句：仅在这一80格展示网格内，46个未定义差值均涉及至少一个零分子解模型在留出组织内给出常数预测。没有把这一原因推广到全部23组织，也未添加30/8/8分列表、重画图或改变严格宏值NA。

英文增补原句为：

> External responses used the CTRP v2 post-QC sigmoid-curve field area_under_curve, retained on its released scale without additional normalization. For this viability-curve endpoint, larger AUC values were interpreted as lower sensitivity; cross-screen evaluation compared within-drug rankings rather than absolute AUC values.

> In this 80-cell display grid, all 46 undefined contrasts involved constant held-out predictions from at least one zero-molecular solution.

没有修改摘要、结果叙述、讨论、结论、表格、引用题录或作者声明；没有改变性能估计、区间、样本量或统计政策。四张原图像素不变，图4替代文本同步；中文Word的版本提示更新。其他审稿意见保持原状态，不推定作者已接受或拒绝，不改原审稿卡、机械总评或来源顺序路线图。

## 交付与核查

新版仅置于本地Git忽略目录：

- `manuscript_v1/author_review_v4/Manuscript_BMC_Genomics_author_review_v4.md/.docx/.pdf`，39页。
- `manuscript_v1/author_review_v4_zh/Manuscript_BMC_Genomics_author_review_v4_zh.md/.docx/.pdf`，22页。

补充方法和Additional files继续使用未改的author_review_v2配套文件。原Word的排版、表格、媒体和样式沿用，仅修改授权文字节点；并非重新套用模板。56件文件在本轮前后SHA一致。一个更早的中文v2 Word在本轮开始时已与历史清单字节不同；本轮记录并保护其现状，未恢复或覆盖。临时Word锁文件不计作科研交付。

机械检查通过：每语种2个内容块变化、197块英中数字与引文对应、两表XML及四原PNG不变、Word/PDF均包含新增文字、无超出页面的文本或图像。英文39页、中文22页的全部61张PNG已按原始尺寸分工实际查看，未见需继续改动的排版问题；另有只读范围复核确认没有过度矫正。见build_manifest、mechanical_audit、visual_inspection及final_audit。

技能render_docx.py已实际尝试，但本机缺LibreOffice；随后使用Microsoft Word只读导出PDF及工作区捆绑Poppler渲染PNG完成核查。未安装替代软件。此次检查是内容/排版与来源一致性核对，不是原始科学复现。

本轮未训练、重抽样、评分或读取CTRP个体响应/索引，未改数据、冻结方案、模型或作者声明。下一步继续由作者通读；其他修改不在本次范围。
