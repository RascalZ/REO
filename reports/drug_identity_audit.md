# 药物身份审计：GDSC2 → CTRPv2

生成日期：2026-10-10。仅使用药物名称、药物/实验标识、结构、筛选浓度和行数；未查看预测效能，未使用 AUC 数值或表达矩阵选择药物。原始项目文件只读。

## 结论与下游接口

当前有 **64 个 GDSC DRUG_ID 对应 63 个 CTRP 单药**通过本轮保守结构核验，足以支持后续按预先声明的样本覆盖规则选择 8 种主分析化合物，但身份通过不等于样本独立性、终点可比性或预测有效性通过。药理类别仅作描述，未用于事后调整面板。

主接口为 `rebuild_v2/data_processed/drug_identity_map.csv`。下游仅可自动纳入 `identity_eligible=true`，保留 `gdsc_drug_id` 与 `ctrp_compound_id`，不能以药名作为唯一键。`drug_identity_eligible.csv` 是同一接口的合格子集。

这里的 `structure_supported` **不是两个数据库实际药瓶的独立化学鉴定**。GDSC Release 8.5 官方注释只有 DRUG_ID、名称、同义名、靶点/通路，没有可直接比较的 SMILES 或 CID。核验链为：GDSC 官方 DRUG_ID→官方主名称→PubChem 名称解析得到的 CID/完整 Standard InChIKey，与 CTRP 官方 `cpd_smiles` 生成的完整 Standard InChIKey 比较。保留 `gdsc_name_derived_*` 命名，防止误写成 GDSC 直接提供的结构。

## 输入与溯源

| 输入 | 结果/来源 |
|---|---|
| 本地 `1.数据/1.GDSC/GDSC2_fitted_dose_response_27Oct23.xlsx` | 242,036 条拟合记录；读取元数据列，不使用响应值；295 个 DRUG_ID、286 个主名称 |
| 本地 `3.Project/Pan-Drugs-REOs/data/CTRP--Drug--Cellline--AUC.txt` | 545 个药名；原表缺 master_cpd_id；本轮只读取 Drug 和行数 |
| [GDSC 官方 Release 8.5 screened compounds](https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/screened_compounds_rel_8.5.csv) | 621 条药物注释；覆盖本地全部 295 个 GDSC2 ID，主药名全部一致 |
| [CTD2 CTRPv2 ExpandedDataset 原包](https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset.zip) 中 `v20.meta.per_compound.txt` | 545 个 master_cpd_id、545 个药名；含 Broad compound ID、供应商和 SMILES |
| [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) | 对 80 个官方主名称候选查询 CID、SMILES、InChIKey，HTTP 成功；所有返回结果、查询 URL、抓取 UTC 时间已缓存 |

CTRP 原包由另一数据整理模块按 HTTP Range 提取所需成员，其 `rebuild_v2/data_raw/ctrp/download_manifest.json` 记录来源、范围、CRC32 和 SHA256。本模块从该已验证文件复制到 `data_raw/drug_metadata`。该 compound 文件大小为 116,632 字节，原包 CRC32=1614100610。

`data_raw/drug_metadata/drug_source_manifest.json` 包含本模块下载/缓存文件的字节数和 SHA256。CTRP 下载地址通过药敏数据整理维护者的公开配置定位，已检验原包结构及 CRC；不要把该地址发现过程误写成新版本 CTRP 数据发布。

本地另有旧 `approved_druglist.txt`、`info_drugbank.csv`、`full database.xml`、`data.json`。已盘点其文件和 CSV 结构：导出 CSV 以 DrugBank ID/药名/靶点为主，没有 GDSC–CTRP 标识连接。本轮不使用其“approved”标签推断当前监管状态，也不使用未经核验来源/版本的旧别名覆盖官方映射；1.35 GB XML 和 1.16 GB JSON 未全文读取。

## 判定规则

1. 首先以大小写和连续空白规范化后的官方主名称/明确同义名生成候选；不做字符串模糊匹配，不去掉盐型、立体前缀或组合药比例。GDSC 同义名不单独作为身份成立依据。
2. 用 GDSC 官方主名称查询 PubChem；保留所有返回结果。多 CID 不任意选第一条。
3. 用 RDKit **2026.09.1** 从 CTRP 官方 SMILES 生成完整 Standard InChIKey。比较完整 key，不能只比较第一段连接性。
4. 只自动纳入完整 key 一致、单组分、无未指派潜在立体元素的单药。保留结构原文；不去盐、不删碎片、不抹除手性。
5. 49 个带明确比例的组合药全体排除单药主任务；总共有 50 个多组分 SMILES，不能将所有点号一律解释为盐并去掉其余部分。
6. 同一个 GDSC ID 若出现多个同等合格 CTRP ID，全部标记 ambiguous；同药多个 GDSC ID 仍各留一行，不合并响应。
7. 未通过者留在 `drug_identity_ambiguous.csv`，不能删除证据或把连接性一致升级成无条件合格。不存在名称候选者标记 `unmatched_no_candidate`，表示本轮策略未建立证据，不表示证明两个筛选绝无该药。

候选生成仍受原始名称/官方别名覆盖所限；本轮不是对 295×545 所有化学结构的穷尽标准化，也没有 GDSC 原始供应商批次结构可用于彻底解决冲突。扩展候选应单独登记新证据，冻结后不得根据预测好坏回改身份。

## 审计结果

| 候选关系状态 | 行数 | 下游处理 |
|---|---:|---|
| structure_supported | 64 | 可进入样本/浓度/独立性资格审查 |
| ambiguous_stereochemistry_or_protonation | 8 | 暂不进入主任务 |
| ambiguous_unassigned_stereochemistry | 5 | 暂不进入主任务 |
| ambiguous_salt_or_multicomponent | 4 | 暂不进入主任务 |
| ambiguous_structure_conflict | 3 | 暂不进入主任务 |
| ambiguous_multiple_pubchem_cids | 1 | 暂不进入主任务 |
| unmatched_no_candidate | 211 | 本轮未建立映射 |

总计 296 条关系覆盖全部 295 个 GDSC ID，多一行来自 Selumetinib 1736 的一个错误官方别名候选，保留冲突证据而非静默删除。

旧 CTRP 的 545 个名称全部能唯一回连到官方 master_cpd_id，详见 `drug_legacy_ctrp_name_to_id.csv`。这仅恢复身份，不能证明旧 AUC 数值、聚合方法或样本独立性正确；后者由 CTRP 原始响应整理模块另行审计。

### 必须保留的具体风险

- **Docetaxel**：GDSC 1007、1819；CTRP 单药 660364。PubChem CID148124 的完整 key 是 `ZDZOTLJHXYCWBA-VCVYQWHSSA-N`，CTRP SMILES 得到 `ZDZOTLJHXYCWBA-JXHJHTFNSA-N`。连接块相同而立体块不同，保留冲突，当前不应以旧稿常用药为由强行纳入。CTRP 660217 是 docetaxel:tanespimycin (2:1) 组合，绝不能与 660364 混同。
- **Selumetinib**：GDSC 1062、1736 都能结构支持回连 CTRP 348991，但官方 1736 的同义名中还包含 AZD1480。该错误候选回连 CTRP 660306 后结构不符，被标记 `ambiguous_structure_conflict`。靶点或别名“看起来接近”不能补救。
- **Bosutinib、Serdemetan**：同名但当前官方来源结构的连接性不同，不能自动接受。
- **Paclitaxel、Tretinoin、SN-38 等**：出现立体信息不一致或缺失。未按“大家都知道是同一种药”手动覆盖。
- **Oxaliplatin、Sepantronium bromide、BMS-345541、Obatoclax Mesylate**：涉及多组分/盐型、络合结构表示或多个 CID，不使用未经验证的去盐规则。

已明确核验的常见示例：Lapatinib 1558→634309；AZD7762 1022→660777；Methotrexate 1008→30371；Gefitinib 1010→52926；Vorinostat 1012→56554；Olaparib 1017→411867。这些仅为身份示例，不构成最终主药预选。

## 多 DRUG_ID 与浓度范围

[Sanger 官方文档](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/) 明确指出，GDSC2 同一化合物多个 DRUG_ID 与内部追踪有关。当前本地数据有 9 组同名多 ID：Docetaxel、Selumetinib、Oxaliplatin、Fulvestrant、Uprosertib、GSK343、Acetalax、Dactinomycin、Ulixertinib。`gdsc_same_name_ids` 和 `gdsc_duplicate_name` 保留该关系。

有 458 个唯一 DRUG_ID–MIN_CONC–MAX_CONC 组合；219 个 ID 只有一个范围，另 76 个存在 2–7 个范围。例如 Navitoclax 1011 有 6 个范围，Vorinostat 1012 有 3 个范围。药物身份通过后仍应审计筛选浓度；不同浓度区间上的 AUC 不天然相同。`drug_gdsc2_assay_inventory.csv` 完整保留范围及其记录数，以供覆盖规则或敏感性分析使用。最低浓度的近似差别不应先凭肉眼归为等价；主流程需明确容差/范围规则。

63 个合格 CTRP 单药中 Selumetinib 对应两个 GDSC ID。主药数应按独立化合物计算；建议一个月主分析优先选单一 GDSC ID 的药物，若选择 Selumetinib，则先固定一个 ID 的技术理由并将另一个作为独立重现/敏感性分支，不可按效能择优。

## 文件、运行与验证

- `src/curate_drugs.py`：元数据盘点、PubChem 查询缓存、结构映射和验证。
- `data_processed/drug_gdsc2_assay_inventory.csv`：GDSC 标识、靶点、公司/版本、浓度及行数。
- `data_processed/drug_ctrp_compound_identity.csv`：545 个 CTRP ID、结构、组合标志和未指派立体元素数。
- `data_processed/drug_identity_map.csv` / `drug_identity_eligible.csv` / `drug_identity_ambiguous.csv`：全部/合格/待审关系。
- `data_processed/drug_legacy_ctrp_name_inventory.csv` / `drug_legacy_ctrp_name_to_id.csv`：旧名回连。
- `data_processed/drug_identity_summary.json`：机器可读计数和软件版本。

使用项目环境依次运行：

```powershell
& 'rebuild_v2\.venv\Scripts\python.exe' 'rebuild_v2\src\curate_drugs.py' --inventory-only
& 'rebuild_v2\.venv\Scripts\python.exe' 'rebuild_v2\src\curate_drugs.py' --pubchem-candidates
& 'rebuild_v2\.venv\Scripts\python.exe' 'rebuild_v2\src\curate_drugs.py' --build-map
& 'rebuild_v2\.venv\Scripts\python.exe' 'rebuild_v2\src\curate_drugs.py' --self-check
```

前提是两个已下载并校验的官方元数据文件位于 `data_raw/drug_metadata`。PubChem 默认复用缓存，避免重跑时静默换版本；若有意更新，应另立新快照，不在冻结分析中覆盖。

已运行验证通过：全部 GDSC ID 保留、无重复映射边、每个合格 GDSC ID 仅对应一个 CTRP ID、完整 InChIKey 一致、无组合/多组分/未指派立体元素进入合格集，以及 Lapatinib 正确 ID、Docetaxel 冲突不纳入、错误 AZD1480 别名不纳入三项具体回归检查。
