# REO · 药物反应预测研究重构

以可追溯的样本内相对表达顺序（REO）研究癌细胞系药物反应预测。公平比较REO、普通表达模型和组织基线，评估未见细胞系上的跨药敏筛选迁移。首选投稿期刊：**BMC Genomics**。

**唯一持续维护的总计划：[RESEARCH_PLAN.md](RESEARCH_PLAN.md)。** 每次执行前先看§0；错误、方案变更及依据写入变更表。旧计划和旧稿作为历史快照保留。

历史阶段1冻结8种化合物、627个开发和180个主外验模型/生物组；后续癌队列v2修订为626/175。阶段2完成MIRA-1/PLX-4720的开发专用嵌套试跑（REO80、表达EN和组织基线）：462次拟合、约135秒、采样峰值RSS约0.66 GiB。历史pilot与当前正式结果均保留。

历史问题见[阶段2试跑报告](reports/阶段2_开发集嵌套试跑.md)。资格复核发现少量非癌对照；v1/pilot保留原样，正式队列修订见总计划与metadata审查记录。

**固定8药正式评价完成。** [阶段4报告](reports/阶段4_固定八药重复嵌套评价.md)与[科学评估](docs/FORMAL_SCIENTIFIC_ASSESSMENT.md)记录840任务、18,120fits、104,307 OOF及10,000次配对bootstrap。主EN宏ΔSpearman −0.028458，RF −0.036612，当前REO流程未显示强表达对照的优势；85/120主REO外层为零分子解。原673检查点SHA保持，Windows控制器失败记录保留，数值源码和正式lock未改。

**当前：四药表示诊断完成，总计划更新至v2.7。** [已有OOF分解](docs/OOF_DECOMPOSITION.md)保留旧结果；[冻结协议](docs/ORDINAL_DIAGNOSTIC_V1_PROTOCOL.md)的两种候选规则×表达/样本内秩/完整pair三种表示、同政策组织参照及9次全流程工程置换已全部执行。独立核验通过：840个成功外层任务、22,440次设计拟合、173,600条OOF记录，其中真实标签为84任务、2,244次拟合和17,360条OOF；其余属于置换控制。详见[阶段5报告](reports/阶段5_秩与基因对诊断试跑.md)与[科学评估及下一步](docs/ORDINAL_DIAGNOSTIC_ASSESSMENT.md)。

在这四个诊断案例中，组织内响应关联选基因相对MAD的宏Δρ分别为表达+0.020439、秩+0.032689、pair+0.041295，三者的探索性配对95%区间均高于0。但响应关联规则下pair相对表达的宏Δρ仅+0.002021（95% CI −0.010449至+0.014278），尚不支持REO具有可靠预测优势，也不证明等效。真实标签72个分子外层拟合中26个为零分子解；Tamoxifen在响应关联规则下三种表示均退回组织参照。四药由先前结果启发选择，区间以已有OOF为条件且未校正多重比较；9次置换是工程控制，不能用于常规0.05显著性结论。旧八药正式结果保持有效，CTRP继续封存。

诊断科学源固定于`e9b576d`。两次内存中断的失败记录、601个唯一成功检查点SHA及原版本恢复证据保留，见[执行记录](docs/ORDINAL_MEMORY_RECOVERY.md)。包含失败尝试的31次完整拟合事件及1次不预测/不评分的资源探针，总计22,472次完整拟合事件，另有2次中断拟合；控制器各会话合计约2小时59分，首启至完成约3小时10分。这些是本机实际执行资源记录，不是隔离的算法速度基准。

## 快速查看

```sh
git clone https://github.com/RascalZ/REO.git
cd REO
python scripts/verify_public_snapshot.py
python scripts/verify_public_ordinal.py
```

这些入口仅用Python标准库核验公开汇总的哈希、结构及计数，不重新拟合模型或核验个体预测。完整身份与拟合审计需要版本化原始输入，见[数据与许可](docs/DATA_AVAILABILITY.md)、[复现命令及边界](docs/REPRODUCIBILITY.md)。

## 当前材料

- [阶段5报告](reports/阶段5_秩与基因对诊断试跑.md)、[四药诊断科学评估](docs/ORDINAL_DIAGNOSTIC_ASSESSMENT.md)、[公开聚合结果](results/ordinal_diagnostic/)。
- [阶段1报告](reports/阶段1_数据核验与队列锁定.md)：药物ID/样本量/浓度、排除记录、外验限制。
- [主面板](results/main_drugs.csv)、[295药物样本量](results/all_295_GDSC_drug_sample_counts.csv)、[结构支持候选](results/all_structure_supported_drug_coverage.csv)。
- [下载登记](data_manifest/download_register.csv)、[旧代码处置](reports/现有代码处置清单.md)。
- [表达审计](reports/expression_audit.md)、[化合物审计](reports/drug_identity_audit.md)、[CTRP原包审计](reports/ctrp_provenance_audit.md)。
- 本地阅读版：`reports/阶段1_数据核验与队列锁定.html`、`reports/阶段1_队列与下载清单.xlsx`。

## 结构

```text
RESEARCH_PLAN.md            # 唯一维护的总计划与变更台账
src/, scripts/             # 整理/验证实现和顺序执行入口
config/                    # 版本化准入配置；本机路径和完整freeze留本地
docs/                      # 数据获取/许可、复现状态与命令
data_manifest/             # 下载登记；个人路径/原文件盘点不提交
data_raw/, data_processed/ # 本地生成和缓存，Git忽略
results/                   # 可公开聚合计数；逐样本manifest留本地
reports/                   # 科学审计；二进制交付和完整验证留本地
.github/workflows/         # 无大型输入的语法/公开汇总检查
```

## 安装与重放

阶段1使用Python 3.12.14及固定依赖。将`config/input_paths.example.json`复制为`input_paths.local.json`并填写原文件路径，或设置`REO_ORIGINAL_ROOT`；原文件只读。

```sh
python -m venv .venv
# 按操作系统激活环境
python -m pip install -r requirements.txt
python scripts/run_curation.py --download --lock
```

已有锁时，该入口只执行验证。该锁冻结队列及数据整理实现；正式八药与四药诊断各有独立的模型/统计锁及已完成运行。干净环境全流程源数据重放和模型复现仍为发布待验项；CI通过不等于预测有效。

开发试跑依赖另列，保留冻结的阶段1依赖文件：

```sh
python -m pip install --index-url https://pypi.org/simple -r requirements-modeling.txt
python -m unittest discover -s tests -p test_pilot_integrity.py -v
python src/run_development_pilot.py --run-dir runs/my-development-pilot
```

不要覆盖已有run。原执行代码提交为`629b8f5`，后续异常资源释放/检查点修复不改变已记录的数值结果；具体执行源哈希随报告发布。

队列v2/REO200/RF入口（需要本地已冻结v1/v2及原输入）：

```sh
python scripts/validate_cancer_cohort.py
python -m unittest discover -s tests -p "test_pilot*_integrity.py" -v
python src/run_development_pilot.py --config config/pilot_v0.2.json --run-dir runs/my-cancer-v2-pilot
python scripts/validate_pilot_run.py --run-dir runs/my-cancer-v2-pilot
python src/report_reo200_rf_pilot.py --run-dir runs/my-cancer-v2-pilot
```

阶段3成功运行提交为`5d05349`；随后修复的纯组织响应舍入边界未在真实906个fit中触发，实际运行保留，聚合和审计在`results/pilot_reo200_rf/`。完整freeze/个体manifest与预测留本地；公开癌资格来源事实缓存位于`data_raw/cancer_metadata/`并保留来源许可。

启动器自动记录Git HEAD，并在拟合前校验科学代码/配置与该提交一致。先提交新的科学配置或实现再启动；未提交改动会在拟合前被拒绝，失败证据保留。文档/报告的未提交改动不影响科学源快照。

CTRP曾用于历史项目；当前是回顾性重新隔离的已知身份验证。两库共享CCLE表达且剂量区间不同，结论限于跨筛选、跨细胞系的排序迁移，未证明跨表达平台、相同绝对AUC终点或患者效用。

## 引用与许可

软件引用见[CITATION.cff](CITATION.cff)。暂未发表的新稿不提供论文DOI。新写重构代码和说明采用[MIT](LICENSE)；上游数据、旧代码和原始稿件保留各自权利和使用条款。
