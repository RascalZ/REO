# REO · 药物反应预测研究重构

以可追溯的样本内相对表达顺序（REO）研究癌细胞系药物反应预测。公平比较REO、普通表达模型和组织基线，评估未见细胞系上的跨药敏筛选迁移。首选投稿期刊：**BMC Genomics**。

**唯一持续维护的总计划：[RESEARCH_PLAN.md](RESEARCH_PLAN.md)。** 每次执行前先看§0；错误、方案变更及依据写入变更表。旧计划和旧稿作为历史快照保留。

阶段1已冻结8种化合物、627个开发和180个主外验模型/生物组。阶段2已完成MIRA-1/PLX-4720的开发专用嵌套试跑（REO80、表达EN和组织基线）：462次拟合、约135秒、采样峰值RSS约0.66 GiB。分数仅为探索性内部结果；正式REO200、RF、5×3重复CV和CTRP评估尚未执行。结果允许REO不优于强基线。

最新进度与问题见[阶段2试跑报告](reports/阶段2_开发集嵌套试跑.md)。资格复核发现少量非癌对照；当前v1/pilot保留原样，正式队列修订见总计划与metadata审查记录。

## 快速查看

```sh
git clone https://github.com/RascalZ/REO.git
cd REO
python scripts/verify_public_snapshot.py
```

该入口仅用Python标准库核验公开汇总。完整身份核验需要版本化原始输入，见[数据与许可](docs/DATA_AVAILABILITY.md)、[复现命令及边界](docs/REPRODUCIBILITY.md)。

## 当前材料

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

已有锁时，该入口只执行验证。冻结的是队列及数据整理实现；模型/统计配置需开发专用试跑后另行冻结。干净环境全流程源数据重放和模型复现仍为发布待验项；CI通过不等于预测有效。

开发试跑依赖另列，保留冻结的阶段1依赖文件：

```sh
python -m pip install --index-url https://pypi.org/simple -r requirements-modeling.txt
python -m unittest discover -s tests -p test_pilot_integrity.py -v
python src/run_development_pilot.py --run-dir runs/my-development-pilot
```

不要覆盖已有run。原执行代码提交为`629b8f5`，后续异常资源释放/检查点修复不改变已记录的数值结果；具体执行源哈希随报告发布。

CTRP曾用于历史项目；当前是回顾性重新隔离的已知身份验证。两库共享CCLE表达且剂量区间不同，结论限于跨筛选、跨细胞系的排序迁移，未证明跨表达平台、相同绝对AUC终点或患者效用。

## 引用与许可

软件引用见[CITATION.cff](CITATION.cff)。暂未发表的新稿不提供论文DOI。新写重构代码和说明采用[MIT](LICENSE)；上游数据、旧代码和原始稿件保留各自权利和使用条款。
