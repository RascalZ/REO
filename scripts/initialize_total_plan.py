"""One-time reviewed-plan import; refuse to overwrite the maintained successor."""
from pathlib import Path
import argparse,csv
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('reviewed_plan',type=Path);args=ap.parse_args()
    dest=ROOT/'RESEARCH_PLAN.md'
    if dest.exists():raise SystemExit('Total plan already exists; update explicitly with change history.')
    prior=args.reviewed_plan.read_text(encoding='utf-8')
    prior=prior.replace('**评估日期：2026年10月10日；本版补充Frontiers现行范围比较、SMOTE与scDEAL取舍。工作周期：10月10日至11月8日，11月9日至10日保留投稿检查缓冲。**','**计划v2.0，更新日期：2026-10-10。唯一维护版本为本文件。工作周期：10月10日至11月8日，11月9日至10日保留投稿检查缓冲。**')
    prior=prior.replace('方向定义已核实的连续药敏曲线AUC','需在模型阶段完成单位/方向验收的连续药敏曲线AUC')
    prior=prior.replace('先选**8个主药物** ；覆盖与运行时间允许时在冻结前增至最多12个。','本轮按覆盖规则确定**8种主分析化合物**；队列v1锁定后不按模型效果扩大或替换。')
    prior=prior.replace('建议启动目标为每药开发集≥200个独立细胞系、≥4个组织，冻结外验尽量≥50且争取≥80。','按独立生物组计数：开发≥200、组织≥4、外验最低≥50；足够8种均达到≥80时优先在该池排序。本轮58种达到80目标，故主面板全部≥80。')
    prior=prior.replace('完整证据见同目录`code_audit.md`、`code_audit_summary.json`、`manuscript_audit.md`。','历史证据在原计划归档；当前处置见[现有代码处置清单](reports/现有代码处置清单.md)。')
    status='''## 0 执行台账与当前版本（每次工作先读此节）

**首选期刊：BMC Genomics。当前：数据、身份、覆盖和队列构建完成；正式验证/锁定状态待本轮收尾写入。模型拟合为0，尚无新的预测性能。**

本文件是唯一持续维护的总计划。原计划是历史快照，阶段报告是证据；发现错误、版本/范围/方法变化，应先更新本节及变更表，再建立对应配置版本。不得因外验成绩不好更换主药、身份规则或测试样本。

|工作包|状态|依据/下一步|
|---|---|---|
|P0 数据/代码审计|完成|[阶段1报告](reports/阶段1_数据核验与队列锁定.md)、[代码处置](reports/现有代码处置清单.md)|
|P1 下载/表达|必要部分完成|本地TPM 1,019×57,820、float64；补CTRP原包成员、结构/身份/GENCODE；[下载登记](data_manifest/download_register.csv)|
|P2 细胞/药物/准入|完成|64个GDSC IDs对应63种结构支持CTRP单药；全部295 IDs有样本量台账；冲突排除|
|P3 队列v1/验证|待最终验证|627开发、180外验唯一模型/生物组；药物×模型4,975/1,222不能累加当独立N|
|P4 开发专用试跑|下一步|连续端点解析/方向、REO平局、强表达和组织基线、内层fit边界、时间/内存预算|
|P5 模型/统计冻结和正式运行|未开始|外层5折×3重复/内层3折；所有监督筛选仅在内层训练；生物学单位CI/置换|
|P6 冻结外验|未开始|最终模型/统计配置冻结后评估CTRP响应，完整报告8药|
|P7 图表/正文/投稿|未开始|结果决定结论；不承诺旧高准确率或接收|

### 本轮输入与面板

表达为2018 CCLE RSEM TPM；2024 Broad Model和2026 Sanger快照只补身份/同源关系，未混用表达版本。开发来自GDSC2 27Oct23；外验为原始CTRPv2 expanded post-QC、SNP核验通过、在全部GDSC2及已知关联组之外的模型。自然外验池183，8药覆盖180；未启用20%共有模型保留。

|化合物|GDSC DRUG_ID|CTRP master_cpd_id|开发独立组|主外验独立组|
|---|---:|---:|---:|---:|
'''
    with (ROOT/'results/main_drugs.csv').open(encoding='utf-8-sig',newline='') as f:
        for r in csv.DictReader(f):status+=f"|{r['gdsc_drug_name']}|{r['gdsc_drug_id']}|{r['ctrp_compound_id']}|{r['n_development_groups']}|{r['n_external_primary_groups']}|\n"
    status+='''
逐药浓度/组织覆盖见[main_drugs.csv](results/main_drugs.csv)。面板含研究探针，不是8种均获批药物；药理类别只描述，未用于改面板。

### 解释边界与未完成项

- 旧全标签DEG/Fisher/RFE及下游结果须重算。CTRP历史用过，本次为回顾性重新设计的身份隔离验证，不能称前瞻盲测。
- 两库均用CCLE表达，当前检验跨药敏screen及未知细胞系排序迁移，不证明跨表达平台或患者效用。
- 每药GDSC2固定浓度区间，按独立组/模型覆盖及浓度数值排序；仅解决库内混合。跨库积分窗口不同，Spearman也不能消除其改变排序的可能。暂不做跨库绝对RMSE/校准主张。
- 共同浓度窗重积分先查拟合参数可行性；若增加，先设端点版本再计算，不能事后挑更好端点。
- 独立范围只指全部本地GDSC2及已知关联组；未知亲缘仍有限制，不能扩称所有GDSC1/2。
- 原项目兼容输入仍需存在；公开数据独立重放和干净环境全流程安装是发布前待验项。[复现状态](docs/REPRODUCIBILITY.md)写清完成/未完成。

### 计划变更记录

|日期|变更|依据|解释影响|
|---|---|---|---|
|2026-10-10|首选BMC Genomics，排除全部12本已拒期刊|用户确认及范围核查|不再将已拒BMC Bioinformatics/Human Genomics列备选|
|2026-10-10|固定现有2018 TPM，无需重下载全表达|找到完整矩阵、数值/秩序核验|未宣称远端压缩包字节同一性|
|2026-10-10|纳入patient/parent图及SNP matched准入|权威元数据可得|按独立生物组计数，未知亲缘仍保留|
|2026-10-10|优先外验≥80，固定8药|58种达80；未看响应大小或效能|计算前覆盖/精度修订，不按成绩挑药|
|2026-10-10|自然CTRP-only外验，不预留GDSC2的20%|183自然组覆盖足够|开发可用N较多，全局角色隔离|
|2026-10-10|加入库内浓度窗筛选及跨库限制|76个GDSC IDs多窗，跨库上限不同|外验保留剂量设计差异，不以Spearman掩盖|
|2026-10-10|Git仓库与唯一总计划|用户指定RascalZ/REO|原稿/大数据留本地；代码/聚合队列版本化|

---

'''
    prior=prior.replace('## 1 核心判断与发表可行性',status+'## 1 核心判断与发表可行性',1)
    prior=prior.replace('### 4.3 要重新下载什么、保留什么\n','### 4.3 要重新下载什么、保留什么\n\n**执行更新：** 必要补充已下载，完整本地TPM已确认。最新状态以§0及下载登记为准；下表扩展并非本月全部必做。\n',1)
    prior=prior.replace('### 4.4 两种外部测试设计与选择顺序\n','### 4.4 两种外部测试设计与选择顺序\n\n**执行更新：** 队列v1采用自然CTRP-only已知生物组。下面共有细胞保留为未启用备用设计。\n',1)
    dest.write_text(prior,encoding='utf-8');print(dest)

if __name__=='__main__':main()
