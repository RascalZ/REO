"""Publish aggregate pilot evidence; individual predictions/features remain local."""
from pathlib import Path
import argparse,csv,html,json,shutil,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
LABELS={'mean':'训练均值','tissue':'组织Ridge','expression':'表达EN（2000基因）','expression_matched':'表达EN（80基因）','reo':'REO EN（80基因）'}
DRUGS={1931:'MIRA-1',1036:'PLX-4720'}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run-dir',type=Path,default=ROOT/'runs/gdsc2-mira1-plx4720-development-nested-en-reo80');args=ap.parse_args()
    run=args.run_dir;meta=json.loads((run/'run_manifest.json').read_text(encoding='utf-8'))
    assert meta['status']=='completed' and not meta['external_response_accessed'] and not meta['blocked_external_access_attempts']
    dest=ROOT/'results/pilot_development';dest.mkdir(parents=True,exist_ok=True)
    for source,target in [('summary_metrics.csv','summary_metrics.csv'),('fold_metrics.csv','fold_metrics.csv'),('chosen_parameters.json','chosen_parameters.json'),('feature200_benchmarks.json','feature200_benchmarks.json'),('endpoint_audit.json','endpoint_audit.json')]:shutil.copyfile(run/source,dest/target)
    metrics=pd.read_csv(run/'summary_metrics.csv');folds=pd.read_csv(run/'fold_metrics.csv')
    choices=json.loads((run/'chosen_parameters.json').read_text(encoding='utf-8'))
    benchmarks=json.loads((run/'feature200_benchmarks.json').read_text(encoding='utf-8'))
    validation=json.loads((run/'pilot_validation.json').read_text(encoding='utf-8')) if (run/'pilot_validation.json').exists() else {'status':'pending'}
    resources=meta['resources'];rss=resources['sampled_peak_rss_bytes']/1024**3;peak=resources['windows_peak_working_set_bytes']/1024**3
    public={'pilot_version':'v0.1','scientific_status':'engineering pilot; exploratory development OOF only','cohort_version':'v1',
            'code_execution_commit':'629b8f5','successful_execution_code_hashes':meta['source_hashes'],'config_sha256':meta['config_sha256'],
            'cohort_freeze_sha256':meta['cohort_lock']['freeze_sha256'],'n_fits':meta['n_fits'],'n_oof_predictions':meta['n_oof_predictions'],
            'resources':resources,'external_response_accessed':False,'external_response_access_attempts':0,
            'split_warnings':meta['split_warnings'],'runtime_threads':'single process inside threadpool_limits(1); saved threadpool_info was obtained after context exit and may show restored defaults',
            'code_patch_after_successful_run':'source-row validation exception cleanup only; numerical model outputs are from recorded execution commit',
            'integrity_validation':validation,'final_confirmatory_results':False}
    (dest/'pilot_run_summary.json').write_text(json.dumps(public,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    blocks=[]
    def title(s):blocks.append(('h2',s))
    def para(s):blocks.append(('p',s))
    def table(headers,rows):blocks.append(('table',(headers,rows)))
    para(f'日期：2026-10-10。MIRA-1与PLX-4720开发专用试跑已完成：{meta["n_fits"]}次拟合、{meta["n_oof_predictions"]}条逐模型OOF预测。总耗时{resources["wall_seconds"]:.1f}秒，采样峰值RSS {rss:.3f} GiB，Windows报告进程峰值工作集{peak:.3f} GiB。CTRP响应未读取用于试跑，当前禁读路径未出现访问尝试。')
    para('本报告是工程与预算验收，不是正式论文结果。两药根据冻结覆盖排名选择，全部样本及成绩保留；不据此换药、筛掉低成绩或宣称REO已经优于基线。')
    title('试跑设计和防泄漏验收')
    para('626/623个开发模型；同一全局已知生物组划分，外层3折×1次、内层3折。每个候选与每个内层重新拟合检测率、MAD筛选、基因对过滤、标准化、组织编码和目标统计。输入仅按基因注释预定蛋白编码范围，未按全开发或外验表达分布预选。')
    para('组织项不惩罚：先在原TPM构造REO，再训练内做组织投影并恢复组织贡献。未知组织的组织效应为0；组织基线回退全局训练均值，分子模型保留其分子预测。合成独立联合优化复核与FWL实现一致。')
    para('内层以原GDSC AUC尺度的合并OOF RMSE选6个配置；外层标签最后才参与评价。24项独立合成/访问/源行/哈希/异常收尾验收通过；实际运行的fit成员审计状态：'+str(validation.get('status'))+'。训练成员、outer/inner边界及每层特征均有本地记录。')
    title('开发集OOF比较：探索性点估计')
    rows=[]
    for _,r in metrics.iterrows():rows.append([DRUGS[int(r.gdsc_drug_id)],LABELS[r.model],str(int(r.n)), 'NA' if pd.isna(r.spearman) else f'{r.spearman:.3f}',f'{r.rmse:.5f}',f'{r.r2:.3f}',f'{r.wall_seconds:.1f}'])
    table(['药物','模型','N','Spearman','RMSE','R²','模型流程耗时/秒'],rows)
    para('均值基线在每个外折内常数，拼接不同折均值可能制造OOF相关，因此总体Spearman明确记NA；误差和R²仍计算。逐折分数另存fold_metrics.csv，不能把三个相依折当独立重复作普通t检验。各模型用同一测试ID，未裁剪超出0–1的预测。')
    para('MIRA-1的REO80 RMSE高于表达2000与组织基线；PLX-4720的REO80排序点估计略高，RMSE与组织基线非常接近。差异未经重复CV/统计推断，当前没有足够证据宣称REO具有可靠增量价值。')
    title('运行时间、内存和正式预算边界')
    table(['项目','实测'],[['完整试跑（数据准备/拟合；不含解释器启动前时间）',f'{resources["wall_seconds"]:.1f}s'],['累计CPU时间',f'{resources["cpu_seconds"]:.1f}s'],['采样峰值RSS（完整单进程，不是独立模型内存）',f'{rss:.3f} GiB'],['Windows进程峰值工作集',f'{peak:.3f} GiB'],['真实拟合/OOF条数',f'{meta["n_fits"]}/{meta["n_oof_predictions"]}'],['收敛警告/回退次数',f'{int(metrics.convergence_warning_fits.sum())}/{int(metrics.fallback_fits.sum())}']])
    table(['药物','REO200保留对数','训练内构造秒数','int8训练特征大小/MiB'],[[DRUGS[b['gdsc_drug_id']],str(b['retained_pairs']),f'{b["fit_transform_wall_seconds"]:.3f}',f'{b["matrix_bytes"]/1024**2:.2f}'] for b in benchmarks])
    para('REO200只测构造，未拟合/调参，所以不能按这0.6秒估计完整200基因EN训练。BLAS运行时受单线程上下文约束；run_manifest的threadpool_info在退出后采集，可能显示恢复后的默认值。资源CSV为100ms完整进程采样。')
    title('发现的问题与处理')
    table(['发现','处理/下一步'],[
        ['全部18个EN外折选择最强alpha=0.1边界','正式分析前预定包含更强正则的路径，或训练子折内计算alpha_max；仍只用开发数据，不按外验选网格。'],
        ['冻结v1含1个开发/5个外验明确Non-Cancerous模型，1个外验疾病注释缺失','试跑保持原v1；已按metadata记录资格复核，正式主队列v2先设癌模型/非癌控制/未知状态规则，不看响应决定。'],
        ['源行校验失败时Excel生成器未关闭，Windows句柄占用','已补try/finally关闭迭代器与workbook，错误药物/曲线/剂量/AUC测试现通过；旧执行代码保留Git提交629b8f5。'],
        ['runner原先依赖手动核验锁','拟合前已加freeze与实际开发输入哈希核验，该加固包含在本次执行版本中。'],
        ['异常中断可能残留running状态，OOF需等全run结束','已补失败状态/资源收尾与逐外层模型OOF/参数检查点；正式长运行仍需验收自动续跑，本次462fits已正常完成。'],
        ['GDSC发布AUC上游借用其他细胞系曲线信息','明确严格训练内拟合仅针对下游预测，不声称上游端点估计已CV隔离。'],
        ['稀有组织不能每折都有，出现1个未知组织测试样本','记录split warnings及未知组织数；预测按固定回退规则，不复制稀有组跨折。'],
        ['80基因与仅线性试跑尚未完成正式比较','下一步冻结资格v2、扩大正则路径，再测REO200实际求解，加入预设RF和5×3重复评价。']])
    title('可追溯材料和下一步')
    para('唯一总计划：RESEARCH_PLAN.md。配置：config/pilot_v0.1.json。执行代码版本：629b8f5；运行源文件哈希和队列锁哈希存聚合pilot_run_summary.json。原v1队列、8药面板及CTRP封存状态均保留；癌模型资格修订必须另建版本。')
    para('本地runs保留：开发响应源行、全局fold manifest、每次fit样本ID/参数/警告、每外层基因/基因对清单、逐模型OOF预测、资源采样及实际运行验证。仓库同步聚合指标、参数、设计审查、测试和本报告，不上传大型矩阵或个人原文件清单。')
    md=['# 阶段2：开发集嵌套试跑报告\n'];body=[]
    for kind,value in blocks:
        if kind=='h2':md.append('## '+value+'\n');body.append('<h2>'+html.escape(value)+'</h2>')
        elif kind=='p':md.append(value+'\n');body.append('<p>'+html.escape(value)+'</p>')
        else:
            headers,rows=value;md.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(row)+' |' for row in rows)+'\n')
            body.append('<div class="table"><table><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>')
    target=ROOT/'reports/阶段2_开发集嵌套试跑.md';target.write_text('\n'.join(md),encoding='utf-8')
    style='body{background:#f2f5f7;color:#193749;font:16px/1.9 "Microsoft YaHei",sans-serif;margin:0}main{max-width:1200px;margin:30px auto;padding:32px 45px;background:white;border-radius:10px}h1{font-size:28px}h2{font-size:22px;border-bottom:2px solid #d7e6ed;margin-top:36px;padding-bottom:8px}.table{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:14px;margin:16px 0}th,td{padding:10px 12px;text-align:left;border:1px solid #d8e4ea;vertical-align:top}th{background:#e7f0f5}tr:nth-child(even) td{background:#f7fafc}@media(max-width:650px){main{padding:16px;margin:10px}}'
    target.with_suffix('.html').write_text('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>阶段2：开发集嵌套试跑</title><style>'+style+'</style></head><body><main><h1>阶段2：开发集嵌套试跑</h1>'+''.join(body)+'</main></body></html>',encoding='utf-8')
    print(json.dumps({'report':str(target),'fits':meta['n_fits'],'peak_rss_GiB':rss,'validation':validation.get('status')},ensure_ascii=False))

if __name__=='__main__':main()
