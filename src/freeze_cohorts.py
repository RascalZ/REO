"""Choose drug panel and external identities using coverage only; never load response values."""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import hashlib,json,math
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; P=ROOT/'data_processed'; R=ROOT/'results'
def read(name): return pd.read_csv(P/name,dtype=str,keep_default_na=False)
def csv(df,name): df.to_csv(R/name,index=False,encoding='utf-8-sig')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()
def eligible_bool(s): return s.astype(str).str.lower().isin(['true','1'])
def count_models(df): return df.canonical_model_id.nunique()
def count_groups(df): return df.biological_group_id.nunique()

def main():
    lock_path=ROOT/'config/cohort_freeze_v1.lock.json'
    freeze_path=ROOT/'config/cohort_freeze_v1.json'
    if lock_path.exists():
        lock=json.loads(lock_path.read_text(encoding='utf-8'))
        assert sha(freeze_path)==lock['freeze_sha256'], 'Locked freeze manifest changed; create a new cohort version.'
        frozen=json.loads(freeze_path.read_text(encoding='utf-8'))
        for item in frozen['files']:
            p=Path(item['path']) if item.get('scope')=='original_read_only' else ROOT/item['path']
            assert p.is_file() and sha(p)==item['sha256'], f'Locked dependency changed: {p}. Create a new cohort version.'
        print('Cohort v1 is locked; all frozen dependency hashes verified. No files were rewritten.')
        return
    rules=json.loads((ROOT/'config/cohort_rules.json').read_text(encoding='utf-8'))
    expression=read('cell_line_dictionary.csv')
    gm=read('gdsc2_cell_mapping.csv')
    graw=read('gdsc2_observation_index.csv')
    cm=read('ctrp_cell_mapping.csv')
    craw=read('ctrp_curve_index.csv')
    drugs=read('drug_identity_eligible.csv')
    drugs=drugs[eligible_bool(drugs.identity_eligible)].copy()
    assert not drugs.gdsc_drug_id.duplicated().any(), 'More than one accepted compound for GDSC ID'
    modelgroups=read('canonical_model_groups.csv').set_index('canonical_model_id').biological_group_id.to_dict()
    # Exclude all known GDSC models AND related groups, including identity-conflict candidates,
    # from natural external validation. Availability in selected drugs is not used here.
    all_gdsc_ids={a for v in gm.candidate_model_ids for a in v.split('|') if a}
    all_gdsc_groups={modelgroups[a] for a in all_gdsc_ids if a in modelgroups}
    g=graw.merge(gm,on=['sanger_model_id','cosmic_id','cell_line_name','tcga_desc'],how='left',validate='many_to_one')
    g['base_eligible']=g.expression_eligible.eq('1') & g.response_available.eq('1') & g.webrelease.eq('Y') & g.dataset.eq('GDSC2')
    g=g[g.base_eligible].copy()
    assert not g.duplicated(['gdsc_drug_id','canonical_model_id']).any()
    # Dominant exact dose interval avoids mixing different AUC integration ranges within drug.
    dose=g.groupby(['gdsc_drug_id','min_conc','max_conc']).agg(n_models=('canonical_model_id','nunique'),n_groups=('biological_group_id','nunique')).reset_index()
    dose['lo_numeric']=pd.to_numeric(dose.min_conc); dose['hi_numeric']=pd.to_numeric(dose.max_conc)
    dose=dose.sort_values(['gdsc_drug_id','n_groups','n_models','lo_numeric','hi_numeric'],ascending=[True,False,False,True,True])
    dose['primary_dose_interval']=~dose.gdsc_drug_id.duplicated()
    csv(dose,'gdsc_dose_range_decisions.csv')
    primarydose=dose[dose.primary_dose_interval][['gdsc_drug_id','min_conc','max_conc']]
    gprimary=g.merge(primarydose,on=['gdsc_drug_id','min_conc','max_conc'],validate='many_to_one')
    all_gdsc=graw.groupby(['gdsc_drug_id','gdsc_drug_name']).agg(n_raw_models=('sanger_model_id','nunique'),n_raw_curve_records=('source_row_number','size')).reset_index()
    beforedose=g.groupby('gdsc_drug_id').agg(n_identity_expression_models=('canonical_model_id','nunique'),n_identity_expression_groups=('biological_group_id','nunique')).reset_index()
    afterdose=gprimary.groupby('gdsc_drug_id').agg(n_primary_dose_models=('canonical_model_id','nunique'),n_primary_dose_groups=('biological_group_id','nunique')).reset_index()
    identities=read('drug_identity_map.csv').groupby('gdsc_drug_id').agg(identity_statuses=('match_status',lambda x:'|'.join(sorted(set(x))))).reset_index()
    all_gdsc=all_gdsc.merge(beforedose,on='gdsc_drug_id',how='left').merge(afterdose,on='gdsc_drug_id',how='left').merge(identities,on='gdsc_drug_id',how='left')
    csv(all_gdsc,'all_295_GDSC_drug_sample_counts.csv')
    c=craw.merge(cm[['master_ccl_id','canonical_model_id','biological_group_id','ccle_id','tissue','expression_row_index','expression_eligible','mapping_status']],on='master_ccl_id',how='left',validate='many_to_one')
    c['base_eligible']=c.expression_eligible.eq('1') & c.response_available.eq('1')
    c=c[c.base_eligible].copy()
    c['natural_external']=~c.biological_group_id.isin(all_gdsc_groups)
    c['identity_qc_primary']=c.snp_fp_status.eq('SNP-matched-reference')
    cprimary=c[c.identity_qc_primary].copy()

    def coverage(holdout_groups):
        rows=[]
        for _,d in drugs.iterrows():
            gd=gprimary[gprimary.gdsc_drug_id==d.gdsc_drug_id]
            development=gd[~gd.biological_group_id.isin(holdout_groups)]
            ct=cprimary[cprimary.master_cpd_id==d.ctrp_compound_id]
            external=ct[ct.natural_external | ct.biological_group_id.isin(holdout_groups)]
            untreated_qc=c[(c.master_cpd_id==d.ctrp_compound_id) & (c.natural_external | c.biological_group_id.isin(holdout_groups))]
            intervals=primarydose[primarydose.gdsc_drug_id==d.gdsc_drug_id]
            lo=intervals.iloc[0].min_conc if len(intervals) else ''
            hi=intervals.iloc[0].max_conc if len(intervals) else ''
            ng=count_groups(development); ne=count_groups(external); nt=development.tissue.nunique()
            reasons=[]
            if ng<rules['min_development_biological_groups']: reasons.append('development_groups_below_200')
            if ne<rules['min_external_biological_groups']: reasons.append('external_groups_below_50')
            if nt<rules['min_development_tissues']: reasons.append('development_tissues_below_4')
            rows.append({**d.to_dict(),'gdsc_min_conc':lo,'gdsc_max_conc':hi,
                         'n_development_models':count_models(development),'n_development_groups':ng,
                         'n_development_tissues':nt,'n_external_primary_models':count_models(external),
                         'n_external_primary_groups':ne,'n_external_primary_tissues':external.tissue.nunique(),
                         'n_external_primary_curve_records':len(external),
                         'n_external_all_identity_status_models':count_models(untreated_qc),
                         'coverage_eligible':not reasons,'coverage_exclusion_reason':';'.join(reasons),
                         'n_natural_external_models':count_models(external[external.natural_external]),
                         'n_reserved_shared_models':count_models(external[~external.natural_external]),
                         'n_gdsc_expression_models_all_dose_ranges':count_models(g[g.gdsc_drug_id==d.gdsc_drug_id])})
        panel=pd.DataFrame(rows)
        panel['drug_id_numeric']=pd.to_numeric(panel.gdsc_drug_id)
        panel=panel.sort_values(['coverage_eligible','n_development_groups','n_external_primary_groups','n_development_tissues','drug_id_numeric'],ascending=[False,False,False,False,True])
        panel['primary_assay_for_chemical']=~panel.ctrp_inchikey.duplicated()
        return panel

    holdout=set(); cov=coverage(holdout)
    sufficient=int((cov.coverage_eligible & cov.primary_assay_for_chemical).sum())>=rules['primary_drug_count']
    design='natural_CTRP_only_known_biological_groups'
    if not sufficient:
        # Fixed deterministic tissue-stratified group reserve; used only if natural design cannot support eight drugs.
        shared_groups=set(g.biological_group_id)&set(cprimary.biological_group_id)
        by_tissue=defaultdict(list)
        for group in sorted(shared_groups):
            tissues=sorted(set(g.loc[g.biological_group_id==group,'tissue']))
            by_tissue['|'.join(tissues)].append(group)
        for tissue,groups in by_tissue.items():
            groups.sort(key=lambda x:hashlib.sha256(f"{rules['split_seed']}:{x}".encode()).hexdigest())
            n=max(1,int(round(len(groups)*rules['fallback_holdout_fraction']))) if len(groups)>=5 else 0
            holdout.update(groups[:n])
        cov=coverage(holdout); design='natural_CTRP_only_plus_global_group_holdout'
    candidates=cov[cov.coverage_eligible & cov.primary_assay_for_chemical].copy()
    assert len(candidates)>=rules['primary_drug_count'], 'Insufficient eligible drugs: do not freeze an underpowered panel'
    preferred=candidates[candidates.n_external_primary_groups>=rules['preferred_external_biological_groups']]
    precision_preference_applied=len(preferred)>=rules['primary_drug_count']
    if precision_preference_applied: candidates=preferred
    selected=candidates.head(rules['primary_drug_count']).copy()
    selected.insert(0,'primary_drug_rank',range(1,len(selected)+1))
    selected_ids=set(selected.gdsc_drug_id)
    cov['selected_primary']=cov.gdsc_drug_id.isin(selected_ids)
    cov['selection_decision']=cov.apply(lambda x:'selected_primary' if x.selected_primary else 'excluded_coverage:'+x.coverage_exclusion_reason if not x.coverage_eligible else 'alternative_assay_same_chemical' if not x.primary_assay_for_chemical else 'eligible_minimum_but_below_preferred_external_80' if precision_preference_applied and x.n_external_primary_groups<80 else 'eligible_not_selected_by_prespecified_coverage_rank',axis=1)
    csv(cov,'all_structure_supported_drug_coverage.csv')
    csv(selected,'main_drugs.csv')
    csv(cov[['gdsc_drug_id','gdsc_drug_name','ctrp_compound_id','selection_decision']],'drug_selection_decisions.csv')
    develop=[]; external=[]; sensitivity=[]
    for _,d in selected.iterrows():
        gd=gprimary[(gprimary.gdsc_drug_id==d.gdsc_drug_id)&(~gprimary.biological_group_id.isin(holdout))]
        for _,r in gd.iterrows():
            develop.append({'gdsc_drug_id':d.gdsc_drug_id,'drug_name':d.gdsc_drug_name,
                            'canonical_model_id':r.canonical_model_id,'biological_group_id':r.biological_group_id,
                            'ccle_id':r.ccle_id,'expression_row_index':r.expression_row_index,'tissue':r.tissue,
                            'gdsc_source_row_number':r.source_row_number,'nlme_curve_id':r.nlme_curve_id,
                            'min_conc':r.min_conc,'max_conc':r.max_conc,'role':'development'})
        cs=c[(c.master_cpd_id==d.ctrp_compound_id)&(c.natural_external|c.biological_group_id.isin(holdout))]
        for which,target in [(cs[cs.identity_qc_primary],external),(cs,sensitivity)]:
            for model,rows in which.groupby('canonical_model_id'):
                r=rows.iloc[0]
                target.append({'gdsc_drug_id':d.gdsc_drug_id,'drug_name':d.gdsc_drug_name,
                               'ctrp_compound_id':d.ctrp_compound_id,'canonical_model_id':model,
                               'biological_group_id':r.biological_group_id,'ccle_id':r.ccle_id,
                               'expression_row_index':r.expression_row_index,'tissue':r.tissue,
                               'ctrp_master_ccl_ids':'|'.join(sorted(set(rows.master_ccl_id))),
                               'ctrp_experiment_ids':'|'.join(sorted(set(rows.experiment_id),key=int)),
                               'ctrp_source_row_numbers':'|'.join(rows.source_row_number),
                               'n_source_experiments':rows.experiment_id.nunique(),
                               'identity_qc_statuses':'|'.join(sorted(set(rows.snp_fp_status))),
                               'external_subtype':'natural_CTRP_only' if bool(r.natural_external) else 'reserved_shared_group',
                               'future_response_aggregation':'median_of_included_experiments',
                               'role':'external_primary' if target is external else 'external_identity_qc_sensitivity'})
    dev=pd.DataFrame(develop); ext=pd.DataFrame(external); sens=pd.DataFrame(sensitivity)
    csv(dev,'development_sample_manifest.csv'); csv(ext,'external_primary_sample_manifest.csv'); csv(sens,'external_sensitivity_sample_manifest.csv')
    dids=set(dev.canonical_model_id); xids=set(ext.canonical_model_id)
    dgroups=set(dev.biological_group_id); xgroups=set(ext.biological_group_id)
    assert not (dids & xids), 'MODEL LEAKAGE'
    assert not (dgroups & xgroups), 'RELATED-MODEL LEAKAGE'
    assert not ext.duplicated(['gdsc_drug_id','canonical_model_id']).any()
    assert ext.identity_qc_statuses.eq('SNP-matched-reference').all()
    assert not selected.ctrp_inchikey.duplicated().any()
    # Role applies globally across every drug, including unused candidate drugs.
    split=expression.copy()
    ctrp_primary_model_pool=set(cprimary.canonical_model_id)
    ctrp_all_model_pool=set(c.canonical_model_id)
    def role(r):
        if r.expression_eligible!='1': return 'excluded_identity'
        if r.biological_group_id in holdout: return 'reserved_external_group_excluded_from_all_GDSC_development'
        if r.biological_group_id in all_gdsc_groups: return 'GDSC_development_pool_never_external'
        if r.canonical_model_id in ctrp_primary_model_pool: return 'natural_external_primary_pool'
        if r.canonical_model_id in ctrp_all_model_pool: return 'natural_external_sensitivity_only_pool'
        return 'unused_no_eligible_CTRP_or_GDSC_observation'
    split['global_role']=split.apply(role,axis=1)
    split['in_selected_development']=split.canonical_model_id.isin(dids)
    split['in_selected_external_primary']=split.canonical_model_id.isin(xids)
    csv(split,'global_split_manifest.csv')
    qc={'original_gdsc_rows':len(graw),'original_gdsc_model_ids':graw.sanger_model_id.nunique(),
        'gdsc_after_identity_expression_rows':len(g),'gdsc_after_identity_expression_models':count_models(g),
        'gdsc_after_primary_dose_range_rows':len(gprimary),
        'ctrp_original_curve_rows':len(craw),'ctrp_measured_models':craw.master_ccl_id.nunique(),
        'ctrp_after_identity_expression_rows':len(c),'ctrp_after_identity_expression_models':count_models(c),
        'ctrp_after_SNP_matched_rows':len(cprimary),'ctrp_after_SNP_matched_models':count_models(cprimary),
        'natural_external_primary_pool_models':count_models(cprimary[cprimary.natural_external]),
        'natural_external_primary_pool_groups':count_groups(cprimary[cprimary.natural_external]),
        'primary_design':design,'fallback_reserved_groups':len(holdout),
        'identity_supported_assay_ids':len(drugs),'coverage_eligible_unique_chemicals':int((cov.coverage_eligible&cov.primary_assay_for_chemical).sum()),
        'preferred_external_80_pool_chemicals':len(preferred),'preferred_external_80_applied':precision_preference_applied,
        'selected_drugs':len(selected),'development_drug_model_observations':len(dev),
        'external_drug_model_observations':len(ext),'development_unique_models':len(dids),
        'external_unique_models':len(xids),'development_unique_groups':len(dgroups),'external_unique_groups':len(xgroups),
        'model_overlap':len(dids&xids),'biological_group_overlap':len(dgroups&xgroups),
        'response_magnitudes_used_for_selection':False,'models_fitted':0,
        'historical_external_use_disclosed':True,
        'independence_limit':'New analysis identity separation within known metadata; CTRP was used in historical work, not prospectively blind validation.'}
    (R/'cohort_summary.json').write_text(json.dumps(qc,ensure_ascii=False,indent=2),encoding='utf-8')
    for kind,df in [('development',dev),('external_primary',ext)]:
        dist=df.groupby(['gdsc_drug_id','drug_name','tissue']).agg(n_models=('canonical_model_id','nunique'),n_groups=('biological_group_id','nunique')).reset_index()
        csv(dist,kind+'_tissue_counts.csv')
    # Complete dependency snapshot; a separate validation step writes the lock after checks pass.
    # Later training outputs and reports are deliberately outside this cohort/input snapshot.
    files=[ROOT/'config/cohort_rules.json',ROOT/'config/curation_environment.json']
    files.extend(ROOT/name for name in ['requirements.txt','pyproject.toml','config/input_paths.example.json'])
    for directory in ['src','scripts','data_raw','data_processed','data_manifest','results']:
        files.extend(p for p in (ROOT/directory).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
    records=[{'path':str(p.relative_to(ROOT)),'scope':'rebuild_v2','bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(set(files))]
    original=pd.read_csv(ROOT/'data_manifest/original_source_sha256.csv',keep_default_na=False).to_dict('records')
    original+=json.loads((P/'expression_audit.json').read_text(encoding='utf-8'))['source_inputs']
    seen=set()
    for item in original:
        p=Path(item['path'])
        if str(p) in seen: continue
        seen.add(str(p))
        current=sha(p)
        assert current==item['sha256'], f'Original input changed: {p}'
        records.append({'path':str(p),'scope':'original_read_only','bytes':p.stat().st_size,'sha256':current})
    freeze={'freeze_type':'cohort_and_drug_panel_v1','created_utc':datetime.now(timezone.utc).isoformat(),
            'rules':rules,'design':design,'main_drugs':selected[['gdsc_drug_id','gdsc_drug_name','ctrp_compound_id','gdsc_min_conc','gdsc_max_conc']].to_dict('records'),
            'files':records,'status':'snapshot_pending_independent_validation_lock',
            'training_configuration_frozen':False,'model_hyperparameters':'to be fixed after development-only pilot',
            'response_magnitudes_accessed_for_selection':False,'historical_CTRP_use':'yes; disclose retrospective redesign'}
    (ROOT/'config/cohort_freeze_v1.json').write_text(json.dumps(freeze,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(qc,ensure_ascii=False))
    print(selected[['gdsc_drug_id','gdsc_drug_name','n_development_models','n_development_groups','n_external_primary_models','n_external_primary_groups']].to_string(index=False))

if __name__=='__main__': main()
