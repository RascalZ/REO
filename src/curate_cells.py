"""Conservative exact-ID reconciliation and biological grouping; no label access."""
from pathlib import Path
from collections import defaultdict
import csv, json, re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data_processed'
from project_paths import original_input

def read(p, sep=','):
    return pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)

def id_value(x):
    return re.sub(r'\.0$', '', str(x).strip())

def mapping(df, key, value):
    out=defaultdict(set)
    for k,v in zip(df[key],df[value]):
        if k and v: out[id_value(k)].add(v)
    return out

class UnionFind:
    def __init__(self): self.parent={}
    def find(self,x):
        if x not in self.parent: self.parent[x]=x
        if self.parent[x]!=x: self.parent[x]=self.find(self.parent[x])
        return self.parent[x]
    def union(self,a,b):
        x,y=self.find(a),self.find(b)
        if x!=y: self.parent[max(x,y)]=min(x,y)

def main():
    broad=read(ROOT/'data_raw/cell_metadata/DepMap_24Q4_v1_Model.csv')
    legacy=read(original_input('legacy_cell_annotations'),'\t')
    expr=read(OUT/'expression_samples.tsv','\t')
    sid=read(OUT/'sanger_identifiers.csv')
    sm=read(OUT/'sanger_models.csv')
    broad_ids=set(broad.ModelID)
    ccle_ach=defaultdict(set)
    evidence=defaultdict(set)
    for df,k,v,source in [(broad,'CCLEName','ModelID','Broad24Q4_CCLEName'),
                          (legacy,'CCLE_ID','depMapID','CCLE2018_annotation')]:
        for a,b in zip(df[k],df[v]):
            if a and re.fullmatch(r'ACH-\d+',b):
                ccle_ach[a].add(b); evidence[a].add(source)
    sid_ach=defaultdict(set)
    for _,r in broad.iterrows():
        if r.SangerModelID: sid_ach[r.SangerModelID].add(r.ModelID)
    for _,r in sid.iterrows():
        if r.identifier_source in ('BROAD_ID','DepMap_ID (Broad_ID)') and re.fullmatch(r'ACH-\d+',r.identifier):
            sid_ach[r.sanger_model_id].add(r.identifier)
    cosmic_ach=mapping(broad,'COSMICID','ModelID')
    sid_ccle=defaultdict(set)
    for _,r in sid.iterrows():
        if r.identifier_source=='CCLE_ID':
            sid_ccle[r.sanger_model_id].add(r.identifier)
            ccle_ach[r.identifier].update(sid_ach[r.sanger_model_id])
            if sid_ach[r.sanger_model_id]: evidence[r.identifier].add('Sanger_official_CCLE_to_ACH')
        elif r.identifier_source=='COSMIC_ID': cosmic_ach[id_value(r.identifier)].update(sid_ach[r.sanger_model_id])
    # IDs that lack a direct Broad ID can use an exact CCLE accession cross-reference.
    for s,cc in sid_ccle.items():
        for c in cc: sid_ach[s].update(ccle_ach[c])
    for _,r in sid[sid.identifier_source=='COSMIC_ID'].iterrows():
        cosmic_ach[id_value(r.identifier)].update(sid_ach[r.sanger_model_id])

    uf=UnionFind(); edges=[]
    def edge(a,b,kind,source_id):
        if a and b and a!=b:
            uf.union(a,b); edges.append({'model_a':a,'model_b':b,'relation':kind,'source_id':source_id})
    for a in set(broad.ModelID)|set(expr.depmap_id): uf.find(a)
    for patient,g in broad[broad.PatientID!=''].groupby('PatientID'):
        ids=sorted(set(g.ModelID))
        for a in ids[1:]: edge(ids[0],a,'Broad_PatientID_same_tissue_sample',patient)
    # Sanger patient covers samples from same donor; parent links also retained.
    for patient,g in sm[sm.sanger_patient_id!=''].groupby('sanger_patient_id'):
        ids=sorted({a for s in g.sanger_model_id for a in sid_ach[s]})
        for a in ids[1:]: edge(ids[0],a,'Sanger_same_patient',patient)
    for _,r in sm.iterrows():
        ids=sorted(sid_ach[r.sanger_model_id])
        for a in ids[1:]: edge(ids[0],a,'Sanger_model_multiple_Broad_IDs',r.sanger_model_id)
        if r.parent_sanger_model_id:
            for a in ids:
                for b in sid_ach[r.parent_sanger_model_id]: edge(a,b,'Sanger_parent_model',r.parent_sanger_model_id)
    bdict=broad.set_index('ModelID').to_dict('index')
    ldict=legacy.set_index('CCLE_ID').to_dict('index')
    expr_rows=[]
    for _,r in expr.iterrows():
        candidates=ccle_ach[r.ccle_id]
        accepted=len(candidates)==1
        a=next(iter(candidates)) if accepted else ''
        b=bdict.get(a,{})
        l=ldict.get(r.ccle_id,{})
        expr_rows.append({'expression_row_index':r.expression_row_index,'ccle_id':r.ccle_id,
                         'canonical_model_id':a,'candidate_model_ids':'|'.join(sorted(candidates)),
                         'mapping_status':'resolved_exact_accession' if accepted else 'ambiguous' if candidates else 'unmapped',
                         'mapping_evidence':'|'.join(sorted(evidence[r.ccle_id])),
                         'biological_group_id':('BIOGROUP_'+uf.find(a)) if a else '',
                         'broad_patient_id':b.get('PatientID',''),
                         'tissue':b.get('OncotreeLineage','') or l.get('Site_Primary',''),
                         'tissue_source':'Broad24Q4_OncotreeLineage' if b.get('OncotreeLineage','') else 'CCLE2018_Site_Primary',
                         'tcga_code_2018':l.get('tcga_code',''),
                         'sanger_model_id_broad':b.get('SangerModelID','')})
    e=pd.DataFrame(expr_rows)
    duplicated=e.canonical_model_id.ne('') & e.canonical_model_id.duplicated(keep=False)
    e.loc[duplicated,'mapping_status']='duplicate_expression_model_excluded'
    e['expression_eligible']=e.mapping_status.eq('resolved_exact_accession').astype(int)
    e.to_csv(OUT/'cell_line_dictionary.csv',index=False,encoding='utf-8-sig')
    edict=e[e.expression_eligible.eq(1)].set_index('canonical_model_id').to_dict('index')
    g=read(OUT/'gdsc2_observation_index.csv')
    identity_cols=['sanger_model_id','cosmic_id','cell_line_name','tcga_desc']
    unique=g[identity_cols].drop_duplicates()
    grows=[]
    for _,r in unique.iterrows():
        by_sid=sid_ach[r.sanger_model_id]
        by_cosmic=cosmic_ach[id_value(r.cosmic_id)]
        candidates=by_sid|by_cosmic
        a=next(iter(candidates)) if len(candidates)==1 else ''
        x=edict.get(a,{})
        status='resolved_exact_ID' if a else 'conflicting_IDs' if candidates else 'unmapped_ID'
        grows.append({**r.to_dict(),'canonical_model_id':a,'biological_group_id':('BIOGROUP_'+uf.find(a)) if a else '',
                      'ccle_id':x.get('ccle_id',''),'expression_row_index':x.get('expression_row_index',''),
                      'tissue':x.get('tissue',''),'mapping_status':status,
                      'mapping_evidence':'SANGER_MODEL_ID;COSMIC_ID' if by_sid and by_cosmic else 'SANGER_MODEL_ID' if by_sid else 'COSMIC_ID' if by_cosmic else '',
                      'candidate_model_ids':'|'.join(sorted(candidates)), 'expression_eligible':int(bool(x))})
    gm=pd.DataFrame(grows)
    gm.to_csv(OUT/'gdsc2_cell_mapping.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(edges,columns=['model_a','model_b','relation','source_id']).drop_duplicates().to_csv(OUT/'biological_group_edges.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame([{'canonical_model_id':a,'biological_group_id':'BIOGROUP_'+uf.find(a)} for a in sorted(uf.parent)]).to_csv(OUT/'canonical_model_groups.csv',index=False,encoding='utf-8-sig')
    # Export all authoritative accession relationships for CTRP reconciliation.
    cross=[{'ccle_id':c,'canonical_model_id':next(iter(v)) if len(v)==1 else '',
            'candidate_model_ids':'|'.join(sorted(v)), 'mapping_status':'resolved' if len(v)==1 else 'ambiguous'} for c,v in sorted(ccle_ach.items())]
    pd.DataFrame(cross).to_csv(OUT/'ccle_accession_crosswalk.csv',index=False,encoding='utf-8-sig')
    sizes=e[e.expression_eligible.eq(1)].groupby('biological_group_id').size()
    summary={'expression_samples':len(e),'eligible_expression_samples':int(e.expression_eligible.sum()),
             'expression_mapping_status_counts':e.mapping_status.value_counts().to_dict(),
             'expression_biological_groups':len(sizes),'groups_with_multiple_expression_models':int((sizes>1).sum()),
             'gdsc_unique_identity_rows':len(gm),'gdsc_unique_sanger_models':gm.sanger_model_id.nunique(),
             'gdsc_mapping_status_counts':gm.mapping_status.value_counts().to_dict(),
             'gdsc_expression_matched_models':gm[gm.expression_eligible.eq(1)].canonical_model_id.nunique(),
             'relationships_known_only':True,'note':'Unknown biological relationships may remain. Broad PatientID is same tissue sample, augmented by Sanger patient/parent links.'}
    (ROOT/'results'/'cell_identity_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    assert not e[e.expression_eligible.eq(1)].canonical_model_id.duplicated().any()
    assert gm.groupby('sanger_model_id').canonical_model_id.nunique().max() == 1
    print(json.dumps(summary,ensure_ascii=False))

if __name__=='__main__': main()
