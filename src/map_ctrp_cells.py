"""Resolve official CTRP cell names against named accessions, flag ambiguity/site conflict."""
from pathlib import Path
from collections import defaultdict
import json, re
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data_processed'
from project_paths import original_input
def read(p,sep=','): return pd.read_csv(p,sep=sep,dtype=str,keep_default_na=False)
def norm(x): return re.sub('[^A-Z0-9]','',str(x).upper())
def main():
    ctrp=read(ROOT/'data_raw/ctrp/v20.meta.per_cell_line.txt','\t')
    broad=read(ROOT/'data_raw/cell_metadata/DepMap_24Q4_v1_Model.csv')
    legacy=read(original_input('legacy_cell_annotations'),'\t')
    bridge=read(original_input('legacy_ctrp_bridge'),'\t')
    cross=read(OUT/'ccle_accession_crosswalk.csv')
    expr=read(OUT/'cell_line_dictionary.csv')
    groups=read(OUT/'canonical_model_groups.csv').set_index('canonical_model_id').biological_group_id.to_dict()
    ccle_candidates={r.ccle_id:set(r.candidate_model_ids.split('|'))-{''} for _,r in cross.iterrows()}
    lookup=defaultdict(set); sources=defaultdict(set)
    def add(name,ids,source):
        if name:
            lookup[norm(name)].update(ids); sources[norm(name)].add(source)
    for _,r in broad.iterrows():
        for col in ['CellLineName','StrippedCellLineName']:
            add(r[col],{r.ModelID},'Broad24Q4_'+col)
    for _,r in legacy.iterrows():
        add(r.Name,ccle_candidates.get(r.CCLE_ID,set()),'CCLE2018_Name')
    for c,v in ccle_candidates.items():
        add(c.split('_')[0],v,'exact_normalized_CCLE_accession_prefix')
    # A legacy bridge is supplementary evidence, never a fuzzy nearest-name match.
    for _,r in bridge.iterrows():
        add(r['CTRP-Cellline-Name'],ccle_candidates.get(r['CCLE-Cellline-Name'],set()),'legacy_bridge_crosschecked_against_accession')
    sites=defaultdict(set)
    for _,r in legacy.iterrows():
        for a in ccle_candidates.get(r.CCLE_ID,set()): sites[a].add(norm(r.Site_Primary))
    e=expr[expr.expression_eligible=='1'].set_index('canonical_model_id').to_dict('index')
    rows=[]
    for _,r in ctrp.iterrows():
        key=norm(r.ccl_name); ids=lookup[key]
        a=next(iter(ids)) if len(ids)==1 else ''
        has_independent_name_source=any(not x.startswith('legacy_bridge') for x in sources[key])
        known_sites=sites[a]-{'','NA','NS','UNKNOWN'} if a else set()
        source_site=norm(r.ccle_primary_site)
        site_ok=not a or not known_sites or not source_site or source_site in known_sites
        status='resolved_exact_normalized_name_accession' if a and has_independent_name_source and site_ok else 'site_conflict' if a and not site_ok else 'legacy_only_unconfirmed' if a else 'ambiguous' if ids else 'unmapped'
        accepted=status=='resolved_exact_normalized_name_accession'
        x=e.get(a,{}) if accepted else {}
        rows.append({**r.to_dict(),'canonical_model_id':a if accepted else '',
                     'candidate_model_ids':'|'.join(sorted(ids)),'mapping_status':status,
                     'mapping_evidence':'|'.join(sorted(sources[key])),
                     'expected_site_normalized':'|'.join(sorted(sites[a])) if a else '',
                     'site_check':'matched' if a and known_sites and source_site and site_ok else 'conflict' if not site_ok else 'metadata_unavailable',
                     'expression_eligible':int(bool(x)),'ccle_id':x.get('ccle_id',''),
                     'expression_row_index':x.get('expression_row_index',''),
                     'biological_group_id':groups.get(a,'') if accepted else '',
                     'tissue':x.get('tissue','')})
    result=pd.DataFrame(rows)
    result.to_csv(OUT/'ctrp_cell_mapping.csv',index=False,encoding='utf-8-sig')
    result[~result.mapping_status.str.startswith('resolved')].to_csv(OUT/'ctrp_cell_mapping_exclusions.csv',index=False,encoding='utf-8-sig')
    summary={'official_ctrp_models':len(result),'mapping_status_counts':result.mapping_status.value_counts().to_dict(),
             'expression_available_master_ccl_ids':int(result.expression_eligible.sum()),
             'expression_available_canonical_models':result[result.expression_eligible.eq(1)].canonical_model_id.nunique(),
             'expression_available_biological_groups':result[result.expression_eligible.eq(1)].biological_group_id.nunique(),
             'name_matching_note':'Only exact names after case/punctuation normalization and corroborated named accessions; primary site checked when both sources provide it. Missing site metadata is recorded, not treated as a match. No fuzzy matching.'}
    (ROOT/'results/ctrp_cell_mapping_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False))
if __name__=='__main__': main()
