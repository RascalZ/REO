"""Check the published aggregate panel and hashes without private or large inputs."""
from pathlib import Path
import csv,hashlib,json

ROOT=Path(__file__).resolve().parents[1]

def main():
    snapshot=json.loads((ROOT/'results/public_cohort_snapshot.json').read_text(encoding='utf-8'))
    for item in snapshot['files']:
        p=ROOT/item['path']
        # Public text hashes use normalized LF; Windows checkout CRLF is harmless.
        raw=p.read_bytes().replace(b'\r\n',b'\n')
        assert hashlib.sha256(raw).hexdigest()==item['sha256_lf'],f'Public snapshot changed: {p.name}'
    with (ROOT/'results/main_drugs.csv').open(encoding='utf-8-sig',newline='') as f: panel=list(csv.DictReader(f))
    summary=json.loads((ROOT/'results/cohort_summary.json').read_text(encoding='utf-8'))
    assert len(panel)==8 and len({r['ctrp_inchikey'] for r in panel})==8
    assert all(r['gdsc_name_derived_inchikey']==r['ctrp_inchikey'] for r in panel)
    assert all(int(r['n_development_groups'])>=200 and int(r['n_external_primary_groups'])>=80 and int(r['n_development_tissues'])>=4 for r in panel)
    assert sum(int(r['n_development_models']) for r in panel)==summary['development_drug_model_observations']
    assert sum(int(r['n_external_primary_models']) for r in panel)==summary['external_drug_model_observations']
    assert summary['model_overlap']==summary['biological_group_overlap']==0
    assert summary['models_fitted']==0 and snapshot['local_identity_validation_status']=='passed'
    print('PASS: public phase-1 panel/hash/count consistency. Full identity validation requires local inputs.')

if __name__=='__main__':main()
