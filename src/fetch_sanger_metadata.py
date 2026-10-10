"""Fetch public model identity/relationship metadata; no response/expression data."""
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from datetime import datetime, timezone
import json, hashlib, csv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data_raw' / 'sanger_metadata'
OUT.mkdir(parents=True, exist_ok=True)
query = {
    'include': 'identifiers.source,sample.patient,sample.tissue,parent',
    'fields[model]': 'names,identifiers,sample,parent,model_type,model_treatment',
    'fields[sample]': 'patient,tissue',
    'fields[patient]': 'id',
    'fields[tissue]': 'name',
    'page[size]': '0',
}
url = 'https://api.cellmodelpassports.sanger.ac.uk/models?' + urlencode(query)
dest = OUT / 'cell_model_passports_2026-10-10.json'
if not dest.exists():
    with urlopen(Request(url, headers={'User-Agent': 'AcademicDataCuration/1.0'}), timeout=55) as r:
        body = r.read()
        data = json.loads(body)
        assert isinstance(data.get('data'), list) and len(data['data']) > 1000
        dest.write_bytes(body)
        (OUT / 'download_manifest.json').write_text(json.dumps({
            'source_url': url, 'documentation_url': 'https://depmap.sanger.ac.uk/documentation/api/',
            'retrieved_utc': datetime.now(timezone.utc).isoformat(),
            'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest(),
            'models': len(data['data']), 'source_type': 'official_public_api_snapshot',
            'license_note': 'Sanger API public non-commercial research use; retain source attribution.'
        }, indent=2), encoding='utf-8')
data = json.loads(dest.read_text(encoding='utf-8'))
included = {(x['type'], str(x['id'])): x for x in data.get('included', [])}
def related(item, field):
    rel = item.get('relationships', {}).get(field, {}).get('data')
    return included.get((rel['type'], str(rel['id'])), {}) if rel else {}
records = []
identifiers = []
for model in data['data']:
    a = model.get('attributes', {})
    sample = related(model, 'sample')
    sr = sample.get('relationships', {})
    parent = model.get('relationships', {}).get('parent', {}).get('data') or {}
    row = {'sanger_model_id': model['id'], 'names': '|'.join(a.get('names', [])),
           'parent_sanger_model_id': parent.get('id', ''),
           'sanger_sample_id': sample.get('id', ''),
           'sanger_patient_id': (sr.get('patient', {}).get('data') or {}).get('id', ''),
           'sanger_tissue': related(sample, 'tissue').get('attributes', {}).get('name', ''),
           'model_type': a.get('model_type', ''), 'model_treatment': a.get('model_treatment', '')}
    for rel in model.get('relationships', {}).get('identifiers', {}).get('data', []):
        ident = included[(rel['type'], str(rel['id']))]
        source = related(ident, 'source')
        name = source.get('attributes', {}).get('name', str(source.get('id', 'unknown')))
        value = ident.get('attributes', {}).get('identifier', '')
        identifiers.append({'sanger_model_id': model['id'], 'identifier_source': name, 'identifier': value})
    records.append(row)
for name, rows in [('sanger_models.csv', records), ('sanger_identifiers.csv', identifiers)]:
    with (ROOT / 'data_processed' / name).open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print(json.dumps({'models': len(records), 'identifiers': len(identifiers),
                  'sources': sorted({x['identifier_source'] for x in identifiers})}))
