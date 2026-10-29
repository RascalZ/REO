"""Package English submission v2 only after source, mechanical and visual checks."""
from pathlib import Path, PurePosixPath
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'manuscript_v1/submission_final_v2'
AUDIT=ROOT/'docs/scope_diagnostics_revision_v1'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    mechanical=json.loads((AUDIT/'mechanical_audit.json').read_text(encoding='utf-8'))
    visual=json.loads((AUDIT/'visual_review.json').read_text(encoding='utf-8'))
    content=json.loads((AUDIT/'content_audit.json').read_text(encoding='utf-8'))
    assert mechanical['mechanical_checks_passed'] and visual['status']=='passed' and content['passed']
    for kind in ['english','chinese','cover']:
        assert visual['pages_reviewed'][kind]==list(range(1,mechanical['documents'][kind]['pages']+1))
    protected=json.loads((AUDIT/'local_private/preservation_manifest.json').read_text(encoding='utf-8'))['files']
    assert all(digest((ROOT/p).read_bytes())==h for p,h in protected.items())
    paths={
        'Manuscript.docx':'Manuscript_BMC_Genomics_submission_v2.docx',
        'Cover_letter.docx':'Cover_letter_BMC_Genomics.docx',
        'Previews/Manuscript.pdf':'Manuscript_BMC_Genomics_submission_v2.pdf',
        'Previews/Cover_letter.pdf':'Cover_letter_BMC_Genomics.pdf',
        'Supporting_sources/Manuscript.md':'Manuscript_BMC_Genomics_submission_v2.md',
        'Supporting_sources/Cover_letter.txt':'Cover_letter_BMC_Genomics.txt',
        'Supporting_sources/Supplementary_methods.md':'Supporting_sources/Supplementary_methods.md',
        'Supporting_sources/Supplementary_methods.pdf':'Supporting_sources/Supplementary_methods.pdf',
        'Supporting_sources/references.bib':'Supporting_sources/references.bib',
        'Author_only/Author_details_for_submission.csv':'Author_details_for_submission.csv',
        'Author_only/AUTHOR_INFORMATION_CHECK.md':'AUTHOR_INFORMATION_CHECK_作者信息核对.md',
    }
    for n in range(1,15):
        name=f'Additional_files/Additional_file_{n}.'+('docx' if n==14 else 'csv')
        paths[name]=name
    for p in sorted((OUT/'Figures').glob('*.png')):paths['Figures/'+p.name]='Figures/'+p.name
    assert len(paths)==29
    payload={name:(OUT/path).read_bytes() for name,path in paths.items()}
    for kind,name in [('english','Manuscript.docx'),('cover','Cover_letter.docx')]:
        assert digest(payload[name])==mechanical['documents'][kind]['docx_sha256']
        pdfname='Previews/'+('Manuscript.pdf' if kind=='english' else 'Cover_letter.pdf')
        assert digest(payload[pdfname])==mechanical['documents'][kind]['pdf_sha256']
    note=('BMC Genomics formal English submission revision v2 — 29 October 2026\n\n'
          'This package supersedes the local submission_final_v1 reading/upload package.\n'
          'Scope: narrower complete-pair elastic-net claims, clearer primary vs matched\n'
          'comparisons, existing zero-fit diagnostics, Table 2 display, and Figure 4A.\n'
          'No retraining, retuning, resampling, rescoring or individual CTRP access.\n'
          'All performance estimates and 14 supplementary files are retained.\n\n'
          'UPLOAD AFTER EXTRACTION BY CATEGORY:\n'
          'Manuscript.docx; Cover_letter.docx; four Figures/ PNGs; 14 Additional_files/.\n'
          'Previews/ and Supporting_sources/ are reading/source aids, not additional uploads.\n'
          'Author_only/ contains confirmed contact information for system entry only.\n'
          'Do not publish or upload Author_only/ as scientific supplementary material.\n'
          'The separate Chinese Word/PDF is an author reading translation and is not\n'
          'included in this English journal package.\n\n'
          'Existing eight-author order and author-confirmed declarations are retained.\n'
          'The present bounded revision was authorized by the corresponding project user.\n'
          'Actual manuscript submission has not been performed; an author must submit\n'
          'and review the journal-generated PDF.\n\n'
          'Frozen scientific release: https://github.com/RascalZ/REO/releases/tag/v0.1.0\n'
          'Software DOI: https://doi.org/10.5281/zenodo.23681422\n'
          'The immutable release/model assets are unchanged. Subsequent presentation and\n'
          'Stage 7 aggregate-report corrections are documented in the repository.\n'
          'Code MIT; author-owned derived assets CC BY 4.0 within applicable rights;\n'
          'third-party resources retain their original terms.\n').encode('utf-8')
    payload['00_UPLOAD_GUIDE.txt']=note
    payload['SHA256SUMS.txt']=''.join(f'{digest(data)}  {name}\n' for name,data in sorted(payload.items())).encode('ascii')
    for name,data in payload.items():
        assert not PurePosixPath(name).is_absolute() and '..' not in PurePosixPath(name).parts
        if name.startswith('Figures/'):assert len(data)<10*1024*1024
        if name.startswith('Additional_files/'):assert len(data)<20*1024*1024
    package=OUT/'BMC_Genomics_submission_v2_2026-10-29.zip'
    assert not package.exists(),'Preserve earlier package'
    with ZipFile(package,'x',ZIP_DEFLATED,compresslevel=9) as z:
        for name,data in sorted(payload.items()):
            item=ZipInfo(name,date_time=(2026,10,29,0,0,0));item.compress_type=ZIP_DEFLATED
            item.create_system=3;item.external_attr=0o100644<<16
            z.writestr(item,data,compress_type=ZIP_DEFLATED,compresslevel=9)
    with ZipFile(package) as z:
        assert z.testzip() is None and len(z.namelist())==len(set(z.namelist()))==31
        assert sorted(z.namelist())==sorted(payload)
        assert all(digest(z.read(n))==digest(data) for n,data in payload.items())
        for line in z.read('SHA256SUMS.txt').decode().splitlines():
            h,n=line.split('  ',1);assert digest(z.read(n))==h
    (OUT/'00_UPLOAD_GUIDE.txt').write_bytes(note)
    (OUT/'SHA256SUMS.txt').write_bytes(payload['SHA256SUMS.txt'])
    report={'date':'2026-10-29','package':package.relative_to(ROOT).as_posix(),
            'sha256':digest(package.read_bytes()),'bytes':package.stat().st_size,
            'members':[{ 'name':n,'sha256':digest(data),'bytes':len(data)} for n,data in sorted(payload.items())],
            'member_count':31,'upload_file_count':20,'verified_checksum_entries':30,
            'crc_and_all_member_hashes_passed':True,'strict_allowlist_passed':True,
            'prior_sources_preserved':len(protected),'training_or_scoring_performed':False,
            'individual_CTRP_access':False,'scientific_estimates_changed':False,
            'author_contacts_publicly_uploaded':False,'manuscript_submitted':False}
    (AUDIT/'submission_package_manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['package','sha256','bytes','member_count']}))


if __name__=='__main__':main()
