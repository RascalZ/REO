"""Package exactly the final English upload files, previews and author aids.

Preserves v5 and rejects unreviewed renders. This only copies existing
publication artifacts; it does not evaluate models or open response datasets.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'manuscript_v1/submission_final_v1'
DOCS=ROOT/'docs/submission_final_v1'
STEM='Manuscript_BMC_Genomics_submission'
STAMP='2026-10-27'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    mechanical=json.loads((DOCS/'mechanical_audit.json').read_text(encoding='utf-8'))
    visual=json.loads((DOCS/'visual_review.json').read_text(encoding='utf-8'))
    assert visual['status']=='passed'
    assert visual['main_pages_reviewed']==list(range(1,41))
    assert visual['cover_pages_reviewed']==[1]
    preservation=json.loads((DOCS/'local_private/preservation_manifest.json').read_text(encoding='utf-8'))['files']
    assert all(sha((ROOT/path).read_bytes())==value for path,value in preservation.items())
    readme=FOLDER/'00_UPLOAD_GUIDE.txt'
    sums=FOLDER/'SHA256SUMS.txt'
    package=FOLDER/f'BMC_Genomics_formal_submission_{STAMP}.zip'
    manifest=DOCS/'submission_package_manifest.json'
    assert not any(p.exists() for p in [readme,sums,package,manifest]),'Preserve earlier package'
    sources={
        'Manuscript.docx':FOLDER/(STEM+'.docx'),
        'Cover_letter.docx':FOLDER/'Cover_letter_BMC_Genomics.docx',
        'Previews/Manuscript.pdf':FOLDER/(STEM+'.pdf'),
        'Previews/Cover_letter.pdf':FOLDER/'Cover_letter_BMC_Genomics.pdf',
        'Supporting_sources/Manuscript.md':FOLDER/(STEM+'.md'),
        'Supporting_sources/Cover_letter.txt':FOLDER/'Cover_letter_BMC_Genomics.txt',
        'Supporting_sources/Supplementary_methods.md':FOLDER/'Supporting_sources/Supplementary_methods.md',
        'Supporting_sources/Supplementary_methods.pdf':FOLDER/'Supporting_sources/Supplementary_methods.pdf',
        'Supporting_sources/references.bib':FOLDER/'Supporting_sources/references.bib',
        'Author_only/AUTHOR_CONFIRMATIONS_REQUIRED.md':FOLDER/'AUTHOR_CONFIRMATIONS_REQUIRED_作者待确认.md',
        'Author_only/Author_details_for_submission.csv':FOLDER/'Author_details_for_submission.csv',
    }
    for n in range(1,15):
        name=f'Additional_file_{n}.'+('docx' if n==14 else 'csv')
        sources['Additional_files/'+name]=FOLDER/'Additional_files'/name
    for name in ['Figure1_study_design.png','Figure2_development_primary.png','Figure3_representation_and_zero_solutions.png','Figure4_loto_and_external.png']:
        sources['Figures/'+name]=FOLDER/'Figures'/name
    assert len(sources)==29
    payloads={}
    for name,path in sources.items():
        assert path.resolve().is_relative_to(FOLDER.resolve()) and path.is_file() and not path.is_symlink()
        payloads[name]=path.read_bytes()
        if name.startswith('Figures/'):
            assert len(payloads[name])<10*1024*1024
        if name.startswith('Additional_files/'):
            assert len(payloads[name])<20*1024*1024
    for stem,key in [(STEM,'Manuscript.docx'),('Cover_letter_BMC_Genomics','Cover_letter.docx')]:
        assert sha(payloads[key])==mechanical['documents'][stem]['docx_sha256']
    assert sha(payloads['Previews/Manuscript.pdf'])==mechanical['documents'][STEM]['pdf_sha256']
    assert sha(payloads['Previews/Cover_letter.pdf'])==mechanical['documents']['Cover_letter_BMC_Genomics']['pdf_sha256']
    note=(
        'BMC Genomics formal English submission materials\n27 October 2026\n\n'
        'Scientific content, figures, tables, analysis, and references are unchanged from v5.\n'
        'Current author-confirmed affiliation, address, funding, emails and declarations are adopted.\n'
        'Original v5 English/Chinese copies and ZIP remain unchanged.\n\n'
        'UPLOAD BY CATEGORY AFTER EXTRACTING THIS ZIP:\n'
        'Manuscript.docx: editable main manuscript.\n'
        'Cover_letter.docx: cover letter.\n'
        'Figures/: four original independent figures.\n'
        'Additional_files/: 13 aggregate CSV tables and Additional_file_14.docx.\n\n'
        'DO NOT automatically upload Previews/, Supporting_sources/, Author_only/, this guide\n'
        'or SHA256SUMS as publishable supplementary files. The PDFs are reading previews.\n'
        'The author email CSV contains contact information for submission-system entry only.\n\n'
        'CURRENT AUTHOR-ONLY CHECK: consult Author_only/AUTHOR_CONFIRMATIONS_REQUIRED.md.\n'
        'The eight-author manuscript list is preserved. Any unmatched email/name must be\n'
        'resolved by the authors before filling the submission system. No author is added\n'
        'or replaced by inference. Approval and exclusive-submission status were confirmed\n'
        'in this task; recheck them if circumstances change before actual submission.\n\n'
        'Published software release: https://github.com/RascalZ/REO/releases/tag/v0.1.0\n'
        'Software archive DOI: https://doi.org/10.5281/zenodo.23681422\n'
        'Code MIT; author-owned derived assets CC BY 4.0 to the extent of applicable rights.\n'
        'Third-party resources retain source terms; no new dataset redistribution.\n\n'
        'Actual manuscript submission has not been performed. An author must submit.\n'
        'Current guidelines: https://link.springer.com/journal/12864/submission-guidelines\n'
    ).encode('utf-8')
    payloads['00_UPLOAD_GUIDE.txt']=note
    payloads['SHA256SUMS.txt']=''.join(f'{sha(data)}  {name}\n' for name,data in sorted(payloads.items())).encode('ascii')
    assert len(payloads)==31
    for name in payloads:
        pure=PurePosixPath(name)
        assert not pure.is_absolute() and '..' not in pure.parts and '\\' not in name
    readme.write_bytes(note); sums.write_bytes(payloads['SHA256SUMS.txt'])
    with ZipFile(package,'x',ZIP_DEFLATED,compresslevel=9) as z:
        for name,data in sorted(payloads.items()):
            item=ZipInfo(name,date_time=(2026,10,27,0,0,0))
            item.compress_type=ZIP_DEFLATED
            item.create_system=3;item.external_attr=0o100644<<16
            z.writestr(item,data,compress_type=ZIP_DEFLATED,compresslevel=9)
    with ZipFile(package) as z:
        assert z.testzip() is None and sorted(z.namelist())==sorted(payloads)
        assert len(z.namelist())==len(set(z.namelist()))
        assert all(sha(z.read(name))==sha(data) for name,data in payloads.items())
        for line in z.read('SHA256SUMS.txt').decode('ascii').splitlines():
            checksum,name=line.split('  ',1)
            assert checksum==sha(z.read(name))
    assert all(sha((ROOT/path).read_bytes())==value for path,value in preservation.items())
    report={'date':STAMP,'package':package.relative_to(ROOT).as_posix(),
        'bytes':package.stat().st_size,'sha256':sha(package.read_bytes()),'member_count':31,
        'members':[{'name':name,'bytes':len(data),'sha256':sha(data)} for name,data in sorted(payloads.items())],
        'crc_and_member_hash_checks_passed':True,'SHA256SUMS_entries_verified':30,
        'strict_member_allowlist_verified':True,'upload_file_count':20,
        'all_author_approval_confirmed':mechanical['all_author_approval_confirmed'],
        'author_only_items_pending':mechanical['author_only_items_pending'],
        'submission_system_ready':not mechanical['author_only_items_pending'],
        'manuscript_submitted':False,'original_v5_file_count_preserved':len(preservation),
        'scientific_results_changed':False,'training_or_scoring_performed':False,
        'individual_CTRP_access':False,'author_contact_payload_publicly_uploaded':False}
    manifest.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['package','member_count','bytes','sha256','submission_system_ready']}))


if __name__=='__main__':
    main()
