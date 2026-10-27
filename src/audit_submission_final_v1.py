"""Verify the editorial boundary and render already-exported PDFs for review."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess
from zipfile import ZipFile

from lxml import etree as ET
from pypdf import PdfReader
import pdfplumber

from audit_public_access_revision_v1 import font_audit, render_pages
from build_submission_final_v1 import ROOT, DOCS, PRIVATE, SOURCE, DEST, STEM, W, section, text


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(value):
    return re.sub(r'[\s*`|\u200b\u00ad]', '', value)


def main():
    config = json.loads((PRIVATE/'authoring_config.json').read_text(encoding='utf-8'))
    preservation = json.loads((PRIVATE/'preservation_manifest.json').read_text(encoding='utf-8'))['files']
    assert all(sha(ROOT/path)==value for path,value in preservation.items())
    old_md=(SOURCE/'Manuscript_BMC_Genomics_author_review_v5.md').read_text(encoding='utf-8')
    md=(DEST/(STEM+'.md')).read_text(encoding='utf-8')
    science=['Abstract','Background','Methods','Results','Discussion','Conclusions','References']
    assert all(section(old_md,h)==section(md,h) for h in science)
    assert [s for s in old_md.splitlines() if s.startswith('|')]==[s for s in md.splitlines() if s.startswith('|')]
    assert not re.search(r'[\u3400-\u9fff]',md)
    for banned in ['Author reading copy','Draft status','Author-specific contribution roles will','will be finalized','[Author confirmation required']:
        assert banned not in md,banned
    assert 'QC2025C053' not in md
    assert 'LH2022C056' in md
    refs=[x for x in section(md,'References').splitlines() if re.match(r'^\d+\. ',x)]
    assert len(refs)==31
    abstract=section(md,'Abstract')
    abstract_words=len(re.findall(r'\S+',abstract.replace('**','')))
    assert abstract_words<=350,abstract_words

    with ZipFile(SOURCE/'Manuscript_BMC_Genomics_author_review_v5.docx') as old,ZipFile(DEST/(STEM+'.docx')) as new:
        assert old.namelist()==new.namelist()
        changed=[n for n in old.namelist() if old.read(n)!=new.read(n)]
        assert changed==['word/document.xml']
        xo,xn=[ET.fromstring(z.read('word/document.xml')) for z in [old,new]]
        tables=lambda x:[ET.tostring(n,method='c14n') for n in x.findall('.//w:tbl',W)]
        assert tables(xo)==tables(xn) and len(tables(xn))==2
        media=[n for n in old.namelist() if n.startswith('word/media/') and not n.endswith('/')]
        assert len(media)==4 and all(old.read(n)==new.read(n) for n in media)
        rows={row['old']:row for row in config['docx_replacements']}
        expected=[]
        untouched=0
        for p in xo.findall('.//w:body/w:p',W):
            old_text=text(p)
            if old_text in rows:
                if rows[old_text]['new'] is not None:
                    expected.append(rows[old_text]['new'])
            else:
                untouched+=1
                expected.append(old_text)
        actual=[text(p) for p in xn.findall('.//w:body/w:p',W)]
        assert expected==actual,'Unexpected paragraph change'
        assert len(xn.findall('.//w:lnNumType',W))==1
        assert not xn.findall('.//w:br[@w:type="page"]',W)
        assert not xn.findall('.//w:pageBreakBefore',W)
    build=json.loads((DOCS/'build_manifest.json').read_text(encoding='utf-8'))
    assert all(sha(DEST/path)==value for path,value in build['byte_identical_support_files'].items())
    summary={'date':config['date'],'unchanged_scientific_sections':science,
        'unchanged_main_tables':2,'unchanged_embedded_figures':4,'reference_count':31,
        'abstract_word_count':abstract_words,'unchanged_body_paragraphs':untouched,
        'changed_ooxml_members':changed,'original_v5_files_preserved':len(preservation),
        'unchanged_support_file_count':21,'all_author_approval_confirmed':config['all_author_approval_confirmed'],
        'author_only_items_pending':config['author_only_items_pending'],
        'renderer':'Isolated native Microsoft Word read-only PDF export and bundled Poppler PNGs',
        'skill_render_fallback_reason':'LibreOffice soffice.exe unavailable; native Word used after sandbox CScript registry-settings access failed',
        'scientific_results_changed':False,'training_or_scoring_performed':False,
        'individual_CTRP_access':False,'documents':{},'visual_review':'pending'}
    for stem in [STEM,'Cover_letter_BMC_Genomics']:
        pdf=DEST/(stem+'.pdf')
        reader=PdfReader(pdf)
        outside=[]; glyphs=[]; body=[]; figs=[]
        with pdfplumber.open(pdf) as doc:
            for i,p in enumerate(doc.pages,1):
                full=p.extract_text() or ''
                if stem==STEM:
                    assert re.search(r'Page\s+'+str(i)+r'\b',full),(i,'page numbering')
                    body.append(p.crop((70,45,p.width-50,p.height-40)).extract_text() or '')
                else:
                    body.append(full)
                outside.extend([i for c in p.chars if c['x0']<-.5 or c['x1']>p.width+.5 or c['top']<-.5 or c['bottom']>p.height+.5])
                glyphs.extend([i for c in p.chars if c['text'] in ['\ufffd','\u25a1','\u0000']])
                if p.images:
                    figs.append(i)
        assert not outside and not glyphs,(stem,outside,glyphs)
        full='\n'.join(body)
        if stem==STEM:
            for edit in config['docx_replacements']:
                if edit['new']:
                    assert normalize(edit['new']) in normalize(full),edit['id']
            assert len(figs)==4
            word_text=''.join(text(p) for p in xn.findall('.//w:body/w:p',W))
            assert 'Author reading copy' not in word_text
        else:
            for row in config['cover_letter']:
                assert normalize(row['text']) in normalize(full),'Cover paragraph missing in PDF'
            assert 'AUTHOR CONFIRMATION REQUIRED' not in full
        (PRIVATE/(stem+'_extracted.txt')).write_text(full,encoding='utf-8')
        # The shared helper creates full-size pages and contact sheets locally.
        info=render_pages({'output':(DEST/(stem+'.md')).relative_to(ROOT).as_posix()}) if stem==STEM else None
        if stem!=STEM:
            pages=DEST/'qa_cover/pages'; pages.mkdir(parents=True)
            ppm=Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
            subprocess.run([str(ppm),'-scale-to','2000','-png',str(pdf),str(pages/'page')],check=True,capture_output=True)
            info={'pages':len(reader.pages),'page_directory':pages.relative_to(ROOT).as_posix(),'page_size_max_pixels':2000}
        summary['documents'][stem]={'pdf_pages':len(reader.pages),'pdf_sha256':sha(pdf),
            'docx_sha256':sha(DEST/(stem+'.docx')),'text_outside_page':len(outside),
            'suspect_missing_glyphs':len(glyphs),'figure_pages':figs,
            'fonts':font_audit(reader),'render':info}
    with (DOCS/'mechanical_audit.json').open('x',encoding='utf-8') as f:
        json.dump(summary,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'main_pages':summary['documents'][STEM]['pdf_pages'],
        'cover_pages':summary['documents']['Cover_letter_BMC_Genomics']['pdf_pages'],
        'abstract_words':abstract_words,'original_v5_preserved':len(preservation),
        'scientific_results_changed':False,'visual_review':'pending'}))


if __name__=='__main__':
    main()
