"""Verify only the permitted v4-to-v5 editorial edits; render PDF pages for review.

This opens manuscript artifacts and publication receipts only. It never reads
individual drug-response data, predictions, source indexes or fitted checkpoints.
"""
from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import pdfplumber
from pypdf import PdfReader
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / 'docs/public_access_revision_v1'
PRIVATE = AUDIT / 'local_private'
PPM = Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
W = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
SCIENTIFIC_COMMITS = ['5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43',
                      'fedc745940374b6870cb6c4bcb66ec8f731ded27',
                      'ca0f837ac10665ab8d4d6ab97691e58e5bed4577']
RELEASE = 'https://github.com/RascalZ/REO/releases/tag/v0.1.0'
DOI = 'https://doi.org/10.5281/zenodo.23681422'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plain(text: str) -> str:
    return text.replace('**','')


def normalize(text: str) -> str:
    return re.sub(r'[\s*`|\u200b\u00ad]', '', text)


def numeric_tokens(text: str):
    text = text.replace('October','10').replace('September','9')
    return Counter(re.findall(r'(?<![A-Za-z])[+−-]?\d+(?:,\d{3})*(?:\.\d+)?', text))


def undo_word(raw: bytes, edits: list[dict], lang: str, member: str, source_raw: bytes) -> bytes:
    targets = {plain(e['new_block']):plain(e['old_block']) for e in edits}
    if lang == 'en':
        targets['REVISED AUTHOR READING COPY — 27 OCTOBER 2026'] = 'REVISED AUTHOR READING COPY — 26 OCTOBER 2026'
    else:
        targets['英文通读v5的中文直接翻译'] = '英文通读v4的中文直接翻译'
        targets['中文修订通读v5 · 源：英文通读v5 · 2026-10-27'] = '中文修订通读v3 · 源：英文通读v4 · 2026-10-26'
    counts = Counter()
    source_paragraphs = {}
    for match in re.finditer(r'(<w:p(?:\s[^>]*)?>)(.*?)(</w:p>)',source_raw.decode('utf-8'),flags=re.S):
        text=''.join(html.unescape(m.group(2)) for m in
                     re.finditer(r'(<w:t(?:\s[^>]*)?>)([^<]*)(</w:t>)',match.group(2)))
        source_paragraphs.setdefault(text,[]).append(match.group(0))

    def without_text(paragraph_xml):
        return re.sub(r'(<w:t(?:\s[^>]*)?>)([^<]*)(</w:t>)',
                      lambda m:m.group(1)+'TEXT'+m.group(3),paragraph_xml)

    def paragraph(match):
        body = match.group(2)
        nodes = list(re.finditer(r'(<w:t(?:\s[^>]*)?>)([^<]*)(</w:t>)', body))
        texts = [html.unescape(x.group(2)) for x in nodes]
        text = ''.join(texts)
        replacement = targets.get(text)
        if lang == 'zh' and text.startswith('译稿日期：2026年10月27日。'):
            replacement = text.replace('译稿日期：2026年10月27日。','译稿日期：2026年10月26日。',1)
        if replacement is None:
            return match.group(0)
        # Original generators may use &apos;/&quot; entity spellings while the
        # bounded patch uses literal quotes. Compare preserve-only structure,
        # then restore the exact original paragraph bytes for the scope check.
        candidates=source_paragraphs.get(replacement,[])
        assert len(candidates)==1,(lang,member,'source paragraph ambiguity')
        restored=candidates[0]
        assert without_text(restored)==without_text(match.group(0)),(lang,member,'paragraph structure changed')
        counts[text] += 1
        return restored

    restored = re.sub(r'(<w:p(?:\s[^>]*)?>)(.*?)(</w:p>)',paragraph,raw.decode('utf-8'),flags=re.S)
    if member == 'word/document.xml':
        assert all(counts[plain(e['new_block'])] == 1 for e in edits)
        if lang == 'en':
            assert counts['REVISED AUTHOR READING COPY — 27 OCTOBER 2026'] == 1
        else:
            assert counts['英文通读v5的中文直接翻译'] == 1
            assert sum(v for k,v in counts.items() if k.startswith('译稿日期：2026年10月27日。')) == 1
    return restored.encode('utf-8')


def section(text: str, title: str) -> str:
    return text.split('## '+title+'\n',1)[1].split('\n## ',1)[0]


def font_audit(reader):
    result = {}
    for page in reader.pages:
        resources = page.get('/Resources')
        if not resources:
            continue
        fonts = resources.get_object().get('/Font',{})
        for font in fonts.get_object().values() if hasattr(fonts,'get_object') else fonts.values():
            font = font.get_object()
            descendants = font.get('/DescendantFonts')
            descriptor = font.get('/FontDescriptor')
            if descendants:
                descriptor = descendants[0].get_object().get('/FontDescriptor')
            embedded = False
            if descriptor:
                desc = descriptor.get_object()
                embedded = any(k in desc for k in ['/FontFile','/FontFile2','/FontFile3'])
            name = str(font.get('/BaseFont'))
            result[name] = {'embedded':embedded,'unicode_map':bool(font.get('/ToUnicode')),
                            'subtype':str(font.get('/Subtype'))}
    assert result and all(v['embedded'] for v in result.values()), result
    return result


def render_pages(info):
    md = ROOT/info['output']
    pages = md.parent/'qa/pages'
    pages.mkdir(parents=True,exist_ok=True)
    assert not list(pages.glob('page-*.png')), 'Never overwrite earlier render images'
    subprocess.run([str(PPM),'-scale-to','2000','-png',str(md.with_suffix('.pdf')),str(pages/'page')],
                   check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    images = sorted(pages.glob('page-*.png'))
    qa = pages.parent
    sheets = []
    for start in range(0,len(images),12):
        group = images[start:start+12]
        sheet = Image.new('RGB',(1440,3*530),(220,220,220))
        draw = ImageDraw.Draw(sheet)
        for j,path in enumerate(group):
            with Image.open(path) as im:
                im.thumbnail((350,500))
                x,y = (j%4)*360,(j//4)*530
                sheet.paste(im,(x+(350-im.width)//2,y+23))
            draw.text((x+8,y+4),'Page '+str(start+j+1),fill='black')
        target=qa/f'contact-{start+1:02d}-{start+len(group):02d}.png'
        assert not target.exists()
        sheet.save(target)
        sheets.append(target.relative_to(ROOT).as_posix())
    return {'pages':len(images),'page_directory':pages.relative_to(ROOT).as_posix(),
            'page_size_max_pixels':2000,'contact_sheets':sheets}


def main():
    assert not (AUDIT/'mechanical_audit.json').exists(), 'Preserve existing audit'
    build=json.loads((AUDIT/'build_manifest.json').read_text(encoding='utf-8'))
    protection=json.loads((AUDIT/'preservation_manifest.json').read_text(encoding='utf-8'))
    edit_specs=json.loads((PRIVATE/'edit_specifications.json').read_text(encoding='utf-8'))
    assert all(sha(ROOT/name)==value for name,value in protection['files'].items())
    assert all(sha(ROOT/name)==value for name,value in build['preserved_support_files'].items())
    assert len(protection['files']) == 91 and len(build['preserved_support_files']) == 20
    summary={'date':'2026-10-27','scope':'Bounded bilingual access/licensing/citation revision',
        'new_training_resampling_scoring':False,'individual_CTRP_access':False,
        'preserved_original_files':len(protection['files']),'preserved_support_files':20,
        'renderer':'Native Microsoft Word read-only export and bundled Poppler PNG rasterization',
        'skill_renderer_attempt':'render_docx.py failed because LibreOffice soffice.exe was absent; its temporary-directory cleanup also raised PermissionError. Native Word fallback used.',
        'languages':{},'visual_review_status':'pending','audit_script_sha256':sha(Path(__file__))}
    paired={}
    private_pdf_text={}
    for lang,info in build['languages'].items():
        old,new=ROOT/info['source'],ROOT/info['output']
        edits=edit_specs[lang]
        assert sha(new)==info['md_sha256'] and sha(new.with_suffix('.docx'))==info['docx_sha256']
        a,b=old.read_text(encoding='utf-8'),new.read_text(encoding='utf-8')
        undone=b
        for edit in edits:
            assert undone.count(edit['new_block'])==1
            undone=undone.replace(edit['new_block'],edit['old_block'],1)
        assert undone==a
        assert len(a.splitlines())==len(b.splitlines())
        assert sum(x!=y for x,y in zip(a.splitlines(),b.splitlines()))==5
        assert len([x for x in a.splitlines() if x.strip()])-5 == 192
        assert [x for x in a.splitlines() if x.startswith('|')]==[x for x in b.splitlines() if x.startswith('|')]
        titles=['Abstract','Results','Discussion','Conclusions'] if lang=='en' else ['摘要','结果','讨论','结论']
        assert all(section(a,t)==section(b,t) for t in titles)
        oldrefs=section(a,'References' if lang=='en' else '参考文献').splitlines()
        newrefs=section(b,'References' if lang=='en' else '参考文献').splitlines()
        oldlist=[x for x in oldrefs if re.match(r'^\d+\. ',x)]
        newlist=[x for x in newrefs if re.match(r'^\d+\. ',x)]
        assert len(oldlist)==len(newlist)==31 and oldlist[:30]==newlist[:30]
        assert DOI in newlist[30] and 'Zhang N, Zhang H, Huang W, Yang S, Zhang Y, Li J, et al.' in newlist[30]
        assert all(commit in a and commit in b for commit in SCIENTIFIC_COMMITS)
        assert RELEASE in b and DOI in b and 'CC BY 4.0' in b and 'MIT License' in b
        for banned in ['e4820b333655b59e47dd944ca1a2e08a0d68b2bd','currently private','目前为私有']:
            assert banned not in b
        with ZipFile(old.with_suffix('.docx')) as zo,ZipFile(new.with_suffix('.docx')) as zn:
            assert zo.namelist()==zn.namelist()
            changed=[name for name in zo.namelist() if zo.read(name)!=zn.read(name)]
            assert changed==info['changed_ooxml_members']
            assert all(undo_word(zn.read(name),edits,lang,name,zo.read(name))==zo.read(name) for name in changed)
            media=[name for name in zo.namelist() if name.startswith('word/media/') and not name.endswith('/')]
            assert len(media)==4 and all(zo.read(n)==zn.read(n) for n in media)
            xo,xn=[ET.fromstring(z.read('word/document.xml')) for z in [zo,zn]]
            assert [ET.tostring(x) for x in xo.findall('.//w:tbl',W)]==[ET.tostring(x) for x in xn.findall('.//w:tbl',W)]
            assert len(xn.findall('.//w:tbl',W))==2
            word_text=''.join(t.text or '' for t in xn.findall('.//w:t',W))
            assert all(normalize(plain(e['new_block'])) in normalize(word_text) for e in edits)
        pdf=new.with_suffix('.pdf')
        assert pdf.exists()
        reader=PdfReader(pdf)
        body_pages=[]; full_pages=[]; image_pages=[]; outside=[]; glyphs=[]
        with pdfplumber.open(pdf) as document:
            for i,page in enumerate(document.pages,1):
                full=page.extract_text() or ''
                full_pages.append(full)
                if lang=='en':
                    text=page.crop((70,45,page.width-50,page.height-40)).extract_text() or ''
                    assert re.search(r'Page\s+'+str(i)+r'\b',full),(lang,i,'page number')
                else:
                    text=full
                    assert re.search(r'第\s*'+str(i)+r'\s*页\s*/\s*共\s*'+str(len(reader.pages))+r'\s*页',full),(lang,i,'page number')
                body_pages.append(text)
                outside.extend([i for c in page.chars if c['x0']<-.5 or c['x1']>page.width+.5 or c['top']<-.5 or c['bottom']>page.height+.5])
                glyphs.extend([i for c in page.chars if c['text'] in ['\ufffd','\u25a1','\u0000']])
                for img in page.images:
                    assert img['x0']>=0 and img['x1']<=page.width+.5
                    assert img['top']>=0 and img['bottom']<=page.height+.5
                    image_pages.append(i)
        assert not outside and not glyphs and len(image_pages)==4,(lang,outside,glyphs,image_pages)
        pdf_body=normalize('\n'.join(body_pages))
        affected={}
        for edit in edits:
            target=normalize(plain(edit['new_block']))
            assert target in pdf_body,(lang,edit['id'],'new block not found in PDF')
            prefix=target[:min(75,len(target))]
            pages=[i for i,p in enumerate(body_pages,1) if prefix in normalize(p)]
            assert pages,(lang,edit['id'],'page location')
            affected[edit['id']]=pages
        assert all(s in pdf_body for s in SCIENTIFIC_COMMITS)
        assert '2026-10-27' in ''.join(full_pages) if lang=='zh' else True
        paired[lang]=[line for line in b.splitlines() if line.strip()]
        private_pdf_text[lang]=body_pages
        summary['languages'][lang]={'pdf_pages':len(reader.pages),'pdf_sha256':sha(pdf),
            'word_sha256':sha(new.with_suffix('.docx')),'changed_MD_blocks':5,'unchanged_MD_blocks':192,
            'word_package_reversal_byte_identical':True,'unchanged_tables':2,'unchanged_media':4,
            'unchanged_first_references':30,'reference_count':31,'all_scientific_snapshot_hashes_preserved':True,
            'unchanged_scientific_sections':titles,'text_outside_page':0,'suspect_missing_glyphs':0,
            'page_numbers_continuous':True,'fonts':font_audit(reader),'figure_pages':image_pages,
            'updated_block_start_pages':affected,'new_blocks_present_in_MD_DOCX_PDF':True}
    assert len(paired['en'])==len(paired['zh'])==197
    for i,(en,zh) in enumerate(zip(paired['en'],paired['zh']),1):
        assert numeric_tokens(en)==numeric_tokens(zh),(i,'bilingual numeric mismatch')
        assert re.findall(r'\[[\d,]+\]',en)==re.findall(r'\[[\d,]+\]',zh),(i,'bilingual citation mismatch')
    old_bib=(ROOT/'manuscript_v1/author_review_v4/Supporting_sources/references.bib').read_text(encoding='utf-8')
    new_bib=(ROOT/'manuscript_v1/author_review_v5/Supporting_sources/references.bib').read_text(encoding='utf-8')
    anchor='@misc{REOReconstruction2026,'
    assert old_bib.split(anchor)[0]==new_bib.split(anchor)[0]
    assert new_bib.count('@')==old_bib.count('@')==31 and DOI in new_bib and RELEASE in new_bib
    creators=list(csv.DictReader((ROOT/'docs/public_release_v1/CREATORS.csv').open(encoding='utf-8')))
    assert len(creators)==8
    assert all(row['family_name']+', '+row['given_name'] in new_bib for row in creators)
    assert sha(ROOT/'manuscript_v1/author_review_v5/Supporting_sources/references.bib')==build['bibliography']['sha256']
    with ThreadPoolExecutor(max_workers=2) as pool:
        renders=list(pool.map(render_pages,build['languages'].values()))
    for lang,render in zip(build['languages'],renders):
        assert render['pages']==summary['languages'][lang]['pdf_pages']
        summary['languages'][lang]['render']=render
    summary['bilingual_numeric_and_citation_blocks_checked']=197
    summary['bibliography_first_30_unchanged']=True
    summary['software_reference_authors_version_DOI_match']=True
    summary['mechanical_checks_passed']=True
    assert all(sha(ROOT/name)==value for name,value in protection['files'].items())
    with (AUDIT/'mechanical_audit.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(summary,f,ensure_ascii=False,indent=2);f.write('\n')
    with (PRIVATE/'rendered_body_text.json').open('x',encoding='utf-8',newline='\n') as f:
        json.dump(private_pdf_text,f,ensure_ascii=False)
    print(json.dumps({'mechanical_checks':'passed','pages':{k:v['pdf_pages'] for k,v in summary['languages'].items()},
                      'updated_pages':{k:v['updated_block_start_pages'] for k,v in summary['languages'].items()},
                      'visual_review':'pending'},ensure_ascii=False))


if __name__=='__main__':
    main()
