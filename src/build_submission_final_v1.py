"""Prepare a bounded English submission copy from v5 and local source facts.

Configuration and author records remain Git-ignored. No scientific data,
prediction records, model states or external response labels are opened.
Unconfirmed author facts are never converted into affirmative declarations.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
from zipfile import ZipFile

from docx import Document
from docx.shared import Inches, Pt
from lxml import etree as ET

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs/submission_final_v1'
PRIVATE = DOCS / 'local_private'
SOURCE = ROOT / 'manuscript_v1/author_review_v5'
DEST = ROOT / 'manuscript_v1/submission_final_v1'
STEM = 'Manuscript_BMC_Genomics_submission'
W = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, ensure_ascii=False, indent=2)
        f.write('\n')


def text(element):
    return ''.join(t.text or '' for t in element.findall('.//w:t', W))


def replace_text(element, value):
    nodes = element.findall('.//w:t', W)
    assert nodes
    nodes[0].text = value
    nodes[0].set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    for node in nodes[1:]:
        node.text = ''


def section(md, heading):
    return md.split('## '+heading+'\n', 1)[1].split('\n## ', 1)[0]


def main():
    config = json.loads((PRIVATE / 'authoring_config.json').read_text(encoding='utf-8'))
    assert not DEST.exists(), 'Preserve previously generated submission files'
    preserved = {}
    for folder in ['author_review_v5', 'author_review_v5_zh']:
        for path in sorted((ROOT/'manuscript_v1'/folder).rglob('*')):
            if path.is_file() and not path.name.startswith('~$'):
                preserved[path.relative_to(ROOT).as_posix()] = sha(path)
    for item in config['source_records']:
        path = ROOT.parent / item['relative_to_project']
        assert sha(path) == item['sha256'], 'Author source has changed'

    old_md = (SOURCE/'Manuscript_BMC_Genomics_author_review_v5.md').read_text(encoding='utf-8')
    replacements = config['replacements']
    revised = old_md
    for item in replacements:
        assert revised.count(item['old']) == 1, item['id']
        revised = revised.replace(item['old'], item['new'], 1)
    for heading in ['Abstract', 'Background', 'Methods', 'Results', 'Discussion', 'Conclusions']:
        assert section(old_md, heading) == section(revised, heading), heading
    # Reference 31 and availability are not rewritten in this bounded pass.
    assert section(old_md, 'References') == section(revised, 'References')
    DEST.mkdir()
    (DEST/(STEM+'.md')).write_text(revised, encoding='utf-8', newline='\n')

    old_docx = SOURCE/'Manuscript_BMC_Genomics_author_review_v5.docx'
    docx_changes = config['docx_replacements']
    applied = {row['id']: 0 for row in docx_changes}
    changed_members = []
    with ZipFile(old_docx) as source, ZipFile(DEST/(STEM+'.docx'), 'x') as target:
        for member in source.infolist():
            raw = source.read(member.filename)
            new = raw
            if member.filename == 'word/document.xml':
                x = ET.fromstring(raw)
                for p in list(x.findall('.//w:body/w:p', W)):
                    value = text(p)
                    for row in docx_changes:
                        if value != row['old']:
                            continue
                        applied[row['id']] += 1
                        if row['new'] is None:
                            p.getparent().remove(p)
                        else:
                            replace_text(p, row['new'])
                        break
                new = ET.tostring(x, xml_declaration=True, encoding='UTF-8', standalone=True)
            if raw != new:
                changed_members.append(member.filename)
            target.writestr(member, new)
    assert all(v == 1 for v in applied.values()), applied
    assert changed_members == ['word/document.xml']
    supports = {}
    for folder in ['Figures', 'Additional_files', 'Supporting_sources']:
        shutil.copytree(SOURCE/folder, DEST/folder)
        for path in (DEST/folder).rglob('*'):
            if path.is_file():
                supports[path.relative_to(DEST).as_posix()] = sha(path)
                assert sha(SOURCE/path.relative_to(DEST)) == sha(path)
    assert len(supports) == 21

    # New cover letter follows a plain journal correspondence layout.
    cover = Document()
    sec = cover.sections[0]
    sec.page_width, sec.page_height = Inches(8.27), Inches(11.69)
    sec.top_margin = sec.bottom_margin = Inches(0.8)
    sec.left_margin = sec.right_margin = Inches(0.9)
    normal = cover.styles['Normal']
    normal.font.name = 'Times New Roman'
    normal.font.size = Pt(11)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.05
    for row in config['cover_letter']:
        para = cover.add_paragraph(row['text'])
        if row.get('bold'):
            for run in para.runs:
                run.bold = True
    cover.save(DEST/'Cover_letter_BMC_Genomics.docx')
    (DEST/'Cover_letter_BMC_Genomics.txt').write_text(
        '\n\n'.join(row['text'] for row in config['cover_letter'])+'\n',
        encoding='utf-8', newline='\n')
    (DEST/'AUTHOR_CONFIRMATIONS_REQUIRED_作者待确认.md').write_text(
        config['author_checklist'], encoding='utf-8', newline='\n')
    assert all(sha(ROOT/path) == value for path, value in preserved.items())
    dump(PRIVATE/'preservation_manifest.json', {'files': preserved})
    dump(DOCS/'build_manifest.json', {
        'date': config['date'], 'source_version': 'v5',
        'scope': 'Submission packaging and source-backed author information only',
        'output_folder': DEST.relative_to(ROOT).as_posix(),
        'source_record_hashes': [row['sha256'] for row in config['source_records']],
        'source_record_locators': [row['locators'] for row in config['source_records']],
        'changed_ooxml_members': changed_members,
        'edit_ids': [row['id'] for row in replacements],
        'unchanged_scientific_sections': ['Abstract', 'Background', 'Methods', 'Results', 'Discussion', 'Conclusions', 'References'],
        'byte_identical_support_files': supports,
        'preserved_v5_file_count': len(preserved),
        'scientific_results_changed': False, 'training_or_scoring_performed': False,
        'individual_CTRP_access': False,
        'all_author_approval_confirmed': config['all_author_approval_confirmed'],
        'author_only_items_pending': config['author_only_items_pending'],
        'authoritative_details_publicly_copied': False,
    })
    print(json.dumps({'output': DEST.relative_to(ROOT).as_posix(),
                      'unchanged_scientific_sections': 7,
                      'preserved_v5_files': len(preserved),
                      'identical_support_files': len(supports),
                      'pending_author_items': len(config['author_only_items_pending'])}))


if __name__ == '__main__':
    main()
