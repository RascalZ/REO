"""Create bounded bilingual v5 copies from preserved v4; no scientific data reads."""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
import re
import shutil
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / 'docs/public_access_revision_v1'
COMMIT = '643e4c9d0bd8b598cb9d9e98a4003bdd6004ec45'
RELEASE = 'https://github.com/RascalZ/REO/releases/tag/v0.1.0'
DOI = 'https://doi.org/10.5281/zenodo.23681422'
REFERENCE = ('31. Zhang N, Zhang H, Huang W, Yang S, Zhang Y, Li J, et al. '
             'REO: Auditable drug-response research reconstruction. Version 0.1.0. '
             'Zenodo; 2026. '+DOI+'.')


def metadata_replacements(lang):
    if lang == 'en':
        return {'REVISED AUTHOR READING COPY — 26 OCTOBER 2026':
                'REVISED AUTHOR READING COPY — 27 OCTOBER 2026'}
    old_note = '译稿日期：2026年10月26日。以下按原稿逐段翻译，不改写或增删研究内容。中文版独立编页码，反馈请注明“中文PDF第X页”。药物名、作者姓名、参考文献题录及图内标签保留原文，图注和表头译为中文。'
    return {
        '英文通读v4的中文直接翻译': '英文通读v5的中文直接翻译',
        old_note: old_note.replace('2026年10月26日', '2026年10月27日'),
        '中文修订通读v3 · 源：英文通读v4 · 2026-10-26':
            '中文修订通读v5 · 源：英文通读v5 · 2026-10-27',
    }


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, obj) -> None:
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
        f.write('\n')


def plain(text: str) -> str:
    return text.replace('**', '')


def patch_docx(source: Path, target: Path, edits, lang: str):
    applied = {e['id']: 0 for e in edits}
    metadata_hits = []
    changed = []
    with ZipFile(source) as zin, ZipFile(target, 'x') as zout:
        for info in zin.infolist():
            raw = zin.read(info.filename)
            result = raw
            if info.filename in {'word/document.xml', 'word/header1.xml'}:
                xml = raw.decode('utf-8')

                def paragraph(match):
                    body = match.group(2)
                    nodes = list(re.finditer(r'(<w:t(?:\s[^>]*)?>)([^<]*)(</w:t>)', body))
                    texts = [html.unescape(m.group(2)) for m in nodes]
                    text = ''.join(texts)
                    replacement = None
                    for edit in edits:
                        if text == plain(edit['old_block']):
                            replacement = plain(edit['new_block'])
                            applied[edit['id']] += 1
                            break
                    meta = metadata_replacements(lang)
                    if text in meta:
                        replacement = meta[text]
                        metadata_hits.append(info.filename)
                    if replacement is None:
                        return match.group(0)
                    if len(texts) == 1:
                        new_texts = [replacement]
                    else:
                        # Preserve the Chinese bold status label and its run style.
                        assert len(texts) == 2 and texts[0] == '稿件状态：'
                        assert replacement.startswith(texts[0])
                        new_texts = [texts[0], replacement[len(texts[0]):]]
                    iterator = iter(new_texts)
                    body = re.sub(r'(<w:t(?:\s[^>]*)?>)([^<]*)(</w:t>)',
                                  lambda m: m.group(1)+html.escape(next(iterator), quote=False)+m.group(3), body)
                    return match.group(1)+body+match.group(3)

                xml = re.sub(r'(<w:p(?:\s[^>]*)?>)(.*?)(</w:p>)', paragraph, xml, flags=re.S)
                result = xml.encode('utf-8')
            if result != raw:
                changed.append(info.filename)
            zout.writestr(info, result)
    assert all(n == 1 for n in applied.values()), applied
    assert len(metadata_hits) == len(metadata_replacements(lang))
    expected = ['word/document.xml', 'word/header1.xml'] if lang == 'zh' else ['word/document.xml']
    assert changed == expected, changed
    return changed


def specifications(text: str, lang: str):
    if lang == 'en':
        specs = [
            ('STATUS', '**Draft status:', '**Draft status:** Author reading copy, 27 October 2026. Public-access and licensing update to v4; scientific content and other author declarations retained. Final author approval remains pending.'),
            ('GOVERNANCE_ACCESS', 'The repository contains source versions,', None),
            ('DERIVED_ASSET_ACCESS', 'Complete derived identity mappings,',
             'Original fold assignments, 64 fixed numeric models, candidate-gene and coefficient tables, preprocessing states, and minimal frozen identity and drug-membership projections are publicly downloadable from release v0.1.0 [31]. Model exports and fold files retain their original bytes, with checksums and zero-coefficient solutions preserved. The minimal metadata omits source experiment, curve and response pointers; source datasets and wider historical compatibility inputs must be obtained separately. These packages do not constitute all inputs for every historical curation or external rescoring step.'),
            ('CODE_ARCHIVE_LICENSE', 'The reconstruction code,',
             'The source code, protocols, aggregate results, figure/table generators and response-free inference examples are publicly available in the immutable REO v0.1.0 release ('+RELEASE+') [31], archived at '+DOI+'. The release binds source commit '+COMMIT+'. Newly authored code uses the MIT License. Author-owned model parameters and derived research assets use CC BY 4.0 to the extent that the authors hold the applicable rights; third-party datasets, annotations and identifiers retain their source terms and are not relicensed. The package includes a hand-specified toy example and a generated-input example for frozen-model inference; neither is a drug-efficacy evaluation.'),
            ('SOFTWARE_REFERENCE_31', '31. REO reconstruction repository.', REFERENCE),
        ]
        before = 'The REO project repository is currently private; its contents and remaining access limitations are specified in the data and code availability statements.'
        after = 'The repository is now public; versioned code, frozen model exports and minimal identity/fold assets are available through the release and archive specified in the data and code availability statements.'
    else:
        specs = [
            ('STATUS', '**稿件状态：', '**稿件状态：** 作者通读版，2026年10月27日。对v4作公开访问与许可更新；保留科学内容和其他作者声明。最终作者批准仍待完成。'),
            ('GOVERNANCE_ACCESS', '仓库包含源码版本、', None),
            ('DERIVED_ASSET_ACCESS', '完整衍生身份映射、',
             '原始分折定义、64个固定数值模型、候选基因与系数表、预处理状态，以及最小冻结身份和药物成员投影已通过v0.1.0发布版公开下载[31]。模型导出与分折文件保留原字节，并提供校验和及保留零系数解。最小元数据省略源实验、曲线和响应指针；来源数据集及更广泛的历史兼容输入须另行获取。这些资产包不构成每一步历史数据整理或外部重新评分所需的全部输入。'),
            ('CODE_ARCHIVE_LICENSE', '重构代码、',
             '源代码、方案、聚合结果、图表生成程序及不使用响应标签的推理示例已公开于不可变REO v0.1.0发布版（'+RELEASE+'）[31]，并归档于'+DOI+'。发布版绑定源码提交'+COMMIT+'。新编写代码采用MIT License。作者拥有权利的模型参数和派生研究资产，在作者持有适用权利的范围内采用CC BY 4.0；第三方数据集、注释和标识保留其来源条款，不被重新许可。资产包含手工设定的示例和使用生成输入的冻结模型推理示例；二者均不是药效评价。'),
            ('SOFTWARE_REFERENCE_31', '31. REO reconstruction repository.', REFERENCE),
        ]
        before = 'REO项目仓库目前为私有；其内容和剩余访问限制在数据及代码可用性声明中说明。'
        after = '仓库现已公开；版本化代码、冻结模型导出及最小身份/分折资产可通过数据及代码可用性声明中的发布版和归档获取。'
    result = []
    for edit_id, prefix, replacement in specs:
        matches = [line for line in text.splitlines() if line.startswith(prefix)]
        assert len(matches) == 1, (lang, edit_id)
        old = matches[0]
        if replacement is None:
            assert old.count(before) == 1
            replacement = old.replace(before, after)
        result.append({'id': edit_id, 'old_block': old, 'new_block': replacement})
    return result


def main():
    assert not (AUDIT/'build_manifest.json').exists(), 'Preserve prior output'
    for folder in ['author_review_v5', 'author_review_v5_zh']:
        assert not (ROOT/'manuscript_v1'/folder).exists(), 'Never overwrite a reading copy'
    AUDIT.mkdir(parents=True, exist_ok=True)
    private = AUDIT/'local_private'
    private.mkdir(exist_ok=False)
    protected = {}
    for folder in ['author_review_v4', 'author_review_v4_zh']:
        for path in sorted((ROOT/'manuscript_v1'/folder).rglob('*')):
            if path.is_file() and not path.name.startswith('~$'):
                protected[path.relative_to(ROOT).as_posix()] = sha(path)
    dump(AUDIT/'preservation_manifest.json', {'date':'2026-10-27', 'files':protected})
    assert json.loads((ROOT/'docs/public_release_v1/GITHUB_RELEASE_AUDIT.json').read_text(encoding='utf-8'))['immutable']
    assert json.loads((ROOT/'docs/public_release_v1/ZENODO_DOWNLOAD_AUDIT.json').read_text(encoding='utf-8'))['all_11_files_verified']
    manifest = {'date':'2026-10-27', 'scope':'Public access, archival citation and license only',
                'source_release_commit':COMMIT, 'new_training_resampling_or_scoring':False,
                'individual_CTRP_access':False, 'languages':{}, 'preserved_support_files':{}}
    private_edits = {}
    for lang in ['en','zh']:
        suffix = '_zh' if lang == 'zh' else ''
        source = ROOT/f'manuscript_v1/author_review_v4{suffix}/Manuscript_BMC_Genomics_author_review_v4{suffix}.md'
        output = ROOT/f'manuscript_v1/author_review_v5{suffix}/Manuscript_BMC_Genomics_author_review_v5{suffix}.md'
        output.parent.mkdir()
        original = source.read_bytes().decode('utf-8')
        edits = specifications(original, lang)
        revised = original
        for edit in edits:
            assert revised.count(edit['old_block']) == 1
            revised = revised.replace(edit['old_block'], edit['new_block'], 1)
        assert len(original.splitlines()) == len(revised.splitlines())
        assert sum(a != b for a,b in zip(original.splitlines(),revised.splitlines())) == 5
        output.write_bytes(revised.encode('utf-8'))
        changed = patch_docx(source.with_suffix('.docx'), output.with_suffix('.docx'), edits, lang)
        private_edits[lang] = edits
        manifest['languages'][lang] = {'source':source.relative_to(ROOT).as_posix(),
            'output':output.relative_to(ROOT).as_posix(), 'md_sha256':sha(output),
            'docx_sha256':sha(output.with_suffix('.docx')), 'changed_blocks':5,
            'unchanged_nonempty_blocks':len([x for x in revised.splitlines() if x.strip()])-5,
            'changed_ooxml_members':changed,
            'edits':[{'id':e['id'],'source_line':original.splitlines().index(e['old_block'])+1,
                      'old_sha256':hashlib.sha256(e['old_block'].encode()).hexdigest(),
                      'new_sha256':hashlib.sha256(e['new_block'].encode()).hexdigest()} for e in edits]}
    v4, v5 = ROOT/'manuscript_v1/author_review_v4', ROOT/'manuscript_v1/author_review_v5'
    for folder in ['Figures','Additional_files','Supporting_sources']:
        (v5/folder).mkdir()
        for path in sorted((v4/folder).iterdir()):
            if path.name == 'references.bib':
                continue
            assert path.is_file()
            target = v5/folder/path.name
            shutil.copy2(path,target)
            assert sha(path) == sha(target)
            manifest['preserved_support_files'][target.relative_to(ROOT).as_posix()] = sha(target)
    assert len(manifest['preserved_support_files']) == 20
    old_bib = (v4/'Supporting_sources/references.bib').read_text(encoding='utf-8')
    anchor = '@misc{REOReconstruction2026,'
    assert old_bib.count(anchor) == 1
    new_entry = '''@misc{REOReconstruction2026,
  author = {Zhang, Nan and Zhang, Huicheng and Huang, Wen and Yang, Siru and Zhang, Yutong and Li, Jinglin and Liu, Hui and Wang, Shuyuan},
  title = {REO: Auditable drug-response research reconstruction},
  version = {0.1.0},
  year = {2026},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.23681422},
  url = {https://doi.org/10.5281/zenodo.23681422},
  urldate = {2026-10-27},
  note = {Software and research-asset archive. GitHub release: https://github.com/RascalZ/REO/releases/tag/v0.1.0; source commit 643e4c9d0bd8b598cb9d9e98a4003bdd6004ec45. Code MIT; author-owned model and derived assets CC BY 4.0 to the extent applicable rights are held; third-party resources not relicensed.}
}
'''
    (v5/'Supporting_sources/references.bib').write_text(old_bib.split(anchor)[0]+new_entry, encoding='utf-8', newline='\n')
    manifest['bibliography'] = {'unchanged_first_entries':30, 'updated_key':'REOReconstruction2026',
                               'sha256':sha(v5/'Supporting_sources/references.bib')}
    dump(private/'edit_specifications.json', private_edits)
    assert all(sha(ROOT/name) == expected for name,expected in protected.items())
    manifest['preserved_original_file_count'] = len(protected)
    manifest['builder_sha256'] = sha(Path(__file__))
    dump(AUDIT/'build_manifest.json',manifest)
    print(json.dumps({'changed_blocks_per_language':5,'original_files_preserved':len(protected),
                      'support_files_byte_identical':20,'scientific_recomputation':False}))


if __name__ == '__main__':
    main()
