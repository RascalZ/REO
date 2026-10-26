"""Apply only the two author-approved clarifications to new local reading copies.

This edits existing OOXML text nodes without regenerating layout, tables, images,
references or author declarations. No individual response or prediction is read.
"""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
import re
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/minor_clarification_v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8", newline="\n")


def edit_ooxml(source: Path, output: Path, edits: list[dict], language: str) -> list[str]:
    changed = []
    applied = {e["id"]: 0 for e in edits}
    alt_count = 0
    metadata_count = 0
    with ZipFile(source) as zin, ZipFile(output, "w") as zout:
        for member in zin.infolist():
            raw = zin.read(member.filename)
            result = raw
            if member.filename in {"word/document.xml", "word/header1.xml"}:
                xml = raw.decode("utf-8")

                def replace_text(match: re.Match) -> str:
                    nonlocal metadata_count
                    content = html.unescape(match.group(2))
                    original = content
                    for edit in edits:
                        if edit["anchor"] in content:
                            assert content.count(edit["anchor"]) == 1
                            content = content.replace(edit["anchor"], edit["replacement"], 1)
                            applied[edit["id"]] += 1
                    if language == "zh":
                        for before, after in [
                            ("英文通读v3的中文直接翻译", "英文通读v4的中文直接翻译"),
                            ("中文修订通读v2 · 源：英文通读v3", "中文修订通读v3 · 源：英文通读v4"),
                        ]:
                            if before in content:
                                content = content.replace(before, after, 1)
                                metadata_count += 1
                    if content == original:
                        return match.group(0)
                    return match.group(1) + html.escape(content, quote=False) + match.group(3)

                xml = re.sub(r"(<w:t(?:\s[^>]*)?>)([^<]*)(</w:t>)", replace_text, xml)

                def replace_alt(match: re.Match) -> str:
                    nonlocal alt_count
                    content = html.unescape(match.group(1))
                    caption = edits[1]
                    if content == caption["old_block"]:
                        alt_count += 1
                        return 'descr="' + html.escape(caption["new_block"], quote=True) + '"'
                    return match.group(0)

                xml = re.sub(r'descr="([^"]*)"', replace_alt, xml)
                result = xml.encode("utf-8")
            if result != raw:
                changed.append(member.filename)
            zout.writestr(member, result)
    assert all(value == 1 for value in applied.values()), applied
    assert alt_count == 1, alt_count
    assert metadata_count == (2 if language == "zh" else 0), metadata_count
    assert changed == (["word/document.xml", "word/header1.xml"] if language == "zh"
                       else ["word/document.xml"]), changed
    return changed


def main() -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    assert not (AUDIT / "build_manifest.json").exists(), "Preserve the prior build record"
    protected = {}
    historical_differences = []
    ignored_transient_locks = []
    historical = json.loads((ROOT / "docs/author_feedback_r1/preservation_manifest.json")
                            .read_text(encoding="utf-8"))["files"]
    for name, previous_sha in historical.items():
        path = ROOT / name
        if path.name.startswith("~$"):
            ignored_transient_locks.append(name)
            continue
        current_sha = sha(path)
        if current_sha != previous_sha:
            historical_differences.append({"file": name, "historical_sha256": previous_sha,
                                           "pre_edit_sha256": current_sha,
                                           "status": "preexisting_difference_preserved_not_restored"})
        protected[name] = current_sha
    reviewed = json.loads((ROOT / "docs/reviews/skill_panel_v1/intake_manifest.json")
                          .read_text(encoding="utf-8"))["files"]
    assert all(sha(ROOT / name) == expected for name, expected in reviewed.items())
    protected.update(reviewed)
    for folder in ["author_review_v3", "author_review_v3_zh"]:
        for suffix in [".md", ".docx", ".pdf"]:
            for path in (ROOT / "manuscript_v1" / folder).glob("*" + suffix):
                if not path.name.startswith("~$"):
                    protected[path.relative_to(ROOT).as_posix()] = sha(path)
    for folder in ["precommitments", "review_reports"]:
        for path in (ROOT / "docs/reviews/skill_panel_v1" / folder).glob("*.md"):
            protected[path.relative_to(ROOT).as_posix()] = sha(path)
    for name in ["editorial_decision.md", "revision_roadmap.md", "revision_roadmap.json"]:
        path = ROOT / "docs/reviews/skill_panel_v1" / name
        protected[path.relative_to(ROOT).as_posix()] = sha(path)
    assert all(sha(ROOT / name) == expected for name, expected in protected.items())
    write_json(AUDIT / "preservation_manifest.json", {
        "files": protected, "preexisting_historical_differences": historical_differences,
        "ignored_transient_word_locks": ignored_transient_locks})

    specifications = {
        "en": [
            {
                "id": "AUC_ENDPOINT",
                "prefix": "The continuous outcome was",
                "anchor": "lower in-vitro drug sensitivity [26].",
                "addition": " External responses used the CTRP v2 post-QC sigmoid-curve field area_under_curve, retained on its released scale without additional normalization. For this viability-curve endpoint, larger AUC values were interpreted as lower sensitivity; cross-screen evaluation compared within-drug rankings rather than absolute AUC values.",
            },
            {
                "id": "LOTO_CAPTION_REASON",
                "prefix": "Figure 4. Defined leave-one-tissue-out",
                "anchor": "Only 34 of 80 drug-by-tissue contrasts are defined, so the strict equal-drug macro is undefined.",
                "addition": " In this 80-cell display grid, all 46 undefined contrasts involved constant held-out predictions from at least one zero-molecular solution.",
            },
        ],
        "zh": [
            {
                "id": "AUC_ENDPOINT",
                "prefix": "连续结局为来源数据报告的",
                "anchor": "即体外药物敏感性越低[26]。",
                "addition": "外部反应使用CTRP v2发布的经质量控制（post-QC）的sigmoid拟合曲线字段area_under_curve，保留发布尺度，不作额外归一化。对于这一细胞活力曲线端点，本研究将较大AUC解释为较低敏感性；跨筛选评价比较药物内排序，而非AUC绝对值。",
            },
            {
                "id": "LOTO_CAPTION_REASON",
                "prefix": "图4. 有定义的留一组织法",
                "anchor": "80个药物—组织差值中仅34个有定义，因此严格等权药物宏指标未定义。",
                "addition": "在这一80格展示网格中，全部46个未定义差值均涉及至少一个零分子解模型，其在留出组织内的预测为常数。",
            },
        ],
    }
    manifest = {"date": "2026-10-26", "scope": "Two author-approved reporting clarifications only",
                "protected_files": len(protected), "languages": {},
                "new_training_resampling_or_scoring": False, "individual_CTRP_access": False,
                "reviewer_roadmap_modified": False, "other_review_items_adjudicated": False}
    for lang, old_folder, new_folder in [
        ("en", "author_review_v3", "author_review_v4"),
        ("zh", "author_review_v3_zh", "author_review_v4_zh"),
    ]:
        suffix = "_zh" if lang == "zh" else ""
        source = ROOT / f"manuscript_v1/{old_folder}/Manuscript_BMC_Genomics_author_review_v3{suffix}.md"
        output = ROOT / f"manuscript_v1/{new_folder}/Manuscript_BMC_Genomics_author_review_v4{suffix}.md"
        assert not output.parent.exists(), "Never overwrite an existing reading copy"
        output.parent.mkdir()
        original = source.read_bytes().decode("utf-8")
        result = original
        edits = specifications[lang]
        for edit in edits:
            blocks = [line for line in original.splitlines() if line.startswith(edit["prefix"])]
            assert len(blocks) == 1
            edit["old_block"] = blocks[0]
            assert blocks[0].count(edit["anchor"]) == 1
            edit["replacement"] = edit["anchor"] + edit["addition"]
            edit["new_block"] = blocks[0].replace(edit["anchor"], edit["replacement"], 1)
            assert result.count(blocks[0]) == 1
            result = result.replace(blocks[0], edit["new_block"], 1)
        old_lines, new_lines = original.splitlines(), result.splitlines()
        assert len(old_lines) == len(new_lines)
        assert sum(a != b for a, b in zip(old_lines, new_lines)) == 2
        output.write_bytes(result.encode("utf-8"))
        members = edit_ooxml(source.with_suffix(".docx"), output.with_suffix(".docx"), edits, lang)
        manifest["languages"][lang] = {
            "source": source.relative_to(ROOT).as_posix(), "source_sha256": sha(source),
            "output": output.relative_to(ROOT).as_posix(), "md_sha256": sha(output),
            "docx_sha256": sha(output.with_suffix(".docx")), "changed_blocks": 2,
            "unchanged_nonempty_blocks": len([line for line in new_lines if line.strip()]) - 2,
            "changed_ooxml_members": members, "edits": edits,
        }
    assert all(sha(ROOT / name) == expected for name, expected in protected.items())
    manifest["builder_sha256"] = sha(Path(__file__))
    write_json(AUDIT / "build_manifest.json", manifest)
    print(json.dumps({"changed_blocks_per_language": 2, "protected_files": len(protected),
                      "new_fits_or_scoring": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
