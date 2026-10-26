"""Audit paired editorial revisions against the exact prior manuscripts and page renders."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET
import zipfile

import pdfplumber
from PIL import Image, ImageDraw
from audit_chinese_translation_v1 import nums, normalize

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/author_feedback_r1"
W = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(args: list[str]) -> bytes:
    return subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout


def section(text: str, name: str) -> str:
    return text.split(f"## {name}\n", 1)[1].split("\n## ", 1)[0].strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--content-only", action="store_true")
    parser.add_argument("--visual-confirmed", action="store_true")
    args = parser.parse_args()
    build = json.loads((DOCS / "build_manifest.json").read_text(encoding="utf-8"))
    preservation = json.loads((DOCS / "preservation_manifest.json").read_text(encoding="utf-8"))["files"]
    assert all(sha(ROOT / name) == want for name, want in preservation.items())
    summary = {"date": "2026-10-26", "scope": "Focused bilingual editorial and layout audit, not independent scientific replication",
               "protected_prior_files_unchanged": len(preservation), "abstract_words": build["abstract_words"],
               "scientific_recomputation_or_individual_CTRP_access": False, "languages": {},
               "manual_visual_inspection_confirmed": args.visual_confirmed,
               "audit_script_sha256": sha(Path(__file__)),
               "render_script_sha256": sha(ROOT / "scripts/render_author_feedback_r1.vbs")}
    paired = {}
    for lang, info in build["languages"].items():
        old = ROOT / info["source"]
        md = ROOT / info["output"]
        docx = md.with_suffix(".docx")
        assert sha(old) == info["source_sha256"] and sha(md) == info["output_md_sha256"] and sha(docx) == info["docx_sha256"]
        a, b = old.read_text(encoding="utf-8"), md.read_text(encoding="utf-8")
        aa, bb = a.splitlines(), b.splitlines()
        assert len(aa) == len(bb)
        actual_changed = [i + 1 for i, (x, y) in enumerate(zip(aa, bb)) if x != y]
        assert actual_changed == sorted(e["line"] for e in info["changed_blocks"])
        assert len(actual_changed) == 8
        assert section(a, "Declarations" if lang == "en" else "声明") == section(b, "Declarations" if lang == "en" else "声明")
        assert section(a, "References" if lang == "en" else "参考文献") == section(b, "References" if lang == "en" else "参考文献")
        assert [l for l in aa if l.startswith("|")] == [l for l in bb if l.startswith("|")]
        assert "earlier manuscript" not in b and "earlier project manuscript" not in b and "早期稿件" not in b
        paired[lang] = [l for l in bb if l.strip()]
        with zipfile.ZipFile(docx) as z:
            xml = ET.fromstring(z.read("word/document.xml"))
            body = normalize("".join(t.text or "" for t in xml.findall(".//w:t", W)))
            images = [hashlib.sha256(z.read(n)).hexdigest() for n in z.namelist() if n.startswith("word/media/") and n.endswith(".png")]
            original_images = [sha(p) for p in (ROOT / "figures/manuscript_m08_v1").glob("Figure*.png") if "_proof" not in p.stem]
            assert len(images) == 4 and sorted(images) == sorted(original_images)
            assert len(xml.findall(".//w:tbl", W)) == 2
            checked = 0
            for line in paired[lang]:
                if re.fullmatch(r"\|[\s|:-]+", line):
                    continue
                if lang == "en" and (line.startswith("# BMC ") or line in ["## Title page", "## Figure legends"]):
                    continue
                if lang == "zh" and line in ["## 标题页", "## 图注"]:
                    continue
                want = re.sub(r"^#{1,3} ", "", line)
                if lang == "zh" and line.startswith("|"):
                    want = want.replace(" to ", " 至 ")
                if re.match(r"^(Table [12]\. |表[12]\. )", want):
                    stop = re.match(r"^(Table [12]\. [^.]+\.)", want).end() if lang == "en" else want.index("。") + 1
                    assert normalize(want[:stop]) in body and normalize(want[stop:]) in body
                else:
                    assert normalize(want) in body, (lang, checked)
                checked += 1
        item = {"changed_blocks": 8, "unmodified_nonempty_blocks": len(paired[lang]) - 8,
                "word_content_blocks_checked": checked, "tables_unchanged": 2, "original_pngs_byte_identical": 4,
                "references_preserved": 31, "declarations_unchanged": True,
                "md_sha256": sha(md), "docx_sha256": sha(docx)}
        if not args.content_only:
            pdf = md.with_suffix(".pdf")
            qa = md.parent / "qa"
            qa.mkdir(exist_ok=True)
            run(["pdftotext", "-layout", str(pdf), str(qa / "main.txt")])
            text_pages = [p for p in (qa / "main.txt").read_text(encoding="utf-8").split("\f") if p.strip()]
            total = len(text_pages)
            footer = all(re.search(rf"第\s*{i}\s*页\s*/\s*共\s*{total}\s*页", p) if lang == "zh"
                         else re.search(rf"^\s*Page {i}\s*$", p, re.M) for i, p in enumerate(text_pages, 1))
            fonts = run(["pdffonts", str(pdf)]).decode(errors="replace").splitlines()[2:]
            embedded = bool(fonts) and all(re.search(r"\byes\s+(?:yes|no)\s+(?:yes|no)\s+\d+\s+\d+\s*$", s) for s in fonts if s.strip())
            outside, figures = 0, []
            with pdfplumber.open(pdf) as doc:
                for i, p in enumerate(doc.pages, 1):
                    outside += sum(c["x0"] < 0 or c["x1"] > p.width + .5 or c["top"] < 0 or c["bottom"] > p.height + .5 for c in p.chars)
                    figures.extend({"page": i, "within_page": im["x0"] >= 0 and im["x1"] <= p.width and im["top"] >= 0 and im["bottom"] <= p.height}
                                   for im in p.images)
            assert footer and embedded and outside == 0 and len(figures) == 4 and all(f["within_page"] for f in figures)
            assert pdf.stat().st_mtime >= docx.stat().st_mtime
            thumbs = qa / "thumbs"
            thumbs.mkdir(exist_ok=True)
            run(["pdftoppm", "-scale-to", "500", "-png", str(pdf), str(thumbs / "page")])
            paths = sorted(thumbs.glob("page-*.png"))
            assert len(paths) == total
            for start in range(0, total, 12):
                sheet = Image.new("RGB", (1200, 2110), "#cccccc")
                draw = ImageDraw.Draw(sheet)
                for slot, p in enumerate(paths[start:start + 12]):
                    im = Image.open(p).convert("RGB")
                    im.thumbnail((394, 500))
                    x, y = (slot % 3) * 400, (slot // 3) * 525
                    draw.text((x + 8, y + 3), f"{lang} page {start + slot + 1}", fill="black")
                    sheet.paste(im, (x + (400 - im.width) // 2, y + 20))
                sheet.save(qa / f"contact_{start // 12 + 1}.png")
            item.update({"pdf_sha256": sha(pdf), "pages": total, "consecutive_footers": footer,
                         "fonts_embedded": embedded, "text_outside_page": outside, "figures": figures})
        summary["languages"][lang] = item
    assert len(paired["en"]) == len(paired["zh"]) == 197
    for index, (en, zh) in enumerate(zip(paired["en"], paired["zh"]), 1):
        assert nums(en) == nums(zh), (index, "bilingual numeric mismatch")
        assert re.findall(r"\[[\d,]+\]", en) == re.findall(r"\[[\d,]+\]", zh)
    summary["paired_numeric_citation_blocks_passed"] = 197
    summary["all_assessed_checks_passed"] = True
    (DOCS / ("content_audit.json" if args.content_only else "final_audit.json")).write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
