"""Source-fidelity and reading-layout audit of the local Chinese manuscript translation."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET
import zipfile

import pdfplumber
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/translation_zh_v1"
OUT = ROOT / "manuscript_v1/author_review_v2_zh"
EN = ROOT / "manuscript_v1/author_review_v2/Manuscript_BMC_Genomics_author_review_v2.md"
STEM = "Manuscript_BMC_Genomics_author_review_v2_zh"
W = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, data: object) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def normalize(text: str) -> str:
    return re.sub(r"[\s*`|]", "", text)


def nums(text: str) -> Counter:
    text = text.replace("October", "10").replace("September", "9")
    return Counter(re.findall(r"(?<![A-Za-z])[+−-]?\d+(?:,\d{3})*(?:\.\d+)?", text))


def command(args: list[str]) -> bytes:
    return subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--content-only", action="store_true")
    parser.add_argument("--visual-confirmed", action="store_true")
    args = parser.parse_args()
    en_text = EN.read_text(encoding="utf-8")
    zh_path = OUT / f"{STEM}.md"
    zh_text = zh_path.read_text(encoding="utf-8")
    en = [(i, s) for i, s in enumerate(en_text.splitlines(), 1) if s.strip()]
    zh = [(i, s) for i, s in enumerate(zh_text.splitlines(), 1) if s.strip()]
    assert len(en) == len(zh) == 197
    block_checks = []
    for index, ((ei, e), (zi, z)) in enumerate(zip(en, zh), 1):
        nc, zc = nums(e), nums(z)
        citations = re.findall(r"\[[\d,]+\]", e) == re.findall(r"\[[\d,]+\]", z)
        heading = len(e) - len(e.lstrip("#")) == len(z) - len(z.lstrip("#"))
        assert nc == zc, (index, "numeric mismatch", nc - zc, zc - nc)
        assert citations and heading, (index, "citation or structural mismatch")
        block_checks.append({"block": index, "en_line": ei, "zh_line": zi,
                             "en_sha256": hashlib.sha256(e.encode()).hexdigest(),
                             "zh_sha256": hashlib.sha256(z.encode()).hexdigest(),
                             "numeric_tokens_preserved": True, "citations_preserved": True})
    refs_en = en_text.split("## References\n", 1)[1].split("## Additional files", 1)[0].strip()
    refs_zh = zh_text.split("## 参考文献\n", 1)[1].split("## 补充文件", 1)[0].strip()
    assert refs_en == refs_zh and len(refs_en.splitlines()) == 31
    protected = json.loads((DOCS / "preservation_manifest.json").read_text(encoding="utf-8"))["files"]
    assert all(sha(ROOT / name) == value for name, value in protected.items())
    with zipfile.ZipFile(OUT / f"{STEM}.docx") as package:
        xml = ET.fromstring(package.read("word/document.xml"))
        body = normalize("".join(t.text or "" for t in xml.findall(".//w:t", W)))
        paragraphs_checked = 0
        for _, line in zh:
            if line in ["## 标题页", "## 图注"] or re.fullmatch(r"\|[\s|:-]+", line):
                continue
            want = re.sub(r"^#{1,3} ", "", line)
            if line.startswith("|"):
                want = want.replace(" to ", " 至 ")
            if re.match(r"^表[12]\. ", want):
                stop = want.index("。") + 1
                assert normalize(want[:stop]) in body and normalize(want[stop:]) in body
            else:
                assert normalize(want) in body, ("Word text missing", paragraphs_checked + 1)
            paragraphs_checked += 1
        tables = xml.findall(".//w:tbl", W)
        assert len(tables) == 2
        images = [package.read(n) for n in package.namelist() if n.startswith("word/media/") and n.endswith(".png")]
    expected_images = [p.read_bytes() for p in (ROOT / "figures/manuscript_m08_v1").glob("Figure*.png") if "_proof" not in p.name]
    assert len(images) == 4
    assert sorted(hashlib.sha256(b).hexdigest() for b in images) == sorted(hashlib.sha256(b).hexdigest() for b in expected_images)
    content = {"date": "2026-10-26", "scope": "Direct translation fidelity and document layout; not a new scientific review",
               "english_source_sha256": sha(EN), "chinese_md_sha256": sha(zh_path),
               "chinese_docx_sha256": sha(OUT / f"{STEM}.docx"),
               "aligned_nonempty_blocks": 197, "all_block_numeric_and_citation_checks_pass": True,
               "references_preserved_exactly": 31, "word_content_blocks_checked": paragraphs_checked,
               "native_tables": 2, "original_pngs_byte_identical": 4,
               "english_files_unchanged": len(protected), "new_fit_score_resampling_or_CTRP_access": False,
               "manual_translation_review": "Completed by the translating agent; not an independent bilingual review",
               "builder_sha256": sha(ROOT / "manuscript_v1/build_chinese_reading_v1.js"),
               "word_render_script_sha256": sha(ROOT / "scripts/render_chinese_reading_v1.vbs"),
               "audit_script_sha256": sha(Path(__file__)),
               "blocks": block_checks}
    dump(DOCS / "content_summary.json", content)
    if args.content_only:
        print(json.dumps({k: v for k, v in content.items() if k != "blocks"}, ensure_ascii=False))
        return
    pdf = OUT / f"{STEM}.pdf"
    qa = OUT / "qa"
    qa.mkdir(exist_ok=True)
    command(["pdftotext", "-layout", str(pdf), str(qa / "main_zh.txt")])
    pages_text = (qa / "main_zh.txt").read_text(encoding="utf-8").split("\f")
    pages_text = [p for p in pages_text if p.strip()]
    n = len(pages_text)
    footer_ok = all(re.search(rf"第\s*{i}\s*页\s*/\s*共\s*{n}\s*页", p) for i, p in enumerate(pages_text, 1))
    fonts = command(["pdffonts", str(pdf)]).decode(errors="replace").splitlines()[2:]
    embedded = bool(fonts) and all(re.search(r"\byes\s+(?:yes|no)\s+(?:yes|no)\s+\d+\s+\d+\s*$", s) for s in fonts if s.strip())
    figures, outside = [], 0
    with pdfplumber.open(pdf) as document:
        for i, p in enumerate(document.pages, 1):
            outside += sum(c["x0"] < 0 or c["x1"] > p.width + .5 or c["top"] < 0 or c["bottom"] > p.height + .5 for c in p.chars)
            for im in p.images:
                figures.append({"page": i, "width_mm": im["width"] / (72 / 25.4),
                                "within_page": im["x0"] >= 0 and im["x1"] <= p.width and im["top"] >= 0 and im["bottom"] <= p.height})
    assert footer_ok and embedded and outside == 0 and len(figures) == 4 and all(f["within_page"] for f in figures)
    assert pdf.stat().st_mtime >= (OUT / f"{STEM}.docx").stat().st_mtime
    thumbs = qa / "thumbs"
    thumbs.mkdir(exist_ok=True)
    command(["pdftoppm", "-scale-to", "500", "-png", str(pdf), str(thumbs / "page")])
    paths = sorted(thumbs.glob("page-*.png"))
    assert len(paths) == n
    for start in range(0, n, 12):
        sheet = Image.new("RGB", (1200, 2110), "#cccccc")
        draw = ImageDraw.Draw(sheet)
        for slot, p in enumerate(paths[start:start + 12]):
            im = Image.open(p).convert("RGB")
            im.thumbnail((394, 500))
            x, y = (slot % 3) * 400, (slot // 3) * 525
            draw.text((x + 8, y + 3), f"CN page {start + slot + 1}", fill="black")
            sheet.paste(im, (x + (400 - im.width) // 2, y + 20))
        sheet.save(qa / f"contact_{start // 12 + 1}.png")
    layout = {"date": "2026-10-26", "pdf_sha256": sha(pdf), "docx_sha256": sha(OUT / f"{STEM}.docx"),
              "pages": n, "consecutive_chinese_page_footers": footer_ok,
              "all_fonts_embedded": embedded, "text_outside_page": outside, "figures": figures,
              "mechanical_layout_checks_pass": True, "manual_visual_inspection_confirmed": args.visual_confirmed}
    dump(DOCS / "render_summary.json", layout)
    print(json.dumps(layout, ensure_ascii=False))


if __name__ == "__main__":
    main()
