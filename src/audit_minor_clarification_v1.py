"""Check the bounded bilingual edit and render its PDFs for local visual review."""
from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import pdfplumber
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "docs/minor_clarification_v1"
PPM = (Path.home() / ".cache/codex-runtimes/codex-primary-runtime/dependencies/native/"
       "poppler/Library/bin/pdftoppm.exe")
W = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(text: str) -> str:
    return re.sub(r"[\s*`|]", "", text)


def numeric_tokens(text: str) -> Counter:
    text = text.replace("October", "10").replace("September", "9")
    return Counter(re.findall(r"(?<![A-Za-z])[+−-]?\d+(?:,\d{3})*(?:\.\d+)?", text))


def section(text: str, name: str) -> str:
    return text.split(f"## {name}\n", 1)[1].split("\n## ", 1)[0]


def render(info: dict) -> tuple[str, int]:
    md = ROOT / info["output"]
    pages = md.parent / "qa/pages"
    pages.mkdir(parents=True, exist_ok=True)
    assert not list(pages.glob("page-*.png")), "Preserve existing render images"
    subprocess.run([str(PPM), "-scale-to", "1600", "-png", str(md.with_suffix(".pdf")),
                    str(pages / "page")], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return info["output"], len(list(pages.glob("page-*.png")))


def main() -> None:
    build = json.loads((AUDIT / "build_manifest.json").read_text(encoding="utf-8"))
    protection = json.loads((AUDIT / "preservation_manifest.json").read_text(encoding="utf-8"))
    assert all(sha(ROOT / name) == value for name, value in protection["files"].items())
    summary = {"date": "2026-10-26", "scope": "Two clarifications; no scientific recomputation",
               "preserved_files": len(protection["files"]), "languages": {},
               "visual_review_status": "pending", "individual_CTRP_access": False,
               "new_training_resampling_scoring": False,
               "renderer": "Native Microsoft Word read-only PDF export, bundled Poppler page PNGs",
               "skill_renderer_attempt": "render_docx.py failed: no LibreOffice soffice.exe; native Word fallback used",
               "preexisting_historical_differences": protection["preexisting_historical_differences"],
               "audit_sha256": sha(Path(__file__)),
               "render_helper_sha256": sha(ROOT / "scripts/render_minor_clarification_v1.vbs")}
    paired = {}
    for lang, info in build["languages"].items():
        old, new = ROOT / info["source"], ROOT / info["output"]
        assert sha(new) == info["md_sha256"] and sha(new.with_suffix(".docx")) == info["docx_sha256"]
        a, b = old.read_text(encoding="utf-8"), new.read_text(encoding="utf-8")
        undone = b
        for edit in info["edits"]:
            assert undone.count(edit["new_block"]) == 1
            undone = undone.replace(edit["new_block"], edit["old_block"], 1)
        assert undone == a
        assert [line for line in a.splitlines() if line.startswith("|")] == [
            line for line in b.splitlines() if line.startswith("|")]
        for name in (["Abstract", "Results", "Discussion", "Conclusions", "Declarations", "References"]
                     if lang == "en" else ["摘要", "结果", "讨论", "结论", "声明", "参考文献"]):
            assert section(a, name) == section(b, name), (lang, name)
        with ZipFile(old.with_suffix(".docx")) as zo, ZipFile(new.with_suffix(".docx")) as zn:
            assert zo.namelist() == zn.namelist()
            changed = [name for name in zo.namelist() if zo.read(name) != zn.read(name)]
            assert changed == info["changed_ooxml_members"]
            media = [name for name in zo.namelist() if name.startswith("word/media/")]
            assert all(zo.read(name) == zn.read(name) for name in media)
            xo, xn = [ET.fromstring(z.read("word/document.xml")) for z in [zo, zn]]
            assert [ET.tostring(x) for x in xo.findall(".//w:tbl", W)] == [
                ET.tostring(x) for x in xn.findall(".//w:tbl", W)]
            assert len(xn.findall(".//w:tbl", W)) == 2
            word_text = "".join(t.text or "" for t in xn.findall(".//w:t", W))
            assert all(normalize(edit["new_block"]) in normalize(word_text) for edit in info["edits"])
        pdf = new.with_suffix(".pdf")
        reader = PdfReader(pdf)
        body_pages = []
        for page in reader.pages:
            lines = page.extract_text().splitlines()
            lines = [line for line in lines if not re.fullmatch(r"\s*Page \d+\s*", line)
                     and not re.fullmatch(r"\s*第\s*\d+\s*页\s*/\s*共\s*\d+\s*页\s*", line)
                     and not line.startswith("中文修订通读")]
            body_pages.append("\n".join(lines))
        if lang == "en":
            # Native Word line numbers occupy the left margin but pypdf may
            # append them to body lines. Exclude that margin geometrically.
            with pdfplumber.open(pdf) as document:
                body_pages = [page.crop((70, 45, page.width - 50, page.height - 40))
                              .extract_text() or "" for page in document.pages]
        pdf_text = normalize("\n".join(body_pages))
        assert all(normalize(edit["addition"]) in pdf_text for edit in info["edits"]), lang
        out_of_page = 0
        images = []
        with pdfplumber.open(pdf) as document:
            for i, page in enumerate(document.pages, 1):
                out_of_page += sum(c["x0"] < -.5 or c["x1"] > page.width + .5 or
                                   c["top"] < -.5 or c["bottom"] > page.height + .5 for c in page.chars)
                for image in page.images:
                    assert image["x0"] >= 0 and image["x1"] <= page.width + .5
                    assert image["top"] >= 0 and image["bottom"] <= page.height + .5
                    images.append(i)
        assert out_of_page == 0 and len(images) == 4
        paired[lang] = [line for line in b.splitlines() if line.strip()]
        summary["languages"][lang] = {
            "changed_blocks": 2, "unchanged_nonempty_blocks": info["unchanged_nonempty_blocks"],
            "references_abstract_results_discussion_declarations_unchanged": True,
            "tables_identical": 2, "all_media_identical": len(media), "figure_count": len(images),
            "figure_pages": images, "pdf_pages": len(reader.pages), "text_outside_page": 0,
            "new_text_present_in_Word_and_PDF": True, "pdf_sha256": sha(pdf),
        }
    assert len(paired["en"]) == len(paired["zh"]) == 197
    for i, (en, zh) in enumerate(zip(paired["en"], paired["zh"]), 1):
        assert numeric_tokens(en) == numeric_tokens(zh), (i, "numeric mismatch")
        assert re.findall(r"\[[\d,]+\]", en) == re.findall(r"\[[\d,]+\]", zh)
    with ThreadPoolExecutor(max_workers=2) as executor:
        rendered = list(executor.map(render, build["languages"].values()))
    for output, count in rendered:
        lang = next(k for k, v in build["languages"].items() if v["output"] == output)
        assert count == summary["languages"][lang]["pdf_pages"]
        summary["languages"][lang]["page_pngs"] = count
    summary["paired_numeric_citation_blocks_checked"] = 197
    summary["mechanical_checks_passed"] = True
    (AUDIT / "mechanical_audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
                                                encoding="utf-8", newline="\n")
    print(json.dumps({"mechanical_checks": "PASS", "preserved_files": summary["preserved_files"],
                      "pages": {k: v["pdf_pages"] for k, v in summary["languages"].items()},
                      "visual_review": "pending"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
