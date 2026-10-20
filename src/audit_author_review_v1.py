"""Layout/source audit of the LOCAL author reading copy, without study inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import pdfplumber
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript_v1/author_review_v1"
QA = OUT / "qa"
DOCS = ROOT / "docs/author_review_v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args: list[str]) -> bytes:
    return subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout


def dump(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--visual-confirmed", action="store_true")
    args = parser.parse_args()
    QA.mkdir(exist_ok=True)
    preserved = json.loads((DOCS / "local_private/preserved_sources.json").read_text(encoding="utf-8"))["files"]
    unchanged = all(sha(ROOT / path) == value for path, value in preserved.items())
    evidence = json.loads((DOCS / "local_private/source_author_declarations.json").read_text(encoding="utf-8"))
    legacy_unchanged = sha(Path(evidence["source"]["absolute_path"])) == evidence["source"]["sha256"]
    result = {"date": "2026-10-20", "scope": "Source preservation and rendered editorial layout, not scientific replication",
              "original_source_unchanged": legacy_unchanged, "protected_M08_files": len(preserved),
              "all_protected_M08_files_unchanged": unchanged,
              "manual_visual_inspection_confirmed": args.visual_confirmed,
              "render_script_sha256": sha(ROOT / "scripts/render_author_review_v1.ps1"),
              "word_dispatch_script_sha256": sha(ROOT / "scripts/render_author_review_v1.vbs"),
              "audit_script_sha256": sha(Path(__file__)), "documents": {}}
    for kind, stem in [("main", "Manuscript_BMC_Genomics_author_review_v1"),
                       ("supplement", "Supplementary_methods_author_review_v1")]:
        docx, pdf = OUT / f"{stem}.docx", OUT / f"{stem}.pdf"
        text_path = QA / f"{kind}.txt"
        run(["pdftotext", "-layout", str(pdf), str(text_path)])
        pages = [s for s in text_path.read_text(encoding="utf-8").split("\f") if s.strip()]
        footers = [re.findall(r"^\s*Page (\d+)\s*$", s, re.M) for s in pages]
        font_output = run(["pdffonts", str(pdf)]).decode("utf-8", errors="replace")
        font_rows = font_output.splitlines()[2:]
        d = {"docx_sha256": sha(docx), "pdf_sha256": sha(pdf), "pages": len(pages),
             "render_not_older_than_docx": pdf.stat().st_mtime >= docx.stat().st_mtime,
             "sequential_page_numbers": footers == [[str(i)] for i in range(1, len(pages) + 1)],
             "no_blank_pages": all(len(s.strip()) > 30 for s in pages),
             "all_used_fonts_embedded": bool(font_rows) and all(re.search(r"\byes\s+(?:yes|no)\s+(?:yes|no)\s+\d+\s+\d+\s*$", row) for row in font_rows if row.strip()),
             "text_outside_page": 0, "figure_locations": [], "focus_pages": []}
        if kind == "main":
            without_footers = re.sub(r"(?m)^\s*Page \d+\s*$", "", " ".join(pages))
            without_line_numbers = re.sub(r"(?m)^\s*\d+\s+", "", without_footers)
            normalized = re.sub(r"\s+", " ", without_line_numbers)
            d["current_access_wording_rendered"] = ("verified as private on 20 October 2026" in normalized
                and "Private repository; access status verified 20 October 2026." in normalized
                and "are available in the REO repository" not in normalized)
        with pdfplumber.open(pdf) as rendered:
            for n, page in enumerate(rendered.pages, 1):
                d["text_outside_page"] += sum(c["x0"] < 0 or c["x1"] > page.width + .5 or c["top"] < 0 or c["bottom"] > page.height + .5 for c in page.chars)
                for im in page.images:
                    below = [c["bottom"] for c in page.chars if c["top"] >= im["bottom"] and 9.7 <= c["size"] <= 10.3]
                    bottom = max(below, default=im["bottom"])
                    combined = (bottom - im["top"]) / (72 / 25.4)
                    d["figure_locations"].append({"page": n, "width_mm": im["width"] / (72 / 25.4),
                        "combined_height_mm": combined, "within_225mm": combined <= 225,
                        "image_within_body": im["x0"] >= 71.5 and im["x1"] <= page.width - 71.5 and im["top"] >= 71.5 and im["bottom"] <= page.height - 71.5})
        for n, page_text in enumerate(pages, 1):
            if n == 1 or re.search(r"^(?:\s*\d+\s+)?(?:Declarations|Funding|Competing interests|Code availability|Authors' contributions|References|Figure [1-4]\.|Table [12]\.)", page_text, re.M) or any(f["page"] == n for f in d["figure_locations"]):
                d["focus_pages"].append(n)
        thumb = QA / f"{kind}_thumbs"
        thumb.mkdir(exist_ok=True)
        run(["pdftoppm", "-scale-to", "420", "-png", str(pdf), str(thumb / "page")])
        images = sorted(thumb.glob("page-*.png"))
        for start in range(0, len(images), 16):
            sheet = Image.new("RGB", (1200, 1780), "#cccccc")
            draw = ImageDraw.Draw(sheet)
            for slot, path in enumerate(images[start:start + 16]):
                im = Image.open(path).convert("RGB")
                im.thumbnail((294, 420))
                sheet.paste(im, ((slot % 4) * 300 + (300 - im.width) // 2, (slot // 4) * 445 + 20))
                draw.text(((slot % 4) * 300 + 5, (slot // 4) * 445 + 3), f"{kind} {start + slot + 1}", fill="black")
            sheet.save(QA / f"{kind}_contact_{start // 16 + 1}.png")
        result["documents"][kind] = d
    result["mechanical_checks_pass"] = (unchanged and legacy_unchanged and
        all(d["render_not_older_than_docx"] and d["sequential_page_numbers"] and d["no_blank_pages"] and d["all_used_fonts_embedded"] and d["text_outside_page"] == 0
            and all(f["within_225mm"] and f["image_within_body"] for f in d["figure_locations"]) for d in result["documents"].values())
        and result["documents"]["main"]["current_access_wording_rendered"]
        and len(result["documents"]["main"]["figure_locations"]) == 4)
    dump(OUT / "render_validation.json", result)
    # No author names, manuscript prose, page text, or correspondence is copied into this status record.
    dump(DOCS / "render_summary.json", result)
    print(json.dumps({"pass": result["mechanical_checks_pass"], "visual_confirmed": args.visual_confirmed,
                      "documents": result["documents"]}))
    if not result["mechanical_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
