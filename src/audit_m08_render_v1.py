"""Inspect rendered M08 documents and make physical-width proofs, no study data."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
from importlib.metadata import version
from pathlib import Path

import pdfplumber
from PIL import Image, ImageDraw
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph

ROOT = Path(__file__).resolve().parents[1]
QA = ROOT / "manuscript_v1/qa_m08_v1"
DOCS = ROOT / "docs/m08_v1"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(command: list[str]) -> bytes:
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--visual-confirmed", action="store_true",
                        help="Set only after a human/agent actually views the produced pages")
    args = parser.parse_args()
    captions = json.loads((DOCS / "figure_captions.json").read_text(encoding="utf-8"))
    names = ["Figure1_study_design", "Figure2_development_primary",
             "Figure3_representation_and_zero_solutions", "Figure4_loto_and_external"]
    proof = QA / "figure_caption_170mm_proof.pdf"
    cv = canvas.Canvas(str(proof), pagesize=(210 * mm, 297 * mm))
    style = ParagraphStyle("caption", fontName="Times-Roman", fontSize=10, leading=12.5)
    dimensions = []
    for name, caption in zip(names, captions):
        path = ROOT / "figures/manuscript_m08_v1" / f"{name}.png"
        iw, ih = Image.open(path).size
        width, height = 170 * mm, 170 * mm * ih / iw
        paragraph = Paragraph(f"<b>Figure {caption['number']}. {caption['title']}.</b> {caption['legend']}", style)
        _, ch = paragraph.wrap(width, 400 * mm)
        cv.drawImage(str(path), 20 * mm, 277 * mm - height, width=width, height=height)
        paragraph.drawOn(cv, 20 * mm, 277 * mm - height - 3 * mm - ch)
        cv.showPage()
        dimensions.append({"figure": caption["number"], "image_height_mm": height / mm,
                           "caption_height_mm": ch / mm, "gap_mm": 3,
                           "combined_height_mm": (height + ch) / mm + 3,
                           "within_225mm": (height + ch) / mm + 3 <= 225,
                           "caption_font_pt": 10, "caption_leading_pt": 12.5})
    cv.save()
    dim_record = {"proof_width_mm": 170, "dimensions": dimensions,
                  "all_within_225mm": all(x["within_225mm"] for x in dimensions),
                  "scope": "Exact proof layout only; publisher re-typesetting may require further adjustment"}
    (DOCS / "figure_caption_dimensions.json").write_text(json.dumps(dim_record, indent=2) + "\n", encoding="utf-8")
    output = {"scope": "Microsoft Word render and layout QA, not new scientific evaluation",
              "render_engine": "Microsoft Word 2021 COM ExportAsFixedFormat",
              "editorial_runtime_only": {"python": platform.python_version(),
                  **{p: version(p) for p in ["numpy", "pandas", "matplotlib", "Pillow", "pdfplumber", "reportlab"]}},
              "render_script_sha256": sha(ROOT / "scripts/render_m08_word.ps1"),
              "audit_script_sha256": sha(Path(__file__)),
              "manual_visual_inspection_confirmed": args.visual_confirmed,
              "local_full_width_proof": proof.relative_to(ROOT).as_posix(),
              "full_width_proof_sha256": sha(proof),
              "full_width_proof": dim_record, "documents": {}}
    for kind, docx in [("main", ROOT / "manuscript_v1/Manuscript_BMC_Genomics_M08_v1.docx"),
                       ("supplement", ROOT / "manuscript_v1/bmc_additional_files_M08_v1/Additional_file_14.docx")]:
        pdf = QA / f"{kind}_M08.pdf"
        target = QA / f"{kind}_M08.txt"
        run(["pdftotext", "-layout", str(pdf), str(target)])
        pages_text = [p for p in target.read_text(encoding="utf-8").split("\f") if p.strip()]
        footers = [re.findall(r"^\s*Page (\d+)\s*$", p, re.M) for p in pages_text]
        document = {"docx_path": docx.relative_to(ROOT).as_posix(), "docx_sha256": sha(docx),
                    "local_pdf": pdf.relative_to(ROOT).as_posix(), "pdf_sha256": sha(pdf),
                    "pages": len(pages_text), "sequential_page_numbers": footers == [[str(i)] for i in range(1, len(pages_text) + 1)],
                    "no_blank_pages": all(len(p.strip()) > 30 for p in pages_text),
                    "fonts_output": run(["pdffonts", str(pdf)]).decode("utf-8", errors="replace"),
                    "figure_locations": [], "text_outside_page": 0}
        with pdfplumber.open(pdf) as rendered:
            for n, page in enumerate(rendered.pages, 1):
                document["text_outside_page"] += sum(c["x0"] < 0 or c["x1"] > page.width + .5 or c["top"] < 0 or c["bottom"] > page.height + .5 for c in page.chars)
                for im in page.images:
                    caption_chars = [c for c in page.chars if c["top"] >= im["bottom"] and 9.7 <= c["size"] <= 10.3]
                    caption_bottom = max(c["bottom"] for c in caption_chars)
                    document["figure_locations"].append({"page": n, "image_width_mm": im["width"] / mm,
                        "image_height_mm": im["height"] / mm, "combined_height_mm": (caption_bottom - im["top"]) / mm,
                        "within_225mm": (caption_bottom - im["top"]) / mm <= 225,
                        "image_within_body": im["x0"] >= 71.5 and im["x1"] <= page.width - 71.5 and im["top"] >= 71.5 and im["bottom"] <= page.height - 71.5})
        font_rows = document["fonts_output"].splitlines()[2:]
        document["all_used_word_pdf_fonts_embedded"] = bool(font_rows) and all(re.search(r"\byes\s+(?:yes|no)\s+(?:yes|no)\s+\d+\s+\d+\s*$", row) for row in font_rows if row.strip())
        thumb_dir = QA / f"{kind}_thumbs"
        thumb_dir.mkdir(exist_ok=True)
        run(["pdftoppm", "-scale-to", "420", "-png", str(pdf), str(thumb_dir / "page")])
        pages = sorted(thumb_dir.glob("page-*.png"))
        for start in range(0, len(pages), 16):
            sheet = Image.new("RGB", (1200, 1780), "#cccccc")
            draw = ImageDraw.Draw(sheet)
            for slot, path in enumerate(pages[start:start+16]):
                im = Image.open(path).convert("RGB")
                im.thumbnail((294, 420))
                sheet.paste(im, ((slot % 4) * 300 + (300 - im.width) // 2, (slot // 4) * 445 + 20))
                draw.text(((slot % 4) * 300 + 5, (slot // 4) * 445 + 3), f"{kind} page {start + slot + 1}", fill="black")
            sheet.save(QA / f"{kind}_contact_{start // 16 + 1}.png")
        output["documents"][kind] = document
    output["mechanical_checks_pass"] = (dim_record["all_within_225mm"] and
        all(d["sequential_page_numbers"] and d["no_blank_pages"] and d["all_used_word_pdf_fonts_embedded"]
            and d["text_outside_page"] == 0 and all(f["within_225mm"] and f["image_within_body"] for f in d["figure_locations"])
            for d in output["documents"].values()) and len(output["documents"]["main"]["figure_locations"]) == 4)
    (DOCS / "render_validation.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"pass": output["mechanical_checks_pass"], "manual_visual_confirmed": args.visual_confirmed,
                      "pages": {k: v["pages"] for k, v in output["documents"].items()}}))
    if not output["mechanical_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
