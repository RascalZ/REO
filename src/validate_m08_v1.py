"""Bounded text/package/figure audit for M08; never access individual data."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from format_m08_reference_authors import apply_author_replacements
from validate_numeric_references import expand_numeric_group

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/m08_v1"
MS = ROOT / "manuscript_v1"
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def plain(text: str) -> str:
    return re.sub(r"\\([*_])", r"\1", text).replace("*", "").replace("`", "")


def text(node: ET.Element) -> str:
    return "".join(t.text or "" for t in node.findall(".//w:t", NS))


def main() -> None:
    old = (MS / "Manuscript_BMC_Genomics_scientific_revision_v1.md").read_text(encoding="utf-8")
    md = (MS / "Manuscript_BMC_Genomics_M08_v1.md").read_text(encoding="utf-8")
    audit = json.loads((DOCS / "reference_author_audit.json").read_text(encoding="utf-8"))
    strip_status = lambda s: re.sub(r"^\*\*Draft status:\*\* .+$", "", s, flags=re.M)
    checks = {"only_status_and_reference_author_prefixes_changed":
              strip_status(apply_author_replacements(old, audit)) == strip_status(md)}
    preserved = json.loads((DOCS / "preserved_sources.json").read_text(encoding="utf-8"))
    checks["all_37_preserved_sources_byte_identical"] = len(preserved["files"]) == 37 and all(
        sha(ROOT / p) == h for p, h in preserved["files"].items())
    checks["reference_metadata_sources_unchanged"] = all(sha(ROOT / p) == h for p, h in audit["source_sha256"].items())
    first = []
    for group in re.findall(r"\[([^\]]+)\]", md.split("## References")[0]):
        for n in sorted(expand_numeric_group(group, 30) or set()):
            if n not in first:
                first.append(n)
    checks["first_citations_are_1_to_30"] = first == list(range(1, 31))
    citation = json.loads((DOCS / "citation_validation.json").read_text(encoding="utf-8"))
    checks["citation_bibtex_coverage_pass"] = citation["pass"]
    captions = re.findall(r"^Figure ([1-4])\. (.+)$", md, re.M)
    checks["four_caption_titles_and_legends_within_limits"] = len(captions) == 4 and all(
        len(c.split(". ", 1)[0].split()) <= 15 and len(c.split(". ", 1)[1].split()) <= 300 for _, c in captions)
    expected = []
    for line in md.splitlines():
        if not line.strip() or line.startswith(("# BMC ", "## Title page", "## Figure legends", "|")):
            continue
        table_caption = re.match(r"^(Table [12]\. [^.]+\.) (.+)$", line)
        if table_caption:
            expected.extend(plain(s) for s in table_caption.groups())
        else:
            expected.append(plain(re.sub(r"^#{1,3} ", "", line)))
    expected.append("DRAFT — AUTHOR INPUT REQUIRED — NOT FOR SUBMISSION")
    # Figures move next to their Results anchors; compare the full paragraph multiset.
    with zipfile.ZipFile(MS / "Manuscript_BMC_Genomics_M08_v1.docx") as z:
        tree = ET.fromstring(z.read("word/document.xml"))
        body = tree.find("w:body", NS)
        paragraphs = [text(p) for p in body.findall("w:p", NS) if text(p)]
        checks["word_paragraphs_equal_markdown_including_references_and_captions"] = Counter(paragraphs) == Counter(expected)
        checks["one_caption_per_figure"] = all(sum(p.startswith(f"Figure {i}.") for p in paragraphs) == 1 for i in range(1, 5))
        tables = body.findall("w:tbl", NS)
        word_cells = [[[text(c) for c in row.findall("w:tc", NS)] for row in table.findall("w:tr", NS)] for table in tables]
        expected_tables = []
        for block in re.findall(r"(?m)(?:^\|.*\n)+", md):
            rows = [[c.strip() for c in line.strip().strip("|").split("|")] for line in block.splitlines()]
            expected_tables.append([[plain(c) for c in row] for row in rows if not all(re.fullmatch(r":?-+:?", c) for c in row)])
        checks["two_editable_tables_cells_equal_markdown"] = word_cells == expected_tables and len(tables) == 2
        checks["four_drawings"] = len(tree.findall(".//w:drawing", NS)) == 4
        embedded = [hashlib.sha256(z.read(p)).hexdigest() for p in z.namelist()
                    if p.startswith("word/media/") and not p.endswith("/")]
        figdir = ROOT / "figures/manuscript_m08_v1"
        image_names = ["Figure1_study_design", "Figure2_development_primary",
                       "Figure3_representation_and_zero_solutions", "Figure4_loto_and_external"]
        checks["word_images_are_exact_new_png_bytes"] = Counter(embedded) == Counter(sha(figdir / (n + ".png")) for n in image_names)
        checks["docx_zip_integrity"] = z.testzip() is None
        footer = ET.fromstring(z.read("word/footer1.xml"))
        checks["current_page_field_present"] = any("PAGE" in (n.text or "") for n in footer.findall(".//w:instrText", NS))
        checks["line_numbering_present"] = bool(tree.findall(".//w:lnNumType", NS))
    package = json.loads((MS / "bmc_additional_files_M08_v1/package_manifest.json").read_text(encoding="utf-8"))
    checks["14_additional_files_byte_identical_to_R1"] = len(package["files"]) == 14 and all(
        sha(ROOT / r["path"]) == sha(ROOT / r["source"]) == r["sha256"] for r in package["files"])
    checks["package_document_source_and_builder_hashes_current"] = all(sha(ROOT / r["path"]) == r["sha256"]
        for r in package["documents"] + [package[k] for k in ["generator", "builder", "source_manuscript", "figure_manifest", "author_replacements"]])
    fig = json.loads((ROOT / "figures/manuscript_m08_v1/figure_manifest.json").read_text(encoding="utf-8"))
    checks["figure_source_hashes_current"] = all(sha(ROOT / p) == h for p, h in fig["source_sha256"].items())
    checks["all_four_numeric_artists_match_original"] = len(fig["figures"]) == 4 and all(x["same_numerical_artists_as_original"] for x in fig["figures"].values())
    checks["all_four_word_minimum_fonts_at_least_8pt"] = all(x["minimum_font_pt_at_word_590px"] >= 8 for x in fig["figures"].values())
    out = {"scope": "M08 text, figures and package consistency only; not independent study reproduction or submission readiness",
           "checks": checks, "pass": all(checks.values()), "failed_checks": [k for k, v in checks.items() if not v],
           "details": {"paragraphs": len(paragraphs), "table_rows": [len(t) for t in word_cells],
                       "reference_count": 30, "preserved_files": len(preserved["files"])}}
    (DOCS / "validation.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"checks": len(checks), "pass": out["pass"], "failed": out["failed_checks"]}))
    if not out["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
