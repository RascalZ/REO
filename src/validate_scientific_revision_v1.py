"""Coverage-bounded editorial/source audit; never load raw data or model payloads."""
from __future__ import annotations
import csv
import hashlib
import json
import re
import subprocess
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from validate_numeric_references import expand_numeric_group

ROOT = Path(__file__).resolve().parents[1]
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    checks = {}
    details = {}
    path = ROOT / "manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md"
    md = path.read_text(encoding="utf-8")
    body = md.split("## References", 1)[0]
    first = []
    for group in re.findall(r"\[([^\]]+)\]", body):
        expanded = expand_numeric_group(group, 30)
        if expanded:
            for n in sorted(expanded):
                if n not in first:
                    first.append(n)
    checks["reference_first_appearance_1_to_30"] = first == list(range(1, 31))
    citation = json.loads((ROOT / "docs/scientific_revision_v1/citation_validation.json").read_text(encoding="utf-8"))
    checks["numeric_citation_and_bibliography_audit"] = citation["pass"]
    old = ["manuscript_v1/Manuscript_draft.md", "manuscript_v1/Manuscript_BMC_Genomics_draft.md", "manuscript_v1/Manuscript_BMC_Genomics_draft.docx"]
    r = subprocess.run(["git", "diff", "--exit-code", "HEAD", "--", *old], cwd=ROOT, capture_output=True)
    checks["old_manuscripts_unchanged_from_base_HEAD"] = r.returncode == 0
    protected = ["config/", "results/", "figures/", "src/ordinal_repeated_evaluation.py", "src/external_validation_v1.py", "src/run_external_validation_v1.py", "src/synthetic_mechanisms.py", "src/run_synthetic_mechanisms.py", "src/analyze_synthetic_mechanisms.py"]
    changed = subprocess.check_output(["git", "diff", "--name-only", "HEAD", "--", *protected], cwd=ROOT).decode("utf-8").splitlines()
    unchanged = True
    for name in changed:
        # R1 repairs Git storage of existing CRLF aggregates; local recorded
        # bytes are not edited. Only the explicitly documented result namespace
        # may differ from its old Git blob by CRLF/LF alone.
        if not name.startswith("results/synthetic_mechanisms/"):
            unchanged = False
            break
        old_blob = subprocess.check_output(["git", "show", "HEAD:" + name], cwd=ROOT)
        if old_blob.replace(b"\r\n", b"\n") != (ROOT / name).read_bytes().replace(b"\r\n", b"\n"):
            unchanged = False
            break
    checks["protected_scientific_files_unchanged_except_documented_git_EOL_storage_repair"] = unchanged
    ledger = list(csv.DictReader((ROOT / "docs/scientific_revision_v1/methods_source_ledger.csv").read_text(encoding="utf-8-sig").splitlines()))
    failed_sources = [row["id"] for row in ledger if sha(ROOT / row["source_path"]) != row["source_sha256"]]
    checks["methods_ledger_source_hashes"] = not failed_sources
    details["methods_source_rows"] = len(ledger)
    table_manifest = json.loads((ROOT / "tables/scientific_revision_v1/revision_table_manifest.json").read_text(encoding="utf-8"))
    checks["aggregate_table_source_and_output_hashes"] = all(sha(ROOT / p) == v["sha256"] for section in ["sources", "outputs"] for p, v in table_manifest[section].items())
    synth_manifest = json.loads((ROOT / "tables/scientific_revision_v1/synthetic_table_manifest.json").read_text(encoding="utf-8"))
    checks["synthetic_copies_unchanged"] = all((ROOT / row["path"]).read_bytes() == (ROOT / row["source_path"]).read_bytes() for row in synth_manifest["tables"])
    package = json.loads((ROOT / "manuscript_v1/bmc_additional_files_scientific_revision_v1/package_manifest.json").read_text(encoding="utf-8"))
    checks["all_14_additional_files_match_sources"] = len(package["files"]) == 14 and all(sha(ROOT / row["path"]) == sha(ROOT / row["source"]) == row["sha256"] for row in package["files"])
    checks["package_document_hashes"] = all(sha(ROOT / row["path"]) == row["sha256"] for row in package["documents"])
    checks["package_docx_builder_hash"] = sha(ROOT / "manuscript_v1/build_bmc_docx.js") == package["docx_builder_sha256"]
    additional_citations = {int(n) for n in re.findall(r"\*\*Additional file (\d+)\.", md)}
    checks["main_additional_file_list_1_to_14"] = additional_citations == set(range(1, 15))
    details["docx"] = {}
    for kind, name in [("main", "Manuscript_BMC_Genomics_scientific_revision_v1.docx"), ("supplement", "Supplementary_methods_scientific_revision_v1.docx")]:
        with zipfile.ZipFile(ROOT / "manuscript_v1" / name) as z:
            tree = ET.fromstring(z.read("word/document.xml"))
            paragraphs = ["".join(t.text or "" for t in p.findall(".//w:t", NS)) for p in tree.findall(".//w:p", NS)]
            drawings = len(tree.findall(".//w:drawing", NS))
            tables = tree.findall(".//w:tbl", NS)
            cols = [len(t.find("w:tr", NS).findall("w:tc", NS)) for t in tables]
            details["docx"][kind] = {"sha256": sha(ROOT / "manuscript_v1" / name), "drawings": drawings, "table_columns": cols, "paragraphs": len(paragraphs), "zip_integrity": z.testzip() is None}
            checks[kind + "_docx_zip_integrity"] = z.testzip() is None
            footer = ET.fromstring(z.read("word/footer1.xml"))
            checks[kind + "_footer_has_current_page_field"] = any("PAGE" in (n.text or "") for n in footer.findall(".//w:instrText", NS)) and any((n.text or "").startswith("Page ") for n in footer.findall(".//w:t", NS))
            checks[kind + "_docx_draft_status_visible"] = any("NOT FOR SUBMISSION" in p for p in paragraphs)
            if kind == "main":
                checks["main_4_images_2_editable_tables"] = drawings == 4 and cols == [5, 5]
                checks["main_one_caption_per_image"] = all(sum(p.startswith(f"Figure {i}.") for p in paragraphs) == 1 for i in range(1, 5))
                checks["main_30_reference_rows"] = len([p for p in paragraphs if re.match(r"^\d+\. [A-Z]", p)]) == 30
            else:
                checks["supplement_two_tables_correct_columns"] = cols == [2, 3]
                checks["supplement_schema_and_synthetic_sections_present"] = all(any(x in p for p in paragraphs) for x in ["SM4.1.", "SM7.1.", "SM8.6."])
    supp = (ROOT / "manuscript_v1/Supplementary_methods_scientific_revision_v1.md").read_text(encoding="utf-8")
    checks["supplement_no_unrendered_latex_or_reserved_marker"] = not any(x in supp for x in [r"\(", r"\)", "RESERVED_FOR_ROOT"])
    abstract = md.split("## Abstract", 1)[1].split("## Keywords", 1)[0]
    details["abstract_words_whitespace_count"] = len(abstract.split())
    checks["abstract_at_most_350_words"] = len(abstract.split()) <= 350
    out = {"date": "2026-10-16", "scope": "editorial aggregate source checks and package structural validation; not independent raw-data/model reproduction", "no_model_fit_resampling_external_rescoring": True, "checks": checks, "details": details, "failed_checks": [k for k, v in checks.items() if not v], "pass": all(checks.values())}
    target = ROOT / "docs/scientific_revision_v1/revision_validation.json"
    target.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"checks": len(checks), "pass": out["pass"], "failed": out["failed_checks"], "abstract_words": len(abstract.split())}))
    if not out["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
