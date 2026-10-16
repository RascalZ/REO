"""Inspect the submission snapshot without fitting or opening sample responses.

Run from the repository root. The output is a file/format consistency audit,
not a scientific validity score or full computational reproduction.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import re
import struct
import sys
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path.cwd()
OUT = ROOT / "docs/reviews/bmc_pre_submission_v1"
sys.path.insert(0, str(ROOT / "src"))
from validate_numeric_references import expand_numeric_group  # noqa: E402


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def inspect_markdown(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    refs = [int(n) for n in re.findall(r"^(\d+)\. ", text, re.M)]
    body = re.split(r"## (?:\d+\.\s*)?(?:Background|Introduction)\n", text, maxsplit=1)[1].split("## References\n", 1)[0]
    order = []
    for token in re.findall(r"\[([^\]]+)\]", body):
        for number in sorted(expand_numeric_group(token, len(refs)) or []):
            if number not in order:
                order.append(number)
    captions = []
    for number, caption in re.findall(r"^Figure (\d)\. (.+)$", text, re.M):
        captions.append({"number": int(number),
                         "title_words": len(caption.split(". ", 1)[0].split()),
                         "legend_words_including_title": len(caption.split())})
    return {"sha256": sha(path), "reference_count": len(refs),
            "citations_without_reference": sorted(set(order) - set(refs)),
            "first_appearance_order": order,
            "first_appearance_order_is_sequential": order == list(range(1, len(refs) + 1)),
            "author_placeholder_lines": [i for i, line in enumerate(text.splitlines(), 1)
                                         if "AUTHOR INPUT REQUIRED" in line],
            "figure_captions": captions}


def main() -> None:
    manuscript = ROOT / "manuscript_v1"
    tables = ROOT / "tables/manuscript_v1"
    manifests = {}
    for kind in ["figures", "tables"]:
        folder = ROOT / kind / "manuscript_v1"
        manifest = json.loads((folder / ("figure_manifest.json" if kind == "figures"
                                         else "table_manifest.json")).read_text(encoding="utf-8"))
        failures = []
        for file, digest in manifest["source_sha256"].items():
            if sha(ROOT / file) != digest:
                failures.append(file)
        for file, digest in manifest["output_sha256"].items():
            if sha(folder / file) != digest:
                failures.append(file)
        if sha(ROOT / manifest["script"]) != manifest["script_sha256"]:
            failures.append(manifest["script"])
        manifests[kind] = {"source_count": len(manifest["source_sha256"]),
                           "output_count": len(manifest["output_sha256"]),
                           "hash_failures": failures}
    supplements = []
    for i, source in enumerate(sorted(tables.glob("TableS*.csv")), 1):
        target = manuscript / "bmc_additional_files" / f"Additional_file_{i}.csv"
        supplements.append({"number": i, "rows": len(rows(source)),
                            "bytes": target.stat().st_size,
                            "source_matches": sha(source) == sha(target)})
    images = []
    for path in sorted((ROOT / "figures/manuscript_v1").glob("*.png")):
        data = path.read_bytes()
        width, height = struct.unpack(">II", data[16:24])
        images.append({"file": path.name, "width_pixels": width,
                       "height_pixels": height, "bytes": len(data)})
    molecular = [r for r in rows(tables / "TableS6_zero_solution_audit.csv")
                 if r["Model_arm"] != "tissue_null"]
    loto = rows(tables / "TableS1_loto_defined_cells.csv")
    zero_counts = {"outer_molecular_fits": sum(int(r["outer_fits"]) for r in molecular),
                   "outer_zero": sum(int(r["zero_molecular_fits"]) for r in molecular),
                   "final_molecular_arms": len(molecular),
                   "final_zero": sum(int(r["full_development_nonzero_coefficients"]) == 0
                                     for r in molecular)}
    main_path = manuscript / "Manuscript_BMC_Genomics_draft.md"
    main_text = main_path.read_text(encoding="utf-8")
    table2_section = main_text.split("| Compound | Development groups |", 1)[1].split("*Note:*", 1)[0]
    table2_lines = [line for line in table2_section.splitlines()
                    if line.startswith("|") and not line.startswith("| ---")]
    table2_matches = []
    cohort = {r["Compound"]: r for r in rows(tables / "Table1_cohort_by_drug.csv")}
    for source, line in zip(rows(tables / "Table2_primary_effects.csv"), table2_lines):
        cells = [v.strip() for v in line.strip("|").split("|")]
        dev = (f"{float(source['Development_delta_rho_pair_minus_strong_expression']):.5f} "
               f"({float(source['Development_family9_lower']):.5f} to "
               f"{float(source['Development_family9_upper']):.5f})")
        ext = (f"{float(source['CTRPv2_delta_rho_pair_minus_strong_expression']):.5f} "
               f"({float(source['CTRPv2_family9_lower']):.5f} to "
               f"{float(source['CTRPv2_family9_upper']):.5f})")
        macro = source["Compound"] == "Equal-drug macro"
        expected = [source["Compound"], "8-drug mean" if macro else
                    cohort[source["Compound"]]["GDSC2_independent_groups"], dev,
                    "8-drug mean" if macro else source["CTRPv2_shared_groups"], ext]
        table2_matches.append({"compound": source["Compound"], "matches_csv_rounding": cells == expected})
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(manuscript / "Manuscript_BMC_Genomics_draft.docx") as archive:
        doc = ET.fromstring(archive.read("word/document.xml"))
    illustrated_audit_path = manuscript / "main_figure_embedding_audit.json"
    illustrated_audit = (json.loads(illustrated_audit_path.read_text(encoding="utf-8"))
                         if illustrated_audit_path.exists() else None)
    record = {"audited_on": "2026-10-16", "scope": "author-requested snapshot audit",
              "individual_response_values_read": False, "new_model_fits_or_statistics": False,
              "manuscripts": {name: inspect_markdown(manuscript / name) for name in
                              ["Manuscript_BMC_Genomics_draft.md", "Manuscript_draft.md",
                               "Manuscript_BMC_Genomics_with_figures.md"]
                              if (manuscript / name).exists()},
              "manifest_checks": manifests, "supplements": supplements, "figures": images,
              "zero_solution_counts": zero_counts,
              "loto_defined": sum(r["Delta_rho_defined"] == "True" for r in loto),
              "loto_total": len(loto), "table2_matches": table2_matches,
              "docx": {"editable_tables": len(doc.findall(".//w:tbl", ns)),
                       "manual_page_breaks": len(doc.findall('.//w:br[@w:type="page"]', ns)),
                       "page_break_before": len(doc.findall(".//w:pageBreakBefore", ns)),
                       "table_shading_elements": len(doc.findall(".//w:tbl//w:shd", ns))},
              "illustrated_copy_source_current": {
                  file: sha(manuscript / file) == digest
                  for file, digest in illustrated_audit["source_sha256"].items()}
                  if illustrated_audit else {"status": "optional local illustrated copy not present"},
              "notice": "No overall readiness or validity score is calculated. Expected outstanding issues are recorded, not hidden by a single pass flag."}
    (OUT / "snapshot_audit.json").write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"table2_rows_checked": len(table2_matches),
                      "table2_all_match": all(r["matches_csv_rounding"] for r in table2_matches),
                      "manifests_match": all(not r["hash_failures"] for r in manifests.values()),
                      "supplements_match": all(r["source_matches"] for r in supplements),
                      "zero_counts": zero_counts}, ensure_ascii=False))


if __name__ == "__main__":
    main()
