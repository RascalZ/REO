"""Create an author-only Vancouver formatting map for the preserved R1 draft.

This editorial script reads bibliography metadata and manuscript text only. It
never imports analysis/model code or reads study responses. Canonical metadata
and the original manuscript are not overwritten. The six-name truncation is a
consistent implementation illustrated by BMC Genomics' current examples, not a
claim that the guidelines uniquely mandate this cutoff for every source type.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path


SOURCE_MD = "manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md"
SOURCE_BIB = "manuscript_v1/references.bib"
CROSSREF_BIB = "manuscript_v1/references_crossref_normalized.bib"
NUMBER_MAP = "manuscript_v1/reference_number_mapping_scientific_revision_v1.csv"
AIVC_METADATA = "manuscript_v1/aivc_reference_metadata.json"
GUIDELINES = "https://link.springer.com/journal/12864/submission-guidelines/research-article"
GENERAL_GUIDELINES = "https://link.springer.com/journal/12864/submission-guidelines"
ACCESS_DATE = "2026-10-16"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def braced_field(block: str, field: str) -> str:
    """Read the project BibTeX's braced fields, including corporate-name braces."""
    start = re.search(r"(?i)(?<![\w])" + re.escape(field) + r"\s*=\s*\{", block)
    if not start:
        return ""
    depth = 1
    pos = start.end()
    begin = pos
    while pos < len(block):
        char = block[pos]
        if char == "{" and (pos == 0 or block[pos - 1] != "\\"):
            depth += 1
        elif char == "}" and (pos == 0 or block[pos - 1] != "\\"):
            depth -= 1
            if depth == 0:
                return block[begin:pos].strip()
        pos += 1
    raise ValueError(f"Unbalanced {field} field")


def bib_records(path: Path) -> list[dict[str, str]]:
    records = []
    for block in re.split(r"(?m)(?=^@\w+\s*\{)", path.read_text(encoding="utf-8")):
        header = re.match(r"^@(\w+)\s*\{\s*([^,]+),", block)
        if not header:
            continue
        records.append({
            "type": header.group(1).lower(),
            "key": header.group(2).strip(),
            "author": braced_field(block, "author"),
            "doi": braced_field(block, "doi").lower(),
            "url": braced_field(block, "url"),
            "title": braced_field(block, "title"),
        })
    return records


def initials(given: str) -> str:
    parts = []
    for token in given.split():
        segments = token.split("-")
        letters = [next((c.upper() for c in segment if c.isalpha()), "") for segment in segments]
        parts.append("-".join(letter for letter in letters if letter))
    result = "".join(parts)
    if not result:
        raise ValueError(f"Cannot initialize given name: {given!r}")
    return result


def author_list(record: dict[str, str]) -> tuple[list[str], bool]:
    value = record["author"]
    if value.startswith("{") and value.endswith("}"):
        return [value[1:-1]], True
    authors = []
    for raw in value.split(" and "):
        if "," not in raw:
            raise ValueError(f"Missing explicit family/given boundary: {record['key']}: {raw}")
        family, given = raw.split(",", 1)
        authors.append(f"{family.strip()} {initials(given.strip())}")
    # Crossref's author list omits this suffix. The PubMed record and the Nature
    # publisher's author list independently confirm it; preserve the NLM suffix.
    if record["key"] == "Ghandi2019":
        if authors[5] != "McDonald ER":
            raise ValueError("Unexpected sixth Ghandi author; review suffix override")
        authors[5] = "McDonald ER 3rd"
    return authors, False


def author_prefix(reference: str) -> str:
    delimiter = ", et al. "
    if delimiter in reference:
        return reference.split(delimiter, 1)[0] + ", et al."
    return reference.split(". ", 1)[0] + "."


def build_audit(root: Path) -> dict[str, object]:
    records = bib_records(root / SOURCE_BIB)
    metadata = {row["doi"]: row for row in records if row["doi"]}
    resources = {row["key"]: row for row in records if row["type"] == "misc"}
    crossref = {row["doi"]: row for row in bib_records(root / CROSSREF_BIB)}
    manuscript = (root / SOURCE_MD).read_text(encoding="utf-8")
    ref_section = manuscript.split("## References\n", 1)[1].split("\n## ", 1)[0]
    refs = re.findall(r"(?m)^(\d+)\. (.+)$", ref_section)
    if len(refs) != 30:
        raise ValueError("Expected exactly 30 preserved R1 references")
    mapping = list(csv.DictReader((root / NUMBER_MAP).open(encoding="utf-8-sig", newline="")))
    if [int(row["scientific_revision_R1_reference"]) for row in mapping] != list(range(1, 31)):
        raise ValueError("R1 reference-number map is not consecutive")
    rows = []
    for index, (number, reference) in enumerate(refs):
        n = int(number)
        if mapping[index]["reference"] != reference:
            raise ValueError(f"Reference {n} disagrees with preserved R1 number map")
        doi_match = re.search(r"https://doi\.org/(\S+)", reference)
        doi = doi_match.group(1).rstrip(".").lower() if doi_match else ""
        if doi:
            source = metadata[doi]
        elif n == 25:
            source = resources["CellModelPassportsAPI2026"]
        elif n == 26:
            source = resources["GDSCdocs2026"]
        else:
            raise ValueError(f"Unresolved reference identity: {n}")
        authors, corporate = author_list(source)
        displayed = authors if corporate else authors[:6]
        formatted = ", ".join(displayed)
        et_al = not corporate and len(authors) > 6
        formatted += ", et al." if et_al else "."
        old = author_prefix(reference)
        new_reference = formatted + reference[len(old):]
        notes = []
        if corporate:
            notes.append("Corporate author retained verbatim; personal-name truncation does not apply.")
        if n == 1:
            notes.append("Sixth-author 3rd suffix restored from PubMed PMID 31068700 and Nature publisher author list; Crossref omits the suffix.")
        if n == 28:
            notes.append("Conference author metadata cross-checked against saved official NeurIPS BibTeX and Crossref; no journal abbreviation applies.")
        if n == 30:
            notes.append("Six authors are all retained without et al.; ahead-of-print year/DOI and intentionally unassigned volume/pages are unchanged.")
        cr = crossref.get(doi)
        if cr and source["author"] != cr["author"]:
            raise ValueError(f"Canonical/Crossref author mismatch requires review: {n}")
        rows.append({
            "reference_number": n,
            "citation_key": source["key"],
            "entry_type": source["type"],
            "doi": doi,
            "author_count": len(authors),
            "corporate_author": corporate,
            "displayed_author_count": len(displayed),
            "et_al": et_al,
            "old_author_prefix": old,
            "new_author_prefix": formatted,
            "changed": old != formatted,
            "source_author_field": source["author"],
            "metadata_evidence": (
                "canonical BibTeX agrees with saved normalized Crossref author field"
                if cr else "saved canonical corporate/conference metadata"
            ),
            "old_reference": reference,
            "new_reference": new_reference,
            "notes": " ".join(notes),
        })
    return {
        "schema_version": "1.0",
        "scope": "M08 reference author typography only; no numbered citation, title, DOI, year, journal, volume or page changes",
        "accessed_on": ACCESS_DATE,
        "selected_policy": "Surname plus unpunctuated initials; preserve compound surnames, accents and hyphenated given-name initials; list all personal authors up to six, otherwise first six plus et al.; retain corporate authors verbatim.",
        "policy_status": "Consistent Vancouver implementation illustrated by the current BMC Genomics six-author-plus-et-al example; the page does not explicitly state a universal mandatory numeric cutoff.",
        "sources": [
            {"url": GUIDELINES, "purpose": "Current research-article Vancouver examples, including six names before et al. and hyphenated initials in the dataset example", "accessed_on": ACCESS_DATE},
            {"url": GENERAL_GUIDELINES, "purpose": "Current main-text guidelines specify Index Medicus/MEDLINE journal abbreviations", "accessed_on": ACCESS_DATE},
            {"url": "https://pubmed.ncbi.nlm.nih.gov/31068700/", "purpose": "Verify McDonald 3rd suffix in Ghandi 2019 sixth author", "accessed_on": ACCESS_DATE},
            {"url": "https://www.nature.com/articles/s41586-019-1186-3", "purpose": "Publisher author list confirms E. Robert McDonald III", "accessed_on": ACCESS_DATE},
        ],
        "source_sha256": {p: digest(root / p) for p in [SOURCE_MD, SOURCE_BIB, CROSSREF_BIB, NUMBER_MAP, AIVC_METADATA]},
        "summary": {
            "reference_count": len(rows),
            "personal_author_references": sum(not row["corporate_author"] for row in rows),
            "corporate_author_references": sum(row["corporate_author"] for row in rows),
            "author_prefixes_changed": sum(row["changed"] for row in rows),
            "et_al_references_after_formatting": sum(row["et_al"] for row in rows),
            "unresolved_author_metadata": 0,
            "new_scientific_results": 0,
        },
        "replacements": rows,
    }


def apply_author_replacements(manuscript: str, audit: dict[str, object]) -> str:
    before, remainder = manuscript.split("## References\n", 1)
    refs, after = remainder.split("\n## ", 1)
    original_refs = refs
    for row in audit["replacements"]:
        old = f"{row['reference_number']}. {row['old_reference']}"
        new = f"{row['reference_number']}. {row['new_reference']}"
        if refs.splitlines().count(old) != 1:
            raise ValueError(f"Reference {row['reference_number']} is missing or duplicated")
        refs = refs.replace(old, new, 1)
    # Non-author reference tails must remain byte-for-byte identical as text.
    for row in audit["replacements"]:
        assert row["old_reference"][len(row["old_author_prefix"]):] == row["new_reference"][len(row["new_author_prefix"]):]
    assert len(re.findall(r"(?m)^\d+\. ", original_refs)) == len(re.findall(r"(?m)^\d+\. ", refs)) == 30
    return before + "## References\n" + refs + "\n## " + after


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-manuscript", type=Path, help="Optional NEW Markdown output; never the original R1 file")
    args = parser.parse_args()
    root = args.root.resolve()
    audit = build_audit(root)
    directory = root / "docs/m08_v1"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "reference_author_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    with (directory / "reference_author_audit.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit["replacements"][0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit["replacements"])
    summary = audit["summary"]
    report = f"""# M08 reference-author audit

Audit date: {ACCESS_DATE}. This is an editorial author-formatting audit, not a new literature search or scientific evaluation.

The 30 reference identities and numbering remain fixed. The author map changes {summary['author_prefixes_changed']} prefixes; {summary['et_al_references_after_formatting']} references retain six names followed by `et al.`. All authors are named for papers with at most six authors. Three corporate authors remain verbatim. No author metadata are unresolved.

The selected policy preserves family-name particles, compound surnames, accents and hyphenated given-name initials. Names are based on the saved verified bibliography and normalized Crossref records. References 27–30 also have the saved Crossref/PubMed/official NeurIPS evidence in `manuscript_v1/aivc_reference_metadata.json`. The original BibTeX files are retained.

The [current BMC Genomics research-article instructions]({GUIDELINES}) identify Vancouver examples and illustrate six authors before `et al.`. The page does not explicitly prescribe the numeric cutoff in prose; six-name truncation is recorded as this manuscript's consistent implementation, rather than presented as a uniquely mandated rule. Its dataset example also preserves hyphenated given-name initials. The [general instructions]({GENERAL_GUIDELINES}) specify Index Medicus/MEDLINE journal abbreviations; this author-only map leaves the previously normalized journal abbreviations unchanged.

One targeted metadata correction is documented: Ghandi 2019's sixth author is rendered `McDonald ER 3rd`. Saved Crossref metadata lack the suffix, but [PubMed PMID 31068700](https://pubmed.ncbi.nlm.nih.gov/31068700/) and the [Nature author list](https://www.nature.com/articles/s41586-019-1186-3) confirm it. Corporate authors, the conference source type, and the existing ahead-of-print citation are preserved.

`reference_author_audit.json` contains the source hashes and all 30 exact old/new author prefixes and complete reference lines. `format_m08_reference_authors.py` exposes `build_audit(root)` and `apply_author_replacements(text, audit)` for the M08 manuscript builder; default execution writes this audit only. The transform verifies that every reference tail after its author block is unchanged. It never reads study data, fits a model, or scores predictions.
"""
    (directory / "reference_author_audit.md").write_text(report, encoding="utf-8", newline="\n")
    if args.output_manuscript:
        output = args.output_manuscript.resolve()
        if output == (root / SOURCE_MD).resolve() or output.exists():
            raise ValueError("Refusing to overwrite the original or an existing output manuscript")
        original = (root / SOURCE_MD).read_text(encoding="utf-8")
        output.write_text(apply_author_replacements(original, audit), encoding="utf-8", newline="\n")
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
