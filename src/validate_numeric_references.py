"""Validate numbered references in the Markdown manuscript and its BibTeX file.

The bundled citation-management validator recognizes LaTeX citation keys, whereas
the working manuscript uses numeric Vancouver-style citations. This small checker
validates that numeric markers in the Introduction, Methods, and Discussion map
to the numbered reference list and that no numbered reference is left uncited.
It does not verify the truth of a cited claim or resolve DOI URLs.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def expand_numeric_group(group: str, nrefs: int) -> set[int] | None:
    """Expand a numeric citation group; return None for non-citation numbers."""
    group = group.strip()
    # Avoid interpreting sample/row counts such as [1,185] as refs [1, 185].
    if re.fullmatch(r"\d{1,3}(?:,\d{3})+", group):
        return None
    if not re.fullmatch(r"\d+(?:\s*[–-]\s*\d+)?(?:\s*,\s*\d+(?:\s*[–-]\s*\d+)?)*", group):
        return None
    found: set[int] = set()
    for item in group.split(","):
        item = item.strip().replace("–", "-")
        if "-" in item:
            start, end = (int(value.strip()) for value in item.split("-", 1))
            if start > end:
                return None
            found.update(range(start, end + 1))
        else:
            found.add(int(item))
    # An out-of-range item in a bracket is usually a data count, not a citation.
    if not found or min(found) < 1 or max(found) > nrefs:
        return None
    return found


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", default="manuscript_v1/Manuscript_draft.md")
    parser.add_argument("--bib", default="manuscript_v1/references.bib")
    parser.add_argument("--crossref", default="manuscript_v1/references_crossref_normalized.bib")
    parser.add_argument("--output", default="manuscript_v1/citation_validation.json")
    args = parser.parse_args()

    manuscript = Path(args.manuscript).read_text(encoding="utf-8")
    refs_text = manuscript.split("## References", 1)[1].split("## Evidence source map", 1)[0]
    reference_rows = re.findall(r"(?m)^(\d+)\.\s+(.+)$", refs_text)
    reference_numbers = {int(number) for number, _ in reference_rows}
    nrefs = len(reference_numbers)

    body = manuscript.split("## References", 1)[0]
    active_sections: list[str] = []
    section = ""
    for line in body.splitlines():
        if line.startswith("## "):
            section = line[3:].strip()
        if section in {"1. Introduction", "2. Methods", "4. Discussion", "Background", "Methods", "Discussion"}:
            active_sections.append(line)
    citations: set[int] = set()
    ignored_numeric_brackets: list[str] = []
    for group in re.findall(r"\[([^\]]+)\]", "\n".join(active_sections)):
        expanded = expand_numeric_group(group, nrefs)
        if expanded is None:
            if re.fullmatch(r"[\d,–\-\s]+", group):
                ignored_numeric_brackets.append(group)
            continue
        citations.update(expanded)

    bib = Path(args.bib).read_text(encoding="utf-8")
    crossref = Path(args.crossref).read_text(encoding="utf-8")
    bib_blocks = [block for block in re.split(r"(?m)(?=^@\w+\{)", bib) if block.strip()]
    bib_entries = re.findall(r"(?m)^@(\w+)\{([^,]+),", bib)
    bib_keys = [key for _, key in bib_entries]
    bib_article_count = sum(entry_type.lower() == "article" for entry_type, _ in bib_entries)
    bib_doi_by_type: dict[str, set[str]] = {"article": set(), "misc": set()}
    for block in bib_blocks:
        header = re.match(r"(?m)^@(\w+)\{", block)
        doi_match = re.search(r"(?im)^\s*doi\s*=\s*\{([^}]+)\}", block)
        if header and doi_match:
            bib_doi_by_type.setdefault(header.group(1).lower(), set()).add(doi_match.group(1).lower())
    bib_dois = set().union(*bib_doi_by_type.values())
    reference_dois = {
        doi.lower()
        for _, reference in reference_rows
        for doi in re.findall(r"10\.\d{4,}/[^\s.]+(?:\.[^\s.]+)*", reference, flags=re.I)
    }
    crossref_blocks = [block for block in re.split(r"(?m)(?=^@\w+\{)", crossref) if block.strip()]
    crossref_dois = [
        doi.lower()
        for block in crossref_blocks
        if re.match(r"(?m)^@article\{", block)
        # The raw Crossref exporter writes each record as a one-line entry,
        # so the DOI field is not anchored to a line start.
        for doi in re.findall(r"(?i)doi\s*=\s*\{([^}]+)\}", block)
    ]

    result = {
        "manuscript": args.manuscript,
        "bibliography": args.bib,
        "numbered_references": sorted(reference_numbers),
        "citations_found_in_intro_methods_discussion": sorted(citations),
        "uncited_numbered_references": sorted(reference_numbers - citations),
        "citations_without_reference": sorted(citations - reference_numbers),
        "ignored_numeric_brackets_not_treated_as_citations": sorted(set(ignored_numeric_brackets)),
        "bibtex_entry_count": len(bib_entries),
        "bibtex_unique_keys": len(set(bib_keys)),
        "bibtex_key_format_valid": all(re.fullmatch(r"[A-Za-z0-9]+", key) for key in bib_keys),
        "bibtex_article_count": bib_article_count,
        "bibtex_article_dois": sorted(bib_doi_by_type.get("article", set())),
        "bibtex_nonarticle_dois": sorted(bib_dois - bib_doi_by_type.get("article", set())),
        "numbered_reference_doi_count": len(reference_dois),
        "bibtex_doi_count": len(bib_dois),
        "reference_dois_missing_from_bibtex": sorted(reference_dois - bib_dois),
        "bibtex_dois_missing_from_references": sorted(bib_dois - reference_dois),
        "crossref_article_entry_count": len(re.findall(r"(?m)^@article\{", crossref)),
        "crossref_doi_field_count": len(crossref_dois),
        "crossref_unique_doi_count": len(set(doi.lower() for doi in crossref_dois)),
        "pass": (
            bool(nrefs)
            and reference_numbers == set(range(1, nrefs + 1))
            and citations == reference_numbers
            and len(bib_entries) == len(set(bib_keys))
            and all(re.fullmatch(r"[A-Za-z0-9]+", key) for key in bib_keys)
            and len(bib_entries) == nrefs
            and bib_article_count == 23
            and len(bib_dois) == 24
            and reference_dois == bib_dois
            and bib_doi_by_type.get("article", set()) == set(crossref_dois)
            and len([block for block in crossref_blocks if re.match(r"(?m)^@article\{", block)]) == 23
            and len(crossref_dois) == 23
            and len(set(doi.lower() for doi in crossref_dois)) == 23
        ),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
