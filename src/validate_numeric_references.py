"""Validate numbered references in the Markdown manuscript and its BibTeX file.

The bundled citation-management validator recognizes LaTeX citation keys, whereas
the working manuscript uses numeric Vancouver-style citations. This small checker
validates that numeric markers in the Introduction, Methods, and Discussion map
to the numbered reference list and that no numbered reference is left uncited.
It does not verify the truth of a cited claim or resolve DOI URLs.
"""

from __future__ import annotations

import argparse
from collections import Counter
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
    # Zero is not a reference number (for example the model domain [0,1]).
    # Keep positive out-of-range citations so that a missing reference fails.
    if not found or min(found) < 1:
        return None
    return found


def normalize_doi(value: str) -> str:
    """Compare DOI identifiers independently of URL prefix and letter case."""
    value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value.strip(), flags=re.I)
    value = re.sub(r"^doi:\s*", "", value, flags=re.I)
    return value.rstrip(".,;").lower()


def bibtex_records(text: str) -> list[dict[str, object]]:
    """Read entry headers and DOI fields in multiline or one-line BibTeX.

    This is a structural checker for the project's exported bibliography, not a
    general BibTeX evaluator. DOI fields must use braces or double quotes.
    """
    records: list[dict[str, object]] = []
    for block in re.split(r"(?im)(?=^\s*@\w+\s*\{)", text):
        header = re.match(r"(?im)^\s*@(\w+)\s*\{\s*([^,]+),", block)
        if not header:
            continue
        doi_fields = [
            normalize_doi(braced if braced else quoted)
            for braced, quoted in re.findall(
                r'(?i)(?<![\w])doi\s*=\s*(?:\{([^}]*)\}|"([^"]*)")', block
            )
        ]
        records.append({
            "type": header.group(1).lower(),
            "key": header.group(2).strip(),
            "dois": doi_fields,
        })
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", default="manuscript_v1/Manuscript_draft.md")
    parser.add_argument("--bib", default="manuscript_v1/references.bib")
    parser.add_argument("--crossref", default="manuscript_v1/references_crossref_normalized.bib")
    parser.add_argument("--output", default="manuscript_v1/citation_validation.json")
    args = parser.parse_args()

    manuscript = Path(args.manuscript).read_text(encoding="utf-8")
    refs_text = manuscript.split("## References", 1)[1].split("\n## ", 1)[0]
    reference_rows = re.findall(r"(?m)^(\d+)\.\s+(.+)$", refs_text)
    reference_number_counts = Counter(int(number) for number, _ in reference_rows)
    reference_numbers = set(reference_number_counts)
    nrefs = len(reference_rows)

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
    bib_entries = bibtex_records(bib)
    crossref_entries = bibtex_records(crossref)
    bib_keys = [str(entry["key"]) for entry in bib_entries]
    article_entries = [entry for entry in bib_entries if entry["type"] == "article"]
    crossref_articles = [entry for entry in crossref_entries if entry["type"] == "article"]
    bib_article_count = len(article_entries)
    bib_type_counts = Counter(str(entry["type"]) for entry in bib_entries)
    bib_doi_by_type: dict[str, set[str]] = {}
    bib_doi_counts: Counter[str] = Counter()
    for entry in bib_entries:
        dois = entry["dois"]
        bib_doi_by_type.setdefault(str(entry["type"]), set()).update(dois)
        bib_doi_counts.update(dois)
    bib_dois = set().union(*bib_doi_by_type.values())
    reference_doi_rows = [
        (int(number), [normalize_doi(doi) for doi in re.findall(
            r"10\.\d{4,}/[^\s.]+(?:\.[^\s.]+)*", reference, flags=re.I
        )])
        for number, reference in reference_rows
    ]
    reference_doi_counts = Counter(
        doi for _, dois in reference_doi_rows for doi in dois
    )
    reference_dois = set(reference_doi_counts)
    crossref_dois = [
        doi for entry in crossref_articles for doi in entry["dois"]
    ]
    article_dois = bib_doi_by_type.get("article", set())
    valid_doi = re.compile(r"10\.\d{4,9}/[^\s]+", flags=re.I)
    article_doi_invalid = [
        str(entry["key"]) for entry in article_entries
        if len(entry["dois"]) != 1 or not valid_doi.fullmatch(entry["dois"][0])
    ]
    crossref_doi_invalid = [
        str(entry["key"]) for entry in crossref_articles
        if len(entry["dois"]) != 1 or not valid_doi.fullmatch(entry["dois"][0])
    ]
    checks = {
        "reference_list_nonempty": bool(nrefs),
        "reference_numbers_unique": all(count == 1 for count in reference_number_counts.values()),
        "reference_numbers_continuous_and_in_order": (
            [int(number) for number, _ in reference_rows] == list(range(1, nrefs + 1))
        ),
        "all_numbered_references_cited_and_resolved": citations == reference_numbers,
        "bibtex_keys_unique": len(bib_entries) == len(set(bib_keys)),
        "bibtex_key_format_valid": all(re.fullmatch(r"[A-Za-z0-9]+", key) for key in bib_keys),
        "bibtex_count_matches_numbered_references": len(bib_entries) == nrefs,
        "every_article_has_exactly_one_valid_doi": not article_doi_invalid,
        "bibtex_has_no_duplicate_doi_fields_or_identifiers": (
            all(len(entry["dois"]) <= 1 for entry in bib_entries)
            and all(count == 1 for count in bib_doi_counts.values())
        ),
        "bibtex_dois_valid": all(valid_doi.fullmatch(doi) for doi in bib_dois),
        "reference_dois_unique_and_one_per_row": (
            all(count == 1 for count in reference_doi_counts.values())
            and all(len(dois) <= 1 for _, dois in reference_doi_rows)
        ),
        "reference_and_bibtex_doi_sets_match": reference_dois == bib_dois,
        # The Crossref archive is for journal articles only. Conference/resource
        # entries remain in the bibliography but must not masquerade as articles.
        "crossref_archive_contains_journal_articles_only": len(crossref_entries) == len(crossref_articles),
        "every_crossref_article_has_exactly_one_valid_doi": not crossref_doi_invalid,
        "crossref_article_count_matches_bibtex_articles": len(crossref_articles) == bib_article_count,
        "crossref_article_dois_unique": len(crossref_dois) == len(set(crossref_dois)),
        "crossref_and_bibtex_article_doi_sets_match": article_dois == set(crossref_dois),
    }

    result = {
        "manuscript": args.manuscript,
        "bibliography": args.bib,
        "numbered_references": sorted(reference_numbers),
        "numbered_reference_row_count": nrefs,
        "duplicate_numbered_references": sorted(number for number, count in reference_number_counts.items() if count > 1),
        "citations_found_in_intro_methods_discussion": sorted(citations),
        "uncited_numbered_references": sorted(reference_numbers - citations),
        "citations_without_reference": sorted(citations - reference_numbers),
        "ignored_numeric_brackets_not_treated_as_citations": sorted(set(ignored_numeric_brackets)),
        "bibtex_entry_count": len(bib_entries),
        "bibtex_unique_keys": len(set(bib_keys)),
        "bibtex_key_format_valid": all(re.fullmatch(r"[A-Za-z0-9]+", key) for key in bib_keys),
        "bibtex_article_count": bib_article_count,
        "bibtex_inproceedings_count": bib_type_counts.get("inproceedings", 0),
        "bibtex_entry_type_counts": dict(sorted(bib_type_counts.items())),
        "bibtex_articles_without_exactly_one_valid_doi": article_doi_invalid,
        "bibtex_article_dois": sorted(bib_doi_by_type.get("article", set())),
        "bibtex_nonarticle_dois": sorted(bib_dois - bib_doi_by_type.get("article", set())),
        "numbered_reference_doi_count": len(reference_dois),
        "bibtex_doi_count": len(bib_dois),
        "reference_dois_missing_from_bibtex": sorted(reference_dois - bib_dois),
        "bibtex_dois_missing_from_references": sorted(bib_dois - reference_dois),
        "crossref_article_entry_count": len(crossref_articles),
        "crossref_nonarticle_entry_count": len(crossref_entries) - len(crossref_articles),
        "crossref_articles_without_exactly_one_valid_doi": crossref_doi_invalid,
        "crossref_doi_field_count": len(crossref_dois),
        "crossref_unique_doi_count": len(set(doi.lower() for doi in crossref_dois)),
        "crossref_article_dois_missing_from_bibtex": sorted(set(crossref_dois) - article_dois),
        "bibtex_article_dois_missing_from_crossref": sorted(article_dois - set(crossref_dois)),
        "checks": checks,
        "failed_checks": [name for name, passed in checks.items() if not passed],
        "pass": all(checks.values()),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
