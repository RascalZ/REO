"""Merge user-designated legacy author facts into a LOCAL R1 reading copy.

Requires the private, source-bound author extraction. Never guesses missing facts,
publishes personal data, imports modeling code or accesses external test labels.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/author_review_v1"
PRIVATE = DOCS / "local_private"
OUT = ROOT / "manuscript_v1/author_review_v1"
MS = ROOT / "manuscript_v1"
BASE = MS / "Manuscript_BMC_Genomics_M08_v1.md"
COMMIT = "e4820b333655b59e47dd944ca1a2e08a0d68b2bd"
REPO_URL = "https://github.com/RascalZ/REO/tree/" + COMMIT


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def record(path: Path) -> dict:
    return {"path": path.relative_to(ROOT).as_posix(), "sha256": sha(path)}


def section(md: str, title: str) -> str:
    return md.split(f"## {title}\n", 1)[1].split("\n## ", 1)[0]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PACKAGE = OUT / "Additional_files"
    PACKAGE.mkdir(exist_ok=True)
    evidence_path = PRIVATE / "source_author_declarations.json"
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    legacy = Path(evidence["source"]["absolute_path"])
    assert sha(legacy) == evidence["source"]["sha256"], "Legacy source changed"
    source = BASE.read_text(encoding="utf-8")
    md = source
    original_status = re.search(r"^\*\*Draft status:\*\* .+$", md, re.M).group(0)
    md = md.replace(original_status,
        "**Draft status:** Author reading copy, 20 October 2026. Scientific results are unchanged from R1-M08. "
        "Author and funding information has been adopted from the user-designated prior submission. "
        "This copy is for author review; final approval and release completion remain pending.")
    ab = evidence["author_block"]
    author_parts = []
    for author in ab["authors"]:
        assert author["affiliations"] == ["1"]
        author_parts.append(author["name"] + "¹" + ("†" if author["equal_contribution"] else "") + ("*" if author["corresponding"] else ""))
    affiliation = ab["affiliations"][0]["text"].strip()
    affiliation = re.sub(r"^1\s+", "¹ ", affiliation)
    corresponding = [a["name"] for a in ab["authors"] if a["corresponding"]]
    first_equal = [a["name"] for a in ab["authors"] if a["equal_contribution"]]
    md = re.sub(r"^\*\*Authors:\*\* .+$", "**Authors:** " + ", ".join(author_parts), md, flags=re.M)
    md = re.sub(r"^\*\*Affiliations:\*\* .+$", "**Affiliations:** " + affiliation, md, flags=re.M)
    md = re.sub(r"^\*\*Corresponding author:\*\* .+$",
        "**Corresponding authors:** " + " and ".join(corresponding) + ".\n\n† " + " and ".join(first_equal) + " contributed equally to this work.", md, flags=re.M)
    old_availability = re.search(r"^The repository contains source versions,.+$", md, re.M).group(0)
    md = md.replace(old_availability,
        "The repository contains source versions, download/accession provenance, code, aggregate results, "
        "figure/table scripts, checksums, and recorded environment metadata. Sample-level CTRPv2 responses, "
        "group bootstrap weights, individual predictions, and fitted checkpoints are excluded from the version-controlled repository. "
        "Reproduction from source data requires the documented identity/annotation and environment records; "
        "the included checks verify aggregate consistency rather than independently reproducing fitted predictions. "
        "Upstream terms remain applicable and are not replaced by the repository code license. "
        "The REO project repository is currently private; its contents and remaining access limitations are specified in the data and code availability statements.")
    ai_old = re.search(r"^\[AUTHOR INPUT REQUIRED: Confirm the generative-AI disclosure.+$", md, re.M).group(0)
    ai_new = ("OpenAI Codex (OpenAI; used during October 2026) assisted with study-planning discussions, "
        "comparison of methodological options, implementation and debugging of analysis and audit code, "
        "reference checking, interpretation drafts, manuscript drafting and editing, and figure/table and document preparation. "
        "Numerical results were produced by the recorded analysis programs from the specified datasets, "
        "rather than generated as empirical observations by a language model. "
        "AI assistance does not transfer authorship or accountability: the human authors remain responsible for "
        "scientific judgment, verification, and approval of the final manuscript.")
    md = md.replace(ai_old, ai_new)
    declarations = evidence["declarations"]
    competing = declarations["competing_interests"]["verbatim"]
    funding = declarations["funding"]["normalized_body"]
    new_declarations = f"""## Declarations

### Ethics approval and consent to participate

This study was a secondary computational analysis of publicly available cancer-cell-line datasets. No new human participants were recruited, no new human tissue was collected, and no new animal experiments were performed.

### Consent for publication

No identifiable individual participant information is presented in this manuscript.

### Availability of data and materials

The analysis used public GDSC2 drug-response, historical CCLE expression, CTRPv2 drug-response, DepMap identity, Cell Model Passports, and GENCODE annotation resources. The source releases, exact members and checksums, identity rules, and processing provenance are documented in Additional files 8 and 14 and the repository data-availability and input ledgers. Source datasets should be obtained from the original providers under their applicable terms. The code repository does not redistribute the original expression and drug-response matrices or relicense upstream resources.

Complete derived identity mappings, original fold definitions, 64 fixed numeric models, feature coefficients, and preprocessing states have been prepared locally. They are not currently included in the version-controlled repository; this manuscript does not claim that the complete trained-model assets are already publicly downloadable. Their permitted access arrangements and formal release documentation remain to be finalized before submission.

### Code availability

The reconstruction code, protocols, aggregate results, figure/table generators, and a response-free inference example are maintained in the REO project repository [31], which was verified as private on 20 October 2026. Public or reviewer access and an archival release remain to be arranged before submission; this reading copy does not claim that these materials are already publicly accessible. The cited source snapshot is commit {COMMIT}; it is an identifiable Git source snapshot, not a separately deposited archival DOI. The newly authored code is licensed under the MIT License. The included example is hand-specified and is not one of the paper's fitted models.

### Competing interests

{competing}

### Funding

{funding}

### Authors' contributions

{" and ".join(first_equal)} contributed equally to this work. Author-specific contribution roles will be completed by the listed authors before submission.

### Acknowledgements

The acknowledgement statement will be finalized by the authors before submission.

"""
    md = md.split("## Declarations\n", 1)[0] + new_declarations + "## References\n" + md.split("## References\n", 1)[1]
    new_ref = (f"31. REO reconstruction repository. Analysis code and aggregate results. "
               f"GitHub; source snapshot {COMMIT}. {REPO_URL}. Private repository; access status verified 20 October 2026.")
    refs_end = md.index("\n## Additional files")
    md = md[:refs_end].rstrip() + "\n" + new_ref + "\n\n" + md[refs_end:]
    assert "[AUTHOR INPUT REQUIRED" not in md
    # All empirical Results, Discussion, Conclusions, figure captions and table cells stay verbatim.
    for title in ["Abstract", "Keywords", "Background", "Results", "Discussion", "Conclusions", "Additional files", "Figure legends"]:
        assert section(source, title) == section(md, title), f"Scientific section changed: {title}"
    assert re.findall(r"(?m)^\|.*$", source) == re.findall(r"(?m)^\|.*$", md)
    main_md = OUT / "Manuscript_BMC_Genomics_author_review_v1.md"
    main_md.write_text(md, encoding="utf-8", newline="\n")
    bib_source = MS / "references.bib"
    bib = bib_source.read_text(encoding="utf-8") + f"\n@misc{{REOReconstruction2026,\n  title = {{REO reconstruction repository: analysis code and aggregate results}},\n  howpublished = {{GitHub}},\n  url = {{{REPO_URL}}},\n  note = {{Source snapshot {COMMIT}; private repository, access status verified 20 October 2026}}\n}}\n"
    (OUT / "references_author_review_v1.bib").write_text(bib, encoding="utf-8", newline="\n")
    supp_md = OUT / "Supplementary_methods_author_review_v1.md"
    shutil.copyfile(MS / "Supplementary_methods_scientific_revision_v1.md", supp_md)
    node = shutil.which("node") or r"D:\Env\NodeJS\node.exe"
    builder = MS / "build_bmc_docx_author_review_v1.js"
    for path, embed in [(main_md, True), (supp_md, False)]:
        cmd = [node, str(builder), "--input", str(path), "--output", str(path.with_suffix(".docx"))]
        subprocess.run(cmd + (["--embed-figures"] if embed else []), cwd=ROOT, check=True)
    files = []
    original_package = MS / "bmc_additional_files_M08_v1"
    for n in range(1, 15):
        old = original_package / f"Additional_file_{n}{'.docx' if n == 14 else '.csv'}"
        src = supp_md.with_suffix(".docx") if n == 14 else old
        target = PACKAGE / old.name
        shutil.copyfile(src, target)
        assert sha(target) == sha(src)
        files.append({"number": n, **record(target), "source": record(src),
                      "same_bytes_as_M08": sha(target) == sha(old),
                      "scientific_content_unchanged": True})
    manifest = {"date": "2026-10-20", "status": "local_author_reading_copy_not_submission_approved",
                "author_source_sha256": sha(legacy), "source_manuscript": record(BASE),
                "source_bibliography": record(bib_source), "builder": record(builder), "generator": record(Path(__file__)),
                "source_extraction": record(evidence_path),
                "documents": [record(p) for p in [main_md, main_md.with_suffix(".docx"), supp_md, supp_md.with_suffix(".docx")]],
                "additional_files": files, "copied_personal_author_count": len(ab["authors"]),
                "scientific_sections_unchanged": True, "conference_content_added": False,
                "fitting_resampling_scoring_or_test_labels_accessed": False,
                "all_authors_approved_new_text": False, "immutable_archival_doi_created": False,
                "repository_visibility": "private", "repository_public_access_available": False,
                "source_missing_items_are_not_assumed_confirmed": True}
    dump(OUT / "package_manifest.json", manifest)
    # A safe public status record deliberately omits author names and source text.
    dump(DOCS / "migration_summary.json", {"date": "2026-10-20", "status": manifest["status"],
         "source_sha256": sha(legacy), "author_count": len(ab["authors"]), "equal_contributors": len(first_equal),
         "corresponding_authors": len(corresponding), "funding_sources": len(declarations["funding"]["funders"]),
         "competing_interest_statement_reused": True, "reference_count": 31,
         "source_missing_fields": ["correspondence_email", "author_specific_contribution_roles", "formal_ethics_applicability_statement", "acknowledgements", "new_version_approval", "data_model_access_permissions", "archive_identifier", "article_license_selection"],
         "ai_assistance_disclosure": "current_recorded_work_not_legacy_absence",
         "outputs_local_git_ignored": True, "scientific_sections_unchanged": True,
         "repository_visibility": "private", "repository_public_access_available": False,
         "no_model_fitting_scoring_or_test_response_access": True,
         "local_main_docx_sha256": sha(main_md.with_suffix(".docx")),
         "local_package_manifest_sha256": sha(OUT / "package_manifest.json")})
    print(json.dumps({"main": main_md.with_suffix(".docx").relative_to(ROOT).as_posix(),
                      "author_count": len(ab["authors"]), "additional_files": len(files), "scientific_sections_unchanged": True}))


if __name__ == "__main__":
    main()
