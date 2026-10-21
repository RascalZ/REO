"""Build a question-led LOCAL reading revision from frozen manuscript evidence.

No modeling modules, study input matrices, fitted artifacts or individual response
files are loaded. The previous author reading copy is a protected source.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/narrative_v1"
OUT = ROOT / "manuscript_v1/author_review_v2"
OLD = ROOT / "manuscript_v1/author_review_v1"


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p: Path, obj: object) -> None:
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def section(md: str, name: str) -> str:
    return md.split(f"## {name}\n", 1)[1].split("\n## ", 1)[0].strip()


def replace_section(md: str, name: str, text: str) -> str:
    old = section(md, name)
    assert md.count(old) == 1
    return md.replace(old, text.strip(), 1)


def display(text: str) -> str:
    return text.replace("CTRPv2", "CTRP v2")


def revise_results(source: str) -> tuple[str, list[dict]]:
    result = source
    edits: list[dict] = []

    def change(old: str, new: str, purpose: str) -> None:
        nonlocal result
        assert result.count(old) == 1, f"Missing or repeated Results target: {purpose}"
        result = result.replace(old, new, 1)
        edits.append({"purpose": purpose, "old_text_sha256": hashlib.sha256(old.encode()).hexdigest(),
                      "new_text_sha256": hashlib.sha256(new.encode()).hexdigest()})

    for a, b in [
        ("Cohort composition and evaluation coverage", "A common cohort supported paired model comparisons"),
        ("REO-pair versus strong-expression performance in development", "Pair encoding did not establish an average gain over strong expression"),
        ("Matched-candidate representation comparisons", "Matched-gene contrasts varied with the candidate rule"),
        ("Zero molecular solutions and fitted-model status", "Zero molecular solutions retained tissue-reference predictions"),
        ("Leave-one-tissue-out sensitivity", "Tissue-held-out contrasts were frequently not estimable"),
        ("Locked CTRPv2 cross-screen evaluation", "Cross-screen prediction did not establish a pair-specific gain"),
        ("Finite synthetic-control results", "Constructed controls illustrated distinct prediction behaviors"),
    ]:
        change("### " + a, "### " + b, "Finding-led heading: " + a)
    change(
        "The repeated nested development evaluation completed all 960 outer tasks, 26,160 fit events, and 119,208 outer predictions. No solver warnings or fallback fits were recorded. ",
        "", "Remove duplicate execution counts already retained in Methods; solver status retained under zero solutions")
    change(
        "The primary equal-drug macro Spearman correlation was 0.2931",
        "The primary comparison did not establish a ranking benefit from the evaluated REO pipeline. The equal-drug macro Spearman correlation was 0.2931",
        "Lead the primary comparison with the research finding")
    change(
        "The macro difference was negative in all three repeated partitions (−0.00393, −0.00655, and −0.00716), which describes partition-level direction consistency but does not create three independent replications.",
        "The macro difference was negative in all three repeated partitions. This describes consistency across partitions, not independent replication; repeat-level values remain in Additional file 5 (Table S5).",
        "Remove repeated per-partition point values; retain direction and complete Table S5")
    change(
        "The tissue-only macro Spearman was 0.2366;",
        "\n\nThe observed development increment over tissue was shared by both primary molecular models. The tissue-only macro Spearman was 0.2366;",
        "Separate the molecular increment from the primary representation contrast")
    change(
        "The response-association rule produced a descriptive macro difference of +0.05654 over the tissue-only reference, while the strong-expression model differed from the same reference by +0.06242. These results are consistent with molecular predictions adding ranking signal beyond tissue means in some parts of the evaluation, but the average increment was not specific to pair encoding. Candidate selection based on response-associated continuous expression was performed within training folds, yet the contrast remains conditional on the selected outcome-informed workflow.",
        "The matched comparisons concerned both representation and gene selection: the response-associated rule used continuous expression–response association within training tissues. It was nested within training folds, so the observed pattern describes the specified selection-and-learning procedures.",
        "Remove duplicate tissue-increment values; preserve training-only outcome-informed candidate policy")
    change(
        "The rank arm further distinguishes order encoding from the complete pair expansion (Table S10). Under MAD selection, rank minus expression was −0.01244 (existing nominal conditional 95% interval, −0.01879 to −0.00622), and pair minus rank was −0.00521 (−0.01222 to +0.00181). Under response-associated selection, rank minus expression was +0.00441 (−0.00150 to +0.01030), and pair minus rank was −0.00666 (−0.01456 to +0.00103). These rank-related intervals are unadjusted descriptive summaries outside the two adjusted inferential families. They do not identify causal shares of amplitude loss, dimensional expansion, collinearity, or penalization.",
        "The rank arm provided a further comparison of ways to encode gene order (Table S10). Pair minus rank was −0.00521 (nominal conditional 95% interval, −0.01222 to +0.00181) under MAD selection and −0.00666 (−0.01456 to +0.00103) under response-associated selection. Rank-minus-expression point estimates changed from negative under MAD to positive under response-associated selection; only the former interval excluded zero. The full rank contrasts remain in Table S10. These intervals are unadjusted descriptive summaries outside the two adjusted inferential families and do not allocate causal shares to amplitude loss, dimensional expansion, collinearity, or penalization.",
        "Make the role of the rank arm explicit; retain all values in unchanged Table S10")
    change(
        "The full-development response-association pair model had an equal-drug external Spearman correlation difference of −0.000372",
        "The external comparison did not establish a relative REO benefit despite positive observed absolute ranking associations. The full-development response-association pair model had an equal-drug external Spearman correlation difference of −0.000372",
        "Distinguish external absolute association from representation advantage")
    change(
        "Drug-specific point differences were heterogeneous: −0.0237 for Telomerase Inhibitor IX, +0.0303 for BMS-536924, −0.0372 for Gefitinib, +0.0270 for Cediranib, +0.0006 for Bortezomib, and exactly zero for MIRA-1, PLX-4720, and Tamoxifen.",
        "Drug-specific differences varied in direction, with exact zeros for MIRA-1, PLX-4720, and Tamoxifen; every estimate is retained in Table 2.",
        "Avoid duplicating drug-specific point values already displayed in Table 2")
    change(
        "Absolute external ranking performance was descriptively non-zero:",
        "Observed mean external correlations were positive:",
        "Use an explicitly descriptive statement for external absolute point estimates")
    change(
        "All 40 generated cohorts and 320 cohort-by-arm summaries were retained.",
        "The constructed controls illustrated how prediction could arise under different specified signals, without assigning a cause to the real-drug results. All 40 generated cohorts and 320 cohort-by-arm summaries were retained.",
        "Introduce the scientific role of the finite controls")
    change(
        "In the tissue-only control, all arms had mean OOF correlation 0.9997 despite every molecular outer model being zero; predictions reconstructed tissue outcomes to floating-point precision. Under the weak signal, molecular mean correlations of 0.1367–0.1451 were close to tissue's 0.1373. With common sample amplitude, every rank/pair outer fit was zero and mean correlation was 0.1541, versus 0.9830 for each continuous arm. That scenario generated all 267 recorded warning occurrences (89 per continuous arm), with no recorded fallback; warning type is not recoverable from the aggregate counts. In the deliberately constructed pair-order scenario, response-associated expression/rank/pair means were 0.7573/0.7605/0.9055; MAD means were 0.2147/0.2069/0.2112, versus tissue's 0.2147 and strong 220-gene expression's 0.7572. Both driver genes were retained in all response-associated and no MAD outer candidate sets. These controls illustrate possible pipeline behavior under specified signals; they do not establish that these mechanisms caused the real-drug findings. Full metrics, generation invariants, warnings, resource observations, and provenance limits are retained in Tables S11–S12 and Additional file 14.",
        "In the tissue-only control, predictions reconstructed tissue outcomes to floating-point precision despite every molecular outer model being zero. Under the weak molecular signal, performance remained close to the tissue reference. The common-amplitude control exposed a clear difference: every rank/pair outer fit was zero, with mean correlation 0.1541 versus 0.9830 for each continuous arm. This scenario generated all 267 recorded warning occurrences (89 per continuous arm), with no recorded fallback; warning type is not recoverable from the aggregate counts.\n\nThe deliberately constructed pair-order scenario favored the response-associated pair model, whose mean correlation was 0.9055 versus 0.7573 for expression and 0.7605 for ranks. Both driver genes were retained in all response-associated and no MAD outer candidate sets; the MAD-selected arms remained close to the tissue reference. These controls illustrate possible pipeline behavior under specified signals without identifying the mechanisms of the real-drug findings. Every cohort-by-arm value, generation invariant, warning count, resource observation and provenance limit remains in Tables S11–S12 and Additional file 14.",
        "Focus synthetic prose on the constructed contrasts; keep all multi-arm metrics in unchanged Table S11 and Supplement")
    return result, edits


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "Additional_files").mkdir(exist_ok=True)
    preservation = json.loads((DOCS / "preservation_manifest.json").read_text(encoding="utf-8"))
    assert all(sha(ROOT / p) == digest for p, digest in preservation["files"].items()), "Protected source changed"
    base = OLD / "Manuscript_BMC_Genomics_author_review_v1.md"
    original = base.read_text(encoding="utf-8")
    md = original
    ab = (DOCS / "drafts/abstract_background.md").read_text(encoding="utf-8")
    dc = (DOCS / "drafts/discussion_conclusions.md").read_text(encoding="utf-8")
    for name, draft in [("Abstract", ab), ("Background", ab), ("Discussion", dc), ("Conclusions", dc)]:
        md = replace_section(md, name, section(draft, name))
    results, edits = revise_results(section(original, "Results"))
    (DOCS / "drafts/results.md").write_text("## Results\n\n" + results + "\n", encoding="utf-8", newline="\n")
    md = replace_section(md, "Results", results)
    md = re.sub(r"(?m)^\*\*Title:\*\* .+$", "**Title:** Evaluating the incremental value of relative expression ordering for drug-response prediction across eight compounds", md)
    md = re.sub(r"(?m)^\*\*Draft status:\*\* .+$", "**Draft status:** Author reading copy, 21 October 2026. Question-led narrative revision of the fixed eight-compound analysis; results and author declarations retained. Final author approval and release completion remain pending.", md)
    md = display(md)
    for name in ["Methods", "Declarations", "References", "Keywords", "Additional files", "Figure legends"]:
        assert section(md, name) == display(section(original, name)), f"Unexpected protected-section edit: {name}"
    old_tables = re.findall(r"(?m)^\|.*$", original)
    assert re.findall(r"(?m)^\|.*$", md) == [display(x) for x in old_tables]
    for label in ["Authors", "Affiliations", "Corresponding authors"]:
        assert re.search(rf"(?m)^\*\*{label}:\*\* .+$", md).group() == re.search(rf"(?m)^\*\*{label}:\*\* .+$", original).group()
    abstract = section(md, "Abstract")
    assert len(abstract.split()) <= 350
    assert re.findall(r"\*\*(Background|Results|Conclusions):\*\*", abstract) == ["Background", "Results", "Conclusions"]
    assert not re.search(r"\[\d+\]", abstract)
    main_md = OUT / "Manuscript_BMC_Genomics_author_review_v2.md"
    main_md.write_text(md, encoding="utf-8", newline="\n")
    supp_md = OUT / "Supplementary_methods_author_review_v2.md"
    shutil.copyfile(OLD / "Supplementary_methods_author_review_v1.md", supp_md)
    bib = OUT / "references_author_review_v2.bib"
    shutil.copyfile(OLD / "references_author_review_v1.bib", bib)
    builder = ROOT / "manuscript_v1/build_bmc_docx_narrative_v1.js"
    node = shutil.which("node") or r"D:\Env\NodeJS\node.exe"
    for p, embed in [(main_md, True), (supp_md, False)]:
        subprocess.run([node, str(builder), "--input", str(p), "--output", str(p.with_suffix('.docx'))] + (["--embed-figures"] if embed else []), cwd=ROOT, check=True)
    additional = []
    for n in range(1, 15):
        ext = ".docx" if n == 14 else ".csv"
        old = OLD / "Additional_files" / f"Additional_file_{n}{ext}"
        src = supp_md.with_suffix(".docx") if n == 14 else old
        dest = OUT / "Additional_files" / old.name
        shutil.copyfile(src, dest)
        additional.append({"number": n, "sha256": sha(dest), "matches_v1_bytes": sha(dest) == sha(old),
                           "science_unchanged": True})
    manifest = {"date": "2026-10-21", "status": "author_reading_copy_question_led_narrative_revision",
        "source_md_sha256": sha(base), "generator_sha256": sha(Path(__file__)), "builder_sha256": sha(builder),
        "draft_inputs": {p.name: sha(p) for p in (DOCS / "drafts").glob("*.md")},
        "documents": [{"path": p.relative_to(ROOT).as_posix(), "sha256": sha(p)} for p in [main_md, main_md.with_suffix('.docx'), supp_md, supp_md.with_suffix('.docx'), bib]],
        "abstract_whitespace_words": len(abstract.split()), "additional_files": additional,
        "methods_unchanged_except_CTRP_display": True, "declarations_unchanged_except_CTRP_display": True,
        "table_cells_unchanged_except_CTRP_display": True, "references_unchanged": True,
        "original_figure_pixels_unchanged": True, "new_scientific_analysis_or_individual_test_access": False,
        "all_authors_approved_new_copy": False, "repository_visibility_change": False,
        "results_edits": edits}
    dump(OUT / "package_manifest.json", manifest)
    # Only status, hashes, and edit purposes enter the tracked report.
    dump(DOCS / "build_summary.json", manifest)
    print(json.dumps({"main": main_md.with_suffix('.docx').relative_to(ROOT).as_posix(), "abstract_words": len(abstract.split()), "results_edit_count": len(edits)}))


if __name__ == "__main__":
    main()
