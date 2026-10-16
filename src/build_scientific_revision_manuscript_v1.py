"""Build R1 from the preserved BMC working draft; no scientific recomputation.

This is a deterministic editorial transform, not a model/evaluation entry point.
Read access is limited to the draft, numeric-citation helper and journal map.
The separate aggregate-table builder records the new tables' scientific sources.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import re
from pathlib import Path

from validate_numeric_references import expand_numeric_group

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "manuscript_v1/Manuscript_BMC_Genomics_draft.md"
OUTPUT = ROOT / "manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md"
AUDIT = ROOT / "docs/scientific_revision_v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    baseline_hash = sha(BASE)
    text = BASE.read_text(encoding="utf-8")

    def replace_paragraph(prefix: str, replacement: str) -> None:
        nonlocal text
        matches = [p for p in text.split("\n\n") if p.startswith(prefix)]
        if len(matches) != 1:
            raise ValueError(f"Expected one paragraph: {prefix!r}; found {len(matches)}")
        text = text.replace(matches[0], replacement, 1)

    text = text.replace("# BMC Genomics Research Article — working draft", "# BMC Genomics Research Article — scientific revision R1", 1)
    replace_paragraph("**Draft status:**", "**Draft status:** Scientific revision R1, 16 October 2026. This is the current content-review copy; the preceding evidence and BMC drafts remain unchanged. Author declarations, formal archiving, and upstream-license confirmation are deferred until scientific content is near final. Author-input placeholders remain; this is not a submission-ready manuscript.")
    replace_paragraph("**Background:**", "**Background:** The incremental value of within-sample relative expression ordering (REO) for continuous drug-response prediction remains uncertain. We compared a tissue reference, a strong 2,000-gene expression model, and matched 200-gene expression, rank, and full-pair models under two training-fold candidate rules. Eight fixed compounds and 626 development cell-line groups underwent five-fold outer cross-validation repeated three times with three-fold inner tuning. Full-development models underwent one locked retrospective CTRPv2 evaluation in 175 non-overlapping groups. Four finite synthetic controls examined specified implementation behavior.")
    replace_paragraph("**Results:**", "**Results:** In development, the predefined REO-pair minus strong-expression difference in equal-drug mean Spearman correlation was −0.00588 (family-9 adjusted conditional bootstrap interval, −0.01871 to 0.00661). The RMSE improvement, strong expression minus pair, was −0.00075 (nominal conditional 95% interval, −0.00124 to −0.00028), favoring expression. In the locked CTRPv2 evaluation, mean correlations were 0.2284 for tissue, 0.2602 for strong expression, and 0.2599 for REO pair; the primary difference was −0.000372 (family-9 adjusted interval, −0.04217 to 0.04350; 1,185 drug–group responses). Matched-budget comparisons did not establish an average pair advantage. Zero molecular coefficients occurred in 430 of 840 outer fits and 29 of 56 final molecular arms. The strict leave-one-tissue-out macro was undefined: only 34 of 80 prespecified contrasts had defined correlations. Synthetic results illustrated tissue-only prediction, weak-signal shrinkage, order-invariant amplitude loss, and a constructed pair-order signal; they did not identify real-data mechanisms.")
    replace_paragraph("**Conclusions:**", "**Conclusions:** The evaluated REO pipeline did not show a reliable average ranking improvement over strong expression in this eight-compound setting. This does not establish equivalence or non-inferiority or exclude other REO formulations. Reporting tissue contribution, candidate policy, representation contrasts, zero molecular solutions, and assay context is necessary to delimit such comparisons.")
    # Clarify the contribution relative to existing evaluation literature.
    old = "These concerns are especially relevant when feature selection itself is response-associated."
    text = text.replace(old, old + " Prior benchmarks already motivate careful baselines and explicit prediction targets [5,6]. Our narrower contribution is a matched candidate-policy-by-representation comparison under the same tissue-adjusted learning rule, with audits of when a nominal molecular predictor reduces to its covariate reference.", 1)
    replace_paragraph("Machine-readable supporting tables", "Machine-readable supporting tables are supplied as Additional files 1–13 (Tables S1–S13). Additional file 14 provides supplementary methods and complete descriptive synthetic-control results.")
    expr = "Expression gene records were annotated against the CCLE-matched gencode.v19 reference using complete versioned Ensembl identifiers: 56,238 of 57,820 source records matched. Only 19,836 explicitly protein-coding records entered modeling. The terminal numeric version was removed for stable feature identifiers without collapsing or averaging records; no duplicate resulting identifiers were found. Unmatched records were annotation-unresolved, not assigned a noncoding biotype. The reference filename, checksum, and annotation checks are in Additional file 14. The local source audit did not re-download the complete historical TPM matrix to establish byte-for-byte remote identity."
    text = text.replace("retained on the source scale without dichotomization, thresholding, or outcome transformation.", "analyzed on the source AUC scale without dichotomization or response thresholds. Training-only target standardization was inverted for evaluation.", 1)
    text = text.replace("Eight compounds were fixed by their GDSC2 identifiers:", expr + "\n\nEight compounds were fixed by their GDSC2 identifiers:", 1)
    replace_paragraph("Cell-line identifiers were reconciled", "Cell-line identities were reconciled by authoritative exact CCLE, Broad/DepMap, Sanger, and COSMIC cross-references; ambiguous identities were not resolved by name similarity. Common biological groups used Broad same-tissue-sample PatientID, Sanger patient links, and known parent/derivative relationships. Broad PatientID does not establish a complete patient genealogy. The cancer-focused cohort-v2 rules excluded explicit Non-Cancerous annotations and quarantined missing or Unknown disease annotations unless exact-ID cancer provenance was established. One exception, HCC1588 (ACH-001078; CVCL_A351; SIDM01820), was retained with concordant cancer-cell-line evidence. This metadata amendment followed a development-only pilot and preceded external response access by the reconstructed workflow. Tissue names were canonicalized by exact case-insensitive vocabulary matching.\n\nThe primary CTRPv2 cohort required SNP-matched-reference identity status and excluded all known GDSC2 and related biological groups, including candidate identities in unresolved mappings, rather than only groups observed for the eight drugs. A separate list of unconfirmed/untested SNP identities does not constitute an evaluated sensitivity cohort or evidence of misidentification. The fixed panel met cohort-v2 feasibility minima of 200 development and 80 primary external groups per drug; these are not power guarantees. The original cohort-v1 minimum was 50 external groups with a preferred 80-group pool, from which all eight were selected. Cohort-v2 yielded 626 development groups and 4,967 drug–group responses, and 175 external groups and 1,185 responses after prespecified within-compound/group median aggregation. Counts are in Table 1; exact identity and exclusion rules are in Additional file 14.")
    text = text.replace("Counts are in Table 1; exact identity and exclusion rules are in Additional file 14.", "Counts are in Table 1; exact identity and exclusion rules, historical coverage ordering, and per-compound dose-window selection are in Additional file 14.", 1)
    replace_paragraph("For each molecular arm, features", "For each molecular arm, features were constructed before standardization and tissue adjustment. A training-only standard scaler was fit to molecular features, and the current training outcome was standardized. Tissue one-hot effects were unpenalized. The molecular design and outcome were projected onto the orthogonal complement of the training tissue design by Frisch–Waugh–Lovell residualization; an elastic net without intercept was fit to these residuals. Prediction restored the training target mean and scale after combining the molecular and tissue components. For an unseen tissue the one-hot vector was zero: prediction was the training overall mean plus any retained molecular contribution, not a zero AUC. With zero molecular coefficients, this matched the tissue reference's observed-tissue mean/unseen-tissue overall-mean policy. The prediction identity is given in Additional file 14.")
    replace_paragraph("The development evaluation used five outer folds", "Development used five outer folds repeated three times, with three inner folds. StratifiedGroupKFold(shuffle=True) used tissue strata and biological-group identities on the canonical-model-sorted 626-group table. Outer seeds were 20261013, 20261014, and 20261015; inner seeds were repeat_seed + 100 + zero-based outer fold. The global assignments were intersected with each compound's eligible records and shared by all arms. Group separation and exactly-once validation coverage were checked. Rare strata of one or two groups produced three outer and fifteen inner split-preparation warnings; frozen partitions were retained without score-based exclusions or reassignment. These are distinct from solver warnings. Candidate selection, feature construction, scaling, tissue adjustment, and alpha_max were rebuilt in every inner-training fit. The lowest pooled inner-OOF RMSE selected the hyperparameters, followed by one refit on all outer-training data and one held-out evaluation. The tissue reference required one fit and no tuning.")
    replace_paragraph("After the development comparison was complete", "After the development comparison, 64 full-development models were trained with a shared three-fold StratifiedGroupKFold assignment on the sorted 626-group table (seed 2026101501), intersected with each drug's records. One rare-stratum warning was recorded. The same grid, training-only feature pipeline, and pooled inner-OOF RMSE objective were retained. The 1,744 fit events included inner candidates and all-development refits. Artifacts were reloaded in an independent process and checked against locked source and prediction audits. Response-free CTRPv2 predictions were generated and locked before the reconstructed evaluation workflow accessed external response values. Structural schema mismatches were hard failures; invalid individual molecular inputs triggered abstention rather than imputation. Exact schema, missingness, denominator, and abort rules are in Additional file 14.")
    # Append split/version history without deleting the original LOTO estimand description.
    loto = next(p for p in text.split("\n\n") if p.startswith("Each of the 23 metadata-defined tissues"))
    text = text.replace(loto, loto + "\n\nHeld-out tissue names were sorted by Python Unicode order and indexed from zero; three-fold inner tuning used seed 2026101401 + tissue_index on each drug's remaining groups. Across 184 drug–tissue contexts, all emitted rare-stratum warnings (176 with one group; eight with two); group separation and validation coverage were verified. Before fitting, execution amendment v1.1 resolved a prose/config seed conflict in favor of the frozen machine value 2026101401. Subsequent logging and row-identity repairs produced completed execution v1.3 without changing models, seeds, grid, estimands, or scoring. The supplementary methods preserve the failed/repaired execution history.", 1)
    # Definedness rules are stage-specific; do not silently impose external criteria on development.
    stats_end = "No new full-pipeline permutation test was run in this eight-drug analysis."
    text = text.replace("Each reported estimand had 20,000 valid draws.", "Each of the 27 prespecified primary and matched pair-minus-expression contrast estimands had 20,000 valid draws.", 1)
    text = text.replace(stats_end, stats_end + "\n\nDevelopment used shared integer group-multiplicity draws from Generator(PCG64), seed 2026101602; tied midranks were recomputed under multiplicities, requiring total multiplicity at least two and nonconstant weighted ranks. External resampling used seed 2026101502 and required at least three distinct represented groups and nonconstant response/prediction for each correlation. Intervals required a defined point and at least 99% valid draws; no invalid draw was replaced. Percentiles used NumPy quantile(method='linear'), with adjusted tail probability 0.05/(2K). LOTO used a separate descriptive helper, requiring at least two held-out samples and nonconstant response and prediction. No undefined component was replaced by zero or omitted from a strict macro. Additional file 14 specifies the full conditional estimands and provenance limits.\n\nAdditional files 9 and 10 consolidate existing absolute external metrics and matched representation contrasts. They add no resampling or inferential family. Existing absolute-rho intervals are available for five molecular arms at the drug level only. Development rank contrasts retain existing nominal, unadjusted descriptive intervals; external rank contrasts and absolute macros have no precomputed intervals. Missing intervals remain NA, including tissue-reference and rank-arm absolute intervals. Their absence does not mean zero uncertainty.", 1)
    replace_paragraph("Analyses were governed by versioned protocols", "Analyses were governed by versioned protocols, configurations, source hashes, cohort/execution locks, and result manifests. Scientific source snapshots were 5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43 (repeated development), fedc745940374b6870cb6c4bcb66ec8f731ded27 (completed LOTO v1.3), and ca0f837ac10665ab8d4d6ab97691e58e5bed4577 (final models/external execution); these bind scientific files rather than necessarily the mutable dispatch HEAD. All three fitting stages recorded Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, SciPy 1.16.2, scikit-learn 1.7.2, and psutil 7.1.0. Stage-specific recorded and unrecorded fields are explicit in Table S13. In particular, the development-bootstrap analyzer environment and actual synthetic-run package snapshot were not separately logged and are not inferred from the plotting or current audit environment.\n\nThe repository contains source versions, download/accession provenance, code, aggregate results, figure/table scripts, checksums, and recorded environment metadata. Sample-level CTRPv2 responses, group bootstrap weights, individual predictions, and fitted checkpoints are excluded from the public repository. Reproduction from source data requires the documented identity/annotation and environment records; current public checks verify aggregate consistency rather than independently reproducing fitted predictions. Upstream terms remain applicable and are not replaced by the repository code license. Final dataset access statements, release identifiers, and permitted data/model-access arrangements remain to be completed before submission.")
    synth_methods = "### Finite synthetic implementation controls\n\nFour prespecified data-generating processes—tissue-only outcome, weak molecular signal, common sample amplitude, and predefined pair order—were each evaluated in ten separately generated cohorts. Each cohort had 240 independent synthetic groups, eight tissues of 30 groups, and 220 genes. The common-amplitude design held all within-sample orders identical; the pair-order design deliberately assigned signal to two lower-variance driver genes. The same eight arms used three outer folds and three inner folds, with training-only preprocessing and the same relative-alpha grid. Matched arms used 200 genes, while the nominal 2,000-gene strong arm actually used all 220 available genes. All 40 cohorts were retained. Cohort metrics and generation checks are in Tables S11 and S12; equations, RNG rules, execution provenance, and all scenario results are in Additional file 14. These were descriptive controls without inferential intervals, not a representation-neutral benchmark, real-data mechanism identification, or a substitute for full-pipeline real-label permutation. Actual synthetic runtime package versions were not separately recorded."
    text = text.replace("### Reproducibility and data governance", synth_methods + "\n\n### Reproducibility and data governance", 1)
    # Main Results: retain original adjusted families; explicitly show ranks too.
    rank_dev = "The rank arm further distinguishes order encoding from the complete pair expansion (Table S10). Under MAD selection, rank minus expression was −0.01244 (existing nominal conditional 95% interval, −0.01879 to −0.00622), and pair minus rank was −0.00521 (−0.01222 to +0.00181). Under response-associated selection, rank minus expression was +0.00441 (−0.00150 to +0.01030), and pair minus rank was −0.00666 (−0.01456 to +0.00103). These rank-related intervals are unadjusted descriptive summaries outside the two adjusted inferential families. They do not identify causal shares of amplitude loss, dimensional expansion, collinearity, or penalization."
    text = text.replace("### Zero molecular solutions and fitted-model status", rank_dev + "\n\n### Zero molecular solutions and fitted-model status", 1)
    absolute_ext = "Absolute external ranking performance was descriptively non-zero: the strict equal-eight-drug mean correlation was 0.2284 for tissue, 0.2602 for strong expression, and 0.2599 for response-associated pair (Table S9). The near-zero primary contrast therefore does not show that all transferable predictive association disappeared. These absolute macros have no precomputed intervals, and the point difference from tissue is not a new inferential molecular-increment result. Table S9 includes every drug and arm, correlation, sample counts, and final zero-coefficient status; available drug-level absolute-rho intervals are transcribed from the original five-arm comparisons. No interval was derived by subtracting endpoints or borrowed for identical predictions. Uncalibrated cross-screen RMSE and R-squared were not evaluated because the assay AUC scales and dose windows differ.\n\nExternal rank minus expression was +0.01658 under MAD and −0.00872 under response-associated selection; pair minus rank was −0.00876 and +0.00200, respectively. These are descriptive strict macro points without precomputed intervals (Table S10). Together with the adjusted pair-minus-expression results, they do not establish general pair superiority, equivalence, or non-inferiority."
    text = text.replace("Table 2. Primary REO-pair", absolute_ext + "\n\nTable 2. Primary REO-pair", 1)
    synth_results = "### Finite synthetic-control results\n\nAll 40 generated cohorts and 320 cohort-by-arm summaries were retained. Across the 840 selected molecular outer refits, 615 were zero solutions: 210/210 in tissue-only, 198/210 in weak molecular, 120/210 in common amplitude, and 87/210 in predefined pair order. This denominator excludes inner-grid fits, including deliberately evaluated analytic null candidates. Undefined residual correlations remain blank in Table S11 rather than being replaced by zero.\n\nIn the tissue-only control, all arms had mean OOF correlation 0.9997 despite every molecular outer model being zero; predictions reconstructed tissue outcomes to floating-point precision. Under the weak signal, molecular mean correlations of 0.1367–0.1451 were close to tissue's 0.1373. With common sample amplitude, every rank/pair outer fit was zero and mean correlation was 0.1541, versus 0.9830 for each continuous arm. That scenario generated all 267 recorded warning occurrences (89 per continuous arm), with no recorded fallback; warning type is not recoverable from the aggregate counts. In the deliberately constructed pair-order scenario, response-associated expression/rank/pair means were 0.7573/0.7605/0.9055; MAD means were 0.2147/0.2069/0.2112, versus tissue's 0.2147 and strong 220-gene expression's 0.7572. Both driver genes were retained in all response-associated and no MAD outer candidate sets. These controls illustrate possible pipeline behavior under specified signals; they do not establish that these mechanisms caused the real-drug findings. Full metrics, generation invariants, warnings, resource observations, and provenance limits are retained in Tables S11–S12 and Additional file 14."
    text = text.replace("## Discussion", synth_results + "\n\n## Discussion", 1)
    text = text.replace("In the locked CTRPv2 evaluation, mean correlations were", "In the locked CTRPv2 evaluation, equal-drug mean correlations were", 1)
    text = text.replace("external rank contrasts and absolute macros have no precomputed intervals", "external rank contrasts and absolute external macros have no precomputed intervals", 1)
    text = text.replace("the original five-arm comparisons", "the original paired comparisons collectively involving five molecular arms", 1)
    old = "Their reported tasks and estimands differ from the present continuous in-vitro response-ranking task, so their performance values are not direct comparators."
    text = text.replace(old, old + " In particular, sparse top-scoring-pair classification and patient ex-vivo response signatures [21] should not be treated as the same model as a complete 19,900-pair elastic-net regression.", 1)
    text = text.replace("earlier synthetic mechanisms test software behavior under defined scenarios but do not provide an empirical null distribution for these drug outcomes.", "the reported finite synthetic controls illustrate software behavior under intentionally specified scenarios but do not provide an empirical null distribution or identify the mechanisms of these drug outcomes. The common-amplitude control retained 267 warning occurrences, and the original synthetic runtime package snapshot was not separately logged.", 1)
    text = text.replace("add the archived analysis-code release DOI and exact repository tag after release", "add a persistent archived analysis-code release identifier and exact repository tag after release", 1)
    text = text.replace("The repository URL alone is not a versioned archive DOI.", "The repository URL alone does not identify an immutable archived version.", 1)
    extra = """**Additional file 9.** `Additional_file_9.csv` (.csv); Table S9: Absolute locked external ranking performance for all eight drugs and eight arms, plus strict eight-drug macros (72 rows). Includes Spearman, counts and final coefficient status. Existing drug-level absolute-rho intervals are available for five molecular arms; all other uncomputed intervals remain NA. Uncalibrated cross-screen RMSE and R-squared are not included.

**Additional file 10.** `Additional_file_10.csv` (.csv); Table S10: Expression200, rank200 and complete-pair comparisons under both candidate rules, for development and external evaluation (108 rows). Only pair-minus-expression belongs to the existing family-18 adjustment. Development rank contrasts have existing nominal descriptive intervals; external rank contrasts have point estimates only. No new inference or intervals were computed.

**Additional file 11.** `Additional_file_11.csv` (.csv); Table S11: Complete finite synthetic cohort-by-arm results (320 rows). Byte-identical to the recorded cohort metrics; blank undefined/not-applicable residual correlations remain blank. All four scenarios and all cohorts are retained; no inferential intervals are reported.

**Additional file 12.** `Additional_file_12.csv` (.csv); Table S12: All 40 synthetic generation checks, including RNG seeds, dimensions, ranges, scenario-specific invariants, and cohort/split hashes. The original log2_tpm column names refer to L=log2(TPM+1). Invariant fields outside their applicable scenarios remain blank.

**Additional file 13.** `Additional_file_13.csv` (.csv); Table S13: Recorded execution versions for repeated development, LOTO, full-development fitting and the locked external gateway. Missing original environment fields are explicitly marked not_recorded; scientific source bindings are distinguished from dispatch HEAD.

**Additional file 14.** `Additional_file_14.docx` (.docx); Supplementary methods and synthetic-control results. Identity/annotation rules, group splits and warnings, prediction identity, external missingness/abort policy, conditional interval implementation, execution amendments, and finite synthetic equations, provenance, warnings, resource records and all scenario results.

"""
    text = text.replace("## Figure legends", extra + "## Figure legends", 1)
    # Renumber only this new manuscript in first-appearance order.
    body, remainder = text.split("## References\n", 1)
    old_refs, tail = remainder.split("\n## Additional files", 1)
    references = {int(n): r for n, r in re.findall(r"(?m)^(\d+)\.\s+(.+)$", old_refs)}
    order: list[int] = []
    for group in re.findall(r"\[([^\]]+)\]", body):
        expanded = expand_numeric_group(group, len(references))
        if expanded:
            for n in sorted(expanded):
                if n not in references:
                    raise ValueError(f"Unresolved reference {n}")
                if n not in order:
                    order.append(n)
    if set(order) != set(references):
        raise ValueError("Uncited references")
    mapping = {old: new for new, old in enumerate(order, 1)}

    def recite(match: re.Match[str]) -> str:
        expanded = expand_numeric_group(match.group(1), len(references))
        if not expanded:
            return match.group(0)
        return "[" + ",".join(str(mapping[n]) for n in sorted(expanded, key=lambda n: mapping[n])) + "]"

    body = re.sub(r"\[([^\]]+)\]", recite, body)
    refs = "\n".join(f"{mapping[n]}. {references[n]}" for n in order)
    text = body + "## References\n\n" + refs + "\n\n\n## Additional files" + tail
    OUTPUT.write_text(text, encoding="utf-8", newline="\n")
    AUDIT.mkdir(parents=True, exist_ok=True)
    map_path = ROOT / "manuscript_v1/reference_number_mapping_scientific_revision_v1.csv"
    with map_path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["previous_BMC_reference", "scientific_revision_R1_reference", "reference"])
        writer.writeheader()
        writer.writerows({"previous_BMC_reference": n, "scientific_revision_R1_reference": mapping[n], "reference": references[n]} for n in order)
    old_nlm = ROOT / "manuscript_v1/nlm_journal_abbreviations.csv"
    rows = list(csv.DictReader(old_nlm.read_text(encoding="utf-8-sig").splitlines()))
    for row in rows:
        row["reference_number"] = str(mapping[int(row["reference_number"])])
    rows.sort(key=lambda r: int(r["reference_number"]))
    new_nlm = ROOT / "manuscript_v1/nlm_journal_abbreviations_scientific_revision_v1.csv"
    with new_nlm.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    assert sha(BASE) == baseline_hash
    manifest = {"revision": "scientific_R1", "date": "2026-10-16", "operation": "editorial_transform_only_no_model_evaluation_or_resampling", "baseline": {"path": str(BASE.relative_to(ROOT)), "sha256": baseline_hash, "preserved": True}, "generator_sha256": sha(Path(__file__)), "outputs": [{"path": str(p.relative_to(ROOT)), "sha256": sha(p)} for p in [OUTPUT, map_path, new_nlm]], "references": len(references), "first_appearance_mapping": mapping}
    (AUDIT / "manuscript_build_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT.relative_to(ROOT)), "words": len(text.split()), "references": len(references), "baseline_preserved": True}))


if __name__ == "__main__":
    main()
