# Narrative revision: scientific-fidelity review checklist

Prepared: 21 October 2026. Planned comparison: preserved `manuscript_v1/author_review_v1/` versus new `manuscript_v1/author_review_v2/`. Purpose: support the approved research-question/evaluation-led rewrite; not a new experiment or independent replication. This checklist is prospective for the editorial comparison and does not mean v2 has passed review.

## 1. Authorized editorial scope

- [ ] Keep v1 files intact. Create v2 as a separate local author-reading package; personal author details remain outside tracked safe documentation.
- [ ] Reframe Abstract, Background, Results narration, Discussion and Conclusions around the research questions and completed comparisons. Clarity and interpretation may change; results, study design and inferential hierarchy may not.
- [ ] Normalize the displayed resource name to `CTRP v2`. This may change prose, captions and heading text only; do not alter identifiers, accession URLs, file names, checksums, source code or frozen artifact keys containing `CTRPv2`.
- [ ] Correct order-invariance wording to a common **strictly increasing** transformation applied to genes within a sample. A strictly decreasing transform reverses signs; arbitrary gene-specific transformations do not in general preserve order. Do not broaden this mathematical property to platform/batch immunity.
- [ ] Methods remains semantically identical except for the approved mathematical clarification and display-name/nonsemantic formatting changes. Any other wording difference requires line-level justification and parent review.
- [ ] Do not inspect response files, individual external predictions, model objects or checkpoints; do not fit, score, bootstrap, permute or change scientific assets.
- [ ] Preserve unrelated pre-existing worktree edits; do not use a broad cleanup or restore operation.

## 2. Study units, identity and fixed panel

- [ ] Fixed eight compounds and GDSC2 IDs remain MIRA-1 (1931), PLX-4720 (1036), Bortezomib (1191), Telomerase Inhibitor IX (1930), Tamoxifen (1199), BMS-536924 (1091), Gefitinib (1010), Cediranib (1922).
- [ ] Development: 626 independent biological groups, 4,967 compound–group responses; per-drug counts 625, 622, then six counts of 620.
- [ ] External: 175 non-overlapping groups, 1,185 compound–group responses after median aggregation of 1,201 curves. Per-drug groups: 108, 167, 169, 168, 165, 165, 161, 82. Source curves: 109, 169, 172, 170, 167, 168, 163, 83.
- [ ] Counts over drugs are not presented as additive unique cell lines. Shared CCLE expression source does not mean the same biological groups occur in training and test.
- [ ] Cohort-v2 cancer eligibility, grouped identities/related models, exact matching, SNP-reference qualification and known limits of patient genealogy remain as recorded.
- [ ] Historical CCLE 2018 TPM: 1,019 models × 57,820 source gene records; 56,238 annotation matches and 19,836 unique protein-coding genes. No DESeq2/count-based DEG pipeline is imported from the old submission.
- [ ] Outcome remains continuous fitted source AUC, without sensitive/resistant thresholds. Higher viability-curve AUC means lower in-vitro sensitivity; AUC is not ROC-AUC and not a clinical efficacy endpoint.

## 3. Model hierarchy and candidate-policy interpretation

- [ ] Primary model comparison remains response-associated 200-gene complete-pair elastic net minus the strong up-to-2,000-MAD-gene expression elastic net. Do not relabel this as a pure encoding effect, because candidate policies and budgets differ.
- [ ] Matched comparisons remain two training-only policies (MAD and tissue-centered continuous-response association) × expression200, within-sample rank200 and complete pair representations, plus tissue reference and strong expression: eight arms total.
- [ ] Response-associated candidate screening uses continuous log expression, not a representation-neutral or sequence-only rule. Any benefit conditional on this policy cannot establish that REO alone discovered the signal.
- [ ] The complete pair expansion has at most 19,900 features for 200 genes and exact-tie zero encoding. Rank and complete pair encode the same within-set order, while linear predictors and penalties differ; equal information does not require equal learned performance.
- [ ] Candidate choice, scaling, tissue projection, alpha path and tuning all remain training-only, including each inner-training fit. No SMOTE, RFE, RF, deep learning or new comparator may appear as part of the final evaluated family.
- [ ] Tissue effects remain unpenalized. Tissue-only uses observed-tissue training means and unseen-tissue overall mean; a zero molecular solution produces the same-policy tissue prediction. A zero solution is neither a failed fit nor proof of no biological signal.
- [ ] Elastic-net grid, analytic zero endpoint, pooled inner-OOF RMSE selection and tie policy remain unchanged. Selection by RMSE does not become tuning by Spearman in narrative prose.
- [ ] Three outer repeats do not become independent cohort replications. Information-flow controls do not justify claiming complete absence of all possible bias or confounding.

## 4. Development results and uncertainty

- [ ] Five outer folds × three repeats, three inner folds; 960 outer tasks, 26,160 fit events, 119,208 OOF predictions. Retain recorded rare-stratum versus solver-warning distinction.
- [ ] Primary macro Spearman: pair 0.2931, strong expression 0.2990; pair−strong −0.00588 (full-precision project value −0.005882); family-9 conditional interval −0.01871 to +0.00661.
- [ ] Six negative and two positive per-drug primary point contrasts, all adjusted drug intervals include zero. Repeat macro directions remain −0.00393, −0.00655, −0.00716; not three independent confirmations.
- [ ] RMSE values: strong 0.06217, pair 0.06292; improvement defined as strong−pair −0.00075, nominal conditional 95% interval −0.00124 to −0.00028. Preserve direction favoring expression and its secondary, unadjusted status.
- [ ] Tissue macro rho 0.2366; pair−tissue +0.05654 (+0.04299 to +0.07004) and strong−tissue +0.06242 (+0.04926 to +0.07577), nominal conditional/descriptive. These increments are shared by molecular approaches, not REO-exclusive.
- [ ] Within-policy pair−expression macros remain MAD −0.01765 (family-18 −0.03117 to −0.00406); response association −0.00225 (−0.01438 to +0.00949).
- [ ] Rank-related development descriptive intervals remain unadjusted: MAD rank−expression −0.01244 [−0.01879, −0.00622], pair−rank −0.00521 [−0.01222, +0.00181]; response rank−expression +0.00441 [−0.00150, +0.01030], pair−rank −0.00666 [−0.01456, +0.00103].
- [ ] Changes in candidate-policy contrast may be described as policy-dependent representation performance; do not claim a newly tested interaction, causal decomposition, or quantitative attribution to information loss/collinearity/penalties.
- [ ] Primary family has nine contrasts; matched-pair family has 18; separate Bonferroni families are not study-wide simultaneous coverage. Keep conditional bootstrap meaning and point/CI rounding consistent across text/tables.
- [ ] Development uses 10,000 shared biological-group resamples, external 20,000. Intervals condition on saved OOF or fixed models and do not repeat feature selection, tuning or fitting. NA, validity thresholds and no replacement of invalid draws are retained.
- [ ] No non-inferiority margin or equivalence/non-inferiority test has been adopted. Non-significance and algebraic zero intervals are not equivalence, information preservation, or non-inferiority evidence.

## 5. Zero solutions and LOTO

- [ ] Real development molecular zeros remain 430/840 outer fits; final molecular zeros 29/56, within 64 total final arms. Do not mix outer-fit frequency with final-model status or include explicit nulls as molecular zeros.
- [ ] Final fit events 1,744; final 64 models are unchanged. Nonzero coefficients alone do not establish predictive utility, mechanism or molecular resistance signatures.
- [ ] LOTO: 23 tissues, 184 drug–tissue contexts; primary display ten tissues × eight drugs, only 34/80 contrasts defined. Strict eight-drug macro remains NA, not zero or the mean of 34 available cells.
- [ ] Only two ten-tissue strict drug means are defined: Telomerase Inhibitor IX −0.0124 and BMS-536924 −0.0526. Do not interpret missing correlations as demonstrated failure to transfer or as successful generalization.
- [ ] LOTO remains descriptive without a new inferential family/CI. Failed/repaired execution history and rare-stratum warnings remain reported in supplementary records.

## 6. Locked external evaluation and claims

- [ ] One completed locked retrospective GDSC2→CTRP v2 evaluation, not an unused future holdout. Prediction lock preceded response access by the reconstructed workflow; historical manuscript exposure remains disclosed.
- [ ] All 9,480 planned predictions finite, zero abstentions; this is coverage, not predictive efficacy. Structural mismatches remain hard failures; missing/invalid molecular input policy is unchanged.
- [ ] Primary external macro pair−strong −0.000372 (full precision −0.00037164); family-9 conditional interval rounded −0.04217 to +0.04350 (full precision −0.042166 to +0.043505).
- [ ] Per-drug primary contrasts retain sign/magnitude and Table 2: approximately −0.02369 Telomerase, +0.03033 BMS, −0.03721 Gefitinib, +0.02700 Cediranib, +0.00060 Bortezomib; exact zero MIRA-1, PLX-4720, Tamoxifen.
- [ ] Those three exact zeros arise from identical same-policy tissue predictions after both molecular models select zero coefficients; they are not evidence of model equivalence under a statistical test.
- [ ] External matched pair−expression: MAD +0.00782 [−0.01884, +0.03450]; response association −0.00671 [−0.02562, +0.01369], family-18 conditional intervals.
- [ ] Absolute external macro correlations remain tissue 0.2284, strong 0.2602, pair 0.2599. No new interval exists for these macros; a near-zero pair−strong contrast is not disappearance of all transferred association or a formally established external molecular increment.
- [ ] External rank contrasts remain descriptive points only: rank−expression +0.01658 MAD/−0.00872 response; pair−rank −0.00876 MAD/+0.00200 response. No interval is invented, borrowed or obtained by subtracting endpoints.
- [ ] Same CCLE measurement source, different assay/dose windows; target is within-drug ranking across screens, not expression-platform robustness, absolute AUC calibration, prospective blind prediction or clinical translation.
- [ ] New narrative must not turn preserved negative results into a claim that REO is universally useless, nor turn finite observed increments into broad REO superiority.

## 7. Synthetic controls and adjacent-field context

- [ ] All four designed scenarios × ten generated cohorts, 40 cohorts and 320 cohort-arm records remain. Each cohort has 240 groups, eight tissues × 30 and 220 genes.
- [ ] Synthetic strong arm uses 220 genes in fact, not 2,000; matched arms use 200. Common amplitude rescales TPM+1, preserving within-sample order by construction.
- [ ] Synthetic zeros 615/840 selected outer molecular refits: 210/210 tissue-only, 198/210 weak, 120/210 amplitude, 87/210 predefined order. Do not include intentionally evaluated inner analytic null candidates in this denominator.
- [ ] Keep all synthetic scenario values reported in v1, including tissue-only near-perfect 0.9997 despite zero molecular models; amplitude order models 0.1541 versus continuous 0.9830; response-associated constructed pair-order expression/rank/pair 0.7573/0.7605/0.9055.
- [ ] All 267 warning occurrences in amplitude (89 per continuous arm), unknown warning type from aggregates, no recorded fallback, undefined residuals left blank and original synthetic environment gaps remain.
- [ ] These are finite, intentionally specified implementation illustrations; not causal identification of real-drug zero mechanisms, a representation-neutral benchmark, a power analysis or replacement for the unperformed new eight-drug full-pipeline label permutation.
- [ ] AIVC/single-cell perturbation context stays short and balanced, acknowledging competitive linear baselines and gains on calibrated metrics. No direct test of a specific virtual-cell model, universal collapse mechanism, novelty-first or proven governance-template claim is introduced.

## 8. Author facts, access, references and assets

- [ ] Title-page author roster/order, equal-contribution/correspondence markers and affiliation remain exactly as v1; do not place these personal details in this safe checklist or tracked review.
- [ ] Declaration facts remain v1: interests and both funding sources/identifiers, known absence of new recruitment/material collection, actual AI involvement and human responsibility. Do not invent emails, author roles, acknowledgement permissions, ethics waivers, license approval, archive DOI or all-author final approval.
- [ ] Project repository remains private as verified 20 October 2026, distinct from public upstream resources. Code/aggregate/example tracking is not verified anonymous public/reviewer access; complete derived model/mapping assets remain local. No promise of access or permission is added.
- [ ] Existing 31 reference identities, numeric mapping, DOIs/URLs and author formatting remain unless parent explicitly approves a separately verified change. Narrative citations must still support their new adjacent claims; no citation becomes decoration for an untested mechanism.
- [ ] Tables 1/2 retain all numeric cells, directions, intervals, denominators and strict macros. Pure `CTRPv2`→`CTRP v2` display normalization in headers/captions is allowed; no numeric cell changes.
- [ ] Existing four PNG/SVG figure assets remain byte-identical, as do their plotted values/matrices. Figure 3B is final-model status, not outer-fold zero frequency; Figure 4A encodes definedness, not effect sign or magnitude.
- [ ] Additional files 1–13 remain byte-identical, including blanks and NA. Additional file 14 / Supplementary methods content must be unchanged except approved display-name/math/format normalization; if copied unchanged, its file bytes should remain identical. Rebuilt DOCX ZIP bytes can change only when rebuild is declared and source-content comparison is documented.
- [ ] Main figure legends, Additional-file descriptions and table notes remain consistent with source assets. New Word embedded image bytes match the unchanged PNGs exactly.
- [ ] Markdown, Word and rendered PDF belong to the same v2 source. Rebuild and actually inspect rendered figure/table/reference pages; verify page sequence, clipping, readable figures and no duplicate captions. Structural checks alone do not establish visual completion.

## 9. Review record to append after v2 generation

Record exact v1/v2 source hashes; sections changed; the allowed Methods/name/math diff; numeric-token and table-cell comparisons; unchanged figure/Additional-file hashes; reference/citation consistency; author/declaration preservation without copying personal data; concrete unresolved scientific/editorial points. Mark results as passed, failed or not assessed with evidence locators. Do not use a passed count as a certificate of scientific validity, independent reproduction, journal compliance or submission authorization.

**Current status: checklist prepared; v2 independent comparison pending notification from the coordinating agent.**
