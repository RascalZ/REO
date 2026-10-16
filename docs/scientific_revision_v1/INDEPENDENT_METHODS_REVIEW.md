# Second-pass Methods alignment review — scientific revision R1

Date: 2026-10-16. Author-side working review of the scientific revision manuscript against Supplementary Methods SM1–SM6, the 34-entry Methods source ledger, frozen configuration/source and recorded execution metadata. The assistant conducting this pass also helped prepare the earlier Methods amendments; the separate comparison pass is not an independent laboratory replication, an external peer-review decision, or evidence of independent reviewer error processes.

No fitting, resampling, external scoring, sample-level response/prediction/weight access or model-pickle loading was performed. Only source/configuration, metadata and existing aggregate provenance were inspected. Main manuscript and generators were not modified by this reviewer. The authorized reporting repairs below were confined to the supplementary Methods file.

## Outcome

The inspected Methods now align with their recorded implementation and Supplementary Methods. No new result-changing numerical or information-flow error was identified within this pass. Four reporting/wording gaps identified during comparison have been closed as described below; they did not alter data, folds, models, frozen results or confidence intervals. No additional analysis is requested by this review.

## Claim-to-implementation checks

Main-manuscript line numbers below refer to the snapshot inspected during this pass; later editorial changes may move them.

| Topic | Main manuscript | Supplementary location | Source-ledger IDs / implementation | Finding |
|---|---|---|---|---|
| Cancer qualification and exact-ID exception | 63–65 | SM1 | C01–C06; cancer cohort-v2 rules and validation | Explicit Non-Cancerous exclusion, quarantine of unverified disease, HCC1588 provenance, timing after development pilot, and fixed 626/175 cohort counts are accurately bounded. No inference of complete genealogy or performance-based cancer classification. |
| SNP QC and related-group isolation | 63–65 | SM1 | C02–C06; curate_cells.py; freeze_cohorts.py | Primary SNP-matched-reference requirement and exclusion against all known GDSC2 candidate/related groups are stated. Unconfirmed/untested identities are a retained list, not a newly evaluated sensitivity cohort or proof of misidentification. |
| Historical panel and dose-window rules | 45, 55, 65; supplement supplies detail | SM1, added historical-selection paragraph | C02/C04/C05; cohort_rules.json; freeze_cohorts.py | Historical minimum of four development tissues, coverage ordering and exact dose-interval selection are now reported. Original v1 external minimum 50/preferred pool 80 is distinguished from v2 minimum 80. No current-result reranking is implied. |
| Annotation and 19,836 model genes | 51–53 | SM2 | G01–G03 | Full ENSG.version lookup, 56,238/57,820 exact reference matches, protein-coding restriction, suffix removal without averaging/collapse and zero duplicate stable identifiers agree. Unmatched audit records are not silently called noncoding. Remote TPM byte identity is not claimed. |
| Development/final global assignments | 87, 93 | SM3 | S01–S06; make_global_split_records; make_global_inner_manifest | Shuffled StratifiedGroupKFold on the sorted global 626-group table, followed by compound-specific intersection and shared arm assignments, is reported correctly. Development inner rule is repeat_seed + 100 + zero-based outer fold; final seed is 2026101501. |
| LOTO inner partitions and version history | 99–101 | SM3 | S07–S11; loto_execution_v1_3.json; compute_split | Inner partitions use each compound's remaining training rows, not the development global inner manifest. Seed is 2026101401 + sorted tissue index. The original prose/config conflict and v1.1–v1.3 execution repairs remain explicit. |
| Rare tissue warnings | 87, 93, 101 | SM3 | S03/S06/S08 and recorded warnings | Three development outer and fifteen inner warnings, one final-split warning, and 184 LOTO contexts (176 one-group/eight two-group) agree. Zero solver-warning wording is now scoped to relevant real-data stages, excluding synthetic warnings. |
| FWL and unseen tissue in original units | 81–83 | SM4 | F01–F02 | The prediction identity restores target mean/SD; unseen tissue has zero one-hot effect, not zero AUC. With beta=0, prediction equals the training tissue mean or unseen-tissue overall mean. No claim that nonzero coefficients prove utility. |
| External schema, abstention and coverage | 75, 93, 95 | SM4.1, added during review | external_validation_v1_1.json::external_gateway; external_validation_v1.py; run_external_validation_v1.py | Structural failures, full-vector molecular abstention, tissue-only exception, invalid tissue/identity, valid zero solutions, eligible/response-valid/predicted/common denominators and hard-abort rules are now locatable. |
| Definedness and conditional intervals | 105–113 | SM5 | T01–T06 | Development multiplicity/variance criterion is distinguished from the external minimum of three distinct represented groups and LOTO's two-sample descriptive helper. Shared PCG64 draws, 99% validity, linear quantiles, separate families, strict NA propagation and conditional scope align. |
| Execution environments and versions | 121–123 | SM6; Table S13 | V01–V06; stage-specific metadata/locks | Full scientific commits and recorded versions agree. Snapshot commits are not misrepresented as dispatch HEAD. LOTO threadpoolctl/joblib and the separate analyzer environment remain unrecorded rather than inferred. Cohort schema 2.0 is distinguished from LOTO execution v1.3/schema 1.0. |

## Reporting repairs completed during this pass

1. **Closed — promised external gateway detail was missing.** Main line 93 referred readers to Additional file 14 for exact schema/missingness/denominator/abort rules, but those rules were initially absent. New SM4.1 now states that missing/extra/duplicated/reordered gene columns and other structural/identity/hash failures abort evaluation. With a valid schema, invalid full-TPM values or a missing expression row trigger all seven molecular arms to abstain, including zero-coefficient arms; a valid tissue reference can still predict without molecular values. Invalid group/tissue prevents every arm. Frozen eligible, response-valid, predicted and paired-common counts are separate; invalid source experiment values invalidate the group response instead of permitting favorable-subset aggregation. No prediction or label was recalculated.
2. **Closed — real/synthetic warning scopes could be confused.** SM3 previously referred broadly to zero solver warnings in completed fitting runs. Its wording now limits that statement to the relevant real-data fitting stages and points to SM8.4 for synthetic warnings. The 267 synthetic warning occurrences remain reported; no warning was removed or retrospectively reclassified.
3. **Closed — historical metadata selection was incompletely summarized.** SM1 now records the original four-tissue minimum, descending development/external group and tissue coverage followed by ascending numeric drug-ID ordering, and within-drug exact-dose-interval selection by group coverage, model coverage and numeric concentration tie-breaks. These are existing rules, not new thresholds or current-result selections. The v2 overlay's unchanged panel/order/windows remain explicit.
4. **Closed by root integrator — outcome-standardization wording.** The earlier phrase “without ... outcome transformation” could conflict with training-only target standardization. The observed updated main line 59 now states source-scale analysis without dichotomization/response thresholds and explicitly notes inversion of training-only target standardization for evaluation. This is a wording correction, not a change to the endpoint or fitting procedure.

## Checks and boundaries

- All 34 Methods source-ledger entries matched their recorded source SHA-256 values at this pass. This is provenance consistency, not independent scientific reproduction.
- The original external NumPy/SeedSequence logging omission remains explicitly retrospective in SM5; no new field is fabricated in frozen output. Missing runtime package records remain missing, with source/environment evidence distinguished by stage.
- Conditional-OOF coverage limits, historically informed development hypotheses, prior CTRPv2 use, dose-window differences, and strict LOTO NA remain known limits. They are not relabeled as newly discovered calculation errors.
- S9 external-metric and rank-result presentation checks are handled by the root integration/result review; this Methods review does not independently validate those effect values or new graphical rendering.
- Data/model accessibility, author declarations, archive identifiers and upstream permissions remain the separately deferred publication-preparation work. Alignment of scientific Methods does not certify that the complete submission package is ready.

## Snapshot identification

| File | SHA-256 / scope |
|---|---|
| Main scientific revision at the observed wording check | `34945047a6b5ae6fad7fc4664586d3be552de15db4f5f4dddb72c38bb0cd3743` — 326 lines; root may subsequently update cross-references/layout |
| Supplement after all authorized repairs in this pass | `cd3844ae2bed77b680906988fe62302055324b2b5de1229b22b69d0363d35f13` — 208 lines |
| Methods source ledger | `4db15e3e6fba7b633a69960b76d95ad291b863de258e844a7d31ff914c449c76` — 34 source entries |
| Table S13 execution versions | `ad3ca1edc66b46544d15412cbeed7c6277a24d17d5d6d024d36b0b2255547aa8` — four stage rows |

After these reporting repairs, no unresolved substantive Methods mismatch remains in the inspected scope. Final document rendering, refreshed file hashes and cross-reference integration remain the root integrator's responsibility; no CTRPv2 rerun or new model/statistical analysis follows from this review.
