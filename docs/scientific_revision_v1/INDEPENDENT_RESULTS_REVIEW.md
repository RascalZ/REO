# R1 results source-to-text review — author-side working record

Date: 2026-10-16 (Asia/Shanghai).

This is a separate source-to-text consistency pass within the same author-authorized AI-assisted workflow. The reviewer also prepared the aggregate-table generator; it is not evidence of an independent statistical reproduction, independently fitted models, or an independent error process. The scope is the R1 abstract's nonsynthetic numerical claims, the added S9/S10 Results and Methods text, and Additional files 9/10 descriptions. Synthetic results are outside this subreview and require their separate audit. No editorial outcome is assigned.

Only the Markdown manuscript, existing public aggregate JSON/CSV/JSONL files, and the added aggregate tables were inspected. No individual response, individual prediction, model artifact, bootstrap weight/draw file, new fitting, or resampling was accessed or executed. No manuscript or generator was edited by this review task.

## Reviewed snapshot and evidence

Final checked manuscript: `manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md`, SHA-256 `24520de0dd5adcda29dc45960189c860ce4b517dbce819ffe90cd0eec54e2674`.

| Evidence file | SHA-256 | Use |
| --- | --- | --- |
| `tables/scientific_revision_v1/TableS9_external_absolute_performance.csv` | `f0ebd23139f5ef05f7516e8f77c202f7e7431e2aa79751d268929348b5c62706` | 72 absolute external Spearman rows, coverage and final zero status |
| `tables/scientific_revision_v1/TableS10_representation_comparisons.csv` | `d44a5fd9d86cba79b70a36aaf70b7d9e9095a29fbde78efd8e516ff23d677a0d` | 108 development/external matched-representation rows |
| `results/p6_external/evaluation_summary.json` | `2d3d1803fff9a631643fb0e2abff165d35fa6f88e4c33215f2366d7c7bc11acf` | Original per-arm absolute rho, common-group paired effects and existing intervals |
| `results/ordinal_repeated/paired_effect_intervals.csv` | `9b36a9ffb874e5f0a8903d72ad12a4cef51049ebdd137f6d6cc4be3ace4c8f8d` | Original development differences, interval endpoints and family labels |
| `results/ordinal_repeated/model_intervals.csv` | `6cd8c249a2774f751b5b0843ae78085966979b4841c88dd266e520df0e85853e` | Original repeat-averaged absolute development metrics |
| `results/p6_preflight/final_model_manifest.jsonl` | `0491c6faca2b79c74f5c056981ddec0263683600cbc3bb5d24c0f0bd160a86ad` | Existing final coefficient counts and audited model identities; model artifacts were not loaded |

The source/output hashes agree with `tables/scientific_revision_v1/revision_table_manifest.json`. Direct scalar comparisons used absolute tolerance 1e-12; this is a transcription/arithmetic check, not recomputation from predictions.

## Findings and resolutions

No unresolved must-fix numerical, interval-scope, or inferential-family mismatch remains in the checked final snapshot. Three issues observed during the combined review were corrected by the root agent and re-read:

| ID | Initial location | Observation and consequence | Required correction | Final verification |
| --- | --- | --- | --- | --- |
| R1-RESULT-01 | Additional file 9, line 306 | The earlier description listed RMSE, MAE and R-squared although S9 contains only external Spearman performance. Those absolute cross-screen error metrics were not evaluated. This was a substantive table-description mismatch. | Describe Spearman, counts, available existing rho intervals and final coefficient status only. | Corrected: line 306 now says absolute external **ranking** performance and Spearman; RMSE/R-squared are expressly excluded. The CSV contains no RMSE/MAE/R-squared fields. |
| R1-RESULT-02 | External rank comparison, line 181 | Response-associated external rank-minus-expression was initially displayed as -0.00871. The existing precise macro is -0.008715471658316619 and rounds to -0.00872 at five decimal places. This was a low-severity rounding error, not a result-direction change. | Use -0.00872 consistently. | Corrected and verified at line 181. |
| R1-RESULT-03 | External bootstrap description, line 109 | “Each reported estimand had 20,000 valid draws” became overbroad after S9 absolute macro and S10 external rank point-only summaries were added. These added summaries have no precomputed intervals or bootstrap counts. | Restrict the valid-draw statement to the original 27 primary/matched pair-minus-expression contrast estimands. | Corrected: line 109 explicitly names those 27 estimands. Lines 113 and 179–181 retain the no-precomputed-interval limits. |

The pre-correction snapshot used for identifying R1-RESULT-02/03 had SHA-256 `a8c85de2d4dacf030c27900d5557994d182e9609bfd84a962a07c699f31813c4`; all final status statements refer to the later checked hash above.

## Numerical and interpretation checks

### Absolute external performance, S9 and abstract

The 64 drug-arm rho values/counts match `evaluation.per_arm_descriptive`. The eight strict macro points are arithmetic means over all eight existing drug-specific rho values; none omits a drug. Macro group N=175 is the fixed external universe and the separate drug-group record N=1,185 is not an independent-subject count.

| Model | Exact strict macro rho | Abstract/Results display | Status |
| --- | ---: | ---: | --- |
| Tissue reference | 0.22840016540287353 | 0.2284 | Matches |
| Strong expression | 0.2602377014399612 | 0.2602 | Matches |
| Response-associated pair | 0.25986606179222094 | 0.2599 | Matches |

These are descriptive points, not newly tested molecular increments. The revised Results at line 179 explicitly retain that boundary, and correctly prevent a near-zero representation contrast from being described as loss of all predictive association. The three identical zero-molecular primary comparisons retain tissue predictions and remain non-informative about pair-specific molecular utility; they are not equivalence or non-inferiority tests.

All 40 nonmissing absolute-rho intervals in S9 match the original `rho_a_interval`/`rho_b_interval` entries, including their unadjusted status and original point. They cover five molecular arms at the drug level. Tissue reference, both rank arms, and absolute macro intervals remain NA; no interval is borrowed from an identical prediction or inferred from another interval. Tissue-null coefficient status is marked not applicable, not counted as a fitted molecular zero solution.

### Matched representations, S10 and Results

All 54 development rows match the original Spearman contrast point, conditional nominal endpoints, adjusted endpoints where present, interval scope and family label. All 48 external drug-level representation points match the original common-group paired correlations/differences; the six external strict macros retain all eight drugs. The original external pair-minus-expression intervals match the aggregate source. Other external representation rows have no interval.

| Candidate rule | Development contrast | Exact macro difference | Nominal conditional 95% interval | Main-text five-decimal display |
| --- | --- | ---: | --- | --- |
| MAD | Rank minus expression | -0.012438995536362463 | [-0.018786, -0.006221] | -0.01244 [-0.01879, -0.00622] |
| MAD | Pair minus rank | -0.0052062510078559705 | [-0.012215, 0.001811] | -0.00521 [-0.01222, 0.00181] |
| Response association | Rank minus expression | 0.004406144705917732 | [-0.001504, 0.010299] | 0.00441 [-0.00150, 0.01030] |
| Response association | Pair minus rank | -0.006658539602295325 | [-0.014560, 0.001026] | -0.00666 [-0.01456, 0.00103] |

The displays and the line-163 statement that these rank-related intervals are unadjusted/descriptive match S10. The existing MAD and response-associated pair-minus-expression contrasts alone retain family-18 adjustment. S10 contains 72 rows with existing nominal intervals, 36 rows with existing family-18 adjusted intervals, and 36 external rank-related point-only rows. No new interval, significance test or family was added.

| Candidate rule | External contrast | Exact strict macro difference | Correct five-decimal display |
| --- | --- | ---: | ---: |
| MAD | Rank minus expression | 0.016579198040519993 | 0.01658 |
| Response association | Rank minus expression | -0.008715471658316619 | -0.00872 |
| MAD | Pair minus rank | -0.008757215042150115 | -0.00876 |
| Response association | Pair minus rank | 0.0020030100423668526 | 0.00200 |

The revised line 181 matches these points and clearly says they have no precomputed intervals. It does not promote their point ordering to statistical superiority, equivalence or non-inferiority.

### Other nonsynthetic abstract results

- The primary development difference -0.00588 and adjusted interval [-0.01871, 0.00661] match the original development primary macro (-0.005882123341).
- The RMSE improvement -0.00075 and nominal conditional interval [-0.00124, -0.00028] match the original direction strong-expression RMSE minus pair RMSE (-0.000749305469). No direction reversal was introduced.
- The primary external difference -0.000372 and adjusted interval [-0.04217, 0.04350] match the original external primary macro (-0.000371639648).
- The 626 development groups, 175 external groups, 1,185 external drug-group responses, 430/840 molecular outer zero solutions, 29/56 final molecular zero solutions and 34/80 defined LOTO contrasts remain consistent with the existing public cohort and audit tables checked in the prior Methods pass. They were not re-estimated from individual data in this review.
- Synthetic-control abstract values/interpretations are deliberately outside this results subreview; their separate scientific audit is required.

## Remaining optional wording improvements

1. At line 113, “absolute macros” could be made “absolute external macros” to remove any possible confusion with the development absolute macros in S4, which do have existing conditional intervals.
2. At line 179, “original five-arm comparisons” could be “original comparisons collectively involving five molecular arms.” The original external inference contains three paired families, collectively covering five arms, rather than a single five-arm test.
3. The abstract's external “mean correlations” could explicitly say “equal-drug mean correlations,” consistent with the strict macro definition. The Methods and Results already define the weighting correctly.

These are clarity suggestions, not evidence of an incorrect point estimate or an unsupported new inferential claim. They do not require new fitting, resampling, or external evaluation.
