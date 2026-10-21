## Results

### A common cohort supported paired model comparisons

The development cohort comprised 626 unique biological groups and 4,967 drug–group response observations. The per-compound development sample counts ranged from 620 to 625 groups. Counts across compounds are not additive unique-cell-line counts because the same biological group can contribute responses to multiple compounds. The CTRPv2 cohort contained 175 non-overlapping external groups and 1,185 valid drug–group responses after median aggregation of 1,201 source curve records. The number of common external groups for the primary paired comparison ranged from 82 for Cediranib to 169 for Bortezomib; every drug met the prespecified minimum of 80. The complete drug identity and cohort counts are in Table 1 and Additional file 8 (Table S8).

Table 1. Drug identities and cell-line group counts in the development and external screens. Per-drug counts refer to eligible GDSC2 development groups and CTRPv2 groups with valid aggregated drug-response measurements. Eligible CTRPv2 curve records are shown separately from the per-cell-line response count because replicate curves were aggregated by median. The final row reports across-drug unique cell-line groups separately from drug–model response records; per-drug counts must not be summed and interpreted as unique cell lines.

| Compound | GDSC2 ID | Development groups / responses | CTRPv2 ID | External groups / responses / source curves |
| --- | --- | --- | --- | --- |
| MIRA-1 | 1931 | 625 / 625 | 580922 | 108 / 108 / 109 |
| PLX-4720 | 1036 | 622 / 622 | 606256 | 167 / 167 / 169 |
| Bortezomib | 1191 | 620 / 620 | 411737 | 169 / 169 / 172 |
| Telomerase Inhibitor IX | 1930 | 620 / 620 | 411865 | 168 / 168 / 170 |
| Tamoxifen | 1199 | 620 / 620 | 26972 | 165 / 165 / 167 |
| BMS-536924 | 1091 | 620 / 620 | 411859 | 165 / 165 / 168 |
| Gefitinib | 1010 | 620 / 620 | 52926 | 161 / 161 / 163 |
| Cediranib | 1922 | 620 / 620 | 606136 | 82 / 82 / 83 |
| All eight compounds | — | 626 / 4967 | — | 175 / 1185 / 1201 |

*Note:* Development groups and external groups are independent biological groups; response-record counts are drug–model observations. Per-drug counts overlap across compounds and must not be summed as unique cell lines. CTRPv2 source curves are shown separately from median-aggregated responses.

The CTRPv2 evaluation generated finite predictions for all 9,480 planned compound-by-model rows, with zero abstentions. The external comparison therefore had complete prediction coverage under the frozen eligibility rules.

### Pair encoding did not establish an average gain over strong expression

The primary comparison did not establish a ranking benefit from the evaluated REO pipeline. The equal-drug macro Spearman correlation was 0.2931 for the response-association pair model and 0.2990 for the strong 2,000-gene expression model. The prespecified difference, pair minus strong expression, was −0.00588 (family-9 Bonferroni-adjusted conditional interval, −0.01871 to 0.00661). Six of eight drug-specific point estimates were negative and two were positive; every drug-specific family-adjusted interval included zero. The macro difference was negative in all three repeated partitions. This describes consistency across partitions, not independent replication; repeat-level values remain in Additional file 5 (Table S5). The paired drug-wise and equal-drug estimates are shown in Figure 2.

For the macro RMSE, the strong-expression model had 0.06217 and the REO-pair model had 0.06292. The paired improvement, defined as strong-expression RMSE minus pair RMSE, was −0.00075 (unadjusted conditional 95% interval, −0.00124 to −0.00028), favoring the strong-expression model on this error measure. 

The observed development increment over tissue was shared by both primary molecular models. The tissue-only macro Spearman was 0.2366; the descriptive response-association pair minus tissue-only difference was +0.05654 (unadjusted conditional interval, +0.04299 to +0.07004). The strong-expression minus tissue-only difference was +0.06242 (+0.04926 to +0.07577). These tissue-reference contrasts indicate that both molecular modeling pipelines can add ranking information over the tissue reference in this development evaluation, but they do not establish an REO-specific advantage. Drug-level results and all absolute metrics are reported in Table 2 and Additional files 3–5 (Tables S3–S5).

### Matched-gene contrasts varied with the candidate rule

Within the MAD-selected 200-gene set, the equal-drug pair-minus-expression difference was −0.01765 (family-18 adjusted conditional interval, −0.03117 to −0.00406). Within the response-association-selected 200-gene set, the corresponding difference was −0.00225 (−0.01438 to +0.00949). Thus, in the fixed comparisons, full pair encoding did not outperform continuous expression on average within the same candidate-gene budget. The response-association pair-minus-expression interval included both negative and positive values; it does not establish equality or non-inferiority. The rank representation, all drug-specific estimates, and the other candidate-by-representation contrasts are retained in Additional files 3–5 (Tables S3–S5). The matched-candidate contrasts are summarized in Figure 3A.

The matched comparisons concerned both representation and gene selection: the response-associated rule used continuous expression–response association within training tissues. It was nested within training folds, so the observed pattern describes the specified selection-and-learning procedures.

The rank arm provided a further comparison of ways to encode gene order (Table S10). Pair minus rank was −0.00521 (nominal conditional 95% interval, −0.01222 to +0.00181) under MAD selection and −0.00666 (−0.01456 to +0.00103) under response-associated selection. Rank-minus-expression point estimates changed from negative under MAD to positive under response-associated selection; only the former interval excluded zero. The full rank contrasts remain in Table S10. These intervals are unadjusted descriptive summaries outside the two adjusted inferential families and do not allocate causal shares to amplitude loss, dimensional expansion, collinearity, or penalization.

### Zero molecular solutions retained tissue-reference predictions

Among 840 molecular outer-fold fits, 430 returned the exact zero-molecular-coefficient solution. These fits were valid outputs of the prespecified elastic-net path, not failed computations. Across the 56 molecular full-development models, 29 had zero molecular coefficients. The resulting predictions retained the tissue component; they must not be interpreted as evidence that the model learned a molecular resistance signature. The remaining models had non-zero coefficients, which alone does not establish predictive utility or biological mechanism. No solver warning or fallback was recorded in the repeated development evaluation. The per-drug and per-arm audit is reported in Additional file 6 (Table S6). Figure 3B shows the molecular coefficient status of the full-development fits.

### Tissue-held-out contrasts were frequently not estimable

Among the 80 prespecified drug-by-tissue cells in the ten-tissue display grid, 34 had a defined within-tissue Spearman difference for response-association pair minus strong expression. Only Telomerase Inhibitor IX and BMS-536924 had defined strict means across all ten tissues (−0.0124 and −0.0526, respectively). Because the other six drugs did not have fully defined ten-tissue means, the prespecified strict eight-drug macro was NA. We did not replace it with a mean over the 34 available cells. The defined and undefined cells are shown in Figure 4A. Full results for all 23 held-out tissues, including undefined correlations and model status, are in Additional files 1 and 7 (Tables S1 and S7).

### Cross-screen prediction did not establish a pair-specific gain

The external comparison did not establish a relative REO benefit despite positive observed absolute ranking associations. The full-development response-association pair model had an equal-drug external Spearman correlation difference of −0.000372 relative to the strong-expression model (family-9 adjusted conditional interval, −0.04217 to +0.04350). Drug-specific differences varied in direction, with exact zeros for MIRA-1, PLX-4720, and Tamoxifen; every estimate is retained in Table 2. Every non-identical drug-specific family-adjusted interval included zero; the strict macro interval also included zero. The three exact zero differences occurred because both full-development models had zero molecular coefficients for those compounds, so their predictions were identical under the frozen implementation. These are not equivalence or non-inferiority results.

Within the same 200-gene candidate set, the external MAD pair-minus-expression macro difference was +0.00782 (family-18 adjusted conditional interval, −0.01884 to +0.03450), and the response-association pair-minus-expression macro difference was −0.00671 (−0.02562 to +0.01369). Neither comparison showed a reliable average pair-encoding advantage. External intervals are conditional on the fixed models and sampled external groups; they do not account for model refitting. As GDSC2 and CTRPv2 share the CCLE expression source and differ in screening protocols and dose windows, these results describe retrospective cross-screen within-drug ranking transportability rather than independent transcriptomic-platform transfer, absolute AUC calibration, prospective validation, or patient-level clinical performance. Full estimates and common group counts are in Table 2 and Additional file 2 (Table S2). Figure 4B displays the external drug-wise and macro contrasts.

Observed mean external correlations were positive: the strict equal-eight-drug mean correlation was 0.2284 for tissue, 0.2602 for strong expression, and 0.2599 for response-associated pair (Table S9). The near-zero primary contrast therefore does not show that all transferable predictive association disappeared. These absolute macros have no precomputed intervals, and the point difference from tissue is not a new inferential molecular-increment result. Table S9 includes every drug and arm, correlation, sample counts, and final zero-coefficient status; available drug-level absolute-rho intervals are transcribed from the original paired comparisons collectively involving five molecular arms. No interval was derived by subtracting endpoints or borrowed for identical predictions. Uncalibrated cross-screen RMSE and R-squared were not evaluated because the assay AUC scales and dose windows differ.

External rank minus expression was +0.01658 under MAD and −0.00872 under response-associated selection; pair minus rank was −0.00876 and +0.00200, respectively. These are descriptive strict macro points without precomputed intervals (Table S10). Together with the adjusted pair-minus-expression results, they do not establish general pair superiority, equivalence, or non-inferiority.

Table 2. Primary REO-pair minus strong-expression differences in drug-wise Spearman correlation. Development estimates use repeated nested out-of-fold predictions; external estimates use locked full-development models evaluated once against CTRPv2. Confidence limits are family-9 adjusted conditional 95% bootstrap intervals. The development bootstrap is conditional on recorded out-of-fold predictions; the external bootstrap is conditional on fixed trained models and the observed external groups. Neither interval includes uncertainty from repeating the complete model-development procedure. A zero molecular solution retains the tissue component without a molecular coefficient. These intervals do not establish equivalence or non-inferiority.

| Compound | Development groups | Development Δρ (adjusted 95% conditional CI) | Shared external groups | External Δρ (adjusted 95% conditional CI) |
| --- | --- | --- | --- | --- |
| MIRA-1 | 625 | 0.02499 (-0.01522 to 0.06353) | 108 | 0.00000 (0.00000 to 0.00000) |
| PLX-4720 | 622 | -0.01053 (-0.03341 to 0.01177) | 167 | 0.00000 (0.00000 to 0.00000) |
| Bortezomib | 620 | -0.01706 (-0.04072 to 0.00553) | 169 | 0.00060 (-0.09204 to 0.10248) |
| Telomerase Inhibitor IX | 620 | -0.00389 (-0.04047 to 0.03030) | 168 | -0.02369 (-0.10873 to 0.05833) |
| Tamoxifen | 620 | 0.00132 (-0.02128 to 0.02269) | 165 | 0.00000 (0.00000 to 0.00000) |
| BMS-536924 | 620 | -0.03302 (-0.07615 to 0.00946) | 165 | 0.03033 (-0.09321 to 0.16706) |
| Gefitinib | 620 | -0.00281 (-0.04802 to 0.04344) | 161 | -0.03721 (-0.12540 to 0.05041) |
| Cediranib | 620 | -0.00605 (-0.04618 to 0.03510) | 82 | 0.02700 (-0.23477 to 0.28333) |
| Equal-drug macro | 8-drug mean | -0.00588 (-0.01871 to 0.00661) | 8-drug mean | -0.00037 (-0.04217 to 0.04350) |

*Note:* Δρ is the response-association full-pair model minus the strong 2,000-gene expression model. Intervals are family-9 adjusted conditional 95% bootstrap intervals. The development interval conditions on saved out-of-fold predictions; the external interval conditions on the locked models and observed external groups. Full-development coefficient counts and zero-solution status are detailed in Additional file 6 (Table S6).

### Constructed controls illustrated distinct prediction behaviors

The constructed controls illustrated how prediction could arise under different specified signals, without assigning a cause to the real-drug results. All 40 generated cohorts and 320 cohort-by-arm summaries were retained. Across the 840 selected molecular outer refits, 615 were zero solutions: 210/210 in tissue-only, 198/210 in weak molecular, 120/210 in common amplitude, and 87/210 in predefined pair order. This denominator excludes inner-grid fits, including deliberately evaluated analytic null candidates. Undefined residual correlations remain blank in Table S11 rather than being replaced by zero.

In the tissue-only control, predictions reconstructed tissue outcomes to floating-point precision despite every molecular outer model being zero. Under the weak molecular signal, performance remained close to the tissue reference. The common-amplitude control exposed a clear difference: every rank/pair outer fit was zero, with mean correlation 0.1541 versus 0.9830 for each continuous arm. This scenario generated all 267 recorded warning occurrences (89 per continuous arm), with no recorded fallback; warning type is not recoverable from the aggregate counts.

The deliberately constructed pair-order scenario favored the response-associated pair model, whose mean correlation was 0.9055 versus 0.7573 for expression and 0.7605 for ranks. Both driver genes were retained in all response-associated and no MAD outer candidate sets; the MAD-selected arms remained close to the tissue reference. These controls illustrate possible pipeline behavior under specified signals without identifying the mechanisms of the real-drug findings. Every cohort-by-arm value, generation invariant, warning count, resource observation and provenance limit remains in Tables S11–S12 and Additional file 14.
