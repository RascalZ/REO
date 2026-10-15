# BMC Genomics Research Article — working draft

**Draft status:** BMC-structured working copy derived from the evidence-bounded manuscript. Author, institutional, ethics, funding, competing-interest, contribution, data-license, code-archive, and generative-AI disclosure details remain placeholders; this file is not ready for submission. The original evidence draft and all historical submissions remain unchanged.

## Title page

**Title:** Assessing the incremental value of relative expression ordering for continuous drug-response prediction: a repeated nested cross-validation and cross-screen study across eight anticancer compounds

**Authors:** [AUTHOR INPUT REQUIRED: confirm complete author list and order]

**Affiliations:** [AUTHOR INPUT REQUIRED: provide full institutional addresses and map each author]

**Corresponding author:** [AUTHOR INPUT REQUIRED: name, postal address, and email]

## Abstract

**Background:** Within-sample relative expression ordering (REO) has been proposed as a robust representation for cancer prediction, but its incremental value over expression-based models for continuous drug-response prediction remains uncertain. We compared a tissue-only reference, a strong 2,000-gene expression model, and matched 200-gene expression, within-sample rank, and gene-pair representations under two training-fold candidate-selection rules. Evaluation used five-fold outer cross-validation repeated three times with three-fold inner tuning across eight fixed compounds and 626 development cell-line groups. A locked cross-screen evaluation assessed full-development models in 175 non-overlapping cell-line groups with eligible CTRPv2 responses. The external analysis was executed once after model and statistical rules were frozen.

**Results:** In repeated nested development evaluation, the predefined REO-pair minus strong-expression difference in equal-drug mean Spearman correlation was −0.00588 (family-9 Bonferroni-adjusted conditional bootstrap interval, −0.01871 to 0.00661). The corresponding macro improvement in RMSE, defined as strong-expression RMSE minus REO-pair RMSE, was −0.00075 (unadjusted conditional 95% interval, −0.00124 to −0.00028), favoring the expression model. In the locked CTRPv2 evaluation, the corresponding equal-drug difference in Spearman correlation was −0.000372 (family-9 adjusted conditional interval, −0.04217 to 0.04350; 1,185 drug–model responses). Twenty-nine of 56 final molecular model arms had zero molecular coefficients. The strict leave-one-tissue-out macro comparison was undefined because only 34 of 80 prespecified drug–tissue contrasts had defined correlations.

**Conclusions:** The evaluated REO pipeline did not show a reliable average ranking improvement over the strong expression baseline in the prespecified eight-compound setting. These findings do not establish equivalence or non-inferiority and do not rule out other REO formulations. They show why tissue contribution, zero molecular solutions, candidate-selection policy, and assay context need to be reported when evaluating molecular drug-response models.

## Keywords

Cancer cell lines; Drug response prediction; Pharmacogenomics; Relative expression ordering; Nested cross-validation; Elastic net

## Background

Large cancer-cell-line pharmacogenomic resources have enabled systematic comparisons between molecular profiles and drug-response phenotypes [1–3,22,23]. These screens provide useful experimental systems for studying response-associated patterns, but their dose ranges, assay procedures, and response summaries are not interchangeable. Cross-study analyses have documented that performance can change when models move between screens, while recent benchmarking studies emphasize evaluation design, response metrics, and comparisons with simple baselines [3–6,18,19]. Accordingly, a model that predicts a cell-line screen outcome should be interpreted as a predictor of that experimental endpoint unless tested against an independent biological or clinical outcome [12].

Within-sample relative expression ordering (REO) represents the sign of the expression difference between two genes in the same sample. Because a common strictly monotone transformation preserves order, pairwise order can be less sensitive than absolute values to some sample-wide scale changes. REO-based methods and qualitative signatures have been applied to differential-expression inference and treatment-response or resistance classification [13–17,20,21]. Those studies motivate evaluation of order-based representations, but their settings include differential-expression tasks and patient-level categorical signatures; they do not establish that pair features improve continuous drug-response prediction in cell lines. In addition, encoding all pairs among 200 genes produces up to 19,900 features and discards expression magnitude. The resulting model may therefore behave differently from either a 200-dimensional rank model or a continuous-expression model, even though a rank vector and its complete pair expansion encode the same within-sample ordering.

High-dimensional biological prediction can be biased when outcome-informed feature selection, preprocessing, or tuning uses information from held-out observations. Guidance and empirical studies recommend nesting all learned steps within training partitions, respecting related sample identities, reporting uncertainty appropriate to the independent unit, and making the information flow auditable [7–11]. Nested validation alone does not establish useful generalization: finite-sample variability, tissue composition, endpoint definition, model complexity, and the choice of baseline also matter [5,6,9,18,19]. These concerns are especially relevant when feature selection itself is response-associated.

Here, we evaluated whether a response-associated full gene-pair model added drug-ranking information beyond a strong 2,000-gene expression baseline in a fixed panel of eight compounds. We also compared expression, within-sample rank, and complete pair encodings within the same 200-gene candidate sets selected by either a training-only variability rule or a training-only response-association rule. A tissue-only reference was used to describe the contribution of tissue means; repeated nested cross-validation, leave-one-tissue-out analysis, zero-solution auditing, and one locked GDSC2-to-CTRPv2 cross-screen evaluation were used to characterize the scope and limitations of the comparisons. The study was designed as a comparative evaluation of a specified learning pipeline, not as a search for the best-performing drug subset or as a clinical prediction study.

## Methods

Machine-readable supporting tables are supplied as Additional files 1–8, corresponding to Tables S1–S8, respectively.

### Study design and estimands

This was a retrospective computational comparison of continuous drug-response prediction models in cancer cell lines. The primary question was whether a training-fold response-associated gene-pair model improved drug-wise response ranking over a strong expression baseline across a fixed panel of eight compounds. Secondary questions compared expression, within-sample rank, and full pair encodings within matched 200-gene candidate sets; quantified performance relative to a tissue-only reference; described zero molecular solutions and tissue-held-out behavior; and evaluated transfer of fixed models from GDSC2 to CTRPv2.

The eight-compound panel was carried forward from the project's identity-, coverage-, and dose-metadata-frozen panel. The panel was not altered according to results of the current repeated evaluation or CTRPv2 evaluation. However, the formulation of the candidate-selection and representation hypotheses was informed by earlier project analyses. Accordingly, the current repeated evaluation is a prespecified, outcome-informed development comparison, not an independent confirmatory study. The CTRPv2 evaluation was locked before its response values were opened by the new evaluation workflow and executed once; CTRPv2 had nevertheless appeared in an earlier manuscript, so the evaluation is not prospectively blind in the history of the project.


The cohorts, nested partitions, model arms, and evaluation boundaries are summarized in Figure 1.
### Data sources, expression measurements, and compound mapping

Drug-response measurements for model development were taken from the frozen GDSC2 fitted dose-response release dated 27 October 2023 [23]. Gene expression was taken from the historical Cancer Cell Line Encyclopedia (CCLE) RSEM gene TPM matrix dated 29 September 2018 [1]. The expression matrix contained 1,019 models and 57,820 measured gene records before the project-specific identity and annotation mapping. The analysis used 19,836 uniquely mapped protein-coding stable Ensembl gene identifiers. Gene-expression values were used as TPM; no raw-count differential-expression analysis or DESeq2 step was part of the reconstructed pipeline. DepMap 24Q4 Model/README and Cell Model Passports metadata snapshots were used for identity and tissue mapping only [25,26]; their expression data were not mixed into the 2018 CCLE matrix. CTRPv2 response records used only in the locked external evaluation derive from the original CTRP screen [22].

Eight compounds were fixed by their GDSC2 identifiers: MIRA-1 (1931), PLX-4720 (1036), Bortezomib (1191), Telomerase Inhibitor IX (1930), Tamoxifen (1199), BMS-536924 (1091), Gefitinib (1010), and Cediranib (1922). Their CTRPv2 compound identifiers and chemical-identity evidence are reported in Additional file 8 (Table S8). Name and structure cross-references, including Standard InChIKey matches where available, supported the compound mapping; this metadata match does not chemically authenticate the material used in either assay. GDSC2 and CTRPv2 concentration ranges differed and were not harmonized post hoc.

### Response outcome and analysis units

The continuous outcome was the source-reported fitted dose-response area under the curve (AUC), retained on the source scale without dichotomization, thresholding, or outcome transformation. Under the GDSC documentation's viability-curve convention, a larger AUC corresponds to greater residual viability across tested concentrations and thus lower in-vitro drug sensitivity [24]. The analysis unit was a canonical cell-line biological group for one compound. Related or duplicate model identities were grouped before data splitting so that a biological group could not occur in both training and test partitions. Each eligible group contributed at most one response record per compound after the frozen cohort curation rules.

### Cohort curation and eligibility

Cell-line identifiers were reconciled across expression, response, and metadata sources using the versioned cohort-v2 mapping, including canonical model identifiers, biological-group identifiers, and tissue labels. Drug identity was mapped using the frozen compound table and chemical evidence. Eligibility was determined from identity, expression coverage, response validity, and the frozen cohort rules, not from current model performance. The resulting development cohort contained 626 unique biological groups and 4,967 compound–group response records across the eight drugs. Per-drug group counts and the distinction between unique groups and compound–group records are shown in Table 1. The locked CTRPv2 cohort contained 175 unique external groups with 1,185 valid compound–group response records after repeated curve records for a compound and canonical model were aggregated by the prespecified median rule.

### Candidate-gene selection and molecular representations

Every learned operation, including detection filtering, candidate selection, representation construction, scaling, tissue adjustment, regularization-path calculation, and hyperparameter selection, was repeated using only the training rows of the relevant inner fit. The outer test fold was used once for evaluation and did not contribute to feature selection or tuning.

Within each training fit, a gene was considered detected if its raw TPM was at least 1 in at least 10% of training samples. The median-absolute-deviation (MAD) rule selected genes with non-zero MAD in log2(TPM+1), ordered by decreasing MAD and then stable Ensembl identifier for ties. The response-association rule first centered log2(TPM+1) and the continuous AUC outcome within each training tissue. Among detected genes with non-zero within-tissue expression residual variance, genes were ranked by the absolute Pearson correlation between their tissue-centered expression and tissue-centered response. Stable Ensembl identifiers broke exact ties. No second MAD ranking was applied. If the training response was exactly constant within every tissue, all eligible genes received a score of zero and were ordered by stable identifier. Both candidate rules selected up to 200 genes for the matched-budget comparisons. The strong-expression baseline selected up to 2,000 genes by the same training-only detection and MAD policy.

For each selected 200-gene set, three representations were evaluated. The expression representation used log2(TPM+1). The rank representation assigned each sample average midranks among the selected genes, using the raw TPM ordering and exact ties. The full pair representation included every unordered pair of selected genes, with orientation determined by ascending stable Ensembl identifier; each feature was +1 when the first gene's raw TPM exceeded the second, −1 when it was lower, and 0 for an exact tie. At 200 genes this yields at most 19,900 pair features. Constant and high-tie pairs were retained. The rank vector and complete pair vector encode the same within-set order information, but their linear predictors and elastic-net penalties differ because the pair encoding expands that information into a larger, correlated feature space.

All molecular inputs had to be complete, finite, non-negative float64 TPM values. No missing-value imputation or post hoc clipping was used. A structural expression-schema mismatch was treated as a hard failure; a sample with invalid molecular TPM was not silently imputed. No such molecular abstentions occurred in the locked external evaluation.

### Prediction models and hyperparameter selection

Eight model arms were evaluated: (1) a tissue-only reference, which predicted the mean training AUC for the observed tissue and the overall training mean for an unseen tissue; (2) a strong expression elastic-net model using up to 2,000 MAD-selected genes; and (3–8) the two candidate-selection rules crossed with the three matched 200-gene representations.

For each molecular arm, features were constructed before standardization and tissue adjustment. A training-only standard scaler was fit to the molecular features, and the continuous outcome was standardized using the current training data. Tissue was represented by one-hot fixed effects that were not penalized. The molecular design and outcome were projected on the orthogonal complement of the training tissue design using a Frisch–Waugh–Lovell residualization, and an elastic net without an intercept was fit to the residualized data. For a tissue absent from training, the fitted tissue contribution was zero while the molecular component was retained.

The elastic-net grid contained five alpha ratios (1, 0.3, 0.1, 0.03, and 0.01) crossed with L1 ratios 0.5 and 1. For each individual training fit, alpha_max was recomputed from that fit's tissue-residualized features and outcome as max|Zrᵀyr|/(n × L1 ratio); the candidate alpha was the selected ratio multiplied by alpha_max. Models used cyclic coordinate descent with a maximum of 4,000 iterations and tolerance 0.001. Hyperparameters were selected by the lowest pooled inner out-of-fold RMSE on the original AUC scale; exact ties were resolved by the first candidate in the frozen grid order. The ratio-1 endpoint was the analytic zero-molecular-coefficient solution and was retained as a legitimate candidate. Zero molecular coefficients were not considered a fit failure: the tissue component remained in the prediction. Predictions were not clipped to the [0,1] interval.

### Repeated nested cross-validation

The development evaluation used five outer folds repeated three times, with three inner folds for tuning. Outer splits were tissue-stratified group splits based on the 626 frozen biological groups, using seeds 20261013, 20261014, and 20261015. Drug-specific records inherited the global group assignments and were restricted to eligible groups for that drug. The inner folds were generated within each outer-training partition, with the frozen repeat/fold-specific seed rule. All eight model arms for a given drug shared the same outer and inner partitions. The 200-gene selection, all feature construction and preprocessing, alpha_max calculation, and model selection were rerun for each inner-training fit. After selecting a hyperparameter set, each molecular model was refit on the entire outer-training partition and evaluated once on the held-out outer fold. The tissue-only model required one training fit per outer context and no hyperparameter search.

This design generated 960 outer tasks (eight compounds × eight model arms × 15 outer folds), 26,160 successful fit events, and 119,208 outer-fold predictions. Repeats reused the same biological groups under different partitions and were not treated as additional independent samples.

### Full-development models and locked CTRPv2 evaluation

After the development comparison was complete, 64 full-development models (eight compounds × eight arms) were trained from the full eligible development data using a frozen grouped three-fold inner tuning procedure. The prespecified grid, training-only feature pipeline, and RMSE objective were unchanged. The expected total was 1,744 fit events, including the inner candidate fits and the final all-development refits. Model artifacts were reloaded in an independent process and checked against locked source and prediction audits. Response-free CTRPv2 predictions were generated and locked before the new evaluation workflow accessed external response values.

The single locked external evaluation matched predictions to CTRPv2 responses by compound and canonical cell-line group. For each compound, repeated CTRPv2 curve records for the same canonical group were collapsed by their median, yielding 1,185 valid drug–group responses from 1,201 source curve records. External comparisons were restricted to common groups for the two models in each paired contrast. The primary estimand was the within-drug Spearman correlation difference, response-association pair minus strong expression; the eight-drug macro was the strict equal-weight mean over the fixed eight compounds. A strict macro was not recalculated after omitting undefined drugs. All 9,480 planned compound-by-model predictions were finite and covered, with no abstentions. This analysis evaluates retrospective within-drug ranking transport from GDSC2 response labels to CTRPv2 response labels using the same CCLE expression source; it is neither cross-expression-platform validation nor clinical validation.

### Leave-one-tissue-out sensitivity analysis

Each of the 23 metadata-defined tissues was held out in turn. For each held-out tissue, all groups from that tissue were excluded from model development; the same drug panel, model arms, training-only feature pipeline, three-fold inner tuning, and RMSE selection rule were used on the remaining tissues. The held-out tissue was treated as unseen by the model. The primary display grid comprised the ten tissues with at least 20 development groups for every compound, yielding 80 prespecified drug-by-tissue contrasts. A drug-specific strict mean required a defined contrast in all ten tissues, and the strict eight-drug macro required a defined strict mean for all eight drugs. Undefined correlations were preserved as NA; no available-cell mean was substituted for the strict macro. This was a sensitivity analysis without a new inferential family or confidence intervals.

### Performance metrics and statistical uncertainty

For each drug, model, and repeat, Spearman correlation, RMSE, mean absolute error, and R² were calculated from the complete outer out-of-fold predictions on the original response scale. The three repeat-specific metric values were then averaged. Predictions were not averaged across repeats before ranking, and the repeated predictions were not concatenated to inflate the sample size. Drug-macro values were equal-weight averages over all eight drugs; if any prespecified drug or repeat component was undefined, the strict macro was reported as NA with the reason retained.

The primary inferential family comprised nine Spearman contrasts: the primary REO-pair minus strong-expression contrast for each of eight drugs and the strict equal-drug macro. A separate 18-estimand family comprised pair-minus-expression contrasts within the MAD and response-association candidate rules, again for each drug and the equal-drug macro. Development intervals were percentile bootstrap intervals based on 10,000 shared resamples of the 626 biological groups, with the same group multiplicities applied across drugs, arms, and repeats. The family-adjusted intervals used Bonferroni-adjusted two-sided percentile limits at nominal family alpha 0.05. These intervals are conditional on the recorded out-of-fold predictions; they do not include uncertainty from repeating feature selection, tuning, or model fitting and are not guaranteed full-algorithm cross-validation coverage. The inferential families are adjusted separately and do not provide simultaneous coverage over every analysis in the study. Other model contrasts and performance metrics are descriptive with unadjusted conditional intervals.

For the external evaluation, 20,000 shared bootstrap resamples were drawn from the 175 external biological groups. The same group weights were used across all drugs, arms, and paired contrasts. The primary family remained nine estimands and the matched-candidate representation family remained 18. Each reported estimand had 20,000 valid draws. External intervals are conditional on the fitted full-development models and the observed external groups; they do not incorporate model-refitting uncertainty. Zero intervals caused by two identical zero-molecular predictions are algebraic identities and do not constitute evidence of equivalence or non-inferiority. No non-inferiority margin was specified, and no equivalence or non-inferiority test was performed. No new full-pipeline permutation test was run in this eight-drug analysis.

### Reproducibility and data governance

Analyses were governed by versioned protocols, configuration files, source hashes, cohort locks, execution locks, and result manifests stored in the project repository. The manuscript-table/figure generation environment used Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, SciPy 1.16.2, scikit-learn 1.7.2, Matplotlib 3.10.7, and psutil 7.1.0; the modeling and locked-external environments are separately bound in their execution records and must be reported from those records rather than inferred from the current plotting environment. Sample-level CTRPv2 response values, external group bootstrap weights, fitted model checkpoints, and individual predictions are not included in the public repository. The repository records source versions, accession/download provenance, code, aggregate tables, figure/table generation scripts, output hashes, and software environment metadata. Upstream data terms remain applicable; the repository license does not relicense source datasets. Dataset citations, repository accessions, and the final public release/archive DOI must be verified and inserted before submission.

[AUTHOR INPUT REQUIRED: Confirm the generative-AI disclosure required by the journal. OpenAI Codex assisted with restructuring and drafting/editing this working manuscript and with reference and file-format checks. Authors must confirm the accurate scope, whether code-development assistance also requires disclosure, the required tool/model/version details, and final wording. Do not submit this placeholder.]

## Results

### Cohort composition and evaluation coverage

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

The repeated nested development evaluation completed all 960 outer tasks, 26,160 fit events, and 119,208 outer predictions. No solver warnings or fallback fits were recorded. The CTRPv2 evaluation generated finite predictions for all 9,480 planned compound-by-model rows, with zero abstentions. The external comparison therefore had complete prediction coverage under the frozen eligibility rules.

### REO-pair versus strong-expression performance in development

The primary equal-drug macro Spearman correlation was 0.2931 for the response-association pair model and 0.2990 for the strong 2,000-gene expression model. The prespecified difference, pair minus strong expression, was −0.00588 (family-9 Bonferroni-adjusted conditional interval, −0.01871 to 0.00661). Six of eight drug-specific point estimates were negative and two were positive; every drug-specific family-adjusted interval included zero. The macro difference was negative in all three repeated partitions (−0.00393, −0.00655, and −0.00716), which describes partition-level direction consistency but does not create three independent replications. The paired drug-wise and equal-drug estimates are shown in Figure 2.

For the macro RMSE, the strong-expression model had 0.06217 and the REO-pair model had 0.06292. The paired improvement, defined as strong-expression RMSE minus pair RMSE, was −0.00075 (unadjusted conditional 95% interval, −0.00124 to −0.00028), favoring the strong-expression model on this error measure. The tissue-only macro Spearman was 0.2366; the descriptive response-association pair minus tissue-only difference was +0.05654 (unadjusted conditional interval, +0.04299 to +0.07004). The strong-expression minus tissue-only difference was +0.06242 (+0.04926 to +0.07577). These tissue-reference contrasts indicate that both molecular modeling pipelines can add ranking information over the tissue reference in this development evaluation, but they do not establish an REO-specific advantage. Drug-level results and all absolute metrics are reported in Table 2 and Additional files 3–5 (Tables S3–S5).

### Matched-candidate representation comparisons

Within the MAD-selected 200-gene set, the equal-drug pair-minus-expression difference was −0.01765 (family-18 adjusted conditional interval, −0.03117 to −0.00406). Within the response-association-selected 200-gene set, the corresponding difference was −0.00225 (−0.01438 to +0.00949). Thus, in the fixed comparisons, full pair encoding did not outperform continuous expression on average within the same candidate-gene budget. The response-association pair-minus-expression interval included both negative and positive values; it does not establish equality or non-inferiority. The rank representation, all drug-specific estimates, and the other candidate-by-representation contrasts are retained in Additional files 3–5 (Tables S3–S5). The matched-candidate contrasts are summarized in Figure 3A.

The response-association rule produced a descriptive macro difference of +0.05654 over the tissue-only reference, while the strong-expression model differed from the same reference by +0.06242. These results are consistent with molecular predictions adding ranking signal beyond tissue means in some parts of the evaluation, but the average increment was not specific to pair encoding. Candidate selection based on response-associated continuous expression was performed within training folds, yet the contrast remains conditional on the selected outcome-informed workflow.

### Zero molecular solutions and fitted-model status

Among 840 molecular outer-fold fits, 430 returned the exact zero-molecular-coefficient solution. These fits were valid outputs of the prespecified elastic-net path, not failed computations. Across the 56 molecular full-development models, 29 had zero molecular coefficients. The resulting predictions retained the tissue component; they must not be interpreted as evidence that the model learned a molecular resistance signature. The remaining models had non-zero coefficients, which alone does not establish predictive utility or biological mechanism. No solver warning or fallback was recorded in the repeated development evaluation. The per-drug and per-arm audit is reported in Additional file 6 (Table S6). Figure 3B shows the molecular coefficient status of the full-development fits.

### Leave-one-tissue-out sensitivity

Among the 80 prespecified drug-by-tissue cells in the ten-tissue display grid, 34 had a defined within-tissue Spearman difference for response-association pair minus strong expression. Only Telomerase Inhibitor IX and BMS-536924 had defined strict means across all ten tissues (−0.0124 and −0.0526, respectively). Because the other six drugs did not have fully defined ten-tissue means, the prespecified strict eight-drug macro was NA. We did not replace it with a mean over the 34 available cells. The defined and undefined cells are shown in Figure 4A. Full results for all 23 held-out tissues, including undefined correlations and model status, are in Additional files 1 and 7 (Tables S1 and S7).

### Locked CTRPv2 cross-screen evaluation

The full-development response-association pair model had an equal-drug external Spearman correlation difference of −0.000372 relative to the strong-expression model (family-9 adjusted conditional interval, −0.04217 to +0.04350). Drug-specific point differences were heterogeneous: −0.0237 for Telomerase Inhibitor IX, +0.0303 for BMS-536924, −0.0372 for Gefitinib, +0.0270 for Cediranib, +0.0006 for Bortezomib, and exactly zero for MIRA-1, PLX-4720, and Tamoxifen. Every non-identical drug-specific family-adjusted interval included zero; the strict macro interval also included zero. The three exact zero differences occurred because both full-development models had zero molecular coefficients for those compounds, so their predictions were identical under the frozen implementation. These are not equivalence or non-inferiority results.

Within the same 200-gene candidate set, the external MAD pair-minus-expression macro difference was +0.00782 (family-18 adjusted conditional interval, −0.01884 to +0.03450), and the response-association pair-minus-expression macro difference was −0.00671 (−0.02562 to +0.01369). Neither comparison showed a reliable average pair-encoding advantage. External intervals are conditional on the fixed models and sampled external groups; they do not account for model refitting. As GDSC2 and CTRPv2 share the CCLE expression source and differ in screening protocols and dose windows, these results describe retrospective cross-screen within-drug ranking transportability rather than independent transcriptomic-platform transfer, absolute AUC calibration, prospective validation, or patient-level clinical performance. Full estimates and common group counts are in Table 2 and Additional file 2 (Table S2). Figure 4B displays the external drug-wise and macro contrasts.

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

## Discussion

In the prespecified eight-compound comparison, the response-associated REO-pair pipeline did not show a reliable average improvement in drug-wise ranking over the strong expression baseline in either repeated nested development evaluation or the locked CTRPv2 evaluation. The internal Spearman difference was close to zero and the internal RMSE contrast slightly favored expression; the external macro ranking difference was also close to zero with an interval that allowed both positive and negative effects. This is a result about the evaluated candidate rules, elastic-net path, response endpoint, compounds, and cell-line cohorts. It neither establishes equivalence or non-inferiority nor proves that other REO models cannot be useful. The result is consistent with recent benchmarking work showing that complexity or feature transformation does not guarantee improvement over well-defined simpler predictors [5,6,18,19].

Earlier REO studies provide a rationale for evaluating within-sample order and have reported useful qualitative signatures in particular data settings [13–17,20,21]. Their reported tasks and estimands differ from the present continuous in-vitro response-ranking task, so their performance values are not direct comparators. In the present matched-candidate analyses, the full pair encoding did not outperform continuous expression on average. The response-association candidate rule was based on continuous expression–response association within training tissues; it was nested inside each training fit, but it does not constitute a selector neutral to expression magnitude. Thus, the results do not isolate a universal “order-only” effect. Moreover, rank vectors and complete pair vectors carry the same ordering information but define different linear-elastic-net hypothesis classes and penalty geometry. The observed differences cannot be attributed uniquely to lost amplitude, expanded dimensionality, collinearity, or regularization.

Tissue reference performance and molecular zero solutions further constrain interpretation. Both the response-associated pair pipeline and the strong-expression model improved on the tissue-only reference in the development evaluation, but the difference between those molecular models did not favor REO on average. Nearly half of the outer-fold molecular fits and more than half of the final molecular model arms selected the exact zero-coefficient endpoint. In those fits, prediction retains the tissue component and the associated correlation does not demonstrate that a molecular pattern was learned. The three exact-zero external primary differences are an especially clear example: both compared models had zero molecular coefficients, so identical predictions yielded a zero contrast by construction.

The cross-screen evaluation is best interpreted as retrospective within-drug ranking transport from GDSC2 responses to CTRPv2 responses under a shared CCLE expression source. Cross-study analyses and recent evaluation frameworks show that assay protocols, dose ranges, response summaries, and aggregation choices can affect agreement and model transport [3,4,6,22,23]. The current evaluation does not test transfer to a different expression-measurement platform, absolute calibration between source AUC values, prospective prediction, or patient outcome. CTRPv2 had been used in an earlier project manuscript, so the present locked execution should not be described as wholly blind in the project's history. The zero-prediction-abstention coverage is complete for the locked eligible set, but complete coverage alone does not imply strong predictive value.

Several limitations affect the strength and scope of inference. First, although supervised feature selection, representation construction, scaling, tissue adjustment, and tuning were repeated inside training partitions, the eight-compound panel and comparison hypotheses were informed by earlier project analyses. The present run is therefore a prespecified development comparison rather than an independent confirmatory study. Second, the bootstrap intervals condition on saved outer-fold predictions or fixed final models and do not include uncertainty from repeating the entire feature-selection, tuning, and fitting procedure. The repeated partitions reuse the same biological groups. Third, this eight-drug analysis did not include a new full-pipeline label-permutation experiment; earlier synthetic mechanisms test software behavior under defined scenarios but do not provide an empirical null distribution for these drug outcomes. Fourth, the leave-one-tissue-out strict macro was undefined because many prespecified within-tissue correlations were not estimable; the defined-cell average cannot replace the planned macro. Finally, only eight compounds and one external screen were studied, and the AUC dose windows differed across screens. Larger and independently designed screen comparisons would be needed to test broader generality.

## Conclusions

In this fixed eight-compound cell-line study, the evaluated REO-pair pipeline did not show a reliable average improvement in drug-response ranking over the strong expression baseline in repeated nested development evaluation or the locked retrospective CTRPv2 cross-screen comparison. The matched representation analyses do not establish a pair-specific advantage, and frequent zero molecular solutions show that some fitted predictors reduce to their tissue component. These findings delimit the performance of the specified candidate-selection rules, elastic-net models, response endpoint, compounds, and cohorts; they do not establish equivalence or non-inferiority and should not be generalized to all REO methods or clinical prediction.

## Declarations

### Ethics approval and consent to participate

[AUTHOR INPUT REQUIRED: verify the ethics statement applicable to secondary analysis of public cancer-cell-line datasets; provide the approving/waiving body and identifier if relevant. Do not submit this placeholder.]

### Consent for publication

[AUTHOR INPUT REQUIRED: confirm the appropriate statement for the datasets and manuscript.]

### Availability of data and materials

[AUTHOR INPUT REQUIRED: specify persistent links and versions for the public GDSC2, CCLE, CTRPv2, DepMap, and Cell Model Passports resources used; add the archived analysis-code release DOI and exact repository tag after release; state any upstream license/access restrictions. The repository URL alone is not a versioned archive DOI.]

### Competing interests

[AUTHOR INPUT REQUIRED: declare all financial and non-financial interests, or confirm the journal's no-competing-interests wording.]

### Funding

[AUTHOR INPUT REQUIRED: list funding sources and funder roles, or confirm whether the study received no specific funding.]

### Authors' contributions

[AUTHOR INPUT REQUIRED: provide author initials and contribution roles; obtain approval from all authors.]

### Acknowledgements

[AUTHOR INPUT REQUIRED: identify contributors and obtain permission, or confirm that none apply.]

## References

1. Ghandi M, et al. Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature*. 2019;569:503–508. https://doi.org/10.1038/s41586-019-1186-3.
2. Feizi N, et al. PharmacoDB 2.0: improving scalability and transparency of in vitro pharmacogenomics analysis. *Nucleic Acids Res*. 2022;50(D1):D1348–D1357. https://doi.org/10.1093/nar/gkab1084.
3. Xia F, et al. A cross-study analysis of drug response prediction in cancer cell lines. *Brief Bioinform*. 2022;23(1):bbab356. https://doi.org/10.1093/bib/bbab356.
4. Bayer FP, Gander M, Kuster B, The M. CurveCurator: a recalibrated F-statistic to assess, classify, and explore significance of dose-response curves. *Nat Commun*. 2023;14:7902. https://doi.org/10.1038/s41467-023-43696-z.
5. Bernett J, et al. Critical evaluation of drug response prediction models with DrEval. *Nat Commun*. 2026;17:4238. https://doi.org/10.1038/s41467-026-72903-w.
6. Codicè F, et al. The specification game: rethinking the evaluation of drug response prediction for precision oncology. *J Cheminform*. 2025;17:33. https://doi.org/10.1186/s13321-025-00972-y.
7. Kapoor S, Narayanan A. Leakage and the reproducibility crisis in machine-learning-based science. *Patterns (N Y)*. 2023;4(9):100804. https://doi.org/10.1016/j.patter.2023.100804.
8. Bernett J, et al. Guiding questions to avoid data leakage in biological machine learning applications. *Nat Methods*. 2024;21:1444–1453. https://doi.org/10.1038/s41592-024-02362-y.
9. Vabalas A, Gowen E, Poliakoff E, Casson AJ. Machine learning algorithm validation with a limited sample size. *PLoS One*. 2019;14(11):e0224365. https://doi.org/10.1371/journal.pone.0224365.
10. Lewis MJ, et al. nestedcv: an R package for fast implementation of nested cross-validation with embedded feature selection designed for transcriptomics and high-dimensional data. *Bioinform Adv*. 2023;3(1):vbad048. https://doi.org/10.1093/bioadv/vbad048.
11. Heil BJ, et al. Reproducibility standards for machine learning in the life sciences. *Nat Methods*. 2021;18:1132–1135. https://doi.org/10.1038/s41592-021-01256-7.
12. Maeser D, Gruener RF, Huang RS. oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data. *Brief Bioinform*. 2021;22(6):bbab260. https://doi.org/10.1093/bib/bbab260.
13. He J, et al. Qualitative transcriptional signature for predicting pathological response of colorectal cancer to FOLFOX therapy. *Cancer Sci*. 2020;111(1):253–265. https://doi.org/10.1111/cas.14263.
14. Tercan B. Robust predictors for drug response of patients with acute myeloid leukemia. *PLoS One*. 2026;21(2):e0343422. https://doi.org/10.1371/journal.pone.0343422.
15. Cai H, et al. Identifying differentially expressed genes from cross-site integrated data based on relative expression orderings. *Int J Biol Sci*. 2018;14(8):892–900. https://doi.org/10.7150/ijbs.24548.
16. Cai H, et al. Identification and characterization of genes with absolute mRNA abundances changes in tumor cells with varied transcriptome sizes. *BMC Genomics*. 2019;20:134. https://doi.org/10.1186/s12864-019-5502-y.
17. Liu K, et al. Systematic exploration of the underlying mechanism of gemcitabine resistance in pancreatic adenocarcinoma. *Mol Oncol*. 2022;16(16):3034–3051. https://doi.org/10.1002/1878-0261.13279.
18. Eckhart L, Lenhof K, Rolli L-M, Lenhof H-P. A comprehensive benchmarking of machine learning algorithms and dimensionality reduction methods for drug sensitivity prediction. *Brief Bioinform*. 2024;25(4):bbae242. https://doi.org/10.1093/bib/bbae242.
19. Shen B, et al. A systematic assessment of deep learning methods for drug response prediction: from in vitro to clinical applications. *Brief Bioinform*. 2023;24(1):bbac605. https://doi.org/10.1093/bib/bbac605.
20. Zheng X, et al. Less is more: relative rank is more informative than absolute abundance for compositional NGS data. *Brief Funct Genomics*. 2025;24:elae045. https://doi.org/10.1093/bfgp/elae045.
21. Chen Y, et al. A qualitative transcriptional signature for predicting extreme resistance of ER-negative breast cancer to paclitaxel, doxorubicin, and cyclophosphamide neoadjuvant chemotherapy. *Front Mol Biosci*. 2020;7:34. https://doi.org/10.3389/fmolb.2020.00034.
22. Seashore-Ludlow B, et al. Harnessing connectivity in a large-scale small-molecule sensitivity dataset. *Cancer Discov*. 2015;5(11):1210–1223. https://doi.org/10.1158/2159-8290.CD-15-0235.
23. Iorio F, et al. A landscape of pharmacogenomic interactions in cancer. *Cell*. 2016;166(3):740–754. https://doi.org/10.1016/j.cell.2016.06.017.
24. Genomics of Drug Sensitivity in Cancer. Drug sensitivity datasets: fitted dose–response AUC documentation. Accessed 15 October 2026. https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/.
25. DepMap, Broad Institute. DepMap Public 24Q4. Figshare+ dataset, version 1. 2024. https://doi.org/10.25452/figshare.plus.27993248.v1.
26. Cell Model Passports. API documentation. Accessed 10 October 2026. https://depmap.sanger.ac.uk/documentation/api/.


## Additional files

**Additional file 1.** `Additional_file_1.csv` (.csv); Table S1: Primary leave-one-tissue-out results across prespecified drug-by-tissue cells. Includes sample-group counts, correlations, paired differences, definedness, and zero-solution status for ten protocol-defined display tissues. Undefined values remain NA.

**Additional file 2.** `Additional_file_2.csv` (.csv); Table S2: All prespecified CTRPv2 contrast estimands. Drug-specific and equal-drug macro estimates for the primary and two matched-candidate comparison families; intervals condition on locked models and shared external-group resampling.

**Additional file 3.** `Additional_file_3.csv` (.csv); Table S3: All inferential contrasts from repeated nested development evaluation. Nine primary-family and 18 matched-candidate representation estimands; intervals condition on saved out-of-fold predictions.

**Additional file 4.** `Additional_file_4.csv` (.csv); Table S4: Absolute performance of each development model arm. Drug-specific and macro Spearman, RMSE, MAE, and R-squared results with conditional group-bootstrap intervals.

**Additional file 5.** `Additional_file_5.csv` (.csv); Table S5: Variation across the three outer cross-validation repeats. Repeat-level estimates and descriptive ranges/standard deviations; repeats are not independent replications.

**Additional file 6.** `Additional_file_6.csv` (.csv); Table S6: Molecular zero-solution and fit audit. Outer-fold zero/nonzero fits, warnings, fallbacks, active-feature ranges, and full-development molecular coefficient counts.

**Additional file 7.** `Additional_file_7.csv` (.csv); Table S7: Complete leave-one-tissue-out sensitivity results. All 184 drug-by-tissue cells across 23 tissues, retaining undefined statistics and model status.

**Additional file 8.** `Additional_file_8.csv` (.csv); Table S8: Compound identity mapping and assay coverage. GDSC2/CTRPv2 identifiers, identity evidence, chemical identifiers, concentration ranges, and cohort coverage for eight compounds.


## Figure legends

Figure 1. Study design and evaluation boundaries. GDSC2 continuous drug-response AUC and CCLE transcript-per-million expression were used for development across eight prespecified compounds. Cell-line identities were grouped when constructing repeated nested cross-validation folds; all supervised candidate selection, representation construction, and model tuning were restricted to the relevant training folds. Eight model arms comprised a tissue reference, a strong 2,000-gene expression model, and two candidate-selection rules crossed with expression, within-sample rank, and gene-pair representations. Full-development models were audited and generated CTRPv2 predictions before external response values were accessed. The development and external evaluations share the CCLE expression source; the external analysis is retrospective cross-screen ranking, not cross-expression-platform or clinical validation.

Figure 2. Primary REO-pair versus strong-expression contrast in repeated nested development evaluation. Points show the drug-specific and equal-drug macro difference in Spearman correlation, defined as the response-association gene-pair model minus the strong 2,000-gene expression model, using out-of-fold predictions from five outer folds repeated three times. Bars show family-9 adjusted conditional bootstrap intervals based on the recorded out-of-fold predictions. Positive values favor the gene-pair model. Intervals are conditional on the recorded predictions and do not include a full repetition of feature selection, tuning, and model fitting.

Figure 3. Candidate-set representation contrasts and full-development molecular solutions. (A) Within each candidate-selection rule, the difference in Spearman correlation between the full gene-pair representation and expression using the same 200 candidate genes; “MAD” denotes selection by median absolute deviation and “response” denotes training-fold response-associated selection. Bars are family-18 adjusted conditional intervals from the repeated nested evaluation. (B) For each drug and molecular model arm fitted on the full development set, gray denotes zero nonzero molecular coefficients and teal denotes at least one nonzero molecular coefficient. This panel describes the final fitted models, not the frequency of zero solutions across outer folds; nonzero coefficients alone do not establish predictive utility or biological mechanism.

Figure 4. Defined leave-one-tissue-out contrasts and locked CTRPv2 cross-screen evaluation. (A) For each drug and each of the ten prespecified held-out tissues, teal marks a defined within-tissue Spearman difference between the response-association gene-pair model and the strong-expression model; gray marks an undefined difference. Counts at the left report the number of defined contrasts out of ten. Only 34 of 80 drug-by-tissue contrasts are defined, so the strict equal-drug macro is undefined. (B) Drug-specific and equal-drug macro differences in external Spearman correlation for the full-development response-association gene-pair model versus the strong-expression model; right-side labels give the number of shared external cell-line groups for each drug. Bars are family-9 adjusted conditional bootstrap intervals based on the locked models and external groups. Three drug-specific differences are exactly zero because both models have zero molecular coefficients. CTRPv2 is a retrospective external screen with the same CCLE expression input, not an independent expression-platform or clinical validation.
