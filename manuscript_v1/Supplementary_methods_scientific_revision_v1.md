# Supplementary methods and synthetic-control results

Prepared from existing locked protocols, source and run metadata on 16 October 2026. This document is Additional file 14 in the scientific revision package. The supporting tables are Table S9 (Additional file 9), external absolute performance; Table S10 (Additional file 10), representation comparisons; Table S11 (Additional file 11), complete synthetic cohort results; Table S12 (Additional file 12), synthetic generation checks; and Table S13 (Additional file 13), execution versions. No model refitting, new resampling or external scoring was performed in preparing this document. The scientific revision manuscript is maintained separately from the earlier working drafts.

## SM1. Identity, cancer eligibility and external identity QC

Canonical identities were reconciled using authoritative exact CCLE, Broad/DepMap, Sanger and COSMIC identifiers. Ambiguous or conflicting mappings were excluded instead of repaired with fuzzy cell-line names. Biological groups joined known canonical identities, Broad PatientID links for models from the same tissue sample, Sanger same-patient links and parent/derivative links. Broad PatientID was not assumed to encode complete patient genealogy. For the natural external cohort, exclusions used all known GDSC2 candidate identities and their related groups, including unresolved-mapping candidates, rather than only the eight selected drug cohorts. The original cohort construction needed no fallback shared-group reserve.

The cohort-v2 cancer overlay retained nonempty Broad disease annotations other than explicit Non-Cancerous or Unknown. Explicit Non-Cancerous models were excluded. Missing, empty or Unknown disease annotations were quarantined unless authoritative cancer provenance for an exact identity was verified. HCC1588 (ACH-001078; CVCL_A351; SIDM01820; KCLB 71588) was retained using concordant cancer-cell-line provenance, including lung squamous cell carcinoma classification. This metadata decision followed a development-only pilot and preceded external response access by the reconstructed workflow. It did not use response magnitudes or model performance, resolve identity conflicts, or alter the eight compounds and dose windows.

Primary CTRPv2 eligibility retained published curves-post-QC membership, response presence, expression availability and the SNP-matched-reference identity status. SNP-unconfirmed and SNP-not-tested records were reserved in an identity-QC sensitivity list and were not thereby declared misidentified. A retained sensitivity list does not imply that a sensitivity model evaluation was performed. The original cohort-v1 feasibility rule used a 50-group external minimum and preferred compounds with at least 80 external groups when at least eight such chemicals were available; all eight were selected from that preferred pool. Cancer cohort-v2 explicitly retained minima of 200 development and 80 primary external groups. These are project feasibility criteria, not a prospective power guarantee. The validated cancer cohort contains 626 development groups (4,967 drug-group responses) and 175 primary external groups (1,185 drug-group responses). Metadata tissue spellings were canonicalized by exact case-insensitive equivalence to the Broad lineage vocabulary and then the Sanger vocabulary, without fuzzy or histology-based merging.

Historical cohort-v1 panel selection additionally required at least four development tissues. Eligible assay candidates were ordered by descending development biological-group count, descending primary external biological-group count, descending development tissue count, and increasing numeric GDSC drug identifier; one primary assay per chemical identity was retained within this coverage ordering. For each GDSC drug identifier, one exact MIN_CONC/MAX_CONC interval was selected by descending eligible biological-group coverage, then descending model coverage, then increasing numeric MIN_CONC and MAX_CONC to resolve ties. These were metadata/presence-based coverage decisions, without response-magnitude or model-performance selection. The cancer cohort-v2 overlay preserved the eight original compound identities, panel order and exact dose windows; it did not rerank or replace drugs or windows according to the current model results.

Sources: `config/cohort_rules.json`; `config/cancer_cohort_v2_rules.json`; `src/curate_cells.py`; `src/freeze_cohorts.py`; `src/freeze_cancer_cohort.py`; `reports/cancer_cohort_v2_validation.json`.

## SM2. GENCODE19 gene mapping and precision

The expression audit used `gencode.v19.genes.v7_model.patched_contigs.gtf.gz` from the Broad CCLE historical directory (SHA-256 `225017a79843a51d3428810bc9656e6353dcce570c816acf02ead8f7d23e8956`). GTF annotations were matched to the complete versioned expression gene identifier. Of 57,820 original unique gene records, 56,238 had an exact full-version GTF match. Unmatched records remained in the source audit with unresolved annotation and were not automatically called noncoding. The model worker selected only records explicitly annotated as protein_coding, giving 19,836 model input genes.

The terminal numeric Ensembl version suffix was removed for stable feature names and deterministic tie-breaking; no gene-symbol aggregation, version-based averaging or expression-column collapse was performed. There were no duplicate unversioned identifiers in the audited source. Original float64 TPM was retained, with round-trip floating-point parsing and exact ties available to the rank/pair builder. The expression source audit checked local format, headers, precision and annotation provenance but did not redownload the complete historical TPM for byte-for-byte comparison with its remote counterpart.

Sources: `src/audit_expression.py`; `data_processed/expression_audit.json`; `data_raw/expression/expression_provenance.json`; `src/run_formal_evaluation.py::worker_init`.

## SM3. Split generation, seeds and rare tissue strata

Every grouped stratified split used scikit-learn StratifiedGroupKFold with shuffle=True, tissue strata and biological-group labels. The repeated development run first sorted the 626-group canonical identity table, generated five outer folds under each repeat seed, and then generated three inner folds in each outer-training context. The global assignment was intersected with each compound's records; all arms shared the same resulting assignment. The final full-development fit similarly generated a single global three-fold assignment before the compound intersection. For LOTO, inner folds were instead generated within each compound's current outer-training rows after the complete held-out tissue was removed. Group separation and exactly-once inner validation coverage were checked in all settings.

| Component | Seed or fixed rule |
|---|---|
| Repeated development outer folds | 20261013, 20261014, 20261015, in repeat order 0, 1, 2 |
| Repeated development inner folds | repeat_seed + 100 + zero-based outer_fold (outer_fold 0–4) |
| Final full-development three-fold tuning | 2026101501 |
| LOTO three-fold tuning | 2026101401 + zero-based index in the tissue order below |
| Development shared group bootstrap | 2026101602; 10,000 draws |
| External shared group bootstrap | 2026101502; 20,000 draws |

| LOTO tissue index | Held-out tissue | Inner seed |
|---:|---|---:|
|0|Biliary Tract|2026101401|
|1|Bladder/Urinary Tract|2026101402|
|2|Bone|2026101403|
|3|Bowel|2026101404|
|4|Breast|2026101405|
|5|CNS/Brain|2026101406|
|6|Cervix|2026101407|
|7|Esophagus/Stomach|2026101408|
|8|Head and Neck|2026101409|
|9|Kidney|2026101410|
|10|Liver|2026101411|
|11|Lung|2026101412|
|12|Lymphoid|2026101413|
|13|Myeloid|2026101414|
|14|Ovary/Fallopian Tube|2026101415|
|15|Pancreas|2026101416|
|16|Peripheral Nervous System|2026101417|
|17|Pleura|2026101418|
|18|Prostate|2026101419|
|19|Skin|2026101420|
|20|Soft Tissue|2026101421|
|21|Thyroid|2026101422|
|22|Uterus|2026101423|

Some tissue strata contained only one or two groups. The repeated development split preparation recorded three warnings for the five-fold outer splits and fifteen warnings for the three-fold inner splits. The final global three-fold preparation recorded one warning for a one-group stratum. The LOTO preflight recorded a warning in each of the 184 drug-by-held-out-tissue contexts: one-group strata in 176 contexts and two-group strata in eight. The fixed grouped folds were retained without removing rare tissues or changing assignments on the basis of performance. These are splitter warnings, distinct from the zero solver warnings and fallbacks recorded for the relevant completed real-data fitting stages; they limit the balance and tuning evidence available for rare tissues. Synthetic fit warnings are reported separately in SM8.4.

The original LOTO protocol prose recorded a conflicting seed base of 2026101801, whereas the machine configuration recorded 2026101401. Before any LOTO model fitting, execution amendment v1.1 selected the explicit machine value and published the tissue-specific map. Execution v1.2 repaired a resource-logging CSV schema defect in v1.1, which had launched no fit task. Execution v1.3 added the already locked inner_seed field to OOF identity records after the first v1.2 tissue-null task was rejected by the checkpoint contract. That preserved v1.2 attempt contained one tissue-null fit and no molecular fit. The completed sensitivity analysis used v1.3; the logging and row-identity repairs did not change its model choices, seeds, split rules, input data, estimands or metrics. The cohort JSON format's schema_version 2.0 is a separate metadata-overlay identifier, not a LOTO execution version.

Sources: `src/modeling/nested_cv.py::make_grouped_splits`; `src/run_formal_evaluation.py::make_global_split_records`; `src/final_model_v1.py::make_global_inner_manifest`; `src/run_tissue_loto_v1_3.py::compute_split`; `config/loto_execution_v1_3.json`; `docs/LOTO_EXECUTION_AMENDMENT_V1_1.md`; `docs/LOTO_EXECUTION_AMENDMENT_V1_2.md`; `docs/LOTO_EXECUTION_AMENDMENT_V1_3.md`; stage-specific run metadata.

## SM4. Original-scale prediction and unseen tissues

Let z be the molecular feature vector after the training-fitted standard scaler, d the tissue one-hot row, M the training tissue-specific means of scaled molecular features, a the training tissue-specific means of the standardized response, beta the fitted molecular coefficients, and mu_y and s_y the training response mean and population standard deviation. The implemented prediction is

`y_hat = mu_y + s_y * [(z - d M) beta + d a]`.

The molecular coefficients are fitted on the training tissue-residualized molecular features and response. The tissue encoder uses handle_unknown='ignore', so an unseen tissue has d=0. Its original-scale prediction is therefore the overall training response mean plus the molecular contribution, not zero AUC. If beta=0, prediction equals the training tissue mean for an observed tissue and the overall training mean for an unseen tissue, exactly matching the same-policy tissue reference. A retained nonzero molecular contribution alone establishes neither predictive utility nor biological mechanism.

Source: `src/modeling/nested_cv.py::_fit_elastic_net` and `predict`; `src/ordinal_diagnostic_models.py::DiagnosticRegressor`.

### SM4.1. External input schema, abstention, hard failures and scoring denominators

The external gateway used the locked primary roster and the same versioned full float64 TPM input, with all 19,836 stable-gene columns in the exact frozen order. It checked the full input schema, rather than only the selected genes of a particular model. Missing, extra, duplicated or reordered gene columns were structural failures; duplicate or unexpected expression-row identifiers, duplicate/conflicting canonical identities, an identity or biological-group overlap with the forbidden training/candidate pool, or an inconsistent source/input/model hash also caused a hard failure. The frozen contract additionally treated wrong compound mapping, an unexpected endpoint field or dose-window convention, and nonfinite output from an otherwise valid fitted-model prediction as hard failures. These failures stopped evaluation rather than being repaired with test-set information.

With a valid column schema, a missing expression row or any missing, nonnumeric, nonfinite or negative value anywhere in a sample's full raw TPM vector caused all seven molecular arms to abstain. This policy also applied to molecular arms whose fitted coefficients were zero and to invalid values in genes that a particular arm had not selected. No missing value was converted into a tie, imputed, or rescued according to model performance. The tissue-only reference could still predict when canonical identity, biological-group identity and a nonmissing valid tissue label were available, because its prediction did not depend on raw molecular values. Any internal molecular placeholders passed to that reference were ignored, not treated as observed expression. An invalid biological-group identifier or missing/invalid tissue label instead prevented prediction by every arm. An otherwise valid, nonblank tissue label absent from the individual drug model's training tissue map remained evaluable under the unseen-tissue policy in SM4. A legitimate zero molecular solution was a valid prediction state, not an abstention or a failed fit.

Each frozen compound-by-model-by-arm row remained in the eligible prediction denominator, with a status, reason and zero-solution/unknown-tissue flags, regardless of whether a prediction was possible. Response validity was tracked separately. For the response aggregation, all selected source experiments for a compound and canonical model had to contain finite parseable AUC values; if any selected value was invalid, that group-level response was marked invalid rather than computed from a favorable subset. Otherwise, the median of all selected eligible experiment values supplied its one response label. Source experiment/compound identities and row pointers were checked against the frozen roster.

For each contrast, both correlations were calculated using the same intersection of response-valid canonical groups and finite, predicted rows from both arms. Coverage records retained frozen eligible N, response-valid N, predicted N over the frozen roster, each arm's predicted N among response-valid groups, paired common-ID/group N, and noncoverage reasons. Missing/invalid responses and abstentions did not silently shrink the original eligibility denominator or the fixed eight-drug panel. Fewer than 80 paired common groups was a low-coverage flag, not a replacement drug or an automatic redefinition of the estimand; correlation definedness followed SM5. If a prespecified drug contrast was undefined, the strict eight-drug macro remained undefined. In the completed locked evaluation all 1,185 planned drug-group responses were valid and all 9,480 arm predictions were finite, with no abstentions; these observations do not remove the need to report the frozen rules.

Sources: `config/external_validation_v1_1.json::external_gateway`; `src/external_validation_v1.py::prepare_gateway`, `validate_model_contract`, `predict_external_arm`, and `paired_common_table`; `src/run_external_validation_v1.py::read_frozen_curve_responses`, `_arm_descriptive`, and `_paired_even_if_no_valid_response`.

## SM5. Correlation definedness and conditional bootstrap intervals

The development bootstrap generator sorts the frozen biological-group labels, then calls NumPy Generator(PCG64(seed)).multinomial with equal group probabilities to produce a shared matrix of integer multiplicities. All 626 groups are included in the development sampling universe; the same draw applies across drugs, arms and all three repeats. Within each drug/repeat the implementation recomputes exact tied midranks under these multiplicities. A development draw's Spearman value requires total sampled multiplicity at least two and positive weighted rank variance of both response and prediction. Paired contrasts and repeat/drug means propagate undefined components. This implementation does not apply the additional external minimum of three distinct represented groups.

External point correlations require aligned finite values, unique biological-group identities, at least three independent common groups, and nonconstant responses and predictions. External draws use the 175-group universe, the same PCG64-multinomial generator and shared multiplicities, with ranks recomputed as above. Every arm correlation in a draw additionally requires at least three distinct groups with positive multiplicity. Reasons for invalid components are retained; any undefined prespecified drug component makes that draw's strict eight-drug macro undefined. The LOTO descriptive helper separately uses NA for fewer than two held-out samples, constant response, constant prediction or nonfinite correlation, and supplies no LOTO confidence interval.

Intervals require a defined point and at least 99% valid draws: 9,900/10,000 for development and 19,800/20,000 for external evaluation. Undefined or insufficient-validity intervals remain NA. No extra draws replace invalid draws. Percentile endpoints use `numpy.quantile(..., method='linear')`; unadjusted limits are 0.025 and 0.975, and each family-adjusted two-sided limit uses tails 0.05/(2K), with K=9 for the primary family or K=18 for the matched-candidate family. Family adjustments remain separate and do not imply study-wide coverage. All 27 prespecified development estimands and 27 prespecified external estimands in the existing results have all requested draws valid. No group bootstrap is presented as full feature-selection, tuning or training uncertainty; no equivalence or non-inferiority test was performed.

The external protocol requested exact NumPy and SeedSequence fields in its bootstrap output. The original completed summary instead retained the seed, group order and multiplicity hash without those explicit fields. NumPy 2.3.5 is identified by the enforced pre-fit environment lock, and the frozen source uses Generator(PCG64(seed)); stage 14 records a later seed-sequence reconstruction as retrospective. It must not be described as an original runtime log field, and the completed summary is not rewritten to suggest otherwise.

Sources: `src/analyze_formal_results.py::generate_global_bootstrap_weights`, `weighted_metrics`, `interval_summary`; `src/analyze_ordinal_repeated.py`; `src/external_validation_v1.py::spearman_checked`, `_interval`, `bootstrap_paired_rho`; `src/analyze_tissue_loto_v1_3.py::spearman`; `config/ordinal_repeated_v1.json`; `config/external_validation_v1_1.json`; `reports/阶段14_P6一次性CTRP外验与研究方向判定.md`.

## SM6. Scientific execution versions and recorded environments

The stage-specific scientific commits, recorded software versions, metadata hashes and execution-lock hashes are given in proposed Table S13 (`tables/scientific_revision_v1/TableS13_execution_versions.csv`). They identify the source snapshots checked by the corresponding execution contracts, not a claim that the working documentation HEAD never changed during execution. Development fitting and the separate conditional-bootstrap analysis share the recorded scientific commit, but the bootstrap provenance does not separately log its interpreter/package versions. The fitting environment is therefore not misrepresented as an independently captured analyzer environment.

The completed LOTO metadata records Python, NumPy, pandas, SciPy, scikit-learn and psutil versions. It does not separately record threadpoolctl or joblib versions; the table explicitly marks those fields not_recorded rather than borrowing versions from the development or final-model stage. The final-model stage bound the environment in its first lock and compared installed versions with that lock; the external gateway repeats that validation before response access. The one-time external summary is bound to the second lock and frozen source snapshot, rather than carrying a newly invented stand-alone environment log.

## SM7. Finite synthetic implementation controls: methods

### SM7.1. Purpose and execution provenance

Four prespecified data-generating processes (DGPs) were used to illustrate the behavior of the evaluated candidate-selection and representation pipelines under explicitly constructed signals. The controls were not designed to estimate the prevalence or causes of zero molecular solutions in real drug-response data, establish a general advantage of REO, or replace label-permutation controls for the eight real compounds. Ten separately generated cohorts were evaluated for each DGP, giving 40 cohorts. Each cohort contained 240 independent synthetic sample groups, eight tissue labels with 30 samples each, and 220 synthetic genes. These controls were frozen under execution protocol P5e-synthetic-v1.1 before any synthetic model fit. The original v1 configuration was retained; the v1.1 amendment fixed the support of random variates and clarified the outer/inner split seeds. The executed source commit was `d151b668f916fda643e90ba1a26035cd103b9036`, with execution-lock SHA-256 `c149bb2d41d3d1a13c874e99ac73c087b25f12d4d5684c32582ff54d461043fd`.

### SM7.2. Data generation

Let t(i) be the tissue index of sample i, taking values 0–7; let h(i) = [t(i) − 3.5]/3.5. The gene index j takes values 0–219. Define the baseline b(j) = 3 + 0.008 × j, and let ℓ(j) be 220 equally spaced values from −1 to 1. The generated expression-scale variable is L(i,j) = log₂[X(i,j) + 1], with synthetic TPM input X(i,j) = 2 raised to L(i,j), minus 1. The original implementation uses the variable name `log2_tpm`; its mathematical meaning follows this conversion and is therefore log₂(TPM + 1).

Independent standard-normal draws used for ε, η, and a were rejection-truncated to [−3,3]. Within a cohort, draws followed the frozen NumPy Generator(PCG64(seed)) stream: arrays were generated in C row-major order and out-of-support elements were redrawn in increasing flattened-index order until all were in support. All outcomes were clipped to [0.01,0.99]; expression was not subsequently clipped or imputed.

The four DGPs were defined as follows.

**Tissue-only signal.**

L(i,j) = b(j) + 0.20 × h(i) × ℓ(j) + 0.35 × ε(i,j).

y(i) = 0.50 + 0.28 × h(i).

The outcome is exactly constant within each tissue; expression noise does not enter the outcome.

**Weak molecular signal.**

L(i,j) = b(j) + 0.15 × h(i) × ℓ(j) + 0.35 × ε(i,j).

z(i) = √5 × mean[ε(i,0), ε(i,1), ε(i,2), ε(i,3), ε(i,4)].

y(i) = 0.50 + 0.05 × h(i) + 0.04 × z(i) + 0.14 × η(i).

**Common sample amplitude.**

L(i,j) = b(j) + 0.20 × a(i).

y(i) = 0.50 + 0.04 × h(i) + 0.12 × a(i) + 0.02 × η(i).

Equivalently, X(i,j) + 1 = 2^[b(j)] × 2^[0.20 × a(i)], where the caret denotes exponentiation of the complete bracketed term. The sample-specific factor therefore scales TPM + 1, rather than TPM itself. The within-sample gene order, and hence both rank and complete pair representations, is identical across samples by construction. Generation checks recorded this exact invariant for all 220 genes; this scenario is an intentional information-removal control.

**Predefined pair-order signal.**

Initially, L(i,j) = b(j) + ε(i,j). The first two genes were then replaced with:

L(i,0) = 3.5 + 0.35 × ε(i,0).

L(i,1) = 3.5 + 0.35 × ε(i,1).

y(i) = 0.50 + 0.04 × h(i) + 0.15 × sign[L(i,0) − L(i,1)] + 0.03 × η(i).

Here sign returns −1, 0 or +1 according to its argument. This deliberately places a response component in the ordering of two lower-variance driver genes. It is a constructed selector/representation stress case, not a representation-neutral benchmark.

Stable synthetic gene identifiers were `SIMGENE_0001` through `SIMGENE_0220`. For scenario indices s = 0, 1, 2, 3 in the order above and replicate indices r = 0–9, the cohort seed was 2026101400 + 100 × s + r. The exact 40 seeds, ranges, dimensions, and cohort/split file hashes are retained in Table S12. Blank invariant fields outside the common-amplitude scenario mean not applicable, not failed.

### SM7.3. Model comparison and nested evaluation

Eight arms were retained: a tissue-only reference, the strong expression elastic-net baseline, and two candidate-gene rules (MAD and within-training-tissue response association) crossed with continuous expression, within-sample rank, and complete three-valued pair encoding. The same training-only detection, candidate construction, scaling, unpenalized tissue projection, and relative-alpha elastic-net policies as in the development comparison were used. The matched arms selected 200 genes and encoded 200 continuous or rank features, or 19,900 pair features. Although the strong-expression budget was 2,000, only 220 genes were generated; all 220 were retained in its outer training fits. It must therefore be interpreted as the 220-gene strong-expression arm in this control experiment. No random-forest arm was fitted.

Each cohort used one three-fold outer evaluation and three-fold inner tuning. Splits used shuffled StratifiedGroupKFold, stratified by tissue, with one unique biological group per synthetic sample. The outer seed was cohort_seed + 1000; the inner seed was cohort_seed + 1100 + zero_based_outer_fold on the outer-training samples. Every arm shared the cohort's splits. Each outer test fold contained 80 samples, with 160 in its training set. Feature selection and all fitted preprocessing were repeated within the relevant training partition. The ten inner candidates combined relative alpha values {1, 0.3, 0.1, 0.03, 0.01} with L1 ratios {0.5, 1}; the pooled inner OOF RMSE selected the candidate, retaining the analytical zero endpoint and the frozen exact-tie rule. Each molecular outer task contained 30 inner fits and one outer refit; each tissue-reference task required one fit. The complete design consequently contained 960 outer tasks, 26,160 fit events, and 76,800 sample-by-arm OOF predictions. The ten DGP cohorts per scenario are generated cohorts, not ten repeats on a single cohort.

### SM7.4. Descriptive endpoints, undefined results and operational records

Each cohort-by-arm row used its complete 240-sample OOF predictions to calculate Spearman correlation, RMSE, MAE, and R-squared. Molecular increments were additionally characterized using the OOF tissue reference:

e(i) = y(i) − tissue_reference_prediction(i).

q(i) = model_prediction(i) − tissue_reference_prediction(i).

Spearman correlation between e and q was recorded both pooled and within each tissue. A within-tissue correlation was undefined if either vector was constant. The strict eight-tissue mean required all eight tissue correlations to be defined; the available-tissue mean and its coverage were recorded separately. RMSE gain was tissue-reference RMSE minus model RMSE, so positive gain favors the model. No undefined correlation was replaced by zero or removed from Table S11.

Existing scenario summaries report the mean, median, minimum, maximum, and number of defined values across ten cohorts. They are descriptive; no p-values, inferential confidence intervals, power claims, or new label-permutation estimates were calculated. Selected outer-refit zero solutions are distinguished from zero coefficients among all inner-grid and outer fits: the latter denominator includes deliberately evaluated analytical-null endpoints. All warnings, fallbacks, failed-task status, and operational resources were retained in the original execution records. Execution used one isolated task process at a time and one numerical-library thread; runtime and sampled memory are operational observations, not isolated algorithm-efficiency comparisons. Actual-runtime package versions were not separately recorded in the synthetic run metadata/lock; the reproducible source is identified by its execution commit and hashes, and repository dependency pins should not be presented as a newly verified runtime snapshot.

## SM8. Finite synthetic implementation controls: results

### SM8.1. Coverage and retained undefined values

All four prespecified DGPs and all ten generated cohorts per DGP were retained. The complete design comprised 40 cohorts, 960 outer tasks, 26,160 fits, and 76,800 recorded OOF predictions. Existing execution metadata recorded one completed session, 960 completed tasks, no archived failed attempts, and no fallback. The revision checks confirmed complete 320-row cohort-by-arm and 960-row outer-task grids, with no omitted model or cohort. Table S11 retains every cohort-by-arm result, including 235 undefined values in each of its pooled, available-tissue, and strict tissue-macro residual-correlation columns. Point Spearman, RMSE, MAE, and R-squared were defined in all 320 rows. These counts describe the recorded run and do not replace direct validation of individual predictions.

### SM8.2. Tissue-only signal

In the tissue-only scenario, the response was exactly constant within each tissue in all ten generation checks. All 210 molecular outer refits selected zero molecular coefficients, and their OOF predictions retained the tissue component. Mean OOF Spearman was 0.9997 for every arm; the maximum cohort RMSE was 1.962615573354719 × 10⁻¹⁷. Thus the tissue outcome was reconstructed to floating-point precision, while exact recorded rank correlations should not all be replaced by one. Molecular-increment correlations were undefined, as expected for constant outcome/molecular residuals. This is a control illustrating that an extremely high overall correlation can coexist with no learned molecular coefficient.

### SM8.3. Weak molecular signal

In the weak-molecular scenario, the mean OOF Spearman ranged from 0.1367 to 0.1451 across molecular arms, compared with 0.1373 for the tissue reference; mean RMSE ranged from 0.1479 to 0.1485, compared with 0.1479 for the reference. The molecular arms selected zero outer solutions in 198 of 210 refits. Both pair arms selected zero in all 30 outer refits each. These descriptive outcomes illustrate limited recovery under this specific weak-signal, sample-size, and regularization setting; they do not establish an empirical null distribution or prove that corresponding real-drug data have the same mechanism.

### SM8.4. Common sample amplitude

The common-amplitude DGP preserved the exact within-sample rank and complete pair signs across samples in all ten cohorts, by design. All four ordinal arms selected zero molecular coefficients in all 30 outer refits per arm and matched the tissue reference (mean Spearman 0.1541, mean RMSE 0.1216). Each continuous-expression arm used a nonzero molecular solution in all 30 outer refits and had mean Spearman 0.9830 and RMSE 0.0207. The original fit diagnostics recorded 267 warning occurrences, all in this DGP: 89 each for the strong-expression, MAD-expression, and response-associated-expression arms. All outcomes and warning counts are retained; no warning type or failure was inferred from the aggregate count. This deliberately order-invariant DGP illustrates a type of signal unavailable to rank/pair representations; it does not quantify amplitude information loss in real drug responses.

### SM8.5. Predefined pair-order signal

In the predefined pair-order DGP, both driver genes were retained in all 30 outer candidate sets under response-associated selection and in none of the 30 MAD candidate sets; the 220-gene strong-expression arm also retained both drivers in all 30. With response-associated candidates, mean Spearman was 0.7573 for expression, 0.7605 for ranks, and 0.9055 for pairs, with mean RMSE 0.1045, 0.1032, and 0.0324, respectively. With MAD candidates, the corresponding mean correlations were 0.2147, 0.2069, and 0.2112, close to the tissue reference's 0.2147. The strong-expression arm had mean Spearman 0.7572 and RMSE 0.1045. This is a constructed case in which candidate selection determines access to two lower-variance order drivers; the apparent pair advantage is specific to that DGP and model class and is not independent biological validation or a causal decomposition of real-data representation differences.

### SM8.6. Zero-solution denominators, resources and interpretation

Across all scenarios, 615 of 840 selected molecular outer refits were zero solutions. By scenario, the counts were 210/210 (tissue only), 198/210 (weak molecular), 120/210 (common amplitude), and 87/210 (predefined pair order). The distinct all-fit denominator contained 13,575 zero molecular fits among 26,040 molecular grid/refit events, including the deliberately evaluated analytical-null candidates. Those all-fit fractions must not be described as selected-model zero rates. A separate existing candidate-set audit found identical 200-gene sets for expression/rank/pair arms in all 240 same-selector outer contexts; the strong-expression arm used all 220 available genes.

The original controller session lasted 3333.73129409994 seconds; sampled controller/worker-tree peak RSS was 454250496 bytes and peak private memory was 490225664 bytes. These measurements include process management and other execution overhead. All conclusions above are descriptive illustrations of specified synthetic mechanisms. They do not identify the cause of real-data zero solutions, establish a universal representation ranking, evaluate single-cell foundation models, or substitute for real-label full-pipeline permutation controls.

All 32 scenario-by-arm summaries, 32 fit-diagnostic rows, 32 gene-selection rows, 960 outer-task rows and 240 candidate-set checks remain in the repository. Their source paths and hashes are listed in `synthetic_reportability_audit.json`. No scenario was omitted.
