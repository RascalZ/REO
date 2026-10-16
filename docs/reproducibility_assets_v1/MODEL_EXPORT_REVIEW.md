# Existing final-model export: source-only implementation review

Prepared: 2026-10-16. Scope: inference/export design, not model fitting or efficacy evaluation.

This review read frozen Python source, the existing public final-model manifest, the expression gene annotation header/order, and the external protocol configuration. It did **not** deserialize any model, open CTRPv2 responses or response indexes, inspect individual external predictions, train, predict, resample, or score. Consequently, the field names and arithmetic below are source-verified expectations; the exporting process must validate the actual object fields and dimensions before claiming that all 64 payloads conform.

## 1. What the final panel actually contains

`results/p6_preflight/final_model_manifest.jsonl` contains 64 rows: eight drugs × the following eight arms. All 64 rows report fresh-process reload verification, and no row has a recorded fallback reason. Its SHA-256 at this review is `0491c6faca2b79c74f5c056981ddec0263683600cbc3bb5d24c0f0bd160a86ad`.

| Arm | Selection policy | Actual representation | Candidate budget |
|---|---|---|---:|
| tissue_null | none | tissue_null | 0 |
| expression_strong | mad | expression | 2,000 |
| mad_expression | mad | expression | 200 |
| mad_rank | mad | rank | 200 |
| mad_pair_full | mad | pair_full | 200 |
| response_association_expression | response_association | expression | 200 |
| response_association_rank | response_association | rank | 200 |
| response_association_pair_full | response_association | pair_full | 200 |

These are `DiagnosticRegressor` models using the inherited FWL elastic-net arithmetic, **not** RF artifacts. The manifest reports 29/56 molecular models with exactly zero molecular coefficients and 27/56 with nonzero molecular coefficients. The eight explicit tissue-null arms are a separate denominator; they must not inflate the molecular zero-solution rate. Reported fitted kinds are expression 24, rank 16, pair_full 16, tissue_null 8. Preserve all 64 models and the zero state of every one.

Sources: `src/ordinal_repeated_models.py:model_spec`, `src/final_model_v1.py:fit_full_development_model`, `src/ordinal_diagnostic_models.py:DiagnosticRegressor`, and the manifest above.

## 2. Artifact and fitted-state field inventory

`src/run_final_model_v1.py` creates a dictionary artifact (lines 463 onward) with the keys below. Export only SHA-bound trusted local artifacts; `joblib` is not a safe interchange format for untrusted downloads. The independent public inference path should load portable JSON, not pickle/joblib or sklearn classes.

| Existing artifact key | Export requirement |
|---|---|
| schema_version | Record the original version and a distinct portable-export schema version. |
| task | Preserve drug ID, arm, selection_rule, representation, n_genes. |
| freeze_sha256 | Preserve original frozen scientific provenance. |
| model | Extract the numeric/static fitted state listed below; do not export a live Python object. |
| stable_gene_ids | Preserve all original 19,836 IDs in exact input order. |
| tissue_code | Preserve the exact full-development tissue-label → numeric-code dictionary. |
| expression_genes_sha256 | Preserve the source annotation hash. |
| expression_tpm_sha256 | Preserve input provenance; no need to load this expression matrix for synthetic-input checks. |
| development_prediction_sha256 | Retain as historical provenance; do not recompute it on real observations in this round. |

Relevant model attributes (the names include trailing underscores):

| Existing attribute | Required meaning and expected shape |
|---|---|
| model.model_kind / model.selection_rule | Constructor declaration, distinct from fitted_kind_. |
| model.fitted_kind_ | prediction branch: expression, rank, pair_full, tissue_null; fail closed on unhandled branches. |
| model.n_features_in_ | 19,838 = 19,836 genes + tissue code + audit ID. |
| model.stable_gene_ids | Same input-gene tuple as artifact.stable_gene_ids. |
| model.target_mean_ / model.target_sd_ | Training-response mean and population SD (`np.std`, ddof=0); scalars. |
| model.encoder_.categories_[0] | Exact ordered per-drug observed tissue codes; size T. |
| model.training_tissue_codes_ | Should agree exactly with categories_[0]. |
| model.target_tissue_means_ | Mean **standardized response** by observed tissue; shape (T,). Not raw-AUC tissue means. |
| model.builder_.kind / selection_rule | Must agree with the task representation/policy. |
| model.builder_.n_features_in_ | 19,836 full TPM input columns. |
| model.builder_.gene_indices_ | Selected positions in the full TPM schema; length G. |
| model.builder_.gene_ids_ | Selected stable IDs in output order; length G. |
| model.builder_.selection_ranked_gene_ids_ | Score-ranked candidate list before ID-based output ordering. Useful selection provenance, not inference order. |
| model.builder_.selected_scores_ | Training selection scores, aligned to ID-sorted gene output order. |
| model.builder_.selected_signed_associations_ | Signed training associations, same order; MAD arms retain zero entries. |
| model.builder_.pair_indices_ | Shape (P,2), full-input gene positions; **not** positions in selected G genes. |
| model.builder_.feature_names_ | Ordered encoded-feature names; length P for pairs or G for expression/rank. |
| model.builder_.n_output_features_ | Encoded feature count p. |
| model.builder_.audit_metadata_ | Retain selection, encoding, no-pair-filter and tie rules. |
| model.scaler_.mean_ | Training encoded-feature means; shape (p,). |
| model.scaler_.scale_ | Exact fitted training scale factors; shape (p,), including the scaler's handling of constant/near-constant columns. |
| model.scaler_.var_ | Retain for audit; shape (p,). Use scale_, not a freshly recomputed sqrt(var_), at inference. |
| model.feature_tissue_means_ | Means of **standardized encoded features** by tissue; shape (T,p). This is not a single p-vector. |
| model.regressor_.coef_ | Standardized-feature/standardized-target EN coefficients; shape (p,), including every exact zero. |
| model.alpha_ratio / model.l1_ratio | Selected constructor hyperparameters. |
| model.actual_alpha_ / alpha_max_ / path_ratio_ | Realized full-development regularization state; retain without recomputing it. |
| model.null_endpoint_ / solver_reason_ | Distinguish the analytic alpha-max endpoint from other zero reasons. |
| model.n_nonzero_molecular_coeffs_ | Must equal count(coef_ != 0.0), with no tolerance. |
| model.fallback_reason_ / warning_messages_ / exact_tissue_constant_target_ | Preserve scientific outcome and numerical records. |

Explicit tissue-null models have `builder_ is None` and no molecular scaler, feature projection, or regressor state. Do not invent empty arrays and imply that they were learned. Store absent molecular state as null and preserve the declared same-policy tissue-null branch.

The generic source also supports a `mean` branch after a constant-training-response fallback and a tissue-null fallback after an empty encoded feature set. Neither is reported in the current 64-row manifest. Export validation should either implement those separately with an explicit state or reject any unexpected object; it must not silently reinterpret them as the ordinary EN branch.

## 3. Exact input schema and encoding

The frozen loader restricts the original annotated expression matrix to `gene_type == 'protein_coding'`, retaining that original row/column order and using `ensembl_id` (stable ENSG IDs). It does not remap by gene symbol, merge duplicate IDs, average probes, or infer missing genes. At review there are exactly 19,836 unique stable IDs in this list. They happen to be lexicographically sorted in the current annotation; this is a fact to verify against the exported schema, **not** permission to independently re-sort arbitrary input and silently accept it.

For a public input interface, require the exact expected ID list and declared order. Missing, extra, duplicated, or reordered columns are a hard schema failure under the frozen protocol. After complete schema verification, check **all** 19,836 raw TPM values in every molecular input row, including unselected genes: finite, numeric, nonnegative, raw TPM, float64. Do not infer that log expression is TPM and do not treat missing values as zero or exact ties.

Selected candidate genes are sorted by stable ID for **all three final representations**, including the 2,000-gene strong expression arm. This differs from the earlier generic `FoldFeatureBuilder` expression arm, which retains decreasing-MAD order. Use the final `DiagnosticFeatureBuilder` state, not the old generic builder convention.

Encoding is fixed, sample-local, and precedes feature standardization and tissue projection:

1. **expression:** `F = log2(raw_TPM[:, gene_indices_] + 1.0)`; shape (n,G), float64.
2. **rank:** `F = rankdata(raw_TPM[:, gene_indices_], axis=1, method='average')`; shape (n,G), float64. Ascending, one-based ranks; exact equal TPMs receive their average midrank. Rank within each sample's **entire selected candidate set**, not across samples, all 19,836 genes, or only genes with active coefficients.
3. **pair_full:** ordered upper-triangular pair list over selected ID-sorted genes. For each exported `(a,b)` full-schema pair, `F_j = (TPM_a > TPM_b) - (TPM_a < TPM_b)`, giving exact -1/0/+1. The first ID is lexically smaller; `first > second` gives +1. Compare, do not subtract TPM values or apply an epsilon. An exact tie gives 0. All pairs, including training-constant and high-tie pairs, are retained. With G=200, p=19,900.

Do not apply the old `FoldFeatureBuilder`'s constant/high-tie pair filters to final `pair_full` models. Do not rerun detection, MAD ranking, response association, or any candidate selection at inference. Selection metadata are provenance, not a new adaptive decision.

## 4. FWL prediction order and tissue handling

For n rows, p encoded features, T observed training tissues, let:

- F be the encoded matrix;
- c = scaler.mean_, s = scaler.scale_;
- Z = (F-c)/s, exactly in feature order;
- D be n×T full one-hot matching the exported ordered encoder categories;
- M = feature_tissue_means_, shape T×p;
- a = target_tissue_means_, shape T;
- beta = regressor.coef_, shape p;
- mu = target_mean_, sigma = target_sd_.

The existing EN prediction computes, in this order:

```text
Z = (F - c) / s
Z_res = Z - D @ M
standardized_prediction = Z_res @ beta + D @ a
prediction = mu + sigma * standardized_prediction
```

There is no separately fitted molecular intercept (`ElasticNet(fit_intercept=False)`), no prediction clipping, no external standardization, and no test-derived tissue adjustment. Do not collapse this to a different intercept formula if the objective is close numerical agreement with the original implementation: mathematically equivalent rearrangements can change floating-point reduction order.

For a valid observed tissue, its one-hot column is 1. For a valid unseen tissue, **all** D columns are zero: the prediction is `mu + sigma*(Z @ beta)`. For an exact molecular zero solution, this becomes the training global response mean for unseen tissue, and `mu + sigma*a_t` for an observed tissue. A standardized unseen-tissue contribution of zero does **not** mean a prediction of zero AUC.

The artifact's `tissue_code` dictionary covers the development union, but each drug's encoder contains only tissues observed for that drug. Determine the unknown-tissue flag using the **per-model encoder categories**, not merely membership in the global dictionary. A valid global tissue absent from one drug is unseen to that model. A truly new valid tissue must receive a noncolliding unknown code or be translated directly into an all-zero D row. A blank/missing/invalid tissue is an abstention, not an unseen tissue to rescue. The interface must distinguish a validated new label from a typo or malformed value; the documented vocabulary/identity contract supplies that validation.

For `tissue_null`, predict only:

```text
prediction = mu + sigma * (D @ a)
```

Tissue-null uses no raw TPM in its arithmetic. With valid identity and tissue, a sample that fails molecular TPM checks can still receive its tissue-null prediction. Conversely, all seven molecular arms—including models whose beta is entirely zero—must abstain on a row with any invalid full TPM value under the frozen gateway policy. Do not add a zero-model exemption. Zero is a valid fitted state, not abstention, failed fit, evidence of drug inactivity, or equivalence/noninferiority.

The low-level estimator's predict method is not the entire frozen gateway: in particular, the tissue-null branch does not validate TPM, and it does not independently implement all identity/missingness rules. The new portable API must explicitly implement the gateway instead of describing the low-level method as already sufficient.

## 5. Export validation and numerical-probe checklist

Before reading a payload, bind its path and SHA-256 to the existing 64-row manifest and original lock. Fail on unlisted paths, duplicate drug/arm keys, missing models, unexpected branches, or source/schema mismatch. Do not call a training controller as an exporter: `run_final_model_v1.py` has training, response loading, checkpoint mutation and legacy verification paths; a new prediction-only exporter should not invoke its main/worker functions.

For each payload validate all dimension/order invariants, finite scalars/arrays, positive scale factors and positive target SD. In particular:

- artifact.stable_gene_ids and model.stable_gene_ids agree, and their length is 19,836;
- model.n_features_in_ is 19,838 and molecular builder.n_features_in_ is 19,836;
- selected indices are unique/in-range and resolve exactly to gene_ids_;
- selected IDs and feature names reproduce existing NUL-separated SHA-256 values;
- pair indices map to the correct exported names and cover the full upper triangle without reorientation;
- p equals feature-name count and coefficient/scaler vector lengths; M has shape (T,p), a has length T;
- encoder categories equal training_tissue_codes_, and exact numeric codes resolve to the dictionary's tissue labels;
- coefficient zero count/state agrees with the manifest, rather than being thresholded or rounded;
- preserve scaler.scale_ exactly, even when var_ is zero or very small;
- JSON must preserve finite float64 values with round-trip precision, without rounding for readability.

Use generated TPM inputs only, not development sample rows or external expression/response. Include deterministic cases for ordinary positive values, zeros, ties, near ties, common amplitude changes, valid observed and unseen tissues, and the whole-panel zero-model subset. Compute predictions through both the existing trusted loaded model and the independent JSON implementation, with no labels, fit calls, metrics, refit, or scientific result update. A tolerance such as absolute 1e-12 / relative 0 follows the original sentinel contract; report any failure instead of relaxing tolerance after observing it. Distinguish a finite input-generation check from the numerical prediction agreement itself.

The original `make_reload_sentinels` deliberately copied three **development** rows, inserted one unseen tissue, and made one all-zero TPM row. It was a historical reload test; do not reuse it in this round when promising entirely synthetic input and no original observations. Independent synthetic inputs need not reproduce its saved digest.

Separately test gateway behavior (schema mismatch, nonfinite/negative full TPM, malformed tissue, missing identity) without calling the original predictor on invalid input. Preserve row IDs, reason codes, predicted/abstained states and unknown/zero flags. An input validation test does not require labels. Do not compute Spearman, RMSE, MAE, R² or reuse the external scorer.

Do not publish sample-level inference outputs from genuine CTRP rows as an example. A small deterministic synthetic-input fixture and generated prediction CSV can demonstrate use without exposing or rereading the sealed labels. State explicitly that the values are software checks, not new biological findings.

## 6. Practical risks and reporting limits

- **Export size:** full pair-model M matrices dominate the portable files (up to about 23×19,900 numbers per model). Compression may be appropriate; deterministic gzip requires a fixed mtime, and the manifest should bind exact delivered bytes plus the canonical uncompressed JSON where used. Do not silently discard full projection state or inactive candidate genes to reduce size.
- **Rank sparse-coefficient trap:** coefficients may be sparse, but every selected candidate participates in each sample's ranks. Re-ranking only active genes changes the model.
- **Units and score meaning:** predictions retain original GDSC2 AUC units/direction. They are not calibrated CTRPv2 AUC, clinical sensitivity probabilities, or per-cell drug-response measurements.
- **Selection vs inference:** portable replay demonstrates that the existing frozen model can be evaluated from its state. It does not independently reproduce model selection, nested-CV uncertainty or historical results; those require the documented development inputs, folds, source and original environment.
- **Floating point:** sklearn/scipy-free inference can be independent while still requiring NumPy float64 linear algebra. Exact platform-independent bitwise identity should not be promised; report the measured maximum discrepancy and tested environment separately from original fitting versions.
- **No preservation overclaim:** a current payload hash match and synthetic prediction agreement establish this export's integrity, not the absence of all historical unsaved attempts or a new prospective blind validation.
- **Scientific lock:** original scientific source, checkpoints, joblib artifacts, configurations, manifests and results remain unchanged. Author declarations, archive publication and license decisions remain pending for the agreed final stage.

## Source pointers

- `src/ordinal_repeated_models.py:13`: final arm declarations and gene budgets.
- `src/final_model_v1.py:fit_full_development_model`: full-development grid/refit contract.
- `src/ordinal_diagnostic_models.py:144`: candidate/feature fitted state.
- `src/ordinal_diagnostic_models.py:186`: final representation transform.
- `src/ordinal_diagnostic_models.py:230`: DiagnosticRegressor fitted state and target/tissue policies.
- `src/ordinal_diagnostic_models.py:300`: same-policy tissue-null prediction.
- `src/modeling/nested_cv.py:121`: StandardScaler/FWL/EN fitted arithmetic.
- `src/modeling/nested_cv.py:184`: inherited molecular prediction arithmetic.
- `src/modeling/features.py:_encode_pair_block`: exact three-value raw-TPM comparison.
- `src/run_final_model_v1.py:358`: original complete input assembly.
- `src/run_final_model_v1.py:392`: historical reload sentinels.
- `src/run_final_model_v1.py:463`: original serialized artifact dictionary.
- `config/external_validation_v1_1.json:external_gateway`: complete schema, molecular abstention, tissue and zero-state policies.

