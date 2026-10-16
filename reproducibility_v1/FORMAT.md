# Numeric model and input formats, version 1.0

The executable schema validator is `src/repro_inference_v1.py::validate_bundle`.
The checked-in toy model is a complete small example. Actual fixed models have
19,836 input genes; the demonstration has five. Native joblib serialization is
not used by the public inference interface.

## Portable model: `reo-portable-model-1.0`

| Field | Type / shape | Meaning |
|---|---|---|
| `task`, `source` | objects | Drug/arm declaration and original artifact/lock/input hashes; no labels or observed predictions |
| `input_gene_ids` | G unique strings | Exact complete raw-TPM input order; G=19,836 for actual models |
| `tissue_code` | canonical label → integer | Global label vocabulary; it is not this drug's list of observed tissues |
| `training_tissue_codes` | T ordered numbers | The fitted encoder's categories, in their original order |
| `model_kind`, `fitted_kind`, `selection_rule` | strings | Declared representation/policy and actual prediction branch; retain fallback distinctions |
| `candidate_gene_indices`, `candidate_gene_ids` | K integers / strings | Positions in the complete input and stable-ID-sorted output order; actual K=200 or 2,000 |
| `feature_names` | p unique strings | `log2TPM:gene`, `sample_rank:gene`, or `REO:left>right` in fitted order |
| `pair_indices` | p×2 integers, or empty | Full-input positions, not indices in a reduced active set; complete upper triangle in candidate order |
| `feature_mean`, `feature_scale` | p float64-compatible numbers | Original training StandardScaler state; population variance, positive scale; constant features have recorded scale 1 |
| `feature_tissue_means` | T×p numbers | Training means of standardized molecular features, before FWL residualization |
| `target_tissue_means` | T numbers | Training means of the standardized target; not original-AUC means |
| `target_mean`, `target_sd` | finite scalars, SD>0 | Training target inverse transformation; population SD (ddof=0), with original constant-response fallback policy |
| `coefficients` | p numbers, zeros included | Complete fitted coefficient vector in standardized-feature/standardized-target units |
| `n_nonzero_molecular_coeffs` | integer | Exact count of entries unequal to floating-point zero; no threshold pruning |
| `declared_molecular_arm` | boolean | Determines whether full valid molecular input is required, even if coefficients are zero |
| `fallback_reason`, `solver_reason`, `chosen_params` | strings / object | Original fitted state and chosen relative-alpha/L1 parameters, without rerunning selection |

The eight explicit tissue references have p=0 and no candidate genes; their
T×0 feature matrix can be represented as T empty rows. A molecular empty-feature
fallback retains `declared_molecular_arm=true`. The exporter currently observes
no final-model fallback, but the interface distinguishes that branch and a
constant-target `mean` branch in its contract tests.

Models are UTF-8 canonical JSON (float64 round-trip decimal values), optionally
gzip compressed with mtime=0. The export manifest records hashes of compressed
bytes and uncompressed canonical JSON separately. Coefficients and complete
tissue matrices are retained; no sparse-only approximation is called complete.

This is complete numeric prediction state, not a substitute for the original
training configuration, eligibility/filter audit, selected hyperparameter path
history, fold membership or fit-event logs. Those are separately inventoried.

## Input: `reo-inference-input-1.0`

Exactly these keys are accepted: `schema_version`, `sample_ids`,
`biological_group_ids`, `tissues`, `gene_ids`, `raw_tpm`. No response key is
accepted. The aligned row lists must be nonempty, and valid sample/group IDs
must be unique. `gene_ids` must equal the complete model order: no automatic
sorting, suffix removal, symbol mapping, normalization or imputation occurs.

Each TPM row contains G scalar values, or is null for an absent expression row.
Values are parsed as float64. Missing, nonnumeric, nonfinite or negative values
cause molecular abstention. Wrong row dimensions or nested entries are hard
failures. Missing identity/tissue abstains for every arm. Exact raw-TPM ties are
preserved; do not pre-round or supply float32 values if reproducing original
rank/pair behavior. This JSON transport does not recover precision previously
lost upstream.

## Output: `reo-inference-output-1.0`

The output retains aligned `sample_ids`, `predictions` (number or null),
`statuses` (`predicted` or `abstained`), `reasons`, `unknown_tissue`, and the
task. `zero_molecular_solution` is a model-level boolean for a declared molecular
arm and null for the explicit tissue reference. `no_response_labels_used=true`
states the interface contract; it is not content authentication of renamed
files or proof of cohort eligibility. A nonfinite computed prediction is a
hard failure, not a value to clip or replace.

## Tables

- `candidate_genes.csv`: every candidate's original column index, stable ID,
  training selection score and signed association; no new selection occurs.
- `all_coefficients.csv.gz`: all 340,800 encoded coefficients, including zeros,
  with feature name/index, scaler mean/scale and unit conversion. No coefficient
  is a causal gene effect. Complete tissue projection remains in model JSON.
- `tissue_effects.csv`: all 1,472 drug-arm-tissue rows, standardized-target
  means and original-scale tissue-reference values; these are derived fitted
  summaries, not individual response observations.

All these actual fitted-parameter tables remain local until the planned release
scope and terms are confirmed. Hash/count summaries are public; that does not
mean the full fitted tables are currently downloadable from GitHub.
