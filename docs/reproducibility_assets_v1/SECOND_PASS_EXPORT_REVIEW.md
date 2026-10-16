# Second-pass review of prediction-only model export

Review date: 2026-10-16. Review mode: read-only source/contract inspection. This review did not execute an export, train or score a model, load a serialized model/portable parameter file, or read real expression rows, response labels, CTRP response indexes, original sentinel predictions, individual external predictions or bootstrap files. It is a second implementation review, not an independent reproduction of the scientific study or a security certification.

Files reviewed:

| Source | SHA-256 at review |
|---|---|
| `scripts/export_repro_models_v1.py` | `d94b072e5b81363f57c1e9ff18f93ca107bd03aa3958fdc544245f06f59602a6` |
| `src/repro_inference_v1.py` | `2815e700ea073deb3a01e0763e6959790fcb1ea356f628bb285a5cf28cb8509b` |

The corresponding original implementation was checked in `src/ordinal_diagnostic_models.py` and `src/modeling/nested_cv.py`; the two original locks were read only for their metadata fields. Both lock execution-commit fields match the export's recorded scientific commit, `ca0f837ac10665ab8d4d6ab97691e58e5bed4577`.

## Review outcome

No blocking defect was identified in the inspected export/inference path for the existing fixed 64-model collection. This outcome is conditional on the original trusted locks and artifact hashes and on successful completion of the export's own checks. At review time `model_export_public_summary.json` had not yet been written, so this review does not attest that all 64 exports or 384 numerical comparisons have completed. The root process must accept the completed manifest, count/hash checks and no-label example separately.

## 1. No-fit, no-score and no-CTRP data path

The export reads the two existing locks, frozen public model manifest, worker context, final-run metadata, bound original scientific source files and the 64 hash-bound native artifacts. It does not call the training runner, response extraction, source-curation functions, cohort validator, original external gateway or external scoring code. The only expression values created by the export are six generated TPM rows under a fixed PCG64 seed. The generated row/group IDs and an unseen tissue label are synthetic.

Before native artifact loading, the script disables `fit` on the diagnostic regressor/builder, inherited train-only regressor, ElasticNet, StandardScaler and OneHotEncoder, and disables the diagnostic regressor's `score`. The inspected code contains no other training or predictive-performance call. It does call the already fitted native model's `predict` and the new numerical implementation's prediction function to compare arithmetic on generated inputs. This is an implementation-agreement check, not an accuracy, AUC, Spearman, RMSE, model-selection or external-validation score.

The export audit hook allows reads of an explicit data-asset list and code modules, restricts writes to the requested local export namespace/public-summary path, and denies network/subprocess operations after hook installation. This is a useful guard for the inspected script. It is not a general sandbox: code-extension reads remain allowed, imports occur before hook installation, and a trusted pickle/joblib load can execute code. The design is therefore appropriate only for the project's known hash-bound artifacts; it must not be described as safe loading of arbitrary third-party joblib files.

The standalone predictor reads explicit JSON/JSON.gz files, rejects known sealed CTRP locations/response basenames before opening, uses an exact response-free sample-key set and does no implicit source lookup. Its path check is documented as a convenience guard, not a scanner that can detect labels renamed into an arbitrary JSON file. The generated example and its explicit schema provide the substantive no-label guarantee for this workflow.

## 2. Binding to existing artifacts and keeping originals

The export requires the frozen manifest to contain 64 rows and the second lock to record 64 models. It verifies the current first-lock hash against the second lock and the current manifest hash against the second lock before reaching artifact loading. For **each** native artifact, SHA-256 must match both the row's expected artifact hash and the second lock's artifact dictionary before `joblib.load`. It then verifies task identity, first-lock binding and tissue coding in the loaded artifact. Worker-context bytes and final metadata are checked against their existing records, and bound original scientific-source files must still match the first lock.

This is a chain of consistency checks anchored in the original trusted locks; the new export script does not independently establish the authenticity of those historical locks or repeat their historical scientific validation. The 64-row task coverage is inherited from the exact frozen manifest whose hash is checked, rather than re-proven from a newly generated factorial task list.

Every original artifact is rehashed after its prediction check. The new script has no joblib dump/write-back operation and opens original artifacts only for reading; output creation uses exclusive modes and refuses an existing export. Native prediction changes in-memory audit fields such as prediction-call counters, elapsed prediction times and synthetic prediction sample IDs. Those changes are not serialized back to the original artifacts. A completed manifest plus unchanged original hashes are the appropriate evidence; do not claim that the in-memory model objects are never mutated.

## 3. Prediction state and numerical order

The portable bundle contains the full ordered 19,836-gene input schema, candidate positions/IDs, feature order, full oriented pair list where applicable, fitted scaler means/scales, feature-by-tissue means, target-by-tissue means, training tissue codes, target mean/SD and the complete coefficient vector including exact zeros. No fitting is required for prediction. Candidate scores and signed associations are also exported to a separate local candidate table, while the coefficient table retains inactive coefficients rather than selecting only apparent successes.

The new calculation follows the original fitted prediction order:

1. Expression: `log2(TPM + 1)` on the frozen candidate positions; rank: average midranks within the same candidate set using exact raw TPM ties; pair: all oriented candidate comparisons coded −1/0/+1.
2. Subtract the fitted feature means, then divide by fitted scales.
3. Subtract the training tissue projection from the already encoded/standardized features; multiply by the frozen coefficient vector; add the fitted standardized target tissue means.
4. Multiply by the fitted target SD and add the target mean. Do not clip the response.

For a legal unseen tissue, the one-hot row is zero, so the model uses the global training target mean plus any retained molecular component. This agrees with the original encoder policy. The whole ordered gene schema and row validity are checked before prediction. Missing/invalid identity and tissue cause explicit row abstention; structural schema faults are hard failures. The predictor does not infer canonical identity, external cohort eligibility, compound identity, genealogy or tissue normalization itself; these remain upstream responsibilities.

The export JSON round-trip is checked against the native model on six generated cases per model with the stated absolute/relative tolerance of 1e−12 and one computational thread. Those cases include ties, all-zero expression, a very small numerical range and an unseen tissue. Such finite probes can reveal many transcription/order errors but are not an exhaustive proof over every possible input. Numerical tolerance agreement must not be relabeled byte-identical prediction equality.

The bundle is complete **for frozen prediction arithmetic**. It is not a replacement for the full training audit: original feature eligibility/filter settings, the entire tuning grid, solver iterations/dual gaps/alpha maxima, original warnings, fit membership and checkpoints remain in the frozen configuration/logs. The original-response-per-feature coefficient exported in the CSV is a conditional coefficient conversion; the CSV alone omits the centering/tissue components needed to reconstruct predictions. The full bundle supplies those components. Coefficients are predictive model parameters, not causal effects.

## 4. Molecular zero versus explicit tissue reference

The exporter separately records `declared_molecular_arm`, exact nonzero coefficient count and `zero_molecular_solution`; the explicit tissue reference gets `null` for the last field. Thus the eight tissue references do not inflate the 29/56 molecular zero count.

For declared molecular arms, the predictor still requires the complete finite nonnegative TPM row even when all molecular coefficients are exactly zero. Their tissue-only numerical predictions remain legitimate frozen model outputs, but the input policy is not relaxed after seeing the zero solution. The explicit tissue-null arm may predict without molecular values while preserving the same identity/tissue/schema requirements. This distinction matches the intended missingness/abstention policy. It does not turn a zero solution into a failed fit or prove absence of underlying molecular biology.

## 5. Public summary and publication limits

The inspected public-summary writer publishes the same **summary object** as the local export manifest: counts, paths/hashes, package versions, per-drug/arm feature/nonzero counts, fallback reasons, synthetic-case numerical agreement, and read-access/guard metadata. It does not include coefficient/scaler/tissue-mean vectors, real cell/model/group identities, real expression vectors, training/test labels, individual original predictions or synthetic raw prediction values. Those parameter tables and portable bundles remain in the ignored local-assets directory pending the release/access decision.

The summary does include real **drug assay IDs** and relative artifact/source paths. Consequently, the accurate description is “no individual cell/group IDs or labels and no fitted parameter vectors,” not “no real IDs of any kind.” Package versions, native-source hashes and exported-file hashes are provenance, not predictive evidence. Publishing only this summary does not make the local portable models or complete identity manifests accessible to an independent researcher.

This agrees with `INPUTS_AND_IDENTITIES.md`: the new interface separates generated-input inference, portable-model replay and full historical data-curation reproduction. Final acceptance should explicitly report which local assets are prepared and which remain unpublished, retain a clear synthetic/example label, and keep author declarations, hosting/archiving and upstream redistribution decisions at the deferred final publication gate.
