# LOTO execution amendment v1.1

Date: 2026-10-14. This amendment resolves one execution-parameter conflict before any LOTO model is fitted. The original P5e v1 document and lock remain unchanged.

## Conflict and resolution

The frozen prose in `docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md`, section C, states an inner-fold seed of `2026101801 + sorted_tissue_index`. The machine-readable frozen file `config/tissue_sensitivity_loto_v1.json` instead records `inner_seed_base: 2026101401`. Both files were covered by the P5e v1 lock. They cannot both define the executed split.

Execution v1.1 follows the machine-readable configuration value, because it was the explicit executable seed field and preserves the locked task design. For each held-out tissue, sort the frozen 23 tissue names by Python/Unicode code-point ordering, assign indexes 0–22, and use `random_state = 2026101401 + index` for the inner, training-only, tissue-stratified grouped 3-fold split. The resulting per-tissue seed map is stored in `config/loto_execution_v1_1.json` and in the prepared run manifest. The original protocol and its lock are not edited.

No model scores, external data, or outcomes informed this resolution. It changes only the deterministic inner split seed; all drugs, tissues, models, outer leave-one-tissue splits, inner-fold policy, feature selection, tuning grid, missing-value handling, sample universe, task budget, metrics, and reporting rules remain as frozen in P5e v1.

## Locked scope and gate

- Run the fixed eight drugs, all 23 tissues, and all eight model arms under a shared tissue-level outer split.
- The held-out tissue is absent from every inner and outer training set. Inner folds are generated only among the other tissues using the frozen grouped stratified splitter.
- Expected workload: 1,472 tasks, 40,112 fits, and 39,736 outer predictions.
- Preserve constant-prediction Spearman as NA; for the held-out tissue, the same-policy tissue-null score is expected to be constant and therefore undefined.
- Do not omit failed tasks, zero solutions, or undefined statistics. Do not use held-out tissues to tune or choose models.
- CTRP response, expression, indexes, and outcome-derived artifacts remain inaccessible until the separate final-model/external-statistics gate.

The runner, validator, analyzer, source configuration, amendment, and input hashes are frozen in `config/loto_execution_v1_1.lock.json`. The prepared split manifest is generated and audited before model execution; its hash and all task/test memberships are then bound in run metadata.

## No-fit split preflight finding

The complete preflight found the expected scikit-learn stratification warning in all 184 drug × held-out-tissue contexts: after withholding one tissue, the least-populated remaining tissue has one group in 176 contexts and two groups in 8 contexts, fewer than three inner folds. The grouped split still returns three complete partitions; each outer-training sample appears in exactly one inner validation fold and no biological group crosses inner train/validation. Because the protocol fixes the current three-fold grouped tuning policy, these warnings are retained and reported rather than changing the splitter after inspecting any model score. The limitation is that tuning performance for the one/two-group tissue strata is intrinsically unstable and its validation sample can be absent from that fold's inner training subset.
