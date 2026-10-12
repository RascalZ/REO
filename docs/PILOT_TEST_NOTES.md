# Development pilot integrity acceptance

This acceptance suite uses synthetic expression and response values. It does not read CTRP response files and does not assess model performance on the external cohort. Its purpose is to detect the full-data feature-selection and tuning errors identified in the original analysis.

## Required invariants

1. **REO ties and direction.** A pair is directed by stable gene identity, never by the order of the input columns or by the response labels. Equal expression produces the third value `0`; greater and smaller expression produce opposite signs. This test includes exact zero ties and nonzero ties.
2. **Train-only feature eligibility and ranking.** Holding the training rows fixed while radically changing the held-out expression distribution must leave the selected genes, selected pairs and transformed training matrix unchanged. A fitted builder must not refit when transforming held-out rows.
3. **Train-only supervised selection.** If any selector uses response labels, changes restricted to held-out labels must not change the selected features or hyperparameters. The present pilot uses label-independent training-fold expression eligibility/variance selection; future supervised selectors must satisfy the same boundary.
4. **Group separation.** Models sharing one biological group remain together in every outer and inner split. This remains true when one group has several models or samples. Random row splitting is not accepted.
5. **Inner refitting.** Every inner training partition must fit its own feature builder. One outer-training or full-data builder reused in all inner partitions is a failure even if no outer-test rows were seen.
6. **Outer-label isolation.** Perturbing only the test labels of one fixed outer fold must leave that fold's hyperparameter choice and predictions unchanged. Its evaluation scores may change. This separates evaluation from selection.
7. **Development-only data access.** The runner may load GDSC responses for the frozen development rows, the frozen development manifest and baseline-expression matrices/identities. It must not open historical or current CTRP response tables during the pilot. The runner's development preflight checks the lock and actually used development files without opening the external-response descriptors.
8. **Frozen assignment integrity.** Altering a development assignment or the locked cohort description must fail before model fitting. Selected GDSC rows must match the frozen drug ID, curve ID and exact concentration interval. Unselected response rows cannot enter the extracted pilot response table.
9. **Failed-fit state and audit columns.** An invalid refit must clear an earlier fold's fitted feature state. Changing audit-only sample IDs cannot affect hyperparameter selection or predictions.

## Why these are substantive tests

The tests change information that would be available to a leaking procedure while keeping the legitimate training inputs unchanged. They then compare training artifacts or predictions, rather than testing that an implementation returns its own advertised constants. Group and refit tests record actual training memberships, so a passing score alone cannot conceal leakage.

## Scope and limitations

Passing this suite validates the tested synthetic cases and the inspected pilot data-loading paths. It does not prove the absence of every possible leakage route, validate clinical utility, or justify selecting drugs based on pilot performance. The main panel and CTRP cohort remain locked under the phase-1 protocol. Any later selector, imputation step, feature-ranking rule or model must be covered by an equivalent fold-local test.

## Execution

From the repository root, run the repository Python environment:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -p test_pilot_integrity.py -v
```

## Verified execution (2026-10-10)

`24 tests` passed in the repository environment; the latest complete execution took `2.208 seconds` of unittest time (process startup is excluded). These tests fit only synthetic models and use temporary mock workbooks/files. Their fits are not counted as fits in the real two-drug pilot.

| Test family | Tests | What was checked |
|---|---:|---|
| Fold feature construction | 9 | Exact three-value REO, canonical direction, column-order invariance, training-only gene/pair filters, label independence, invalid-input handling, failed-refit state |
| Nested-CV boundaries | 7 | Group isolation/coverage/reproducibility, actual inner-fit memberships, outer-label and outer-expression perturbations, audit-ID independence, unseen-tissue handling |
| External-response access guard | 1 | Six prohibited mock read attempts rejected, including byte-path `os.open`; three development mock inputs readable |
| Frozen response-row extraction | 4 | Only selected rows extracted; wrong drug/curve/concentration, duplicate/missing rows and invalid AUC rejected |
| Frozen input verification | 2 | Changed assignment/lock refused; used development hashes verified without reading the deliberately absent external-response file |
| Failure evidence and resource shutdown | 1 | Simulated runtime failure and interruption recorded as failed; resources saved, monitor thread stopped and partial checkpoint preserved |

The unseen-tissue baseline returns the current training target mean. The expression model assigns unseen tissues zero tissue effect while retaining its fitted gene component. The tests verify that changing one unseen tissue code to another changes neither prediction nor fitted category state.

## Source-path review

The inspected runner loads the frozen development manifest and panel, GDSC2 source responses, protein-coding annotation and baseline TPM. Its metadata-based global splits are intersected with drug-specific observations; inner indices are relative to the local outer-training rows. Target scaling, feature filters, REO construction, feature scaling, tissue encoding and tissue adjustment are refitted on each training partition. Outer-test labels are used for metrics after tuning and prediction.

The current guard rejects the raw CTRP directory, historical CTRP response table and three known processed response-index files. Behavioral acceptance checks those paths. This is protection against accidental access to the current files; it is not a claim that every possible future external-response file name is automatically forbidden. Any new loader requires renewed review.

The mean predictor has no sample-specific ranking within an outer fold. The runner records pooled mean-baseline Spearman as unavailable rather than interpreting differences among fold means as ranking performance.

## Problems identified and repaired before the real pilot

- The initial runner recorded the cohort lock but did not independently verify the development assignment and expression-input hashes. It now verifies the lock and eight actually used development files, so an edited expression-row assignment cannot silently pass the curve-ID checks.
- Validation failures in the workbook response reader initially left the input workbook/iterator open on Windows. The failure fixtures detected file-handle retention through `WinError 32`. The reader now closes both in a `finally` block; all rejection fixtures and their cleanup pass.

No unresolved leakage or file-access failure was found in the tested cases. Run failures, convergence problems, performance estimates and resource measurements from the actual pilot are reported separately in its run artifacts and master research plan.

## Independent audit of the actual pilot artifacts

Run:

```powershell
.\.venv\Scripts\python.exe scripts/validate_pilot_run.py --run-dir runs/gdsc2-mira1-plx4720-development-nested-en-reo80 --source-revision 629b8f5
```

The validator passed all recorded boundaries: `462` actual model fits, including `342` gene-feature fits; each actual training sample-ID set matched its intended inner-training partition or final outer-training refit. Current validation/test identities and known biological groups were absent from those fit sets. Every drug/model/identity had exactly one outer OOF prediction: `626` for MIRA-1 and `623` for PLX-4720, across five models, giving `6,245` OOF records.

The audit independently reconstructed target means/scales from the corresponding training observations in run artifacts, checked fixed-grid selections against logged inner scores, and reconstructed pooled outer RMSE, MAE, R² and Spearman. All four recorded configuration/source hashes matched immutable execution revision `629b8f5`; the later workbook-failure/resource-cleanup changes are explicitly distinguishable from the executed source. No original GDSC response workbook, CTRP response table or expression matrix was reopened by this validator.

The result is `pilot_validation.json` inside the run directory. Its paths are repository-relative and its counts are aggregate. This audit verifies the recorded fit memberships; the synthetic acceptance suite separately inspects actual builder inputs and perturbs held-out information. Inner OOF predictions were not exported, so the validator checks hyperparameter choice against recorded inner scores rather than independently reconstructing those inner scores.

## Cancer-cohort and expanded-model acceptance (2026-10-12)

The v2 suite adds `11` synthetic tests to the original `24`, giving `35 tests`, all passing in `3.071 seconds` of unittest time. The added cases independently reconstruct each training fit's EN `alpha_max` from standardized, within-tissue residuals; inspect actual RF selector inputs; change outer-test labels and expression; verify stable within-tissue Pearson pair ranks; and check that the forest receives original molecular values plus training-fitted one-hot tissue. Tissue and audit-ID columns are excluded from gene/pair construction. A zero molecular EN `alpha_max` preserves training tissue effects. Zero within-tissue response association gives RF zero molecular pairs with an explicit tissue-only status.

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_pilot*_integrity.py" -v
```

The independent metadata validator is:

```powershell
.\.venv\Scripts\python.exe scripts/validate_cancer_cohort.py --lock
```

It locked cohort v2 after verifying metadata-only filtering of all three v1 sample lists, unchanged drug rank/IDs/structures/dose windows and biological identities, every retained/excluded source row, frozen-v1 metadata hashes, and separation from **all known GDSC identity-candidate groups**, including ambiguous candidates. The resulting main cohort contains `626` development and `175` primary-external identities/groups, with `4,967` and `1,185` drug/model observations. Every drug retains at least `200` development and `80` external groups; the smallest external count is `82`.

HCC1588 cancer eligibility is corroborated by the exact `ACH-001078` / `CVCL_A351` / `SIDM01820` / `KCLB 71588` chain and cached official cancer classification. The other unverified disease annotation remains excluded. Only capitalization-equivalent tissue changes are accepted, and every original tissue value remains recorded.

This validation opens metadata/manifests/code only. It rejects attempts to open the raw CTRP directory, known response indices and expression arrays. Inherited expression-input hashes are compared with the locked v1 descriptors; their large-file bytes are deliberately deferred to the development runner's input verification. A successful recheck preserves the exact v2 lock and validation-report bytes; it does not rewrite the historical snapshot. The v2 amendment was made after the development pilot and before external-response access, as disclosed in the research plan.
