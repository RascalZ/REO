# Development pilot integrity acceptance

This acceptance suite uses synthetic expression and response values. It does not read CTRP response files and does not assess model performance on the external cohort. Its purpose is to detect the full-data feature-selection and tuning errors identified in the original analysis.

## Required invariants

1. **REO ties and direction.** A pair is directed by stable gene identity, never by the order of the input columns or by the response labels. Equal expression produces the third value `0`; greater and smaller expression produce opposite signs. This test includes exact zero ties and nonzero ties.
2. **Train-only feature eligibility and ranking.** Holding the training rows fixed while radically changing the held-out expression distribution must leave the selected genes, selected pairs and transformed training matrix unchanged. A fitted builder must not refit when transforming held-out rows.
3. **Train-only supervised selection.** If any selector uses response labels, changes restricted to held-out labels must not change the selected features or hyperparameters. The present pilot uses label-independent training-fold expression eligibility/variance selection; future supervised selectors must satisfy the same boundary.
4. **Group separation.** Models sharing one biological group remain together in every outer and inner split. This remains true when one group has several models or samples. Random row splitting is not accepted.
5. **Inner refitting.** Every inner training partition must fit its own feature builder. One outer-training or full-data builder reused in all inner partitions is a failure even if no outer-test rows were seen.
6. **Outer-label isolation.** Perturbing only the test labels of one fixed outer fold must leave that fold's hyperparameter choice and predictions unchanged. Its evaluation scores may change. This separates evaluation from selection.
7. **Development-only data access.** The runner may load GDSC responses for the frozen development rows, the frozen development manifest and baseline-expression matrices/identities. It must not open historical or current CTRP response tables during the pilot. Merely validating a CTRP file's existing provenance hash is distinct from loading its response values.

## Why these are substantive tests

The tests change information that would be available to a leaking procedure while keeping the legitimate training inputs unchanged. They then compare training artifacts or predictions, rather than testing that an implementation returns its own advertised constants. Group and refit tests record actual training memberships, so a passing score alone cannot conceal leakage.

## Scope and limitations

Passing this suite validates the tested synthetic cases and the inspected pilot data-loading paths. It does not prove the absence of every possible leakage route, validate clinical utility, or justify selecting drugs based on pilot performance. The main panel and CTRP cohort remain locked under the phase-1 protocol. Any later selector, imputation step, feature-ranking rule or model must be covered by an equivalent fold-local test.

## Execution

From the repository root, run the repository Python environment:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -p test_pilot_integrity.py -v
```

The test execution result and any remaining issues will be recorded here after the pilot interfaces are finalized.
