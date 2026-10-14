# LOTO execution amendment v1.3

Date: 2026-10-14. This amendment repairs a checkpoint/OOF row identity omission found in the first v1.2 task. The frozen v1.2 preflight had passed. Its first task (drug 1931, held-out Biliary Tract, tissue-null arm) executed and wrote a completed checkpoint with one OOF row and one tissue-null fit event, but the controller rejected it during checkpoint validation; the formal v1.2 run therefore stopped at 0 accepted tasks. No molecular model task was run.

## Observed contract failure and correction

The checkpoint contains the required `payload_sha256`; checksum verification succeeds. The actual failure is that the OOF row omitted `inner_seed`, although the frozen task includes it. The shared checkpoint validator checks every expected task field against every OOF row, so it rejected the row with a missing-key error. The task, split, fit and prediction were otherwise completed and remain preserved as a failed/unaccepted v1.2 artifact.

Execution v1.3 adds the locked `inner_seed` value to each per-sample OOF record and to the aggregate OOF CSV field list; the completed-run validator independently checks every aggregate row against the locked seed for its held-out tissue. The corrected row will carry every task identity field. A no-model-fit unit regression test uses a minimal synthetic checkpoint fixture to verify that omitting the seed is rejected, the correct seed passes, and an incorrect seed is rejected. Separately, the preserved v1.2 checkpoint will be checked in memory with the seed added to a copy and a recomputed payload digest; its on-disk content remains unchanged.

## Scope and status

This is an output-provenance/contract repair only. It does not alter the 23 tissues, eight drugs, eight arms, response inputs, sample universe, seeds, splits, fitting, feature selection, tuning grid, estimands, metrics, NA/missingness rules, or task budget. v1.3 receives its own config, code snapshot, immutable lock, and run directory. The v1.2 run and checkpoint are retained locally as a failed implementation attempt and will not be resumed or imported into v1.3.

The expected v1.3 contract remains 1,472 tasks, 40,112 fit events, and 39,736 outer predictions. CTRP remains sealed; no score or external outcome informed this correction.
