# LOTO execution amendment v1.2

Date: 2026-10-14. This amendment repairs a controller logging-schema defect encountered on the first v1.1 startup. The v1.1 preflight had passed; its first run stopped before launching any task worker or fitting any model. The failed v1.1 run record is retained locally and must not be resumed under changed code.

## Defect and correction

On Windows, the frozen `process_tree_memory()` helper returns `system_commit_bytes` and `system_commit_limit_bytes` in addition to the memory fields already listed in the LOTO controller resource CSV schema. The v1.1 controller attempted to write the complete resource object into a CSV row without declaring those two fields. It raised a `ValueError` on the first resource sample, before worker startup. The v1.2 controller adds both returned fields to the resource CSV schema.

The retained v1.1 run metadata records 0 of 1,472 tasks completed, 0 completed-task wall seconds, no checkpoints, no fit-event file, and no OOF file. No model fit or external response access occurred. v1.2 uses a new run directory and does not modify or resume v1.1 evidence.

## Scope

This is an execution-observability repair only. It does not alter the frozen 23 tissues, eight drugs, eight model arms, group boundaries, split assignments, seeds, feature selection, tuning grid, estimands, missingness/NA handling, task budget, or resource thresholds. The v1.2 config is copied from the frozen v1.1 config, with only the protocol/amendment/run identifiers and this code-change note updated. The expected contract remains 1,472 tasks, 40,112 fits, and 39,736 outer predictions.

The v1.2 runner, validator, analyzer, config, this amendment, and source/input hashes must be locked in a new immutable execution lock before a new no-fit preflight. CTRP remains inaccessible. No score, model result, or external outcome informed this correction.
