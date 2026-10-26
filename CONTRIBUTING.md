# Contributing

This repository records a completed, frozen computational study. Contributions should separate proposed new work from the recorded analysis.

- Describe the problem, affected command and evidence. Use synthetic inputs or aggregate examples; do not attach source response matrices, individual predictions, credentials or unpublished author documents.
- Preserve frozen configurations, completed results, source hashes, negative results and NA/zero-solution behavior. Propose a new analysis in a separate version rather than overwriting the study.
- Document dependency changes and run checks appropriate to the modified component. Prediction-only changes use `tests/test_repro_inference_v1.py`; aggregate-contract changes use the `verify_public_*` scripts.
- Keep data terms separate from code licensing. MIT covers the new reconstruction code/documentation, not upstream datasets, original manuscripts or legacy code. Do not assign authorship or confirm permission on another person's behalf.

See [REPRODUCING.md](REPRODUCING.md) for commands and limits, and [RESEARCH_PLAN.md](RESEARCH_PLAN.md) for the maintained study record.
