# Minimal reproduction technical audit — 2026-10-16

This is the M04 local-asset step of master plan v3.36, completed into v3.37.
No scientific result, model fit, split assignment, external score or manuscript
was changed. Complete derived mappings/models remain local pending access and
rights decisions; this directory publishes counts, hashes and technical checks.

## Entry points

- [Runnable package guide](../../reproducibility_v1/README.md) and
  [numeric formats](../../reproducibility_v1/FORMAT.md).
- `INPUTS_AND_IDENTITIES.md` / `input_asset_inventory.csv`: 69 source/dependency
  records, with previously recorded hashes; distinguish prediction, final-fit
  numeric inputs, freeze verification and historical source-data reconstruction.
- `FOLD_DEFINITIONS.md` / `fold_sources.csv`: exact recorded memberships,
  global-versus-drug-specific split policies, seeds, rare strata and amendment
  history. No split was regenerated.
- `MODEL_EXPORT_REVIEW.md` / `SECOND_PASS_EXPORT_REVIEW.md`: source-only,
  role-separated checks; neither review deserialized the native payloads.
- `metadata_export_public_summary.json`: local extraction counts/hashes,
  without individual cell/group identities or response values.
- `model_export_public_summary.json`: all 64 models' source/portable hashes,
  dimensions, zero status and generated-input numerical agreement. Real drug
  IDs are retained; individual cell/group IDs and fitted vectors are absent.
- `asset_bundle_manifest.json`: hashes of 87 complete local derived assets;
  listing their hashes does not mean they are publicly downloadable.
- `public_example_isolation_audit.json`: a fresh minimal copy of public files
  worked with no old pipeline, local model package, sklearn/joblib import or
  CTRP files. It reused the installed pinned environment, not a fresh install.
- `validation.json`: integration checks, code hashes, original-file protection,
  fit/score/access boundaries and remaining limitations.

## What was actually prepared

Local ignored assets contain 626 development and 175 external identities;
4,967/1,185 corresponding label-free drug memberships; four byte-identical
recorded split CSVs; 64 numeric model bundles; 25,600 candidate rows; 340,800
encoded coefficient rows including zeros; 1,472 fitted tissue rows; complete
gene/tissue schema; and generated-only inference examples. All original native
models retained their SHA-256 values. The four split CSVs contain 31,926,
1,878, 14,398 and 109,274 rows, respectively.

The native exporter loaded only trusted artifacts after checking both the
public manifest and the existing second-lock hashes. The original first-lock
environment/scientific sources were checked; `.fit` methods and estimator
scoring were disabled during export. No real expression observation, original
reload sentinel, CTRP label/index, original individual prediction or bootstrap
weight was opened. Every model was checked on six generated inputs; all 384
native-versus-portable outputs agreed exactly in the recorded environment.
This is numerical implementation agreement, not a drug-response score or an
independent replication of the study.

The public inference module has no sklearn, joblib, estimator `.fit` or scoring
dependency. The toy fixture is hand-specified, not fitted and not a paper model.
Its seven rows retain four predictions and three abstentions. Sixteen contract
tests cover formulas, ties, full candidate ranks, gene schemas, zero/null
distinctions, missingness, unseen tissues, path checks and overwrite refusal.

## Verification commands

```sh
python scripts/verify_repro_assets_v1.py
python -m unittest discover -s tests -p test_repro_inference_v1.py -v
```

The first uses the standard library and public records. The second requires
the two pinned numerical dependencies. Locally, `python
scripts/verify_repro_assets_v1.py --local` additionally checks exported numeric
bundle hashes, complete gene order, zero counts, table rows and nine metadata
outputs. It never scores an outcome or opens a source CTRP file.

The export scripts intentionally refuse an existing output. Do not rerun them
into the completed directory or invoke historical curation/model validators to
test this package. Source-data access and the eventual public model/mapping
route still need to be closed alongside the deferred author/archive/rights
work. M04's technical preparation is complete; M04's public-access requirement
is not yet fully closed.
