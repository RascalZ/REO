# REO v0.1.0 frozen model assets

This package supplies 64 existing portable JSON.gz models (eight drugs by eight
methods), three parameter/candidate tables, the complete ordered gene schema,
the tissue vocabulary and a generated full-gene input. Model bytes are copied
unchanged from the prepared local export; their SHA-256 values are preserved.
No model was fitted, selected, rescored or changed while building this release.

Use `scripts/predict_repro_v1.py` from the matching REO repository version with
an explicit model and a response-free input JSON. Required input keys are
`schema_version`, `sample_ids`, `biological_group_ids`, `tissues`, `gene_ids`,
and `raw_tpm`. Gene IDs must match the complete ordered model schema. The
packaged synthetic input contains generated numbers, not expression
observations or a performance benchmark. No real response labels or individual
external predictions are included.

Target means/scales, per-tissue target means, feature means/scales, per-tissue
feature means and coefficients are existing fitted prediction parameters.
They must not be interpreted as raw individual response arrays. Zero molecular
solutions remain valid recorded solutions; this release changes no branch or
abstention rule.

`ASSET_MANIFEST.json` inventories all other files in this archive. The outer
release file manifest also hashes this archive's manifest. See
MODEL_ASSETS_LICENSE.md and THIRD_PARTY_ATTRIBUTION.md for the limited author
rights grant and upstream attribution. Obtain original upstream inputs
independently. This asset package supports fixed-model numerical replay; it
does not by itself reproduce the historical evaluation or supply training data.
