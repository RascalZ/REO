# Minimal reproduction assets and response-free inference

This package separates numerical inference, reconstruction of recorded splits,
and full historical source-data reproduction. It does not retrain or score a
model, reopen CTRPv2 labels, or establish independent scientific replication.
The existing manuscript and all frozen analyses remain unchanged.

## Current release access — 27 October 2026

The public [GitHub Release v0.1.0](https://github.com/RascalZ/REO/releases/tag/v0.1.0)
is published and immutable. Its tag resolves to
`643e4c9d0bd8b598cb9d9e98a4003bdd6004ec45`. The same 11 release files are
published on [Zenodo](https://zenodo.org/records/23681422), with software DOI
[10.5281/zenodo.23681422](https://doi.org/10.5281/zenodo.23681422).
This is a software/research-asset identifier, not a paper DOI. All 11 GitHub
attachments were downloaded without credentials or cookies and checked for
exact size/SHA-256 equality with the published local files; the
[release audit](../docs/public_release_v1/GITHUB_RELEASE_AUDIT.json) also records
the anonymous tag and immutable-state checks.

| Published file / direct download | Bytes | SHA-256 |
|---|---:|---|
| [reo-code-v0.1.0.zip](https://github.com/RascalZ/REO/releases/download/v0.1.0/reo-code-v0.1.0.zip) | 3,645,168 | `693e46cc8a62f2608715021a3322a9c49a0b94eee9ce038d9c9bd7a2ecef1d40` |
| [reo-model-assets-v0.1.0.zip](https://github.com/RascalZ/REO/releases/download/v0.1.0/reo-model-assets-v0.1.0.zip) | 72,937,103 | `85ee348b5632ca20cda6218efb874e034c482deab84484b43eecf4cb5dfe1334` |
| [reo-metadata-assets-v0.1.0.zip](https://github.com/RascalZ/REO/releases/download/v0.1.0/reo-metadata-assets-v0.1.0.zip) | 1,077,796 | `df5562a5cc9a07ac3ac70563e91b125690e5dc296b340941c5e1c243a1b3a49f` |

The [public summary](../docs/public_release_v1/PAYLOAD_SUMMARY.json) and
[per-file catalog](../docs/public_release_v1/PAYLOAD_ASSET_CATALOG.json) describe
the published asset bytes. All 64 model payloads and four fold CSVs preserve
their original bytes. Four new identity/membership projections preserve every
row and retained field string while omitting source experiment/curve/response
pointers, QC records and dose fields. They are not complete inputs for the old
external runner's extraction or rescoring. Original upstream data and wider
historical curation inputs must be obtained independently.

Newly authored code is MIT. The packages include a scoped CC BY 4.0 grant for
author-owned model/derived assets, including project assignments and
author-created identity decisions; third-party data, identifiers and
annotations are not relicensed. Read the packaged license and attribution
notices. Public access does not replace upstream access or redistribution terms.

## Immediately runnable public example

Use Python 3.12 with the two pinned numerical dependencies in
`requirements-repro-inference.txt`. From the repository root:

```sh
python -m pip install -r requirements-repro-inference.txt
python scripts/predict_repro_v1.py --model examples/repro_inference_v1/toy_rank_model.json --input examples/repro_inference_v1/toy_input.json --output toy_predictions_new.json
```

The output path must be new. All five gene names, tissues, TPM values and model
parameters in this example were invented by hand, with no fitting. It is a
software demonstration, not a paper model or a drug-efficacy result. See the
[worked arithmetic](../examples/repro_inference_v1/README.md).

## Preserved complete local technical exports

`local_assets/` remains Git-ignored as the preserved complete technical export.
The published release above uses an audited subset/projection of this export,
not a wholesale copy of its source-pointer metadata. The local export is not
an additional public download. Historical summaries below retain their
pre-release scope; the current release catalog above identifies the new packages.

| Asset | Local location | What is preserved |
|---|---|---|
| Development identities and drug membership | `local_assets/metadata/` | 626 groups, 4,967 drug–group rows; exact canonical, CCLE, group and tissue fields |
| External identity metadata, without labels | `local_assets/metadata/` | 175 groups, 1,185 drug–group rows; existing identities/QC/source pointers only |
| Four recorded split manifests | `local_assets/metadata/` | Original repeated development, final inner, LOTO outer and LOTO drug-inner bytes; no new splits |
| Complete 64 fixed numeric models | `local_assets/models/*.json.gz` | Full ordered input schema, candidates, all zero/nonzero coefficients, scalers, tissue projection and target scaling |
| Candidate and coefficient tables | `local_assets/tables/` | Every candidate and encoded coefficient, including zeros; unpenalized tissue targets are separate |
| Generated full-gene input | `local_assets/synthetic_inputs/` | Invented numerical values, exact 19,836-gene schema and observed/unseen tissues; no source expression row |

The original eight-compound identity and dose-window table already exists as
[Table S8](../tables/manuscript_v1/TableS8_compound_identity_mapping.csv);
the current manuscript-package copy is
[Additional file 8](../manuscript_v1/bmc_additional_files_scientific_revision_v1/Additional_file_8.csv).
Do not infer a new chemical match or harmonized dose window from inference.

## Frozen-model inference using the published release

Download the model archive above or from Zenodo and check its SHA-256 before
extraction. Use the versioned code ZIP or the v0.1.0 checkout. This
standard-library check must succeed:

```sh
python -c "import hashlib,pathlib; p=pathlib.Path('reo-model-assets-v0.1.0.zip'); assert hashlib.sha256(p.read_bytes()).hexdigest() == '85ee348b5632ca20cda6218efb874e034c482deab84484b43eecf4cb5dfe1334'; print('model archive SHA-256 verified')"
python -m zipfile -e reo-model-assets-v0.1.0.zip release_assets/model
python scripts/predict_repro_v1.py --model release_assets/model/models/d1930_response_association_pair_full.json.gz --input release_assets/model/synthetic_inputs/full_gene_order.synthetic.json --output released_model_synthetic_predictions_new.json
```

This uses a real frozen molecular model with explicitly generated input.
No metadata CSV or response file is opened, and no scientific performance
measure is calculated. The output must be new. The public download audit
establishes file access and integrity; it is separate from model inference or
scientific performance validation.

## Preserved local inference route

After the one-time local technical export has been prepared, this command uses
an existing molecular model and the generated input, with no response file:

```sh
python scripts/predict_repro_v1.py --model reproducibility_v1/local_assets/models/d1930_response_association_pair_full.json.gz --input reproducibility_v1/local_assets/synthetic_inputs/full_gene_order.synthetic.json --output local_synthetic_predictions_new.json
```

No native joblib or sklearn object is loaded by this inference CLI. Exact gene
column order and complete float64-compatible raw TPM are required. It does not
repeat detection filtering, gene selection, scaling fitting or hyperparameter
selection on new samples. Output rows retain abstentions and their reasons.

The numerical module expects identities and tissues that have already been
canonicalized. It does **not** authenticate sample ancestry, map chemical IDs,
apply cancer eligibility, prove external-group independence or replace the
frozen external gateway. Known sealed paths are rejected before opening, but
that convenience check is not an OS sandbox or a detector of renamed content.

## Parameter and zero-solution interpretation

Let `x` denote an encoded feature row, `m,s` its recorded training scaler,
`z=(x−m)/s`, `d` the one-hot vector over this drug model's observed tissues,
`M` the training tissue means of standardized features, `a` the training
standardized-target tissue means, and `β` the complete molecular vector.

`prediction = target_mean + target_sd × ((z − dM)β + da)`.

Feature construction comes first: log2(TPM+1), within-selected-set average
midrank, or every oriented raw-TPM pair sign −1/0/+1. Full pair expansion and
ranks preserve the same selected-set order but have different linear models
and penalties. Every candidate is needed to compute a rank, even when its own
coefficient is zero. Complete pair columns, including constant/high-tie ones,
and the full tissue matrix are retained in the local model; no intercept
folding or active-only approximation changes the numerical sequence.

For an unseen tissue, `d=0`. For `β=0`, a molecular arm predicts its training
tissue mean, or the overall training mean if unseen. This is a valid selected
model, not a fitting failure, a proof of no biological signal, or an equivalence
result. All 29/56 molecular zero solutions are retained, separate from the eight
explicit tissue references. Zero-coefficient molecular arms still abstain on
any invalid full-TPM value, including an unselected gene; only the explicit
tissue reference can predict without valid molecular values when identity and
tissue are available. No prediction is clipped to [0,1].

`coefficient_original_response_per_encoded_feature` is target-SD × β / feature
scale. It is conditional on the tissue projection and the other model terms;
it is not a raw-TPM effect or a causal resistance mechanism. The complete
preprocessing state is needed for inference, rather than a feature list alone.

## What remains outside this package

- Actual training-response and large expression matrices are inventoried with
  versions and previously recorded hashes, not recopied or downloaded here.
- Original source-data curation and controller lock verification retain wider
  metadata/legacy dependencies. The final numeric worker does not directly use
  those legacy tables, but that does not remove historical audit dependencies.
- The numeric export covers the 64 final models. Historical outer-fold model
  checkpoints, original OOF predictions and bootstrap objects are not newly
  exported or replayed here; original scientific records remain in place.
- Software and model/metadata publication is complete; author declarations and
  final manuscript approval remain separate closing steps.
- Synthetic numerical agreement is not predictive scoring, an efficacy result,
  a new validation cohort, or a reproduction of the manuscript's correlations.

See [input/identity audit](../docs/reproducibility_assets_v1/INPUTS_AND_IDENTITIES.md),
[fold definitions](../docs/reproducibility_assets_v1/FOLD_DEFINITIONS.md), and
[model export review](../docs/reproducibility_assets_v1/MODEL_EXPORT_REVIEW.md).
Do not run historical curation, validators that open response-bearing sources,
training entry points or external scoring to test this inference example.
