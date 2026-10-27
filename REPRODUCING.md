# Reproduction entry points

Use Python 3.12. These commands distinguish software checks from reproducing the scientific analysis. A passing check does not imply a clean dependency installation or end-to-end scientific retraining.

## Current access — 27 October 2026

The repository is public. The immutable [v0.1.0 release](https://github.com/RascalZ/REO/releases/tag/v0.1.0)
and its 11 attachments are published, with the same files archived on
[Zenodo](https://zenodo.org/records/23681422) under software DOI
[10.5281/zenodo.23681422](https://doi.org/10.5281/zenodo.23681422).
This is not a paper DOI. Anonymous GitHub metadata, tag and all attachment
downloads were verified in the [release audit](docs/public_release_v1/GITHUB_RELEASE_AUDIT.json).
The tag resolves to commit `643e4c9d0bd8b598cb9d9e98a4003bdd6004ec45`.
See the [asset guide](reproducibility_v1/README.md) for filenames, SHA-256 values
and the response-free frozen-model command.
The immediately runnable toy below is already tracked in this public repository.

## Dependency layers

| File | Purpose |
|---|---|
| `requirements-repro-inference.txt` | NumPy/SciPy prediction-only implementation and handwritten-fixture tests |
| `requirements.txt` / `pyproject.toml` | Frozen original curation layer, not the full modeling environment |
| `requirements-modeling.txt` | Modeling and scientific contract tests |
| `requirements-figures.txt` | Figure generation |
| `requirements-manuscript-qa.txt` | Document and layout helpers |

Install only the relevant layer in an isolated environment. For the toy:

```sh
python -m pip install -r requirements-repro-inference.txt
python scripts/predict_repro_v1.py --model examples/repro_inference_v1/toy_rank_model.json --input examples/repro_inference_v1/toy_input.json --output toy_predictions_new.json
python -m unittest discover -s tests -p "test_repro_inference_v1.py" -v
```

The seven invented rows include predicted and abstained cases. They test the numerical input contract, retained zero solutions and unknown-tissue behavior. They do not evaluate drug response or demonstrate paper-model availability.

## Recorded aggregate checks

These standard-library checks do not require source matrices or response labels:

```sh
python scripts/verify_public_snapshot.py
python scripts/verify_public_ordinal.py
python scripts/verify_public_ordinal_repeated.py
python scripts/verify_repro_assets_v1.py
```

The last command defaults to recorded summary/fixture contracts. Its optional `--local` mode needs separately prepared local assets. See [the asset guide](reproducibility_v1/README.md). A summary or inventory hash is not an accessible trained model.

## Frozen-model inference with the published release

Download the model ZIP from the release or Zenodo archive, verify its SHA-256
against the asset guide, then run from the v0.1.0 checkout or matching code ZIP:

```sh
python -m zipfile -e reo-model-assets-v0.1.0.zip release_assets/model
python scripts/predict_repro_v1.py --model release_assets/model/models/d1930_response_association_pair_full.json.gz --input release_assets/model/synthetic_inputs/full_gene_order.synthetic.json --output released_model_synthetic_predictions_new.json
```

Use a new output path. This replays an existing fitted model on generated
numbers; it reads no drug-response labels, fits no parameters and does not
recompute a scientific performance metric. Metadata CSVs are not required for
this inference CLI. The metadata package supplies frozen project identities
and assignments for interpretation, not the full source-pointer inputs needed
by the historical external runner. Preserve all zero molecular solutions and
abstention rules.

## Scientific reconstruction

Full reconstruction needs the precise upstream versions, identity/group mappings, folds, source-access permissions and environments described in [data availability](docs/DATA_AVAILABILITY.md) and [reproduction status](docs/REPRODUCIBILITY.md). Obtain upstream expression/response data independently under their own terms. The published model and minimal metadata packages do not include raw matrices, original external extraction pointers or all legacy curation dependencies, so a checkout plus these packages alone cannot replay every historical audit or all paper results. The model/derived-asset CC BY 4.0 grant covers only rights the authors own; it does not relicense upstream material.

Authoritative designs are the [development protocol](docs/ORDINAL_REPEATED_V1_PROTOCOL.md), [LOTO protocol](docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md) and [external protocol](docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md). Historical commands remain in the reproduction-status document; use a separate authorized workspace for a new reconstruction. They are not instructions to alter the completed run.

The unique CTRP external evaluation is complete. Individual response files, indices and predictions remain sealed during publication work. These checks never repeat external scoring or use its result for tuning. Conditional intervals exclude repeating the entire feature-selection/tuning/fitting procedure; NA and zero-model results must remain.
