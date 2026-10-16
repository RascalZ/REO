# Label-free inference example

These five gene names, tissue names, TPM values and model numbers were manually invented. This fixture was **not fitted**, is **not a paper model**, and contains no development or CTRP observations, responses, sample predictions, or efficacy metrics.

Install `requirements-repro-inference.txt` in an appropriate Python environment, then run from the repository root:

```powershell
python scripts/predict_repro_v1.py --model examples/repro_inference_v1/toy_rank_model.json --input examples/repro_inference_v1/toy_input.json --output toy_predictions_new.json
```

The output path must be new. The checked-in `toy_expected_output.json` is a deterministic software-check result, not evidence about drug sensitivity. The interface accepts no response-label field and offers no fitting or scoring command. Only JSON/JSON.gz numeric model/input files are loaded; sklearn and joblib are unnecessary.

The example has seven rows: two observed tissues; a globally known tissue absent from this model's training categories; missing TPM; missing tissue; a negative TPM in an unselected gene; and a new tissue label. The last is numerically unseen, since this interface expects already canonicalized exact labels and performs no fuzzy matching. It does not establish that a newly supplied biological model is eligible for external validation.

`toy_1` has tied candidate-gene TPM values. Only four selected genes determine its average ranks; the fifth gene is still part of the complete schema and validity checks. Its first-gene rank is 2.5, so the handwritten FWL prediction is `0.5 + 0.2 * (((2.5 - 2.5) - 0.1) * 0.25 - 0.5) = 0.395`. `toy_2` predicts 0.68; the two unseen-tissue cases predict 0.5 and 0.575. Three invalid rows remain present with null predictions and reason codes.

Input keys are exactly `schema_version`, `sample_ids`, `biological_group_ids`, `tissues`, `gene_ids`, `raw_tpm`. IDs must already denote canonical samples and biological groups; duplicate valid IDs are hard failures. This module does not infer parent/derivative relationships, compound mappings, cancer qualification, cohort eligibility, or train/external overlap. Those are separate documented curation gates. Missing/invalid identity or tissue abstains in every arm. A true tissue-null arm may predict despite invalid molecular TPM; a zero-coefficient molecular arm must abstain on invalid TPM.

The sealed-path check blocks known CTRP/raw/processed/external-run/legacy-response paths before file opening. It is a local guard, not an OS access-control boundary or a guarantee about arbitrarily renamed file contents. Use only the explicit synthetic fixture or documented response-free inputs.
