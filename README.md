# REO drug-response evaluation

A leakage-controlled comparison of expression, within-sample ranks and relative-expression gene pairs for continuous cancer-cell-line drug-response prediction. The study separates tissue-reference performance, molecular increments and representation differences under matched candidate rules.

**Access:** this repository is private (verified 26 October 2026). Tracked code and aggregate results are available to authorized collaborators, not anonymously. The 64 fitted models and complete derived identity/fold assets are prepared locally but are not downloadable here. See [data and access conditions](docs/DATA_AVAILABILITY.md).

## Study and results

- Development: CCLE expression and GDSC2 responses, 626 biological groups, eight fixed compounds, 5 outer folds × 3 repeats and 3 inner folds.
- Primary development pair-minus-strong-expression macro ΔSpearman: −0.00588, family-adjusted conditional interval [−0.01871, +0.00661].
- Locked external evaluation: 175 non-overlapping groups, 1,185 drug–group responses; primary macro ΔSpearman −0.000372 [−0.042166, +0.043505].
- No average pair advantage was established. These intervals establish neither equivalence nor non-inferiority. Both screens use the same CCLE expression source: this is retrospective cross-screen ranking transport, not cross-expression-platform or clinical validation.
- Zero molecular solutions and undefined correlations remain: 430/840 molecular outer fits were zero; the strict leave-one-tissue-out macro is NA (34/80 display contrasts defined).

![Study design and evaluation boundaries](figures/manuscript_m08_v1/Figure1_study_design.png)

The manuscript is unpublished. The current local reading copy is English v4 with a paired Chinese copy; author packages are Git-ignored. [Stage 31](reports/阶段31_Minor1与2最小补充.md) records the two latest clarifications. There is no paper DOI or archived release yet. The sole maintained plan is [RESEARCH_PLAN.md](RESEARCH_PLAN.md).

## Quick start

Use Python 3.12 in an isolated environment. From the repository root:

```sh
python -m pip install -r requirements-repro-inference.txt
python scripts/predict_repro_v1.py --model examples/repro_inference_v1/toy_rank_model.json --input examples/repro_inference_v1/toy_input.json --output toy_predictions_new.json
```

The output path must be new. Parameters and inputs are handwritten; this is not a paper model or a drug-efficacy result. It reads no response labels and performs no fitting. See the [worked example](examples/repro_inference_v1/README.md).

## Reproduction checks

```sh
python scripts/verify_public_snapshot.py
python scripts/verify_public_ordinal.py
python scripts/verify_public_ordinal_repeated.py
python scripts/verify_repro_assets_v1.py
python -m unittest discover -s tests -p "test_repro_inference_v1.py" -v
```

The commands named `verify_public_*` inspect release-suitable aggregates; their names do not imply public repository access. They verify recorded contracts and toy inference, not scientific retraining. [REPRODUCING.md](REPRODUCING.md) explains dependency layers, training protocols, model-access limits and the completed external-evaluation boundary. CI checks aggregates and inference without opening sealed CTRP responses.

## Data and model access

| Material | Current location and access |
|---|---|
| Code, configurations, aggregate tables and figures | Tracked here; private repository access required |
| Upstream CCLE, GDSC2, CTRP v2 and annotations | Versioned sources and applicable terms in [data availability](docs/DATA_AVAILABILITY.md) |
| Actual fitted models, preprocessing and coefficients | Local, Git-ignored; permitted reviewer/public access still to be arranged |
| Complete derived identities and fold definitions | Local, Git-ignored; see [minimal reproduction assets](reproducibility_v1/README.md) |
| Raw matrices, individual responses/predictions and checkpoints | Not distributed in this repository |
| Manuscript and submission materials | Local author package; no public release or all-author submission approval claimed |

## Project map

| Path | Purpose |
|---|---|
| `src/`, `scripts/`, `tests/` | Curation, modeling, inference, analysis and checks |
| `config/` | Frozen development and external contracts |
| `data_manifest/` | Source versions and provenance |
| `results/` | Aggregate outputs supporting the paper |
| `figures/manuscript_m08_v1/` | Four current figures and editable sources |
| `reproducibility_v1/` | Numeric model contract and toy inference |
| `docs/`, `reports/` | Protocols, evidence, audits and history |

Start with the [development protocol](docs/ORDINAL_REPEATED_V1_PROTOCOL.md), [external protocol](docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md), [statistical assessment](docs/ORDINAL_REPEATED_ASSESSMENT.md) and [repository audit](docs/repository_release_v1/REPOSITORY_AUDIT.md). Detailed history remains in stage reports and the [preserved preceding README](docs/repository_release_v1/README_before_cleanup_20261026.md).

## Citation, contribution and license

Use [CITATION.cff](CITATION.cff), identifying the exact commit used. The development version is not an immutable release; publication/archival metadata will be added only when it exists. Contributor citation metadata remains collective pending author-side confirmation.

[CONTRIBUTING.md](CONTRIBUTING.md) describes contribution boundaries. Newly authored reconstruction code and documentation use [MIT](LICENSE); upstream datasets, original manuscripts and legacy code are not relicensed.
