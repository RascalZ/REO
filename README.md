# REO drug-response prediction: auditable reconstruction

This repository contains a leakage-controlled comparison of relative expression ordering (REO), expression, rank, and tissue-only models for continuous drug-response prediction in cancer cell lines. The analysis is complete through its locked external evaluation. Manuscript reconstruction is in progress.

## Current status

As of 2026-10-15, the single maintained research plan is [RESEARCH_PLAN.md](RESEARCH_PLAN.md), version 3.27. Read its current-status section before using any older stage notes. The current scientific readouts are:

- Fixed eight-drug nested development evaluation: 626 independent development groups, 4,967 drug–model observations, 5 outer folds × 3 repeats and 3 inner folds. The predefined REO-pair minus strong-expression macro ΔSpearman was −0.00588, with family-adjusted conditional interval [−0.01871, +0.00661]. This does not establish REO superiority, equivalence, or non-inferiority.
- Tissue sensitivity: the strict prespecified leave-one-tissue-out macro was undefined; 34/80 drug–tissue cells had defined correlations.
- Final model audit: 64/64 models and 1,744 fit events passed independent reload/source audit. Of 56 molecular arms, 29 selected a zero molecular solution.
- One-time CTRPv2 evaluation: 175 non-overlapping external groups and 1,185 drug–model response units, with 9,480/9,480 predictions covered. The primary equal-drug REO-pair minus strong-expression ΔSpearman was −0.000372, family-9 adjusted conditional interval [−0.042166, +0.043505]. The result is a retrospective cross-screen transfer evaluation using the same CCLE expression source, not cross-expression-platform or clinical validation.
- The original manuscript’s “pan-cancer REO atlas” claims, binary classification results, pathway/PPI claims, GEO fine-tuning and single-cell predictions are not supported by the rebuilt analysis. See [Stage 15 manuscript reconstruction blueprint](reports/阶段15_投稿稿件重构蓝图.md).

## Research records

- [Total research plan](RESEARCH_PLAN.md): sole maintained plan, decision gates, history and change log.
- [Stage 14 external evaluation and scientific interpretation](reports/阶段14_P6一次性CTRP外验与研究方向判定.md): coverage, effect estimates, intervals and limitations.
- [Stage 15 manuscript reconstruction blueprint](reports/阶段15_投稿稿件重构蓝图.md): old-claim disposition, proposed research framing, abstract draft and figure/table map.
- [Nested evaluation protocol](docs/ORDINAL_REPEATED_V1_PROTOCOL.md): frozen development estimands and training-fold rules.
- [External validation protocol](docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md): locked cross-screen evaluation rules.
- [Data availability and provenance](docs/DATA_AVAILABILITY.md) and [reproducibility status](docs/REPRODUCIBILITY.md).

## Key artifacts

| Artifact | Purpose |
|---|---|
| `config/external_validation_v1_1.lock.json` | First lock, binding final-model scientific sources and environment |
| `config/external_execution_v1_1.lock.json` | Second lock, binding final models and response-free external predictions before CTRPv2 response access |
| `results/formal_development/` | Aggregate results for the original fixed-eight development evaluation |
| `results/ordinal_repeated/` | Aggregate results for the two candidate rules × three representations evaluation |
| `results/tissue_loto/` | Aggregate tissue leave-out results, including undefined cells |
| `results/p6_external/evaluation_summary.json` | Public aggregate result from the unique completed external attempt; contains no individual response values |
| `results/p6_external/attempt_0001_completed.json` | Completion marker linked to the aggregate summary hash |
| `reports/` | Stage-by-stage scientific and operational records |

The original manuscript and submission packages remain outside the repository and are treated as read-only source material. Raw/processed expression matrices, individual sample-level predictions, fit checkpoints, local paths and response files are not committed. Upstream dataset terms remain in force; the repository license does not relicense those data.

## Inspecting public aggregates

With Python 3.12, the aggregate verifiers can be run from the repository root:

```sh
python scripts/verify_public_snapshot.py
python scripts/verify_public_ordinal.py
python scripts/verify_public_ordinal_repeated.py
```

These verify published cohort/model aggregate snapshots. They do not refit models or reproduce sample-level predictions. Full locked-model verification depends on the versioned input data and local environment described in the data-availability record.

The second-lock check can be run without opening CTRPv2 response values:

```sh
python scripts/freeze_external_execution_v1.py --verify-only
```

The CTRPv2 result was executed once under the verified second lock. Do not rerun the external evaluation as a tuning or model-selection step. No source code in this repository should be used to change the locked drug panel, predictions, estimands, or analysis after seeing the external results.

## Statistical interpretation

The external primary interval is conditional on the fitted models and observed external groups. It does not include the uncertainty from repeating feature selection, tuning and model fitting. The internal bootstrap intervals are likewise conditional on the recorded out-of-fold predictions. Repeated folds are not additional independent biological samples. A zero molecular solution means the penalized model retained the tissue component without a molecular coefficient; a correlation from that model does not demonstrate a learned molecular resistance signal. Undefined correlations and strict macro values remain NA.

The cross-screen assessment uses GDSC2 training responses and CTRPv2 evaluation responses, while expression comes from CCLE in both settings. CTRPv2 had been used in the historical manuscript, so the present analysis is not prospective or fully blind. Dose windows and screening protocols differ; the external target is within-drug ranking transportability, not absolute AUC calibration or patient response.

## Project structure

```text
RESEARCH_PLAN.md            # sole maintained plan and change log
config/                     # frozen model, analysis and execution configurations
data_manifest/              # source registry and data versions
docs/                       # protocols, data availability and reproducibility notes
reports/                    # scientific audits and stage reports
results/                    # public aggregate outputs; no sample-level labels
scripts/, src/              # curation, model, validation and analysis code
tests/                      # integrity and contract checks
```

The next work package is manuscript reconstruction from the locked evidence: review the new claim, build replacement figures and tables, write Methods/Results, then update the journal-fit decision. BMC Genomics remains a candidate for reassessment; publication timing and acceptance cannot be inferred from these results.
