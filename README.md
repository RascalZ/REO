# REO drug-response prediction: auditable reconstruction

This repository contains a leakage-controlled comparison of relative expression ordering (REO), expression, rank, and tissue-only models for continuous drug-response prediction in cancer cell lines. The analysis is complete through its locked external evaluation. Manuscript reconstruction is in progress.

## Current status

As of 2026-10-20, the single maintained research plan is [RESEARCH_PLAN.md](RESEARCH_PLAN.md), version 3.39. Read its current-status section before using any older stage notes. The current scientific readouts are:

- Fixed eight-drug nested development evaluation: 626 unique development biological groups, 4,967 drug–model observations, 5 outer folds × 3 repeats and 3 inner folds. The predefined REO-pair minus strong-expression macro ΔSpearman was −0.00588, with family-adjusted conditional interval [−0.01871, +0.00661]. This does not establish REO superiority, equivalence, or non-inferiority.
- Tissue sensitivity: the strict prespecified leave-one-tissue-out macro was undefined; 34/80 drug–tissue cells had defined correlations.
- Final model audit: 64/64 models and 1,744 fit events passed independent reload/source audit. Of 56 molecular arms, 29 selected a zero molecular solution.
- One-time CTRPv2 evaluation: 175 non-overlapping external groups and 1,185 drug–model response units, with 9,480/9,480 predictions covered. The primary equal-drug REO-pair minus strong-expression ΔSpearman was −0.000372, family-9 adjusted conditional interval [−0.042166, +0.043505]. The result is a retrospective cross-screen transfer evaluation using the same CCLE expression source, not cross-expression-platform or clinical validation.
- The original manuscript’s “pan-cancer REO atlas” claims, binary classification results, pathway/PPI claims, GEO fine-tuning and single-cell predictions are not supported by the rebuilt analysis. See [Stage 15 manuscript reconstruction blueprint](reports/阶段15_投稿稿件重构蓝图.md).
- The complete current author-reading package is local and Git-ignored at `manuscript_v1/author_review_v1/` (main Markdown/Word/PDF, supplementary methods and Additional files). It incorporates the user-designated previous submission’s recorded author order, affiliation, equal-first/correspondence markers, competing-interest declaration and two funding sources. Author personal details stay in that local package. Missing emails, individual contributions, acknowledgements, permissions and all-author approval of the reconstructed manuscript are not invented. See [Stage 26](reports/阶段26_作者资料迁入与最终通读稿.md); 26 content checks passed; the 39-page main and 19-page supplement were rendered by Word and visually inspected, with layout/source checks recorded there. The preserved [R1-M08 scientific manuscript](manuscript_v1/Manuscript_BMC_Genomics_M08_v1.md), [Word](manuscript_v1/Manuscript_BMC_Genomics_M08_v1.docx), [supplementary methods](manuscript_v1/Supplementary_methods_scientific_revision_v1.md), and [Additional files 1–14](manuscript_v1/bmc_additional_files_M08_v1/) remain unchanged. The old local `with_figures` derivative is historical and must not be circulated as current.
- R1 adds all eight drugs' external absolute correlations (strict macros: tissue 0.2284, strong expression 0.2602, response-associated pair 0.2599), and full matched expression200/rank200/pair contrasts from saved aggregates. A near-zero pair-minus-strong difference does not imply disappearance of all predictive association. Uncomputed intervals remain NA; external rank contrasts are descriptive points only. Cross-screen uncalibrated RMSE/R-squared are not evaluated.
- All four finite synthetic controls are reported, including all 40 cohorts, 615/840 selected molecular zero outer fits, 267 warning occurrences, and undefined statistics. The synthetic strong-expression arm actually uses 220 genes, not 2,000. These constructed illustrations do not identify real-data mechanisms or substitute for eight-drug full-pipeline label permutation, which was not run.
- The [Stage 23 R1 record](reports/阶段23_科学稿第一轮修订.md) and [revision audit directory](docs/scientific_revision_v1/) document source-to-text reviews, Methods additions, first-appearance reference renumbering and package checks. These are coverage-bounded audits, not independent raw-data/model reproduction. M01–M03 are addressed for R1; M04 local technical assets and response-free inference are now prepared; actual public model/mapping access remains open until the planned rights/release closing step. M08 figure/reference style is complete for the current review copy; see the [Stage 25 render record](reports/阶段25_M08图文样式与实际渲染验收.md) and [bounded checks](docs/m08_v1/). Recorded author/declaration facts have now been incorporated into the local reading package; genuinely absent facts, actual access arrangements, formal archiving, rights and approval of the revised manuscript remain closing tasks. The next gate is the user’s read-through and discussion; do not request facts already supplied in the designated historical manuscript again. The manuscript is not submission-ready. The unique locked CTRPv2 evaluation is complete; saved aggregates remain usable in writing, while individual responses/indices stay sealed. This editorial step performs no label access, retraining, resampling or rescoring.

## Actual repository access — 20 October 2026

Authenticated GitHub metadata confirms that `RascalZ/REO` is **private**; the anonymous repository endpoint is inaccessible. Code, aggregate tables, provenance and invented-input examples tracked here are not thereby publicly or anonymously accessible. Earlier “public” asset descriptions denote the intended release scope, not a verified current access route. Public/editor/reviewer access, actual model/mapping availability and formal archiving remain to be completed. This editorial work does not change repository visibility or release local assets. Public upstream resource availability is separate from project-repository access; see the [access correction](docs/DATA_AVAILABILITY.md).

The local author-review package is for the user’s read-through, not a record that all authors have approved submission. Its source manuscript and personal author data are not included in the tracked safe report. Scientific results, locked models/statistics and the M08 originals are preserved.

## Current manuscript style build

The [M08 figures](figures/manuscript_m08_v1/) use 170 mm canvases and 600 dpi PNGs, with editable SVG sources. Figures 3/4 have vertical panels; minimum font is 9.5 pt at 170 mm (about 8.72 pt in the embedded Word copy). All plotted points, intervals and matrices agree with the preceding figures. References retain their 30 identities and numbering; 19 author prefixes were normalized. Scientific prose, table cells, captions and all 14 Additional files remain unchanged.

Optional editorial dependencies are pinned in `requirements-manuscript-qa.txt` (Python 3.12.14); they do not redefine the historical modeling environment. Use that environment to run `python src/format_m08_reference_authors.py`, `python src/plot_manuscript_figures_m08_v1.py`, `python src/build_m08_package_v1.py`, then the numeric-reference validator targeting `Manuscript_BMC_Genomics_M08_v1.md` and `docs/m08_v1/citation_validation.json`, and `python src/validate_m08_v1.py`. Actual Word rendering uses `scripts/render_m08_word.ps1` on Windows with Word installed; `python src/audit_m08_render_v1.py` inspects its output. Set `--visual-confirmed` only after actually viewing the pages. Local PDFs and page PNGs are ignored by Git. Rebuilding changes DOCX ZIP metadata, so refresh provenance and render checks before circulating a rebuilt copy. These commands perform no model training or external scoring.

## Response-free inference and local reproduction assets

The [package guide](reproducibility_v1/README.md) and [numeric formats](reproducibility_v1/FORMAT.md) document minimal inputs, exact identity/fold definitions, full preprocessing, feature coefficients and zero/null policies. The public toy is hand-specified and is not a paper model. With Python 3.12 and the two pinned dependencies:

```sh
python -m pip install -r requirements-repro-inference.txt
python scripts/predict_repro_v1.py --model examples/repro_inference_v1/toy_rank_model.json --input examples/repro_inference_v1/toy_input.json --output toy_predictions_new.json
```

The output path must be new. The public example works without sklearn/joblib, old pipeline sources, local full assets or CTRP files. It produces four predictions and three explicit abstentions from seven invented rows. No label/scoring input is accepted. See the [Stage 24 record](reports/阶段24_最小复现资产与无标签推理.md) and [technical audit](docs/reproducibility_assets_v1/).

Complete local assets preserve 64 fixed numeric models, all 29/56 molecular zero solutions, 25,600 candidate rows, 340,800 coefficient rows, 1,472 tissue rows, 626/175 identities, label-free drug memberships and four original split CSVs. Native/portable outputs agreed on 384 generated cases; this is a numerical implementation check, not predictive efficacy or independent study replication. Actual derived models and mappings remain Git-ignored pending access/rights decisions; GitHub contains code, invented examples and hash/count summaries. CTRPv2 labels remain sealed.

## Research records

- [Total research plan](RESEARCH_PLAN.md): sole maintained plan, decision gates, history and change log.
- [Stage 14 external evaluation and scientific interpretation](reports/阶段14_P6一次性CTRP外验与研究方向判定.md): coverage, effect estimates, intervals and limitations.
- [Stage 15 manuscript reconstruction blueprint](reports/阶段15_投稿稿件重构蓝图.md): old-claim disposition, proposed research framing, abstract draft and figure/table map.
- [Stage 16 captions and evidence ledger](reports/阶段16_图注与稿件证据台账.md): English draft captions/table notes and claim-level evidence boundaries.
- [First-draft Figures 1–4](figures/manuscript_v1/): editable SVG and 600 dpi PNG; inputs, outputs and software versions are recorded in `figure_manifest.json`. Figure 4 includes the LOTO defined-cell panel and locked CTRPv2 comparison.
- [First-draft Tables 1–2 and S1–S8](tables/manuscript_v1/): 10 CSVs generated from frozen aggregates, with source/output hashes and row counts in `table_manifest.json`.
- [Manuscript working draft](manuscript_v1/Manuscript_draft.md): first-pass English text grounded in frozen protocols, implementation and aggregate results.
- Rebuild the R1 documents/package with `python src/build_scientific_revision_manuscript_v1.py` then `python src/build_scientific_revision_package_v1.py` (Node.js and `manuscript_v1` npm dependencies required). The old `npm run build` default targets the preserved preceding draft; do not use it to overwrite that historical copy during R1 review.
- Run `python src/validate_numeric_references.py --manuscript manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md --output docs/scientific_revision_v1/citation_validation.json` then `python src/validate_scientific_revision_v1.py` for R1 coverage/structure checks. These do not read external individual responses or fit models.
- [Bibliography and citation audit](manuscript_v1/references.bib): formatted reference metadata. Run `python src/validate_numeric_references.py` to validate numbered citation coverage against the manuscript and Crossref export.
- [NLM ISO journal-abbreviation map](manuscript_v1/nlm_journal_abbreviations.csv): source title, abbreviation, DOI, and field used for all 26 journal-article references.
- [Stage 18 literature and claim audit](reports/阶段18_引文元数据与主张核验.md): reference provenance, year/article-number reconciliation, and evidence boundaries.
- [Stage 19 full manuscript and guideline audit](reports/阶段19_全文主张代码与期刊规范二次核验.md) and [Stage 20 BMC structure/file package audit](reports/阶段20_BMC结构稿与文件包审计.md): evidence alignment, journal requirements, and submission gates.
- [Stage 21 bounded AIVC context and citation audit](reports/阶段21_AIVC讨论边界与引文整合.md): claim limits, balanced sources, and manuscript verification; [new-reference provenance](manuscript_v1/aivc_reference_metadata.json) records Crossref, PubMed and official conference metadata.
- Historical preceding [BMC Genomics draft (Markdown)](manuscript_v1/Manuscript_BMC_Genomics_draft.md), [editable Word draft](manuscript_v1/Manuscript_BMC_Genomics_draft.docx), and [Additional files 1–8](manuscript_v1/bmc_additional_files/). These author-input-required working files are not a final submission package.
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

The next work package is manuscript reconstruction from the locked evidence: manually audit the first-draft figures and tables, write Methods/Results with evidence links, then update the journal-fit decision. BMC Genomics remains a candidate for reassessment; publication timing and acceptance cannot be inferred from these results.
