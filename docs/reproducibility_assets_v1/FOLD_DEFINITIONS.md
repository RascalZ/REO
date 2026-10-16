# Frozen fold definitions and existing-manifest audit

Prepared 2026-10-16 for the v3.36 reproducibility-asset task. This document describes existing scientific assignments; it does not create new folds. The audit used only development eligibility/identity manifests, saved fold manifests, configuration, source code, metadata and already published aggregate LOTO results. No response table, individual prediction, bootstrap payload, model payload or CTRPv2 response/index file was opened. No splitter, fitting routine or scoring routine was called.

`fold_sources.csv` gives the relative paths, current SHA-256, sizes, CSV row counts and schema, together with available metadata-bound hashes. This is a local asset inventory; it does not decide redistribution permissions or certify a complete remote data archive.

## Audit result and direct reuse

All four saved fold CSVs match their own completed-run metadata's `prepared_input_hashes` byte for byte. They can be copied as existing evidence without rerunning scikit-learn or reconstructing random assignments.

| Stage and existing file | Data rows, excluding header | SHA-256 |
|---|---:|---|
| Repeated development: `global_fold_manifest.csv` | 31,926 | `7b519fc04b713ca69cf013a90e7fe9e31479944ad3f7b355f098e6ece895dcf4` |
| Full-development final fitting: `global_inner_fold_manifest.csv` | 1,878 | `68fa4a7b6274a40ef066f97f0f44d8f75bfd51ecbf940755b7fc2ddbcea5fdc9` |
| LOTO v1.3: `global_fold_manifest.csv` | 14,398 | `73c6b6aeecfd3d3c9528758d51526abdd3d87cb8002c4abc3ee469a8e99f38db` |
| LOTO v1.3: `loto_inner_fold_manifest.csv` | 109,274 | `0ef218bfec2cc9d0e3012f1f4bd3aab61867d2b3c10fa93385ad032ad35e96b8` |

The 4,967-row development eligibility manifest identifies 626 unique canonical models, 626 biological groups and 23 tissues. The audit confirmed a one-to-one canonical-model/group relationship in this cohort and consistent tissue/group attributes for every row of all four fold manifests. Biological groups encode the frozen known-relationship policy; this is not proof that every unknown cell-line relationship has been discovered.

By inspecting the stored assignments as sets, the audit confirmed:

- Development: complete disjoint train/test membership in all 15 outer contexts; complete disjoint inner train/validation membership in all 45 inner contexts; no outer-test model in an inner context; exactly one outer-test membership per model per repeat and exactly one inner-validation membership per outer-training model per outer context.
- Final fitting: three complete disjoint train/validation partitions and exactly one validation membership per global model.
- LOTO: 23 complete tissue outer partitions, each global model held out once; 184 drug-by-held-out-tissue inner contexts; each context contains exactly that drug's outer-training models, with one validation assignment per model, all three fold indexes present, and no held-out-tissue model included.
- Intersecting the existing development and final global assignments with each drug's existing eligibility IDs preserves disjointness and complete coverage. This verification performed set operations on stored membership, not fresh splitting.

These checks do not independently revalidate historical fitted checkpoints, coefficients or predictions. The completed LOTO validation record reports checkpoint validation from the original execution; that claim remains attributed to that record.

## 1. Canonical identities and per-drug eligibility

The source identity/eligibility table is `results/cohort_v2/development_sample_manifest.csv` (4,967 rows, SHA-256 `21efc28d5ee21603f2466b740c4aa29eb24fbfd6be4fec15f6392aadb661e608`). It contains no numeric response column. Its required linking fields are:

`gdsc_drug_id`, `drug_name`, `canonical_model_id`, `biological_group_id`, `ccle_id`, `expression_row_index`, `tissue`, `role`, `cohort_version`.

The original also retains source curve/row, dose-window, eligibility and evidence fields. For exact reproduction preserve the full existing file locally; a future thin public projection must document its omitted fields and original hash. `expression_row_index` is an offset into the project-specific prepared expression matrix, not a universally portable CCLE index; the matrix/sample ordering must be supplied separately.

The repeated and final splitters first deduplicate this table by canonical ID and sort canonical IDs. They partition the global 626-group universe, then intersect the same membership with each drug's eligibility set. The arm is not part of the fold key: all eight model arms share each applicable drug partition. Neither response ranks nor prediction performance defines these assignments.

| GDSC drug ID | Drug | Eligible development groups | Repeated outer-test group range | Repeated inner-validation group range | Final inner-validation counts, folds 0/1/2 |
|---|---|---:|---:|---:|---|
| 1931 | MIRA-1 | 625 | 124–126 | 165–168 | 209 / 208 / 208 |
| 1036 | PLX-4720 | 622 | 123–126 | 163–167 | 208 / 207 / 207 |
| 1191 | Bortezomib | 620 | 122–125 | 162–167 | 208 / 207 / 205 |
| 1930 | Telomerase Inhibitor IX | 620 | 122–125 | 162–167 | 208 / 207 / 205 |
| 1199 | Tamoxifen | 620 | 122–125 | 162–167 | 208 / 207 / 205 |
| 1091 | BMS-536924 | 620 | 122–125 | 162–167 | 208 / 207 / 205 |
| 1010 | Gefitinib | 620 | 122–125 | 162–167 | 208 / 207 / 205 |
| 1922 | Cediranib | 620 | 122–125 | 162–167 | 208 / 207 / 205 |

The repeated ranges summarize the actual stored 15 outer contexts and 45 inner contexts after intersection; they are not target sample sizes.

## 2. Repeated development: five outer folds × three repeats, three inner folds

Scientific source commit: `5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43`.

Configuration/lock: `config/ordinal_repeated_v1.json` and `config/ordinal_repeated_v1.lock.json`. Execution source: `src/run_ordinal_repeated.py`, with the fold preparation/intersection functions in `src/run_formal_evaluation.py` and grouped splitter in `src/modeling/nested_cv.py`.

Existing fold path:

`runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1/global_fold_manifest.csv`

The frozen algorithm was `StratifiedGroupKFold(n_splits=..., shuffle=True, random_state=seed)`, tissue as the stratification label and biological group as the grouping unit. `repeat` and `outer_fold` are zero based. The three outer seeds were 20261013, 20261014 and 20261015. Inner seeds were `repeat_seed + 100 + outer_fold`, evaluated only on the current global outer-training partition. Fold preparation used the recorded scikit-learn 1.7.2 environment; exact saved memberships take precedence over trying to reproduce assignments with a newer library.

CSV schema:

`canonical_model_id, biological_group_id, tissue, repeat, outer_fold, inner_fold, role`

For `inner_fold=-1`, role is `train` or `test` and covers all 626 global IDs once in each outer context. For `inner_fold=0,1,2`, role is `train` or `validation` and covers only that outer-training set once per inner fold. Do not treat all rows as distinct observations. The data-row count is 9,390 outer-membership rows plus 22,536 inner-membership rows. Global outer-test counts are 125–126; global inner-validation counts are 166–168.

The recorded preparation retained 18 stratification warning occurrences: three outer warnings and fifteen inner warnings, because the least-populated tissue had one or two groups. There is one global Biliary Tract group and two Cervix groups. Tissue stratification cannot guarantee every rare tissue appears in every train or validation partition. No tissue was dropped or fold reassigned in response to scores. These warnings are not model fitting failures or elastic-net solver warnings.

The stored completed-run contract is 960 model/outer tasks, 26,160 fit events and 119,208 OOF rows. Those are metadata counts, not newly regenerated outputs. The 15 outer contexts are not independent datasets; repeated observations of a biological group must not be counted as independent experimental units.

## 3. Full-development final fitting: shared three-fold tuning

Scientific source commit: `ca0f837ac10665ab8d4d6ab97691e58e5bed4577`.

Configuration: `config/external_validation_v1_1.json` (`final_training`). Relevant source: `src/final_model_v1.py` (`make_global_inner_manifest`, `local_inner_splits`) and `src/run_final_model_v1.py`.

Existing fold path:

`runs/gdsc2-cancer-v2-fixed8-full-development-final-v1_1/global_inner_fold_manifest.csv`

The same sorted 626-group identity universe received one shared three-fold tissue-stratified grouped assignment using seed 2026101501. The three global validation sets contain 209, 209 and 208 groups. The assignment is intersected with each drug and shared across arms. Once tuning selected an existing hyperparameter configuration, the historical execution fitted its final model on all eligible development records; the inner folds are tuning partitions, not a final model test set.

CSV schema:

`canonical_model_id, biological_group_id, tissue, inner_fold, role`

`inner_fold` is zero based, roles are `train` and `validation`. Each model occurs in all three fold partitions but is validation exactly once. The recorded preparation retained one warning: the least-populated tissue had one group, fewer than three folds. The completed-run metadata records 64 final artifacts and 1,744 fit events. This asset-preparation task does not reopen model payloads or repeat that fitting.

## 4. LOTO: hold out a complete tissue and tune within the remaining drug cohort

Scientific source commit: `fedc745940374b6870cb6c4bcb66ec8f731ded27`; completed execution version: P5e-LOTO-v1.3.

Configuration/lock: `config/loto_execution_v1_3.json` and `config/loto_execution_v1_3.lock.json`, retaining the earlier `config/tissue_sensitivity_loto_v1.json` scientific definitions. Relevant source: `src/run_tissue_loto_v1_3.py`.

Existing paths:

- `runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_3/global_fold_manifest.csv`
- `runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_3/loto_inner_fold_manifest.csv`

The 23 frozen tissue names are ordered by Python Unicode code-point order and indexed 0–22. Every group of the current tissue is held out; the global outer assignment is shared by all drugs and model arms. The inner assignment differs from the repeated-development/final policy: it is made on **each drug's available outer-training records**, not on a global 626-group inner split followed by intersection. All arms for that drug-by-held-out-tissue context share its saved assignment.

The inner seed is `2026101401 + tissue_index`. The original frozen prose had instead stated base 2026101801; amendment v1.1 resolved that conflict before LOTO model fitting in favor of the frozen machine-config value 2026101401. v1.2/v1.3 preserved this resolution; v1.3 included the already frozen `inner_seed` in prediction identity fields after a checkpoint-contract repair. Do not erase that history or label the execution as an unchanged, internally consistent v1 protocol.

Outer CSV schema:

`canonical_model_id, biological_group_id, tissue, repeat, outer_fold, inner_fold, role, heldout_tissue`

`repeat=0`, `outer_fold=tissue_index`, `inner_fold=-1`; role is `test` exactly when `tissue==heldout_tissue`, otherwise `train`. There are 626 × 23 = 14,398 rows.

Inner CSV schema:

`drug_id, tissue_index, heldout_tissue, canonical_model_id, biological_group_id, tissue, validation_fold, inner_seed`

This is a compact validation-assignment table: each drug's outer-training model has one `validation_fold` value 0/1/2. For a given inner fold j, rows with `validation_fold==j` form validation and the other rows in the **same drug-by-held-out-tissue context** form training. It does not list outer-test models. It contains 109,274 assignments over 184 contexts; context validation-set sizes are 162–208. Drug-by-held-out-tissue outer tests contain 1–132 groups.

The frozen preflight recorded one rare-stratum warning in every one of the 184 contexts: 176 had a smallest remaining tissue with one group; eight had two groups. The three grouped folds were retained, and the present set audit confirms complete, separated saved memberships. This does not remove the instability of tuning for rare tissue strata.

All 23 tissues and all eight arms were retained in the historical run (1,472 tasks, 40,112 fits, 39,736 OOF rows). The ten-tissue display/strict macro set was prespecified by at least 20 groups per drug: Bowel, Breast, CNS/Brain, Esophagus/Stomach, Lung, Lymphoid, Myeloid, Ovary/Fallopian Tube, Pancreas and Skin. It is not the set on which folds were fitted; the other 13 tissues remain in the outputs.

The LOTO descriptive Spearman helper returns NA for fewer than two held-out observations, constant response or constant prediction. An unseen tissue has an all-zero tissue one-hot vector, so an otherwise zero molecular solution predicts the training overall mean and has undefined within-tissue Spearman. Preserve the NA and its reason. Any missing component makes the ten-tissue drug macro NA, and any undefined drug macro makes the strict eight-drug macro NA. The saved aggregate retains 34/80 defined primary contrasts and an NA strict macro. Available-cell means and concatenated cross-tissue correlations are separate descriptive quantities and cannot replace that result. No LOTO confidence-interval/testing family was defined.

## 5. Minimal local package and future release schema

For this round, direct byte copies of the four fold CSVs are sufficient as local membership assets. Use separate unambiguous names such as `development_repeated_global_folds.csv`, `final_development_inner_folds.csv`, `loto_global_outer_folds.csv` and `loto_drug_inner_assignments.csv`; retain original relative path and SHA-256 in a package manifest. Preserve UTF-8 BOM, line endings, column order and row order if declaring a byte-identical copy. No `make_grouped_splits`, `StratifiedGroupKFold`, random generator or frozen runner needs to be invoked.

The minimal interpretation package also needs: the 4,967-row development eligibility manifest (or a documented ID-only projection), a 626-row identity/tissue inventory derived without altering labels, the frozen configurations/locks or source hashes, and the source ledger. Identity-only projections are derived convenience views and must not be mislabeled as original frozen files. Numeric development responses and expression inputs belong to separately inventoried assets; they are not needed to interpret fold membership or demonstrate label-free inference.

A future release manifest should include `asset_id`, `stage`, `original_relative_path`, `package_relative_path`, `sha256`, `bytes`, `data_rows`, `columns`, `execution_commit`, `metadata_relative_path`, `metadata_bound_sha256`, `response_columns_present`, and `redistribution_status`. For these fold CSVs `response_columns_present=false`; redistribution status remains pending author/license review. CTRPv2 identities, labels and individual external predictions are unnecessary for interpreting development folds. The separate local metadata export in section 6 retains existing external identities and source pointers, without response labels or individual predictions.

Reproducibility of fold assignment does not itself establish absence of all other leakage. Train-only candidate selection, preprocessing, tuning and model serialization are separate audited components; exact membership merely supplies the boundary on which those claims can be checked.

## 6. Executed local metadata export

`scripts/export_repro_metadata_v1.py` implements the limited export described above. It was executed with `.venv/Scripts/python.exe -X utf8 scripts/export_repro_metadata_v1.py`; output is `reproducibility_v1/local_assets/metadata`. The runtime audit guard permits input opens only for the two existing cohort manifests, four fold CSVs and three run metadata JSONs. It prohibits other file reads, network/process calls, replacement/deletion, and file writes outside the new output directory. Output files use exclusive creation; a second invocation into the existing directory was checked and refused before any export read or write.

The external cohort manifest was first inspected through an unbuffered header-only read. Its columns contain source identity, experiment/row references, identity QC and aggregation-policy metadata, without a response-value column. Only after passing that positive field whitelist was its body read. Response-presence fields, if present in a permitted source schema, are excluded from all output projections. No numeric CTRPv2 label or individual prediction was read.

The completed package contains 626 development identities/4,967 drug memberships and 175 external identities/1,185 drug memberships, with no canonical-model or known-biological-group overlap. Exact existing CCLE IDs and source identifiers are retained; no new name match, inferred cross-reference or download is introduced. All four folds were copied byte for byte with their original metadata-bound hashes. `mapping_schema.json` specifies keys, field meanings and the pending redistribution status; `metadata_export_manifest.json` records source/output hashes and columns. The public copy `docs/reproducibility_assets_v1/metadata_export_public_summary.json` contains counts and hashes, without individual model IDs.

The initial launch stopped before reading inputs or writing output files because the guard did not permit creation of the missing `reproducibility_v1` ancestor directory. The directory-only rule was narrowed to the chosen output's repository-local ancestors while retaining exclusive file writes inside output. The subsequent export completed with zero blocked access attempts. No scientific file, fold or sample selection changed during this implementation repair.
