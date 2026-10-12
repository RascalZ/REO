# Reproduction status and commands

## Verified versus pending

- Implemented: original input hashes; expression audit; authoritative identity/relatedness graph; compound support and conflicts; observation indices; prespecified coverage panel; group-isolated external manifests; independent source-row validation and immutable local lock.
- Implemented in stage 2: training-only feature/scaling/tissue/target states, grouped 3-outer/3-inner development CV, exact GDSC source-row/dose/response checks, two-drug EN/tissue/mean comparisons and resource audit. See the stage-2 report; these are exploratory engineering outputs, not confirmatory performance estimates.
- Pending: primary cancer-model eligibility revision, full REO200/RF/5x3 repeated analysis, expanded regularization path, uncertainty/permutation inference, final external predictions and manuscript results. There are no final trained-model checkpoints or external results yet.
- Local environment: Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, openpyxl 3.1.5, RDKit 2026.9.1. Expression mmap files are two float64 arrays of 1,019 × 57,820, about 471 MB each. No GPU is needed for phase 1. Model runtime/RAM must be measured in the next development pilot.

## Quick inspection without original data

```sh
git clone https://github.com/RascalZ/REO.git
cd REO
python scripts/verify_public_snapshot.py
```

This verifies published aggregate hashes/counts only, not individual identity or prediction performance.

This quick inspection and syntax compilation were tested on 2026-10-10 in a fresh Git clone using isolated Python standard-library execution (`-I -S`). Full source replay and a newly installed dependency environment have not yet been verified.

## Full local phase-1 reconstruction

Use Python 3.12 and a fresh checkout/workspace; do not overwrite a locked workspace. Prepare the original versioned inputs listed in DATA_AVAILABILITY. Set `REO_ORIGINAL_ROOT` to the parent original-data directory or copy `config/input_paths.example.json` to `config/input_paths.local.json` and edit the file paths. Only the latter is machine-specific and Git-ignored.

```sh
python -m venv .venv
# Activate .venv using your shell's normal activation command.
python -m pip install -r requirements.txt
python scripts/run_curation.py --download --lock
```

The orchestrator runs inventory → expression/reference audit → Sanger/CTRP acquisition → drug identity → cell identity → coverage and freeze → independent validation. Required original project inputs are described explicitly; a clean checkout installation/source replay must be tested before a publication release is claimed. Mutable Sanger metadata is cached as a dated local snapshot; replay against a later API snapshot constitutes a new cohort audit, not exact reproduction of cohort v1. `--download` does not overwrite present CTRP members, and source downloads should never overwrite original inputs.

For the already locked local workspace:

```sh
python src/validate_cohorts.py
python src/freeze_cohorts.py
```

The second command verifies recorded dependencies and exits without rewriting a locked cohort. Changes to frozen curation code, mappings, rules or inputs require a new cohort version and an entry in the total research plan. Adding independent model scripts does not authorize changes to cohort v1.

## Outputs and scientific boundaries

Development pilot entry points:

```sh
python -m pip install --index-url https://pypi.org/simple -r requirements-modeling.txt
python -m unittest discover -s tests -p test_pilot_integrity.py -v
python src/run_development_pilot.py --run-dir runs/new-development-pilot
python scripts/validate_pilot_run.py --run-dir runs/new-development-pilot
python src/report_development_pilot.py --run-dir runs/new-development-pilot
```

Execution commit `629b8f5` produced the recorded v0.1 results. The current reader includes a subsequent validation-failure handle cleanup, and the launcher now records failure/resource state and partial outer-fold checkpoints. These repairs do not change the successful v0.1 numerical execution. Full run inputs, individual predictions and fit-membership traces stay local; aggregate evidence is under `results/pilot_development/`.

The authoritative phase-1 outputs are the local `config/cohort_freeze_v1.json`, `.lock.json`, `reports/validation_summary.json`, `results/*sample_manifest.csv` and aggregate `main_drugs.csv`. Selection does not use response magnitudes. External response values are evaluated only after development/model/statistical protocol freeze.

The historical v1 external pool had 180 models. Cancer-focused v2 excludes five explicit non-cancerous controls and verifies HCC1588 through exact authoritative IDs, leaving 175. It preserves separation from all **GDSC2** identities and known related groups, including conflicting identity candidates. Unknown genealogy remains a limitation. CTRP was used historically, so this is a retrospective redesign. GDSC2 and CTRP share the CCLE expression source and differ in dose windows: the planned external question concerns cross-screen drug-response rank transportability, not cross-expression-platform performance or equal absolute AUC calibration.

## Cancer-v2 / REO200 and RF pilot

```sh
python scripts/validate_cancer_cohort.py
python -m unittest discover -s tests -p "test_pilot*_integrity.py" -v
python src/run_development_pilot.py --config config/pilot_v0.2.json --run-dir runs/new-cancer-v2-pilot
python scripts/validate_pilot_run.py --run-dir runs/new-cancer-v2-pilot
python src/report_reo200_rf_pilot.py --run-dir runs/new-cancer-v2-pilot
```

Commit scientific source/configuration changes before fitting. The launcher automatically records and checks Git HEAD (or an explicitly provided REO_EXECUTION_COMMIT) so immutable source validation does not rely on a forgotten environment variable. The completed 906-fit execution uses `5d05349`. Later exact-constant tissue-target maintenance applies only to a degenerate case absent in every actual training fit, and does not overwrite historical predictions. Full independent audit and aggregate results are under the local run directory and `results/pilot_reo200_rf/`; the source cache for cancer eligibility is versioned under `data_raw/cancer_metadata/`.

Both primary REO EN drug models selected a zero molecular endpoint in all three outer folds. This pilot supports no claim of reliable REO molecular advantage; full eight-drug repeated analysis/statistics and sealed CTRP evaluation remain pending.

## Fixed eight-drug formal development evaluation

The formal v1 lock binds scientific model/statistics configuration, protocol,
execution and analysis code, and cancer-v2 identity inputs. Prior pilots were seen;
this is a pre-formal-run lock, not prospective registration of the entire study.

```sh
python -m unittest discover -s tests -p "test_formal*.py" -v
python scripts/freeze_formal_analysis.py --verify-only
python src/run_formal_evaluation.py
python scripts/validate_formal_run.py --run-dir runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1
python src/analyze_formal_results.py --run-dir runs/gdsc2-cancer-v2-fixed8-nested5x3-inner3-formal-v1
```

Use `--resume` explicitly after interruption, with `REO_EXECUTION_COMMIT` set to
the run's recorded commit if HEAD has advanced. Completed checkpoint checksums,
registered file hashes, bound prepared inputs and frozen scientific versions must
match. Previous failed attempts are archived; they are never silently erased.
The resource scheduler caps two workers, with one BLAS/RF thread each, and admits
a second task only with at least 3 GiB available RAM. Dispatch does not change
folds, seeds, candidates or estimands.

840 outer tasks produce 18,120 fits and 104,307 development OOF rows. Three repeats
retain 626 independent development groups; rows are not independent extra N.
Paired bootstrap intervals are conditional on recorded OOF outcomes/predictions,
without train-refit uncertainty. Individual records, complete features, bootstrap
group IDs and count matrices stay local. Aggregate tables and reports can be
versioned. CTRP response access remains forbidden throughout this work package.
