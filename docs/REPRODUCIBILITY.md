# Reproduction status and commands

## Verified versus pending

- Implemented: original input hashes; expression audit; authoritative identity/relatedness graph; compound support and conflicts; observation indices; prespecified coverage panel; group-isolated external manifests; independent source-row validation and immutable local lock.
- Pending: leakage-free feature/model engine, nested grouped CV, continuous response parsing/direction checks, baselines, statistical uncertainty, final external predictions and manuscript results. There are no trained-model checkpoints or valid new performance numbers to reproduce yet.
- Local environment: Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, openpyxl 3.1.5, RDKit 2026.9.1. Expression mmap files are two float64 arrays of 1,019 × 57,820, about 471 MB each. No GPU is needed for phase 1. Model runtime/RAM must be measured in the next development pilot.

## Quick inspection without original data

```sh
git clone https://github.com/RascalZ/REO.git
cd REO
python scripts/verify_public_snapshot.py
```

This verifies published aggregate hashes/counts only, not individual identity or prediction performance.

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

The authoritative phase-1 outputs are the local `config/cohort_freeze_v1.json`, `.lock.json`, `reports/validation_summary.json`, `results/*sample_manifest.csv` and aggregate `main_drugs.csv`. Selection does not use response magnitudes. External response values are evaluated only after development/model/statistical protocol freeze.

The 180 external models are separated from all **GDSC2** identities and known related groups, including conflicting identity candidates. Unknown genealogy remains a limitation. CTRP was used historically, so this is a retrospective redesign. GDSC2 and CTRP share the CCLE expression source and differ in dose windows: the planned external question concerns cross-screen drug-response rank transportability, not cross-expression-platform performance or equal absolute AUC calibration.
