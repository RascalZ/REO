# Reproduction status and commands

## Verified versus pending

- Implemented: original input hashes; expression audit; authoritative identity/relatedness graph; compound support and conflicts; observation indices; prespecified coverage panel; group-isolated external manifests; independent source-row validation and immutable local lock.
- Implemented in stage 2: training-only feature/scaling/tissue/target states, grouped 3-outer/3-inner development CV, exact GDSC source-row/dose/response checks, two-drug EN/tissue/mean comparisons and resource audit. See the stage-2 report; these are exploratory engineering outputs, not confirmatory performance estimates.
- Completed: cancer-model eligibility v2 (626 development groups and 175 sealed external groups), expanded regularization path and REO200/RF pilot, fixed eight-drug 5-outer × 3-repeat / 3-inner evaluation, and conditional paired-bootstrap analysis. Formal v1 remains unchanged.
- Completed: four-drug 3-outer × 1-repeat / 3-inner diagnostic with two candidate rules, three molecular representations, the same-policy tissue reference and nine full-pipeline within-tissue permutations. Independent validation, frozen statistics and reporting passed. The successful design contains 840 outer tasks, 22,440 fits and 173,600 OOF rows; operational overhead is counted separately below.
- Pending: final all-development models and sealed external predictions, final manuscript results and a clean-environment full source/model replay. The fixed-eight repeated candidate/representation evaluation is now complete. Local outer-fit audit checkpoints exist; these are not final deployment models or external validation results.
- Local environment: Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, openpyxl 3.1.5, RDKit 2026.9.1. Expression mmap files are two float64 arrays of 1,019 × 57,820, about 471 MB each. No GPU is needed for phase 1. Pilot/formal/diagnostic runtime and sampled memory measurements are reported with their execution boundaries; they are not isolated algorithm benchmarks.

## Quick inspection without original data

```sh
git clone https://github.com/RascalZ/REO.git
cd REO
python scripts/verify_public_snapshot.py
python scripts/verify_public_ordinal.py
python -m unittest discover -s tests -p "test_public_ordinal_snapshot.py" -v
```

These standard-library checks verify published aggregate hashes, schemas and counts only. They do not read private response/OOF/checkpoint files, fit models, or independently establish predictive validity. CI also runs synthetic integrity contracts after installing modeling dependencies; it does not launch the private full-run audit.

The original cohort quick inspection and syntax compilation were tested on 2026-10-10 in a fresh Git clone using isolated Python standard-library execution (`-I -S`). The later ordinal aggregate verifier and its four tests passed locally; this is not a claim of a new clean-clone model reproduction. Full source replay and a newly installed dependency environment have not yet been verified.

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

Both primary REO EN drug models selected a zero molecular endpoint in all three outer folds of this historical pilot. It supports no claim of reliable REO molecular advantage. The subsequent fixed eight-drug repeated analysis and statistics are complete as described below; CTRP evaluation remains sealed and pending.

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

Windows execution repair: the first formal session stopped after 673 successful
tasks when atomic metadata replacement raised WinError5. All 673 checkpoint
checksums and actual fit memberships passed independent recovery audit; no model
failure was discarded. Live readers can block replacement on this filesystem,
even with a shared-delete handle. Do not read live formal_metadata.json or run
partial validation while fits are running. Progress reads append-only logs and
resource CSV. The sequential launcher waits for a successful training-process
exit before reading metadata for QA/statistics/reporting:

```sh
python scripts/run_formal_stage.py --resume --execution-commit db31a35
python scripts/formal_progress.py
```

This repair leaves frozen numerical sources and the formal scientific lock
unchanged. Original failed-session/controller/postprocess records and the 673
checkpoint SHA recovery snapshot remain local, alongside the resumed session.


## Four-drug candidate-rule and ordinal-representation diagnostic

The completed additive diagnostic is frozen at scientific execution commit
`e9b576d060f8ff59fed7f0946b7f48d5ec94fa9a`. It does not replace formal v1.
See [the frozen protocol](ORDINAL_DIAGNOSTIC_V1_PROTOCOL.md) for
3 outer folds × 1 repeat, 3 inner folds, both selectors, three molecular
representations, the identical-policy tissue reference and 9 full-pipeline nulls.
Independent validation confirmed 840 successful tasks, 22,440 fits and 173,600
OOF predictions. True-label evaluation accounts for 84 tasks, 2,244 fits and
17,360 OOF rows; the nine permuted-label evaluations account for the remainder.
These four outcome-informed cases and nine nulls are exploratory engineering
controls; the null count does not support a conventional 0.05 significance test.
The 10,000 paired group-bootstrap draws condition on recorded OOF predictions,
do not include train-refit uncertainty, and yield unadjusted exploratory intervals.
The four-case macro-average must not be reported as a repeat of the original
fixed eight-drug study.

```sh
python scripts/freeze_ordinal_diagnostic.py --verify-only
python -m unittest discover -s tests -p "test_ordinal*.py" -v
python scripts/run_ordinal_stage.py --execution-commit e9b576d
# Explicitly resume the SAME scientific version after an interruption:
python scripts/run_ordinal_stage.py --execution-commit e9b576d --resume
```

These are the recorded execution commands, not instructions to overwrite the
completed run. Use `--resume` only for an interrupted run with exactly matching
scientific source, inputs and checkpoints. An independent replay needs a separate
output location and its own provenance; never delete completed output to bypass
the existing-run protections.

The launcher applies six numerical-library single-thread environment variables
BEFORE child-process imports as well as the frozen runtime thread limit. On this
Windows host, setting runtime limits alone left substantial startup commit
reservations. Two allocation failures and their original checkpoint hashes are
retained locally. Their completed fit events total 31 (4 + 27); one resource probe
adds one completed fit without prediction or scoring. Thus the 22,440 successful
pipeline fits correspond to 22,472 completed fit events overall, with two further
aborted fits making the attempted-fit lower bound 22,474. Recovery verified
unchanged hashes for 601 unique successful checkpoints; overlapping recovery
snapshots must not be added as independent checkpoints. Scientific source and
the scientific lock remained unchanged. See [the recovery record](ORDINAL_MEMORY_RECOVERY.md).

Successful-task wall time sums to 10,614.268 seconds, failed-attempt wall time to
34.041 seconds, and probe time to 6.050 seconds. Controller sessions sum to
10,762.477 seconds; first start to final completion spans 11,419.146 seconds,
including 656.657 seconds between sessions for investigation, the probe and
recovery. These are different measurement boundaries, not additive estimates of
isolated compute cost. The largest sampled worker RSS was 854,458,368 bytes and
process-tree RSS 1,024,815,104 bytes; shared pages can be counted more than once
in process sums. The audit records these values and their provenance in local
`operational_overhead.json` and `ordinal_validation.json`; public aggregates and
the report expose the relevant summaries.

Do not inspect
live `ordinal_metadata.json` or run the validator concurrently with training.
The launcher waits for training exit before independent validation, frozen
statistics and report rendering, logging each stage separately.

```sh
python scripts/validate_ordinal_diagnostic.py
python src/analyze_ordinal_diagnostic.py
python src/report_ordinal_diagnostic.py
```

The analyzer refuses an incomplete or unverifiable run and existing outputs.
Public summaries contain no sample-level responses, predictions or permutation
donor IDs. Private audit files and bootstrap weights remain under ignored `runs/`.
CTRP access stays forbidden. Aggregate-only old-OOF decomposition is available in
`results/oof_decomposition/`; paired intervals share the original group resamples.

In the completed four-case diagnostic, response-association selection improved
macro Spearman relative to MAD for expression (+0.020439), rank (+0.032689) and
pair (+0.041295), with all three exploratory paired intervals above zero. Under
that rule, pair versus expression was +0.002021 (95% CI −0.010449 to +0.014278),
so neither reliable pair superiority nor equivalence is established. There were
26 zero molecular solutions among 72 true-label molecular outer fits, versus
637 among 648 permuted-label molecular fits. The combined count of 663 excludes
tissue-only fits and must not be presented as a true-label result. Tamoxifen's
three response-associated representations all equaled the tissue reference;
it is a prespecified diagnostic comparator, not proof that leakage is impossible.
Feature stability remains limited, and no mechanistic biomarker claim follows.
See [the stage-5 report](../reports/阶段5_秩与基因对诊断试跑.md),
[scientific assessment](ORDINAL_DIAGNOSTIC_ASSESSMENT.md), and the sole maintained
[research plan](../RESEARCH_PLAN.md) for interpretation and the next decision.


## Fixed-eight repeated candidate/representation evaluation

The newly authorized run uses `config/ordinal_repeated_v1.json` and
`docs/ORDINAL_REPEATED_V1_PROTOCOL.md`, preserving both historical scientific
locks. It compares the six 200-gene rule/representation arms, same-policy tissue
means and MAD2000 expression EN on shared 5-outer × 3-repeat / 3-inner splits.
The expected observed-label contract is 960 tasks, 26,160 fits and 119,208 OOF
rows from 626 global biological groups and 4,967 drug-model observations.
There are no new permutation fits in this version; the completed four-case nine
null controls remain engineering evidence. CTRP response/index access is forbidden.

Single-task child processes exit after each outer task. Numerical thread limits
are set before imports. The controller checks available physical RAM (512 MiB)
and commit (1 GiB), polling every five seconds for up to 300 seconds before a
checkpoint-preserving `resource_paused` exit 75. Repeating a stopped run requires
explicit same-lock resume; no score-dependent arm or feature changes are allowed.
Resource observations include preparation/process startup/waits and must not be
interpreted as isolated representation efficiency benchmarks.

The primary Spearman comparison is response-associated pair versus strong
expression EN, with eight drugs plus equal-drug macro as a nominal family of nine.
Within-rule pair versus expression200 comparisons form a separate nominal family
of eighteen. Conditional OOF group-bootstrap intervals do not cover complete
train-refit uncertainty. Noninferiority/equivalence remains a pending scientific
question: no numerical margin, including −0.02, is adopted, and no verdict is made.

Implementation and freeze are complete. Execution began on 2026-10-13 at
18:14 Asia/Shanghai using scientific commit `5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43`.
Training completed on 2026-10-13 at 22:48 Asia/Shanghai. Independent QA, statistics and reporting then passed; numerical conclusions and limits are recorded below.

Pre-fit acceptance on 2026-10-13 passed 19 synthetic integrity and six statistical
contracts, all entry points compiled, and both historical locks verified unchanged.
The new repeated lock SHA256 is `ef8f588b99c42cb54db73b495127bc41882fa3e5b450e62167c6b46c53d583e7`.
Each fresh worker binds source/config/protocol hashes before fitting and uses
round-trip response parsing. Recorded commands for this separate run are:

```sh
python scripts/freeze_ordinal_repeated.py --verify-only
python scripts/run_ordinal_repeated_stage.py --execution-commit 5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43
# Only if an interrupted run has stopped:
python scripts/run_ordinal_repeated_stage.py --execution-commit 5e54a693fe3a047c53b3d9e2dfeb7638ff1b9b43 --resume
```

Progress reads append-only logs with `python scripts/ordinal_repeated_progress.py`.
Do not audit live metadata. Actual scientific completion requires the post-exit
validator and statistics; a preflight test pass is not a prediction result.


The completed repeated run has 960 successful tasks, 26,160 fits and 119,208 OOF
rows. Independent QA verified 960 fresh process identities, 1,440 shared-gene
contexts, the original/source float labels and exact aggregate outputs; warning
and fallback fits are zero. There were no failed/paused sessions or extra probes.
Controller wall time was 16,431.152 seconds, successful task wall sum 14,340.372
seconds, and child wall sum 16,312.905 seconds. These overlapping timings cannot
be added. The process-tree sampled peak RSS was 1,253,429,248 bytes; private peak
was 1,114,628,096 bytes. Task samples exclude interpreter/import startup, include
worker initialization, and are operational observations rather than an isolated
complexity benchmark. Resource-policy check/wait timing includes immediate
checks, not only actual pressure-induced waiting.

Aggregate-only verification is now available without original inputs:

```sh
python scripts/verify_public_ordinal_repeated.py
python -m unittest discover -s tests -p "test_public_ordinal_repeated_snapshot.py" -v
```

Fourteen products (thirteen CSVs plus operational_overhead.json) passed hashes
and coverage/algebra checks. The paired table has 900 effects from 25 contrasts
and four metrics, with separate Spearman families of nine and eighteen. Main
association-pair versus strong-expression macro difference is -0.005882 with
family-nine nominal conditional interval [-0.018714, 0.006614]; it does not
establish superiority, noninferiority, equivalence or full information retention.
The corresponding RMSE improvement is negative. Full results and methodological
limits are in the stage-6 report and ORDINAL_REPEATED_ASSESSMENT.md. No new null
fits were performed; four-case historical nulls retain their engineering scope.
CTRP remains sealed and final models have not yet been trained.
