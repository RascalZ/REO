# Scientific revision R1 evidence package — 2026-10-16

The authoritative content-review manuscript is `manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md`, with its same-source illustrated DOCX. Supplementary Methods/Results and Additional files 1–14 are in the corresponding `scientific_revision_v1` paths. Earlier manuscripts are preserved. This package is not submission-ready; author facts, formal archive/rights confirmation and remaining reproducibility-access work are not inferred or declared complete.

## Scope

This revision transcribes existing public aggregates and existing synthetic controls. It does not fit models, read individual CTRPv2 responses or predictions, open model payloads/bootstrap weights, resample, perform new inference or rerun external scoring. Source-to-text checks are not independent reproduction of the scientific analysis.

## Records

- `SYNTHETIC_REPORTABILITY.md` and `synthetic_reportability_audit.json`: coverage, equations, source/execution hashes, all-scenario diagnostics and caveats.
- `methods_source_ledger.csv`: 34 source-bound Methods facts, with hashes. The Methods review additionally identifies the frozen gateway sources used for SM4.1.
- `INDEPENDENT_RESULTS_REVIEW.md`, `INDEPENDENT_METHODS_REVIEW.md`, `INDEPENDENT_SCIENTIFIC_REVIEW.md`: role-separated second-pass checks within the same AI-assisted workflow; these are not statistically independent reviews or raw-data replication.
- `REVIEW_RESOLUTION.md`: integration decisions and resolution of the reported issues.
- `manuscript_build_manifest.json`, `citation_validation.json`, `revision_validation.json`, and `document_qa.json`: editorial transformation, reference coverage, package/structure checks and final Word-render provenance.
- `staged_byte_validation.json`: 47 hash-bound source/output files matched the Git index byte-for-byte before commit.
- `git_eol_storage_audit.json`: scoped Git storage repair for 16 existing aggregate/supplement files; differences are CRLF/LF only, preserving the originally recorded source bytes without data/statistical edits. New hash-bound artifacts also opt out of text normalization.
- Table provenance: `tables/scientific_revision_v1/revision_table_manifest.json` and `synthetic_table_manifest.json`.
- Distribution manifest: `manuscript_v1/bmc_additional_files_scientific_revision_v1/package_manifest.json`.

`ADDED_RESULTS_DRAFT.md`, `METHODS_AMENDMENTS_DRAFT.md`, `synthetic_methods_draft.md` and `synthetic_results_draft.md` are retained integration inputs. They are not additional authoritative manuscripts; numbering, rounding, final wording and package status are governed by R1 and its manifests.

## Rebuild

From the repository root, with the existing Node.js/npm document dependencies:

```sh
python src/build_scientific_revision_tables_v1.py
python src/build_scientific_revision_manuscript_v1.py
python src/build_scientific_revision_package_v1.py
python src/validate_numeric_references.py --manuscript manuscript_v1/Manuscript_BMC_Genomics_scientific_revision_v1.md --output docs/scientific_revision_v1/citation_validation.json
python src/validate_scientific_revision_v1.py
```

The table builder reads four explicitly allowed public aggregate files. The manuscript builder applies a deterministic editorial transform to the preserved preceding BMC draft. DOCX core timestamps can change on rebuilding; refresh the package and render audit after an intended rebuild. Do not run these editorial scripts as a substitute for model reproduction or invoke an external-scoring entry point.
