# Internal evidence map and editorial checklist for the BMC Genomics draft

## Evidence source map for this draft

| Draft section | Primary source |
| --- | --- |
| Study design, data, and nested evaluation | `docs/ORDINAL_REPEATED_V1_PROTOCOL.md`; `config/ordinal_repeated_v1.json`; `src/ordinal_diagnostic_models.py`; `src/ordinal_repeated_models.py`; `src/modeling/nested_cv.py`; `src/run_ordinal_repeated.py` |
| Data provenance and cohort sizes | `docs/DATA_AVAILABILITY.md`; `results/cohort_v2/public_cohort_snapshot.json`; `tables/manuscript_v1/Table1_cohort_by_drug.csv`; `tables/manuscript_v1/TableS8_compound_identity_mapping.csv` |
| Development estimates and zero solutions | `reports/阶段6_固定八药表示与候选规则重复评价.md`; `tables/manuscript_v1/Table2_primary_effects.csv`; `TableS3_development_27_estimands.csv`; `TableS4_development_model_performance.csv`; `TableS6_zero_solution_audit.csv` |
| Leave-one-tissue-out results | `docs/TISSUE_SENSITIVITY_LOTO_PROTOCOL_V1.md`; `reports/阶段10_LOTO组织外推敏感性.md`; `tables/manuscript_v1/TableS1_loto_defined_cells.csv`; `TableS7_loto_all_23_tissues.csv` |
| Locked CTRPv2 analysis | `docs/EXTERNAL_VALIDATION_PROTOCOL_V1_1.md`; `reports/阶段14_P6一次性CTRP外验与研究方向判定.md`; `results/p6_external/evaluation_summary.json`; `tables/manuscript_v1/TableS2_external_27_estimands.csv` |
| Citation metadata and claim support | `references.bib`; `reference_dois.txt`; `reports/阶段18_引文元数据与主张核验.md` |

## Author/editorial checks still required

1. Verify all text against the linked protocol, code and aggregate table, especially units, the response direction, the meaning of independent biological groups, and each interval's conditioning scope.
2. Confirm the author-preferred reference style and retain the claim-to-source audit; bibliography metadata and DOIs have been reconciled against Crossref plus publisher/PubMed records, with article identifiers used where journals publish article numbers rather than page ranges.
3. Add a complete data/code availability statement, software versions and release/archive identifiers; review upstream licenses and metadata redistribution conditions.
4. Replace this working file with the journal-formatted manuscript only after author details, funding, conflicts, contributions and journal-fit review are confirmed. Do not edit or overwrite historical submissions.
