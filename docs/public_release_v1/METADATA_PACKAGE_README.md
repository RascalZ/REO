# REO v0.1.0 minimal identity and assignment metadata

The four identity/membership CSVs are new column projections of the existing
prepared metadata. They preserve all rows, row order and retained field strings.
They are derived views, not byte-identical replacements for the old full
metadata exports. Identity tables retain canonical model ID, biological group
ID, CCLE ID, prepared expression-row offset, frozen tissue and cohort role;
cohort version and project eligibility status are retained as context. Drug
membership tables additionally retain GDSC drug ID and drug name.

The four fold CSVs are unchanged byte copies of the existing project
assignments. Their original hashes remain binding. These files contain no
source experiment/curve/response pointers and are not newly generated splits.
Use the matching repository's frozen fold definitions to interpret repeated,
final-fit and leave-one-tissue-out assignments. Fold rows are assignments, not
independent observations.

The projection omits CTRPv2 compound/master/experiment/source-row pointers,
experiment counts, external QC/subtype records, GDSC curve/source-row/dose
fields and source-evidence strings. It contains no response values or real
individual predictions. This is **not a complete input for the old external
runner's rescoring or raw-response extraction**: that runner requires the
original source pointers and aggregation/QC contract. Original upstream inputs
must be obtained independently under their own terms.

`expression_row_index` is a zero-based offset into the project's historical
prepared expression matrix. It is not a portable DepMap row number. Preserve
the frozen biological groups and tissues; do not infer them anew. DepMap 24Q4
alone cannot reconstruct every frozen historical external identity. See
THIRD_PARTY_ATTRIBUTION.md. The separate frozen-model replay accepts explicit
identity/tissue strings and does not read these CSV files.

`mapping_schema.json` describes the published columns. `ASSET_MANIFEST.json`
hashes all other files in this archive. The outer release file manifest also
hashes this manifest. Third-party identifiers are not relicensed by
MODEL_ASSETS_LICENSE.md; that notice grants only rights the authors own.
