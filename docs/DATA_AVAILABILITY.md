# Data availability and provenance

## Current access — 2026-10-27

`RascalZ/REO` is now **public**. Its tracked code, frozen configurations,
aggregate results, source/provenance records and invented-input example are
accessible anonymously. This supersedes the historical private-repository
status below; it does not change upstream dataset terms.

Two v0.1.0 packages have been prepared and hash-checked: 64 byte-preserved
numeric models with parameter/feature tables and generated input, and minimal
identity/membership projections with four unchanged project fold files. The
[planned GitHub Release](https://github.com/RascalZ/REO/releases/tag/v0.1.0)
and actual package downloads are **pending publication and download verification**.
See the [asset guide](../reproducibility_v1/README.md),
[prepared-byte summary](public_release_v1/PAYLOAD_SUMMARY.json) and
[per-file catalog](public_release_v1/PAYLOAD_ASSET_CATALOG.json).
Zenodo has genuinely reserved `10.5281/zenodo.23681422` for the software archive,
but it is **not yet registered/published**. This is not a paper DOI and does not
yet establish an accessible archived record. Repository release immutability
is enabled for future releases; no published immutable release is asserted here.

The metadata projection preserves all 626/175 frozen development/external
identities and 4,967/1,185 drug memberships but excludes source experiment,
curve, response-index, QC and dose-pointer fields. It is not a complete input
for historical external extraction or rescoring. DepMap24Q4 alone cannot
reconstruct every frozen historical identity. Original expression/response
matrices and wider legacy curation dependencies must be obtained independently
under their source terms; a public checkout and these minimal packages alone
do not reproduce every historical curation audit or scientific result.

Newly authored code/repository documentation use MIT. The separate asset
CC BY 4.0 grant applies only to rights the authors own; it does not relicense
third-party identifiers, annotations or data. See the
[asset notice](public_release_v1/MODEL_ASSETS_LICENSE.md) and
[upstream attribution](public_release_v1/THIRD_PARTY_ATTRIBUTION.md).
Individual CTRPv2 responses, response indices and real individual external
predictions remain sealed; release preparation performs no new training,
resampling or scientific scoring. The local v4 manuscript/submission package
is preserved and is separate from the software release.

## Historical access correction — 2026-10-20

At that date, authenticated GitHub repository metadata confirmed that `RascalZ/REO` was **private**, and the anonymous endpoint was inaccessible. Tracked code, aggregates, source records and invented-input examples were then private. Their presence and older descriptions as “public” did not establish public or anonymous editor/reviewer access. This paragraph records the historical correction; the current access section above now governs.

Publicly accessible GDSC2, CCLE, CTRPv2, DepMap and Cell Model Passports resources have their own terms; public source access does not itself authorize redistribution. At that stage complete derived models/mappings remained local and Git-ignored, and publication/archive arrangements were pending. The 20 October editorial step itself changed no visibility or asset access.

The current main manuscript is in local Git-ignored `manuscript_v1/author_review_v4/`, with its paired Chinese copy in `manuscript_v1/author_review_v4_zh/`. The two clarifications in [Stage 31](../reports/阶段31_Minor1与2最小补充.md) change no analysis, rights or access permission. The unchanged supplementary files are now self-contained in `manuscript_v1/author_review_v4/Additional_files/` and `Supporting_sources/`, copied byte-for-byte before the old v1/v2 reading folders were archived and cleaned. Prior non-QA materials are recoverable from the local archive; v3 remains in place. See [Stage 32](../reports/阶段32_仓库规范投稿打包与旧副本清理.md) and [history restoration](repository_release_v1/LOCAL_HISTORY.md). The Stage 29 editorial revision changes no analysis or access permission. The 21 October 2026 question-led narrative revision changes no analysis or access permission. Source-recorded author/declaration facts were adopted, but absent contact/contribution facts, formal publication permission and approval by all authors of the revised manuscript are not fabricated. See [Stage 26](../reports/阶段26_作者资料迁入与最终通读稿.md) and plan v3.46; the narrative revision is recorded in Stage 27. User read-through and discussion are the next gate; current-package rendering/check results are recorded only after verification. The original M08 manuscript and scientific files are preserved.

The unique locked CTRPv2 evaluation is already complete. Saved aggregate results may continue to support the manuscript; individual responses, response indices and individual predictions remain sealed during this editorial work. No new label access, fitting, resampling or scoring is performed.

The repository contains newly authored curation, modeling, audit and figure/table code, public aggregate results and manuscript drafts. Large expression/response matrices, raw downloaded metadata, old manuscripts, individual local paths and original-file inventories are not committed. Source terms continue to apply; the repository MIT license does not relicense upstream data.

## Historical publication status — 2026-10-16

Development, final model fitting and the single retrospective CTRPv2 evaluation are complete; do not repeat external scoring as part of publication preparation. The source descriptions below preserve the historical cohort-v1 input provenance. The paper uses cancer-model cohort v2 with 626 development and 175 external groups, and 19,836 uniquely mapped protein-coding genes from the GENCODE19-based audit. Keeping unmatched records in the source audit does not mean they were all included in modeling.

Source/model accessibility gaps at that stage were recorded in stage 22 and `docs/reviews/bmc_pre_submission_v1/issue_register.csv` (M04/M06). Legacy compatibility dependencies were not silently replaced, and that stage did not establish a complete public-input model replay or public model access. A release needs precise source/member/version/hash information and actual artifact-access routes or justified restrictions. Source terms remain separate from MIT for new code. Current release progress is recorded above rather than retroactively claimed by this stage.

## Historical prepared derived assets — stage 24

The [input/identity ledger](reproducibility_assets_v1/INPUTS_AND_IDENTITIES.md) now inventories 69 current and historical dependencies with recorded versions/hashes. Local, Git-ignored technical assets preserve development/external metadata projections (626/175 identities; 4,967/1,185 drug memberships), four original split CSVs, 64 full numeric models and complete candidate/coefficient/tissue tables. No source CTRPv2 label/index or individual external prediction was read, and no model was retrained/scored. Metadata projections were taken only from positively whitelisted response-free manifests.

The guide, format definition, hand-invented fixture and hash/count summaries were committed at that stage. That step did not release actual fitted vectors, full mappings or original expression/response matrices. It prepared technical assets without establishing public download access. Current release packaging and access are described above; the MIT license has never implied permission to redistribute all upstream material.

The final numeric worker uses frozen development responses, expression and fold assignments, without direct legacy-table reads. Original controller verification and source-data curation retain wider metadata/legacy dependencies. Those have not been silently removed, replaced or rerun. Numerical inference readiness is therefore distinct from complete public-data-only historical reconstruction.

## Inputs used for cohort v1

| Input | Frozen source | Availability and use |
|---|---|---|
| GDSC2 fitted response | `GDSC2_fitted_dose_response_27Oct23.xlsx` | Existing local file; [GDSC download portal](https://www.cancerrxgene.org/downloads/bulk_download). SHA256 `2f9e5aab648a6785ba4a915f04d254a1443c891cbb0f376d1a25855d399bf864`. Retrieve exactly this version for reproduction. |
| CCLE gene TPM | RSEM gene TPM 20180929 | [Broad historical directory](https://data.broadinstitute.org/ccle/). Local uncompressed SHA256 `1d46dd8d77f635b7716be75c503623229efc761982247022b37ee749ce906dbd`. Internal checks passed; remote byte identity was not verified. |
| CCLE counts/annotations | 20180929 / 20181226 | Historical source as above. Counts support the input audit; the primary REO matrix uses TPM, preserving cross-gene length normalization. |
| DepMap Model + README | 24Q4 Public v1 | [Versioned Figshare DOI](https://doi.org/10.25452/figshare.plus.27993248.v1), CC BY 4.0, official MD5 verified. Identity/related tissue-sample metadata only; no 24Q4 expression is mixed into the 2018 matrix. |
| Cell Model Passports | API snapshot 2026-10-10 | [Official API documentation](https://depmap.sanger.ac.uk/documentation/api/). Public non-commercial research terms; retain attribution and check redistribution terms before releasing full raw metadata. |
| GDSC compound annotations | release 8.5 | [Official screened compound table](https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/screened_compounds_rel_8.5.csv). Contains official names/IDs/targets; not native measured chemical structures. |
| CTRPv2 expanded curves and metadata | 2015 expanded archive | [Archive URL](https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset.zip), found through [curator configuration](https://github.com/BHKLAB-DataProcessing/CTRPv2-Pharmacoset_Snakemake). Required members were range-downloaded and checked against ZIP CRC and supplied MD5 where present. The entire archive/raw-well data were not downloaded. |
| PubChem identity candidates | cached queries 2026-10-10 | [PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest). Full Standard InChIKey support from official name resolution versus CTRP SMILES; does not chemically authenticate assay bottles. |
| Gene annotation | matching GENCODE19 GTF | [Broad historical directory](https://data.broadinstitute.org/ccle/). Unmatched genes are retained and not automatically called noncoding. |

Detailed local URLs, times, hashes and status are recorded under `data_raw/*/*manifest*.json`, `data_processed/expression_audit.json` and [download_register.csv](../data_manifest/download_register.csv). The complete freeze and individual source rows remain in the local workspace; `results/public_cohort_snapshot.json` records their validation hash without publishing machine paths.

## Original project compatibility inputs

The current identity audit additionally reads the legacy annotation/CTRP name bridge and legacy four-column CTRP table to audit historical identity coverage. These are original project inputs and are not silently synthesized from another source. Configure exact files using `config/input_paths.local.json`; see the example. A third-party user without these files can inspect aggregate outputs and code, but cannot yet reproduce every historical audit step. A clean public-data-only alternative is a pending task and must preserve/revalidate cohort identity before it is used.

## Historical optional data-work discussion

At the early timing-pilot phase no further download was required for that pilot. This paragraph is historical and does not resolve the current publication-access gaps above. Common-dose-window re-integration needs a separate feasibility assessment of fit parameters and possibly dose-level files. Independent expression measurements/PRISM or GEO are optional extensions with separate selection and licensing checks; they are not evidence already produced by this stage. Clinical patient data and full FASTQ reprocessing are outside this month's primary analysis.
