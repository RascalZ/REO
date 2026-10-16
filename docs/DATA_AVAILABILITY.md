# Data availability and provenance

The repository contains newly authored curation, modeling, audit and figure/table code, public aggregate results and manuscript drafts. Large expression/response matrices, raw downloaded metadata, old manuscripts, individual local paths and original-file inventories are not committed. Source terms continue to apply; the repository MIT license does not relicense upstream data.

## Current publication status — 2026-10-16

Development, final model fitting and the single retrospective CTRPv2 evaluation are complete; do not repeat external scoring as part of publication preparation. The source descriptions below preserve the historical cohort-v1 input provenance. The paper uses cancer-model cohort v2 with 626 development and 175 external groups, and 19,836 uniquely mapped protein-coding genes from the GENCODE19-based audit. Keeping unmatched records in the source audit does not mean they were all included in modeling.

The current source/model accessibility gaps are recorded in stage 22 and `docs/reviews/bmc_pre_submission_v1/issue_register.csv` (M04/M06). In particular, the legacy compatibility dependencies below have not been silently replaced, and no complete public-input model replay or public availability of all trained models is claimed. A paper release must supply precise source/member/version/hash information and an actual artifact-access route or justified restrictions. The MIT license is already selected for new code; exact upstream release terms remain separate. An archived release needs a real persistent identifier; a DOI is one recommended choice, not a fabricated or universally exclusive requirement.

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
