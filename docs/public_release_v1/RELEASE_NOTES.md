# REO v0.1.0 — release notes

**Prepared release; publication and archive verification are pending.** This
document is not evidence that a GitHub Release or a Zenodo DOI already exists.
After publication, cite the verified software version and archive identifier;
a software archive DOI is not the unpublished manuscript's DOI.

This version packages the existing research software and frozen model exports
for an auditable comparison of expression, within-sample ranks and
relative-expression gene pairs, separating tissue-reference performance,
molecular increments and representation differences under matched candidate
rules. The study did not establish an average predictive advantage for REO.
Negative, undefined and zero-molecular-solution results remain part of the
record.

## Included assets

- 64 portable numeric JSON.gz models: eight fixed drugs by eight methods,
  including eight tissue references and 56 molecular models. All model payloads
  are byte-identical to the existing export. The 29 zero and 27 nonzero final
  molecular solutions are retained.
- Candidate genes, coefficients, tissue parameters, complete ordered gene and
  tissue schemas, and a generated full-gene input for response-free inference.
- Minimal frozen identity/membership projections for 626 development and 175
  external groups (4,967 and 1,185 drug memberships), with four original fold
  assignment files preserved byte-for-byte.
- Per-file SHA-256 inventories, archive checksums, projection audits and model
  preservation evidence. The published asset catalog records their exact scope.

The minimal metadata deliberately omits source experiment/curve/response
pointers and is not a complete input for historical external rescoring.
Original upstream data must be obtained independently under their own terms.
The frozen-model inference example uses generated expression numbers and no
response labels. It is not an efficacy evaluation or a new performance result.
No training, selection or scientific scoring was performed to build these
packages. Individual CTRPv2 responses, response indices and real individual
external predictions remain sealed and are not included.

## Authorship and licensing

Release creators follow the author-confirmed v4 order: Nan Zhang, Huicheng
Zhang, Wen Huang, Siru Yang, Yutong Zhang, Jinglin Li, Hui Liu and Shuyuan Wang.
Their confirmed affiliation is the College of Bioinformatics Science and
Technology, Harbin Medical University, Harbin 150081, China. No ORCID has been
invented.

Newly authored repository code uses MIT. Original model parameters, derived
tables, project fold assignments, author-created identity decisions, generated
examples and asset documentation use CC BY 4.0 only to the extent the authors
own the applicable rights. This grant does not relicense third-party data,
identifiers or annotations. Read each package's `MODEL_ASSETS_LICENSE.md` and
`THIRD_PARTY_ATTRIBUTION.md`; source-specific rights and attribution remain.

The manuscript and author submission package are not part of this software
release. Publication of these research assets does not imply all-author
approval of a journal submission.
