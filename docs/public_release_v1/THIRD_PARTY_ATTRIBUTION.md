# Third-party attribution and upstream inputs

This package contains author-produced frozen model state, parameter summaries,
project assignments and a minimal identity projection. It does not distribute
the upstream expression or drug-response matrices. Obtain upstream inputs
independently under their own terms. Public availability of this release does
not change those terms.

## DepMap identity metadata

DepMap 24Q4 Public, version 1, Broad Institute / DepMap:
https://doi.org/10.25452/figshare.plus.27993248.v1 .
Official release: https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248/1 .
Model.csv: https://ndownloader.figshare.com/files/51065297 .
The retained official metadata snapshot records CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/), Model.csv size 645696 bytes and
MD5 675210d17675f3517b0ce39a3c274f16. The local source file matched that size
and MD5. Its SHA-256 is
`b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa`.

Canonical model and CCLE identifiers are third-party identifiers; their
inclusion does not relicense them. Project biological groups, frozen tissue
decisions, cohort membership and prepared expression-row offsets must retain
their existing interpretation. They cannot all be reconstructed from this
DepMap release alone: all 626 development identities have exact ModelID and
CCLEName matches; among 175 external identities, 174 ModelIDs are present and
172 of those have exact CCLEName matches. Frozen identities were not replaced
or rematched for this release.

## Additional identity evidence

The historical eligibility/identity audit also retained Cellosaurus evidence,
SIB Swiss Institute of Bioinformatics, with its recorded CC BY 4.0 attribution:
https://www.cellosaurus.org/ and https://www.cellosaurus.org/description.html .
Historical Sanger and KCLB identity evidence is described in the repository's
source records. The minimal release omits their experimental/source pointers
and does not purport to relicense those sources.

## GDSC, CTRPv2 and historical CCLE expression

GDSC drug identifiers/names and historical cohort membership identify the
frozen project panel. CTRPv2 response values, source row indexes, experiment
identifiers, QC records and individual predictions are not included. GDSC
source curve IDs, source row indexes, concentration windows and response values
are not included. Historical CCLE expression observations are not included;
the packaged complete-gene input is explicitly generated numerical data.
Model parameters and candidate-selection summaries are fitted/derived state,
not raw source observations. No new rights in upstream datasets or identifiers
are granted by MODEL_ASSETS_LICENSE.md.
