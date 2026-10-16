# Added results for scientific revision R1 — author-side working draft

Generated only from the four public aggregate sources listed in the table manifest. No model fit, sample-level response/prediction access, bootstrap draw/weight access, resampling, new interval, or new inferential family was performed.

## Proposed Results addition: external absolute performance

Absolute ranking associations in the locked CTRPv2 evaluation were descriptively non-zero. The strict equal-eight-drug mean Spearman correlation was 0.228400 for the tissue-only reference, 0.260238 for the strong-expression model, and 0.259866 for the response-associated pair model (Table S9). These descriptive macro points are strict arithmetic means of the eight existing drug-specific correlations; no absolute-macro intervals were precomputed. A near-zero pair-minus-expression contrast therefore does not imply absence of all transferable predictive association. The three drugs with identical zero-molecular primary models retain tissue-based predictions and provide no pair-specific molecular comparison.

This is an association/ranking statement, not a test of causal molecular information. Existing absolute-rho conditional intervals are transcribed for the five molecular arms represented in the original inferential comparisons (40 drug–arm rows). Tissue-only, rank-arm, and absolute-macro intervals were not precomputed and remain NA. No interval was inferred from equal predictions or from subtracting interval endpoints.

## Proposed Results addition: complete matched-candidate representation contrasts

Table S10 makes the prespecified representation comparisons explicit for both development and external evaluation. Only pair-minus-expression carries the existing 18-estimand family adjustment; development rank contrasts retain their existing unadjusted descriptive intervals, and external rank contrasts remain descriptive point estimates without new intervals.

| Stage | Candidate rule | Contrast | Δρ | Existing nominal conditional 95% interval | Existing family-adjusted interval |
| --- | --- | --- | ---: | --- | --- |
| Repeated_nested_development | mad | rank − expression | -0.012439 | [-0.018786, -0.006221] | NA (not precomputed) |
| Repeated_nested_development | mad | pair_full − expression | -0.017645 | [-0.026595, -0.008907] | [-0.031170, -0.004064] |
| Repeated_nested_development | mad | pair_full − rank | -0.005206 | [-0.012215, +0.001811] | NA (not precomputed) |
| Repeated_nested_development | response_association | rank − expression | +0.004406 | [-0.001504, +0.010299] | NA (not precomputed) |
| Repeated_nested_development | response_association | pair_full − expression | -0.002252 | [-0.010382, +0.005824] | [-0.014384, +0.009491] |
| Repeated_nested_development | response_association | pair_full − rank | -0.006659 | [-0.014560, +0.001026] | NA (not precomputed) |
| Locked_CTRPv2 | mad | rank − expression | +0.016579 | NA (not precomputed) | NA (not precomputed) |
| Locked_CTRPv2 | mad | pair_full − expression | +0.007822 | [-0.009541, +0.025451] | [-0.018837, +0.034499] |
| Locked_CTRPv2 | mad | pair_full − rank | -0.008757 | NA (not precomputed) | NA (not precomputed) |
| Locked_CTRPv2 | response_association | rank − expression | -0.008715 | NA (not precomputed) | NA (not precomputed) |
| Locked_CTRPv2 | response_association | pair_full − expression | -0.006712 | [-0.019105, +0.006186] | [-0.025624, +0.013687] |
| Locked_CTRPv2 | response_association | pair_full − rank | +0.002003 | NA (not precomputed) | NA (not precomputed) |

The rank-arm comparisons do not yield a general REO advantage. In development, the response-associated rank-minus-expression macro interval includes zero; the response-associated pair-minus-rank macro point estimate is negative and its existing descriptive interval also includes zero. External rank contrasts have no precomputed uncertainty intervals and cannot support a new significance conclusion. These results neither establish equivalence/non-inferiority nor identify causal shares of amplitude loss, feature expansion, collinearity, or elastic-net regularization.

## Table interpretation and boundaries

- Table S9 has 72 rows: 64 drug–arm points and eight strict equal-drug macro points. N_unique_external_groups for a macro is the fixed external universe of 175, not 8 independent subjects; the 1,185 drug–group response records are a different denominator.
- Table S10 has 108 rows: two evaluation stages × two candidate rules × three representation contrasts × nine drug/macro estimands. Spearman only; positive Δρ favors the candidate representation.
- NA is retained for unavailable/undefined quantities. Tissue-only coefficient-zero status is not treated as a fitted molecular zero solution. Macro coefficient counts are not summed as if they were one model.
- Development intervals condition on saved OOF predictions and omit training/tuning/refitting uncertainty; external intervals condition on fixed fitted models and observed external groups. The two study families remain separately adjusted and do not guarantee whole-study simultaneous coverage.
- No performance-retention margin was chosen and no equivalence or non-inferiority test was added. The negative/near-zero primary REO advantage conclusion remains unchanged.
