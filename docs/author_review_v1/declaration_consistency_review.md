# Author-review copy: declaration consistency review

Review date: 20 October 2026. Scope: read-only comparison of the final Human Genomics submission manuscript with the current M08 manuscript and local reproducibility documentation. No model training, scoring, or access to CTRPv2 individual responses, response indices, or predictions was performed. This report deliberately omits author identities, institutional addresses, and contact details.

## Sources and limits

- Historical source: the user-designated `12-Human Genomics/Manuscript.docx`, SHA-256 `bdfaa1fee9ad68f5fc47e1c7d109f1fda42721ab1dd333fa6b767f903ed47add`.
- The source has 152 top-level body paragraphs (zero-based indices 0–151). All Word XML parts were additionally searched for contact, contribution, ethics/consent, availability, AI, permission, and approval statements. No such additional declaration was found in headers, footers, or tables. A search not finding a statement is not evidence that the underlying fact is false.
- Current comparison: `manuscript_v1/Manuscript_BMC_Genomics_M08_v1.md`, especially title-page placeholders, Methods “Reproducibility and data governance,” and Declarations; `reproducibility_v1/README.md`; repository `README.md`.
- Repository visibility update supplied by the coordinating agent: the authenticated GitHub connector `list_repositories` reports `RascalZ/REO` as private; the anonymous repository endpoint returns 404. Accordingly, repository tracking or synchronization must not be described as verified public access. The anonymous 404 alone would not establish the reason for inaccessibility; the authenticated visibility result supplies that distinction.
- User instruction authorizes direct adoption of recorded author/declaration facts into the local author-reading copy. It does not establish facts absent from that source and does not mean all listed authors have approved the substantially revised manuscript for submission.

## Carry forward from the historical source

| Item | Source locator | Recommendation |
|---|---|---|
| Author sequence, affiliation, equal-first-author notation, and correspondence markers | Body paragraphs 1–4 | Preserve exactly, as instructed. The recorded title-page roster has eight authors; one affiliation; first two marked equal contribution; final two marked corresponding authors. Do not infer extra author roles from these markers. |
| Competing interests | Paragraphs 120–121 | Carry forward the recorded no-competing-interest declaration; it may be condensed to BMC wording without broadening the claim. |
| Funding | Paragraph 122 | Preserve both named funding sources and grant identifiers exactly. Correct the missing space between the heading and sentence. The source does not describe funder roles; do not add “the funders had no role” as an established fact. |

The author roster in this specific user-designated final manuscript is the source of authority for this copying operation; do not reconstruct it from the earlier submission-history note.

## Facts the historical source does not supply

| Item | Evidence status | Treatment in the author-reading copy |
|---|---|---|
| Corresponding email addresses | Not present in the inspected DOCX text | Do not invent or obtain by matching people on the internet. Use a concise separate author-review note, rather than guessed addresses. |
| Individual author roles / CRediT | Only equal-contribution notation is present | Preserve the equal-contribution fact; specific conceptualization, coding, analysis, funding-acquisition, supervision, and writing roles remain unverified. Do not allocate these by author order. |
| Final approval of the reconstructed manuscript | Not established; the user is requesting a reading copy | Do not state that all authors have read and approved this new final manuscript. Describe the deliverable as a final copy for author review, not as submission-authorized. |
| Ethics committee approval or formal exemption identifier | No statement found | Do not fabricate an approval, waiver, committee, or identifier. A factual description of secondary public-cell-line data analysis, with no new recruitment or material collection, can be based on the current study design. It must not purport to be a documented institutional exemption. |
| Consent for publication | No statement found | “Not applicable” may be drafted from the current manuscript’s absence of identifiable individual-person case data; this is a present-design assessment, not a copied historical declaration or an all-author approval statement. |
| Acknowledgements and permissions to name contributors | No statement found | Do not invent acknowledged individuals, their permissions, or a confirmed absence of contributions. List the unresolved item in the accompanying review note. |
| Specific upstream redistribution permission or institutional release restrictions | No statement found | Public resource use does not establish unrestricted redistribution. Do not silently remove existing rights/access limitations. |
| Archived code release / dataset DOI | No study-specific archived release identifier found | Do not fabricate DOI, release tag, public archive, or an immutable accession. Upstream article DOIs in the bibliography are not a study-code archival DOI. |

## Availability, inference assets, and model preprocessing

Scientific code, aggregate results, provenance, and an invented-input response-free inference example are tracked in the version-controlled project repository. The repository is private, so none of those tracked assets has verified public accessibility through that URL; editor/reviewer access arrangements also remain to be established. Actual fitted models, mappings, sample-level responses/predictions, and bootstrap objects are not publicly released. The later M04 step prepared local, label-free reproducibility assets. Local preparation, Git tracking, synchronization, and numerical implementation checks are not public availability or independent scientific replication. The current availability wording must distinguish these states instead of describing the project repository or its example as public.

Historical claims of broad robustness, local retraining, clinical transfer, and the previous analysis pipeline must not replace current frozen Methods. Current preprocessing, candidate selection, identity grouping, elastic-net tissue adjustment, zero-coefficient semantics, and abstention/missingness rules belong to the reconstructed analysis and remain unchanged. A title-page/declaration merge is not an authorization to import old scientific claims or run new prediction checks.

Suggested availability scope for an author-reading draft: distinguish the publicly accessible upstream source resources from the private project repository. State that code, aggregate results, provenance, and a synthetic response-free example are maintained in the private repository; complete derived model/mapping assets have been prepared locally. Public release, anonymous editor/reviewer access, formal archiving, and upstream-rights arrangements remain incomplete. Do not promise “available on reasonable request” unless the named corresponding authors can actually provide the specified assets on those terms. Do not change repository visibility or expose local assets as part of this editorial review.

## AI disclosure

The historical manuscript contains no AI-use statement. Its silence cannot justify omitting assistance during the subsequent reconstruction. The conversation and project records support disclosure of OpenAI Codex assistance with research-plan refinement, analysis-code development/checks, reference checks, interpretation/manuscript drafting and revision, and document formatting. Calling this only language polishing would understate the recorded use.

A draft may state these supported purposes and that human authors retain responsibility for the scientific content. It must not claim that all affected material has already received final author verification while the user is about to begin that review. Tool-provider/product names can be recorded; an unverified uniform model/version for the whole project must not be invented.

## Current official policy check

Checked 20 October 2026 without transmitting source manuscript content:

- [BMC Genomics Research Article instructions](https://link.springer.com/journal/12864/submission-guidelines/research-article), title page and Declarations: supply corresponding contact, author contributions, funding, interests, data access, and applicable ethics/consent statements. A missing historical statement is not automatically “Not applicable.”
- [BMC editorial policies](https://link.springer.com/brands/bmc/editorial-policies), Authorship, AI, and Software/code: human accountability and transparent AI disclosure apply; submitting implies all-author agreement; actual access conditions and an archived code identifier must be represented accurately.

These checks support preparation of an honest author-reading copy. They do not certify policy compliance or submission readiness.

## Recommended handoff

Deliver the new author-reading manuscript with the verified historical author/title-page, funding, and interest facts incorporated. Supply a short separate note for the genuinely absent contact/contribution facts and release-state limitations. Leave scientific results and all prior files untouched. Final submission approval, any new publication-license agreement, and any archive creation remain separate actions after the requested read-through.
