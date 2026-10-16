# M08 reference-author audit

Audit date: 2026-10-16. This is an editorial author-formatting audit, not a new literature search or scientific evaluation.

The 30 reference identities and numbering remain fixed. The author map changes 19 prefixes; 14 references retain six names followed by `et al.`. All authors are named for papers with at most six authors. Three corporate authors remain verbatim. No author metadata are unresolved.

The selected policy preserves family-name particles, compound surnames, accents and hyphenated given-name initials. Names are based on the saved verified bibliography and normalized Crossref records. References 27–30 also have the saved Crossref/PubMed/official NeurIPS evidence in `manuscript_v1/aivc_reference_metadata.json`. The original BibTeX files are retained.

The [current BMC Genomics research-article instructions](https://link.springer.com/journal/12864/submission-guidelines/research-article) identify Vancouver examples and illustrate six authors before `et al.`. The page does not explicitly prescribe the numeric cutoff in prose; six-name truncation is recorded as this manuscript's consistent implementation, rather than presented as a uniquely mandated rule. Its dataset example also preserves hyphenated given-name initials. The [general instructions](https://link.springer.com/journal/12864/submission-guidelines) specify Index Medicus/MEDLINE journal abbreviations; this author-only map leaves the previously normalized journal abbreviations unchanged.

One targeted metadata correction is documented: Ghandi 2019's sixth author is rendered `McDonald ER 3rd`. Saved Crossref metadata lack the suffix, but [PubMed PMID 31068700](https://pubmed.ncbi.nlm.nih.gov/31068700/) and the [Nature author list](https://www.nature.com/articles/s41586-019-1186-3) confirm it. Corporate authors, the conference source type, and the existing ahead-of-print citation are preserved.

`reference_author_audit.json` contains the source hashes and all 30 exact old/new author prefixes and complete reference lines. `format_m08_reference_authors.py` exposes `build_audit(root)` and `apply_author_replacements(text, audit)` for the M08 manuscript builder; default execution writes this audit only. The transform verifies that every reference tail after its author block is unchanged. It never reads study data, fits a model, or scores predictions.
