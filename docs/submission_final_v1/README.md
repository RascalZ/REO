# Formal English submission preparation

The local manuscript folder is `manuscript_v1/submission_final_v1/` and is Git-ignored. It contains the English main manuscript, a new BMC Genomics cover letter, the unchanged four figures and 14 additional files, reading PDFs, and separate author-only contact/check files. Original English and Chinese v5 files and the old ZIP are preserved.

## Editorial and verification boundary

Seven bounded Markdown edits remove the reading-copy heading/status and replace correspondence, contribution, acknowledgement, affiliation and funding blocks. Eight historical author roles were obtained from the separate contribution record accompanying the last Human Genomics submission. Current author replies confirmed no additional acknowledgements and supplied updated institution/address, grants and corresponding emails. The eight-author manuscript order has not been changed by inference.

Abstract, Background, Methods, Results, Discussion, Conclusions and all 31 references remain identical to v5. Two native tables, four embedded images and 21 supporting-source files retain their contents. No model training, parameter tuning, resampling, scoring or individual CTRP response access was performed.

`render_docx.py` was attempted but LibreOffice was unavailable. The isolated native Microsoft Word fallback exported the main manuscript and cover letter read-only. All 40 main pages and the one-page cover letter were rendered to 2000-pixel PNGs and visually inspected. The exact-file recheck did not confirm an initial batched-image concern about line numbers; no numbering edit was made.

The ZIP has 31 strictly allowlisted members; CRC, all member hashes and 30 SHA256SUMS entries passed. Only 20 are journal-upload files: main DOCX, cover DOCX, four figures and 14 additional files. Previews, sources, author contact/check files and checksums are explicitly separated from publishable supplements.

## Author gate

Current-version approval and exclusive/unpublished submission status were confirmed by the user for the version under review. An author-list/email discrepancy remains under clarification; the authors must supply an unambiguous final list and any missing system contact information before actual submission. No additional author is inferred from an unmatched email. A changed list requires corresponding contribution and approval confirmation.

The manuscript was not submitted. Author names, emails, private source text, the manuscript and the ZIP are not uploaded to this public audit folder. Published source/research-asset release v0.1.0 and Zenodo DOI 10.5281/zenodo.23681422 are unchanged.

## Records

- `build_manifest.json`: source hashes/locators, permitted block IDs and supporting-file hashes.
- `mechanical_audit.json`: section equality, original preservation, document and render checks.
- `visual_review.json`: actual page coverage and observed layout findings.
- `submission_package_manifest.json`: exact ZIP members and hash verification.
- Local-only `local_private/`: source text/configuration, full preservation hashes, author replies and email mapping.

Guidelines checked 27 October 2026: [Research article](https://link.springer.com/journal/12864/submission-guidelines/research-article), [submission preparation](https://link.springer.com/journal/12864/submission-guidelines). The former specifies declaration sections and individual contributions; the latter specifies cover-letter approval/exclusive-submission statements, author emails and upload formats.
