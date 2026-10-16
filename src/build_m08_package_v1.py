"""Build a style-only R1 derivative from existing text and aggregate figures.

No science pipeline, models, sample-level data or scoring module is imported.
The preceding manuscripts, builders, tables and figures are never overwritten.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path
from format_m08_reference_authors import apply_author_replacements

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs/m08_v1"
MS = ROOT / "manuscript_v1"
BASE = MS / "Manuscript_BMC_Genomics_scientific_revision_v1.md"
DEST = MS / "Manuscript_BMC_Genomics_M08_v1.md"
PACKAGE = MS / "bmc_additional_files_M08_v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path: Path) -> dict:
    return {"path": path.relative_to(ROOT).as_posix(), "sha256": sha(path)}


def preserved_check() -> None:
    frozen = json.loads((DOCS / "preserved_sources.json").read_text(encoding="utf-8"))
    changed = [name for name, value in frozen["files"].items() if sha(ROOT / name) != value]
    if changed:
        raise RuntimeError(f"Preserved source changed: {changed}")


def main() -> None:
    preserved_check()
    source = BASE.read_text(encoding="utf-8")
    replacements = json.loads((DOCS / "reference_author_audit.json").read_text(encoding="utf-8"))
    md = apply_author_replacements(source, replacements)
    old_status = re.search(r"^\*\*Draft status:\*\* .+$", md, re.M).group(0)
    new_status = ("**Draft status:** Scientific revision R1, M08 visual/style copy, 16 October 2026. "
                  "Scientific content is unchanged from R1. The preceding files remain preserved. "
                  "Author declarations, formal archiving, and upstream-license confirmation remain deferred; "
                  "author-input placeholders remain and this is not a submission-ready manuscript.")
    md = md.replace(old_status, new_status)
    before_body = re.sub(r"^\*\*Draft status:\*\* .+$", "", source.split("## References")[0], flags=re.M)
    after_body = re.sub(r"^\*\*Draft status:\*\* .+$", "", md.split("## References")[0], flags=re.M)
    assert before_body == after_body, "M08 must not change scientific prose or table cells"
    # Captions and Additional-file descriptions also remain verbatim.
    assert source.split("## Additional files")[1] == md.split("## Additional files")[1]
    DEST.write_text(md, encoding="utf-8", newline="\n")
    captions = []
    for n, caption in re.findall(r"^Figure ([1-4])\. (.+)$", md, re.M):
        title, legend = caption.split(". ", 1)
        captions.append({"number": int(n), "title": title, "title_words": len(title.split()),
                         "legend": legend, "legend_words": len(legend.split()),
                         "exact_caption": f"Figure {n}. {caption}"})
    assert len(captions) == 4 and all(x["title_words"] <= 15 and x["legend_words"] <= 300 for x in captions)
    (DOCS / "figure_captions.json").write_text(json.dumps(captions, indent=2) + "\n", encoding="utf-8")
    node = shutil.which("node") or r"D:\Env\NodeJS\node.exe"
    subprocess.run([node, str(MS / "build_bmc_docx_m08_v1.js"), "--embed-figures"], cwd=ROOT, check=True)
    PACKAGE.mkdir(exist_ok=True)
    original_package = MS / "bmc_additional_files_scientific_revision_v1"
    entries = []
    for n in range(1, 15):
        suffix = ".docx" if n == 14 else ".csv"
        source_file = original_package / f"Additional_file_{n}{suffix}"
        target = PACKAGE / source_file.name
        shutil.copyfile(source_file, target)
        assert sha(target) == sha(source_file)
        entries.append({"number": n, **record(target), "source": source_file.relative_to(ROOT).as_posix(),
                        "byte_identical_to_scientific_R1": True})
    figure_manifest = ROOT / "figures/manuscript_m08_v1/figure_manifest.json"
    manifest = {
        "revision": "scientific_R1_M08_style_v1", "date": "2026-10-16",
        "status": "content_review_author_access_rights_archive_pending_not_submission_ready",
        "operation": "editorial_style_and_existing_aggregate_graphics_only",
        "retrained_or_rescored": False, "individual_external_labels_or_predictions_read": False,
        "generator": record(Path(__file__)), "builder": record(MS / "build_bmc_docx_m08_v1.js"),
        "source_manuscript": record(BASE), "figure_manifest": record(figure_manifest),
        "author_replacements": record(DOCS / "reference_author_audit.json"),
        "documents": [record(p) for p in [DEST, DEST.with_suffix(".docx"),
                       MS / "Supplementary_methods_scientific_revision_v1.md",
                       MS / "Supplementary_methods_scientific_revision_v1.docx"]],
        "files": entries,
        "note": "Supplementary content and all Additional files are unchanged byte copies. "
                "Zero, negative and NA results retain their original scientific meaning."
    }
    (PACKAGE / "package_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    preserved_check()
    print(json.dumps({"additional_files": 14, "captions": 4, "scientific_text_unchanged": True,
                      "old_files_unchanged": True}))


if __name__ == "__main__":
    main()
