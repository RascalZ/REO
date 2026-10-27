"""Package the rendered English v5 manuscript with an exact 26-file allowlist.

Run from the repository root after both new Word/PDF copies have been rendered
and reviewed. Original v4 files are verified, never edited. This script does
not train, score, read response data, clean folders, or create history archives.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


ROOT = Path(__file__).resolve().parents[1]
V4 = ROOT / "manuscript_v1/author_review_v4"
V5 = ROOT / "manuscript_v1/author_review_v5"
DOCS = ROOT / "docs/public_access_revision_v1"
STAMP = "2026-10-27"
RELEASE_URL = "https://github.com/RascalZ/REO/releases/tag/v0.1.0"
DOI_URL = "https://doi.org/10.5281/zenodo.23681422"
FIGURES = (
    "Figure1_study_design.png",
    "Figure2_development_primary.png",
    "Figure3_representation_and_zero_solutions.png",
    "Figure4_loto_and_external.png",
)
NOTE_NAME = "00_READ_BEFORE_SUBMISSION.txt"
SUMS_NAME = "SHA256SUMS.txt"
PACKAGE_NAME = f"BMC_Genomics_v5_submission_bundle_{STAMP}.zip"


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def checked_file(path: Path) -> bytes:
    """Read only an explicitly selected regular file inside this repository."""
    if not path.resolve().is_relative_to(ROOT.resolve()):
        raise ValueError("Input escapes the repository workspace")
    if path.is_symlink() or not path.is_file():
        raise ValueError("Expected an existing regular publication file: " + path.name)
    return path.read_bytes()


def json_bytes(value: dict) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode("utf-8")


def original_v4_hashes() -> dict[str, str]:
    """Check the six delivered v4 manuscript files and the original v4 ZIP."""
    delivery = json.loads(checked_file(ROOT / "docs/minor_clarification_v1/delivery_manifest.json"))
    originals = delivery["files"]
    if len(originals) != 6:
        raise ValueError("Expected exactly six preserved bilingual v4 manuscript files")
    for relative, expected in originals.items():
        if not relative.startswith(("manuscript_v1/author_review_v4/", "manuscript_v1/author_review_v4_zh/")):
            raise ValueError("Unexpected original-manuscript preservation path")
        if sha_bytes(checked_file(ROOT / relative)) != expected:
            raise ValueError("Original v4 manuscript hash changed: " + Path(relative).name)
    old_package = json.loads(checked_file(ROOT / "docs/repository_release_v1/submission_package_manifest.json"))
    old_path = old_package["package"]
    if old_path != "manuscript_v1/author_review_v4/BMC_Genomics_v4_submission_bundle_2026-10-26.zip":
        raise ValueError("Unexpected original v4 package path")
    if sha_bytes(checked_file(ROOT / old_path)) != old_package["sha256"]:
        raise ValueError("Original v4 submission ZIP changed")
    return {**originals, old_path: old_package["sha256"]}


def input_members() -> dict[str, Path]:
    """Return the 24 publication inputs; generated note/checksums make 26."""
    stem = V5 / "Manuscript_BMC_Genomics_author_review_v5"
    members = {
        "Manuscript.docx": stem.with_suffix(".docx"),
        "Manuscript.pdf": stem.with_suffix(".pdf"),
        "Supporting_sources/Manuscript.md": stem.with_suffix(".md"),
        "Supporting_sources/Supplementary_methods.md": V5 / "Supporting_sources/Supplementary_methods.md",
        "Supporting_sources/Supplementary_methods.pdf": V5 / "Supporting_sources/Supplementary_methods.pdf",
        "Supporting_sources/references.bib": V5 / "Supporting_sources/references.bib",
    }
    for number in range(1, 15):
        name = f"Additional_file_{number}." + ("docx" if number == 14 else "csv")
        members["Additional_files/" + name] = V5 / "Additional_files" / name
    for name in FIGURES:
        members["Figures/" + name] = V5 / "Figures" / name
    if len(members) != 24:
        raise ValueError("Input-member allowlist is not exactly 24 files")
    return members


def note_bytes() -> bytes:
    return (
        "BMC Genomics v5 submission materials\n"
        "27 October 2026\n\n"
        "This package updates the v4 manuscript's access, software-archive citation and\n"
        "license statements following the actual public release. Scientific methods,\n"
        "results, tables and figure data are unchanged. Original v4 files are preserved.\n"
        "This is an author-review package, not evidence of final approval by every author.\n\n"
        "Manuscript.docx: editable English main text with two tables and four embedded figures.\n"
        "Manuscript.pdf: the corresponding reading preview. Select the main-text format\n"
        "required by the submission system; do not upload both as duplicate main manuscripts.\n"
        "Figures/: four separate original 600 dpi PNG figures.\n"
        "Additional_files/: 13 CSV supplementary tables and Additional_file_14.docx methods.\n"
        "Supporting_sources/: English manuscript Markdown, unchanged supplementary methods\n"
        "Markdown/PDF and the updated bibliography, for author retention and use if requested.\n"
        "SHA256SUMS.txt: SHA-256 hashes for the other 25 package members.\n"
        "Unzip this package and assign files to the submission system's requested categories.\n"
        "Chinese reading copies, QA previews and temporary files are not included.\n\n"
        "Published research software and asset release:\n" + RELEASE_URL + "\n"
        "Published software/research-asset archive DOI (not a manuscript DOI):\n" + DOI_URL + "\n\n"
        "Newly authored code is MIT-licensed. Author-owned model and derived assets carry\n"
        "a CC BY 4.0 grant only to the extent that the authors hold the applicable rights.\n"
        "Third-party datasets, annotations and identifiers are not relicensed. Original\n"
        "upstream inputs must be obtained under their own terms. Minimal released metadata\n"
        "and fixed-model inference do not supply every historical curation/rescoring input.\n\n"
        "Before journal submission, the authors must confirm the complete final manuscript,\n"
        "any remaining missing declaration/contact/contribution facts, and final approval\n"
        "by all authors. This package does not invent these facts or submit the manuscript.\n"
        "Existing draft/author-reading status is retained for this review stage.\n\n"
        "Journal instructions:\n"
        "https://link.springer.com/journal/12864/submission-guidelines/research-article\n"
    ).encode("utf-8")


def main() -> None:
    note_path, sums_path = V5 / NOTE_NAME, V5 / SUMS_NAME
    package_path = V5 / PACKAGE_NAME
    manifest_path = DOCS / "submission_package_manifest.json"
    outputs = (note_path, sums_path, package_path, manifest_path)
    if any(path.exists() for path in outputs):
        raise FileExistsError("Refusing to overwrite an existing v5 package output")
    if not V5.is_dir() or V5.is_symlink():
        raise ValueError("The rendered v5 English folder must already exist")

    originals = original_v4_hashes()
    status = json.loads(checked_file(ROOT / "docs/public_release_v1/PUBLICATION_STATUS.json"))
    if (status.get("github_release") != RELEASE_URL or status.get("github_release_state") != "published"
            or status.get("zenodo_doi") != "10.5281/zenodo.23681422" or status.get("zenodo_state") != "published"):
        raise ValueError("Actual release publication evidence does not match the package note")

    sources = input_members()
    payloads = {name: checked_file(path) for name, path in sources.items()}
    main_md = payloads["Supporting_sources/Manuscript.md"].decode("utf-8")
    if re.search(r"[\u3400-\u9fff]", main_md):
        raise ValueError("The English manuscript Markdown contains Chinese text; inspect before packaging")
    if RELEASE_URL not in main_md or "10.5281/zenodo.23681422" not in main_md:
        raise ValueError("English manuscript does not contain both formal release and DOI")
    if "10.5281/zenodo.23681422" not in payloads["Supporting_sources/references.bib"].decode("utf-8"):
        raise ValueError("Updated bibliography does not contain the formal software DOI")
    if not payloads["Manuscript.pdf"].startswith(b"%PDF-"):
        raise ValueError("Rendered main PDF signature is invalid")

    preserved_supports = []
    for name, data in payloads.items():
        if name.startswith(("Figures/", "Additional_files/")) or name in {
            "Supporting_sources/Supplementary_methods.md", "Supporting_sources/Supplementary_methods.pdf"
        }:
            before = checked_file(V4 / name)
            if before != data:
                raise ValueError("Scientific supplementary/figure bytes changed: " + name)
            preserved_supports.append({"name": name, "sha256": sha_bytes(data)})
    if len(preserved_supports) != 20:
        raise ValueError("Expected 20 unchanged scientific supporting files")

    payloads[NOTE_NAME] = note_bytes()
    sums = "".join(f"{sha_bytes(data)}  {name}\n" for name, data in sorted(payloads.items()))
    payloads[SUMS_NAME] = sums.encode("ascii")
    if len(payloads) != 26:
        raise ValueError("Final package must have exactly 26 allowlisted members")
    for name in payloads:
        path = PurePosixPath(name)
        if path.is_absolute() or ".." in path.parts or "\\" in name:
            raise ValueError("Unsafe package member path")

    # Every input and original-preservation check passed before the first write.
    DOCS.mkdir(parents=True, exist_ok=True)
    with note_path.open("xb") as handle:
        handle.write(payloads[NOTE_NAME])
    with sums_path.open("xb") as handle:
        handle.write(payloads[SUMS_NAME])
    with ZipFile(package_path, "x", ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(payloads.items()):
            info = ZipInfo(name, date_time=(2026, 10, 27, 0, 0, 0))
            info.compress_type = ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data, compress_type=ZIP_DEFLATED, compresslevel=9)
    with ZipFile(package_path) as archive:
        if archive.testzip() is not None or sorted(archive.namelist()) != sorted(payloads):
            raise ValueError("ZIP CRC or exact member allowlist verification failed")
        if len(archive.namelist()) != len(set(archive.namelist())):
            raise ValueError("Duplicate ZIP member")
        for name, data in payloads.items():
            if sha_bytes(archive.read(name)) != sha_bytes(data):
                raise ValueError("ZIP member SHA-256 mismatch: " + name)
        for line in archive.read(SUMS_NAME).decode("ascii").splitlines():
            expected, name = line.split("  ", 1)
            if sha_bytes(archive.read(name)) != expected:
                raise ValueError("SHA256SUMS mismatch: " + name)
    if original_v4_hashes() != originals:
        raise ValueError("Original v4 preservation verification changed during packaging")

    sources.update({NOTE_NAME: note_path, SUMS_NAME: sums_path})
    manifest = {
        "date": STAMP,
        "package": package_path.relative_to(ROOT).as_posix(),
        "sha256": sha_bytes(package_path.read_bytes()),
        "bytes": package_path.stat().st_size,
        "member_count": len(payloads),
        "members": [{"name": name, "source": sources[name].relative_to(ROOT).as_posix(),
                     "bytes": len(data), "sha256": sha_bytes(data)} for name, data in sorted(payloads.items())],
        "preserved_v4_originals": originals,
        "unchanged_scientific_support_file_count": len(preserved_supports),
        "unchanged_scientific_support_files": sorted(preserved_supports, key=lambda row: row["name"]),
        "archive_CRC_and_all_member_hashes_verified": True,
        "SHA256SUMS_entries_verified": 25,
        "strict_member_allowlist_verified": True,
        "Chinese_reading_copies_or_QA_files_included": False,
        "scientific_results_changed": False,
        "main_manuscript_scientific_equivalence_scope": "Separate revision audit required; this packager does not edit manuscript content and verifies supporting-file byte preservation only",
        "rendering_review_performed_by_packager": False,
        "training_scoring_or_response_access_performed": False,
        "submission_approval_claimed": False,
        "submission_ready_claimed": False,
        "publication_scope": "English author-review submission materials; only access, licensing and software-archive citation text changed",
    }
    with manifest_path.open("xb") as handle:
        handle.write(json_bytes(manifest))
    print(json.dumps({"package": manifest["package"], "member_count": 26,
                      "bytes": manifest["bytes"], "sha256": manifest["sha256"],
                      "preserved_scientific_support_files": 20}, ensure_ascii=False))


if __name__ == "__main__":
    main()
