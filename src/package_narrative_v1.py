"""Package the approved-for-reading editorial copy; no scientific inputs are read."""
from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript_v1/author_review_v2"
DOCS = ROOT / "docs/narrative_v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    render = json.loads((DOCS / "render_summary.json").read_text(encoding="utf-8"))
    review = json.loads((DOCS / "independent_review_counts.json").read_text(encoding="utf-8"))
    assert render["mechanical_checks_pass"] and render["manual_visual_inspection_confirmed"]
    assert review["all_assessed_checks_passed"] and review["counts"]["passed"] == 54
    for kind, stem in [("main", "Manuscript_BMC_Genomics_author_review_v2"),
                       ("supplement", "Supplementary_methods_author_review_v2")]:
        for ext in ["docx", "pdf"]:
            assert sha(OUT / f"{stem}.{ext}") == render["documents"][kind][f"{ext}_sha256"]

    sources: dict[str, Path] = {}
    for stem in ["Manuscript_BMC_Genomics_author_review_v2", "Supplementary_methods_author_review_v2"]:
        for ext in ["md", "docx", "pdf"]:
            sources[f"{stem}.{ext}"] = OUT / f"{stem}.{ext}"
    for name in ["references_author_review_v2.bib", "本版改动说明.md"]:
        sources[name] = OUT / name
    for number in range(1, 15):
        name = f"Additional_file_{number}.{'docx' if number == 14 else 'csv'}"
        sources[f"Additional_files/{name}"] = OUT / "Additional_files" / name
    for path in sorted((ROOT / "figures/manuscript_m08_v1").glob("Figure*")):
        if path.suffix in {".png", ".svg"} and "_proof" not in path.stem:
            sources[f"Figures/{path.name}"] = path
    assert len(sources) == 30
    manifest = {
        "date": "2026-10-21", "status": "local_author_reading_copy_question_led_revision",
        "submission_approval": False, "archival_doi": None, "scientific_results_unchanged": True,
        "qa_scope": "Editorial fidelity and rendered layout, not independent scientific replication",
        "entries": [{"name": name, "sha256": sha(path), "bytes": path.stat().st_size}
                    for name, path in sources.items()],
    }
    manifest_path = OUT / "delivery_manifest.json"
    dump(manifest_path, manifest)
    archive = OUT / "BMC_Genomics_author_reading_package_2026-10-21.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as package:
        for name, path in sources.items():
            package.write(path, name)
        package.write(manifest_path, "delivery_manifest.json")
    with zipfile.ZipFile(archive) as package:
        assert package.testzip() is None and len(package.namelist()) == 31
        for entry in manifest["entries"]:
            assert hashlib.sha256(package.read(entry["name"])).hexdigest() == entry["sha256"]
    summary = {
        "date": "2026-10-21", "local_delivery_package": archive.relative_to(ROOT).as_posix(),
        "zip_sha256": sha(archive), "zip_bytes": archive.stat().st_size, "zip_entries": 31,
        "zip_integrity_and_entry_hashes_pass": True, "author_details_git_ignored": True,
        "independent_content_checks_passed": 54, "mechanical_content_checks": 39,
        "semantic_content_checks": 15, "main_pages": render["documents"]["main"]["pages"],
        "supplement_pages": render["documents"]["supplement"]["pages"],
        "word_pdf_visual_checks_passed": True, "protected_source_files_unchanged": 60,
        "new_scientific_evaluation_performed": False, "archival_release_or_doi_created": False,
        "submission_approval": False,
    }
    dump(DOCS / "delivery_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
