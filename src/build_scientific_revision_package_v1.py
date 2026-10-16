"""Package existing tables and render editorial R1 DOCX copies; never score data."""
from __future__ import annotations
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript_v1"
DEST = MANUSCRIPT / "bmc_additional_files_scientific_revision_v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    DEST.mkdir(exist_ok=True)
    subprocess.run(["node", str(MANUSCRIPT / "build_bmc_docx.js"), "--input", "Manuscript_BMC_Genomics_scientific_revision_v1.md", "--output", "Manuscript_BMC_Genomics_scientific_revision_v1.docx", "--embed-figures"], cwd=ROOT, check=True)
    subprocess.run(["node", str(MANUSCRIPT / "build_bmc_docx.js"), "--input", "Supplementary_methods_scientific_revision_v1.md", "--output", "Supplementary_methods_scientific_revision_v1.docx"], cwd=ROOT, check=True)
    sources = {i: MANUSCRIPT / f"bmc_additional_files/Additional_file_{i}.csv" for i in range(1, 9)}
    tables = ROOT / "tables/scientific_revision_v1"
    for i, name in enumerate(["TableS9_external_absolute_performance.csv", "TableS10_representation_comparisons.csv", "TableS11_synthetic_cohort_results.csv", "TableS12_synthetic_generation_checks.csv", "TableS13_execution_versions.csv"], 9):
        sources[i] = tables / name
    sources[14] = MANUSCRIPT / "Supplementary_methods_scientific_revision_v1.docx"
    entries = []
    for i, source in sources.items():
        target = DEST / f"Additional_file_{i}{source.suffix}"
        shutil.copyfile(source, target)
        assert sha(source) == sha(target)
        row = {"number": i, "source": source.relative_to(ROOT).as_posix(), "path": target.relative_to(ROOT).as_posix(), "sha256": sha(target), "byte_identical_to_source": True}
        if source.suffix == ".csv":
            with source.open(encoding="utf-8-sig", newline="") as f:
                records = list(csv.reader(f))
            row.update(rows=len(records)-1, columns=len(records[0]))
        entries.append(row)
    manifest = {
        "revision": "scientific_R1", "date": "2026-10-16", "author_declarations_archive_rights": "deferred_not_submission_ready",
        "operation": "copy_existing_aggregates_and_build_documents_only_no_response_access_or_scoring",
        "generator_sha256": sha(Path(__file__)), "docx_builder_sha256": sha(MANUSCRIPT / "build_bmc_docx.js"),
        "documents": [{"path": (MANUSCRIPT / n).relative_to(ROOT).as_posix(), "sha256": sha(MANUSCRIPT / n)} for n in ["Manuscript_BMC_Genomics_scientific_revision_v1.md", "Manuscript_BMC_Genomics_scientific_revision_v1.docx", "Supplementary_methods_scientific_revision_v1.md", "Supplementary_methods_scientific_revision_v1.docx"]],
        "files": entries,
        "note": "S1-S8 preserved byte-for-byte; missing precomputed CIs remain NA; original blank synthetic undefined values remain blank; external endpoints restricted to Spearman ranking, not uncalibrated cross-screen error."
    }
    (DEST / "package_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"additional_files": len(entries), "documents": 2, "all_source_copies_identical": True}))


if __name__ == "__main__":
    main()
