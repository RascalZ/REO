"""Assemble existing v4 materials and a recoverable old-reading-copy archive.

Only publication documents, figures and aggregate supplementary files are read.
This command never edits manuscripts, trains, scores or deletes anything.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[1]
V4 = ROOT / "manuscript_v1/author_review_v4"
DOCS = ROOT / "docs/repository_release_v1"
STAMP = "2026-10-26"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")


def checked_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        assert sha(source) == sha(destination), f"Preserve differing existing file: {destination.name}"
    else:
        shutil.copy2(source, destination)
    assert sha(source) == sha(destination)


def main() -> None:
    DOCS.mkdir(parents=True, exist_ok=True)
    assert V4.is_dir()
    stem = V4 / "Manuscript_BMC_Genomics_author_review_v4"
    prior = json.loads((ROOT / "docs/minor_clarification_v1/delivery_manifest.json").read_text(encoding="utf-8"))
    assert all(sha(ROOT / name) == value for name, value in prior["files"].items())
    origins = {}
    for i in range(1, 15):
        name = f"Additional_file_{i}." + ("docx" if i == 14 else "csv")
        target = V4 / "Additional_files" / name
        source = target if target.exists() else ROOT / "manuscript_v1/author_review_v2/Additional_files" / name
        checked_copy(source, target)
        origins[target.relative_to(ROOT).as_posix()] = source.relative_to(ROOT).as_posix()
    figure_names = ["Figure1_study_design.png", "Figure2_development_primary.png",
                    "Figure3_representation_and_zero_solutions.png", "Figure4_loto_and_external.png"]
    for name in figure_names:
        source, target = ROOT / "figures/manuscript_m08_v1" / name, V4 / "Figures" / name
        checked_copy(source, target)
        origins[target.relative_to(ROOT).as_posix()] = source.relative_to(ROOT).as_posix()
    for name, original in [
        ("Supplementary_methods.md", "Supplementary_methods_author_review_v2.md"),
        ("Supplementary_methods.pdf", "Supplementary_methods_author_review_v2.pdf"),
        ("references.bib", "references_author_review_v2.bib"),
    ]:
        target = V4 / "Supporting_sources" / name
        source = target if target.exists() else ROOT / "manuscript_v1/author_review_v2" / original
        checked_copy(source, target)
        origins[target.relative_to(ROOT).as_posix()] = source.relative_to(ROOT).as_posix()
    note = V4 / "00_READ_BEFORE_SUBMISSION.txt"
    assert not note.exists()
    note.write_text(
        "BMC Genomics v4 投稿材料包 / Submission materials\n"
        "2026-10-26\n\n"
        "本ZIP保留已确认的v4正文原件，未改科学内容或作者声明。\n"
        "它是上传前整理包，不代表所有作者已完成最终提交批准。\n\n"
        "Manuscript.docx：可编辑英文正文，含两张表和四幅嵌入图。\n"
        "Manuscript.pdf：同版阅读预览，不必与Word重复上传为两份主稿。\n"
        "Figures/：四张独立600 dpi PNG，不含proof或QA图片。\n"
        "Additional_files/：13个CSV表和Additional_file_14.docx补充方法。\n"
        "Supporting_sources/：主稿MD、补充方法MD/PDF及Bib，仅供作者保管，按系统需要选取。\n"
        "SHA256SUMS.txt：包内文件校验和。请先解压，按投稿系统文件类型分别上传。\n\n"
        "仍待真实确认：通讯联系方式（如适用）、逐位作者贡献与致谢、全作者最终批准；\n"
        "编辑/审稿人实际代码与模型访问、再分发许可及正式归档安排。\n"
        "正文仍保留author-reading/draft标识，正式提交前应与最终声明一起收尾。\n"
        "本包未虚构这些事实，也未自动投递稿件或变更仓库可见性。\n\n"
        "指南：https://link.springer.com/journal/12864/submission-guidelines/research-article\n",
        encoding="utf-8", newline="\n")
    members = {"Manuscript.docx": stem.with_suffix(".docx"),
               "Manuscript.pdf": stem.with_suffix(".pdf"),
               "Supporting_sources/Manuscript.md": stem.with_suffix(".md"),
               "00_READ_BEFORE_SUBMISSION.txt": note}
    for folder in ["Figures", "Additional_files", "Supporting_sources"]:
        members.update({f"{folder}/{p.name}": p for p in sorted((V4 / folder).iterdir()) if p.is_file()})
    sums = "".join(f"{sha(path)}  {name}\n" for name, path in sorted(members.items()))
    sums_path = V4 / "SHA256SUMS.txt"
    assert not sums_path.exists()
    sums_path.write_text(sums, encoding="utf-8", newline="\n")
    members["SHA256SUMS.txt"] = sums_path
    package = V4 / f"BMC_Genomics_v4_submission_bundle_{STAMP}.zip"
    assert not package.exists()
    with ZipFile(package, "w", ZIP_DEFLATED, compresslevel=9) as z:
        for name, path in sorted(members.items()):
            z.write(path, name)
    with ZipFile(package) as z:
        assert z.testzip() is None
        assert sorted(z.namelist()) == sorted(members)
        assert all(hashlib.sha256(z.read(name)).hexdigest() == sha(path) for name, path in members.items())
    dump(DOCS / "submission_package_manifest.json", {
        "date": STAMP, "package": package.relative_to(ROOT).as_posix(), "sha256": sha(package),
        "bytes": package.stat().st_size, "member_count": len(members),
        "members": [{"name": name, "source": path.relative_to(ROOT).as_posix(), "sha256": sha(path)}
                    for name, path in sorted(members.items())],
        "migration_origins": origins, "archive_CRC_and_all_member_hashes_verified": True,
        "article_and_declarations_modified": False, "submission_approval_claimed": False})

    archive_dir = ROOT / "local_archives"
    archive_dir.mkdir(exist_ok=True)
    archive = archive_dir / f"author_review_v1_v2_history_{STAMP}.zip"
    assert not archive.exists()
    old_folders = [ROOT / "manuscript_v1" / name for name in
                   ["author_review_v1", "author_review_v2", "author_review_v2_zh"]]
    archived, disposable, targets = [], [], []
    for directory in old_folders:
        assert directory.resolve().parent == (ROOT / "manuscript_v1").resolve()
        files = sorted(p for p in directory.rglob("*") if p.is_file())
        assert not any(p.name.startswith("~$") for p in files), "An old reading copy is open in Word"
        assert not any(p.is_symlink() for p in directory.rglob("*")), "Do not traverse links"
        targets.append({"directory": directory.relative_to(ROOT).as_posix(),
                        "files": len(files), "bytes": sum(p.stat().st_size for p in files)})
        for p in files:
            relative = p.relative_to(ROOT).as_posix()
            row = {"path": relative, "sha256": sha(p), "bytes": p.stat().st_size}
            (disposable if "qa" in p.relative_to(directory).parts else archived).append(row)
    backup_manifest = {"date": STAMP, "preserved_non_QA": archived, "discardable_QA": disposable,
                       "scope": "Only three superseded local author reading folders; no scientific data or results"}
    with ZipFile(archive, "w", ZIP_DEFLATED, compresslevel=9) as z:
        for row in archived:
            z.write(ROOT / row["path"], row["path"])
        z.writestr("ARCHIVE_MANIFEST.json", json.dumps(backup_manifest, ensure_ascii=False, indent=2))
    with ZipFile(archive) as z:
        assert z.testzip() is None
        assert all(hashlib.sha256(z.read(row["path"])).hexdigest() == row["sha256"] for row in archived)
    dump(DOCS / "cleanup_plan.json", {
        "date": STAMP, "targets": targets, "archive": archive.relative_to(ROOT).as_posix(),
        "archive_sha256": sha(archive), "archive_bytes": archive.stat().st_size,
        "archived_file_count": len(archived), "discardable_QA_files": len(disposable),
        "archive_CRC_and_hashes_verified": True, "v4_supplements_migrated_and_verified": True,
        "deleted": False, "scientific_data_models_configs_code_history_in_scope": False,
        "all_files_before_cleanup": archived + disposable})
    assert all(sha(ROOT / name) == value for name, value in prior["files"].items())
    print(json.dumps({"submission_zip": str(package), "members": len(members),
                      "package_bytes": package.stat().st_size, "history_archive": str(archive),
                      "cleanup_targets": targets, "files_deleted": 0}, ensure_ascii=False))


if __name__ == "__main__":
    main()
