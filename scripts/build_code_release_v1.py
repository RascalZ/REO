"""Export a deterministic software archive from one explicit Git commit.

This reads selected tracked blobs with Git plumbing. It never opens source
files from the working tree, traverses ignored files, checks out a commit, or
executes training/prediction/evaluation. Manuscripts and author packages are
excluded before any blob is requested. Outputs require a fresh private folder.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import zipfile


ROOT = Path(__file__).resolve().parents[1]
PRIVATE = ROOT / "docs/public_release_v1/local_private"
VERSION = "0.1.0"
ARCHIVE_NAME = f"reo-code-v{VERSION}.zip"
ALLOWED_DIRECTORIES = {
    ".github", "src", "scripts", "config", "configs", "results", "docs",
    "figures", "data_manifest", "examples", "tests", "reproducibility_v1",
    "tables", "reports",
}
PRIVATE_SEGMENTS = {
    ".git", ".codex", ".agents", ".aws", ".venv", "__pycache__",
    ".pytest_cache", "node_modules", "local_private", "local_assets",
    "local_archives", "checkpoints", "runs", "data_raw", "data_processed",
    "legacy", "_qa_checkout",
}
BINARY_PACKAGE_SUFFIXES = {
    ".zip", ".gz", ".tar", ".tgz", ".7z", ".rar", ".doc", ".docx",
    ".dotx", ".ppt", ".pptx", ".potx", ".xls", ".xlsx", ".xlsm",
    ".npy", ".npz", ".parquet", ".rds", ".rdata", ".pkl", ".pickle",
    ".joblib", ".bin", ".pyc", ".pyo",
}
ROOT_FILE_PATTERN = re.compile(
    r"(?i)^(?:readme(?:[._-].*)?|license(?:[._-].*)?|copying(?:[._-].*)?|"
    r"notice(?:[._-].*)?|citation(?:[._-].*)?|requirements[^/]*\.txt|"
    r"pyproject\.toml|research_plan\.md|reproducing\.md|contributing\.md|setup\.(?:py|cfg)|makefile|dockerfile|"
    r"environment[^/]*\.ya?ml|\.gitignore|\.gitattributes|\.editorconfig)$")
OBSERVATION_NAME_PATTERN = re.compile(
    r"(?i)(?:sample_manifest|global_split_manifest|eligibility_decisions|"
    r"expression_samples|individual|predictions?|responses(?:[_.-]|$)|response_auc|observed_auc|"
    r"(?:^|[_.-])oof(?:[_.-]|$)|bootstrap_group|source_index|"
    r"original_file_inventory|original_source_sha256)")
NOTICE_NAME = "CODE_ARCHIVE_NOTICE.md"
MANIFEST_NAME = "CODE_ARCHIVE_MANIFEST.json"
NOTICE = """# REO v0.1.0 software archive notice

This archive contains selected tracked software and public supporting material
from the exact Git commit recorded in CODE_ARCHIVE_MANIFEST.json. Files were
read from Git blobs, without working-tree substitution. Manuscripts, author
packages, private/ignored assets, raw data, source response indexes and model
payload archives are excluded. Inclusion here does not execute or validate an
experiment or reproduce historical numerical results.

The original repository LICENSE is preserved. Its MIT grant applies to the
software covered by that license; this notice does not replace or broaden it.
Third-party identifiers, annotations, datasets and other third-party material
retain their existing terms and attribution requirements and are not relicensed
by the software license. Existing attribution notices remain in the selected
documentation. Frozen model/metadata assets are separate release attachments
with their own limited author-rights notice and upstream attribution. Obtain
upstream inputs independently under their own terms.

The manifest records included paths, Git blob IDs, file SHA-256 values and every
tracked exclusion with its reason. Untracked/ignored working-tree files are
never enumerated or read. ZIP member order and timestamps are fixed.
""".encode("utf-8")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def encode(value) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def exclusion_reason(path: str, mode: str = "100644", kind: str = "blob"):
    """Decide from tracked names/modes before requesting any file content."""
    parsed = PurePosixPath(path)
    parts = parsed.parts
    lowered = [part.casefold() for part in parts]
    if not parts or parsed.is_absolute() or ".." in parts or "\\" in path:
        return "unsafe_or_nonrelative_path"
    if any(part == "manuscript_v1" or part.startswith("manuscript") for part in lowered):
        return "manuscript_not_part_of_software_archive"
    if any(part in PRIVATE_SEGMENTS for part in lowered):
        return "private_raw_or_nonrelease_tree"
    if any(part.startswith(("author_", "author-")) or "author_package" in part or
           "author-package" in part for part in lowered):
        return "author_package_or_author_material"
    if any(part == ".env" or part.startswith(".env.") for part in lowered):
        return "local_environment_or_secret_file"
    if mode not in {"100644", "100755"} or kind != "blob":
        return "symlink_submodule_or_nonregular_entry"
    if parsed.suffix.casefold() in BINARY_PACKAGE_SUFFIXES:
        return "prebuilt_package_estimator_or_nonsource_binary"
    if parsed.suffix.casefold() == ".pdf" and lowered[0] != "figures":
        return "document_binary_outside_public_figures"
    if len(parts) == 1:
        return None if ROOT_FILE_PATTERN.fullmatch(path) else "root_file_outside_software_allowlist"
    if parts[0] not in ALLOWED_DIRECTORIES:
        return "top_level_tree_outside_software_allowlist"
    if lowered[0] in {"results", "tables", "reports", "data_manifest"} and \
            OBSERVATION_NAME_PATTERN.search(parsed.name):
        return "individual_observation_or_source_index_name"
    return None


def git(*args, input_bytes=None) -> bytes:
    result = subprocess.run(["git", *args], cwd=ROOT, input=input_bytes,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode:
        # Do not print arbitrary Git object contents or host-specific paths.
        raise RuntimeError("Read-only Git command failed: " + args[0])
    return result.stdout


def tracked_entries(commit: str):
    entries = []
    for record in git("ls-tree", "-r", "-z", "--full-tree", commit).split(b"\0"):
        if not record:
            continue
        header, path = record.split(b"\t", 1)
        mode, kind, oid = header.decode("ascii").split(" ")
        entries.append({"path": path.decode("utf-8"), "git_mode": mode,
                        "git_type": kind, "blob_sha": oid})
    return sorted(entries, key=lambda entry: entry["path"])


def blob_contents(entries):
    """Read only already-selected object IDs using one Git batch process."""
    object_ids = list(dict.fromkeys(entry["blob_sha"] for entry in entries))
    batch = git("cat-file", "--batch", input_bytes="".join(oid + "\n" for oid in object_ids).encode("ascii"))
    contents = {}
    offset = 0
    for requested in object_ids:
        end = batch.index(b"\n", offset)
        oid, kind, length = batch[offset:end].decode("ascii").split(" ")
        if oid != requested or kind != "blob":
            raise ValueError("Git batch returned an unexpected object")
        start = end + 1
        stop = start + int(length)
        data = batch[start:stop]
        if len(data) != int(length) or batch[stop:stop + 1] != b"\n":
            raise ValueError("Truncated Git blob stream")
        contents[requested] = data
        offset = stop + 1
    if offset != len(batch):
        raise ValueError("Unexpected trailing Git batch bytes")
    return contents


def zip_bytes(files, modes):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as handle:
        for name, data in sorted(files.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = (0o100755 if modes.get(name) == "100755" else 0o100644) << 16
            handle.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    return buffer.getvalue()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--commit", required=True, help="Explicit complete Git commit object ID; no branch/HEAD/tag")
    parser.add_argument("--out", type=Path, required=True,
                        help="New directory inside docs/public_release_v1/local_private")
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[0-9a-fA-F]{40}|[0-9a-fA-F]{64}", args.commit):
        raise ValueError("--commit must be a complete literal Git commit object ID")
    commit = args.commit.lower()
    out = (ROOT / args.out).resolve() if not args.out.is_absolute() else args.out.resolve()
    if not out.is_relative_to(PRIVATE.resolve()) or out == PRIVATE.resolve():
        raise ValueError("Output must be a new subdirectory of docs/public_release_v1/local_private")
    if out.exists():
        raise FileExistsError("Output already exists; choose a new --out directory")
    if git("cat-file", "-t", commit).strip() != b"commit":
        raise ValueError("The explicit object is not a Git commit")
    tree_sha = git("rev-parse", commit + "^{tree}").decode("ascii").strip()
    tracked = tracked_entries(commit)
    included, excluded = [], []
    for entry in tracked:
        reason = exclusion_reason(entry["path"], entry["git_mode"], entry["git_type"])
        if reason:
            excluded.append({**entry, "reason": reason})
        else:
            included.append(entry)
    if not included or not any(entry["path"] == "LICENSE" for entry in included):
        raise ValueError("The selected commit must include the original root LICENSE")
    if any(entry["path"] in {NOTICE_NAME, MANIFEST_NAME} for entry in included):
        raise ValueError("Reserved archive notice/manifest path already exists in the commit")
    contents = blob_contents(included)
    files = {entry["path"]: contents[entry["blob_sha"]] for entry in included}
    if b"MIT License" not in files["LICENSE"] or b"Permission is hereby granted" not in files["LICENSE"]:
        raise ValueError("Original LICENSE is not the expected MIT software license")
    modes = {entry["path"]: entry["git_mode"] for entry in included}
    records = [{**entry, "bytes": len(files[entry["path"]]), "sha256": sha(files[entry["path"]])}
               for entry in included]
    manifest = {
        "schema_version": "reo-code-archive-1.0", "release_version": VERSION,
        "git_commit": commit, "git_tree_sha": tree_sha,
        "source_contract": "Only selected tracked Git blobs from the explicit commit; no working-tree source reads",
        "software_license": "Original MIT LICENSE preserved; no relicensing of third-party material",
        "manuscript_v1_included": False, "ignored_or_untracked_files_enumerated": False,
        "working_tree_source_files_read": False, "training_scoring_or_prediction_executed": False,
        "selection_scope": "Tracked software, configuration, public aggregates, docs, figures, source manifests, fixtures and tests",
        "included": records, "excluded": excluded,
        "generated_files": [{"path": NOTICE_NAME, "bytes": len(NOTICE), "sha256": sha(NOTICE)}],
        "manifest_self_hash_scope": "This in-archive manifest excludes its own hash; the outer manifest hashes it",
    }
    files[NOTICE_NAME] = NOTICE
    manifest_data = encode(manifest)
    files[MANIFEST_NAME] = manifest_data
    archive_data = zip_bytes(files, modes)
    outer = {**manifest, "archive": {"path": ARCHIVE_NAME, "bytes": len(archive_data),
                                     "sha256": sha(archive_data), "members": len(files)},
             "archive_manifest": {"path": MANIFEST_NAME, "bytes": len(manifest_data),
                                  "sha256": sha(manifest_data)}}
    outputs = {ARCHIVE_NAME: archive_data, "code_manifest.json": encode(outer)}
    outputs["SHA256SUMS"] = "".join(f"{sha(data)}  {name}\n" for name, data in sorted(outputs.items())).encode("ascii")
    out.mkdir(parents=True, exist_ok=False)
    for name, data in sorted(outputs.items()):
        with (out / name).open("xb") as handle:
            handle.write(data)
    print(json.dumps({"release_version": VERSION, "git_commit": commit,
                      "output_relative_to_repo": out.relative_to(ROOT).as_posix(),
                      "included_tracked_blobs": len(included), "excluded_tracked_entries": len(excluded),
                      "archive": outer["archive"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
