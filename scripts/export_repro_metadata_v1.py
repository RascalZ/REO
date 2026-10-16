"""Export existing label-free identity metadata and byte-identical frozen folds.

No training, scoring, random splitting, expression/response loading or network
access is implemented here. Run once into a new directory; existing output is
never overwritten. The default output is intentionally a local, non-release
asset package. Public redistribution permissions remain a separate author task.
"""
from __future__ import annotations

import argparse
import codecs
import collections
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
LOCAL_ROOT = ROOT / "reproducibility_v1/local_assets"
DEFAULT_OUT = LOCAL_ROOT / "metadata"
DEV_RUN = "runs/gdsc2-cancer-v2-fixed8-mad-association-expression-rank-pair-nested5x3-inner3-repeated-v1"
FINAL_RUN = "runs/gdsc2-cancer-v2-fixed8-full-development-final-v1_1"
LOTO_RUN = "runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_3"
DEV_INPUT = "results/cohort_v2/development_sample_manifest.csv"
EXTERNAL_INPUT = "results/cohort_v2/external_primary_sample_manifest.csv"
FOLD_INPUTS = [
    (f"{DEV_RUN}/global_fold_manifest.csv", f"{DEV_RUN}/ordinal_repeated_metadata.json", "development_repeated_global_folds.csv", 31926),
    (f"{FINAL_RUN}/global_inner_fold_manifest.csv", f"{FINAL_RUN}/final_model_metadata.json", "final_development_inner_folds.csv", 1878),
    (f"{LOTO_RUN}/global_fold_manifest.csv", f"{LOTO_RUN}/loto_metadata.json", "loto_global_outer_folds.csv", 14398),
    (f"{LOTO_RUN}/loto_inner_fold_manifest.csv", f"{LOTO_RUN}/loto_metadata.json", "loto_drug_inner_assignments.csv", 109274),
]
ALLOWED_INPUTS = {
    DEV_INPUT, EXTERNAL_INPUT,
    *(item[0] for item in FOLD_INPUTS),
    *(item[1] for item in FOLD_INPUTS),
}
PRESENCE_ONLY_FIELDS = {"has_response", "response_present", "response_available", "has_ctrp_response"}
DEV_COLUMNS = {
    "gdsc_drug_id", "drug_name", "canonical_model_id", "biological_group_id", "ccle_id",
    "expression_row_index", "tissue", "gdsc_source_row_number", "nlme_curve_id", "min_conc",
    "max_conc", "role", "v1_tissue", "eligibility_status", "evidence_source", "cohort_version",
}
EXTERNAL_COLUMNS = {
    "gdsc_drug_id", "drug_name", "ctrp_compound_id", "canonical_model_id", "biological_group_id",
    "ccle_id", "expression_row_index", "tissue", "ctrp_master_ccl_ids", "ctrp_experiment_ids",
    "ctrp_source_row_numbers", "n_source_experiments", "identity_qc_statuses", "external_subtype",
    "future_response_aggregation", "role", "v1_tissue", "eligibility_status", "evidence_source",
    "cohort_version",
}
IDENTITY_REQUIRED = ["canonical_model_id", "biological_group_id", "ccle_id", "expression_row_index", "tissue"]
IDENTITY_OPTIONAL = ["role", "cohort_version", "eligibility_status", "evidence_source"]
DEV_MEMBERSHIP = [
    "gdsc_drug_id", "drug_name", "canonical_model_id", "biological_group_id", "ccle_id",
    "expression_row_index", "tissue", "gdsc_source_row_number", "nlme_curve_id", "min_conc",
    "max_conc", "role", "cohort_version", "eligibility_status", "evidence_source",
]
EXTERNAL_MEMBERSHIP = [
    "gdsc_drug_id", "drug_name", "ctrp_compound_id", "canonical_model_id", "biological_group_id",
    "ccle_id", "expression_row_index", "tissue", "ctrp_master_ccl_ids", "ctrp_experiment_ids",
    "ctrp_source_row_numbers", "n_source_experiments", "identity_qc_statuses", "external_subtype",
    "role", "cohort_version", "eligibility_status", "evidence_source",
]
EXPECTED_DRUGS = {"1931", "1036", "1191", "1930", "1199", "1091", "1010", "1922"}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class MetadataAuditGuard:
    """Deny every non-whitelisted file read and every file write outside out.

    The stdlib imports and codec lookup finish before installation. The only
    permitted data inputs are six CSV files and three existing metadata JSONs.
    No subprocess or network operation is permitted. New repository-local
    ancestor directories may be created solely to establish the chosen output.
    File writes require exclusive-create semantics inside the fresh out tree.
    """

    def __init__(self, out: Path):
        self.out = out
        self.allowed = {(ROOT / item).resolve() for item in ALLOWED_INPUTS}
        self.read_paths: set[str] = set()
        self.write_paths: set[str] = set()
        self.blocked_attempts = 0

    def reject(self, message: str):
        self.blocked_attempts += 1
        raise PermissionError(message)

    def __call__(self, event, args):
        if event == "open":
            name, mode, flags = args
            if isinstance(name, int):
                self.reject("Opening an existing file descriptor is outside this export contract")
            path = Path(os.fsdecode(name)).resolve()
            is_write = bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND))
            if is_write:
                if not path.is_relative_to(self.out) or not (flags & os.O_EXCL):
                    self.reject("Only exclusive new files inside the export output are writable")
                self.write_paths.add(path.relative_to(self.out).as_posix())
            else:
                if path not in self.allowed:
                    self.reject("Input file is outside the six-CSV/three-JSON metadata whitelist")
                self.read_paths.add(path.relative_to(ROOT).as_posix())
        elif event == "os.mkdir":
            path = Path(os.fsdecode(args[0])).resolve()
            if not (path.is_relative_to(self.out) or
                    (self.out.is_relative_to(path) and path.is_relative_to(ROOT))):
                self.reject("Directory creation is outside the local export tree")
        elif event in {"os.remove", "os.rmdir", "os.rename", "os.link", "os.symlink", "os.chmod", "os.truncate"}:
            self.reject("Deletion, replacement and link operations are not part of this export")
        elif event.startswith(("socket.", "subprocess.")) or event in {"os.system", "os.exec", "os.posix_spawn"}:
            self.reject("Network and process execution are not part of this metadata export")


def header_only(path: Path) -> list[str]:
    # Unbuffered readline stops at the header, so a rejected schema does not
    # require reading/hash-computing the remaining file body.
    with path.open("rb", buffering=0) as handle:
        first_line = handle.readline()
    header = next(csv.reader([first_line.decode("utf-8-sig")]))
    if not header or len(header) != len(set(header)) or any(not value for value in header):
        raise ValueError("Empty or duplicate CSV header fields")
    return header


def csv_bytes(data: bytes) -> tuple[list[str], list[dict[str, str]]]:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("Malformed CSV row width")
    return list(reader.fieldnames or []), rows


def encode_csv(rows: list[dict], fields: list[str]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8-sig")


def encode_json(value) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def write_new(out: Path, relative: str, data: bytes):
    path = out / relative
    with path.open("xb") as handle:
        handle.write(data)


def safe_manifest_csv(relative: str, allowed_fields: set[str], sources: list[dict]):
    path = ROOT / relative
    if not path.exists():
        sources.append({"relative_path": relative, "status": "missing", "sha256": None})
        return None
    header = header_only(path)
    rejected = sorted(set(header) - allowed_fields - PRESENCE_ONLY_FIELDS)
    if rejected:
        sources.append({"relative_path": relative, "status": "header_rejected_body_not_read",
                        "sha256": None, "header_columns": header, "rejected_fields": rejected})
        return None
    data = path.read_bytes()
    columns, rows = csv_bytes(data)
    if columns != header:
        raise ValueError("Header changed between schema check and body read")
    sources.append({"relative_path": relative, "status": "label_free_metadata_read", "sha256": digest(data),
                    "bytes": len(data), "data_rows": len(rows), "columns": columns,
                    "presence_only_fields_not_exported": sorted(set(header) & PRESENCE_ONLY_FIELDS),
                    "response_value_fields_present": False})
    return rows


def project_cohort(rows: list[dict], *, cohort: str, expected_rows: int, expected_ids: int):
    if len(rows) != expected_rows or {row["gdsc_drug_id"] for row in rows} != EXPECTED_DRUGS:
        raise ValueError(f"{cohort}: frozen drug membership coverage changed")
    memberships = DEV_MEMBERSHIP if cohort == "development" else EXTERNAL_MEMBERSHIP
    required = set(IDENTITY_REQUIRED + memberships)
    if not all(required.issubset(row) for row in rows):
        raise ValueError(f"{cohort}: required metadata fields missing")
    ids: dict[str, list[dict]] = collections.defaultdict(list)
    membership_keys = set()
    for row in rows:
        if any(not row[field] for field in IDENTITY_REQUIRED):
            raise ValueError(f"{cohort}: an identity field is empty")
        if cohort == "development" and row["role"] != "development":
            raise ValueError("Non-development record in development manifest")
        key = (row["gdsc_drug_id"], row["canonical_model_id"])
        if key in membership_keys:
            raise ValueError(f"{cohort}: duplicate drug/model metadata record")
        membership_keys.add(key)
        ids[row["canonical_model_id"]].append(row)
    if len(ids) != expected_ids:
        raise ValueError(f"{cohort}: identity inventory changed")
    for group_rows in ids.values():
        for field in IDENTITY_REQUIRED:
            if len({row[field] for row in group_rows}) != 1:
                raise ValueError(f"{cohort}: inconsistent required identity metadata")
    fields = list(IDENTITY_REQUIRED)
    omitted = []
    candidates = IDENTITY_OPTIONAL + ([] if cohort == "development" else ["ctrp_master_ccl_ids", "identity_qc_statuses", "external_subtype"])
    for field in candidates:
        if all(field in row for row in rows) and all(len({row[field] for row in group_rows}) == 1 for group_rows in ids.values()):
            fields.append(field)
        else:
            omitted.append(field)
    identities = [{field: ids[key][0][field] for field in fields} for key in sorted(ids)]
    if len({row["biological_group_id"] for row in identities}) != expected_ids:
        raise ValueError(f"{cohort}: one-model-per-group condition changed")
    projected = [{field: row[field] for field in memberships} for row in rows]
    counts = {
        "identities": len(identities), "biological_groups": len({row["biological_group_id"] for row in identities}),
        "drug_memberships": len(projected), "drugs": len(EXPECTED_DRUGS),
        "tissues": len({row["tissue"] for row in identities}),
        "drug_membership_counts": dict(sorted(collections.Counter(row["gdsc_drug_id"] for row in projected).items())),
        "optional_identity_fields_omitted_because_absent_or_variable": omitted,
    }
    return identities, fields, projected, memberships, counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT,
                        help="New directory under reproducibility_v1/local_assets; never overwrite")
    args = parser.parse_args()
    out = args.out.resolve()
    if not out.is_relative_to(LOCAL_ROOT) or out == LOCAL_ROOT:
        raise ValueError("Output must be a new subdirectory of reproducibility_v1/local_assets")
    if out.exists():
        raise FileExistsError("Output already exists; preserve it and select a new --out directory")
    codecs.lookup("utf-8-sig")
    guard = MetadataAuditGuard(out)
    sys.addaudithook(guard)
    out.mkdir(parents=True, exist_ok=False)
    sources: list[dict] = []
    outputs: list[dict] = []
    counts = {}
    cohort_identity_rows = {}

    def add_output(relative, data, *, source_paths, rows=None, columns=None, byte_identical=False):
        write_new(out, relative, data)
        outputs.append({"relative_path": relative, "sha256": digest(data), "bytes": len(data),
                        "data_rows": rows, "columns": columns,
                        "source_relative_paths": source_paths, "byte_identical_source_copy": byte_identical,
                        "response_value_fields_present": False, "redistribution_status": "pending_author_and_license_review"})

    for cohort, path, allowed, expected_rows, expected_ids in [
        ("development", DEV_INPUT, DEV_COLUMNS, 4967, 626),
        ("external", EXTERNAL_INPUT, EXTERNAL_COLUMNS, 1185, 175),
    ]:
        rows = safe_manifest_csv(path, allowed, sources)
        if rows is None:
            counts[cohort] = {"status": "not_exported_missing_or_header_rejected"}
            continue
        identities, fields, memberships, membership_fields, stats = project_cohort(
            rows, cohort=cohort, expected_rows=expected_rows, expected_ids=expected_ids)
        cohort_identity_rows[cohort] = identities
        counts[cohort] = {"status": "exported", **stats}
        add_output(f"{cohort}_identities.csv", encode_csv(identities, fields), source_paths=[path], rows=len(identities), columns=fields)
        add_output(f"{cohort}_drug_membership.csv", encode_csv(memberships, membership_fields), source_paths=[path], rows=len(memberships), columns=membership_fields)

    if set(cohort_identity_rows) == {"development", "external"}:
        dev_ids = {row["canonical_model_id"] for row in cohort_identity_rows["development"]}
        ext_ids = {row["canonical_model_id"] for row in cohort_identity_rows["external"]}
        dev_groups = {row["biological_group_id"] for row in cohort_identity_rows["development"]}
        ext_groups = {row["biological_group_id"] for row in cohort_identity_rows["external"]}
        counts["development_external_overlap"] = {"canonical_models": len(dev_ids & ext_ids), "biological_groups": len(dev_groups & ext_groups)}
        if dev_ids & ext_ids or dev_groups & ext_groups:
            raise ValueError("Development and external metadata overlap")

    metadata_cache = {}
    for source, metadata_path, target, expected_rows in FOLD_INPUTS:
        if metadata_path not in metadata_cache:
            raw = (ROOT / metadata_path).read_bytes()
            metadata_cache[metadata_path] = json.loads(raw.decode("utf-8-sig"))
            sources.append({"relative_path": metadata_path, "status": "run_metadata_read", "sha256": digest(raw), "bytes": len(raw)})
        raw = (ROOT / source).read_bytes()
        expected_hash = metadata_cache[metadata_path]["prepared_input_hashes"][Path(source).name]
        if digest(raw) != expected_hash:
            raise ValueError("Existing fold file differs from its bound prepared hash")
        fields, fold_rows = csv_bytes(raw)
        if len(fold_rows) != expected_rows:
            raise ValueError("Existing fold row count differs from frozen inventory")
        allowed_fold_fields = {"canonical_model_id", "biological_group_id", "tissue", "repeat", "outer_fold",
                               "inner_fold", "role", "heldout_tissue", "drug_id", "tissue_index", "validation_fold", "inner_seed"}
        if not set(fields).issubset(allowed_fold_fields):
            raise ValueError("Unexpected fold fields: output was not copied")
        sources.append({"relative_path": source, "status": "existing_folds_hash_verified", "sha256": digest(raw),
                        "bytes": len(raw), "data_rows": len(fold_rows), "columns": fields,
                        "metadata_relative_path": metadata_path, "metadata_bound_sha256": expected_hash,
                        "source_hash_binding_match": True, "response_value_fields_present": False})
        add_output(target, raw, source_paths=[source], rows=len(fold_rows), columns=fields, byte_identical=True)

    schema = {
        "schema_version": "1.0", "scope": "Existing source metadata projections and byte-identical fold copies",
        "identity_key": ["canonical_model_id"], "group_key": ["biological_group_id"],
        "drug_membership_key": ["gdsc_drug_id", "canonical_model_id"],
        "identity_fields": {field: "Exact string from the frozen manifest; not a newly inferred mapping" for field in IDENTITY_REQUIRED},
        "expression_row_index_semantics": "Zero-based offset in this project's prepared expression matrix; requires matching sample ordering",
        "source_id_semantics": {
            "ccle_id": "Exact CCLE name from the frozen manifest",
            "ctrp_master_ccl_ids": "Existing source field, retained verbatim where consistent; no new aggregation or identity inference",
            "ctrp_experiment_ids": "Existing source experiment identifier string, not an external label",
            "ctrp_source_row_numbers": "Existing source metadata row identifiers, not an external label",
            "gdsc_source_row_number": "Existing source row identifier",
            "nlme_curve_id": "Existing GDSC curve identifier",
        },
        "presence_only_fields": "May be inspected only as existing source metadata; never exported in projections",
        "fold_interpretation": "docs/reproducibility_assets_v1/FOLD_DEFINITIONS.md",
        "output_schema": [{"relative_path": item["relative_path"], "columns": item["columns"], "data_rows": item["data_rows"]} for item in outputs],
        "redistribution_status": "Local preparation only; source data permissions, author declarations, formal archive and license remain pending",
        "not_a_complete_retraining_package": "No expression matrix, response labels, trained model payload or bootstrap draws are included",
    }
    add_output("mapping_schema.json", encode_json(schema), source_paths=[DEV_INPUT, EXTERNAL_INPUT])
    gates = {
        "no_ctrp_response_labels_read": True, "no_response_tables_read": True,
        "no_individual_predictions_read": True, "no_ctrp_raw_or_index_read": True,
        "no_model_payload_read": True, "no_training": True, "no_scoring": True,
        "no_random_generation_or_fold_regeneration": True, "existing_outputs_overwritten": False,
        "all_four_fold_source_hash_bindings_match": True,
        "blocked_access_attempts": guard.blocked_attempts,
        "data_read_whitelist_files": len(ALLOWED_INPUTS), "actual_input_paths_read": sorted(guard.read_paths),
    }
    manifest = {
        "schema_version": "1.0", "exporter": "scripts/export_repro_metadata_v1.py",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "output_relative_to_repo": out.relative_to(ROOT).as_posix(),
        "sources": sources, "outputs": outputs, "counts": counts, "gates": gates,
        "redistribution_status": "pending_author_and_license_review",
    }
    manifest_bytes = encode_json(manifest)
    write_new(out, "metadata_export_manifest.json", manifest_bytes)
    public_summary = {
        "schema_version": "1.0", "exporter": "scripts/export_repro_metadata_v1.py",
        "manifest_sha256": digest(manifest_bytes), "counts": counts,
        "source_hashes": {item["relative_path"]: item.get("sha256") for item in sources},
        "output_hashes": {item["relative_path"]: item["sha256"] for item in outputs},
        "gates": {key: value for key, value in gates.items() if key != "actual_input_paths_read"},
        "contains_individual_ids": False,
        "redistribution_status": "Local metadata assets prepared; public redistribution permissions remain pending",
    }
    write_new(out, "metadata_export_public_summary.json", encode_json(public_summary))
    print(json.dumps({"status": "completed", "output": out.relative_to(ROOT).as_posix(), "counts": counts,
                      "fold_copies": 4, "manifest_sha256": digest(manifest_bytes),
                      "blocked_access_attempts": guard.blocked_attempts}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
