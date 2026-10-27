"""Build deterministic, outcome-free v0.1.0 release assets from frozen exports.

Only fixed allowlisted files under local_assets and the release notice templates
are read. No source response/index/prediction file, estimator, network, training
or scoring operation is used. Existing output directories are never replaced.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import zipfile


ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "reproducibility_v1/local_assets"
DOCS = ROOT / "docs/public_release_v1"
PRIVATE = DOCS / "local_private"
VERSION = "0.1.0"
MODEL_ARCHIVE = f"reo-model-assets-v{VERSION}.zip"
METADATA_ARCHIVE = f"reo-metadata-assets-v{VERSION}.zip"
DRUGS = {1010, 1036, 1091, 1191, 1199, 1922, 1930, 1931}
METHODS = {"tissue_null", "expression_strong", "mad_expression", "mad_rank",
           "mad_pair_full", "response_association_expression",
           "response_association_rank", "response_association_pair_full"}
IDENTITY_FIELDS = ["canonical_model_id", "biological_group_id", "ccle_id",
                   "expression_row_index", "tissue", "role", "cohort_version",
                   "eligibility_status"]
MEMBERSHIP_FIELDS = ["gdsc_drug_id", "drug_name", *IDENTITY_FIELDS]
FOLDS = {
    "development_repeated_global_folds.csv": 31926,
    "final_development_inner_folds.csv": 1878,
    "loto_global_outer_folds.csv": 14398,
    "loto_drug_inner_assignments.csv": 109274,
}
FOLD_FIELDS = {"canonical_model_id", "biological_group_id", "tissue", "repeat",
               "outer_fold", "inner_fold", "role", "heldout_tissue", "drug_id",
               "tissue_index", "validation_fold", "inner_seed"}
TABLES = {"tables/candidate_genes.csv": 25600,
          "tables/all_coefficients.csv.gz": 340800,
          "tables/tissue_effects.csv": 1472}
MODEL_FILES = {f"models/d{drug}_{method}.json.gz"
               for drug in DRUGS for method in METHODS}
METADATA_FILES = {f"metadata/{cohort}_{kind}.csv" for cohort in
                  ("development", "external") for kind in
                  ("identities", "drug_membership")} | {f"metadata/{name}" for name in FOLDS}
SCHEMA_FILES = {"schemas/input_gene_ids.json", "schemas/tissue_code.json"}
SYNTHETIC = "synthetic_inputs/full_gene_order.synthetic.json"
ALLOWLIST = MODEL_FILES | METADATA_FILES | set(TABLES) | SCHEMA_FILES | {
    SYNTHETIC, "assets_manifest_v1.json", "model_export_manifest.json"}
FORBIDDEN_FIELDS = {
    "response_auc", "observed_auc", "predicted_auc", "response", "responses",
    "label", "labels", "prediction", "predictions", "y", "y_true", "y_pred",
    "ctrp_compound_id", "ctrp_master_ccl_ids", "ctrp_experiment_ids",
    "ctrp_source_row_numbers", "n_source_experiments", "identity_qc_statuses",
    "external_subtype", "future_response_aggregation", "gdsc_source_row_number",
    "nlme_curve_id", "min_conc", "max_conc", "evidence_source",
}
PRIVATE_PATTERN = re.compile(
    r"(?i)(?:(?<![a-z0-9])[a-z]:[\\/]|/(?:home|users|mnt)/|\b(?:ghp|github_pat)_[a-z0-9_]+|"
    r"\b(?:api[_-]?key|password|secret[_-]?key)\s*[:=]|\bbearer\s+[a-z0-9._-]+)")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def encode(value) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode("utf-8")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows_of(data: bytes):
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    fields = reader.fieldnames or []
    require(fields and len(fields) == len(set(fields)), "Invalid CSV header")
    rows = list(reader)
    require(all(None not in row and all(v is not None for v in row.values())
                for row in rows), "Malformed CSV row width")
    return fields, rows


def csv_encode(fields, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8-sig")


def field_audit(value, *, allow_generated_tpm=False):
    """Check keys/strings only; do not expose fitted vectors or identity rows."""
    if isinstance(value, dict):
        require(not (set(value) & FORBIDDEN_FIELDS), "Forbidden observation/source field")
        require(allow_generated_tpm or "raw_tpm" not in value, "Observed expression field")
        for child in value.values():
            field_audit(child, allow_generated_tpm=allow_generated_tpm)
    elif isinstance(value, list):
        for child in value:
            if isinstance(child, (str, dict, list)):
                field_audit(child, allow_generated_tpm=allow_generated_tpm)
    elif isinstance(value, str):
        require(not PRIVATE_PATTERN.search(value), "Machine path or credential pattern")


def snapshot():
    return {p.relative_to(LOCAL).as_posix(): (p.stat().st_size, p.stat().st_mtime_ns)
            for p in LOCAL.rglob("*") if p.is_file()}


def archive(entries):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as handle:
        for name, data in sorted(entries.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            handle.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    return buffer.getvalue()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=Path("docs/public_release_v1/local_private/artifacts_v0.1.0"))
    args = parser.parse_args(argv)
    out = (ROOT / args.out).resolve() if not args.out.is_absolute() else args.out.resolve()
    require(out.is_relative_to(PRIVATE.resolve()) and out != PRIVATE.resolve(),
            "Output must be a new subdirectory of docs/public_release_v1/local_private")
    if out.exists():
        raise FileExistsError("Output already exists; select a new --out directory")
    before = snapshot()
    opened = {}

    def read(relative):
        require(relative in ALLOWLIST, "Input is outside the fixed release allowlist")
        path = (LOCAL / relative).resolve()
        require(path.is_relative_to(LOCAL.resolve()), "Input escapes local_assets")
        data = path.read_bytes()
        opened[relative] = {"bytes": len(data), "sha256": sha(data)}
        return data

    inventory = json.loads(read("assets_manifest_v1.json"))
    inventory_by_path = {entry["path"]: entry for entry in inventory["files"]}
    export = json.loads(read("model_export_manifest.json"))
    require(export["model_count"] == 64 and len(export["models"]) == 64,
            "Frozen model count differs")
    model_records = {entry["path"]: entry for entry in export["models"]}
    require(set(model_records) == MODEL_FILES, "Frozen drug/method model grid differs")

    def frozen(relative):
        data = read(relative)
        bound = inventory_by_path[relative]
        require(sha(data) == bound["sha256"] and len(data) == bound["bytes"],
                "Prepared asset differs from its frozen inventory: " + relative)
        return data

    model_entries = {}
    metadata_entries = {}
    catalog = []
    preservation = []
    projection = []
    memberships = {}
    cohort_ids = {}
    fold_audits = []
    sensitive = {"forbidden_field_names": sorted(FORBIDDEN_FIELDS),
                 "machine_path_or_credential_pattern_matches": 0,
                 "raw_response_or_real_prediction_payloads_included": False,
                 "models_checked": 0, "tables_checked": 0,
                 "metadata_csvs_checked": 0,
                 "generated_expression_input_only": True}

    def add(target, package, name, data, *, source=None, columns=None, count=None,
            byte_identical=False):
        require(not Path(name).is_absolute() and ".." not in Path(name).parts,
                "Archive member path must be relative")
        target[name] = data
        entry = {"archive": package, "path": name, "sha256": sha(data), "bytes": len(data),
                 "byte_identical_source_copy": byte_identical}
        if source:
            entry.update(source_relative_path="reproducibility_v1/local_assets/" + source,
                         source_sha256=opened[source]["sha256"])
        if columns is not None:
            entry["columns"] = columns
        if count is not None:
            entry["data_rows"] = count
        catalog.append(entry)

    for relative in sorted(MODEL_FILES):
        data = frozen(relative)
        expected = model_records[relative]
        require(sha(data) == expected["portable_sha256"] and len(data) == expected["portable_bytes"],
                "Portable model differs from model export manifest")
        plain = gzip.decompress(data)
        require(sha(plain) == expected["canonical_json_sha256"], "Canonical model hash differs")
        field_audit(json.loads(plain))
        add(model_entries, MODEL_ARCHIVE, relative, data, source=relative, byte_identical=True)
        preservation.append({"path": relative, "source_sha256": sha(data),
                             "packaged_sha256": sha(data), "byte_identical": True})
        sensitive["models_checked"] += 1

    for relative, expected_count in sorted(TABLES.items()):
        data = frozen(relative)
        plain = gzip.decompress(data) if relative.endswith(".gz") else data
        reader = csv.DictReader(io.StringIO(plain.decode("utf-8-sig"), newline=""))
        columns = reader.fieldnames or []
        require(not set(columns) & FORBIDDEN_FIELDS, "Forbidden table column")
        count = 0
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()), "Malformed parameter table")
            field_audit(row)
            count += 1
        require(count == expected_count, "Parameter table row count differs")
        add(model_entries, MODEL_ARCHIVE, relative, data, source=relative, columns=columns,
            count=count, byte_identical=True)
        sensitive["tables_checked"] += 1

    for relative in sorted(SCHEMA_FILES):
        data = frozen(relative)
        field_audit(json.loads(data))
        add(model_entries, MODEL_ARCHIVE, relative, data, source=relative, byte_identical=True)
    data = frozen(SYNTHETIC)
    generated = json.loads(data)
    field_audit(generated, allow_generated_tpm=True)
    gene_order = json.loads(model_entries["schemas/input_gene_ids.json"])
    require(generated["gene_ids"] == gene_order and len(gene_order) == 19836,
            "Generated input complete gene order differs")
    n = len(generated["sample_ids"])
    require(n == 6 and all(len(generated[k]) == n for k in
                           ["biological_group_ids", "tissues", "raw_tpm"]) and
            all(len(row) == len(gene_order) for row in generated["raw_tpm"]),
            "Generated full-gene input dimensions differ")
    add(model_entries, MODEL_ARCHIVE, SYNTHETIC, data, source=SYNTHETIC, count=n,
        byte_identical=True)

    for cohort, expected_ids, expected_members in [("development", 626, 4967), ("external", 175, 1185)]:
        for kind, wanted, expected_count in [("identities", IDENTITY_FIELDS, expected_ids),
                                              ("drug_membership", MEMBERSHIP_FIELDS, expected_members)]:
            relative = f"metadata/{cohort}_{kind}.csv"
            original = frozen(relative)
            fields, rows = rows_of(original)
            require(len(rows) == expected_count and set(wanted).issubset(fields),
                    "Frozen metadata coverage/schema differs")
            projected = [{name: row[name] for name in wanted} for row in rows]
            require(all(all(row[name] for name in IDENTITY_FIELDS) for row in projected),
                    "Retained identity/context field is empty")
            for row in projected:
                field_audit(row)
            encoded = csv_encode(wanted, projected)
            check_fields, check_rows = rows_of(encoded)
            require(check_fields == wanted and check_rows == projected,
                    "Projection changed retained field strings or row order")
            if kind == "identities":
                require(len({row["canonical_model_id"] for row in rows}) == expected_ids and
                        len({row["biological_group_id"] for row in rows}) == expected_ids,
                        "Identity/group uniqueness differs")
                cohort_ids[cohort] = {row["canonical_model_id"]: row for row in projected}
            else:
                require(len({(row["gdsc_drug_id"], row["canonical_model_id"]) for row in rows}) == expected_members,
                        "Duplicate frozen drug membership")
                require({int(row["gdsc_drug_id"]) for row in rows} == DRUGS,
                        "Drug panel differs")
                require(set(row["canonical_model_id"] for row in rows) == set(cohort_ids[cohort]),
                        "Membership and identity universes differ")
                require(all(all(row[name] == cohort_ids[cohort][row["canonical_model_id"]][name]
                                for name in IDENTITY_FIELDS) for row in projected),
                        "Membership/identity retained fields disagree")
                memberships[cohort] = {"identities": expected_ids, "memberships": len(rows),
                                       "drug_membership_counts": {str(drug): sum(int(row["gdsc_drug_id"]) == drug
                                                                            for row in rows) for drug in sorted(DRUGS)}}
            projection.append({"file": relative, "source_sha256": sha(original),
                               "projection_sha256": sha(encoded), "rows_before": len(rows),
                               "rows_after": len(check_rows), "retained_columns": wanted,
                               "omitted_columns": [name for name in fields if name not in wanted],
                               "row_order_and_retained_strings_preserved": True})
            add(metadata_entries, METADATA_ARCHIVE, relative.removeprefix("metadata/"),
                encoded, source=relative, columns=wanted, count=len(rows))
            sensitive["metadata_csvs_checked"] += 1

    require(not set(cohort_ids["development"]) & set(cohort_ids["external"]),
            "Development/external canonical overlap")
    require(not {row["biological_group_id"] for row in cohort_ids["development"].values()} &
            {row["biological_group_id"] for row in cohort_ids["external"].values()},
            "Development/external group overlap")
    for name, expected_count in FOLDS.items():
        relative = "metadata/" + name
        data = frozen(relative)
        fields, rows = rows_of(data)
        require(set(fields).issubset(FOLD_FIELDS) and not set(fields) & FORBIDDEN_FIELDS,
                "Frozen assignment has an unexpected/source pointer field")
        require(len(rows) == expected_count, "Frozen assignment row count differs")
        for row in rows:
            field_audit(row)
            require(row["canonical_model_id"] in cohort_ids["development"] and
                    all(row[key] == cohort_ids["development"][row["canonical_model_id"]][key]
                        for key in ["biological_group_id", "tissue"]),
                    "Assignment identity/group/tissue disagrees with frozen development roster")
        add(metadata_entries, METADATA_ARCHIVE, name, data, source=relative, columns=fields,
            count=len(rows), byte_identical=True)
        fold_audits.append({"file": name, "rows": len(rows), "sha256": sha(data),
                            "byte_identical": True, "source_pointer_columns_present": False})
        sensitive["metadata_csvs_checked"] += 1

    schema = {"schema_version": "reo-public-minimal-metadata-1.0", "release_version": VERSION,
              "identity_key": ["canonical_model_id"], "group_key": ["biological_group_id"],
              "drug_membership_key": ["gdsc_drug_id", "canonical_model_id"],
              "identity_columns": IDENTITY_FIELDS, "membership_columns": MEMBERSHIP_FIELDS,
              "expression_row_index_semantics": "Zero-based project-prepared matrix offset; not a portable DepMap row number",
              "frozen_assignments_regenerated": False,
              "complete_input_for_historical_external_rescoring": False,
              "upstream_inputs_must_be_obtained_independently": True,
              "historical_identities_fully_reconstructible_from_DepMap_24Q4": False,
              "third_party_identifiers_relicensed": False}
    add(metadata_entries, METADATA_ARCHIVE, "mapping_schema.json", encode(schema))

    for entries, package, readme in [(model_entries, MODEL_ARCHIVE, "MODEL_PACKAGE_README.md"),
                                     (metadata_entries, METADATA_ARCHIVE, "METADATA_PACKAGE_README.md")]:
        for source, target in [(readme, "README.md"), ("MODEL_ASSETS_LICENSE.md", "MODEL_ASSETS_LICENSE.md"),
                               ("THIRD_PARTY_ATTRIBUTION.md", "THIRD_PARTY_ATTRIBUTION.md")]:
            data = (DOCS / source).read_bytes()
            require(not PRIVATE_PATTERN.search(data.decode("utf-8")), "Notice has a private pattern")
            add(entries, package, target, data)
        included = [row for row in catalog if row["archive"] == package]
        add(entries, package, "ASSET_MANIFEST.json", encode({"release_version": VERSION,
            "manifest_scope": "All other entries in this archive; outer release manifest hashes this file",
            "files": included}))

    archives = {MODEL_ARCHIVE: archive(model_entries), METADATA_ARCHIVE: archive(metadata_entries)}
    archive_records = [{"file": name, "bytes": len(data), "sha256": sha(data),
                        "entries": len(model_entries if name == MODEL_ARCHIVE else metadata_entries)}
                       for name, data in sorted(archives.items())]
    require(snapshot() == before, "A local_assets source file changed during construction")
    audit = {"release_version": VERSION, "projections": projection, "folds": fold_audits,
             "cohorts": memberships, "canonical_overlap": 0, "biological_group_overlap": 0,
             "local_assets_file_size_and_mtime_snapshot_unchanged": True,
             "read_allowlist_enforced": True, "read_files": opened,
             "no_training_no_scoring_no_prediction_execution": True}
    summary = {"release_version": VERSION, "archives": archive_records,
               "model_count": 64, "model_payloads_byte_identical": True,
               "molecular_models": export["molecular_models"],
               "zero_molecular_models": export["zero_molecular_models"],
               "nonzero_molecular_models": export["nonzero_molecular_models"],
               "tissue_reference_models": export["tissue_reference_models"],
               "cohorts": memberships, "frozen_assignment_files": 4,
               "raw_response_values_included": False, "real_individual_predictions_included": False,
               "original_expression_observations_included": False,
               "complete_input_for_historical_external_rescoring": False,
               "third_party_identifiers_relicensed": False,
               "license_scope": "CC BY 4.0 only to the extent the authors own the relevant rights; upstream terms remain",
               "source_assets_modified": False, "training_scoring_or_prediction_executed": False}
    catalog_object = {"release_version": VERSION, "archives": archive_records,
                      "contains_identity_rows_or_fitted_vectors": False,
                      "files": sorted(catalog, key=lambda row: (row["archive"], row["path"]))}
    output = {**archives, "public_release_summary.json": encode(summary),
              "asset_catalog.json": encode(catalog_object),
              "release_file_manifest.json": encode({"release_version": VERSION,
                  "archives": archive_records, "archive_entries": catalog_object["files"]}),
              "projection_membership_audit.json": encode(audit),
              "sensitive_pattern_audit.json": encode(sensitive),
              "model_hash_preservation.json": encode({"release_version": VERSION,
                  "model_count": 64, "all_byte_identical": True, "files": preservation})}
    output["SHA256SUMS"] = "".join(f"{sha(data)}  {name}\n" for name, data in sorted(output.items())).encode("ascii")
    out.mkdir(parents=True, exist_ok=False)
    for name, data in sorted(output.items()):
        with (out / name).open("xb") as handle:
            handle.write(data)
    require(snapshot() == before, "A local_assets source file changed during output writing")
    print(json.dumps({"release_version": VERSION, "output_relative_to_repo": out.relative_to(ROOT).as_posix(),
                      "archives": archive_records, "model_hashes_preserved": 64,
                      "metadata_projection_files": 4, "unchanged_fold_files": 4,
                      "source_assets_modified": False}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
