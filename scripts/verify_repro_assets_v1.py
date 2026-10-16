"""Inspect inference asset hashes/contracts without fitting or outcome scoring.

Default: public counts/fixtures only. --local also checks explicitly exported
local numeric bundles and metadata, never source labels/indices/predictions.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).absolute().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def check(condition, message):
    if not condition:
        raise ValueError(message)


def no_parameter_or_observation_fields(obj):
    forbidden = {"input_gene_ids", "candidate_gene_ids", "candidate_gene_indices", "coefficients", "feature_mean", "feature_scale", "feature_tissue_means", "target_tissue_means", "fit_sample_ids", "predictions", "responses", "labels", "raw_tpm"}
    if isinstance(obj, dict):
        return not (set(obj) & forbidden) and all(no_parameter_or_observation_fields(v) for v in obj.values())
    if isinstance(obj, list):
        return all(no_parameter_or_observation_fields(v) for v in obj)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local", action="store_true")
    args = parser.parse_args()
    docs = ROOT / "docs/reproducibility_assets_v1"
    models = read(docs / "model_export_public_summary.json")
    meta = read(docs / "metadata_export_public_summary.json")
    inventory = read(docs / "asset_bundle_manifest.json")
    check(inventory["models"] == 64 and inventory["complete_local_assets_git_ignored"] is True, "Local asset availability status differs")
    check(models["model_count"] == 64 and models["zero_molecular_models"] == 29 and models["nonzero_molecular_models"] == 27 and models["tissue_reference_models"] == 8, "Final model counts differ")
    check(models["synthetic_numeric_checks"] == 384 and models["max_absolute_native_portable_difference"] == 0.0, "Recorded synthetic numerical agreement differs")
    check(models["new_training"] is False and models["metric_or_model_selection_scoring"] is False and models["CTRP_labels_indices_or_individual_predictions_read"] is False and not models["blocked_accesses"], "Model export access/count contract differs")
    check(no_parameter_or_observation_fields(models) and no_parameter_or_observation_fields(meta), "Public summaries must not contain fitted vectors or observation fields")
    check(meta["counts"]["development"]["identities"] == 626 and meta["counts"]["development"]["drug_memberships"] == 4967 and meta["counts"]["external"]["identities"] == 175 and meta["counts"]["external"]["drug_memberships"] == 1185, "Metadata counts differ")
    expected = read(ROOT / "examples/repro_inference_v1/toy_expected_output.json")
    check(expected["statuses"].count("predicted") == 4 and expected["statuses"].count("abstained") == 3, "Toy recorded cases differ")
    verified_models = verified_metadata = verified_tables = 0
    if args.local:
        base = (ROOT / "reproducibility_v1/local_assets").resolve()
        check(read(base / "assets_manifest_v1.json") == inventory, "Local/public bundle hash inventories differ")
        for row in inventory["files"]:
            p = (base / row["path"]).resolve()
            check(p.is_relative_to(base) and digest(p) == row["sha256"], "Exported local bundle asset SHA differs")
        check(read(base / "model_export_manifest.json") == models, "Local/public model-export manifests differ")
        gene_order = read(base / "schemas/input_gene_ids.json")
        for row in models["models"]:
            path = (base / row["path"]).resolve()
            check(path.is_relative_to(base / "models") and path.name.endswith(".json.gz"), "Unexpected model path")
            check(digest(path) == row["portable_sha256"], "Compressed model SHA differs")
            data = gzip.decompress(path.read_bytes())
            check(hashlib.sha256(data).hexdigest() == row["canonical_json_sha256"], "Canonical model SHA differs")
            bundle = json.loads(data)
            check(bundle["input_gene_ids"] == gene_order and len(gene_order) == 19836, "Full gene order differs")
            check(bundle["n_nonzero_molecular_coeffs"] == row["nonzero_coefficients"] == sum(v != 0.0 for v in bundle["coefficients"]), "Zero/nonzero coefficient counts differ")
            verified_models += 1
        for row in models["tables"]:
            path = (base / row["path"]).resolve()
            check(path.is_relative_to(base / "tables"), "Unexpected table path")
            check(digest(path) == row["sha256"], "Fitted table SHA differs")
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", encoding="utf-8", newline="") as f:
                count = sum(1 for _ in csv.reader(f)) - 1
            check(count == row["rows"], "Fitted table row count differs")
            verified_tables += 1
        md = base / "metadata"
        check(digest(md / "metadata_export_manifest.json") == meta["manifest_sha256"], "Local metadata manifest SHA differs")
        manifest = read(md / "metadata_export_manifest.json")
        for row in manifest["outputs"]:
            path = (md / row["relative_path"]).resolve()
            check(path.is_relative_to(md) and digest(path) == row["sha256"], "Local metadata output SHA differs")
            check(row.get("response_value_fields_present") is not True, "Unexpected response-bearing metadata output")
            verified_metadata += 1
    print(json.dumps({"status": "passed", "mode": "local_exported_assets" if args.local else "public_summary_and_fixture_contracts", "verified_numeric_model_bundles": verified_models, "verified_local_metadata_outputs": verified_metadata, "verified_fitted_tables": verified_tables, "new_training_or_scoring": False, "source_CTRP_responses_indices_or_individual_predictions_read": False, "scope": "hash/schema/count checks; recorded generated-input agreement is not recomputed efficacy or independent study replication"}))


if __name__ == "__main__":
    main()
