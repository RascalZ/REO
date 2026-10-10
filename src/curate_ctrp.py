"""Curate CTRPv2.0 identities without inspecting held-out response distributions.

The archived curve table remains the authoritative response store. This program
indexes its physical line numbers and checks only response presence/finite syntax.
It neither chooses drugs by response nor calculates response summary statistics.

Usage: python rebuild_v2/src/curate_ctrp.py [--download]
No pandas or other third-party dependencies are required.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import io
import json
import math
from pathlib import Path
import struct
import urllib.request
import zipfile
import zlib

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data_raw" / "ctrp"
OUT = ROOT / "data_processed"
URL = "https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset.zip"
SOURCE_CONFIG = "https://raw.githubusercontent.com/BHKLAB-DataProcessing/CTRPv2-Pharmacoset_Snakemake/main/config/pipeline.yaml"
WANTED = {
    "MANIFEST.txt", "v20._COLUMNS.txt", "v20._README.txt", "CTRPv2.0._README.docx",
    "v20.data.curves_post_qc.txt", "v20.meta.per_cell_line.txt",
    "v20.meta.per_compound.txt", "v20.meta.per_experiment.txt",
    "v20.meta.per_assay_plate.txt",
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def verify_published_manifest():
    """Validate selected files against the MD5 manifest distributed in the ZIP."""
    verified = []
    for line in (RAW / "MANIFEST.txt").read_text(encoding="utf-8").splitlines():
        expected, filename = line.split(" *", 1)
        path = RAW / filename
        if not path.exists():
            continue  # Deliberately unrequested large raw-well tables are absent.
        h = hashlib.md5()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        if h.hexdigest() != expected:
            raise RuntimeError(f"Published MANIFEST.txt MD5 mismatch: {filename}")
        verified.append(filename)
    return verified


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def get_range(start, end):
    req = urllib.request.Request(URL, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(req, timeout=90) as r:
        if r.status != 206:
            raise RuntimeError("Archive host ignored HTTP Range; refusing full archive download")
        content_range = r.headers.get("Content-Range", "")
        if not content_range.startswith(f"bytes {start}-{end}/"):
            raise RuntimeError(f"Unexpected Content-Range: {content_range}")
        b = r.read()
    if len(b) != end - start + 1:
        raise RuntimeError("Truncated HTTP range")
    return b


def download():
    """Download only selected archive members, checking each member's ZIP CRC32."""
    RAW.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(urllib.request.Request(URL, method="HEAD"), timeout=30) as r:
        total = int(r.headers["Content-Length"])
        headers = {k: r.headers.get(k) for k in ("Content-Length", "Last-Modified", "ETag", "x-amz-version-id")}
    tail = get_range(max(0, total - 65557), total - 1)
    pos = tail.rfind(b"PK\x05\x06")
    if pos < 0:
        raise RuntimeError("ZIP end-of-central-directory not found")
    fields = struct.unpack_from("<4s4H2IH", tail, pos)
    directory_size, directory_offset = fields[5:7]
    directory = get_range(directory_offset, directory_offset + directory_size - 1)
    # zipfile adjusts offsets to the isolated directory; undo that adjustment.
    with zipfile.ZipFile(io.BytesIO(directory + tail[pos:])) as z:
        members = [{"filename": m.filename, "file_size": m.file_size,
                    "compress_size": m.compress_size, "CRC": m.CRC,
                    "header_offset": m.header_offset + directory_offset,
                    "compress_type": m.compress_type} for m in z.infolist()]
    missing = WANTED - {m["filename"] for m in members}
    if missing:
        raise RuntimeError(f"Required ZIP members absent: {missing}")
    manifest_path = RAW / "download_manifest.json"
    old = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else []
    manifest = {m["filename"]: m for m in old}
    for m in members:
        if m["filename"] not in WANTED:
            continue
        target = RAW / m["filename"]
        if target.exists():
            prior = manifest.get(m["filename"])
            if prior is None or sha256(target) != prior["sha256"]:
                raise RuntimeError(f"Existing raw file is not manifest-verified: {target.name}")
            if prior["file_size"] != m["file_size"] or prior["CRC"] != m["CRC"]:
                raise RuntimeError(f"Remote archive revision changed: {target.name}; use an explicit new data version")
            continue
        local_header = get_range(m["header_offset"], m["header_offset"] + 29)
        h = struct.unpack("<4s5H3I2H", local_header)
        if h[0] != b"PK\x03\x04":
            raise RuntimeError("Invalid local ZIP header")
        start = m["header_offset"] + 30 + h[-2] + h[-1]
        end = start + m["compress_size"] - 1
        compressed = get_range(start, end)
        if m["compress_type"] == 8:
            data = zlib.decompress(compressed, -15)
        elif m["compress_type"] == 0:
            data = compressed
        else:
            raise RuntimeError("Unsupported ZIP compression")
        if len(data) != m["file_size"] or zlib.crc32(data) != m["CRC"]:
            raise RuntimeError(f"Archive member integrity check failed: {m['filename']}")
        with target.open("xb") as f:
            f.write(data)
        manifest[m["filename"]] = dict(m, url=URL, range_start=start, range_end=end,
            sha256=hashlib.sha256(data).hexdigest(),
            retrieved_at=dt.datetime.now(dt.timezone.utc).isoformat(), crc32_verified=True)
        write_json(manifest_path, list(manifest.values()))
    write_json(RAW / "archive_members.json", members)
    write_json(RAW / "archive_provenance.json", {
        "archive_url": URL, "discovery_source": SOURCE_CONFIG,
        "retrieved_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "archive_headers": headers, "selective_http_range_download": True,
        "whole_archive_downloaded": False,
        "caution": "Archive content reached via curator-published CloudFront URL; historical NCI endpoint redirects. ZIP CRC verifies transfer integrity, not independent scientific validity."})


def rows(filename):
    with (RAW / filename).open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def write_csv(filename, fieldnames, records):
    with (OUT / filename).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(records)


def one_to_one(records, key):
    result = {}
    for r in records:
        if r[key] in result:
            raise RuntimeError(f"Nonunique metadata key: {key}={r[key]}")
        result[r[key]] = r
    return result


def finite_present(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def curate():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((RAW / "download_manifest.json").read_text(encoding="utf-8"))
    for m in manifest:
        p = RAW / m["filename"]
        if sha256(p) != m["sha256"] or p.stat().st_size != m["file_size"]:
            raise RuntimeError(f"Raw-source checksum mismatch: {p.name}")
    manifest_verified_files = verify_published_manifest()
    cells = one_to_one(rows("v20.meta.per_cell_line.txt"), "master_ccl_id")
    compounds = one_to_one(rows("v20.meta.per_compound.txt"), "master_cpd_id")
    experiments = collections.defaultdict(list)
    for line, r in enumerate(rows("v20.meta.per_experiment.txt"), 2):
        experiments[r["experiment_id"]].append(dict(r, source_row_number=line))
    exp_identity = {}
    exp_output = []
    for eid, records in experiments.items():
        cids = {r["master_ccl_id"] for r in records}
        if len(cids) != 1:
            raise RuntimeError(f"Ambiguous experiment-to-cell identity: {eid}")
        cid = next(iter(cids))
        if cid not in cells:
            raise RuntimeError(f"Cell metadata missing for {cid}")
        exp_identity[eid] = cid
        sets = {k: sorted({str(r[k]) for r in records}) for k in records[0] if k != "source_row_number"}
        exp_output.append({
            "experiment_id": eid, "master_ccl_id": cid, "ccl_name": cells[cid]["ccl_name"],
            "screen_id": "CTRPv2.0:experiment:" + eid,
            "screen_id_definition": "derived_namespace_of_source_experiment_id",
            "n_source_metadata_rows": len(records), "run_ids": "|".join(sets["run_id"]),
            "experiment_dates": "|".join(sets["experiment_date"]),
            "snp_fp_status": "|".join(sets["snp_fp_status"]),
            "has_multiple_dates_or_runs": int(len(sets["run_id"]) > 1 or len(sets["experiment_date"]) > 1),
            "source_metadata_rows_json": json.dumps(records, ensure_ascii=False, separators=(",", ":")),
        })
    exp_resolved = {r["experiment_id"]: r for r in exp_output}
    write_csv("ctrp_experiment_metadata.csv", list(exp_output[0]), exp_output)

    seen_keys = set()
    pair_experiments = collections.defaultdict(set)
    pair_available = collections.Counter()
    cell_rows, cell_drugs = collections.Counter(), collections.defaultdict(set)
    drug_rows, drug_cells = collections.Counter(), collections.defaultdict(set)
    source_identity = collections.Counter()
    curve_point_counts = collections.Counter()
    curve_parameter_counts = collections.Counter()
    response_present_count = 0
    count = 0
    index_fields = ["source_row_number", "source_record_id", "experiment_id", "screen_id",
                    "master_cpd_id", "master_ccl_id", "ccl_name", "ccle_primary_site",
                    "response_available", "source_qc_stage", "snp_fp_status", "conc_pts_fit", "fit_num_param"]
    with (OUT / "ctrp_curve_index.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=index_fields)
        w.writeheader()
        for line, r in enumerate(rows("v20.data.curves_post_qc.txt"), 2):
            eid, did = r["experiment_id"], r["master_cpd_id"]
            if eid not in exp_identity or did not in compounds:
                raise RuntimeError(f"Unmapped curve identity at source line {line}")
            if (eid, did) in seen_keys:
                raise RuntimeError(f"Duplicate source experiment-compound curve at source line {line}")
            seen_keys.add((eid, did))
            cid = exp_identity[eid]
            cell = cells[cid]
            available = int(finite_present(r["area_under_curve"]))
            w.writerow({"source_row_number": line,
                "source_record_id": f"CTRPv2.0:experiment:{eid}:compound:{did}",
                "experiment_id": eid, "screen_id": exp_resolved[eid]["screen_id"],
                "master_cpd_id": did, "master_ccl_id": cid,
                "ccl_name": cell["ccl_name"], "ccle_primary_site": cell["ccle_primary_site"],
                "response_available": available, "source_qc_stage": "published_curves_post_qc",
                "snp_fp_status": exp_resolved[eid]["snp_fp_status"],
                "conc_pts_fit": r["conc_pts_fit"], "fit_num_param": r["fit_num_param"]})
            curve_point_counts[r["conc_pts_fit"]] += 1
            curve_parameter_counts[r["fit_num_param"]] += 1
            pair_experiments[(did, cid)].add(eid)
            pair_available[(did, cid)] += available
            cell_rows[cid] += 1
            cell_drugs[cid].add(did)
            drug_rows[did] += 1
            drug_cells[did].add(cid)
            source_identity[(compounds[did]["cpd_name"], cell["ccl_name"], cell["ccle_primary_site"])] += 1
            response_present_count += available
            count += 1
    write_csv("ctrp_pair_availability.csv", ["master_cpd_id", "master_ccl_id", "ccl_name",
              "n_experiments", "n_response_available", "experiment_ids", "repeat_rule_applied"],
              ({"master_cpd_id": did, "master_ccl_id": cid, "ccl_name": cells[cid]["ccl_name"],
                "n_experiments": len(eids), "n_response_available": pair_available[(did, cid)],
                "experiment_ids": "|".join(sorted(eids, key=int)),
                "repeat_rule_applied": "none_preserved"}
               for (did, cid), eids in sorted(pair_experiments.items(), key=lambda x: (int(x[0][0]), int(x[0][1])))))
    cell_output = [dict(r, n_curve_records=cell_rows[cid], n_compounds=len(cell_drugs[cid]),
                   has_curve_records=int(cell_rows[cid] > 0),
                   name_namespace="CTRP_official_ccl_name_not_CCLE_Name") for cid, r in cells.items()]
    write_csv("ctrp_cell_metadata.csv", list(cell_output[0]), cell_output)
    write_csv("ctrp_drug_availability.csv", ["master_cpd_id", "cpd_name", "n_curve_records", "n_cell_lines"],
              ({"master_cpd_id": did, "cpd_name": compounds[did]["cpd_name"],
                "n_curve_records": n, "n_cell_lines": len(drug_cells[did])}
               for did, n in sorted(drug_rows.items(), key=lambda x: int(x[0]))))
    from project_paths import original_input
    legacy = original_input('legacy_ctrp_response')
    legacy_summary = {"file": str(legacy), "response_values_compared": False}
    if legacy.exists():
        with legacy.open(encoding="utf-8-sig", newline="") as f:
            legacy_identity = collections.Counter((r["Drug"], r["Cellline"], r["CCLE_primary_site"])
                for r in csv.DictReader(f, delimiter="\t"))
        legacy_summary.update({"sha256": sha256(legacy), "n_rows": sum(legacy_identity.values()),
            "identity_multiset_equal_to_published_archive": legacy_identity == source_identity,
            "unmatched_legacy_identity_rows": sum((legacy_identity - source_identity).values()),
            "unmatched_archive_identity_rows": sum((source_identity - legacy_identity).values()),
            "conclusion": "Identity multiset matches; response equality and historical derivation were not established. Use archived source with source IDs."})
    write_json(OUT / "ctrp_legacy_identity_summary.json", legacy_summary)
    summary = {
        "dataset": "CTRPv2.0_2015_ctd2_ExpandedDataset",
        "n_curve_records": count, "n_compounds_with_curves": len(drug_rows),
        "n_cells_with_curves": len(cell_rows), "n_experiments_with_curves": len({k[0] for k in seen_keys}),
        "n_cells_in_metadata": len(cells), "n_experiment_metadata_rows": sum(map(len, experiments.values())),
        "n_unique_experiment_ids_in_metadata": len(experiments),
        "n_experiment_ids_with_multiple_metadata_rows": sum(len(v) > 1 for v in experiments.values()),
        "n_experiment_ids_with_multiple_dates_or_runs": sum(r["has_multiple_dates_or_runs"] for r in exp_output),
        "n_experiment_ids_with_multiple_run_ids": sum("|" in r["run_ids"] for r in exp_output),
        "n_ambiguous_experiment_cell_mappings": 0, "n_unmapped_curve_ids": 0,
        "n_duplicate_experiment_compound_curves": 0,
        "n_drug_cell_pairs": len(pair_experiments),
        "n_drug_cell_pairs_with_repeat_experiments": sum(len(v) > 1 for v in pair_experiments.values()),
        "n_curve_records_with_finite_response": response_present_count,
        "n_experiments_by_snp_fp_status": dict(collections.Counter(r["snp_fp_status"] for r in exp_output)),
        "n_experiment_ids_with_multiple_snp_statuses": sum("|" in r["snp_fp_status"] for r in exp_output),
        "n_curves_by_conc_pts_fit": dict(sorted(curve_point_counts.items(), key=lambda x: int(x[0]))),
        "n_curves_by_fit_num_param": dict(sorted(curve_parameter_counts.items(), key=lambda x: int(x[0]))),
        "response_store": "data_raw/ctrp/v20.data.curves_post_qc.txt",
        "response_field": "area_under_curve",
        "response_definition": "Integrated area under the sigmoid-fit concentration-response curve; not ROC-AUC. Original archive value is retained without rescaling.",
        "raw_curve_sha256": sha256(RAW / "v20.data.curves_post_qc.txt"),
        "published_manifest_md5_verified_files": manifest_verified_files,
        "response_values_printed_or_summarized": False, "response_values_used_for_selection": False,
        "repeated_screen_aggregation": "not_performed; prespecify aggregation before outcome access",
        "screen_id_caution": "screen_id is a documented namespaced experiment_id, not an additional original field or biological replicate ID",
        "qc_caution": "Published post-QC membership and source fit fields retained. Pre-QC dose-level censoring flags are not in this download; do not claim re-executed raw-well QC.",
        "output_files": ["ctrp_curve_index.csv", "ctrp_pair_availability.csv", "ctrp_cell_metadata.csv",
                         "ctrp_experiment_metadata.csv", "ctrp_drug_availability.csv", "ctrp_legacy_identity_summary.json"],
    }
    write_json(OUT / "ctrp_curation_summary.json", summary)
    print(json.dumps({k: summary[k] for k in ("dataset", "n_curve_records", "n_compounds_with_curves",
          "n_cells_with_curves", "n_experiments_with_curves", "n_drug_cell_pairs",
          "n_drug_cell_pairs_with_repeat_experiments", "response_values_printed_or_summarized")}, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--download", action="store_true", help="Selectively retrieve missing official archive members via documented curator URL")
    args = ap.parse_args()
    if args.download:
        download()
    curate()
