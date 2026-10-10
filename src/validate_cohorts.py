"""Independently validate the frozen identity-only cohort; never inspect responses.

Usage: python rebuild_v2/src/validate_cohorts.py [--lock]

Reads original inputs only to hash bytes, and CTRP curve rows only to compare
identifiers. No response values are printed, summarized, selected, or evaluated.
The deterministic report is bound to the freeze file's hash. An existing lock or
its validation report is never replaced with inconsistent contents.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "data_processed"
R = ROOT / "results"
FREEZE = ROOT / "config" / "cohort_freeze_v1.json"
LOCK = ROOT / "config" / "cohort_freeze_v1.lock.json"
REPORT = ROOT / "reports" / "validation_summary.json"
HASH_CACHE = {}
CHECKS = []


class ValidationFailure(RuntimeError):
    pass


def require(condition, name, **evidence):
    CHECKS.append({"check": name, "passed": bool(condition), **evidence})
    if not condition:
        raise ValidationFailure(name)


def digest(path):
    path = Path(path).resolve()
    stat = path.stat()
    key = (str(path), stat.st_size, stat.st_mtime_ns)
    if key not in HASH_CACHE:
        h = hashlib.sha256()
        with path.open("rb") as f:
            for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
                h.update(b)
        HASH_CACHE[key] = h.hexdigest()
    return HASH_CACHE[key]


def json_bytes(obj):
    return (json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def iter_csv(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f, delimiter=delimiter)


def table(path, delimiter=","):
    return list(iter_csv(path, delimiter))


def unique_map(records, key, label):
    keys = [r[key] for r in records]
    require(all(keys) and len(keys) == len(set(keys)), label + "_unique_nonempty_key", n_rows=len(keys))
    return {r[key]: r for r in records}


def as_set(value):
    return set(value.split("|")) - {""}


def yes(value):
    return str(value).lower() in {"1", "true"}


def verify_file(path, item, category):
    require(path.is_file(), category + "_exists", path=str(path))
    require(digest(path) == item["sha256"], category + "_sha256", path=str(path))
    if "bytes" in item:
        require(path.stat().st_size == int(item["bytes"]), category + "_bytes", path=str(path))


def no_response_columns(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        fields = next(csv.reader(f))
    forbidden = {"auc", "areaundercurve", "lnic50", "ic50", "aac", "response",
                 "responsevalue", "responsemagnitude", "observedresponse", "prediction",
                 "predictedresponse", "y", "ytrue", "ypred", "label", "classlabel",
                 "resistant", "sensitive", "sensitivityscore"}
    found = [x for x in fields if re.sub(r"[^a-z0-9]", "", x.lower()) in forbidden]
    require(not found, "identity_artifact_has_no_response_value_columns", path=str(path.relative_to(ROOT)), prohibited_headers=found)


def validate(freeze):
    require(freeze.get("response_magnitudes_accessed_for_selection") is False,
            "freeze_declares_no_response_magnitude_selection")
    require(freeze.get("training_configuration_frozen") is False,
            "cohort_lock_does_not_claim_model_hyperparameter_freeze")
    require(freeze["design"] == "natural_CTRP_only_known_biological_groups",
            "validator_targets_natural_GDSC2_disjoint_design")

    frozen_paths = set()
    for item in freeze["files"]:
        if item.get("scope") == "original_read_only":
            path = Path(item["path"])
            require(path.is_absolute(), "original_read_only_path_is_absolute", path=str(path))
        else:
            path = ROOT / item["path"]
            require(path.resolve().is_relative_to(ROOT.resolve()), "frozen_relative_dependency_stays_in_root", path=item["path"])
        require(str(path.resolve()) not in frozen_paths, "freeze_dependency_not_duplicated", path=str(path))
        frozen_paths.add(str(path.resolve()))
        verify_file(path, item, "frozen_dependency")
    # Ensure the expanded freeze covers its implementation and all substantive inputs.
    must_cover = [Path(__file__).resolve(), ROOT / "src/curate_cells.py", ROOT / "src/map_ctrp_cells.py",
                  ROOT / "src/freeze_cohorts.py", P / "gdsc2_cell_mapping.csv", P / "ctrp_cell_mapping.csv",
                  P / "canonical_model_groups.csv", P / "biological_group_edges.csv", P / "drug_identity_eligible.csv",
                  P / "expression_audit.json", ROOT / "data_manifest/original_source_sha256.csv",
                  ROOT / "data_raw/ctrp/v20.data.curves_post_qc.txt"]
    require(all(str(p.resolve()) in frozen_paths for p in must_cover), "freeze_covers_critical_code_and_mapping_dependencies")

    originals = table(ROOT / "data_manifest/original_source_sha256.csv")
    expression_audit = load_json(P / "expression_audit.json")
    for category, records in [("original_manifest", originals),
                              ("expression_original_input", expression_audit["source_inputs"]),
                              ("expression_reference", expression_audit.get("downloaded_references", [])),
                              ("expression_derived_file", expression_audit.get("derived_files", []))]:
        for item in records:
            path = Path(item["path"])
            require(path.is_absolute(), category + "_path_is_absolute")
            verify_file(path, item, category)
            require(str(path.resolve()) in frozen_paths, category + "_included_in_freeze", path=str(path))
    for item in load_json(ROOT / "data_raw/ctrp/download_manifest.json"):
        path = ROOT / "data_raw/ctrp" / item["filename"]
        verify_file(path, item, "CTRP_download_member")
        require(path.stat().st_size == int(item["file_size"]), "CTRP_download_member_bytes", filename=item["filename"])

    groups = unique_map(table(P / "canonical_model_groups.csv"), "canonical_model_id", "canonical_group_map")
    expression = table(P / "cell_line_dictionary.csv")
    expr = unique_map([r for r in expression if r["expression_eligible"] == "1"], "canonical_model_id", "eligible_expression_model")
    sample_rows = unique_map(table(P / "expression_samples.tsv", "\t"), "expression_row_index", "expression_matrix_row")
    require(len(sample_rows) == int(expression_audit["matrix_shape"][0]), "expression_sample_table_matches_matrix_shape")
    for a, x in expr.items():
        require(a in groups and x["biological_group_id"] == groups[a]["biological_group_id"], "expression_group_matches_canonical_map", model_id=a)
        sr = sample_rows.get(x["expression_row_index"])
        require(sr is not None and sr["ccle_id"] == x["ccle_id"] and sr["depmap_id"] == a,
                "expression_matrix_row_accession_and_ACH_match", model_id=a)

    gm_rows = table(P / "gdsc2_cell_mapping.csv")
    gm_keys = [tuple(r[k] for k in ("sanger_model_id", "cosmic_id", "cell_line_name", "tcga_desc")) for r in gm_rows]
    require(len(gm_keys) == len(set(gm_keys)), "GDSC_mapping_join_key_unique")
    gm = dict(zip(gm_keys, gm_rows))
    cm = unique_map(table(P / "ctrp_cell_mapping.csv"), "master_ccl_id", "CTRP_mapping_join_key")
    all_gdsc_ids = {a for r in gm_rows for a in as_set(r["candidate_model_ids"])}
    require(all(a in groups for a in all_gdsc_ids), "all_GDSC2_candidates_have_group")
    all_gdsc_groups = {groups[a]["biological_group_id"] for a in all_gdsc_ids}
    require(len(all_gdsc_ids) == 976 and len(all_gdsc_groups) == 970,
            "expected_v1_all_GDSC2_candidate_identity_universe", n_ACH=len(all_gdsc_ids), n_groups=len(all_gdsc_groups))

    main = table(R / "main_drugs.csv")
    drugs = unique_map(main, "gdsc_drug_id", "primary_drug_assay")
    require(len(drugs) == 8 and len({r["ctrp_inchikey"] for r in main}) == 8, "eight_unique_primary_chemicals")
    require(sorted(int(r["primary_drug_rank"]) for r in main) == list(range(1, 9)), "primary_drug_ranks_1_to_8")
    eligible = unique_map(table(P / "drug_identity_eligible.csv"), "gdsc_drug_id", "eligible_drug_assay")
    chemical = unique_map(table(P / "drug_ctrp_compound_identity.csv"), "ctrp_compound_id", "CTRP_compound_identity")
    for d in main:
        did, cid = d["gdsc_drug_id"], d["ctrp_compound_id"]
        require(did in eligible and cid in chemical, "selected_drug_has_identity_evidence", gdsc_drug_id=did)
        key = d["ctrp_inchikey"]
        require(bool(re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", key)) and key == d["gdsc_name_derived_inchikey"]
                and key == chemical[cid]["ctrp_inchikey"] and key == eligible[did]["ctrp_inchikey"]
                and cid == eligible[did]["ctrp_compound_id"] and yes(eligible[did]["identity_eligible"])
                and not yes(chemical[cid]["is_combination"]),
                "selected_full_InChIKeys_agree_and_are_single_compound", gdsc_drug_id=did)
    frozen_panel = [{k: d[k] for k in ("gdsc_drug_id", "gdsc_drug_name", "ctrp_compound_id", "gdsc_min_conc", "gdsc_max_conc")} for d in main]
    require(frozen_panel == freeze["main_drugs"], "main_drug_panel_matches_freeze")

    datasets = {"development": table(R / "development_sample_manifest.csv"),
                "external_primary": table(R / "external_primary_sample_manifest.csv"),
                "external_sensitivity": table(R / "external_sensitivity_sample_manifest.csv")}
    for name in ("development_sample_manifest.csv", "external_primary_sample_manifest.csv", "external_sensitivity_sample_manifest.csv",
                 "main_drugs.csv", "global_split_manifest.csv", "all_295_GDSC_drug_sample_counts.csv"):
        no_response_columns(R / name)
    no_response_columns(P / "gdsc2_observation_index.csv")
    no_response_columns(P / "ctrp_curve_index.csv")
    for role, records in datasets.items():
        require(bool(records), role + "_not_empty")
        keys = [(r["gdsc_drug_id"], r["canonical_model_id"]) for r in records]
        require(len(keys) == len(set(keys)), role + "_unique_drug_model")
        require({r["gdsc_drug_id"] for r in records} == set(drugs), role + "_contains_exact_primary_drug_panel")
        for r in records:
            a = r["canonical_model_id"]
            require(a in expr and all(r[k] == expr[a][k] for k in ("ccle_id", "expression_row_index", "biological_group_id", "tissue")),
                    role + "_expression_identity_consistent", model_id=a, gdsc_drug_id=r["gdsc_drug_id"])
    dev, ext, sens = (datasets[k] for k in ("development", "external_primary", "external_sensitivity"))
    dids, dgroups = {r["canonical_model_id"] for r in dev}, {r["biological_group_id"] for r in dev}
    for label, records in (("primary", ext), ("sensitivity", sens)):
        ids = {r["canonical_model_id"] for r in records}
        gs = {r["biological_group_id"] for r in records}
        require(not (ids & all_gdsc_ids) and not (gs & all_gdsc_groups), label + "_external_absent_from_ALL_GDSC2_candidate_models_and_groups")
        require(not (ids & dids) and not (gs & dgroups), label + "_external_disjoint_from_selected_development")
        require(all(r["external_subtype"] == "natural_CTRP_only" for r in records), label + "_external_subtype_natural")
    require(all(r["identity_qc_statuses"] == "SNP-matched-reference" for r in ext), "primary_external_SNP_reference_matched_only")

    per_drug = []
    for did, d in drugs.items():
        dv = [r for r in dev if r["gdsc_drug_id"] == did]
        ex = [r for r in ext if r["gdsc_drug_id"] == did]
        ng = len({r["biological_group_id"] for r in dv})
        ne = len({r["biological_group_id"] for r in ex})
        nt = len({r["tissue"] for r in dv} - {""})
        require(ng >= 200 and ne >= 80 and nt >= 4, "per_drug_minimum_independent_coverage", gdsc_drug_id=did, development_groups=ng, external_groups=ne, development_tissues=nt)
        observed = {"n_development_models": len({r["canonical_model_id"] for r in dv}),
                    "n_development_groups": ng, "n_development_tissues": nt,
                    "n_external_primary_models": len({r["canonical_model_id"] for r in ex}),
                    "n_external_primary_groups": ne,
                    "n_external_primary_tissues": len({r["tissue"] for r in ex} - {""}),
                    "n_external_primary_curve_records": sum(len(as_set(r["ctrp_source_row_numbers"])) for r in ex)}
        require(all(int(d[k]) == v for k, v in observed.items()), "main_drug_counts_reconcile_with_manifests", gdsc_drug_id=did)
        require(all((r["min_conc"], r["max_conc"]) == (d["gdsc_min_conc"], d["gdsc_max_conc"]) for r in dv),
                "development_uses_one_exact_frozen_dose_interval", gdsc_drug_id=did)
        per_drug.append({"gdsc_drug_id": did, "drug_name": d["gdsc_drug_name"], **observed})

    # Validate GDSC mother-table joins and every selected source-row pointer.
    dev_by_source = unique_map(dev, "gdsc_source_row_number", "development_source_row")
    gdsc_drug_ids, gdsc_selected_found, gdsc_rows, unmatched = set(), set(), 0, 0
    for r in iter_csv(P / "gdsc2_observation_index.csv"):
        gdsc_rows += 1
        gdsc_drug_ids.add(r["gdsc_drug_id"])
        key = tuple(r[k] for k in ("sanger_model_id", "cosmic_id", "cell_line_name", "tcga_desc"))
        m = gm.get(key)
        unmatched += m is None
        d = dev_by_source.get(r["source_row_number"])
        if d is not None:
            require(m is not None and r["gdsc_drug_id"] == d["gdsc_drug_id"] and r["nlme_curve_id"] == d["nlme_curve_id"]
                    and r["min_conc"] == d["min_conc"] and r["max_conc"] == d["max_conc"]
                    and m["canonical_model_id"] == d["canonical_model_id"] and m["expression_eligible"] == "1"
                    and r["response_available"] == "1" and r["webrelease"] == "Y" and r["dataset"] == "GDSC2",
                    "development_record_traces_to_eligible_GDSC_index_row", source_row_number=r["source_row_number"])
            gdsc_selected_found.add(r["source_row_number"])
    inventory = table(R / "all_295_GDSC_drug_sample_counts.csv")
    require(len(inventory) == 295 and len({r["gdsc_drug_id"] for r in inventory}) == 295
            and {r["gdsc_drug_id"] for r in inventory} == gdsc_drug_ids, "GDSC_inventory_covers_295_unique_assay_IDs")
    require(gdsc_rows == 242036 and unmatched == 0, "GDSC_mother_join_preserves_all_rows", before_rows=gdsc_rows, after_rows=gdsc_rows, unmatched_rows=unmatched)
    require(gdsc_selected_found == set(dev_by_source), "all_development_source_row_pointers_found")

    # CTRP validation accesses only identifiers and metadata, never AUC fields.
    refs = defaultdict(list)
    for label, records in (("primary", ext), ("sensitivity", sens)):
        for r in records:
            raw_refs = r["ctrp_source_row_numbers"].split("|")
            require(len(raw_refs) == len(set(raw_refs)) and all(v.isdigit() for v in raw_refs), label + "_source_pointer_list_unique_numeric")
            require(int(r["n_source_experiments"]) == len(as_set(r["ctrp_experiment_ids"])), label + "_experiment_count_matches_ID_list")
            for line in raw_refs:
                refs[line].append((label, r))
    exp_cids, exp_status = defaultdict(set), defaultdict(set)
    for r in iter_csv(ROOT / "data_raw/ctrp/v20.meta.per_experiment.txt", "\t"):
        exp_cids[r["experiment_id"]].add(r["master_ccl_id"])
        exp_status[r["experiment_id"]].add(r["snp_fp_status"])
    require(all(len(v) == 1 for v in exp_cids.values()) and all(len(v) == 1 for v in exp_status.values()), "CTRP_experiment_cell_and_SNP_status_unambiguous")
    official_cells = unique_map(table(ROOT / "data_raw/ctrp/v20.meta.per_cell_line.txt", "\t"), "master_ccl_id", "CTRP_official_cell")
    index_selected = {}
    ctrp_rows, ctrp_unmatched = 0, 0
    raw_selected = {}
    seen_curve_keys = set()
    index_iter = iter_csv(P / "ctrp_curve_index.csv")
    raw_iter = iter_csv(ROOT / "data_raw/ctrp/v20.data.curves_post_qc.txt", "\t")
    for line, raw in enumerate(raw_iter, 2):
        idx = next(index_iter, None)
        require(idx is not None, "CTRP_index_not_shorter_than_source")
        ctrp_rows += 1
        eid, cpd = raw["experiment_id"], raw["master_cpd_id"]
        require(eid in exp_cids, "CTRP_source_experiment_exists")
        cid = next(iter(exp_cids[eid]))
        status = next(iter(exp_status[eid]))
        require((eid, cpd) not in seen_curve_keys, "CTRP_source_experiment_compound_unique")
        seen_curve_keys.add((eid, cpd))
        require(idx["source_row_number"] == str(line) and idx["experiment_id"] == eid
                and idx["master_cpd_id"] == cpd and idx["master_ccl_id"] == cid
                and idx["snp_fp_status"] == status and idx["ccl_name"] == official_cells[cid]["ccl_name"],
                "CTRP_index_identity_matches_physical_source_row")
        ctrp_unmatched += cid not in cm
        if str(line) in refs:
            m = cm[cid]
            index_selected[str(line)] = idx
            raw_selected[str(line)] = (eid, cpd, cid, status)
            for label, rec in refs[str(line)]:
                require(cpd == rec["ctrp_compound_id"] and cpd == drugs[rec["gdsc_drug_id"]]["ctrp_compound_id"]
                        and m["canonical_model_id"] == rec["canonical_model_id"]
                        and cid in as_set(rec["ctrp_master_ccl_ids"]) and eid in as_set(rec["ctrp_experiment_ids"])
                        and status in as_set(rec["identity_qc_statuses"]) and m["expression_eligible"] == "1"
                        and idx["response_available"] == "1" and (label != "primary" or status == "SNP-matched-reference"),
                        label + "_external_source_row_compound_cell_SNP_consistent", source_row_number=str(line))
    require(next(index_iter, None) is None, "CTRP_index_not_longer_than_source")
    require(ctrp_rows == 395263 and ctrp_unmatched == 0, "CTRP_mother_join_preserves_all_rows", before_rows=ctrp_rows, after_rows=ctrp_rows, unmatched_rows=ctrp_unmatched)
    require(set(raw_selected) == set(refs), "all_external_source_row_pointers_found")
    for label, records in (("primary", ext), ("sensitivity", sens)):
        for rec in records:
            source = [raw_selected[x] for x in rec["ctrp_source_row_numbers"].split("|")]
            require({x[0] for x in source} == as_set(rec["ctrp_experiment_ids"])
                    and {x[2] for x in source} == as_set(rec["ctrp_master_ccl_ids"])
                    and {x[3] for x in source} == as_set(rec["identity_qc_statuses"]), label + "_source_metadata_lists_exact")

    # Independently check direct patient and parent links from authority metadata,
    # rather than relying solely on the derived biological-group labels.
    broad = table(ROOT / "data_raw/cell_metadata/DepMap_24Q4_v1_Model.csv")
    patients, sanger_ach = defaultdict(set), defaultdict(set)
    for r in broad:
        if r["PatientID"]:
            patients[r["PatientID"]].add(r["ModelID"])
        if r["SangerModelID"]:
            sanger_ach[r["SangerModelID"]].add(r["ModelID"])
    cross = {r["ccle_id"]: as_set(r["candidate_model_ids"]) for r in table(P / "ccle_accession_crosswalk.csv")}
    for r in table(P / "sanger_identifiers.csv"):
        if r["identifier_source"] in ("BROAD_ID", "DepMap_ID (Broad_ID)") and re.fullmatch(r"ACH-\d+", r["identifier"]):
            sanger_ach[r["sanger_model_id"]].add(r["identifier"])
        elif r["identifier_source"] == "CCLE_ID":
            sanger_ach[r["sanger_model_id"]].update(cross.get(r["identifier"], set()))
    spatient, parents = defaultdict(set), []
    for r in table(P / "sanger_models.csv"):
        if r["sanger_patient_id"]:
            spatient[r["sanger_patient_id"]].update(sanger_ach[r["sanger_model_id"]])
        if r["parent_sanger_model_id"]:
            parents.append((sanger_ach[r["sanger_model_id"]], sanger_ach[r["parent_sanger_model_id"]]))
    for label, records in (("primary", ext), ("sensitivity", sens)):
        xid = {r["canonical_model_id"] for r in records}
        require(not any(v & dids and v & xid for v in patients.values()), label + "_raw_Broad_patient_links_do_not_cross_split")
        require(not any(v & dids and v & xid for v in spatient.values()), label + "_raw_Sanger_patient_links_do_not_cross_split")
        require(not any((a & dids and b & xid) or (a & xid and b & dids) for a, b in parents), label + "_raw_Sanger_parent_links_do_not_cross_split")

    summary = load_json(R / "cohort_summary.json")
    recomputed = {"selected_drugs": len(main), "development_drug_model_observations": len(dev),
                  "external_drug_model_observations": len(ext), "development_unique_models": len(dids),
                  "external_unique_models": len({r["canonical_model_id"] for r in ext}),
                  "development_unique_groups": len(dgroups), "external_unique_groups": len({r["biological_group_id"] for r in ext})}
    require(all(summary[k] == v for k, v in recomputed.items()), "cohort_summary_reconciles_with_independent_counts")
    require(summary["models_fitted"] == 0 and summary["response_magnitudes_used_for_selection"] is False,
            "cohort_stage_declares_no_model_fitting_or_response_selection")
    return {"cohort_counts": recomputed, "per_drug_identity_coverage": per_drug,
            "all_GDSC2_candidate_ACH": len(all_gdsc_ids), "all_GDSC2_candidate_groups": len(all_gdsc_groups),
            "GDSC_mother_rows": gdsc_rows, "CTRP_mother_rows": ctrp_rows,
            "external_source_rows_validated": len(raw_selected),
            "independence_scope": "Known GDSC2 identities and related groups; not a claim about all GDSC1 or unknown genealogy."}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lock", action="store_true", help="Create immutable cohort lock only after every validation passes")
    args = ap.parse_args()
    freeze_hash = digest(FREEZE)
    existing = load_json(LOCK) if LOCK.exists() else None
    # A changed locked artifact must not overwrite the evidence of its earlier validation.
    if existing:
        if existing.get("freeze_sha256") != freeze_hash:
            raise ValidationFailure("Existing lock refers to a different freeze; do not overwrite")
        if not REPORT.is_file() or digest(REPORT) != existing.get("validation_report_sha256"):
            raise ValidationFailure("Existing locked validation report is absent or changed; do not overwrite")
    freeze = load_json(FREEZE)
    result = {"schema_version": "1.0", "freeze_sha256": freeze_hash,
              "freeze_created_utc": freeze["created_utc"],
              "response_values_read_for_analysis": False, "response_values_reported": False,
              "validation_scope": "Identity, lineage, source traceability, hashes, cohort coverage and isolation; no statistical performance validation."}
    error = None
    try:
        result.update(validate(freeze))
        require(digest(FREEZE) == freeze_hash, "freeze_unchanged_during_validation")
        result["status"] = "passed"
    except Exception as exc:
        error = exc
        result["status"] = "failed"
        result["failure"] = str(exc)
    # Compact repeated row-level checks; preserve every distinct check and its counts.
    grouped = {}
    for c in CHECKS:
        item = grouped.setdefault(c["check"], {"check": c["check"], "passed": True, "n_checks": 0})
        item["n_checks"] += 1
        item["passed"] = item["passed"] and c["passed"]
        if not c["passed"]:
            item["failure_evidence"] = {k: v for k, v in c.items() if k not in ("check", "passed")}
    result["checks"] = list(grouped.values())
    result["n_assertions"] = len(CHECKS)
    payload = json_bytes(result)
    report_hash = hashlib.sha256(payload).hexdigest()
    if existing:
        if result["status"] != "passed" or report_hash != existing["validation_report_sha256"]:
            raise ValidationFailure("Current validation disagrees with immutable lock/report; no files overwritten")
    else:
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_bytes(payload)
    if error is not None:
        raise ValidationFailure("Validation failed: " + str(error))
    if args.lock and not existing:
        lock = {"schema_version": "1.0", "freeze_sha256": freeze_hash,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "validation_report_sha256": report_hash,
                "validation_report": "reports/validation_summary.json", "status": "validated_cohort_lock",
                "training_configuration_frozen": False}
        with LOCK.open("xb") as f:
            f.write(json_bytes(lock))
    print(json.dumps({"status": "passed", "report": str(REPORT), "lock_exists": LOCK.exists(),
                      "n_assertions": len(CHECKS), **result["cohort_counts"]}, ensure_ascii=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("VALIDATION FAILED: " + str(exc), file=sys.stderr)
        raise SystemExit(1)
