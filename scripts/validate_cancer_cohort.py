"""Independent metadata-only cohort-v2 validation and immutable lock creation.

Does not load drug-response tables/indices or expression-array bytes. Large
scientific inputs inherited from v1 retain their recorded v1 hash descriptors;
the development runner must verify those bytes before their eventual use.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.parse import urlparse

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "results/cohort_v2"
REPORT = ROOT / "reports/cancer_cohort_v2_validation.json"
LOCK = ROOT / "config/cohort_freeze_v2.lock.json"


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_csv(path):
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    forbidden = {"auc", "response_auc", "observed_auc", "predicted_auc", "area_under_curve", "ln_ic50", "z_score"}
    require(not forbidden & {c.lower() for c in frame.columns}, "Response values found in a metadata-only manifest")
    return frame


def truth(value):
    require(str(value).lower() in {"true", "false", "1", "0"}, "Invalid metadata Boolean")
    return str(value).lower() in {"true", "1"}


def relative_path(name):
    path = (ROOT / str(name).replace("\\", "/")).resolve()
    require(path.is_relative_to(ROOT), "Scientific dependency escapes the repository")
    return path


def install_metadata_guard():
    """Reject accidental response/index/array access by this validator."""
    forbidden_indices = {"ctrp_curve_index.csv", "ctrp_pair_availability.csv", "ctrp_drug_availability.csv"}
    blocked = []
    def hook(event, args):
        if event != "open" or not args or not isinstance(args[0], (str, bytes, Path)):
            return
        try:
            path = Path(args[0].decode() if isinstance(args[0], bytes) else args[0]).resolve()
        except (OSError, ValueError):
            return
        if (path.is_relative_to(ROOT / "data_raw/ctrp") or path.name in forbidden_indices
                or path.name.casefold() == "ctrp--drug--cellline--auc.txt"
                or path.suffix.casefold() in {".npy", ".npz", ".gct"}):
            blocked.append(path.name)
            raise PermissionError("Metadata-only validation cannot open response/index/array inputs")
    sys.addaudithook(hook)
    return blocked


def validate_override(freeze):
    evidence = read_json(ROOT / "data_raw/cancer_metadata/evidence_manifest.json")
    for source in evidence["sources"]:
        path = relative_path(source["path"])
        require(path.is_relative_to(ROOT / "data_raw/cancer_metadata"), "Override source outside evidence cache")
        require(sha(path) == source["sha256"], "Cached override evidence hash mismatch")
    api = urlparse(evidence["cellosaurus_api_url"])
    require(api.scheme == "https" and api.hostname == "api.cellosaurus.org" and api.path == "/cell-line/CVCL_A351", "Override lacks official exact-accession API provenance")
    record = read_json(ROOT / "data_raw/cancer_metadata/cellosaurus_CVCL_A351.json")["Cellosaurus"]["cell-line-list"]
    require(len(record) == 1, "Override cache has more than one cell-line record")
    record = record[0]
    require({a["value"] for a in record["accession-list"] if a["type"] == "primary"} == {"CVCL_A351"}, "Cellosaurus primary accession mismatch")
    require(record["category"] == "Cancer cell line", "Cellosaurus does not classify the override as a cancer cell line")
    require(any(d["database"] == "NCIt" and d["accession"] == "C3493" and d["label"] == "Lung squamous cell carcinoma" for d in record["disease-list"]), "Authoritative cancer diagnosis missing")
    crossrefs = {(r["database"], r["accession"]) for r in record["xref-list"]}
    require({("DepMap", "ACH-001078"), ("Cell_Model_Passport", "SIDM01820"), ("KCLB", "71588")} <= crossrefs, "Cellosaurus exact-ID identity chain is incomplete")
    sanger = read_csv(ROOT / "data_processed/sanger_identifiers.csv")
    official = sanger[sanger.sanger_model_id.eq("SIDM01820")]
    official_pairs = set(zip(official.identifier_source, official.identifier))
    required_pairs = {("BROAD_ID", "ACH-001078"), ("CCLE_ID", "HCC1588_LUNG"), ("Cellosaurus (RRID)", "CVCL_A351"), ("KCLB", "71588")}
    require(required_pairs <= official_pairs, "Frozen Sanger exact crossrefs do not corroborate the override")
    cached_sanger = read_json(ROOT / "data_raw/cancer_metadata/sanger_SIDM01820_identity_evidence.json")
    cached_pairs = {(r["identifier_source"], r["identifier"]) for r in cached_sanger["identifiers"]}
    require(cached_sanger["sanger_model_id"] == "SIDM01820" and cached_pairs == official_pairs, "Derived Sanger cache differs from frozen authoritative identity metadata")
    provider = read_json(ROOT / "data_raw/cancer_metadata/kclb_71588_eligibility_facts.json")
    require(urlparse(provider["source_url"]).hostname == "cellbank.snu.ac.kr" and provider["kclb_catalog_id"] == "71588", "Provider source is not the exact authoritative catalog record")
    require(provider["direct_http_snapshot_status"] == "retrieved" and re.fullmatch(r"[0-9a-f]{64}", provider["source_body_sha256"]), "Provider cache lacks successful retrieval/checksum provenance")
    require(provider["origin"].casefold() == "lung" and provider["histopathology"] == "squamous cell carcinoma", "Provider eligibility facts conflict with Cellosaurus")
    license_record = read_json(ROOT / "data_raw/cancer_metadata/cellosaurus_license.json")
    require(urlparse(license_record["source_url"]).hostname == "www.cellosaurus.org" and license_record["license_url"] == "https://creativecommons.org/licenses/by/4.0/", "Cellosaurus source/license provenance missing")
    expected = freeze["rules"]["authoritative_override"]
    require(expected["canonical_model_id"] == "ACH-001078" and expected["cellosaurus_id"] == "CVCL_A351" and expected["sanger_model_id"] == "SIDM01820", "Override rule identity mismatch")
    return {"canonical_model_id": "ACH-001078", "ccle_id": "HCC1588_LUNG", "cellosaurus_id": "CVCL_A351", "sanger_model_id": "SIDM01820", "kclb_id": "71588", "category": "Cancer cell line", "disease": "Lung squamous cell carcinoma", "exact_identity_chain_verified": True}


def validate():
    freeze_path = ROOT / "config/cohort_freeze_v2.json"
    freeze = read_json(freeze_path)
    parent_path = ROOT / "config/cohort_freeze_v1.json"
    parent_lock_path = ROOT / "config/cohort_freeze_v1.lock.json"
    parent = read_json(parent_path)
    parent_lock = read_json(parent_lock_path)
    require(sha(parent_path) == freeze["parent_cohort_freeze_sha256"] == parent_lock["freeze_sha256"], "Historical v1 freeze SHA changed")
    require(sha(parent_lock_path) == freeze["parent_lock_sha256"], "Historical v1 lock SHA changed")
    require(freeze["cohort_version"] == "v2" and freeze["external_response_accessed"] is False and freeze["v1_files_overwritten"] is False, "Unexpected cohort version/access/overwrite declaration")
    rules = read_json(ROOT / "config/cancer_cohort_v2_rules.json")
    require(rules == freeze["rules"] and rules["selection_uses_response_values"] is False and rules["selection_uses_model_performance"] is False, "Eligibility rules differ from freeze or use outcomes")
    require(rules["minimum_development_groups"] >= 200 and rules["minimum_external_primary_groups"] >= 80, "Eligibility weakened declared minimum sample counts")
    inherited = {str(r["path"]).replace("\\", "/"): r for r in parent["files"]}
    dependency_checks = []
    for source in freeze["files"]:
        name = str(source["path"]).replace("\\", "/")
        path = relative_path(name)
        require(not name.startswith(("reports/", "runs/")) and "RESEARCH_PLAN" not in name and name not in {"config/cohort_freeze_v2.json", "config/cohort_freeze_v2.lock.json"}, "Mutable or circular dependency included in scientific freeze")
        if source["scope"] == "inherited_v1_scientific_input":
            require(name in inherited and source["sha256"] == inherited[name]["sha256"] and source["bytes"] == inherited[name]["bytes"], "Inherited scientific input descriptor differs from locked v1")
            dependency_checks.append({"path": name, "verification": "descriptor_identical_to_locked_v1; bytes deferred to development runner"})
            continue
        require(source["scope"] == "v2_direct_dependency", "Unknown dependency scope")
        require(path.suffix.lower() in {".json", ".csv", ".py"} and not name.startswith("data_raw/ctrp/") and path.name not in {"ctrp_curve_index.csv", "ctrp_pair_availability.csv", "ctrp_drug_availability.csv"}, "Non-metadata dependency would breach validation scope")
        require(path.stat().st_size == source["bytes"] and sha(path) == source["sha256"], "Direct scientific dependency changed: " + name)
        if name in inherited:
            require(source["sha256"] == inherited[name]["sha256"], "Historical v1 scientific metadata overwritten: " + name)
        dependency_checks.append({"path": name, "verification": "actual_metadata_or_code_bytes_match_freeze"})
    override = validate_override(freeze)

    old_global = read_csv(ROOT / "results/global_split_manifest.csv").set_index("ccle_id", drop=False)
    new_global = read_csv(V2 / "global_split_manifest.csv").set_index("ccle_id", drop=False)
    decisions = read_csv(V2 / "eligibility_decisions.csv").set_index("ccle_id", drop=False)
    require(old_global.index.is_unique and new_global.index.is_unique and decisions.index.is_unique and set(old_global.index) == set(new_global.index) == set(decisions.index), "Global identity/decision rows were added, removed or duplicated")
    broad = read_csv(ROOT / "data_raw/cell_metadata/DepMap_24Q4_v1_Model.csv").set_index("ModelID")
    expected_eligible = {}
    statuses = {}
    for ccle, original in old_global.iterrows():
        decision, current = decisions.loc[ccle], new_global.loc[ccle]
        for column in ("expression_row_index", "canonical_model_id", "biological_group_id"):
            require(decision[column] == original[column] == current[column], "Canonical identity or expression-row assignment changed")
        require(decision.v1_tissue == current.v1_tissue == original.tissue and decision.tissue == current.tissue, "Original tissue value was not preserved")
        require(decision.tissue.casefold() == original.tissue.casefold(), "Tissue change exceeds capitalization-only rule")
        require(current.v1_global_role == original.global_role, "Original global role provenance lost")
        identifier = original.canonical_model_id
        disease = broad.loc[identifier, "OncotreePrimaryDisease"] if identifier in broad.index else ""
        if not identifier or not truth(original.expression_eligible):
            eligible, status = False, "excluded_unresolved_identity"
        elif disease == "Non-Cancerous":
            eligible, status = False, "excluded_non_cancerous"
        elif disease.strip() and disease.casefold() != "unknown":
            eligible, status = True, "included_broad_disease_annotation"
            require(decision.disease_annotation == disease, "Broad disease annotation changed")
        elif identifier == override["canonical_model_id"]:
            require(ccle == override["ccle_id"], "Override canonical ID linked to the wrong CCLE accession")
            eligible, status = True, "included_authoritative_cancer_override"
        else:
            eligible, status = False, "excluded_unverified_disease"
        require(truth(decision.included_v2) == truth(current.cancer_eligible_v2) == eligible and decision.eligibility_status == status, "Decision does not follow independent authoritative-metadata classification")
        require(decision.exclusion_reason == "" if eligible else decision.exclusion_reason != "", "Eligibility/exclusion reason missing or contradictory")
        if not eligible:
            require(not truth(current.in_selected_development) and not truth(current.in_selected_external_primary), "Excluded identity remains in a selected cohort flag")
            require(current.global_role.startswith("excluded"), "Excluded identity retains a modelling role")
        else:
            require(current.global_role == original.global_role, "Eligible global role was reassigned")
        expected_eligible[ccle] = eligible
        statuses[status] = statuses.get(status, 0) + 1

    manifests = {}
    removed_counts = {}
    excluded = read_csv(V2 / "excluded_sample_manifest.csv")
    for kind in ("development", "external_primary", "external_sensitivity"):
        filename = kind + "_sample_manifest.csv"
        original = read_csv(ROOT / "results" / filename)
        current = read_csv(V2 / filename)
        keys = ["gdsc_drug_id", "canonical_model_id"]
        require(not original.duplicated(keys).any() and not current.duplicated(keys).any(), "Duplicated drug/model manifest rows")
        expected = original[original.ccle_id.map(expected_eligible)].set_index(keys).sort_index()
        retained = current.set_index(keys).sort_index()
        require(list(expected.index) == list(retained.index), "Cohort change exceeds metadata-only eligibility filtering")
        for column in original.columns:
            if column in keys:
                continue
            if column == "tissue":
                require(retained.v1_tissue.equals(expected.tissue) and retained.tissue.map(str.casefold).equals(expected.tissue.map(str.casefold)), "Manifest tissue change is not capitalization only")
            else:
                require(retained[column].equals(expected[column]), "Drug, dose, source record or identity changed: " + column)
        expected_role = "external_identity_qc_sensitivity" if kind == "external_sensitivity" else kind
        require(set(current.role) == {expected_role} and current.ccle_id.map(expected_eligible).all(), "Incorrect retained role or eligibility")
        removed = original[~original.ccle_id.map(expected_eligible)]
        recorded = excluded[excluded.source_manifest.eq(kind)]
        require(set(map(tuple, removed[keys].to_numpy())) == set(map(tuple, recorded[keys].to_numpy())) and len(removed) == len(recorded), "Excluded-source ledger incomplete or duplicated")
        manifests[kind] = current
        removed_counts[kind] = {"unique_models": removed.canonical_model_id.nunique(), "drug_model_observations": len(removed)}

    dev, ext, sensitivity = (manifests[k] for k in ("development", "external_primary", "external_sensitivity"))
    require(set(ext.identity_qc_statuses) == {"SNP-matched-reference"}, "Primary external identity QC requirement changed")
    require(not set(dev.canonical_model_id) & set(ext.canonical_model_id) and not set(dev.biological_group_id) & set(ext.biological_group_id), "Development/external identity or group overlap")
    gdsc = read_csv(ROOT / "data_processed/gdsc2_cell_mapping.csv")
    model_groups = read_csv(ROOT / "data_processed/canonical_model_groups.csv")
    groups_for = dict(zip(model_groups.canonical_model_id, model_groups.biological_group_id))
    all_gdsc_candidates = {identifier for value in gdsc.candidate_model_ids for identifier in value.split("|") if identifier}
    require(all_gdsc_candidates <= set(groups_for), "Unrepresented candidate group prevents complete GDSC exclusion verification")
    all_gdsc_groups = {groups_for[i] for i in all_gdsc_candidates}
    require(not set(ext.biological_group_id) & all_gdsc_groups and not set(sensitivity.biological_group_id) & all_gdsc_groups, "External identity-QC pool overlaps any known GDSC candidate group")
    selected_development = set(dev.ccle_id)
    selected_external = set(ext.ccle_id)
    require({ccle for ccle, row in new_global.iterrows() if truth(row.in_selected_development)} == selected_development and {ccle for ccle, row in new_global.iterrows() if truth(row.in_selected_external_primary)} == selected_external, "Global selected flags differ from actual manifests")

    old_panel = read_csv(ROOT / "results/main_drugs.csv")
    panel = read_csv(V2 / "main_drugs.csv")
    require(len(panel) == len(old_panel) == 8 and panel.gdsc_drug_id.tolist() == old_panel.gdsc_drug_id.tolist(), "Frozen panel/rank changed")
    mutable_counts = {c for c in old_panel.columns if c.startswith("n_")}
    for column in old_panel.columns:
        if column not in mutable_counts:
            require(panel[column].equals(old_panel[column]), "Frozen drug identity/dose/structure/rank changed: " + column)
    per_drug = []
    for row in panel.itertuples(index=False):
        development = dev[dev.gdsc_drug_id.eq(row.gdsc_drug_id)]
        external = ext[ext.gdsc_drug_id.eq(row.gdsc_drug_id)]
        extra = sensitivity[sensitivity.gdsc_drug_id.eq(row.gdsc_drug_id)]
        ndev, next_ = development.biological_group_id.nunique(), external.biological_group_id.nunique()
        require(ndev >= rules["minimum_development_groups"] and next_ >= rules["minimum_external_primary_groups"] and development.tissue.nunique() >= 4, "A frozen drug fails sample-size/tissue gates")
        counts = {"n_development_models": development.canonical_model_id.nunique(), "n_development_groups": ndev,
                  "n_development_tissues": development.tissue.nunique(), "n_external_primary_models": external.canonical_model_id.nunique(),
                  "n_external_primary_groups": next_, "n_external_primary_tissues": external.tissue.nunique(),
                  "n_external_primary_curve_records": int(external.n_source_experiments.astype(int).sum()),
                  "n_external_all_identity_status_models": extra.canonical_model_id.nunique()}
        for name, count in counts.items():
            require(int(getattr(row, name)) == count, "Panel count does not match retained manifests: " + name)
        require(set(development.min_conc) == {row.gdsc_min_conc} and set(development.max_conc) == {row.gdsc_max_conc}, "Retained development dose window changed")
        require(set(external.ctrp_compound_id) == {row.ctrp_compound_id}, "Retained external compound ID changed")
        per_drug.append({"gdsc_drug_id": row.gdsc_drug_id, "drug_name": row.gdsc_drug_name, "development_groups": ndev, "external_primary_groups": next_})
    summary = read_json(V2 / "cohort_summary.json")
    summary_counts = {"development_unique_models": dev.canonical_model_id.nunique(), "development_unique_groups": dev.biological_group_id.nunique(),
                      "external_unique_models": ext.canonical_model_id.nunique(), "external_unique_groups": ext.biological_group_id.nunique(),
                      "development_drug_model_observations": len(dev), "external_drug_model_observations": len(ext),
                      "external_sensitivity_unique_models": sensitivity.canonical_model_id.nunique(), "external_sensitivity_drug_model_observations": len(sensitivity)}
    for name, count in summary_counts.items():
        require(summary[name] == count, "Aggregate cohort summary differs from independent counts: " + name)
    require(summary["response_values_read"] is False and summary["response_indices_read"] is False and summary["expression_distributions_read"] is False, "Metadata-only scope contradicted by cohort summary")
    return {"status": "passed", "schema_version": "1.0", "validated_utc": datetime.now(timezone.utc).isoformat(),
            "cohort_version": "v2", "freeze_sha256": sha(freeze_path), "parent_cohort_freeze_sha256": sha(parent_path),
            "parent_lock_sha256": sha(parent_lock_path), "original_response_or_index_or_expression_array_read": False,
            "metadata_only_qualification_verified": True, "v1_freeze_and_used_scientific_metadata_unchanged": True,
            "inherited_array_byte_validation_deferred_to_development_runner": True,
            "drug_panel_identity_rank_structures_and_dose_windows_unchanged": True,
            "all_known_GDSC_candidate_group_exclusion_verified": True,
            "authoritative_override": override, "global_eligibility_status_counts": statuses,
            "cohort_counts": summary_counts, "removed": removed_counts, "per_drug_counts": per_drug,
            "scientific_dependency_checks": dependency_checks,
            "limitations": ["Eligibility follows the exact frozen Broad disease annotations and one exact-ID authoritative override; no new histology inference is made.",
                            "Large expression input bytes are deliberately not reopened during metadata-only validation; inherited SHA descriptors are matched to v1 and will be verified on development access.",
                            "Known group separation cannot rule out undocumented relatedness; CTRP was historically used and this is a retrospective redesign."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", action="store_true", help="Create the immutable v2 lock only after successful independent validation")
    args = parser.parse_args()
    blocked = install_metadata_guard()
    try:
        report = validate()
        require(not blocked, "A forbidden input access was attempted")
        if LOCK.exists():
            lock = read_json(LOCK)
            require(lock["freeze_sha256"] == report["freeze_sha256"] and REPORT.exists() and sha(REPORT) == lock["validation_report_sha256"], "Existing v2 lock/report was changed")
            previous = read_json(REPORT)
            current = {k: v for k, v in report.items() if k != "validated_utc"}
            previous = {k: v for k, v in previous.items() if k != "validated_utc"}
            require(current == previous, "Revalidation differs from the locked validation report")
            # Preserve the exact locked validation bytes on successful rechecks.
            print(json.dumps({"status": "passed", "lock_status": "existing_immutable_lock_verified", "cohort_counts": report["cohort_counts"]}, ensure_ascii=False))
            return 0
        REPORT.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if args.lock:
            lock = {"status": "validated_cancer_cohort_lock", "schema_version": "2.0", "cohort_version": "v2",
                    "created_utc": datetime.now(timezone.utc).isoformat(), "freeze_sha256": report["freeze_sha256"],
                    "parent_cohort_freeze_sha256": report["parent_cohort_freeze_sha256"], "parent_lock_sha256": report["parent_lock_sha256"],
                    "validation_report": "reports/cancer_cohort_v2_validation.json", "validation_report_sha256": sha(REPORT),
                    "training_configuration_frozen": False, "external_response_accessed": False}
            with LOCK.open("x", encoding="utf-8") as stream:
                json.dump(lock, stream, ensure_ascii=False, indent=2, allow_nan=False)
                stream.write("\n")
        print(json.dumps({"status": "passed", "lock_status": "created" if args.lock else "not_requested", "cohort_counts": report["cohort_counts"]}, ensure_ascii=False))
        return 0
    except Exception as exc:
        failure = {"status": "failed", "validated_utc": datetime.now(timezone.utc).isoformat(), "error_type": type(exc).__name__,
                   "error": str(exc).replace(str(ROOT), "<repository>"), "blocked_input_accesses": blocked}
        target = REPORT if not LOCK.exists() else ROOT / "reports/cancer_cohort_v2_revalidation_failure.json"
        target.write_text(json.dumps(failure, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(failure, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    sys.exit(main())
