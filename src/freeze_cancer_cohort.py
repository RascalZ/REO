"""Create a metadata-only cancer eligibility overlay of locked cohort v1.

No response values, response indices, or expression matrices are opened. The
eight-drug panel, rank, assay IDs, dose windows and source-row identities remain
fixed. v1 and its historical pilot are preserved. Independent validation must
create the v2 lock; this script creates only a pending snapshot.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import re
import urllib.error
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/cohort_v2"
CACHE = ROOT / "data_raw/cancer_metadata"
RULE_PATH = ROOT / "config/cancer_cohort_v2_rules.json"
FREEZE_PATH = ROOT / "config/cohort_freeze_v2.json"
LOCK_PATH = ROOT / "config/cohort_freeze_v2.lock.json"
CELL_API = "https://api.cellosaurus.org/cell-line/CVCL_A351?fields=id,ac,ca,di,dr&format=json"
CELL_PAGE = "https://www.cellosaurus.org/CVCL_A351"
CELL_LICENSE = "https://www.cellosaurus.org/description.html"
KCLB_PAGE = "https://cellbank.snu.ac.kr/cellline/search?q=71588&sc=y"
SANGER_LICENSE = "https://depmap.sanger.ac.uk/documentation/data-usage-policy/"
PARENT_INPUTS = [
    "results/main_drugs.csv", "results/development_sample_manifest.csv",
    "results/external_primary_sample_manifest.csv",
    "results/external_sensitivity_sample_manifest.csv",
    "results/global_split_manifest.csv", "data_processed/cell_line_dictionary.csv",
    "data_processed/canonical_model_groups.csv", "data_processed/gdsc2_cell_mapping.csv",
    "data_raw/cell_metadata/DepMap_24Q4_v1_Model.csv",
    "data_processed/sanger_models.csv", "data_processed/sanger_identifiers.csv",
]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read(path):
    return pd.read_csv(ROOT / path, dtype=str, keep_default_na=False)


def csv(frame, name):
    frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")


def download(url):
    request = urllib.request.Request(url, headers={"User-Agent": "REO-research-metadata-curation/2.0"})
    with urllib.request.urlopen(request, timeout=45) as response:
        return response.read(), response.geturl()


def fetch_evidence():
    """Explicit optional cache step; never fetch changing sources during freeze."""
    if (CACHE / "evidence_manifest.json").exists():
        raise FileExistsError("Preserve cached evidence. Use a new evidence/cohort version to refresh it.")
    CACHE.mkdir(parents=True, exist_ok=True)
    retrieved = datetime.now(timezone.utc).isoformat()
    if (CACHE / "cellosaurus_CVCL_A351.json").exists():
        raw = (CACHE / "cellosaurus_CVCL_A351.json").read_bytes()
        final_url = CELL_API
    else:
        raw, final_url = download(CELL_API)
    cell = json.loads(raw)["Cellosaurus"]["cell-line-list"][0]
    assert cell["category"] == "Cancer cell line"
    assert {x["value"] for x in cell["accession-list"]} == {"CVCL_A351"}
    assert {x["value"] for x in cell["name-list"]} == {"HCC1588"}
    cross = {(x["database"], x["accession"]) for x in cell["xref-list"]}
    assert {("DepMap", "ACH-001078"), ("Cell_Model_Passport", "SIDM01820"), ("KCLB", "71588")} <= cross
    assert any(x["accession"] == "C3493" and x["label"] == "Lung squamous cell carcinoma" for x in cell["disease-list"])
    if not (CACHE / "cellosaurus_CVCL_A351.json").exists():
        (CACHE / "cellosaurus_CVCL_A351.json").write_bytes(raw)
    if not (CACHE / "cellosaurus_license.json").exists():
        license_raw, license_url = download(CELL_LICENSE)
        license_text = " ".join(html.unescape(re.sub(r"<[^>]+>", " ", license_raw.decode("utf-8"))).split())
        assert "CC BY 4.0" in license_text
        dump(CACHE / "cellosaurus_license.json", {
            "source_url": CELL_LICENSE, "resolved_url": license_url, "retrieved_utc": retrieved,
            "source_body_sha256": hashlib.sha256(license_raw).hexdigest(),
            "license": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
            "license_url": "https://creativecommons.org/licenses/by/4.0/",
            "cache_kind": "license fact derived from official description; not full HTML",
        })
    try:
        kclb_raw, kclb_url = download(KCLB_PAGE)
    except (urllib.error.URLError, TimeoutError) as exc:
        # The exact-ID Cellosaurus record plus Sanger crosswalk is sufficient.
        # Do not fabricate a source-body hash or independent bank snapshot.
        kclb_facts = {"source_url": KCLB_PAGE, "retrieved_utc": retrieved,
            "direct_http_snapshot_status": "unavailable", "transport_error_type": type(exc).__name__,
            "not_required_for_override": True,
            "kclb_catalog_id": "71588", "catalog_id_evidence": "Cellosaurus CVCL_A351 exact KCLB cross-reference",
            "cache_kind": "source availability record only; no bank data snapshot asserted"}
    else:
        kclb_text = " ".join(html.unescape(re.sub(r"<[^>]+>", " ", kclb_raw.decode("utf-8", errors="replace"))).split())
        assert all(term in kclb_text for term in ["HCC-1588", "71588", "lung", "squamous cell carcinoma"])
        kclb_facts = {"source_url": KCLB_PAGE, "resolved_url": kclb_url, "retrieved_utc": retrieved,
            "source_body_sha256": hashlib.sha256(kclb_raw).hexdigest(), "direct_http_snapshot_status": "retrieved",
            "cell_line_name": "HCC-1588", "kclb_catalog_id": "71588", "origin": "lung", "histopathology": "squamous cell carcinoma",
            "cache_kind": "minimal source facts, not redistribution of full copyrighted page",
            "rights_note": "KCLB page states all rights reserved; no open-data license inferred."}
    dump(CACHE / "kclb_71588_eligibility_facts.json", kclb_facts)
    identifiers = read("data_processed/sanger_identifiers.csv")
    identifiers = identifiers[identifiers.sanger_model_id.eq("SIDM01820")]
    models = read("data_processed/sanger_models.csv")
    models = models[models.sanger_model_id.eq("SIDM01820")]
    assert len(models) == 1
    values = set(zip(identifiers.identifier_source, identifiers.identifier))
    assert {("BROAD_ID", "ACH-001078"), ("CCLE_ID", "HCC1588_LUNG"), ("Cellosaurus (RRID)", "CVCL_A351"), ("KCLB", "71588")} <= values
    dump(CACHE / "sanger_SIDM01820_identity_evidence.json", {
        "source_url": "https://cellmodelpassports.sanger.ac.uk/passports/SIDM01820",
        "source_snapshot": "data_raw/sanger_metadata/cell_model_passports_2026-10-10.json",
        "derived_inputs": ["data_processed/sanger_models.csv", "data_processed/sanger_identifiers.csv"],
        "sanger_model_id": "SIDM01820", "names": models.iloc[0]["names"],
        "tissue": models.iloc[0]["sanger_tissue"],
        "identifiers": identifiers[["identifier_source", "identifier"]].to_dict("records"),
        "license_source": SANGER_LICENSE,
        "license_note": "Official API permits individual non-commercial research use; no unrestricted redistribution license inferred.",
    })
    files = ["cellosaurus_CVCL_A351.json", "cellosaurus_license.json", "kclb_71588_eligibility_facts.json", "sanger_SIDM01820_identity_evidence.json"]
    dump(CACHE / "evidence_manifest.json", {
        "schema_version": "1.0", "retrieved_utc": retrieved,
        "scope": "Malignancy eligibility and exact identity only; no response values or expression distributions.",
        "cellosaurus_api_url": CELL_API, "cellosaurus_page_url": CELL_PAGE,
        "cellosaurus_resolved_api_url": final_url,
        "sources": [{"path": (CACHE / name).relative_to(ROOT).as_posix(), "sha256": sha(CACHE / name)} for name in files],
        "resolved_identity": {"canonical_model_id": "ACH-001078", "ccle_id": "HCC1588_LUNG", "cellosaurus_id": "CVCL_A351", "sanger_model_id": "SIDM01820", "kclb_catalog_id": "71588"},
        "eligibility_resolution": "Cancer cell line; lung squamous cell carcinoma, directly linked by exact DepMap/Sanger/KCLB accessions.",
        "license_attribution": "Cellosaurus, SIB Swiss Institute of Bioinformatics, CC BY 4.0; official API eligibility record retained without modification.",
    })


def verify_parent():
    parent_path = ROOT / "config/cohort_freeze_v1.json"
    lock_path = ROOT / "config/cohort_freeze_v1.lock.json"
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    assert sha(parent_path) == lock["freeze_sha256"], "v1 freeze/lock mismatch"
    recorded = {Path(row["path"]).as_posix(): row for row in parent["files"] if row.get("scope") != "original_read_only"}
    for name in PARENT_INPUTS:
        assert name in recorded, f"Required parent scientific input is not frozen: {name}"
        assert sha(ROOT / name) == recorded[name]["sha256"], f"Parent input changed: {name}"
    return parent, lock, recorded


def build_rules():
    return {
        "schema_version": "2.0", "cohort_version": "v2", "parent_cohort": "v1",
        "purpose": "Metadata-only cancer-model eligibility overlay for formal analysis; preserve historical cohort v1 and pilot.",
        "selection_uses_response_values": False, "selection_uses_model_performance": False,
        "timing_disclosure": "Amendment identified after development-only pilot, before CTRP response access; not a claim of pre-pilot preregistration.",
        "explicit_non_cancerous": "Exclude exact Broad OncotreePrimaryDisease=Non-Cancerous from cancer-focused analysis; retain in excluded-control manifest.",
        "unknown_disease": "Quarantine absent/empty/Unknown primary-disease annotation unless an authoritative exact-ID cancer classification is verified and cached.",
        "other_broad_disease_annotation": "Retain nonempty disease annotations other than explicit Non-Cancerous/Unknown; preserve exact source term without inferring further histology.",
        "authoritative_override": {"canonical_model_id": "ACH-001078", "ccle_id": "HCC1588_LUNG", "category": "Cancer cell line", "disease": "Lung squamous cell carcinoma", "disease_id": "NCIt:C3493", "cellosaurus_id": "CVCL_A351", "sanger_model_id": "SIDM01820", "kclb_id": "71588", "evidence_manifest": "data_raw/cancer_metadata/evidence_manifest.json"},
        "identity": "Unresolved/ambiguous canonical identity remains excluded; eligibility does not repair identity conflicts.",
        "related_groups": "Preserve v1 biological-group IDs and exclusion against all known GDSC groups, including conflicting identity candidates.",
        "panel": "Keep exactly the eight v1 drugs in the original frozen rank, with unchanged GDSC/CTRP IDs, structures, and exact dose windows. No reranking/replacement.",
        "tissue_rule": "Exact casefold equivalence to authoritative Broad OncotreeLineage vocabulary, then Sanger tissue vocabulary; only capitalization canonicalization, no fuzzy/name/histology merges. Preserve v1_tissue on every copied sample/global manifest.",
        "primary_external_qc": "Preserve v1 SNP-matched-reference requirement and all source records/aggregation; no external response or expression distribution read.",
        "sensitivity": "Cancer-eligible v1 external identity-QC sensitivity list retained separately; explicit non-cancerous/unknown excluded rather than pooled as cancer sensitivity.",
        "minimum_development_groups": 200, "minimum_external_primary_groups": 80,
        "formal_model_config": "Not frozen by this cohort snapshot; independent model/statistical/resource gates remain.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch-metadata", action="store_true")
    args = parser.parse_args()
    if args.fetch_metadata:
        fetch_evidence()
        print("Authoritative metadata/rights evidence cached. No cohort files generated.")
        return
    if LOCK_PATH.exists():
        raise FileExistsError("cohort v2 is locked. Preserve it; create a new version for amendments.")
    if FREEZE_PATH.exists() or OUT.exists() or RULE_PATH.exists():
        raise FileExistsError("Preserve existing v2 snapshot. Use a new version rather than silently overwrite.")
    parent, parent_lock, recorded = verify_parent()
    evidence = json.loads((CACHE / "evidence_manifest.json").read_text(encoding="utf-8"))
    for row in evidence["sources"]:
        assert sha(ROOT / row["path"]) == row["sha256"], "Cached metadata evidence changed"
    cell = json.loads((CACHE / "cellosaurus_CVCL_A351.json").read_text(encoding="utf-8"))["Cellosaurus"]["cell-line-list"][0]
    assert cell["category"] == "Cancer cell line"
    assert ("DepMap", "ACH-001078") in {(x["database"], x["accession"]) for x in cell["xref-list"]}
    assert any(x["accession"] == "C3493" for x in cell["disease-list"])
    rules = build_rules()
    model = read("data_raw/cell_metadata/DepMap_24Q4_v1_Model.csv")
    assert not model.ModelID.duplicated().any()
    dictionary = read("results/global_split_manifest.csv")
    annotation_columns = ["ModelID", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeLineage"]
    decisions = dictionary.merge(model[annotation_columns], left_on="canonical_model_id", right_on="ModelID", how="left", validate="many_to_one").fillna("")
    decisions["v1_tissue"] = decisions.tissue
    decisions["eligibility_status"] = "included_broad_disease_annotation"
    decisions["exclusion_reason"] = ""
    decisions["evidence_source"] = "Broad24Q4_Model_OncotreePrimaryDisease"
    decisions["disease_annotation"] = decisions.OncotreePrimaryDisease
    unknown = decisions.OncotreePrimaryDisease.str.strip().str.casefold().isin(["", "unknown", "not reported"])
    decisions.loc[unknown, ["eligibility_status", "exclusion_reason"]] = ["excluded_unverified_disease", "no_verified_cancer_classification"]
    noncancer = decisions.OncotreePrimaryDisease.eq("Non-Cancerous")
    decisions.loc[noncancer, ["eligibility_status", "exclusion_reason"]] = ["excluded_non_cancerous", "explicit_Broad_Non-Cancerous"]
    override = decisions.canonical_model_id.eq("ACH-001078")
    assert override.sum() == 1 and decisions.loc[override, "ccle_id"].iloc[0] == "HCC1588_LUNG"
    decisions.loc[override, ["eligibility_status", "exclusion_reason", "evidence_source", "disease_annotation"]] = ["included_authoritative_cancer_override", "", "Cellosaurus_CVCL_A351_exact_DepMap_and_Sanger_crossrefs_plus_KCLB71588", "Lung squamous cell carcinoma (NCIt:C3493)"]
    unresolved = decisions.canonical_model_id.eq("") | ~decisions.expression_eligible.eq("1")
    decisions.loc[unresolved, ["eligibility_status", "exclusion_reason"]] = ["excluded_unresolved_identity", "v1_unresolved_or_expression_ineligible_identity"]
    decisions["included_v2"] = decisions.eligibility_status.str.startswith("included_")
    vocab = {}
    for value in sorted(set(model.OncotreeLineage)):
        if value:
            assert value.casefold() not in vocab or vocab[value.casefold()] == value
            vocab[value.casefold()] = value
    sanger_models = read("data_processed/sanger_models.csv")
    for value in sorted(set(sanger_models.sanger_tissue)):
        if value:
            vocab.setdefault(value.casefold(), value)
    decisions["tissue"] = decisions.tissue.map(lambda value: vocab.get(value.casefold(), value))
    decisions["tissue_casefold_normalized"] = decisions.tissue.ne(decisions.v1_tissue)
    decisions["cohort_version"] = "v2"
    assert not decisions.expression_row_index.duplicated().any()
    eligible_by_id = decisions[decisions.included_v2].set_index("canonical_model_id")
    assert eligible_by_id.index.is_unique

    def overlay(name):
        original = read("results/" + name)
        included = original[original.canonical_model_id.isin(eligible_by_id.index)].copy()
        excluded = original[~original.canonical_model_id.isin(eligible_by_id.index)].copy()
        included["v1_tissue"] = included.tissue
        included["tissue"] = included.canonical_model_id.map(eligible_by_id.tissue)
        included["eligibility_status"] = included.canonical_model_id.map(eligible_by_id.eligibility_status)
        included["evidence_source"] = included.canonical_model_id.map(eligible_by_id.evidence_source)
        included["cohort_version"] = "v2"
        excluded["exclusion_reason"] = excluded.canonical_model_id.map(decisions.set_index("canonical_model_id").exclusion_reason)
        return original, included, excluded

    dev0, dev, dev_excluded = overlay("development_sample_manifest.csv")
    ext0, ext, ext_excluded = overlay("external_primary_sample_manifest.csv")
    sens0, sensitivity, sens_excluded = overlay("external_sensitivity_sample_manifest.csv")
    assert dev.canonical_model_id.nunique() == 626 and ext.canonical_model_id.nunique() == 175
    assert not (set(dev.biological_group_id) & set(ext.biological_group_id))
    groups = read("data_processed/canonical_model_groups.csv").set_index("canonical_model_id").biological_group_id.to_dict()
    gdsc = read("data_processed/gdsc2_cell_mapping.csv")
    all_gdsc_models = {value for values in gdsc.candidate_model_ids for value in values.split("|") if value}
    all_gdsc_groups = {groups[value] for value in all_gdsc_models if value in groups}
    assert not (all_gdsc_models & set(ext.canonical_model_id))
    assert not (all_gdsc_groups & set(ext.biological_group_id))
    assert not (all_gdsc_groups & set(sensitivity.biological_group_id))
    panel = read("results/main_drugs.csv")
    assert panel.gdsc_drug_id.tolist() == ["1931", "1036", "1191", "1930", "1199", "1091", "1010", "1922"]
    count_updates = {}
    for _, row in panel.iterrows():
        drug = row.gdsc_drug_id
        d = dev[dev.gdsc_drug_id.eq(drug)]
        e = ext[ext.gdsc_drug_id.eq(drug)]
        s = sensitivity[sensitivity.gdsc_drug_id.eq(drug)]
        assert d.biological_group_id.nunique() >= 200 and e.biological_group_id.nunique() >= 80
        assert set(d.min_conc) == {row.gdsc_min_conc} and set(d.max_conc) == {row.gdsc_max_conc}
        count_updates[drug] = {
            "n_development_models": d.canonical_model_id.nunique(), "n_development_groups": d.biological_group_id.nunique(), "n_development_tissues": d.tissue.nunique(),
            "n_external_primary_models": e.canonical_model_id.nunique(), "n_external_primary_groups": e.biological_group_id.nunique(), "n_external_primary_tissues": e.tissue.nunique(),
            "n_external_primary_curve_records": int(pd.to_numeric(e.n_source_experiments).sum()),
            "n_external_all_identity_status_models": s.canonical_model_id.nunique(),
            "n_natural_external_models": e[e.external_subtype.eq("natural_CTRP_only")].canonical_model_id.nunique(),
            "n_reserved_shared_models": e[~e.external_subtype.eq("natural_CTRP_only")].canonical_model_id.nunique(),
        }
    for column in next(iter(count_updates.values())):
        panel["v1_" + column] = panel[column]
        panel[column] = panel.gdsc_drug_id.map(lambda drug: count_updates[drug][column])
    panel["coverage_count_scope"] = "cancer_eligible_v2_selected_manifest"
    panel["cohort_version"] = "v2"
    global_split = decisions[dictionary.columns.tolist() + ["v1_tissue", "eligibility_status", "evidence_source", "disease_annotation", "exclusion_reason", "included_v2", "tissue_casefold_normalized", "cohort_version"]].copy()
    global_split["v1_global_role"] = global_split.global_role
    global_split["v1_in_selected_development"] = global_split.in_selected_development
    global_split["v1_in_selected_external_primary"] = global_split.in_selected_external_primary
    global_split.loc[~global_split.included_v2, "global_role"] = global_split.loc[~global_split.included_v2, "eligibility_status"] + "_v2"
    global_split["in_selected_development"] = global_split.canonical_model_id.isin(set(dev.canonical_model_id))
    global_split["in_selected_external_primary"] = global_split.canonical_model_id.isin(set(ext.canonical_model_id))
    global_split["cancer_eligible_v2"] = global_split.included_v2
    exclusions = pd.concat([dev_excluded.assign(source_manifest="development"), ext_excluded.assign(source_manifest="external_primary"), sens_excluded.assign(source_manifest="external_sensitivity")], ignore_index=True)
    summary = {
        "schema_version": "2.0", "cohort_version": "v2", "parent_cohort_version": "v1",
        "status": "snapshot_pending_independent_validation_lock",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "drug_panel_and_rank_unchanged": True, "drug_and_assay_identity_unchanged": True, "dose_windows_unchanged": True,
        "development_unique_models": dev.canonical_model_id.nunique(), "development_unique_groups": dev.biological_group_id.nunique(),
        "external_unique_models": ext.canonical_model_id.nunique(), "external_unique_groups": ext.biological_group_id.nunique(),
        "development_drug_model_observations": len(dev), "external_drug_model_observations": len(ext),
        "external_sensitivity_unique_models": sensitivity.canonical_model_id.nunique(), "external_sensitivity_drug_model_observations": len(sensitivity),
        "removed_development_unique_models": dev_excluded.canonical_model_id.nunique(), "removed_external_primary_unique_models": ext_excluded.canonical_model_id.nunique(), "removed_external_sensitivity_unique_models": sens_excluded.canonical_model_id.nunique(),
        "unknown_primary_external_models_remaining": 0,
        "verified_authoritative_override": {"canonical_model_id": "ACH-001078", "ccle_id": "HCC1588_LUNG", "cellosaurus_id": "CVCL_A351", "category": "Cancer cell line", "disease": "Lung squamous cell carcinoma", "status": "included_v2"},
        "model_overlap": len(set(dev.canonical_model_id) & set(ext.canonical_model_id)), "biological_group_overlap": len(set(dev.biological_group_id) & set(ext.biological_group_id)),
        "external_overlap_with_all_known_GDSC_candidate_groups": len(set(ext.biological_group_id) & all_gdsc_groups),
        "minimum_external_primary_groups_across_eight_drugs": min(v["n_external_primary_groups"] for v in count_updates.values()),
        "tissue_capitalization_normalized_global_rows": int(decisions.tissue_casefold_normalized.sum()),
        "response_values_read": False, "response_indices_read": False, "expression_distributions_read": False,
        "historical_CTRP_use_disclosed": True, "model_performance_used_for_eligibility": False,
        "formal_model_configuration_frozen": False,
    }
    OUT.mkdir(parents=True)
    dump(RULE_PATH, rules)
    csv(panel, "main_drugs.csv")
    csv(dev, "development_sample_manifest.csv")
    csv(ext, "external_primary_sample_manifest.csv")
    csv(sensitivity, "external_sensitivity_sample_manifest.csv")
    csv(global_split, "global_split_manifest.csv")
    columns = ["expression_row_index", "ccle_id", "canonical_model_id", "biological_group_id", "v1_tissue", "tissue", "OncotreePrimaryDisease", "OncotreeSubtype", "eligibility_status", "included_v2", "exclusion_reason", "evidence_source", "disease_annotation", "tissue_casefold_normalized", "cohort_version"]
    csv(decisions[columns], "eligibility_decisions.csv")
    csv(exclusions, "excluded_sample_manifest.csv")
    for name, frame in [("development", dev), ("external_primary", ext)]:
        tissue = frame.groupby(["gdsc_drug_id", "drug_name", "tissue"], sort=True).agg(n_models=("canonical_model_id", "nunique"), n_groups=("biological_group_id", "nunique")).reset_index()
        csv(tissue, name + "_tissue_counts.csv")
    dump(OUT / "cohort_summary.json", summary)
    direct = PARENT_INPUTS + ["config/cohort_freeze_v1.json", "config/cohort_freeze_v1.lock.json", "src/freeze_cancer_cohort.py", "config/cancer_cohort_v2_rules.json"]
    direct += [path.relative_to(ROOT).as_posix() for path in sorted(CACHE.iterdir()) if path.is_file()]
    direct += [path.relative_to(ROOT).as_posix() for path in sorted(OUT.iterdir()) if path.is_file()]
    files = [{"path": name, "scope": "v2_direct_dependency", "bytes": (ROOT / name).stat().st_size, "sha256": sha(ROOT / name)} for name in direct]
    inherited_names = ["data_processed/expression_tpm.npy", "data_processed/expression_genes.tsv", "data_processed/expression_samples.tsv", "data_manifest/original_source_sha256.csv", "src/project_paths.py"]
    inherited = [{**recorded[name], "scope": "inherited_v1_scientific_input", "verified_by_this_builder": False, "reference_note": "Inherited from locked v1; not opened by metadata-only v2 builder. Validate hashes when development data are accessed."} for name in inherited_names]
    dump(FREEZE_PATH, {
        "schema_version": "2.0", "freeze_type": "cancer_eligibility_cohort_v2", "cohort_version": "v2",
        "created_utc": summary["created_utc"], "status": "snapshot_pending_independent_validation_lock",
        "parent_cohort_freeze": "config/cohort_freeze_v1.json", "parent_cohort_freeze_sha256": parent_lock["freeze_sha256"],
        "parent_lock": "config/cohort_freeze_v1.lock.json", "parent_lock_sha256": sha(ROOT / "config/cohort_freeze_v1.lock.json"),
        "rules": rules, "summary": summary, "files": files + inherited,
        "dependency_scope": "Only scientific metadata, inherited locked scientific inputs, cohort-generation logic and output manifests; mutable plans/reports excluded.",
        "training_configuration_frozen": False, "external_response_accessed": False,
        "v1_files_overwritten": False, "independent_validation_required": True,
    })
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
