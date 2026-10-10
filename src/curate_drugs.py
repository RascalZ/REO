"""Audit drug identities without reading response values or expression data.

Sources are immutable inputs. Downloads and outputs stay under rebuild_v2.
This module is intentionally independent of predictive model performance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import openpyxl
import pandas as pd

from project_paths import REPOSITORY, ORIGINAL_ROOT, original_input, original_label
ROOT = ORIGINAL_ROOT
BUILD = REPOSITORY
RAW = BUILD / "data_raw" / "drug_metadata"
OUT = BUILD / "data_processed"
GDSC_XLSX = original_input('gdsc2')
LEGACY_CTRP = original_input('legacy_ctrp_response')
GDSC_URL = "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/screened_compounds_rel_8.5.csv"
GDSC_DOC = "https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/"
CTRP_URL = "https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset.zip"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def exact_name(name: str) -> str:
    """Case and whitespace normalization only; no punctuation/salt stripping."""
    return " ".join(str(name).casefold().split())


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not fields:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def extract_local_inventory():
    """Only drug identifiers and assay metadata are accessed; AUC is never used."""
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.load_workbook(GDSC_XLSX, read_only=True, data_only=True)
    it = iter(wb.active.values)
    header = next(it)
    wanted = ["DATASET", "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME", "COMPANY_ID", "WEBRELEASE", "MIN_CONC", "MAX_CONC"]
    indices = [header.index(n) for n in wanted]
    seen = defaultdict(int)
    for values in it:
        seen[tuple(values[i] for i in indices)] += 1
    wb.close()
    rows = [dict(zip([s.lower() for s in wanted], k), n_assayed_rows=v) for k, v in seen.items()]
    rows.sort(key=lambda r: (r["drug_id"], str(r["min_conc"]), str(r["max_conc"])))
    write_csv(OUT / "drug_gdsc2_assay_inventory.csv", rows)
    ctrp_counts = defaultdict(int)
    # csv.DictReader does not convert or inspect response values.
    with LEGACY_CTRP.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ctrp_counts[row["Drug"]] += 1
    write_csv(OUT / "drug_legacy_ctrp_name_inventory.csv", [{"legacy_drug_name": k, "n_assayed_rows": v} for k, v in sorted(ctrp_counts.items())])
    report = {"gdsc_unique_drug_ids": len({r["drug_id"] for r in rows}), "gdsc_unique_drug_names": len({r["drug_name"] for r in rows}), "gdsc_metadata_rows": len(rows), "legacy_ctrp_unique_names": len(ctrp_counts), "response_values_used": False}
    print(json.dumps(report, ensure_ascii=False))
    return report


def get_json(url: str, cache_name: str):
    """Cache public official metadata; never reinterpret HTTP failures as absence."""
    RAW.mkdir(parents=True, exist_ok=True)
    path = RAW / cache_name
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "ResearchDrugIdentityAudit/1.0"}), timeout=35) as response:
                result = json.loads(response.read().decode("utf-8"))
            result["_retrieval"] = {"source_url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
            path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
            time.sleep(0.35)
            return result
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                result = {"_retrieval": {"source_url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat()}, "not_found": True}
                path.write_text(json.dumps(result), encoding="utf-8")
                return result
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
        except (TimeoutError, urllib.error.URLError):
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


def query_name_properties(name: str):
    encoded = urllib.parse.quote(name, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/InChIKey,IsomericSMILES,CanonicalSMILES/JSON"
    cache = "pubchem_name_" + hashlib.sha256(name.encode()).hexdigest()[:16] + ".json"
    try:
        return name, get_json(url, cache)
    except Exception as exc:
        return name, {"retrieval_error": str(exc), "source_url": url}


def fetch_candidate_pubchem():
    inv = pd.read_csv(OUT / "drug_gdsc2_assay_inventory.csv")
    meta = pd.read_csv(RAW / "screened_compounds_rel_8.5.csv").fillna("")
    meta = meta[meta.DRUG_ID.isin(inv.drug_id)]
    legacy = pd.read_csv(OUT / "drug_legacy_ctrp_name_inventory.csv").fillna("")
    ctrp_names = {exact_name(v) for v in legacy.legacy_drug_name}
    candidates = []
    for row in meta.to_dict("records"):
        aliases = [row["DRUG_NAME"], *str(row["SYNONYMS"]).split(",")]
        if any(exact_name(a) in ctrp_names for a in aliases if a.strip()):
            candidates.append(row["DRUG_NAME"])
    candidates = sorted(set(candidates))
    results = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for name, result in pool.map(query_name_properties, candidates):
            results[name] = result
    (RAW / "pubchem_candidate_properties_index.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"candidate_names":len(candidates),"retrieval_errors":sum('retrieval_error' in r for r in results.values()),"not_found":sum(r.get('not_found',False) for r in results.values())}))
    return results


def build_identity_map():
    from rdkit import Chem, RDLogger, rdBase
    RDLogger.DisableLog("rdApp.warning")
    inv = pd.read_csv(OUT / "drug_gdsc2_assay_inventory.csv").fillna("")
    meta = pd.read_csv(RAW / "screened_compounds_rel_8.5.csv").fillna("")
    meta = meta[meta.DRUG_ID.isin(inv.drug_id)]
    ctrp = pd.read_csv(RAW / "v20.meta.per_compound.txt", sep="\t").fillna("")
    properties = json.loads((RAW / "pubchem_candidate_properties_index.json").read_text(encoding="utf-8"))
    assert meta.DRUG_ID.is_unique and ctrp.master_cpd_id.is_unique
    assert set(inv.drug_id) == set(meta.DRUG_ID)
    assert (inv[["drug_id", "drug_name"]].drop_duplicates().merge(meta, left_on="drug_id", right_on="DRUG_ID").eval("drug_name == DRUG_NAME")).all()
    compounds = []
    for r in ctrp.to_dict("records"):
        smiles = r["cpd_smiles"]
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        key = Chem.MolToInchiKey(mol) if mol is not None else ""
        unassigned = sum(label == "?" for _, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)) if mol is not None else ""
        unassigned_stereo = sum(str(s.specified) != "Specified" for s in Chem.FindPotentialStereo(mol)) if mol is not None else ""
        nfrag = len(Chem.GetMolFrags(mol)) if mol is not None else ""
        combination = "mol/mol" in r["cpd_name"] or bool(re.search(r"\s\+\s", r["cpd_name"]))
        compounds.append({"ctrp_compound_id": r["master_cpd_id"], "ctrp_compound_name": r["cpd_name"], "ctrp_broad_cpd_id": r["broad_cpd_id"], "ctrp_smiles": smiles,
                          "ctrp_inchikey": key, "structure_parse_status": "valid" if key else "invalid_or_missing", "is_combination": combination, "n_structure_fragments": nfrag,
                          "unassigned_chiral_centres": unassigned, "unassigned_stereo_elements": unassigned_stereo, "ctrp_status_historical": r["cpd_status"], "ctrp_target": r["target_or_activity_of_compound"], "ctrp_top_test_conc_umol": r["top_test_conc_umol"], "source_url": CTRP_URL})
    write_csv(OUT / "drug_ctrp_compound_identity.csv", compounds)
    by_name = defaultdict(list)
    by_key = defaultdict(list)
    for c in compounds:
        by_name[exact_name(c["ctrp_compound_name"])].append(c)
        if c["ctrp_inchikey"] and not c["is_combination"]:
            by_key[c["ctrp_inchikey"]].append(c)
    name_ids = meta.groupby("DRUG_NAME").DRUG_ID.apply(list).to_dict()
    counts = inv.groupby("drug_id").n_assayed_rows.sum().to_dict()
    ranges = inv.groupby("drug_id").size().to_dict()
    rows = []
    for g in meta.to_dict("records"):
        name, gid = g["DRUG_NAME"], g["DRUG_ID"]
        props = properties.get(name, {})
        plist = props.get("PropertyTable", {}).get("Properties", [])
        prop = plist[0] if len(plist) == 1 else {}
        pkey = prop.get("InChIKey", "")
        aliases = [name] + [s.strip() for s in str(g["SYNONYMS"]).split(",") if s.strip()]
        candidates = {}
        for alias in aliases:
            for c in by_name.get(exact_name(alias), []):
                candidates[c["ctrp_compound_id"]] = c
        # Structures augment candidate generation, never override contradictory structures.
        for c in by_key.get(pkey, []):
            candidates[c["ctrp_compound_id"]] = c
        if not candidates:
            candidates = {"": {}}
        for cid, c in candidates.items():
            ckey = c.get("ctrp_inchikey", "")
            evidence = []
            if cid:
                matched_aliases = [a for a in aliases if exact_name(a) == exact_name(c["ctrp_compound_name"])]
                if matched_aliases:
                    evidence.append("official_name_or_synonym:" + "|".join(matched_aliases))
            if not cid:
                status = "unmatched_no_candidate"
            elif c["is_combination"]:
                status = "excluded_combination"
            elif len(plist) > 1:
                status = "ambiguous_multiple_pubchem_cids"
            elif not pkey or not ckey:
                status = "ambiguous_missing_structure"
            elif c["n_structure_fragments"] > 1 or "." in prop.get("SMILES", ""):
                status = "ambiguous_salt_or_multicomponent"
            elif pkey == ckey:
                if c["unassigned_stereo_elements"]:
                    status = "ambiguous_unassigned_stereochemistry"
                else:
                    status = "structure_supported"
                evidence.append("full_standard_inchikey_equal")
            elif pkey.split("-")[0] == ckey.split("-")[0]:
                status = "ambiguous_stereochemistry_or_protonation"
                evidence.append("connectivity_block_equal_full_inchikey_differs")
            else:
                status = "ambiguous_structure_conflict"
                evidence.append("full_inchikey_and_connectivity_differ")
            if pkey:
                evidence.append("gdsc_official_primary_name_to_pubchem_not_gdsc_measured_structure")
            rows.append({"gdsc_drug_id": gid, "gdsc_drug_name": name, "ctrp_compound_id": cid,
                         "ctrp_compound_name": c.get("ctrp_compound_name", ""), "match_status": status,
                         "identity_eligible": status == "structure_supported", "match_evidence": ";".join(evidence),
                         "gdsc_pubchem_cid": prop.get("CID", ""), "gdsc_name_derived_inchikey": pkey,
                         "ctrp_inchikey": ckey, "gdsc_name_derived_smiles": prop.get("SMILES", ""), "ctrp_smiles": c.get("ctrp_smiles", ""),
                         "gdsc_synonyms": g["SYNONYMS"], "gdsc_target": g["TARGET"], "gdsc_target_pathway": g["TARGET_PATHWAY"],
                         "gdsc_same_name_ids": "|".join(map(str, name_ids[name])), "gdsc_duplicate_name": len(name_ids[name]) > 1,
                         "gdsc_n_assayed_rows": counts[gid], "gdsc_metadata_range_rows": ranges[gid],
                         "ctrp_broad_cpd_id": c.get("ctrp_broad_cpd_id", ""), "ctrp_top_test_conc_umol": c.get("ctrp_top_test_conc_umol", ""),
                         "source_url": " | ".join([GDSC_URL, CTRP_URL, props.get("_retrieval", {}).get("source_url", "")]).rstrip(" |")})
    # No arbitrary selection when >1 CTRP IDs share a validated full structure.
    valid_per_gid = defaultdict(list)
    for r in rows:
        if r["identity_eligible"]:
            valid_per_gid[r["gdsc_drug_id"]].append(r)
    for group in valid_per_gid.values():
        if len(group) > 1:
            for r in group:
                r["match_status"] = "ambiguous_multiple_ctrp_ids"
                r["identity_eligible"] = False
    rows.sort(key=lambda r: (r["gdsc_drug_id"], str(r["ctrp_compound_id"])))
    write_csv(OUT / "drug_identity_map.csv", rows)
    write_csv(OUT / "drug_identity_eligible.csv", [r for r in rows if r["identity_eligible"]], list(rows[0]))
    write_csv(OUT / "drug_identity_ambiguous.csv", [r for r in rows if r["match_status"].startswith("ambiguous")], list(rows[0]))
    legacy = pd.read_csv(OUT / "drug_legacy_ctrp_name_inventory.csv").fillna("")
    legacy_rows = []
    for r in legacy.to_dict("records"):
        found = by_name.get(exact_name(r["legacy_drug_name"]), [])
        for c in found or [{}]:
            legacy_rows.append({**r, "ctrp_compound_id": c.get("ctrp_compound_id", ""), "ctrp_compound_name": c.get("ctrp_compound_name", ""), "mapping_status": "exact_name_unique" if len(found) == 1 else "ambiguous_or_unmatched", "is_combination": c.get("is_combination", ""), "source_url": CTRP_URL})
    write_csv(OUT / "drug_legacy_ctrp_name_to_id.csv", legacy_rows)
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(), "rdkit_version": rdBase.rdkitVersion,
               "gdsc_ids": meta.DRUG_ID.nunique(), "gdsc_names": meta.DRUG_NAME.nunique(), "ctrp_ids": len(compounds),
               "ctrp_combinations": sum(c["is_combination"] for c in compounds), "ctrp_multifragment": sum(c["n_structure_fragments"] != "" and c["n_structure_fragments"] > 1 for c in compounds),
               "identity_status_counts": dict(pd.Series([r["match_status"] for r in rows]).value_counts().items()),
               "eligible_gdsc_ids": len({r["gdsc_drug_id"] for r in rows if r["identity_eligible"]}), "eligible_ctrp_ids": len({r["ctrp_compound_id"] for r in rows if r["identity_eligible"]}),
               "legacy_unique_name_mappings": sum(r["mapping_status"] == "exact_name_unique" for r in legacy_rows), "response_values_used": False}
    (OUT / "drug_identity_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest = [{"path": original_label(p), "size_bytes": p.stat().st_size, "sha256": sha256(p)} for p in sorted(RAW.iterdir()) if p.is_file() and p.name != "drug_source_manifest.json"]
    original_inputs = [{"path": original_label(p), "size_bytes": p.stat().st_size, "sha256": sha256(p), "access": "read_only_metadata_columns"} for p in [GDSC_XLSX, LEGACY_CTRP]]
    (RAW / "drug_source_manifest.json").write_text(json.dumps({"created_utc": summary["created_utc"], "gdsc_source_url": GDSC_URL, "ctrp_source_url": CTRP_URL, "ctrp_archive_verification_manifest": "rebuild_v2/data_raw/ctrp/download_manifest.json", "original_inputs": original_inputs, "files": manifest}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))
    return summary


def validate_identity_outputs():
    """Check scientifically material joins and exclusion rules on pinned metadata."""
    inv = pd.read_csv(OUT / "drug_gdsc2_assay_inventory.csv")
    mapping = pd.read_csv(OUT / "drug_identity_map.csv").fillna("")
    eligible = mapping[mapping.identity_eligible.astype(str).str.lower() == "true"]
    compounds = pd.read_csv(OUT / "drug_ctrp_compound_identity.csv").fillna("")
    assert set(mapping.gdsc_drug_id) == set(inv.drug_id), "GDSC ID silently dropped"
    assert not mapping.duplicated(["gdsc_drug_id", "ctrp_compound_id"]).any(), "Duplicate mapping edge"
    assert eligible.gdsc_drug_id.is_unique, "More than one eligible CTRP compound per GDSC drug ID"
    assert (eligible.gdsc_name_derived_inchikey == eligible.ctrp_inchikey).all(), "Eligible structure conflict"
    eligible_ct = eligible.merge(compounds, on="ctrp_compound_id", suffixes=("", "_ct"), validate="many_to_one")
    assert (~eligible_ct.is_combination).all(), "Combination admitted to single-agent cohort"
    assert (eligible_ct.n_structure_fragments == 1).all(), "Salt/mixture admitted"
    assert (eligible_ct.unassigned_stereo_elements == 0).all(), "Undefined stereo admitted"
    # Drug-specific sentinel checks address the known old manuscript identity problem.
    lapatinib = eligible[eligible.gdsc_drug_id == 1558]
    assert len(lapatinib) == 1 and int(lapatinib.iloc[0].ctrp_compound_id) == 634309
    assert not set([1007, 1819]).intersection(eligible.gdsc_drug_id), "Docetaxel stereo conflict silently admitted"
    assert not len(eligible[(eligible.gdsc_drug_id == 1736) & (eligible.ctrp_compound_id == 660306)]), "Erroneous AZD1480 synonym admitted"
    legacy = pd.read_csv(OUT / "drug_legacy_ctrp_name_to_id.csv")
    assert (legacy.mapping_status == "exact_name_unique").all()
    print(json.dumps({"identity_checks": "PASS", "eligible_mapping_edges": len(eligible), "original_response_values_examined": False}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--pubchem-candidates", action="store_true")
    parser.add_argument("--build-map", action="store_true")
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    if args.self_check:
        validate_identity_outputs()
    elif args.build_map:
        build_identity_map()
    elif args.pubchem_candidates:
        fetch_candidate_pubchem()
    else:
        extract_local_inventory()
