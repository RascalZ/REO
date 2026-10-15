"""Build manuscript tables from audited aggregate outputs only.

The script does not read individual response values, expression matrices,
sample-level predictions, or model checkpoints. It records source and output
SHA-256 digests beside the generated CSV tables.
"""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tables" / "manuscript_v1"
DRUGS_PATH = "results/cohort_v2/main_drugs.csv"
COHORT_PATH = "results/cohort_v2/public_cohort_snapshot.json"
DEVELOPMENT_PATH = "results/ordinal_repeated/paired_effect_intervals.csv"
MODEL_PATH = "results/p6_preflight/final_model_manifest.jsonl"
EXTERNAL_PATH = "results/p6_external/evaluation_summary.json"
LOTO_PATH = "results/tissue_loto/contrasts_by_drug_tissue.csv"
LOTO_CONFIG_PATH = "config/tissue_sensitivity_loto_v1.json"
DEVELOPMENT_MODELS_PATH = "results/ordinal_repeated/model_intervals.csv"
REPEAT_STABILITY_PATH = "results/ordinal_repeated/repeat_stability.csv"
ZERO_AUDIT_PATH = "results/ordinal_repeated/molecular_and_failure_summary.csv"
SOURCES = [DRUGS_PATH, COHORT_PATH, DEVELOPMENT_PATH, MODEL_PATH,
           EXTERNAL_PATH, LOTO_PATH, LOTO_CONFIG_PATH,
           DEVELOPMENT_MODELS_PATH, REPEAT_STABILITY_PATH, ZERO_AUDIT_PATH]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_model_manifest() -> pd.DataFrame:
    rows = [json.loads(line) for line in (ROOT / MODEL_PATH).read_text(encoding="utf-8").splitlines() if line.strip()]
    result = pd.DataFrame({
        "gdsc_drug_id": [str(r["task"]["drug_id"]) for r in rows],
        "model": [r["task"]["model"] for r in rows],
        "n_nonzero_molecular_coeffs": [int(r["n_nonzero_molecular_coeffs"]) for r in rows],
    })
    if result.duplicated(["gdsc_drug_id", "model"]).any():
        raise ValueError("Final model manifest contains duplicate drug/model entries")
    return result


def build_table1(drugs: pd.DataFrame, cohort: dict, external: dict) -> pd.DataFrame:
    table = drugs.sort_values("primary_drug_rank").copy()
    table = table.rename(columns={
        "gdsc_drug_id": "GDSC2_compound_ID",
        "gdsc_drug_name": "Compound",
        "ctrp_compound_id": "CTRPv2_compound_ID",
        "n_development_groups": "GDSC2_independent_groups",
        "n_development_models": "GDSC2_drug_model_responses",
        "n_external_primary_groups": "CTRPv2_shared_groups",
        "n_external_primary_models": "CTRPv2_valid_drug_model_responses",
        "n_external_primary_curve_records": "CTRPv2_eligible_curve_records",
    })
    columns = ["GDSC2_compound_ID", "Compound", "CTRPv2_compound_ID",
               "GDSC2_independent_groups", "GDSC2_drug_model_responses",
               "CTRPv2_shared_groups", "CTRPv2_valid_drug_model_responses",
               "CTRPv2_eligible_curve_records"]
    table = table[columns]
    counts = cohort["cohort_counts"]
    total = {column: "Across-drug unique total" for column in columns}
    total.update({
        "GDSC2_compound_ID": "",
        "Compound": "All eight compounds",
        "CTRPv2_compound_ID": "",
        "GDSC2_independent_groups": counts["development_unique_groups"],
        "GDSC2_drug_model_responses": counts["development_drug_model_observations"],
        "CTRPv2_shared_groups": counts["external_unique_groups"],
        "CTRPv2_valid_drug_model_responses": counts["external_drug_model_observations"],
        "CTRPv2_eligible_curve_records": external["response_audit"]["selected_source_rows"],
    })
    return pd.concat([table, pd.DataFrame([total])], ignore_index=True)


def build_table2(drugs: pd.DataFrame, dev: pd.DataFrame, external: dict,
                 model_manifest: pd.DataFrame) -> pd.DataFrame:
    order = drugs.sort_values("primary_drug_rank")
    ids = [str(x) for x in order["gdsc_drug_id"].tolist()]
    names = dict(zip(order["gdsc_drug_id"].astype(str), order["gdsc_drug_name"]))
    dev_rows = dev[(dev["contrast"] == "primary:association_pair-strong_expression") &
                   (dev["metric"] == "spearman") &
                   (dev["family_name"] == "primary_fixed9")].copy()
    dev_rows["gdsc_drug_id"] = dev_rows["gdsc_drug_id"].astype(str)
    dev_by_id = dev_rows.set_index("gdsc_drug_id")
    ext = external["evaluation"]["contrasts"]["primary"]["inference"]
    zero = model_manifest.pivot(index="gdsc_drug_id", columns="model",
                                values="n_nonzero_molecular_coeffs")
    records = []
    for drug_id in ids:
        d = dev_by_id.loc[drug_id]
        e = ext["per_drug"][drug_id]
        interval = e["delta_interval"]
        exp_nonzero = int(zero.loc[drug_id, "expression_strong"])
        pair_nonzero = int(zero.loc[drug_id, "response_association_pair_full"])
        records.append({
            "Compound": names[drug_id], "GDSC2_compound_ID": drug_id,
            "Development_delta_rho_pair_minus_strong_expression": float(d["point"]),
            "Development_family9_lower": float(d["family_adjusted_conditional_lower"]),
            "Development_family9_upper": float(d["family_adjusted_conditional_upper"]),
            "CTRPv2_shared_groups": int(e["n_common_groups"]),
            "CTRPv2_delta_rho_pair_minus_strong_expression": float(e["delta_rho"]),
            "CTRPv2_family9_lower": float(interval["family_adjusted_lower"]),
            "CTRPv2_family9_upper": float(interval["family_adjusted_upper"]),
            "Final_expression_strong_nonzero_coefficients": exp_nonzero,
            "Final_REO_pair_nonzero_coefficients": pair_nonzero,
            "Both_final_models_molecular_zero": exp_nonzero == 0 and pair_nonzero == 0,
        })
    dm = dev_by_id.loc["macro_equal_fixed_eight_drugs"]
    em = ext["strict_macro"]
    eci = em["delta_interval"]
    records.append({
        "Compound": "Equal-drug macro", "GDSC2_compound_ID": "",
        "Development_delta_rho_pair_minus_strong_expression": float(dm["point"]),
        "Development_family9_lower": float(dm["family_adjusted_conditional_lower"]),
        "Development_family9_upper": float(dm["family_adjusted_conditional_upper"]),
        "CTRPv2_shared_groups": "8 drugs",
        "CTRPv2_delta_rho_pair_minus_strong_expression": float(em["delta_rho"]),
        "CTRPv2_family9_lower": float(eci["family_adjusted_lower"]),
        "CTRPv2_family9_upper": float(eci["family_adjusted_upper"]),
        "Final_expression_strong_nonzero_coefficients": "",
        "Final_REO_pair_nonzero_coefficients": "",
        "Both_final_models_molecular_zero": "Not applicable",
    })
    return pd.DataFrame(records)


def build_table_s1(loto: pd.DataFrame, config: dict) -> pd.DataFrame:
    display_tissues = set(config["expected_primary_display_tissues"])
    out = loto[loto["heldout_tissue"].isin(display_tissues)].copy()
    out["drug_id"] = out["drug_id"].astype(str)
    out = out.rename(columns={
        "drug_id": "GDSC2_compound_ID", "drug_name": "Compound",
        "heldout_tissue": "Held_out_tissue", "n_samples": "N_groups",
        "pair_minus_expression_strong_delta_rho": "REO_pair_minus_strong_expression_delta_rho",
        "pair_minus_expression_strong_defined": "Delta_rho_defined",
        "pair_rho": "REO_pair_rho", "expression_strong_rho": "Strong_expression_rho",
        "pair_zero_molecular_solution": "REO_pair_zero_molecular_solution",
    })
    cols = ["GDSC2_compound_ID", "Compound", "Held_out_tissue", "N_groups",
            "REO_pair_rho", "Strong_expression_rho",
            "REO_pair_minus_strong_expression_delta_rho", "Delta_rho_defined",
            "REO_pair_zero_molecular_solution"]
    out = out[cols].sort_values(["GDSC2_compound_ID", "Held_out_tissue"]).reset_index(drop=True)
    if len(out) != 80 or int(out["Delta_rho_defined"].sum()) != 34:
        raise ValueError("LOTO table fails the frozen 80-cell/34-defined-cell contract")
    return out


def build_table_s2(external: dict, drugs: pd.DataFrame, model_manifest: pd.DataFrame) -> pd.DataFrame:
    names = dict(zip(drugs["gdsc_drug_id"].astype(str), drugs["gdsc_drug_name"]))
    zeros = model_manifest.pivot(index="gdsc_drug_id", columns="model",
                                 values="n_nonzero_molecular_coeffs")
    records = []
    contrasts = external["evaluation"]["contrasts"]
    for contrast_name, contrast in contrasts.items():
        inference = contrast["inference"]
        for drug_id, row in inference["per_drug"].items():
            ci = row["delta_interval"]
            candidate = contrast["candidate"]
            reference = contrast["reference"]
            candidate_nonzero = int(zeros.loc[str(drug_id), candidate])
            reference_nonzero = int(zeros.loc[str(drug_id), reference])
            records.append({
                "Contrast_family": contrast_name, "Candidate_model": candidate,
                "Reference_model": reference, "Estimand": "Drug-specific",
                "Compound": names[str(drug_id)], "GDSC2_compound_ID": str(drug_id),
                "N_shared_groups": int(row["n_common_groups"]),
                "Candidate_rho": float(row["rho_a"]), "Reference_rho": float(row["rho_b"]),
                "Delta_rho": float(row["delta_rho"]),
                "Nominal95_lower": float(ci["nominal_95_lower"]),
                "Nominal95_upper": float(ci["nominal_95_upper"]),
                "Family_adjusted_lower": float(ci["family_adjusted_lower"]),
                "Family_adjusted_upper": float(ci["family_adjusted_upper"]),
                "Family_size": int(ci["family_size"]), "Bootstrap_valid": int(ci["bootstrap_valid"]),
                "Candidate_nonzero_coefficients_final_model": candidate_nonzero,
                "Reference_nonzero_coefficients_final_model": reference_nonzero,
                "Both_final_models_molecular_zero": candidate_nonzero == 0 and reference_nonzero == 0,
                "Interval_scope": ci["interval_scope"], "Undefined_reason": row["undefined_reason"],
            })
        macro = inference["strict_macro"]
        ci = macro["delta_interval"]
        records.append({
            "Contrast_family": contrast_name, "Candidate_model": contrast["candidate"],
            "Reference_model": contrast["reference"], "Estimand": "Equal-drug macro",
            "Compound": "8-drug macro", "GDSC2_compound_ID": "", "N_shared_groups": "8 drugs",
            "Candidate_rho": "", "Reference_rho": "", "Delta_rho": float(macro["delta_rho"]),
            "Nominal95_lower": float(ci["nominal_95_lower"]),
            "Nominal95_upper": float(ci["nominal_95_upper"]),
            "Family_adjusted_lower": float(ci["family_adjusted_lower"]),
            "Family_adjusted_upper": float(ci["family_adjusted_upper"]),
            "Family_size": int(ci["family_size"]), "Bootstrap_valid": int(ci["bootstrap_valid"]),
            "Candidate_nonzero_coefficients_final_model": "",
            "Reference_nonzero_coefficients_final_model": "",
            "Both_final_models_molecular_zero": "Not applicable",
            "Interval_scope": ci["interval_scope"], "Undefined_reason": macro["undefined_reason"],
        })
    table = pd.DataFrame(records)
    if len(table) != 27:
        raise ValueError(f"Expected 27 external estimands, observed {len(table)}")
    return table


def build_table_s3(dev: pd.DataFrame, drugs: pd.DataFrame) -> pd.DataFrame:
    names = dict(zip(drugs["gdsc_drug_id"].astype(str), drugs["gdsc_drug_name"]))
    keep = dev[dev["family_name"].isin(["primary_fixed9", "representation_fixed18"]) &
               dev["metric"].eq("spearman")].copy()
    if len(keep) != 27:
        raise ValueError(f"Expected 27 internal inferential estimands, observed {len(keep)}")
    keep["gdsc_drug_id"] = keep["gdsc_drug_id"].astype(str)
    keep["Compound"] = keep["gdsc_drug_id"].map(names).fillna("Equal-drug macro")
    keep["N_groups"] = keep["gdsc_drug_id"].map(
        dict(zip(drugs["gdsc_drug_id"].astype(str), drugs["n_development_groups"].astype(int)))).astype("object")
    keep.loc[keep["gdsc_drug_id"].str.startswith("macro_"), "N_groups"] = "8-drug macro"
    cols = ["family_name", "contrast", "candidate", "reference", "gdsc_drug_id", "Compound",
            "N_groups", "point", "conditional_95_lower", "conditional_95_upper",
            "family_adjusted_conditional_lower", "family_adjusted_conditional_upper",
            "family_size", "bootstrap_valid", "ci_status", "interval_scope", "positive_favors_candidate_or_positive_interaction"]
    return keep[cols].rename(columns={
        "family_name": "Inference_family", "candidate": "Candidate_model", "reference": "Reference_model",
        "gdsc_drug_id": "GDSC2_compound_ID", "point": "Delta_rho",
        "conditional_95_lower": "Nominal95_lower", "conditional_95_upper": "Nominal95_upper",
        "family_adjusted_conditional_lower": "Family_adjusted_lower",
        "family_adjusted_conditional_upper": "Family_adjusted_upper",
    }).sort_values(["Inference_family", "GDSC2_compound_ID", "contrast"]).reset_index(drop=True)


def build_table_s4(model_intervals: pd.DataFrame, drugs: pd.DataFrame) -> pd.DataFrame:
    out = model_intervals.copy()
    out["gdsc_drug_id"] = out["gdsc_drug_id"].astype(str)
    names = dict(zip(drugs["gdsc_drug_id"].astype(str), drugs["gdsc_drug_name"]))
    out["Compound"] = out["gdsc_drug_id"].map(names).fillna("Equal-drug macro")
    out = out.rename(columns={
        "gdsc_drug_id": "GDSC2_compound_ID", "model": "Model_arm", "metric": "Metric",
        "n_unique_biological_groups": "N_groups", "point": "Point_estimate",
        "conditional_95_lower": "Conditional95_lower", "conditional_95_upper": "Conditional95_upper",
        "family_adjusted_conditional_lower": "Family_adjusted_lower",
        "family_adjusted_conditional_upper": "Family_adjusted_upper",
        "bootstrap_valid": "Bootstrap_valid", "interval_scope": "Interval_scope",
    })
    cols = ["GDSC2_compound_ID", "Compound", "Model_arm", "Metric", "N_groups",
            "Point_estimate", "Conditional95_lower", "Conditional95_upper",
            "Family_adjusted_lower", "Family_adjusted_upper", "family_size",
            "Bootstrap_valid", "ci_status", "undefined_reason", "Interval_scope"]
    return out[cols].sort_values(["GDSC2_compound_ID", "Model_arm", "Metric"]).reset_index(drop=True)


def build_table_s5(repeat_stability: pd.DataFrame, drugs: pd.DataFrame) -> pd.DataFrame:
    out = repeat_stability.copy()
    out["gdsc_drug_id"] = out["gdsc_drug_id"].astype(str)
    names = dict(zip(drugs["gdsc_drug_id"].astype(str), drugs["gdsc_drug_name"]))
    out["Compound"] = out["gdsc_drug_id"].map(names).fillna("Equal-drug macro")
    out = out.rename(columns={
        "gdsc_drug_id": "GDSC2_compound_ID", "comparison": "Comparison", "metric": "Metric",
        "repeat_values": "Three_repeat_values", "repeat_mean": "Repeat_mean",
        "repeat_min": "Repeat_min", "repeat_max": "Repeat_max", "repeat_range": "Repeat_range",
        "repeat_sample_sd": "Repeat_sample_SD_descriptive", "positive_repeats": "Positive_repeats",
        "negative_repeats": "Negative_repeats", "zero_repeats": "Zero_repeats",
        "undefined_repeats": "Undefined_repeats", "scope": "Interpretation_scope",
    })
    return out


def build_table_s6(zero_audit: pd.DataFrame, model_manifest: pd.DataFrame, drugs: pd.DataFrame) -> pd.DataFrame:
    out = zero_audit.groupby(["gdsc_drug_id", "model"], as_index=False).agg(
        repeats=("repeat", "nunique"), outer_fits=("outer_fits", "sum"),
        zero_molecular_fits=("zero_molecular_fits", "sum"),
        nonzero_molecular_fits=("nonzero_molecular_fits", "sum"),
        warning_outer_fits=("warning_outer_fits", "sum"), fallback_outer_fits=("fallback_outer_fits", "sum"),
        active_feature_min=("active_feature_min", "min"), active_feature_max=("active_feature_max", "max"),
        unknown_tissue_outer_samples=("outer_unknown_tissue_samples", "sum"),
    )
    out["gdsc_drug_id"] = out["gdsc_drug_id"].astype(str)
    final = model_manifest.rename(columns={"n_nonzero_molecular_coeffs": "full_development_nonzero_coefficients"})
    out = out.merge(final, left_on=["gdsc_drug_id", "model"], right_on=["gdsc_drug_id", "model"],
                    how="left", validate="one_to_one")
    out.loc[out["model"].eq("tissue_null"), "full_development_nonzero_coefficients"] = pd.NA
    names = dict(zip(drugs["gdsc_drug_id"].astype(str), drugs["gdsc_drug_name"]))
    out["Compound"] = out["gdsc_drug_id"].map(names)
    out = out.rename(columns={"gdsc_drug_id": "GDSC2_compound_ID", "model": "Model_arm"})
    cols = ["GDSC2_compound_ID", "Compound", "Model_arm", "repeats", "outer_fits",
            "zero_molecular_fits", "nonzero_molecular_fits", "warning_outer_fits",
            "fallback_outer_fits", "active_feature_min", "active_feature_max",
            "unknown_tissue_outer_samples", "full_development_nonzero_coefficients"]
    return out[cols].sort_values(["GDSC2_compound_ID", "Model_arm"]).reset_index(drop=True)


def build_table_s7(loto: pd.DataFrame) -> pd.DataFrame:
    out = loto.rename(columns={
        "drug_id": "GDSC2_compound_ID", "drug_name": "Compound", "heldout_tissue": "Held_out_tissue",
        "n_samples": "N_groups", "pair_minus_expression_strong_delta_rho": "Delta_rho_pair_minus_strong_expression",
        "pair_minus_expression_strong_defined": "Delta_rho_defined", "rank_minus_expression_strong_delta_rho": "Delta_rho_rank_minus_expression",
        "rmse_gain_pair_vs_tissue_null": "RMSE_gain_pair_vs_tissue_reference",
        "rmse_gain_pair_vs_expression_strong": "RMSE_gain_pair_vs_strong_expression",
        "pair_rho": "Pair_rho", "expression_strong_rho": "Strong_expression_rho",
        "tissue_null_rho": "Tissue_reference_rho", "pair_zero_molecular_solution": "Pair_zero_molecular_solution",
    })
    return out.sort_values(["GDSC2_compound_ID", "Held_out_tissue"]).reset_index(drop=True)


def build_table_s8(drugs: pd.DataFrame) -> pd.DataFrame:
    cols = ["primary_drug_rank", "gdsc_drug_id", "gdsc_drug_name", "ctrp_compound_id",
            "ctrp_compound_name", "ctrp_broad_cpd_id", "match_status", "match_evidence", "gdsc_pubchem_cid",
            "gdsc_name_derived_inchikey", "ctrp_inchikey", "gdsc_min_conc", "gdsc_max_conc",
            "ctrp_top_test_conc_umol", "n_development_groups", "n_external_primary_groups",
            "n_external_primary_curve_records"]
    out = drugs[cols].copy()
    return out.rename(columns={
        "primary_drug_rank": "Panel_order", "gdsc_drug_id": "GDSC2_compound_ID",
        "gdsc_drug_name": "GDSC2_name", "ctrp_compound_id": "CTRPv2_compound_ID",
        "ctrp_compound_name": "CTRPv2_name", "ctrp_broad_cpd_id": "CTRPv2_Broad_CPD_ID",
        "match_status": "Identity_status",
        "match_evidence": "Identity_evidence", "gdsc_pubchem_cid": "GDSC2_PubChem_CID",
        "gdsc_name_derived_inchikey": "GDSC2_name_derived_InChIKey", "ctrp_inchikey": "CTRPv2_InChIKey",
        "gdsc_min_conc": "GDSC2_minimum_concentration_uM", "gdsc_max_conc": "GDSC2_maximum_concentration_uM",
        "ctrp_top_test_conc_umol": "CTRPv2_top_test_concentration_uM",
        "n_development_groups": "GDSC2_development_groups", "n_external_primary_groups": "CTRPv2_external_groups",
        "n_external_primary_curve_records": "CTRPv2_curve_records",
    })


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    drugs = pd.read_csv(ROOT / DRUGS_PATH, dtype={"gdsc_drug_id": str, "ctrp_compound_id": str})
    drugs["gdsc_drug_id"] = drugs["gdsc_drug_id"].astype(str)
    cohort = json.loads((ROOT / COHORT_PATH).read_text(encoding="utf-8"))
    dev = pd.read_csv(ROOT / DEVELOPMENT_PATH, dtype={"gdsc_drug_id": str})
    external = json.loads((ROOT / EXTERNAL_PATH).read_text(encoding="utf-8"))
    loto = pd.read_csv(ROOT / LOTO_PATH)
    loto_config = json.loads((ROOT / LOTO_CONFIG_PATH).read_text(encoding="utf-8"))
    model_intervals = pd.read_csv(ROOT / DEVELOPMENT_MODELS_PATH, dtype={"gdsc_drug_id": str})
    repeat_stability = pd.read_csv(ROOT / REPEAT_STABILITY_PATH, dtype={"gdsc_drug_id": str})
    zero_audit = pd.read_csv(ROOT / ZERO_AUDIT_PATH, dtype={"gdsc_drug_id": str})
    models = read_model_manifest()

    outputs = {
        "Table1_cohort_by_drug.csv": build_table1(drugs, cohort, external),
        "Table2_primary_effects.csv": build_table2(drugs, dev, external, models),
        "TableS1_loto_defined_cells.csv": build_table_s1(loto, loto_config),
        "TableS2_external_27_estimands.csv": build_table_s2(external, drugs, models),
        "TableS3_development_27_estimands.csv": build_table_s3(dev, drugs),
        "TableS4_development_model_performance.csv": build_table_s4(model_intervals, drugs),
        "TableS5_repeat_stability.csv": build_table_s5(repeat_stability, drugs),
        "TableS6_zero_solution_audit.csv": build_table_s6(zero_audit, models, drugs),
        "TableS7_loto_all_23_tissues.csv": build_table_s7(loto),
        "TableS8_compound_identity_mapping.csv": build_table_s8(drugs),
    }
    for filename, frame in outputs.items():
        frame.to_csv(OUT / filename, index=False, na_rep="NA", float_format="%.10g")

    manifest = {
        "schema_version": "1.0",
        "status": "first_draft_tables_from_aggregate_outputs",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": "src/build_manuscript_tables_v1.py",
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "individual_response_values_read": False,
        "sample_level_predictions_read": False,
        "source_sha256": {name: sha256(ROOT / name) for name in SOURCES},
        "script_sha256": sha256(ROOT / "src" / "build_manuscript_tables_v1.py"),
        "output_sha256": {name: sha256(OUT / name) for name in outputs},
        "table_rows": {name: int(len(frame)) for name, frame in outputs.items()},
    }
    (OUT / "table_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "built", "output": str(OUT), "tables": manifest["table_rows"]}))


if __name__ == "__main__":
    main()
