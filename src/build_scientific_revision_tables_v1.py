"""Reformat existing public aggregates for scientific revision R1.

Only the four explicit aggregate inputs below are read. This script does not
load responses, individual predictions, bootstrap draws/weights, or model
artifacts. Strict macro points are arithmetic means of all eight pre-existing
drug points; no new fitting, resampling, intervals, or inferential families are
introduced. Previously published tables are never overwritten.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    "external": "results/p6_external/evaluation_summary.json",
    "final_models": "results/p6_preflight/final_model_manifest.jsonl",
    "development_contrasts": "results/ordinal_repeated/paired_effect_intervals.csv",
    "development_models": "results/ordinal_repeated/model_intervals.csv",
}
DRUGS = (
    (1931, "MIRA-1"), (1036, "PLX-4720"), (1191, "Bortezomib"),
    (1930, "Telomerase Inhibitor IX"), (1199, "Tamoxifen"),
    (1091, "BMS-536924"), (1010, "Gefitinib"), (1922, "Cediranib"),
)
MODELS = (
    "tissue_null", "expression_strong", "mad_expression", "mad_rank",
    "mad_pair_full", "response_association_expression",
    "response_association_rank", "response_association_pair_full",
)
RULES = ("mad", "response_association")
COMPARISONS = (("rank", "expression"), ("pair_full", "expression"),
               ("pair_full", "rank"))
MACRO = "macro_equal_fixed_eight_drugs"
EXTERNAL_SCOPE = "conditional_on_fixed_trained_models_and_external_biological_groups"


def number(value):
    if value in (None, "", "NA", "NaN", "nan"):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def strict_mean(values):
    values = list(values)
    if len(values) != len(DRUGS) or any(value is None for value in values):
        return None
    return math.fsum(values) / len(DRUGS)


def require_close(a, b, label):
    if a is None or b is None:
        if a != b:
            raise ValueError("Definedness mismatch: " + label)
    elif not math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("Aggregate mismatch: " + label)


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def csv_records(text):
    return list(csv.DictReader(text.lstrip("\ufeff").splitlines()))


def write_csv(path, rows):
    if not rows or any(list(row) != list(rows[0]) for row in rows):
        raise ValueError("Empty or inconsistent output schema")
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: "NA" if value is None else value for key, value in row.items()}
                         for row in rows)


def existing_interval(interval, source, *, family="descriptive"):
    """Copy an already computed interval without changing coverage or family."""
    if interval is None:
        return {
            "Conditional95_lower": None, "Conditional95_upper": None,
            "Family_adjusted_lower": None, "Family_adjusted_upper": None,
            "Inference_family": "descriptive", "Family_size": None,
            "Bootstrap_requested": None, "Bootstrap_valid": None,
            "CI_status": "not_precomputed_no_new_interval",
            "Interval_scope": "point_only_no_new_interval_or_test",
            "CI_source": "none_existing",
        }
    # Development and external files use different field names.
    return {
        "Conditional95_lower": number(interval.get("conditional_95_lower", interval.get("nominal_95_lower"))),
        "Conditional95_upper": number(interval.get("conditional_95_upper", interval.get("nominal_95_upper"))),
        "Family_adjusted_lower": number(interval.get("family_adjusted_conditional_lower", interval.get("family_adjusted_lower"))),
        "Family_adjusted_upper": number(interval.get("family_adjusted_conditional_upper", interval.get("family_adjusted_upper"))),
        "Inference_family": family,
        "Family_size": number(interval.get("family_size")),
        "Bootstrap_requested": number(interval.get("bootstrap_requested")),
        "Bootstrap_valid": number(interval.get("bootstrap_valid")),
        "CI_status": interval["ci_status"],
        "Interval_scope": interval["interval_scope"],
        "CI_source": source,
    }


def build_tables(contents):
    external = json.loads(contents["external"])
    evaluation = external["evaluation"]
    manifests = [json.loads(line) for line in contents["final_models"].splitlines() if line.strip()]
    dev_contrasts = csv_records(contents["development_contrasts"])
    dev_models = csv_records(contents["development_models"])
    expected = {(drug, model) for drug, _ in DRUGS for model in MODELS}
    arm_records = evaluation["per_arm_descriptive"]
    arms = {(int(row["drug_id"]), row["model"]): row for row in arm_records}
    final = {(int(row["task"]["drug_id"]), row["task"]["model"]): (index, row)
             for index, row in enumerate(manifests, 1)}
    if len(arm_records) != 64 or set(arms) != expected or len(manifests) != 64 or set(final) != expected:
        raise ValueError("The fixed eight-drug/eight-arm aggregate contract is incomplete")
    if external["status"] != "completed_one_time_external_evaluation":
        raise ValueError("External source is not the completed locked evaluation")
    audit = external["identity_audit"]
    if audit["external_unique_models"] != audit["external_unique_groups"]:
        raise ValueError("Model count cannot be presented as independent-group N")
    if any(audit[key] != 0 for key in ("development_model_overlap", "development_group_overlap", "known_candidate_group_overlap")):
        raise ValueError("External identity overlap violates the locked cohort")

    # Harvest absolute rho intervals only where they were already generated.
    absolute_intervals = {}
    for family, contrast in evaluation["contrasts"].items():
        for drug, estimates in contrast["inference"]["per_drug"].items():
            for suffix, arm in (("a", contrast["candidate"]), ("b", contrast["reference"])):
                key = (int(drug), arm)
                interval = estimates["rho_" + suffix + "_interval"]
                require_close(number(interval["point"]), number(arms[key]["rho"]), "external absolute rho")
                source = INPUTS["external"] + "::evaluation.contrasts." + family + ".inference.per_drug." + drug + ".rho_" + suffix + "_interval"
                if key in absolute_intervals:
                    if interval != absolute_intervals[key][0]:
                        raise ValueError("Repeated existing absolute intervals disagree")
                else:
                    absolute_intervals[key] = (interval, source)

    s9 = []
    for model in MODELS:
        for drug, name in DRUGS:
            row = arms[(drug, model)]
            eligible = int(audit["per_drug_eligible"][str(drug)])
            if any(int(row[key]) != eligible for key in ("frozen_eligible_ids", "response_valid_ids", "predicted_on_response_valid_ids")):
                raise ValueError("Incomplete coverage: do not subtract unpaired full-set correlations")
            index, fitted = final[(drug, model)]
            if not fitted["fresh_process_reload_verified"]:
                raise ValueError("Final-model reload audit is incomplete")
            coefficients = int(fitted["n_nonzero_molecular_coeffs"])
            zero = None if model == "tissue_null" else coefficients == 0
            if model != "tissue_null" and int(row["zero_molecular_solution_eligible"]) != (eligible if zero else 0):
                raise ValueError("Existing prediction zero flags disagree with model manifest")
            interval, source = absolute_intervals.get((drug, model), (None, "none_existing"))
            s9.append({
                "Estimand": "Drug-specific", "GDSC2_compound_ID": drug, "Compound": name,
                "Model_arm": model, "N_unique_external_groups": eligible,
                "N_drug_group_response_records": eligible, "N_compounds": 1,
                "Spearman_rho": number(row["rho"]),
                **existing_interval(interval, source),
                "Undefined_reason": row["undefined_reason"],
                "Final_nonzero_molecular_coefficients": None if model == "tissue_null" else coefficients,
                "Final_zero_molecular_solution": zero,
                "N_final_molecular_zero_drug_models": None if model == "tissue_null" else int(zero),
                "Unknown_tissue_drug_group_records": row["unknown_training_tissue_eligible"],
                "Zero_status_source": INPUTS["final_models"] + ":" + str(index),
                "Point_source": INPUTS["external"] + "::evaluation.per_arm_descriptive",
                "Point_aggregation": "unchanged_existing_drug_specific_rho",
            })
        selected = [row for row in s9 if row["Model_arm"] == model and row["Estimand"] == "Drug-specific"]
        point = strict_mean(row["Spearman_rho"] for row in selected)
        s9.append({
            "Estimand": "Strict equal-eight-drug macro", "GDSC2_compound_ID": MACRO,
            "Compound": "Equal-drug macro", "Model_arm": model,
            "N_unique_external_groups": audit["external_unique_groups"],
            "N_drug_group_response_records": sum(row["N_drug_group_response_records"] for row in selected),
            "N_compounds": 8, "Spearman_rho": point,
            **existing_interval(None, "none_existing"),
            "Undefined_reason": "undefined_constituent_drug" if point is None else "",
            "Final_nonzero_molecular_coefficients": None,
            "Final_zero_molecular_solution": None,
            "N_final_molecular_zero_drug_models": None if model == "tissue_null" else sum(row["N_final_molecular_zero_drug_models"] for row in selected),
            "Unknown_tissue_drug_group_records": sum(row["Unknown_tissue_drug_group_records"] for row in selected),
            "Zero_status_source": INPUTS["final_models"] + "::all_eight_records_for_this_arm",
            "Point_source": INPUTS["external"] + "::evaluation.per_arm_descriptive",
            "Point_aggregation": "strict_arithmetic_mean_of_all_eight_existing_drug_rhos_no_new_resampling",
        })

    dev_lookup = {(row["gdsc_drug_id"], row["contrast"]): (index, row)
                  for index, row in enumerate(dev_contrasts, 2) if row["metric"] == "spearman"}
    dev_absolute = {(row["gdsc_drug_id"], row["model"]): row
                    for row in dev_models if row["metric"] == "spearman"}
    descriptive = {(int(row["drug_id"]), row["candidate"], row["reference"]): row
                   for row in evaluation["paired_descriptive"]}
    s10 = []
    for rule in RULES:
        for candidate_rep, reference_rep in COMPARISONS:
            candidate, reference = rule + "_" + candidate_rep, rule + "_" + reference_rep
            label = rule + ":" + candidate_rep + "-" + reference_rep
            for drug, name in (*DRUGS, (MACRO, "Equal-drug macro")):
                index, row = dev_lookup[(str(drug), label)]
                a, b = dev_absolute[(str(drug), candidate)], dev_absolute[(str(drug), reference)]
                if a["n_unique_biological_groups"] != b["n_unique_biological_groups"]:
                    raise ValueError("Development representation comparison has unequal N")
                require_close(number(a["point"]) - number(b["point"]), number(row["point"]), "development representation delta")
                s10.append({
                    "Evaluation_stage": "Repeated_nested_development", "Candidate_rule": rule,
                    "Contrast": label, "Candidate_model": candidate, "Reference_model": reference,
                    "GDSC2_compound_ID": drug, "Compound": name,
                    "N_unique_groups": int(a["n_unique_biological_groups"]),
                    "N_compounds": 8 if drug == MACRO else 1,
                    "Candidate_rho": number(a["point"]), "Reference_rho": number(b["point"]),
                    "Delta_rho_candidate_minus_reference": number(row["point"]),
                    **existing_interval(row, INPUTS["development_contrasts"] + ":" + str(index), family=row["family_name"]),
                    "Undefined_reason": row["undefined_reason"],
                    "Both_final_models_molecular_zero": None,
                    "Point_source": INPUTS["development_contrasts"] + ":" + str(index),
                    "Point_aggregation": "unchanged_existing_mean_of_three_complete_repeat_metrics",
                })
            inferential = evaluation["contrasts"]["representation_" + rule] if candidate_rep == "pair_full" and reference_rep == "expression" else None
            drug_points, candidate_points, reference_points = [], [], []
            for drug, name in DRUGS:
                if inferential is not None:
                    row = inferential["inference"]["per_drug"][str(drug)]
                    interval = row["delta_interval"]
                    point_source = INPUTS["external"] + "::evaluation.contrasts.representation_" + rule + ".inference.per_drug." + str(drug)
                    family = "representation_fixed18"
                else:
                    row = descriptive[(drug, candidate, reference)]
                    interval = None
                    point_source = INPUTS["external"] + "::evaluation.paired_descriptive"
                    family = "descriptive"
                a, b, point = number(row["rho_a"]), number(row["rho_b"]), number(row["delta_rho"])
                if a is not None and b is not None:
                    require_close(a-b, point, "external paired representation delta")
                require_close(a, number(arms[(drug, candidate)]["rho"]), "paired versus per-arm rho")
                require_close(b, number(arms[(drug, reference)]["rho"]), "paired versus per-arm rho")
                n = int(row["n_common_groups"])
                if n != int(audit["per_drug_eligible"][str(drug)]):
                    raise ValueError("External complete paired-coverage contract changed")
                drug_points.append(point); candidate_points.append(a); reference_points.append(b)
                s10.append({
                    "Evaluation_stage": "Locked_CTRPv2", "Candidate_rule": rule,
                    "Contrast": label, "Candidate_model": candidate, "Reference_model": reference,
                    "GDSC2_compound_ID": drug, "Compound": name, "N_unique_groups": n,
                    "N_compounds": 1, "Candidate_rho": a, "Reference_rho": b,
                    "Delta_rho_candidate_minus_reference": point,
                    **existing_interval(interval, point_source + ".delta_interval", family=family),
                    "Undefined_reason": row["undefined_reason"],
                    "Both_final_models_molecular_zero": all(final[(drug, model)][1]["n_nonzero_molecular_coeffs"] == 0 for model in (candidate, reference)),
                    "Point_source": point_source, "Point_aggregation": "unchanged_existing_paired_drug_specific_delta",
                })
            macro = strict_mean(drug_points)
            if inferential is not None:
                row = inferential["inference"]["strict_macro"]
                require_close(macro, number(row["delta_rho"]), "external strict representation macro")
                interval = row["delta_interval"]
                point_source = INPUTS["external"] + "::evaluation.contrasts.representation_" + rule + ".inference.strict_macro"
            else:
                interval = None
                point_source = INPUTS["external"] + "::evaluation.paired_descriptive"
            s10.append({
                "Evaluation_stage": "Locked_CTRPv2", "Candidate_rule": rule,
                "Contrast": label, "Candidate_model": candidate, "Reference_model": reference,
                "GDSC2_compound_ID": MACRO, "Compound": "Equal-drug macro",
                "N_unique_groups": audit["external_unique_groups"], "N_compounds": 8,
                "Candidate_rho": strict_mean(candidate_points), "Reference_rho": strict_mean(reference_points),
                "Delta_rho_candidate_minus_reference": macro,
                **existing_interval(interval, point_source + ".delta_interval", family="representation_fixed18"),
                "Undefined_reason": "undefined_constituent_drug" if macro is None else "",
                "Both_final_models_molecular_zero": None,
                "Point_source": point_source,
                "Point_aggregation": "strict_arithmetic_mean_of_all_eight_existing_paired_drug_deltas_no_new_resampling",
            })
    return s9, s10


def draft_results(s9, s10):
    absolute = {row["Model_arm"]: row for row in s9 if row["GDSC2_compound_ID"] == MACRO}
    compare = {(row["Evaluation_stage"], row["Contrast"]): row
               for row in s10 if row["GDSC2_compound_ID"] == MACRO}
    lines = [
        "# Added results for scientific revision R1 — author-side working draft",
        "", "Generated only from the four public aggregate sources listed in the table manifest. No model fit, sample-level response/prediction access, bootstrap draw/weight access, resampling, new interval, or new inferential family was performed.",
        "", "## Proposed Results addition: external absolute performance", "",
        ("Absolute ranking associations in the locked CTRPv2 evaluation were descriptively non-zero. The strict equal-eight-drug mean Spearman correlation was "
         f"{absolute['tissue_null']['Spearman_rho']:.6f} for the tissue-only reference, "
         f"{absolute['expression_strong']['Spearman_rho']:.6f} for the strong-expression model, and "
         f"{absolute['response_association_pair_full']['Spearman_rho']:.6f} for the response-associated pair model (Table S9). "
         "These descriptive macro points are strict arithmetic means of the eight existing drug-specific correlations; no absolute-macro intervals were precomputed. A near-zero pair-minus-expression contrast therefore does not imply absence of all transferable predictive association. The three drugs with identical zero-molecular primary models retain tissue-based predictions and provide no pair-specific molecular comparison."),
        "", "This is an association/ranking statement, not a test of causal molecular information. Existing absolute-rho conditional intervals are transcribed for the five molecular arms represented in the original inferential comparisons (40 drug–arm rows). Tissue-only, rank-arm, and absolute-macro intervals were not precomputed and remain NA. No interval was inferred from equal predictions or from subtracting interval endpoints.",
        "", "## Proposed Results addition: complete matched-candidate representation contrasts", "",
        "Table S10 makes the prespecified representation comparisons explicit for both development and external evaluation. Only pair-minus-expression carries the existing 18-estimand family adjustment; development rank contrasts retain their existing unadjusted descriptive intervals, and external rank contrasts remain descriptive point estimates without new intervals.",
        "", "| Stage | Candidate rule | Contrast | Δρ | Existing nominal conditional 95% interval | Existing family-adjusted interval |",
        "| --- | --- | --- | ---: | --- | --- |",
    ]
    for stage in ("Repeated_nested_development", "Locked_CTRPv2"):
        for rule in RULES:
            for candidate, reference in COMPARISONS:
                key = rule + ":" + candidate + "-" + reference
                row = compare[(stage, key)]
                def interval(a, b):
                    return "NA (not precomputed)" if row[a] is None else f"[{row[a]:+.6f}, {row[b]:+.6f}]"
                value = row["Delta_rho_candidate_minus_reference"]
                lines.append(f"| {stage} | {rule} | {candidate} − {reference} | {value:+.6f} | "
                             + interval("Conditional95_lower", "Conditional95_upper") + " | "
                             + interval("Family_adjusted_lower", "Family_adjusted_upper") + " |")
    lines += [
        "", "The rank-arm comparisons do not yield a general REO advantage. In development, the response-associated rank-minus-expression macro interval includes zero; the response-associated pair-minus-rank macro point estimate is negative and its existing descriptive interval also includes zero. External rank contrasts have no precomputed uncertainty intervals and cannot support a new significance conclusion. These results neither establish equivalence/non-inferiority nor identify causal shares of amplitude loss, feature expansion, collinearity, or elastic-net regularization.",
        "", "## Table interpretation and boundaries", "",
        "- Table S9 has 72 rows: 64 drug–arm points and eight strict equal-drug macro points. N_unique_external_groups for a macro is the fixed external universe of 175, not 8 independent subjects; the 1,185 drug–group response records are a different denominator.",
        "- Table S10 has 108 rows: two evaluation stages × two candidate rules × three representation contrasts × nine drug/macro estimands. Spearman only; positive Δρ favors the candidate representation.",
        "- NA is retained for unavailable/undefined quantities. Tissue-only coefficient-zero status is not treated as a fitted molecular zero solution. Macro coefficient counts are not summed as if they were one model.",
        "- Development intervals condition on saved OOF predictions and omit training/tuning/refitting uncertainty; external intervals condition on fixed fitted models and observed external groups. The two study families remain separately adjusted and do not guarantee whole-study simultaneous coverage.",
        "- No performance-retention margin was chosen and no equivalence or non-inferiority test was added. The negative/near-zero primary REO advantage conclusion remains unchanged.",
    ]
    return "\n".join(lines) + "\n"


def main():
    contents, provenance = {}, {}
    for key, relative in INPUTS.items():
        path = ROOT / relative
        raw = path.read_bytes()
        contents[key] = raw.decode("utf-8-sig")
        provenance[relative] = {"sha256": sha_bytes(raw), "bytes": len(raw), "scope": "existing_public_aggregate_only"}
    s9, s10 = build_tables(contents)
    if len(s9) != 72 or len(s10) != 108:
        raise ValueError("Unexpected revision table row count")
    if len({(row["GDSC2_compound_ID"], row["Model_arm"]) for row in s9}) != 72:
        raise ValueError("Duplicate S9 identity")
    if len({(row["Evaluation_stage"], row["GDSC2_compound_ID"], row["Contrast"]) for row in s10}) != 108:
        raise ValueError("Duplicate S10 identity")
    if sum(row["Conditional95_lower"] is not None for row in s9) != 40:
        raise ValueError("Existing external absolute-CI coverage changed")
    if sum(row["Conditional95_lower"] is not None for row in s10) != 72:
        raise ValueError("Existing representation-CI coverage changed")
    if sum(row["Family_adjusted_lower"] is not None for row in s10) != 36:
        raise ValueError("Existing family-18 representation coverage changed")
    for row in s10:
        if row["Both_final_models_molecular_zero"] is True:
            require_close(row["Delta_rho_candidate_minus_reference"], 0.0, "identical zero-molecular comparison")
    directory = ROOT / "tables/scientific_revision_v1"
    directory.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for name, rows in (("TableS9_external_absolute_performance.csv", s9),
                       ("TableS10_representation_comparisons.csv", s10)):
        path = directory / name
        write_csv(path, rows)
        outputs[path.relative_to(ROOT).as_posix()] = {
            "sha256": sha_bytes(path.read_bytes()), "rows": len(rows),
            "columns": list(rows[0]),
        }
    draft = ROOT / "docs/scientific_revision_v1/ADDED_RESULTS_DRAFT.md"
    draft.parent.mkdir(parents=True, exist_ok=True)
    draft.write_text(draft_results(s9, s10), encoding="utf-8")
    outputs[draft.relative_to(ROOT).as_posix()] = {"sha256": sha_bytes(draft.read_bytes())}
    manifest = {
        "schema_version": "scientific_revision_v1", "status": "aggregate_only_tables_validated",
        "generator": "src/build_scientific_revision_tables_v1.py",
        "generator_sha256": sha_bytes(Path(__file__).read_bytes()),
        "sources": provenance, "outputs": outputs,
        "contracts": {
            "fixed_drugs": [drug for drug, _ in DRUGS], "fixed_models": list(MODELS),
            "S9_rows": 72, "S9_existing_absolute_rho_intervals": 40,
            "S9_absolute_macro_intervals_not_precomputed": True,
            "S10_rows": 108, "S10_existing_intervals": 72,
            "S10_existing_family18_adjusted_intervals": 36,
            "S10_external_rank_point_only_rows": 36,
            "strict_macro_requires_all_eight_defined_drugs": True,
            "old_tables_not_overwritten": True, "sample_data_read": False,
            "model_artifacts_loaded": False, "new_fits": 0, "new_resampling_draws": 0,
            "new_confidence_intervals": 0, "new_inferential_families": 0,
        },
    }
    output = directory / "revision_table_manifest.json"
    output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "manifest": output.relative_to(ROOT).as_posix(), "S9_rows": len(s9), "S10_rows": len(s10)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
