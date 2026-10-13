"""Public-only snapshot checks using synthetic aggregate fixtures.

No private OOF, checkpoint, response or expression file is accessed. Fixtures
exercise the complete frozen aggregate schema without generating model fits.
"""
from __future__ import annotations
import csv
import importlib.util
import itertools
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("verify_public_ordinal", ROOT / "scripts/verify_public_ordinal.py")
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


def make_snapshot(output):
    config = json.loads((ROOT / "config/ordinal_diagnostic_v1.json").read_text(encoding="utf-8"))
    lock = json.loads((ROOT / "config/ordinal_diagnostic_v1.lock.json").read_text(encoding="utf-8"))
    drugs = list(map(str, config["drug_ids"]))
    scopes, permutations, models = drugs + [VERIFY.MACRO], list(map(str, range(10))), config["models"]
    folds, metrics = ["0", "1", "2"], VERIFY.METRICS
    frames = {name: [] for name in VERIFY.FILES}

    def interval(point):
        return dict(point=point, bootstrap_requested=10000, bootstrap_valid=10000,
            bootstrap_invalid=0, bootstrap_valid_fraction=1., conditional_95_lower=point,
            conditional_95_upper=point, family_adjusted_conditional_lower="",
            family_adjusted_conditional_upper="", family_size="", ci_status="reported",
            interval_scope="conditional_on_observed_OOF_records_not_full_train_refit_or_guaranteed_algorithmic_CV_coverage")

    values = {}
    for drug, permutation, model in itertools.product(scopes, permutations, models):
        drug_position = drugs.index(drug) if drug in drugs else 1.5
        values[(drug, permutation, model)] = dict(
            spearman=.2 + .05 * drug_position + .01 * int(permutation),
            rmse=.1 + .002 * drug_position + .001 * int(permutation),
            mae=.08 + .001 * drug_position, r2=.1 + .02 * drug_position)
        frames["point_metrics.csv"].append(dict(gdsc_drug_id=drug, permutation_id=permutation,
            model=model, n_bootstrap_global_groups=626,
            n_unique_biological_groups=620 if drug in drugs else 624,
            prediction_rms_difference_from_tissue_null=0., **values[(drug, permutation, model)]))
    for drug, model, metric in itertools.product(scopes, models, metrics):
        frames["model_intervals.csv"].append(dict(gdsc_drug_id=drug, model=model,
            metric=metric, **interval(values[(drug, "0", model)][metric])))
    for drug in scopes:
        for name, candidate, reference in VERIFY.contrasts(config):
            for metric in metrics:
                base = dict(gdsc_drug_id=drug, contrast=name, candidate=candidate,
                    reference=reference, metric=metric, positive_favors_candidate=True)
                frames["paired_effect_intervals.csv"].append(dict(**base, **interval(0.)))
                for permutation in permutations[1:]:
                    frames["null_contrasts.csv"].append(dict(**base, permutation_id=permutation, point=0.))
                frames["null_contrast_summary.csv"].append(dict(**base,
                    null_runs_requested=9, null_runs_recorded=9, null_runs_valid=9,
                    null_min=0., null_median=0., null_max=0.))
    for drug, permutation, model in itertools.product(drugs, permutations, models):
        base = dict(gdsc_drug_id=drug, permutation_id=permutation, model=model)
        tissue = model == "tissue_null"
        frames["molecular_and_failure_summary.csv"].append(dict(**base, outer_fits=3,
            zero_molecular_fits=0 if tissue else 3, nonzero_molecular_fits=0,
            active_feature_total=0, active_feature_min="" if tissue else 0,
            active_feature_max="" if tissue else 0, warning_outer_fits=0, fallback_outer_fits=0))
        frames["all_fit_events_summary.csv"].append(dict(**base,
            actual_fit_events=3 if tissue else 93, all_fit_warning_events=0, all_fit_fallback_events=0))
        frames["resource_summary.csv"].append(dict(**base, outer_tasks=3,
            task_wall_seconds_sum=6., task_cpu_seconds_sum=3., outer_fit_seconds_sum=1.5,
            outer_prediction_seconds_sum=.3, sampled_peak_rss_bytes_max=1000,
            windows_lifetime_peak_working_set_bytes_max=1200))
        for fold in folds:
            frames["chosen_parameters.csv"].append(dict(**base, outer_fold=fold,
                active_features=0, selected_genes="" if tissue else 200,
                candidate_features="" if tissue else 19900 if model.endswith("pair_full") else 200,
                chosen_params=json.dumps({} if tissue else dict(alpha_ratio=1., l1_ratio=.5)),
                warnings="[]", fallback_reason=""))
            frames["resource_tasks.csv"].append(dict(**base, outer_fold=fold,
                task_wall_seconds=2., task_cpu_seconds=1., outer_fit_seconds=.5,
                outer_prediction_seconds=.1, sampled_peak_rss_bytes=1000,
                windows_lifetime_peak_working_set_bytes=1200))
        if not tissue:
            for fold_a, fold_b in itertools.combinations(folds, 2):
                frames["active_feature_stability.csv"].append(dict(**base,
                    fold_a=fold_a, fold_b=fold_b, candidate_gene_jaccard=1.,
                    active_feature_jaccard="", same_sign_active_jaccard="", empty_active_union=True))
    for drug, permutation in itertools.product(drugs, permutations):
        observed = permutation == "0"
        frames["permutation_coverage.csv"].append(dict(drug_id=drug, permutation_id=permutation,
            n_groups=620, singleton_groups=2, identity_donors=620 if observed else 2,
            unchanged_response_count=620 if observed else 2,
            changed_response_fraction=0. if observed else 618 / 620))
    for name, rows in frames.items():
        with (output / name).open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    provenance = dict(status="completed",
        config_sha256=VERIFY.sha(ROOT / "config/ordinal_diagnostic_v1.json"),
        lock_sha256=VERIFY.sha(ROOT / "config/ordinal_diagnostic_v1.lock.json"),
        analyzer_sha256=VERIFY.sha(ROOT / "src/analyze_ordinal_diagnostic.py"),
        source_run="runs/" + config["run_name"], external_response_accessed=False,
        blocked_source_access_attempts=[], no_permutation_p_values=True, new_model_fits=0,
        no_repeat_N_inflation=True, all_null_runs_retained=True,
        bootstrap_shared_across_models_drugs=True, confidence_scope=config["statistics"]["confidence_scope"],
        expected_outer_tasks=lock["contract"]["tasks"], actual_fit_events=lock["contract"]["fits"],
        oof_rows=lock["contract"]["oof_predictions"], n_global_biological_groups=626,
        bootstrap_draws=10000, bootstrap_seed=config["statistics"]["bootstrap_seed"],
        output_sha256={name: VERIFY.sha(output / name) for name in frames})
    (output / "analysis_provenance.json").write_text(json.dumps(provenance), encoding="utf-8")
    return provenance


class PublicOrdinalSnapshotTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name)

    def test_complete_synthetic_snapshot_passes(self):
        make_snapshot(self.output)
        result = VERIFY.verify(self.output)
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["files_verified"], 12)
        self.assertEqual(result["successful_pipeline_fits"], 22440)
        self.assertEqual(result["observed_paired_intervals"], 300)
        self.assertEqual(result["null_contrast_points"], 2700)
        self.assertFalse(result["private_data_opened"])

    def test_missing_snapshot_is_failure_not_skip(self):
        with self.assertRaisesRegex(ValueError, "provenance is missing"):
            VERIFY.verify(self.output)

    def test_byte_tamper_is_rejected(self):
        make_snapshot(self.output)
        path = self.output / "resource_summary.csv"
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            VERIFY.verify(self.output)

    def test_rehashed_fit_budget_tamper_is_still_rejected(self):
        provenance = make_snapshot(self.output)
        path = self.output / "all_fit_events_summary.csv"
        raw = path.read_bytes()
        self.assertIn(b",93,", raw)
        path.write_bytes(raw.replace(b",93,", b",94,", 1))
        provenance["output_sha256"][path.name] = VERIFY.sha(path)
        (self.output / "analysis_provenance.json").write_text(json.dumps(provenance), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "fit count differs"):
            VERIFY.verify(self.output)


if __name__ == "__main__":
    unittest.main()
