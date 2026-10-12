"""Synthetic integrity contracts for resumable formal CV and paired inference.

No original response/expression files are opened. Bootstrap fixtures are
checked against explicit duplicated-row calculations with exact ties.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from run_formal_evaluation import (atomic_json, checkpoint_name, find_pending_tasks,
                                   validate_resume_checkpoint, make_global_split_records,
                                   local_split_indices, task_list)
from analyze_formal_results import (generate_global_bootstrap_weights,
                                    weighted_midrank, weighted_metrics, interval_summary, analyze_oof)


class ResumeCheckpointIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.task = {"drug_id": 1931, "repeat": 0, "outer_fold": 0, "model": "mean"}
        self.freeze = "a" * 64
        rows = [{**self.task, "drug_name": "Synthetic compound", "canonical_model_id": f"SYN-{i}",
                 "biological_group_id": f"GROUP-{i}", "tissue": "Synthetic", "observed_auc": y,
                 "predicted_auc": 0.8} for i, y in ((1, 0.7), (2, 0.9))]
        self.valid = {"schema_version": "1.0", "task": self.task, "freeze_sha256": self.freeze,
                      "status": "completed", "oof_rows": rows, "expected_test_ids": ["SYN-1", "SYN-2"],
                      "expected_training_expression_ids": [10, 11, 12], "inner_training_expression_ids": [],
                      "outer_test_expression_ids": [1, 2],
                      "result": {"chosen_params": {}, "predictions": [0.8, 0.8], "metrics": {"n": 2}},
                      "fit_events": [], "resources": {}, "error": None}

    def write(self, folder, payload=None):
        path = Path(folder) / checkpoint_name(self.task)
        atomic_json(path, self.valid if payload is None else payload)
        return path

    def test_completed_matching_checkpoint_is_not_scheduled_again(self):
        with tempfile.TemporaryDirectory(prefix="reo_resume_fixture_") as folder:
            path = self.write(folder)
            before = path.read_bytes()
            self.assertTrue(validate_resume_checkpoint(path, self.task, self.freeze))
            self.assertEqual(find_pending_tasks([self.task], folder, self.freeze), [])
            self.assertEqual(path.read_bytes(), before)

    def test_failed_checkpoint_remains_pending_and_is_not_silently_deleted(self):
        failed = copy.deepcopy(self.valid)
        failed.update(status="failed", error="Synthetic failure", oof_rows=[])
        with tempfile.TemporaryDirectory(prefix="reo_failed_resume_fixture_") as folder:
            path = self.write(folder, failed)
            before = path.read_bytes()
            self.assertFalse(validate_resume_checkpoint(path, self.task, self.freeze))
            self.assertEqual(find_pending_tasks([self.task], folder, self.freeze), [self.task])
            self.assertEqual(path.read_bytes(), before)

    def test_changed_task_or_freeze_is_a_fatal_integrity_error(self):
        with tempfile.TemporaryDirectory(prefix="reo_wrong_resume_fixture_") as folder:
            path = self.write(folder)
            for expected, freeze in (({**self.task, "repeat": 1}, self.freeze), (self.task, "b" * 64)):
                with self.subTest(expected=expected, freeze=freeze):
                    with self.assertRaises(ValueError):
                        validate_resume_checkpoint(path, expected, freeze)

    def test_duplicate_missing_or_nonfinite_oof_cannot_be_reused(self):
        examples = []
        duplicate = copy.deepcopy(self.valid)
        duplicate["oof_rows"].append(copy.deepcopy(duplicate["oof_rows"][0]))
        examples.append(duplicate)
        missing = copy.deepcopy(self.valid)
        missing["oof_rows"].pop()
        examples.append(missing)
        nonfinite = copy.deepcopy(self.valid)
        nonfinite["oof_rows"][0]["predicted_auc"] = None
        examples.append(nonfinite)
        wrong_context = copy.deepcopy(self.valid)
        wrong_context["oof_rows"][0]["model"] = "reo"
        examples.append(wrong_context)
        for number, payload in enumerate(examples):
            with self.subTest(example=number):
                with tempfile.TemporaryDirectory(prefix="reo_bad_resume_fixture_") as folder:
                    path = self.write(folder, payload)
                    with self.assertRaises((ValueError, TypeError)):
                        validate_resume_checkpoint(path, self.task, self.freeze)

    def test_atomic_writer_leaves_no_tmp_or_partial_checkpoint(self):
        with tempfile.TemporaryDirectory(prefix="reo_atomic_fixture_") as folder:
            path = self.write(folder)
            saved = json.loads(path.read_text(encoding="utf-8"))
            digest = saved.pop("payload_sha256")
            expected = hashlib.sha256(json.dumps(saved, sort_keys=True, ensure_ascii=False, allow_nan=False,
                                                  separators=(",", ":")).encode("utf-8")).hexdigest()
            self.assertEqual(digest, expected)
            self.assertEqual(saved, self.valid)
            self.assertEqual(list(Path(folder).glob("*.tmp")), [])

    def test_finite_prediction_tampering_is_rejected_before_resume(self):
        with tempfile.TemporaryDirectory(prefix="reo_checkpoint_tamper_fixture_") as folder:
            path = self.write(folder)
            altered = json.loads(path.read_text(encoding="utf-8"))
            altered["oof_rows"][0]["predicted_auc"] = 0.123
            path.write_text(json.dumps(altered), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "checksum"):
                validate_resume_checkpoint(path, self.task, self.freeze)


class FormalGroupedRepeatIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.ids = pd.DataFrame({"canonical_model_id": [f"C{i}" for i in range(24)],
                                 "biological_group_id": np.repeat([f"G{i}" for i in range(12)], 2),
                                 "tissue": np.repeat([f"T{i % 3}" for i in range(12)], 2)})
        self.config = {"drug_ids": [1931, 1036], "models": ["mean", "reo"],
                       "outer_repeats": 3, "outer_folds": 3, "inner_folds": 2,
                       "repeat_seeds": [301, 302, 303]}

    def test_global_groups_are_isolated_in_every_repeat_and_inner_partition(self):
        records, _ = make_global_split_records(self.ids, self.config)
        for repeat in range(3):
            test_coverage = []
            for outer in range(3):
                block = records[(records["repeat"] == repeat) & (records.outer_fold == outer)]
                train = block[(block.inner_fold == -1) & block.role.eq("train")]
                test = block[(block.inner_fold == -1) & block.role.eq("test")]
                self.assertFalse(set(train.biological_group_id) & set(test.biological_group_id))
                self.assertEqual(set(train.canonical_model_id) | set(test.canonical_model_id), set(self.ids.canonical_model_id))
                test_coverage.extend(test.canonical_model_id.tolist())
                for inner in range(2):
                    it = block[(block.inner_fold == inner) & block.role.eq("train")]
                    iv = block[(block.inner_fold == inner) & block.role.eq("validation")]
                    self.assertFalse(set(it.biological_group_id) & set(iv.biological_group_id))
                    self.assertEqual(set(it.canonical_model_id) | set(iv.canonical_model_id), set(train.canonical_model_id))
            self.assertEqual(sorted(test_coverage), sorted(self.ids.canonical_model_id))

    def test_drug_intersection_preserves_global_roles_and_relative_inner_indices(self):
        records, _ = make_global_split_records(self.ids, self.config)
        data = self.ids.iloc[np.arange(24) % 5 != 0].iloc[::-1].reset_index(drop=True)
        task = {"drug_id": 1931, "repeat": 1, "outer_fold": 2, "model": "reo"}
        train, test, inner = local_split_indices(data, records, task, 2)
        block = records[(records["repeat"] == 1) & (records.outer_fold == 2)]
        for role, positions in (("train", train), ("test", test)):
            expected = set(block[(block.inner_fold == -1) & block.role.eq(role)].canonical_model_id) & set(data.canonical_model_id)
            self.assertEqual(set(data.iloc[positions].canonical_model_id), expected)
        for number, (it, iv) in enumerate(inner):
            outer = data.iloc[train]
            for role, positions in (("train", it), ("validation", iv)):
                expected = set(block[(block.inner_fold == number) & block.role.eq(role)].canonical_model_id) & set(outer.canonical_model_id)
                self.assertEqual(set(outer.iloc[positions].canonical_model_id), expected)

    def test_formal_tasks_are_unique_cartesian_products(self):
        jobs = task_list(self.config)
        keys = [(j["drug_id"], j["repeat"], j["outer_fold"], j["model"]) for j in jobs]
        self.assertEqual(len(keys), 2 * 3 * 3 * 2)
        self.assertEqual(len(set(keys)), len(keys))


class SharedBiologicalBootstrapIntegrityTests(unittest.TestCase):
    def test_bootstrap_draws_count_groups_not_repeats_or_duplicate_models(self):
        first_order, first = generate_global_bootstrap_weights(["G3", "G1", "G2", "G1"] * 3, 17, 1204)
        second_order, second = generate_global_bootstrap_weights(["G1", "G2", "G3"], 17, 1204)
        np.testing.assert_array_equal(first_order, second_order)
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(first.sum(axis=1), np.full(17, 3))

    def test_weighted_tied_ranks_and_metrics_match_explicit_duplicate_rows(self):
        y = np.array([0.7, 0.7, 0.8, 0.95, 0.95])
        pred = np.array([0.8, 0.75, 0.75, 0.96, 0.96])
        counts = np.array([[2, 0, 1, 3, 0], [0, 2, 0, 1, 2], [1, 1, 1, 1, 1], [0, 0, 0, 0, 0]])
        ranks = weighted_midrank(y, counts)
        actual = weighted_metrics(y, pred, counts)
        for draw, multiplicities in enumerate(counts):
            indices = np.repeat(np.arange(len(y)), multiplicities)
            if not len(indices):
                self.assertTrue(all(np.isnan(actual[name][draw]) for name in actual))
                continue
            yy, pp = y[indices], pred[indices]
            explicit_ranks = rankdata(yy, method="average")
            for index in np.flatnonzero(multiplicities):
                np.testing.assert_array_equal(ranks[draw, index], explicit_ranks[np.flatnonzero(indices == index)[0]])
            self.assertAlmostEqual(actual["spearman"][draw], spearmanr(yy, pp).statistic, places=12)
            self.assertAlmostEqual(actual["rmse"][draw], np.sqrt(np.mean((yy - pp) ** 2)), places=12)
            self.assertAlmostEqual(actual["mae"][draw], np.mean(np.abs(yy - pp)), places=12)
            self.assertAlmostEqual(actual["r2"][draw], 1 - np.sum((yy - pp) ** 2) / np.sum((yy - yy.mean()) ** 2), places=12)

    def test_invalid_bootstrap_fraction_withholds_ci_without_discarding_task(self):
        config = {"bootstrap": {"minimum_valid_draw_fraction": 0.99, "unadjusted_level": 0.95},
                  "multiplicity": {"family_alpha": 0.05}}
        samples = np.arange(100, dtype=float)
        samples[:2] = np.nan
        result = interval_summary(1., samples, config, family_size=9)
        self.assertEqual(result["bootstrap_requested"], 100)
        self.assertEqual(result["bootstrap_valid"], 98)
        self.assertEqual(result["bootstrap_invalid"], 2)
        self.assertEqual(result["ci_status"], "withheld_insufficient_valid_draws")
        self.assertTrue(np.isnan(result["conditional_95_lower"]))

    def test_primary_family_nine_percentile_boundaries_are_prespecified(self):
        config = {"bootstrap": {"minimum_valid_draw_fraction": 0.99, "unadjusted_level": 0.95},
                  "multiplicity": {"family_alpha": 0.05}}
        samples = np.linspace(-0.2, 0.8, 10000)
        result = interval_summary(0.3, samples, config, family_size=9)
        expected = np.quantile(samples, [0.05 / 18., 1 - 0.05 / 18.], method="linear")
        np.testing.assert_allclose([result["family_adjusted_conditional_lower"], result["family_adjusted_conditional_upper"]], expected, atol=1e-15, rtol=0)


class FormalEstimandAndPairingTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((ROOT / "config/formal_statistics_v1.json").read_text(encoding="utf-8"))
        self.config["drug_ids_in_frozen_rank"] = [1, 2]
        self.config["bootstrap"].update(draws=4, computational_block_draws=2)
        self.weights = np.array([[2, 1, 1, 1, 1, 0], [0, 1, 2, 1, 1, 1],
                                 [1, 1, 1, 1, 1, 1], [1, 2, 0, 1, 2, 0]], dtype=int)
        self.manifest = pd.DataFrame([{"gdsc_drug_id": drug, "canonical_model_id": f"C{i}",
                                       "biological_group_id": f"G{i}"}
                                      for drug, n in ((1, 6), (2, 5)) for i in range(n)])
        truth = np.arange(1, 7) / 10.
        reference = np.array([0.1, 0.3, 0.2, 0.4, 0.6, 0.5])
        candidates = [truth, truth[::-1], np.array([0.2, 0.1, 0.3, 0.4, 0.5, 0.6])]
        rows = []
        for drug, n in ((1, 6), (2, 5)):
            for model in self.config["models"]:
                for repeat in range(3):
                    prediction = candidates[repeat] if model in ("reo", "reo_rf") else reference
                    if model == "mean":
                        prediction = np.full(6, 0.35)
                    for i in range(n):
                        rows.append({"drug_id": drug, "canonical_model_id": f"C{i}", "biological_group_id": f"G{i}",
                                     "model": model, "repeat": repeat, "outer_fold": i % 3,
                                     "observed_auc": truth[i], "predicted_auc": prediction[i]})
        self.oof = pd.DataFrame(rows)

    def test_repeat_metric_mean_is_not_metric_of_averaged_predictions_or_fake_n(self):
        result = analyze_oof(self.oof, self.manifest, self.config, weights=self.weights)
        data = self.oof[self.oof.drug_id.eq(1) & self.oof.model.eq("reo")]
        expected = np.mean([spearmanr(part.observed_auc, part.predicted_auc).statistic for _, part in data.groupby("repeat")])
        averaged = data.groupby("canonical_model_id")[['observed_auc', 'predicted_auc']].mean()
        wrong = spearmanr(averaged.observed_auc, averaged.predicted_auc).statistic
        row = result["model_intervals"].query("gdsc_drug_id == 1 and model == 'reo' and metric == 'spearman'").iloc[0]
        self.assertAlmostEqual(row.point, expected, places=12)
        self.assertGreater(abs(row.point - wrong), 0.1)
        self.assertEqual(row.n_unique_biological_groups, 6)
        np.testing.assert_array_equal(result["bootstrap_counts"], self.weights)

    def test_paired_macro_ci_uses_one_identical_draw_across_models_drugs_and_repeats(self):
        result = analyze_oof(self.oof, self.manifest, self.config, weights=self.weights)
        samples = []
        for draw in self.weights:
            differences = []
            for drug, n in ((1, 6), (2, 5)):
                indices = np.repeat(np.arange(n), draw[:n])
                repeat_differences = []
                for repeat in range(3):
                    subset = self.oof[self.oof.drug_id.eq(drug) & self.oof["repeat"].eq(repeat)]
                    a = subset[subset.model.eq("reo")].sort_values("canonical_model_id")
                    r = subset[subset.model.eq("expression")].sort_values("canonical_model_id")
                    repeat_differences.append(spearmanr(a.observed_auc.to_numpy()[indices], a.predicted_auc.to_numpy()[indices]).statistic
                                              - spearmanr(r.observed_auc.to_numpy()[indices], r.predicted_auc.to_numpy()[indices]).statistic)
                differences.append(np.mean(repeat_differences))
            samples.append(np.mean(differences))
        row = result["paired_effect_intervals"].query("gdsc_drug_id == 'macro_equal_fixed_drugs' and comparison == 'primary_EN' and metric == 'spearman'").iloc[0]
        expected = np.quantile(samples, [0.05 / 18, 1 - 0.05 / 18], method="linear")
        np.testing.assert_allclose([row.family_adjusted_conditional_lower, row.family_adjusted_conditional_upper], expected, atol=1e-12, rtol=0)
        self.assertEqual(row.n_unique_biological_groups, 6)
        self.assertEqual(row.fixed_drug_count, 2)

    def test_incomplete_repeat_does_not_silently_nanmean_or_remove_drug_from_macro(self):
        missing = self.oof[~(self.oof.drug_id.eq(1) & self.oof.model.eq("reo") & self.oof["repeat"].eq(2))]
        result = analyze_oof(missing, self.manifest, self.config, weights=self.weights)
        rows = result["paired_effect_intervals"].query("comparison == 'primary_EN' and metric == 'spearman'")
        drug = rows[rows.gdsc_drug_id.eq(1)].iloc[0]
        macro = rows[rows.gdsc_drug_id.eq("macro_equal_fixed_drugs")].iloc[0]
        self.assertTrue(np.isnan(drug.point) and np.isnan(macro.point))
        self.assertEqual(drug.ci_status, "withheld_undefined_point")
        self.assertEqual(macro.fixed_drug_count, 2)
        self.assertIn("undefined_constituent_drugs:1", macro.undefined_reason)


if __name__ == "__main__":
    unittest.main()
