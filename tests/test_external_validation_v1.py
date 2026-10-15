"""Synthetic-only tests: no external response, metadata, or source file access."""
from pathlib import Path
import json
import sys
import unittest

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from external_validation_v1 import (
    SharedBootstrap, bootstrap_paired_rho, make_shared_bootstrap,
    names_sha256, paired_common_table, paired_rho_difference,
    predict_external_arm, prepare_gateway, spearman_checked,
    shared_bootstrap_sha256, strict_fixed_eight_macro, validate_model_contract,
)


GENES = ("ENSG000001", "ENSG000002", "ENSG000003")
DRUGS = (11, 12, 13, 14, 15, 16, 17, 18)


class StubBuilder:
    def __init__(self, representation="expression"):
        self.stable_gene_ids = GENES
        self.gene_ids_ = GENES[:2]
        self.gene_indices_ = np.asarray([0, 1])
        if representation == "expression":
            self.feature_names_ = np.asarray(["log2TPM:" + name for name in self.gene_ids_])
        elif representation == "pair_full":
            self.feature_names_ = np.asarray(["REO:ENSG000001>ENSG000002"])
            self.pair_indices_ = np.asarray([[0, 1]])


class StubModel:
    def __init__(self, representation="expression", *, zero=False):
        self.stable_gene_ids = GENES
        self.n_features_in_ = len(GENES) + 2
        self.model_kind = representation
        self.training_tissue_codes_ = [0, 1]
        self.builder_ = None if representation == "tissue_null" else StubBuilder(representation)
        self.fitted_kind_ = representation
        self.n_nonzero_molecular_coeffs_ = 0 if zero else 1
        self.received = []
        self.zero = zero

    def predict(self, x):
        self.received.append(x.copy())
        tissue = np.where(x[:, -2] == 0, .2, np.where(x[:, -2] == 1, .4, .3))
        return tissue if self.model_kind == "tissue_null" or self.zero else tissue + .01 * x[:, 0]


class GatewayTests(unittest.TestCase):
    def _batch(self):
        identity = pd.DataFrame({"canonical_model_id": ["A", "B", "C", "D", "E", "F"],
                                 "biological_group_id": ["G1", "G2", "G3", "G4", "G5", "G6"],
                                 "tissue": ["Novel", "Known", "Known", "Known", "Known", "Known"]})
        tpm = pd.DataFrame([[1, 2, 3], [1, np.nan, 2], [-1, 2, 3],
                            [1, np.inf, 3], [1, "bad", 3]],
                           index=["A", "B", "C", "D", "E"], columns=GENES)
        return prepare_gateway(identity, tpm, stable_gene_ids=GENES,
                               tissue_to_code={"Known": 0, "Other": 1},
                               training_tissue_codes=[0, 1],
                               expected_gene_sha256=names_sha256(GENES))

    def test_abstain_unknown_tissue_and_zero_solution(self):
        batch = self._batch()
        self.assertEqual(batch.unknown_tissue_code, 2)
        self.assertEqual(batch.molecular_valid.tolist(), [True, False, False, False, False, False])
        self.assertEqual(batch.abstain_reason, ("", "missing_tpm", "negative_tpm", "nonfinite_tpm",
                                                  "nonnumeric_tpm", "missing_expression_row"))
        molecular = StubModel(zero=True)
        result = predict_external_arm(molecular, batch, drug_id=11, arm="mad_expression",
                                      representation="expression")
        self.assertEqual(result.status.tolist(), ["predicted"] + ["abstained"] * 5)
        self.assertTrue(result.zero_molecular_solution.all())
        self.assertEqual(molecular.received[0].shape, (1, len(GENES) + 2))
        self.assertEqual(molecular.received[0][0, -2], 2)  # fresh unknown category
        self.assertAlmostEqual(result.iloc[0].predicted_auc, .3)  # global tissue term
        null = StubModel("tissue_null")
        null_result = predict_external_arm(null, batch, drug_id=11, arm="tissue_null",
                                           representation="tissue_null")
        self.assertEqual(null_result.status.tolist(), ["predicted"] * 6)
        self.assertTrue(np.all(null.received[0][:, :len(GENES)] == 0))
        self.assertTrue(np.isnan(result.iloc[1:].predicted_auc).all())

    def test_exact_gene_feature_and_unknown_code_contracts(self):
        batch = self._batch()
        model = StubModel("pair_full")
        validate_model_contract(model, representation="pair_full", stable_gene_ids=GENES,
                                training_tissue_codes=[0, 1],
                                expected_selected_gene_ids=GENES[:2],
                                expected_feature_names=["REO:ENSG000001>ENSG000002"])
        with self.assertRaisesRegex(ValueError, "gene order"):
            validate_model_contract(model, representation="pair_full",
                                    stable_gene_ids=GENES[::-1], training_tissue_codes=[0, 1])
        with self.assertRaisesRegex(ValueError, "feature order|Feature order"):
            validate_model_contract(model, representation="pair_full", stable_gene_ids=GENES,
                                    training_tissue_codes=[0, 1],
                                    expected_feature_names=["REO:ENSG000002>ENSG000001"])
        model.builder_.feature_names_ = np.asarray(["REO:ENSG000002>ENSG000001"])
        with self.assertRaisesRegex(ValueError, "direction"):
            predict_external_arm(model, batch, drug_id=11, arm="mad_pair_full",
                                 representation="pair_full")
        model.builder_.feature_names_ = np.asarray(["REO:ENSG000001>ENSG000002"])
        model.builder_.pair_indices_ = np.asarray([[1, 0]])
        with self.assertRaisesRegex(ValueError, "Pair input indices"):
            validate_model_contract(model, representation="pair_full", stable_gene_ids=GENES,
                                    training_tissue_codes=[0, 1])
        identity = batch.identity
        reversed_columns = pd.DataFrame([[1, 2, 3]], index=["A"], columns=GENES[::-1])
        with self.assertRaisesRegex(ValueError, "exactly match"):
            prepare_gateway(identity, reversed_columns, stable_gene_ids=GENES,
                            tissue_to_code={"Known": 0, "Other": 1}, training_tissue_codes=[0, 1])
        with self.assertRaisesRegex(ValueError, "collides"):
            prepare_gateway(identity, pd.DataFrame([[1, 2, 3]], index=["A"], columns=GENES),
                            stable_gene_ids=GENES, tissue_to_code={"Known": 0, "Other": 1},
                            training_tissue_codes=[0, 1], unknown_tissue_code=1)
        complete = pd.DataFrame([[1, 2, 3]], index=["A"], columns=GENES)
        with self.assertRaisesRegex(ValueError, "exactly match"):
            prepare_gateway(identity, complete.assign(EXTRA=0), stable_gene_ids=GENES,
                            tissue_to_code={"Known": 0, "Other": 1}, training_tissue_codes=[0, 1])
        with self.assertRaisesRegex(ValueError, "exactly match"):
            prepare_gateway(identity, complete.drop(columns=GENES[-1]), stable_gene_ids=GENES,
                            tissue_to_code={"Known": 0, "Other": 1}, training_tissue_codes=[0, 1])

    def test_nonfinite_model_output_and_duplicate_group_are_hard_failures(self):
        batch = self._batch()
        model = StubModel()
        model.predict = lambda x: np.full(len(x), np.nan)
        with self.assertRaisesRegex(ValueError, "Nonfinite fitted-model prediction"):
            predict_external_arm(model, batch, drug_id=11, arm="mad_expression",
                                 representation="expression")
        identity = batch.identity.copy()
        identity.loc[1, "biological_group_id"] = identity.loc[0, "biological_group_id"]
        with self.assertRaisesRegex(ValueError, "Duplicate biological group"):
            prepare_gateway(identity, pd.DataFrame([[1, 2, 3]], index=["A"], columns=GENES),
                            stable_gene_ids=GENES, tissue_to_code={"Known": 0, "Other": 1},
                            training_tissue_codes=[0, 1])
        with self.assertRaisesRegex(ValueError, "overlaps frozen training"):
            prepare_gateway(batch.identity, pd.DataFrame([[1, 2, 3]], index=["A"], columns=GENES),
                            stable_gene_ids=GENES, tissue_to_code={"Known": 0, "Other": 1},
                            training_tissue_codes=[0, 1], forbidden_training_group_ids=["G4"])
        identity = batch.identity.copy()
        identity.loc[1, "canonical_model_id"] = identity.loc[0, "canonical_model_id"]
        with self.assertRaisesRegex(ValueError, "duplicate IDs"):
            prepare_gateway(identity, pd.DataFrame([[1, 2, 3]], index=["A"], columns=GENES),
                            stable_gene_ids=GENES, tissue_to_code={"Known": 0, "Other": 1},
                            training_tissue_codes=[0, 1])

    def test_missing_tissue_abstains_even_for_null(self):
        identity = pd.DataFrame({"canonical_model_id": ["A", "B"],
                                 "biological_group_id": ["G1", "G2"],
                                 "tissue": ["Known", ""]})
        tpm = pd.DataFrame([[1., 2., 3.], [1., 2., 3.]], index=["A", "B"], columns=GENES)
        batch = prepare_gateway(identity, tpm, stable_gene_ids=GENES,
                                tissue_to_code={"Known": 0, "Other": 1},
                                training_tissue_codes=[0, 1])
        null = predict_external_arm(StubModel("tissue_null"), batch, drug_id=11,
                                    arm="tissue_null", representation="tissue_null")
        self.assertEqual(null.status.tolist(), ["predicted", "abstained"])
        self.assertEqual(null.iloc[1].abstain_reason, "invalid_tissue")

    def test_global_cohort_tissue_absent_for_this_drug_is_unknown(self):
        # The global development codebook could contain B=2, yet this drug was
        # trained only on A=0 and C=5. B must not collide with a trained code.
        identity = pd.DataFrame({"canonical_model_id": ["A", "B"],
                                 "biological_group_id": ["G1", "G2"],
                                 "tissue": ["GloballyKnownButDrugUnseen", "DrugSeenC"]})
        tpm = pd.DataFrame([[1., 2., 3.], [1., 2., 3.]], index=["A", "B"], columns=GENES)
        batch = prepare_gateway(identity, tpm, stable_gene_ids=GENES,
                                tissue_to_code={"DrugSeenA": 0, "DrugSeenC": 5},
                                training_tissue_codes=[0, 5])
        self.assertEqual(batch.tissue_codes.tolist(), [6., 5.])
        self.assertEqual(batch.unknown_tissue.tolist(), [True, False])
        model = StubModel("tissue_null")
        model.training_tissue_codes_ = [0, 5]
        result = predict_external_arm(model, batch, drug_id=11, arm="tissue_null",
                                      representation="tissue_null")
        self.assertEqual(result.status.tolist(), ["predicted", "predicted"])
        self.assertEqual(result.unknown_training_tissue.tolist(), [True, False])

    def test_real_frozen_regressor_interface_with_synthetic_tpm(self):
        from ordinal_diagnostic_models import DiagnosticRegressor
        config = json.loads((Path(__file__).resolve().parents[1] /
                             "config/ordinal_repeated_v1.json").read_text(encoding="utf-8"))
        config["feature"]["n_genes"] = 2
        config["feature"]["diagnostic_genes"] = 2
        raw = np.asarray([[1., 4., 3.], [2., 5., 3.], [3., 1., 2.], [4., 2., 2.],
                          [5., 3., 4.], [6., 4., 4.]], dtype=np.float64)
        tissue = np.asarray([0., 0., 0., 5., 5., 5.])
        fit_x = np.column_stack((raw, tissue, np.arange(6, dtype=float)))
        y = np.asarray([.1, .2, .3, .5, .6, .7])
        model = DiagnosticRegressor(model_kind="pair_full", selection_rule="mad",
                                    config=config, stable_gene_ids=GENES,
                                    alpha_ratio=1.0, l1_ratio=.5).fit(fit_x, y)
        identity = pd.DataFrame({"canonical_model_id": ["X", "Y"],
                                 "biological_group_id": ["GX", "GY"],
                                 "tissue": ["Unseen", "Seen"]})
        tpm = pd.DataFrame([[2., 4., 3.], [np.nan, 4., 3.]], index=["X", "Y"], columns=GENES)
        batch = prepare_gateway(identity, tpm, stable_gene_ids=GENES,
                                tissue_to_code={"Seen": 0, "OtherSeen": 5},
                                training_tissue_codes=model.training_tissue_codes_)
        predicted = predict_external_arm(model, batch, drug_id=11, arm="mad_pair_full",
                                         representation="pair_full")
        self.assertEqual(predicted.status.tolist(), ["predicted", "abstained"])
        direct = model.predict(np.asarray([[2., 4., 3., 6., 0.]], dtype=np.float64))
        self.assertAlmostEqual(predicted.iloc[0].predicted_auc, direct[0])


class StatisticsTests(unittest.TestCase):
    def test_constant_and_small_group_rho_are_explicit_na(self):
        self.assertEqual(spearman_checked([1, 2], [2, 1], ["A", "B"]),
                         (None, "insufficient_independent_groups"))
        self.assertEqual(spearman_checked([1, 1, 1], [1, 2, 3], ["A", "B", "C"]),
                         (None, "constant_response"))
        self.assertEqual(spearman_checked([1, 2, 3], [1, 1, 1], ["A", "B", "C"]),
                         (None, "constant_prediction"))
        self.assertAlmostEqual(spearman_checked([1, 2, 3], [3, 2, 1], ["A", "B", "C"])[0], -1.)

    def test_common_id_coverage_and_no_available_case_macro(self):
        response = pd.DataFrame({"drug_id": 11, "canonical_model_id": ["A", "B", "C", "D"],
                                 "biological_group_id": ["G1", "G2", "G3", "G4"],
                                 "observed_auc": [1., 2., 3., 4.]})
        prediction = pd.DataFrame([
            (11, "pair", "A", "G1", .1, "predicted"),
            (11, "pair", "B", "G2", .2, "predicted"),
            (11, "pair", "C", "G3", np.nan, "abstained"),
            (11, "pair", "D", "G4", .4, "predicted"),
            (11, "expression", "A", "G1", .1, "predicted"),
            (11, "expression", "B", "G2", .3, "predicted"),
            (11, "expression", "C", "G3", .2, "predicted"),
            (11, "expression", "D", "G4", .4, "predicted"),
        ], columns=["drug_id", "model", "canonical_model_id", "biological_group_id", "predicted_auc", "status"])
        pair, coverage = paired_common_table(response, prediction, drug_id=11,
                                             arm_a="pair", arm_b="expression", design_min_groups=80)
        self.assertEqual(pair.canonical_model_id.tolist(), ["A", "B", "D"])
        self.assertEqual((coverage["n_arm_a_predicted"], coverage["n_arm_b_predicted"],
                          coverage["n_common_groups"]), (3, 4, 3))
        self.assertTrue(coverage["below_design_coverage"])
        self.assertEqual(coverage["arm_a_noncoverage_reasons"], {"abstained_unspecified": 1})
        self.assertIsNotNone(paired_rho_difference(pair)["delta_rho"])
        self.assertEqual(strict_fixed_eight_macro({d: .1 for d in DRUGS}, DRUGS), (.1, ""))
        point, reason = strict_fixed_eight_macro({d: (.1 if d != 14 else None) for d in DRUGS}, DRUGS)
        self.assertIsNone(point)
        self.assertEqual(reason, "at_least_one_drug_undefined")
        extra = pd.concat([prediction, pd.DataFrame([
            (11, "pair", "E", "G5", .5, "predicted"),
            (11, "expression", "E", "G5", .5, "predicted"),
        ], columns=prediction.columns)], ignore_index=True)
        _, extended_coverage = paired_common_table(response, extra, drug_id=11,
            arm_a="pair", arm_b="expression", frozen_eligible_ids=["A", "B", "C", "D", "E"])
        self.assertEqual(extended_coverage["n_frozen_eligible_ids"], 5)
        self.assertEqual(extended_coverage["n_response_valid_ids"], 4)
        self.assertEqual(extended_coverage["n_arm_a_predicted_all_frozen"], 4)
        self.assertEqual(extended_coverage["n_common_groups"], 3)
        corrupted = prediction.copy()
        corrupted.loc[0, "biological_group_id"] = "WRONG"
        with self.assertRaisesRegex(ValueError, "group differs"):
            paired_common_table(response, corrupted, drug_id=11, arm_a="pair", arm_b="expression")

    def test_shared_group_exact_midrank_and_valid_draw_policy(self):
        groups = ("G0", "G1", "G2", "G3", "G4")
        y = np.array([1., 2., 3., 4., 5.])
        pa = np.array([1., 1., 3., 4., 5.])
        pb = np.array([4., 3., 2., 1., 5.])
        pair = pd.DataFrame({"biological_group_id": groups, "observed_auc": y,
                             "prediction_a": pa, "prediction_b": pb})
        counts = np.array([[2, 0, 1, 2, 0], [5, 0, 0, 0, 0]], dtype=np.uint16)
        shared = SharedBootstrap(groups, counts, 99)
        result = bootstrap_paired_rho({drug: pair.copy() for drug in DRUGS},
                                      fixed_drug_ids=DRUGS, shared=shared,
                                      min_valid_fraction=.5, family_size=9, chunk_size=1)
        indices = np.repeat(np.arange(5), counts[0])
        manual = float(spearmanr(y[indices], pa[indices]).statistic -
                       spearmanr(y[indices], pb[indices]).statistic)
        for drug in DRUGS:
            item = result["per_drug"][drug]["delta_interval"]
            self.assertEqual(item["bootstrap_valid"], 1)
            self.assertEqual(result["per_drug"][drug]["invalid_paired_draw_reasons"],
                             {"insufficient_groups": 1})
            self.assertAlmostEqual(item["nominal_95_lower"], manual)
            self.assertAlmostEqual(item["nominal_95_upper"], manual)
            self.assertAlmostEqual(item["family_adjusted_lower"], manual)
        self.assertEqual(result["strict_macro"]["delta_interval"]["bootstrap_valid"], 1)
        self.assertEqual(result["shared_multiplicity_sha256"], shared_bootstrap_sha256(shared))
        withheld = bootstrap_paired_rho({drug: pair.copy() for drug in DRUGS},
                                        fixed_drug_ids=DRUGS, shared=shared,
                                        min_valid_fraction=.99, chunk_size=1)
        self.assertEqual(withheld["strict_macro"]["delta_interval"]["ci_status"],
                         "withheld_insufficient_valid_draws")
        one_constant = {drug: pair.copy() for drug in DRUGS}
        one_constant[14]["prediction_a"] = .5
        undefined = bootstrap_paired_rho(one_constant, fixed_drug_ids=DRUGS,
                                         shared=shared, min_valid_fraction=.5, chunk_size=1)
        self.assertIsNone(undefined["strict_macro"]["delta_rho"])
        self.assertEqual(undefined["strict_macro"]["delta_interval"]["bootstrap_valid"], 0)

    def test_global_draw_schedule_reproducible_and_all_groups_retained(self):
        first = make_shared_bootstrap(["G3", "G1", "G2"], draws=12, seed=20, expected_groups=3)
        second = make_shared_bootstrap(["G2", "G3", "G1"], draws=12, seed=20, expected_groups=3)
        self.assertEqual(first.groups, ("G1", "G2", "G3"))
        np.testing.assert_array_equal(first.counts, second.counts)
        np.testing.assert_array_equal(first.counts.sum(axis=1), 3)


if __name__ == "__main__":
    unittest.main()
