"""Historical standalone reproducer for the exact-constant maintenance edge.

The maintained default regression cases are in test_pilot_v2_integrity.py.
This standalone fixture preserves the original post-run diagnosis; execute
it manually to compare the historical execution and repaired implementations.
Only synthetic data is used; no original response/expression data is accessed.
"""
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from modeling.reo_rf import TrainOnlyReoPairSelector
from modeling.nested_cv import TrainOnlyRegressor


class ExactWithinTissueConstantRegression(unittest.TestCase):
    def test_affine_response_standardization_cannot_create_pair_information(self):
        rng = np.random.default_rng(903)
        for counts in ((16, 16, 16), (11, 13, 17), (7, 9, 13)):
            with self.subTest(counts=counts):
                tissue = np.repeat(np.arange(3), counts)
                d = np.eye(3)[tissue]
                z = rng.integers(-1, 2, size=(len(tissue), 8)).astype(np.int8)
                names = tuple(f"REO:A{i}>B{i}" for i in range(8))
                original = tissue.astype(float)
                standardized = (original - original.mean()) / original.std()
                for response in (original, standardized):
                    selector = TrainOnlyReoPairSelector(3).fit(z, response, tissue_design=d, feature_names=names)
                    self.assertEqual(len(selector.selected_indices_), 0)

    def test_tissue_only_response_has_exact_zero_alpha_max_in_en(self):
        rng = np.random.default_rng(902)
        for counts in ((16, 16, 16), (7, 9, 13)):
            tissue = np.repeat(np.arange(3), counts)
            genes = rng.lognormal(1.0, 1.0, (len(tissue), 6)).astype(np.float64)
            x = np.column_stack((genes, tissue, np.arange(len(tissue)))).astype(np.float64)
            y = tissue.astype(float)
            config = {"feature": {"expression_genes": 6, "reo_genes": 6, "min_tpm": 1., "min_detection_fraction": 0.1, "max_pair_tie_fraction": 0.9},
                      "elastic_net": {"regularization_mode": "alpha_max_relative", "alpha_ratios": [1., 0.1], "l1_ratios": [0.5, 1.], "max_iter": 1000, "tol": 0.001}}
            for kind in ("expression", "reo"):
                with self.subTest(counts=counts, kind=kind):
                    model = TrainOnlyRegressor(kind, config, tuple(f"G{i}" for i in range(6)), alpha_ratio=0.1, l1_ratio=0.5).fit(x, y)
                    self.assertEqual(model.audit()["alpha_max"], 0.0)
                    np.testing.assert_array_equal(model.regressor_.coef_, np.zeros_like(model.regressor_.coef_))
                    predictions = model.predict(x)
                    for category in np.unique(tissue):
                        values = predictions[tissue == category]
                        np.testing.assert_array_equal(values, np.full(len(values), values[0]))


if __name__ == "__main__":
    unittest.main()
