"""Nominal covariance storage preserves values and releases full parents."""
import gc
import unittest
import weakref

import numpy as np

from covariance_storage import nominal_copy


class TestCovarianceStorage(unittest.TestCase):
    def test_nominal_storage_does_not_retain_all_weight_coordinates(self):
        for weights in (82, 183):
            original = np.arange(5*13*4*weights, dtype=float).reshape(5, 13, 4, weights).copy()
            original[1, 2, 3, 0] = np.nan
            expected = original[..., 0].copy()
            nominal = nominal_copy(original)
            np.testing.assert_array_equal(nominal, expected)
            self.assertTrue(nominal.flags.owndata)
            self.assertFalse(np.shares_memory(nominal, original))
            self.assertEqual(original.nbytes, weights*nominal.nbytes)
            parent = weakref.ref(original)
            del original
            gc.collect()
            self.assertIsNone(parent())
            np.testing.assert_array_equal(nominal, expected)

    def test_covariance_and_bias_values_are_unchanged(self):
        rng = np.random.default_rng(20260914)
        factor = rng.normal(size=(5, 13, 3, 183))
        old = factor[..., 0].reshape(5, -1)
        new = nominal_copy(factor).reshape(5, -1)
        np.testing.assert_array_equal(new, old)
        # Strided/contiguous BLAS paths may accumulate in different orders.
        covariance = old.T@old
        tolerance = 10*np.finfo(float).eps*np.max(np.abs(covariance))
        np.testing.assert_allclose(new.T@new, covariance, rtol=0., atol=tolerance)
        bias = rng.normal(size=(13, 3, 183))
        np.testing.assert_array_equal(nominal_copy(bias), bias[..., 0])


if __name__ == '__main__':
    unittest.main()
