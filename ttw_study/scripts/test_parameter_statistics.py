"""Analytic covariance and flavour-normalization checks for parameter scans."""
import unittest

import numpy as np

from full_flavour_statistics import FLAVOURS
from parameter_statistics import LABELS, comparison_keys, contrast, estimate, group_key


class TestParameterStatistics(unittest.TestCase):
    @staticmethod
    def fixture(full=False):
        flavours = FLAVOURS if full else FLAVOURS[:2]
        comparisons = [dict(scenario=name, w_treatment='onshell', charge='plus',
                            flavours=[list(f) for f in flavours]) for name in ('as117', 'as119')]
        arrays = {}
        spread = np.asarray([.94, .97, 1., 1.03, 1.06])
        for comparison in comparisons:
            for (origin, variant), keys in comparison_keys(comparison, full).items():
                for i, key in enumerate(keys):
                    if key in arrays:
                        continue
                    amplitude = 1. if i == 0 else 9.+i-1
                    amplitude *= (1.5 if i == 0 else 1.1) if variant == 'Pi' else 1.
                    if origin == 'scan':
                        amplitude *= (1.05 if comparison['scenario'] == 'as117' else .95)
                        if variant == 'Pi': amplitude *= 1.01
                    feature = amplitude*np.asarray([1., .4, .6])
                    weights = np.linspace(.8, 1.2, 82)
                    weights[0] = 1.
                    arrays[key] = spread[:, None, None]*feature[None, :, None]*weights[None, None, :]
        return arrays, comparisons

    def test_flavours_are_summed_before_ratios_and_normalization(self):
        arrays, comparisons = self.fixture()
        result = estimate(arrays, comparisons, full_flavour=False)
        expected = 11.4/10.
        i = LABELS.index('reference_Pi_over_S')
        self.assertAlmostEqual(result['value'][0, i, 0, 0], expected)
        self.assertNotAlmostEqual(result['value'][0, i, 0, 0], (1.5+1.1)/2)
        normalized = estimate(arrays, comparisons, transform=lambda row: row[1:]/row[0], full_flavour=False)
        for label in ('reference_S', 'reference_Pi', 'scan_S', 'scan_Pi'):
            np.testing.assert_allclose(normalized['value'][:, LABELS.index(label), :, 0], [[.4, .6]]*2)
            self.assertLess(float(normalized['mc_error'][:, LABELS.index(label)].max()), 1.e-14)

    def test_shared_baseline_covariance_and_independent_replica_order(self):
        arrays, comparisons = self.fixture()
        result = estimate(arrays, comparisons, full_flavour=False)
        index = LABELS.index('product_shift_change')
        factor = result['covariance_factor'][:, :, index, 0, 0]
        expected = sum(np.var(values[:, 0, 0], ddof=1)/len(values)
                       for key, values in arrays.items() if key[0] == 'reference')
        self.assertAlmostEqual(factor[:, 0]@factor[:, 1], expected, places=12)
        self.assertGreater(expected, 0.)
        permuted = dict(arrays)
        key = next(key for key in arrays if key[0] == 'scan')
        permuted[key] = permuted[key][[4, 0, 2, 1, 3]]
        changed = estimate(permuted, comparisons, full_flavour=False)
        np.testing.assert_allclose(changed['value'], result['value'], atol=1.e-13)
        np.testing.assert_allclose(changed['mc_error'], result['mc_error'], atol=1.e-13)
        self.assertEqual(len(result['ensemble_order']), 12)

    def test_full_direct_flavour_scope_and_zero_denominator_guards(self):
        arrays, comparisons = self.fixture()
        with self.assertRaisesRegex(ValueError, 'all eight'): estimate(arrays, comparisons)
        arrays, comparisons = self.fixture(full=True)
        result = estimate(arrays, comparisons)
        self.assertEqual(len(result['ensemble_order']), 48)
        self.assertTrue(np.isfinite(result['value']).all())
        arrays.pop(next(iter(arrays)))
        with self.assertRaisesRegex(ValueError, 'Missing'): estimate(arrays, comparisons)
        result = contrast(np.array([0.]), np.array([1.]), np.array([1.]), np.array([1.1]))
        self.assertTrue(np.isnan(result[LABELS.index('reference_Pi_over_S'), 0]))
        self.assertTrue(np.isnan(result[LABELS.index('product_ratio_double_ratio'), 0]))


if __name__ == '__main__':
    unittest.main()
