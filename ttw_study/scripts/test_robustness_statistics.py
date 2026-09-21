"""Analytical full-flavour, shared-state covariance and charge checks."""
import copy
import unittest

import numpy as np

from full_flavour_statistics import FLAVOURS
from replica_statistics import ratio
from robustness_statistics import (LABELS, CHARGE_CONTRAST_LABELS, comparison_keys,
    estimate, estimate_charges, group_key, pair_comparisons, select)
from run_robustness_campaign import STATES, comparisons as all_comparisons


def fixture(full=False, both=False, weights=3):
    flavours = FLAVOURS if full else FLAVOURS[:2]
    comparisons = [dict(row, flavours=[list(f) for f in flavours]) for row in all_comparisons()
                   if both or row['charge'] == 'plus']
    arrays = {}
    spread = np.asarray([.8, .9, 1., 1.1, 1.2])
    for c, charge in enumerate(('plus', 'minus') if both else ('plus',)):
        for s, (mode, mass) in enumerate(STATES):
            for v, variant in enumerate(('S', 'Pi')):
                for f, flavour in enumerate(flavours):
                    amplitude = (f+1)*(1.+.1*s)*(1.+v*(.05+.01*f+.02*s))*(1.+.7*c)
                    fraction = .2+.03*f+.04*c
                    features = amplitude*np.asarray([fraction, 1.-fraction, 1.])
                    key = group_key(mode, mass, charge, variant, flavour)
                    arrays[key] = spread[:, None, None]*features[None, :, None]*np.linspace(1., 1.2, weights)[None, None]
    return arrays, comparisons


class TestRobustnessStatistics(unittest.TestCase):
    def test_full_flavour_sum_precedes_ratios_and_normalization(self):
        arrays, comparisons = fixture(full=True, both=True, weights=1)
        result = estimate(arrays, comparisons)
        self.assertEqual(len(result['ensemble_order']), 160)
        index = LABELS.index('reference_Pi_over_S')
        expected = sum((f+1)*(1.05+.01*f) for f in range(8))/sum(range(1, 9))
        self.assertAlmostEqual(result['value'][0, index, -1, 0], expected)
        self.assertNotAlmostEqual(expected, np.mean([1.05+.01*f for f in range(8)]))
        chosen = [comparisons[0]]
        normalized = estimate(select(arrays, chosen), chosen, lambda row: ratio(row[:-1], row[-1]))
        expected_shape = sum((f+1)*(.2+.03*f) for f in range(8))/sum(range(1, 9))
        self.assertAlmostEqual(normalized['value'][0, 0, 0, 0], expected_shape)
        factor = normalized['covariance_factor'][:, 0, 0]
        np.testing.assert_allclose(factor.sum(axis=1), 0., atol=1.e-15)

    def test_shared_intermediate_and_baseline_covariances_and_closure(self):
        arrays, comparisons = fixture()
        result = estimate(arrays, comparisons, full_flavour=False)
        for label in ('product_shift_change', 'product_ratio_change'):
            i = LABELS.index(label)
            np.testing.assert_allclose(result['value'][2, i], result['value'][0, i]+result['value'][1, i], atol=1.e-14)
            np.testing.assert_allclose(result['covariance_factor'][:, 2, i],
                result['covariance_factor'][:, 0, i]+result['covariance_factor'][:, 1, i], atol=2.e-14)
        factor = result['covariance_factor'][:, :, LABELS.index('product_shift_change'), -1, 0]
        def variance(state):
            return sum(np.var(row[:, -1, 0], ddof=1)/len(row) for key, row in arrays.items() if key[:2] == state)
        self.assertAlmostEqual(factor[:, 0]@factor[:, 1], -variance(STATES[1]), places=12)
        self.assertAlmostEqual(factor[:, 0]@factor[:, 3], variance(STATES[0]), places=12)
        self.assertAlmostEqual(factor[:, 1]@factor[:, 4], -variance(STATES[2]), places=12)
        changed = dict(arrays)
        for i, key in enumerate(changed):
            changed[key] = changed[key][np.roll(np.arange(5), i % 5)]
        permuted = estimate(changed, comparisons, full_flavour=False)
        np.testing.assert_allclose(permuted['value'], result['value'], atol=2.e-14)
        np.testing.assert_allclose(permuted['mc_error'], result['mc_error'], atol=2.e-14)

    def test_charge_normalization_and_zero_asymmetry_use_absolute_changes(self):
        arrays, comparisons = fixture(both=True)
        pairs, unpaired = pair_comparisons(comparisons, full_flavour=False)
        self.assertFalse(unpaired)
        selected = [row for row in comparisons if row['name'] == 'all_W']
        pairs = [row for row in pairs if row['name'] == 'all_W']
        data = select(arrays, selected, full_flavour=False)
        normalized = estimate_charges(data, pairs, lambda row: ratio(row[:-1], row[-1]), full_flavour=False)
        means = {}
        for charge in ('plus', 'minus'):
            means[charge] = sum(value.mean(axis=0) for key, value in data.items()
                               if key[:2] == STATES[0] and key[2:4] == (charge, 'S'))
        plus, minus = means['plus'], means['minus']
        expected = (plus+minus)[0, 0]/(plus+minus)[-1, 0]
        self.assertAlmostEqual(normalized['value'][0, 0, 0, 0, 0], expected)
        self.assertNotAlmostEqual(expected, .5*(plus[0, 0]/plus[-1, 0]+minus[0, 0]/minus[-1, 0]))
        # Equal reference charge means give zero A_W without invalidating
        # any absolute change or assigning a ratio to that zero asymmetry.
        changed = dict(data)
        for key in data:
            if key[:2] == STATES[0] and key[2] == 'minus':
                changed[key] = data[key[:2]+('plus',)+key[3:]].copy()
        result = estimate_charges(changed, pairs, full_flavour=False)
        self.assertEqual(len(CHARGE_CONTRAST_LABELS), 9)
        self.assertAlmostEqual(result['value'][0, 2, 0, -1, 0], 0.)
        self.assertTrue(np.isfinite(result['value']).all())
        self.assertTrue(np.isfinite(result['mc_error']).all())

    def test_missing_flavours_retrainings_and_mislabeled_states_are_rejected(self):
        arrays, comparisons = fixture()
        with self.assertRaisesRegex(ValueError, 'all eight'):
            estimate(arrays, comparisons)
        broken = dict(arrays)
        first = next(iter(broken))
        broken[first] = broken[first][:4]
        with self.assertRaisesRegex(ValueError, 'five finite'):
            estimate(broken, comparisons, full_flavour=False)
        broken = copy.deepcopy(comparisons[0])
        broken['target'] = comparisons[1]['target']
        with self.assertRaisesRegex(ValueError, 'name does not match'):
            comparison_keys(broken, full_flavour=False)
        zero = dict(arrays)
        for key in zero:
            if key[:2] == STATES[0] and key[3] == 'S':
                zero[key] = np.zeros_like(zero[key])
        result = estimate(zero, comparisons, full_flavour=False)
        self.assertTrue(np.isnan(result['value'][0, LABELS.index('reference_Pi_over_S')]).all())
        pairs, unpaired = pair_comparisons(comparisons, full_flavour=False)
        self.assertEqual(pairs, [])
        self.assertEqual(len(unpaired), 5)


if __name__ == '__main__':
    unittest.main()
