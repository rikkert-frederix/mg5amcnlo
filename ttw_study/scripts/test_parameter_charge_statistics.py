"""Charge normalization and shared-reference covariance checked analytically."""
import unittest

import numpy as np

from full_flavour_statistics import FLAVOURS
from parameter_charge_statistics import (CHARGE_LABELS, LABELS, charge_prediction,
                                         estimate, pair_comparisons)
from parameter_statistics import comparison_keys


class TestParameterChargeStatistics(unittest.TestCase):
    @staticmethod
    def fixture(full=False):
        flavours = FLAVOURS if full else FLAVOURS[:2]
        comparisons = [dict(scenario=scenario, w_treatment='onshell', charge=charge,
                            flavours=[list(row) for row in flavours])
                       for scenario in ('as117', 'as119') for charge in ('plus', 'minus')]
        ensembles = {}
        spread = np.asarray([.93, .97, 1., 1.03, 1.07])
        weights = np.linspace(.8, 1.2, 82)
        weights[0] = 1.
        for comparison in comparisons:
            plus = comparison['charge'] == 'plus'
            for (origin, variant), keys in comparison_keys(comparison, full).items():
                for f, key in enumerate(keys):
                    if key in ensembles:
                        continue
                    amplitude = ((2. if f == 0 else 8.) if plus else (1. if f == 0 else 2.))
                    if variant == 'Pi':
                        amplitude *= 1.12 if plus else .97
                    if origin == 'scan':
                        amplitude *= ((1.06 if plus else .94) if comparison['scenario'] == 'as117'
                                      else (.95 if plus else 1.03))
                        if variant == 'Pi':
                            amplitude *= 1.02 if plus else .99
                    feature = amplitude*np.asarray([.8, .2, 1.] if f == 0 else [.3, .7, 1.])
                    ensembles[key] = spread[:,None,None]*feature[None,:,None]*weights[None,None,:]
        pairs, unpaired = pair_comparisons(comparisons, full)
        assert not unpaired
        return ensembles, comparisons, pairs

    def test_charge_ratios_and_combined_shapes_use_raw_flavour_sums(self):
        arrays, comparisons, pairs = self.fixture()
        result = estimate(arrays, pairs, full_flavour=False)
        s, pi = (LABELS.index(label) for label in ('reference_S', 'reference_Pi'))
        np.testing.assert_allclose(result['value'][0,:,s,2,0], [13.,10./3.,7./13.])
        self.assertAlmostEqual(result['value'][0,1,pi,2,0], (10.*1.12)/(3.*.97))
        self.assertGreater(result['mc_error'][0,1,s,2,0], 0.)
        normalized = estimate(arrays, pairs, lambda row: row[:-1]/row[-1], False)
        self.assertAlmostEqual(normalized['value'][0,0,s,0,0],5.4/13.)
        self.assertNotAlmostEqual(normalized['value'][0,0,s,0,0],(.4+1.4/3.)/2.)
        self.assertAlmostEqual(normalized['value'][0,1,s,0,0],.4/(1.4/3.))
        self.assertAlmostEqual(normalized['value'][0,2,s,0,0],(.4-1.4/3.)/(.4+1.4/3.))
        # Every charge's bin fractions sum to one; so does the combined shape.
        np.testing.assert_allclose(normalized['value'][:,0,:4].sum(axis=2),1.)
        self.assertEqual(result['value'].shape,(2,3,9,3,82))

    def test_shared_reference_covariance_and_independent_charge_permutations(self):
        arrays, comparisons, pairs = self.fixture()
        result = estimate(arrays,pairs,full_flavour=False)
        change, baseline = (LABELS.index(label) for label in ('product_shift_change', 'reference_Pi_minus_S'))
        factors = result['covariance_factor']
        covariance = np.sum(factors[:,0,:,change]*factors[:,1,:,change],axis=0)
        np.testing.assert_allclose(covariance,result['mc_error'][0,:,baseline]**2,rtol=1.e-11,atol=1.e-14)
        self.assertTrue((covariance>0.).all())
        expected_sum_variance = sum(np.var(value[:,2,0],ddof=1)/len(value)
                                    for key,value in arrays.items() if key[0]=='reference')
        self.assertAlmostEqual(covariance[CHARGE_LABELS.index('sum'),2,0],expected_sum_variance)
        changed = dict(arrays)
        for key,value in arrays.items():
            if key[3]=='minus': changed[key]=value[[4,0,2,1,3]]
        permuted = estimate(changed,pairs,full_flavour=False)
        np.testing.assert_allclose(permuted['value'],result['value'],atol=1.e-13)
        np.testing.assert_allclose(permuted['mc_error'],result['mc_error'],atol=1.e-13)

    def test_zero_and_negative_asymmetries_remain_defined(self):
        arrays, comparisons, pairs = self.fixture()
        arrays = {key:np.ones_like(value) for key,value in arrays.items()}
        result = estimate(arrays,pairs,full_flavour=False)
        np.testing.assert_allclose(result['value'][:,:,LABELS.index('product_shift_change')],0.)
        np.testing.assert_allclose(result['value'][:,2],0.)
        self.assertTrue(np.isfinite(result['mc_error']).all())
        self.assertAlmostEqual(charge_prediction(np.asarray([1.]),np.asarray([2.]))[2,0],-1./3.)
        prediction = charge_prediction(np.asarray([1.]),np.asarray([0.]))
        self.assertTrue(np.isnan(prediction[1,0]))
        self.assertEqual(prediction[2,0],1.)

    def test_pairing_and_full_flavour_guards(self):
        arrays, comparisons, pairs = self.fixture()
        with self.assertRaisesRegex(ValueError,'all eight'): estimate(arrays,pairs)
        paired, missing = pair_comparisons([row for row in comparisons if row['charge']=='plus'],False)
        self.assertFalse(paired)
        self.assertEqual(len(missing),2)
        wrong = [dict(row) for row in comparisons]
        wrong[1]['flavours'] = [list(FLAVOURS[0])]
        with self.assertRaisesRegex(ValueError,'same ordered flavour'): pair_comparisons(wrong,False)
        with self.assertRaisesRegex(ValueError,'Duplicate'): pair_comparisons(comparisons+comparisons[:1],False)
        arrays, comparisons, pairs = self.fixture(full=True)
        result = estimate(arrays,pairs)
        self.assertEqual(len(result['ensemble_order']),96)
        arrays.pop(next(iter(arrays)))
        with self.assertRaisesRegex(ValueError,'Missing'): estimate(arrays,pairs)


if __name__=='__main__':
    unittest.main()
