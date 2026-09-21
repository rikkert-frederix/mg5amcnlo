import unittest
import numpy as np
from narrow_w_statistics import comparisons, estimate, estimate_retrained_reference, width_factors


class TestNarrowWStatistics(unittest.TestCase):
    def test_rate_power_is_applied_once_before_ratios(self):
        factors=[1.,.2,.05]
        values=[np.asarray([2.,3.])/e**3 for e in factors]
        result=comparisons(values,factors)
        np.testing.assert_allclose(result[0],[[2.,3.]]*3)
        np.testing.assert_allclose(result[1],0.,atol=1.e-14)
        np.testing.assert_allclose(result[2],1.)

    def test_independent_samples_share_reference_covariance(self):
        factors=[1.,.2,.05]
        base=np.arange(1.,7.)[:,None]/6
        vectors=[base/e**3 for e in factors]
        strata=[np.zeros(6,dtype=int)]*3
        result=estimate(vectors,strata,factors)
        reordered=estimate([vectors[0],vectors[1][::-1],vectors[2]],strata,factors)
        np.testing.assert_allclose(result['mc_error'],reordered['mc_error'])
        self.assertEqual(result['mc_error'][2,0,0],0.)
        self.assertGreater(result['mc_error'][2,1,0],0.)
        cov=result['covariance_factor'][:,2,:,0]
        self.assertGreater(np.dot(cov[:,1],cov[:,2]),0.)
        np.testing.assert_allclose(result['value'][2],1.)

    def test_normalization_keeps_parent_rate_correlations(self):
        factors=[1.,.2,.05]
        base=np.arange(1.,7.)/6
        vectors=[np.stack([base,2*base,4*base],axis=1)[:,:,None]/e**3 for e in factors]
        result=estimate(vectors,[np.zeros(6,dtype=int)]*3,factors,normalized=True)
        np.testing.assert_allclose(result['value'][0,:,:,0],[[.25,.5]]*3)
        np.testing.assert_allclose(result['mc_error'],0.,atol=1.e-14)

    def test_invalid_factors_and_missing_normalization_are_rejected(self):
        for factors in ([.2,.05],[1.,0.],[1.,-1.],[1.,np.nan],[1.,2.]):
            with self.assertRaises(ValueError):
                width_factors(factors,2)
        with self.assertRaises(ValueError):
            comparisons([np.ones(1),np.ones(1)],[1.,.2],normalized=True)

    def test_retrained_reference_mean_and_uncertainty_are_not_summed(self):
        reference=np.arange(8.,13.)[:,None,None]
        widths=[1.,.2,.05]
        conditional=[np.full((6,1,1),10./6/e**3) for e in widths[1:]]
        result=estimate_retrained_reference(reference,conditional,[np.zeros(6)]*2,widths)
        np.testing.assert_allclose(result['value'][0,:,0,0],[10.,10.,10.])
        expected_variance=np.var(reference[:,0,0],ddof=1)/len(reference)
        self.assertAlmostEqual(result['mc_error'][0,0,0,0]**2,expected_variance)
        np.testing.assert_allclose(result['mc_error_conditional_samples'],0.,atol=1.e-13)
        self.assertAlmostEqual(result['mc_error'][1,1,0,0]**2,expected_variance)
        factors=result['covariance_factor'][:,1,:,0,0]
        self.assertAlmostEqual(np.dot(factors[:,1],factors[:,2]),expected_variance)
        self.assertEqual(result['mc_error'][2,0,0,0],0.)
        self.assertEqual(result['mc_error'][1,0,0,0],0.)

    def test_retraining_normalization_and_permutation_keep_correlations(self):
        reference=np.stack([np.arange(8.,13.),2*np.arange(8.,13.),4*np.arange(8.,13.)],axis=1)[...,None]
        base=np.arange(1.,7.)
        conditional=np.stack([base,2*base,4*base],axis=1)[...,None]/6/.2**3
        result=estimate_retrained_reference(reference,[conditional],[np.zeros(6)],[1.,.2],True)
        np.testing.assert_allclose(result['value'][0,:,:,0],[[.25,.5]]*2)
        np.testing.assert_allclose(result['mc_error'],0.,atol=1.e-14)
        a=estimate_retrained_reference(reference,[conditional],[np.zeros(6)],[1.,.2])
        b=estimate_retrained_reference(reference[::-1],[conditional[::-1]],[np.zeros(6)],[1.,.2])
        np.testing.assert_allclose(a['value'],b['value'])
        np.testing.assert_allclose(a['mc_error'],b['mc_error'])
        np.testing.assert_allclose(a['mc_error']**2,
            a['mc_error_reference_retraining']**2+a['mc_error_conditional_samples']**2)
        with self.assertRaises(ValueError):
            estimate_retrained_reference(reference[:4],[conditional],[np.zeros(6)],[1.,.2])


if __name__=='__main__':
    unittest.main()
