"""Synthetic exact checks; these do not certify physical MC convergence."""
import unittest
import numpy as np

from replica_statistics import (covariance_factor, independent_jackknife,
                                independent_jackknife_many, jackknife, ratio)


class TestReplicaStatistics(unittest.TestCase):
    def test_covariance_of_mean_and_linear_jackknife(self):
        x = np.asarray([[1.,2.], [2.,4.], [3.,6.], [4.,8.], [5.,10.]])
        factor = covariance_factor(x)
        np.testing.assert_allclose(factor.T@factor, np.cov(x.T,ddof=1)/5)
        result = jackknife(x,lambda mean: mean)
        np.testing.assert_allclose(result['value'], [3.,6.])
        np.testing.assert_allclose(result['mc_error'], np.sqrt([.5,2.]))
        np.testing.assert_allclose(result['nonlinear_bias_estimate'],0.,atol=1.e-14)

    def test_normalization_keeps_perfect_bin_rate_correlation(self):
        x = np.asarray([[1.,2.], [2.,4.], [3.,6.], [4.,8.], [5.,10.]])
        result = jackknife(x,lambda mean: ratio(mean[0],mean[1]))
        self.assertEqual(float(result['value']), .5)
        self.assertEqual(float(result['mc_error']), 0.)

    def test_independent_samples_do_not_gain_an_artificial_covariance(self):
        x = np.arange(5,dtype=float)[:,None]
        y = 2*x+10.
        result = independent_jackknife(x,y,lambda a,b: b-a)
        np.testing.assert_allclose(result['mc_error'],np.sqrt([2.5]))
        np.testing.assert_allclose(result['value'],[12.])
        swapped = independent_jackknife(x,y[::-1],lambda a,b: b-a)
        np.testing.assert_allclose(result['mc_error'],swapped['mc_error'])

    def test_missing_replicas_and_unstable_denominator(self):
        with self.assertRaises(ValueError):
            covariance_factor([[1.,2.]]*4)
        self.assertTrue(np.isnan(ratio(1.,0.)))
        self.assertTrue(np.isnan(ratio(1.,-1.)))
        x = np.asarray([[1.,-1.]]*4+[[1.,5.]])
        result = jackknife(x,lambda mean: ratio(mean[0],mean[1]))
        self.assertAlmostEqual(float(result['value']),5.)
        self.assertFalse(result['valid'])
        self.assertTrue(np.isnan(result['mc_error']))

    def test_four_independent_ensembles_for_a_treatment_shift(self):
        x = np.arange(5,dtype=float)[:,None]
        samples = [factor*x+10. for factor in (1.,2.,3.,4.)]
        transform = lambda s,p,sref,pref: (p-s)-(pref-sref)
        result = independent_jackknife_many(samples,transform)
        np.testing.assert_allclose(result['value'],[0.])
        np.testing.assert_allclose(result['mc_error'],[np.sqrt(15.)])
        self.assertEqual(result['covariance_factor'].shape,(20,1))
        samples[1] = samples[1][::-1]
        np.testing.assert_allclose(independent_jackknife_many(samples,transform)['mc_error'],
                                   result['mc_error'])


if __name__ == '__main__':
    unittest.main()
