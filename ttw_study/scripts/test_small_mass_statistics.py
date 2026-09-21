import unittest

import numpy as np

from small_mass_report import estimate


class TestSmallMassStatistics(unittest.TestCase):
    def test_common_massless_control_correlates_mass_differences(self):
        base=np.arange(1.,7.)[:,None,None]/6
        vectors=[scale*base for scale in (1.,2.,3.,4.,5.,6.)]
        strata=[np.zeros(6)]*6
        result=estimate(vectors,strata)
        # The change in Pi-S is (Pim-Sm)-(Pi0-S0); its two mass versions
        # share exactly the independent Pi0 and S0 variances.
        covariance=result['covariance_factor'][:,:,4,0,0]
        expected=sum(np.var(v[:,0,0]*6,ddof=1)/6 for v in vectors[:2])
        self.assertAlmostEqual(np.dot(covariance[:,0],covariance[:,1]),expected)
        shuffled=estimate([v[::-1] if i%2 else v for i,v in enumerate(vectors)],strata)
        np.testing.assert_allclose(result['value'],shuffled['value'])
        np.testing.assert_allclose(result['mc_error'],shuffled['mc_error'])

    def test_bin_parent_correlations_cancel_before_mass_ratios(self):
        base=np.arange(1.,7.)/6
        vectors=[np.stack([scale*base,2*scale*base,4*scale*base],axis=1)[...,None]
                 for scale in (1.,2.,3.,4.,5.,6.)]
        result=estimate(vectors,[np.zeros(6)]*6,normalized=True)
        np.testing.assert_allclose(result['value'][:,:2],1.)
        np.testing.assert_allclose(result['mc_error'],0.,atol=1.e-14)
        with self.assertRaises(ValueError):estimate(vectors[:4],[np.zeros(6)]*4)


if __name__=='__main__':
    unittest.main()
