import unittest

import numpy as np

from joint_comparison import comparison_values,paired_variants
from replica_statistics import ratio


class TestJointVariantComparison(unittest.TestCase):
    def test_independent_variants_are_not_arbitrarily_paired(self):
        first = np.arange(1.,7.)[:,None]
        second = 2*first
        strata = np.zeros(6,dtype=int)
        result = paired_variants(first/6,second/6,strata,strata,comparison_values)
        self.assertAlmostEqual(result['value'][2,0],3.5)
        expected = np.sqrt(np.var(first[:,0],ddof=1)/6+np.var(second[:,0],ddof=1)/6)
        self.assertAlmostEqual(result['mc_error'][2,0],expected)
        self.assertGreater(result['mc_error'][3,0],0.)
        reordered = paired_variants(first/6,second[::-1]/6,strata,strata,comparison_values)
        np.testing.assert_allclose(result['mc_error'],reordered['mc_error'])

    def test_bin_normalization_correlations_remain_within_each_batch(self):
        base = np.arange(1.,7.)
        first = np.stack([base,2*base],axis=1)
        second = np.stack([3*base,4*base],axis=1)
        strata = np.zeros(6,dtype=int)
        result = paired_variants(first/6,second/6,strata,strata,
                                 lambda s,p:comparison_values(ratio(s[0],s[1]),ratio(p[0],p[1])))
        np.testing.assert_allclose(result['value'],[.5,.75,.25,1.5])
        np.testing.assert_allclose(result['mc_error'],0.,atol=1.e-14)


if __name__ == '__main__':
    unittest.main()
