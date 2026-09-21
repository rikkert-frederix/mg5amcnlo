import unittest

import numpy as np

from joint_report import derived_rates, estimate, normalization_bins
from replica_statistics import ratio


class JointReportChecks(unittest.TestCase):
    def test_region_normalization_is_not_truncated_spectrum(self):
        self.assertEqual(normalization_bins(11),[4])
        self.assertEqual(normalization_bins(15),[6,7,8])
        self.assertEqual(normalization_bins(16),[7,8])
        self.assertEqual(normalization_bins(17),[5])
        self.assertEqual(normalization_bins(21),[3])
        values = np.arange(1.,11.)[:,None,None]*np.array([1.,2.,6.])[None,:,None]/5
        result = estimate(values,np.repeat([0,1],5),lambda total: ratio(total[:-1],total[-1]))
        np.testing.assert_allclose(result['value'].ravel(),[1/6,1/3])
        np.testing.assert_allclose(result['mc_error'],0.,atol=1e-15)

    def test_nested_rate_and_jet_partition(self):
        result = derived_rates(np.array([10.,8.,7.,6.,4.,1.,2.,.8,.2])[:,None])
        np.testing.assert_allclose(result.ravel(),[.6,.4,2/3,.25,.5,.2,.05])
        self.assertAlmostEqual(float(result[3:].sum()),1.)


if __name__ == '__main__':
    unittest.main()
