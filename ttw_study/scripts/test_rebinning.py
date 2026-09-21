import unittest

import numpy as np

from joint_report import estimate
from rebinning import rebin, spans


class RebinningChecks(unittest.TestCase):
    def test_exact_sums_for_all_replicas_and_weights(self):
        edges=np.array([[0.,1.],[1.,2.],[2.,3.],[3.,4.]])
        values=np.arange(5*4*3).reshape(5,4,3)
        result=rebin(values,edges,[0.,2.,4.])
        np.testing.assert_array_equal(result[:,0],values[:,:2].sum(axis=1))
        np.testing.assert_array_equal(result.sum(axis=1),values.sum(axis=1))
        for invalid in ([0.,1.5,4.],[1.,2.,4.],[0.,4.,2.],[0.,2.,3.]):
            with self.assertRaises(ValueError):
                spans(edges,invalid)

    def test_signed_bin_correlations_are_retained(self):
        rows=np.arange(1.,11.)
        values=np.stack([rows,20.-rows],axis=1)[:,:,None]/5
        merged=rebin(values,[[0.,1.],[1.,2.]],[0.,2.])
        joint=estimate(merged,np.repeat([0,1],5))
        np.testing.assert_allclose(joint['value'],[[40.]])
        np.testing.assert_allclose(joint['mc_error'],0.,atol=1.e-14)
        unmerged=estimate(values,np.repeat([0,1],5))
        self.assertTrue((unmerged['mc_error'] > 0.).all())


if __name__ == '__main__':
    unittest.main()
