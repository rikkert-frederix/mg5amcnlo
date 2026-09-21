import unittest

import numpy as np

from main_replicas import collapse_batches


class TestMainReplicaVectors(unittest.TestCase):
    def test_complete_signed_vector_and_conditional_variance(self):
        base=np.arange(1.,7.)[:,None,None]/6
        values=np.concatenate([base,-2*base])
        strata=np.repeat([0,1],6)
        mean,variance=collapse_batches(values,strata,values.sum(axis=0))
        self.assertAlmostEqual(mean[0,0],-3.5)
        expected=(1+4)*np.var(np.arange(1.,7.),ddof=1)/6
        self.assertAlmostEqual(variance[0,0],expected)
        # Cancellation-sensitive validation uses absolute contributions, not
        # an unstable relative comparison to a nearly zero final histogram.
        values=np.concatenate([base,-base])
        mean,variance=collapse_batches(values,strata,np.zeros((1,1)))
        np.testing.assert_allclose(mean,0.,atol=1.e-14)
        self.assertGreater(variance[0,0],0.)

    def test_incomplete_or_changed_batches_are_rejected(self):
        values=np.ones((6,2,3))/6
        for changed,strata,final in ((values[:4],np.zeros(4),values[:4].sum(axis=0)),
                (values,np.zeros(6),np.zeros((2,3))),
                (values,np.zeros(5),values.sum(axis=0))):
            with self.assertRaises(ValueError):collapse_batches(changed,strata,final)


if __name__=='__main__':
    unittest.main()
