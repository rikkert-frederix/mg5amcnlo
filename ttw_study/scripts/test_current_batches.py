from pathlib import Path
import unittest

import numpy as np

from campaign import ROOT
from compare_current_batches import absolute_inputs, contrast
from joint_comparison import paired_variants


class TestCurrentBatches(unittest.TestCase):
    def test_physical_input_aliases_are_normalized_without_relabelling_hashes(self):
        path='ttw_study/inputs/benchmark.json'
        self.assertEqual(absolute_inputs({path:'digest'}),{str(ROOT/path):'digest'})
        self.assertEqual(absolute_inputs({str(ROOT/path):'digest'}),{str(ROOT/path):'digest'})
        with self.assertRaises(ValueError):
            absolute_inputs({path:'digest',str(ROOT/path):'digest'})

    def test_independent_strata_preserve_weight_correlation_and_branching_squared(self):
        rows=np.arange(1.,11.)
        native=np.stack([rows,2*rows],axis=1)/5
        decayed=native*.01+np.asarray([.002,.004])
        strata=np.repeat([0,1],5)
        result=paired_variants(native,decayed,strata,strata,lambda a,b:contrast(a,b,.01))
        np.testing.assert_allclose(result['value'][0],native.sum(axis=0)*.01)
        np.testing.assert_allclose(result['value'][1],decayed.sum(axis=0))
        np.testing.assert_allclose(result['value'][2],decayed.sum(axis=0)-native.sum(axis=0)*.01)
        self.assertGreater(result['mc_error'][2,0],0.)
        np.testing.assert_allclose(result['mc_error'][2,1],2*result['mc_error'][2,0])


if __name__=='__main__':
    unittest.main()
