import unittest
import json
from pathlib import Path
import tempfile
from unittest import mock

import numpy as np

from main_convergence import diagnostic, linear_summary, run


class TestMainConvergence(unittest.TestCase):
    def test_independent_linear_contrast_and_degrees_of_freedom(self):
        a=np.arange(1.,6.)[:,None]
        b=2.*a
        result=linear_summary([a,b],[1.,-1.])
        self.assertAlmostEqual(result['value'][0],-3.)
        self.assertAlmostEqual(result['mc_error'][0]**2,2.5)
        self.assertAlmostEqual(result['effective_dof'][0],2.5**2/((.5**2+2.**2)/4.))
        self.assertAlmostEqual(result['max_ensemble_variance_fraction'][0],.8)
        # Permuting runs of an independent ensemble cannot change its covariance.
        reversed_result=linear_summary([a,b[::-1]],[1.,-1.])
        for key in result:np.testing.assert_allclose(result[key],reversed_result[key])

    def test_paired_scale_cancellation_keeps_common_run_noise(self):
        noise=np.arange(1.,6.)
        scales=np.column_stack((noise,noise+0.01*np.array([-2.,-1.,0.,1.,2.])))
        result=linear_summary([(scales[:,1]-scales[:,0])[:,None]])
        self.assertAlmostEqual(result['value'][0],0.)
        self.assertAlmostEqual(result['mc_error'][0],np.sqrt(.00005))

    def test_zero_variance_is_not_a_resolved_deviation(self):
        result=linear_summary([np.ones((5,2))])
        detail=diagnostic(result,['a','b'])
        self.assertEqual(detail['tested_nonzero_error_entries'],0)
        self.assertEqual(detail['entries_above_three_mc_errors'],0)
        self.assertIsNone(detail['maximum_label'])
        self.assertTrue(np.isnan(result['effective_dof']).all())

    def test_insufficient_or_nonfinite_ensembles_fail(self):
        for rows in ([],[np.ones((4,1))],[np.full((5,1),np.nan)],
                     [np.ones((5,1)),np.ones((5,2))]):
            with self.assertRaises(ValueError):linear_summary(rows)

    def test_complete_inspection_keeps_raw_bins_and_does_not_release(self):
        from test_main_report import TestMainReport
        source,layout,vectors,conditional=TestMainReport().fixture(2)
        source.update(samples=[{}]*480,distinct_stage_initializations=960,cross_run_rng_overlap=0)
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'input.json'; path.write_text('{}\n')
            output=Path(folder)/'inspection.json'
            with mock.patch('main_convergence.read_main',return_value=(source,layout,vectors,conditional)):
                run(path,output)
            report=json.loads(output.read_text())
            self.assertFalse(report['publication_certified'])
            self.assertTrue(report['no_runs_discarded'])
            self.assertEqual(len(report['rates']),10)
            self.assertEqual(len(report['inclusive']),10)
            self.assertEqual(len(report['inclusive_charge_sum']),5)
            self.assertEqual(len(report['inclusive_decay_scales']),8)
            self.assertEqual(len(report['retraining_groups']),96)
            self.assertEqual(len(report['original_shape_precision']),20)
            with np.load(output.with_suffix('.npz')) as data:
                self.assertEqual(data['original_shape_2_values'].shape,(11,2))
                self.assertEqual(data['strict_identity_plus_value'].shape[1],183)


if __name__=='__main__':
    unittest.main()
