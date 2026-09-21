import unittest
from unittest.mock import patch

import numpy as np

from inclusive_report import STABLE_VARIANT, compare, decay_scale_residuals, input_bin
from joint_report import estimate
from load_results import POINTS


class TestInclusiveNormalization(unittest.TestCase):
    def test_reference_orders(self):
        self.assertEqual(STABLE_VARIANT,dict(LO='LO',D='LO',PiD='LO',P='P',S='P',Pi='P'))

    def test_decay_differences_at_fixed_production_point(self):
        values = np.ones(183)
        for i,(r,f,t,a) in enumerate(POINTS,1):
            values[i] = 100*r+10*f+2*t+a
        np.testing.assert_allclose(decay_scale_residuals(values),[2*(t-1)+(a-1) for _,_,t,a in POINTS])

    def test_joint_cancellation_preserves_weight_correlations(self):
        batches = np.asarray([np.arange(183)+i for i in range(6)],dtype=float)
        reduced = estimate(batches/6,np.zeros(6,dtype=int),decay_scale_residuals)
        np.testing.assert_allclose(reduced['mc_error'],0.,atol=1.e-12)
        np.testing.assert_allclose(reduced['value'],decay_scale_residuals(batches.mean(axis=0)))

    def test_input_is_before_all_selections_and_agrees_across_cuts(self):
        result = dict(values=np.ones((10,183)),offsets=np.arange(11),
                      titles=['W+ rates: cut%d' % i for i in range(5)]+['other']*5)
        self.assertEqual(input_bin(result),0)
        result['values'][3,100] *= 2
        with self.assertRaises(ValueError):
            input_bin(result)

    def test_branching_normalization_and_mismatch_guards(self):
        manifest = dict(w_treatment='onshell',production_scale='core-w-ht-half',
                        settings=dict(lhaid=331700,iseed=1),w_width=2.,run_name='S')
        result = dict(values=np.full((5,183),.125),errors=np.full(5,.01),offsets=np.arange(6),
                      titles=['W+ rates: cut%d' % i for i in range(5)],
                      report=dict(variant='S',manifest=manifest,path='/fake/process/Events/S/MADatNLO.HwU',
                                  execution={'export_generation':{'args':dict(
                                      charge='plus',flavours=['e','e','mu'])}}))
        stable = dict(values=np.ones(183),error=.01,manifest={'settings':dict(lhaid=331700,iseed=2)},
                      execution=dict(variant='P',charge='plus',branching_e=.5),
                      reference=dict(w_width=2.),parameters={})
        with patch('inclusive_report.production_parameters',return_value={}):
            row = compare(result,stable,.5)
            self.assertEqual(row['ratio']['scale']['combined63']['envelope'],[1.,1.])
            self.assertEqual(row['central_pull'],0.)
            self.assertAlmostEqual(row['central_mc_error_pb'],np.hypot(.01,.125*.01))
            for key,value in (('w_width',3.),('w_treatment','all-bw')):
                original = manifest[key]
                manifest[key] = value
                with self.assertRaises(ValueError):
                    compare(result,stable,.5)
                manifest[key] = original
            manifest['w_treatment'] = 'top-bw'
            with patch('inclusive_report.validate_bw_normalization',return_value={'verified':True}) as proof:
                row = compare(result,stable,.5)
                proof.assert_called_once_with(manifest,.5)
                self.assertEqual(row['central_pull'],0.)
                self.assertTrue(row['bw_current_normalization']['verified'])
            with patch('inclusive_report.validate_bw_normalization',side_effect=ValueError('unverified current')):
                with self.assertRaises(ValueError):
                    compare(result,stable,.5)
            manifest['w_treatment'] = 'onshell'
            stable['manifest']['settings']['iseed'] = 1
            with self.assertRaises(ValueError):
                compare(result,stable,.5)


if __name__ == '__main__':
    unittest.main()
