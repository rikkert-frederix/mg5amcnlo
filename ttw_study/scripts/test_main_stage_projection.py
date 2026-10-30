import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from load_results import POINTS
from main_stage_projection import CENTRAL, PRODUCTION, TERMS, project, run, stage_contribution


class TestMainStageProjection(unittest.TestCase):
    def fixture(self,product):
        ratios={.5:1.11,1.:1.,2.:.91};g=-.08
        q=np.array([3.,-.2,-.3,.04] if product else [3.,-.2,-.3])
        x=np.zeros((5,2,183))
        for j,p in enumerate(POINTS,1):
            t,a=ratios[p[2]],ratios[p[3]]
            y=q[0]+t*q[1]+a*q[2]
            if product:y=(y+t*a*q[3])/((1+g*t)*(1+g*a))
            x[:,:,j]=y
        x[:,:,0]=x[:,:,1+POINTS.index((1.,1.,1.,1.))]
        return x,q,ratios,g

    def test_product_and_strict_projection_at_fixed_production_coordinates(self):
        for product in (True,False):
            x,q,r,g=self.fixture(product)
            result,check=project(x,r,g,product)
            np.testing.assert_allclose(result,np.broadcast_to(q,result.shape),rtol=1e-11,atol=1e-12)
            self.assertLess(check['maximum_absolute_closure_pb'],1e-12)
            self.assertEqual(result.shape,(5,2,len(PRODUCTION),len(q)))

    def test_nonpolynomial_scale_response_is_rejected(self):
        x,_,r,g=self.fixture(True)
        x[:,:,1+POINTS.index((1.,1.,1.,1.))]+=.01
        with self.assertRaisesRegex(ValueError,'do not close'):project(x,r,g,True)

    def test_width_matched_inclusive_cancellation_and_stage_orders(self):
        g=-.08;born=2.;production_correction=.7
        u=np.array([born,g*born,g*born,g*g*born])
        v=u*(born+production_correction)/born
        terms=stage_contribution(v,'Pi',g)+stage_contribution(u,'PiD',g)
        out=dict(zip(TERMS,terms))
        self.assertAlmostEqual(out['A0'],born)
        self.assertAlmostEqual(out['100'],production_correction)
        self.assertAlmostEqual(out['strict_reconstructed'],born+production_correction)
        for name in ('production_decay_second_order','decay_second_order','total_second_order',
                     'delta_reconstructed','higher_remainder'):
            self.assertAlmostEqual(out[name],0.)

    def test_fiducial_second_order_and_exact_shift_close(self):
        u=np.array([2.,-.12,-.17,.007]);v=np.array([2.7,-.23,-.19,.011]);g=-.08
        terms=stage_contribution(v,'Pi',g)+stage_contribution(u,'PiD',g)
        out=dict(zip(TERMS,terms))
        self.assertAlmostEqual(out['product_reconstructed'],sum(v)/(1+g)**2)
        self.assertAlmostEqual(out['strict_reconstructed'],v[0]+u[1]+u[2]-2*g*u[0])
        self.assertAlmostEqual(out['total_second_order'],out['production_decay_second_order']+out['decay_second_order'])
        self.assertAlmostEqual(out['delta_reconstructed'],out['total_second_order']+out['higher_remainder'])

    def test_complete_reporting_keeps_controls_and_does_not_release(self):
        from test_main_report import TestMainReport
        source,layout,vectors,conditional=TestMainReport().fixture(2)
        # Constant strict weights and width-normalized constant numerators
        # satisfy the projection while allowing independent run fluctuations.
        from campaign import STUDY
        from main_stage_projection import width_parameters
        ratios,g=width_parameters(json.loads((STUDY/'inputs/benchmark.json').read_text()))
        for key,x in vectors.items():
            for j,p in enumerate(POINTS,1):
                x[:,:,j]=x[:,:,0]
                if key[1] in ('Pi','PiD'):
                    x[:,:,j]*=(1+g)**2/((1+g*ratios[p[2]])*(1+g*ratios[p[3]]))
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'input.json';path.write_text('{}\n')
            output=Path(folder)/'projection.json'
            with mock.patch('main_stage_projection.read_main',return_value=(source,layout,vectors,conditional)):
                run(path,output)
            report=json.loads(output.read_text())
            self.assertFalse(report['publication_certified'])
            self.assertEqual(report['new_integrations'],0)
            self.assertEqual(len(report['rates']),10)
            self.assertEqual(len(report['closure_checks']),64)
            self.assertEqual(len(report['strict_decay_diagnostics']),4)
            row=next(iter(report['rates'].values()))['input']
            self.assertEqual(set(row['independent_control_residuals']),{'S','P','D','LO'})
            with np.load(output.with_suffix('.npz')) as data:
                self.assertEqual(data['plus__term_values'].shape,(90,9,len(TERMS)))
            with self.assertRaisesRegex(ValueError,'overwrite'):run(path,output)


if __name__=='__main__':
    unittest.main()
