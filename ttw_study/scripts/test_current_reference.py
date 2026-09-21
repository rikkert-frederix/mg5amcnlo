import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from campaign import digest
from run_current_reference import commands, compare, settings


class TestCurrentReference(unittest.TestCase):
    def test_native_current_and_matching_fixed_settings(self):
        command=commands(Path('/tmp/independent_current_test'))
        self.assertIn('generate p p > t t~ mu+ vm QCD=2 QED=2 [QCD]',command)
        self.assertNotIn('t >',command)
        self.assertIn('define p = g u c d s b u~ c~ d~ s~ b~',command)
        reference=dict(settings=dict(lhaid=331700,ebeam1=6500.,ebeam2=6500.,
                                     fixed_ren_scale=False,iseed=1))
        original=copy.deepcopy(reference)
        result=settings(reference,.005,55001)
        self.assertEqual(reference,original)
        self.assertTrue(all(result[k] for k in ('fixed_ren_scale','fixed_fac_scale','fixed_qes_scale')))
        self.assertEqual([result[k] for k in ('mur_ref_fixed','muf_ref_fixed','qes_ref_fixed')],
                         [212.6925]*3)
        self.assertEqual((result['lhaid'],result['iseed'],result['req_acc_fo']),(331700,55001,.005))

    def test_only_two_branching_factors_and_independent_nominal_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'inputs').mkdir()
            benchmark=root/'inputs/benchmark.json'
            benchmark.write_text(json.dumps(dict(W_width=dict(branching_e=.1))))
            native_output=root/'native.HwU'
            native_output.write_text('synthetic native output')
            native_path=root/'execution.json'
            native=dict(status='finished',output=str(native_output),output_sha256=digest(native_output),
                        canonical_weights_pb=[.1]*183,nominal_mc_error_pb=.001)
            native_path.write_text(json.dumps(native))
            options=settings(dict(settings=dict(lhaid=331700)),.005,55001)
            manifest=dict(files={},settings=options,internal_W_width_GeV=2.,no_branching_factor_applied=True,
                          benchmark_sha256=digest(benchmark))
            (root/'manifest.json').write_text(json.dumps(manifest))
            audit=root/'audit.json'
            audit.write_text('synthetic audit')
            dm=dict(variant='LO',w_treatment='all-bw',production_scale='fixed',decay_bottom_mass=0.,
                    run_name='LO',settings=dict(options,iseed=55002),w_width=2.)
            result=dict(report=dict(manifest=dm,path=str(root/'process/Events/LO/result.HwU'),
                                    execution=dict(export_generation=dict(args=dict(
                                        charge='plus',flavours=['e','e','mu'])))),
                        values=np.full((1,183),.00102),errors=np.asarray([.000003]))
            with patch('run_current_reference.STUDY',root),patch('run_current_reference.load',return_value=result),\
                 patch('run_current_reference.input_bin',return_value=0),\
                 patch('run_current_reference.production_parameters',return_value='same'),\
                 patch('run_current_reference.validate_bw_normalization',return_value=dict(proof='fixture')):
                output=root/'comparison.json'
                compare(native_path,audit,output)
                report=json.loads(output.read_text())
                self.assertAlmostEqual(report['branching_factor'],.01)
                self.assertFalse(report['associated_branching_factor_applied'])
                self.assertAlmostEqual(report['expected_pb'][0],.001)
                self.assertAlmostEqual(report['nominal_combined_mc_error_pb'],np.hypot(.000003,.00001))
                with self.assertRaisesRegex(ValueError,'overwrite'):
                    compare(native_path,audit,output)
                dm['settings']['iseed']=55001
                with self.assertRaisesRegex(ValueError,'settings'):
                    compare(native_path,audit,root/'invalid.json')


if __name__=='__main__':
    unittest.main()
