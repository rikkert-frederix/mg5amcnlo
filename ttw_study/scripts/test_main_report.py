"""Synthetic end-to-end reporting test; never creates study physics results."""
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from full_flavour_statistics import FLAVOURS,VARIANTS
from joint_report import normalization_bins
from load_results import WEIGHTS
from main_binning import CONTINUOUS_IDS, layout_signature
from main_report import run


class TestMainReport(unittest.TestCase):
    def fixture(self,continuous_bins=1):
        configs=('R04_b25','R04_b30','R04_b40','R03_b25','R05_b25')
        titles,edges,offsets,locations=[],[],[0],{}
        for cut,config in enumerate(configs):
            for charge,sign in (('plus','W+'),('minus','W-')):
                rates=np.asarray([100.,80.,60.,40.*(1-.05*cut),30.*(1-.1*cut),
                                  15.*(1-.1*cut),8.*(1-.1*cut),5.*(1-.1*cut),2.*(1-.1*cut)])
                start=offsets[-1]
                locations[config,charge]=(start,rates,len(titles))
                for local in range(1,22):
                    length=9 if local==1 else continuous_bins
                    titles.append(config+' '+sign+(' rates: synthetic' if local==1 else ' observable_%d'%local))
                    edges.extend([[i,i+1] for i in range(length)])
                    offsets.append(offsets[-1]+length)
        layout=dict(titles=titles,offsets=np.array(offsets),edges=np.array(edges),weights=WEIGHTS)
        vectors,conditional={},{}
        amplitudes=dict(LO=1.,P=1.2,D=.9,S=1.1,Pi=1.13,PiD=.915)
        variation=np.ones(183)
        variation[82:]=np.linspace(.99,1.01,101)
        variation[82]=1.
        spread=np.asarray([.96,.98,1.,1.02,1.04])
        for charge,charge_factor in (('plus',1.),('minus',.6)):
            base=np.zeros(len(edges))
            for config in configs:
                start,rates,h=locations[config,charge]
                base[start:start+9]=rates
                for local in range(2,22):
                    a,b=layout['offsets'][h+local-1:h+local+1]
                    base[a:b]=.6*sum(rates[i] for i in normalization_bins(local))/continuous_bins
            for variant in VARIANTS:
                for i,flavour in enumerate(FLAVOURS):
                    key=charge,variant,flavour
                    vectors[key]=(spread[:,None,None]*base[None,:,None]*variation[None,None,:]*
                                  amplitudes[variant]*charge_factor*(i+1)/36)
                    conditional[key]=np.full((5,len(edges)),.001)
        source=dict(arrays_sha256='synthetic-fixture',benchmark_sha256='synthetic-fixture')
        return source,layout,vectors,conditional

    def test_all_outputs_keep_physical_ratio_order_and_covariance(self):
        self.check_complete_report(False)

    def test_common_bin_plan_keeps_charge_normalization_and_rate_diagnostics(self):
        self.check_complete_report(True)

    def check_complete_report(self,merge):
        source,layout,vectors,conditional=self.fixture(2 if merge else 1)
        with tempfile.TemporaryDirectory(prefix='synthetic_main_report_') as folder:
            path=Path(folder)/'synthetic_input.json'
            path.write_text('{}\n')
            output=Path(folder)/'synthetic_report.json'
            plan_path=None
            if merge:
                plan_path=Path(folder)/'synthetic_bins.json'
                plan_path.write_text(json.dumps(dict(schema='main_common_binning_v1',
                    layout_sha256=layout_signature(layout),histograms={str(local):[0,2] for local in CONTINUOUS_IDS})))
            with mock.patch('main_report.read_main',return_value=(source,layout,vectors,conditional)):
                run(path,output,plan_path)
            result=json.loads(output.read_text())
            self.assertEqual(len(result['rates']),10)
            self.assertEqual(len(result['spectra']),200)
            self.assertEqual(len(result['charge_spectra']),100)
            rate=result['rates']['R04_b25 W+ rates: synthetic']['input']
            self.assertAlmostEqual(rate['S']['scale']['central'],110.)
            self.assertAlmostEqual(rate['S_minus_P_minus_D_plus_LO']['scale']['central'],0.)
            self.assertGreater(rate['S_minus_P_minus_D_plus_LO']['nominal_retraining_mc_error'],0.)
            charge=result['charge_rates']['R04_b25 W+/W- rates: synthetic']['input']['S']
            self.assertAlmostEqual(charge['plus_over_minus']['scale']['central'],1/.6)
            self.assertAlmostEqual(charge['asymmetry']['scale']['central'],.25)
            h='R04_b25 W+ observable_3'
            detail=json.loads(gzip.decompress(Path(result['spectra'][h]['path']).read_bytes()))
            self.assertAlmostEqual(detail['normalized'][0]['S']['scale']['central'],.6)
            self.assertAlmostEqual(detail['measured_range_coverage']['S']['scale']['central'],.6)
            self.assertNotIn('S_minus_P_minus_D_plus_LO',detail['normalized'][0])
            self.assertEqual(len(detail['normalized'][0]['S']['scale']['points']),81)
            ch='R04_b25 W+/W- observable_3'
            detail=json.loads(gzip.decompress(Path(result['charge_spectra'][ch]['path']).read_bytes()))
            self.assertAlmostEqual(detail['normalized'][0]['S']['combined']['scale']['central'],.6)
            self.assertAlmostEqual(detail['normalized'][0]['S']['plus_shape_over_minus_shape']['scale']['central'],1.)
            with np.load(output.with_suffix('.npz')) as arrays:
                expected_bins=310 if merge else 290
                self.assertEqual(arrays['absolute'].shape,(12,expected_bins,183))
                self.assertEqual(arrays['normalized'].shape,(11,expected_bins,183))
                self.assertIn('rates_0_covariance_factor_nominal',arrays)
                if merge:
                    for bank in range(10):
                        for local in (2,14):
                            a,b=arrays['offsets'][bank*21+local-1:bank*21+local+1]
                            self.assertEqual(b-a,2)
            if merge:
                self.assertEqual(result['binning']['original_bins'],490)
                self.assertEqual(result['binning']['retained_bins'],310)
                self.assertEqual(len(result['binning']['changed_histograms']),180)
                diagnostic=result['retraining_diagnostics']['plus__S__eee']
                np.testing.assert_allclose(diagnostic['mean_conditional_variance_of_combination'],.0002)
            else:
                self.assertIsNone(result['binning'])
            with self.assertRaises(ValueError):run(path,output)


if __name__=='__main__':
    unittest.main()
