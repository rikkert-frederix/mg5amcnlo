"""Complete synthetic parameter-report and main-bin compatibility tests."""
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from campaign import STUDY, digest
from full_flavour_statistics import FLAVOURS
from load_results import POINTS, WEIGHTS
from main_binning import layout_signature
from parameter_charge_statistics import LABELS as CHARGE_CONTRAST_LABELS
from parameter_statistics import LABELS, comparison_keys, select
from parameter_report import CONFIGS, apply_common_binning, run


class TestParameterReport(unittest.TestCase):
    def fixture(self):
        titles, edges, offsets, banks = [], [], [0], {}
        for config in CONFIGS:
            for charge, sign in (('plus', 'W+'), ('minus', 'W-')):
                banks[config, charge] = len(titles)
                for local in range(1, 22):
                    length = 9 if local == 1 else 2
                    titles.append(config+' '+sign+(' rates: synthetic' if local == 1 else ' synthetic_%d' % local))
                    edges.extend([[i, i+1] for i in range(length)])
                    offsets.append(offsets[-1]+length)
        layout = dict(titles=titles, edges=np.asarray(edges), offsets=np.asarray(offsets), weights=WEIGHTS[:82])
        comparisons = [dict(scenario=scenario, w_treatment='onshell', charge=charge,
                            flavours=[list(f) for f in FLAVOURS[:2]])
                       for scenario in ('as117', 'as119') for charge in ('plus', 'minus')]
        ensembles, conditional = {}, {}
        spread = np.asarray([.94, .97, 1., 1.03, 1.06])
        weight = np.r_[1., [1.+.01*np.log(r)+.02*np.log(f)+.03*np.log(t)+.04*np.log(a) for r,f,t,a in POINTS]]
        for comparison in comparisons:
            for (origin, variant), keys in comparison_keys(comparison, False).items():
                for flavour_index, key in enumerate(keys):
                    if key in ensembles: continue
                    amplitude = 1. if flavour_index == 0 else 9.
                    if variant == 'Pi': amplitude *= 1.5 if flavour_index == 0 else 1.1
                    if origin == 'scan':
                        amplitude *= 1.05 if comparison['scenario'] == 'as117' else .95
                        if variant == 'Pi': amplitude *= 1.02 if comparison['scenario'] == 'as117' else .98
                    if comparison['charge'] == 'minus': amplitude *= .6
                    base = np.zeros(len(edges))
                    for cut, config in enumerate(CONFIGS):
                        h = banks[config, comparison['charge']]
                        start = offsets[h]
                        factor = 1.-.1*cut
                        rates = np.asarray([100., 80., 60., 40.*factor, 30.*factor,
                                            20.*factor, 8.*factor, 2.*factor, 0.])
                        if variant == 'S': rates[6:9] = [10.*factor, 0., 0.]
                        base[start:start+9] = rates
                        for local in range(2, 22):
                            a,b = offsets[h+local-1:h+local+1]
                            parent = rates[4] if local in (11,12,13,14) else rates[3]
                            shape = np.asarray([.6,.3]) if flavour_index == 0 else np.asarray([.15,.75])
                            base[a:b] = parent*shape
                    ensembles[key] = spread[:,None,None]*amplitude*base[None,:,None]*weight[None,None,:]
                    conditional[key] = np.full((5,len(edges)), .01)
        source = dict(comparisons=comparisons, arrays_sha256='synthetic fixture',
                      benchmark_sha256=digest(STUDY/'inputs/benchmark.json'))
        return source, layout, ensembles, conditional

    def test_main_bins_apply_without_inventing_pdf_weights_or_merging_categories(self):
        source, layout, values, conditional = self.fixture()
        plan = dict(schema='main_common_binning_v1',
                    layout_sha256=layout_signature(dict(layout, weights=WEIGHTS)), histograms={'11':[0.,2.]})
        updated, rebinned, diagnostics, metadata = apply_common_binning(layout, values, conditional, plan)
        self.assertEqual(updated['weights'], WEIGHTS[:82])
        self.assertEqual(len(updated['edges']), 480)
        self.assertEqual(len(metadata['changed_histograms']), 10)
        for bank in range(10):
            for local in (2,14):
                a,b = updated['offsets'][bank*21+local-1:bank*21+local+1]
                self.assertEqual(b-a, 2)
        key = next(iter(values))
        a,b = updated['offsets'][10:12]
        np.testing.assert_allclose(rebinned[key][:,a:b], values[key][:,layout['offsets'][10]:layout['offsets'][11]].sum(axis=1)[:,None])
        self.assertTrue(np.isnan(diagnostics[key][:,a:b]).all())
        wrong = dict(plan, histograms={'14':[0.,2.]})
        with self.assertRaises(ValueError): apply_common_binning(layout, values, conditional, wrong)
        wrong = dict(plan, layout_sha256='different geometry')
        with self.assertRaises(ValueError): apply_common_binning(layout, values, conditional, wrong)

    def test_complete_selected_report_preserves_shared_covariance_and_normalization(self):
        source, layout, values, conditional = self.fixture()
        with tempfile.TemporaryDirectory(prefix='synthetic_parameter_report_') as folder:
            folder = Path(folder)
            path = folder/'synthetic_input.json'
            path.write_text('{"data_scope":"synthetic report fixture only"}\n')
            binning = folder/'synthetic_main_bins.json'
            binning.write_text(json.dumps(dict(schema='main_common_binning_v1',
                layout_sha256=layout_signature(dict(layout,weights=WEIGHTS)), histograms={'11':[0.,2.]})))
            output = folder/'synthetic_report.json'
            actual_save = np.savez_compressed
            def checked_save(filename, **arrays):
                for name, value in arrays.items():
                    if 'covariance_factor_nominal' in name or 'nonlinear_bias_nominal' in name:
                        self.assertTrue(value.flags.owndata, name)
                actual_save(filename, **arrays)
            with mock.patch('parameter_report.read_vectors', return_value=(source,layout,values,conditional)), \
                 mock.patch('parameter_report.np.savez_compressed', side_effect=checked_save):
                # Explicit subset labels are compulsory for this two-flavour test.
                with self.assertRaisesRegex(ValueError, 'all eight'):
                    run(path,output,binning,locals=(11,14))
                self.assertFalse(output.parent.joinpath(output.stem+'_spectra').exists())
                run(path,output,binning,locals=(11,14),full_flavour=False)
            report = json.loads(output.read_text())
            self.assertEqual(len(report['rates']),4)
            self.assertFalse(report['full_flavour'])
            self.assertEqual(report['binning']['retained_bins'],480)
            rate = report['rates']['as117__onshell__plus']['R04_b25']['input']
            self.assertAlmostEqual(rate['reference_S']['scale']['central'],1000.)
            self.assertAlmostEqual(rate['reference_Pi_over_S']['scale']['central'],1.14)
            self.assertAlmostEqual(rate['product_ratio_change']['scale']['central'],.0228)
            self.assertAlmostEqual(rate['product_ratio_double_ratio']['scale']['central'],1.02)
            self.assertNotIn('pdf',rate['scan_S'])
            minus = report['rates']['as117__onshell__minus']['R04_b25']['input']
            self.assertAlmostEqual(minus['reference_S']['scale']['central'],600.)
            files = report['spectra']['as117__onshell__plus']
            spectrum = json.loads(gzip.decompress(Path(files['R04_b25_h11']['path']).read_bytes()))
            self.assertEqual(len(spectrum['normalized']),1)
            self.assertEqual(spectrum['normalization_rate_bins'],['fiducial_2b'])
            self.assertAlmostEqual(spectrum['normalized'][0]['reference_S']['scale']['central'],.9)
            self.assertAlmostEqual(spectrum['normalized'][0]['scan_Pi']['scale']['central'],.9)
            self.assertAlmostEqual(spectrum['normalized'][0]['product_ratio_change']['scale']['central'],0.)
            self.assertEqual(len(spectrum['normalized'][0]['scan_Pi']['scale']['points']),81)
            categorical = json.loads(gzip.decompress(Path(files['R04_b25_h14']['path']).read_bytes()))
            self.assertEqual(len(categorical['normalized']),2)
            charge = report['charge_observables']
            self.assertEqual(len(charge['pairs']),2)
            self.assertEqual(charge['unpaired'],[])
            rate = charge['rates']['as117__onshell']['R04_b25']['input']
            self.assertAlmostEqual(rate['sum']['reference_S']['scale']['central'],1600.)
            self.assertAlmostEqual(rate['sum']['product_shift_change']['scale']['central'],49.504)
            self.assertAlmostEqual(rate['plus_over_minus']['reference_S']['scale']['central'],5./3.)
            self.assertAlmostEqual(rate['asymmetry']['reference_S']['scale']['central'],.25)
            self.assertAlmostEqual(rate['asymmetry']['product_shift_change']['scale']['central'],0.)
            self.assertGreater(rate['asymmetry']['product_shift_change']['nominal_retraining_mc_error'],0.)
            self.assertEqual(len(rate['asymmetry']['scan_Pi']['scale']['points']),81)
            self.assertNotIn('pdf',rate['sum']['scan_S'])
            target = charge['spectra']['as117__onshell']['R04_b25_h11']['path']
            spectrum = json.loads(gzip.decompress(Path(target).read_bytes()))
            self.assertAlmostEqual(spectrum['normalized'][0]['combined']['reference_S']['scale']['central'],.9)
            self.assertAlmostEqual(spectrum['normalized'][0]['plus_shape_over_minus_shape']['reference_S']['scale']['central'],1.)
            self.assertAlmostEqual(spectrum['normalized'][0]['shape_asymmetry']['reference_S']['scale']['central'],0.)
            self.assertEqual(report['source_hashes']['parameter_charge_statistics.py'],digest(STUDY/'scripts/parameter_charge_statistics.py'))
            with np.load(output.with_suffix('.npz')) as arrays:
                factor = arrays['joint_all_rates_nominal_covariance_factor_nominal']
                label = LABELS.index('product_shift_change')
                self.assertGreater(factor[:,0,label,0]@factor[:,2,label,0],0.)
                self.assertAlmostEqual(factor[:,0,label,0]@factor[:,1,label,0],0.,places=10)
                self.assertEqual(arrays['c0_R04_b25_h11_normalized'].shape,(13,1,82))
                factor = arrays['joint_charge_all_rates_nominal_covariance_factor_nominal']
                delta = CHARGE_CONTRAST_LABELS.index('product_shift_change')
                base = CHARGE_CONTRAST_LABELS.index('reference_Pi_minus_S')
                for q in (0,1,2):
                    self.assertAlmostEqual(factor[:,0,q,delta,0]@factor[:,1,q,delta,0],
                        arrays['joint_charge_all_rates_nominal_mc_errors'][0,q,base,0,0]**2,places=9)
                self.assertEqual(arrays['charge_c0_R04_b25_h11_normalized'].shape,(3,9,1,82))
            with self.assertRaises(ValueError): run(path,output)

    def test_unpaired_selection_is_explicit_and_mismatched_charge_bins_are_rejected(self):
        source, layout, values, conditional = self.fixture()
        with tempfile.TemporaryDirectory(prefix='synthetic_parameter_charge_scope_') as folder:
            folder = Path(folder)
            path = folder/'synthetic_input.json'
            path.write_text('{"data_scope":"synthetic charge scope fixture only"}\n')
            output = folder/'synthetic_rejected.json'
            wrong = dict(layout,edges=layout['edges'].astype(float))
            wrong['edges'][layout['offsets'][21+10],0] += .1
            with mock.patch('parameter_report.read_vectors',return_value=(source,wrong,values,conditional)):
                with self.assertRaisesRegex(ValueError,'different bin boundaries'):
                    run(path,output,locals=(11,),full_flavour=False)
            self.assertFalse(output.parent.joinpath(output.stem+'_spectra').exists())
            source = dict(source,comparisons=[row for row in source['comparisons'] if row['charge']=='plus'])
            values = select(values,source['comparisons'],False)
            conditional = {key:conditional[key] for key in values}
            output = folder/'synthetic_unpaired.json'
            with mock.patch('parameter_report.read_vectors',return_value=(source,layout,values,conditional)):
                run(path,output,locals=(14,),full_flavour=False)
            report = json.loads(output.read_text())
            self.assertEqual(len(report['rates']),2)
            self.assertEqual(report['charge_observables']['pairs'],[])
            self.assertEqual(len(report['charge_observables']['unpaired']),2)
            self.assertFalse(report['charge_observables']['rates'])
            with np.load(output.with_suffix('.npz')) as arrays:
                self.assertFalse(any(key.startswith(('joint_charge_','charge_c')) for key in arrays
                                     if key != 'charge_contrast_labels'))


if __name__=='__main__':
    unittest.main()
