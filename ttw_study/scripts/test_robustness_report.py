"""Synthetic report checks; physical collector validation is separate."""
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
from pdf_statistics import uncertainty
from robustness_report import CONFIGS, DEFAULT_SPECTRA, run
from robustness_statistics import LABELS, group_key
from run_robustness_campaign import STATES, comparisons as all_comparisons


def fixture():
    titles, edges, offsets, banks = [], [], [0], {}
    for config in CONFIGS:
        for charge, sign in (('plus', 'W+'), ('minus', 'W-')):
            banks[config, charge] = len(titles)
            for local in range(1, 22):
                length = 9 if local == 1 else 2
                titles.append(config+' '+sign+(' rates: synthetic' if local == 1 else ' synthetic_%d' % local))
                edges.extend([[i, i+1] for i in range(length)])
                offsets.append(offsets[-1]+length)
    layout = dict(titles=titles, edges=np.asarray(edges, dtype=float),
                  offsets=np.asarray(offsets), weights=WEIGHTS)
    comparisons = [dict(row, flavours=[list(f) for f in FLAVOURS[:2]]) for row in all_comparisons()]
    ensembles, conditional = {}, {}
    spread = np.asarray([.94, .97, 1., 1.03, 1.06])
    scale = [1.+.01*np.log(r)+.02*np.log(f)+.03*np.log(t)+.04*np.log(a) for r, f, t, a in POINTS]
    for c, charge in enumerate(('plus', 'minus')):
        for s, (mode, mass) in enumerate(STATES):
            for v, variant in enumerate(('S', 'Pi')):
                for f, flavour in enumerate(FLAVOURS[:2]):
                    amplitude = (1. if f == 0 else 9.)*(1.+.1*s)*(1. if c == 0 else .6)
                    if v:
                        amplitude *= (1.5 if f == 0 else 1.1)*(1.+.025*s)
                    pdf = 1.+(.015+.004*s+.009*f+.012*v+.003*c)*np.arange(101)/100.
                    weight = np.r_[1., scale, pdf]
                    base = np.zeros(len(edges))
                    for cut, config in enumerate(CONFIGS):
                        h = banks[config, charge]
                        start = offsets[h]
                        factor = (1.-.08*cut)*(1.-.015*cut*f)*(1.+.002*s*cut*(1+v))
                        rates = np.asarray([100., 80., 60., 40.*factor, 30.*factor,
                                            20.*factor, 8.*factor, 2.*factor, 0.])
                        if variant == 'S':
                            rates[6:9] = [10.*factor, 0., 0.]
                        base[start:start+9] = rates
                        for local in range(2, 22):
                            a, b = offsets[h+local-1:h+local+1]
                            parent = rates[4] if local in (11, 12, 13, 14) else rates[3]
                            shape = np.asarray([.6, .3] if f == 0 else [.15, .75])
                            if c:
                                shape += [-.2, .2] if f == 0 else [-.1, .1]
                            base[a:b] = parent*shape
                    key = group_key(mode, mass, charge, variant, flavour)
                    ensembles[key] = spread[:, None, None]*amplitude*base[None, :, None]*weight[None, None, :]
                    conditional[key] = np.full((5, len(edges)), .01)
    source = dict(comparisons=comparisons, arrays_sha256='synthetic fixture',
                  benchmark_sha256=digest(STUDY/'inputs/benchmark.json'))
    return source, layout, ensembles, conditional, banks


def spectrum(section, identifier, name):
    record = section['spectra'][identifier][name]
    if digest(record['path']) != record['sha256']:
        raise AssertionError('Spectrum summary hash mismatch')
    detail = json.loads(gzip.decompress(Path(record['path']).read_bytes()))
    if digest(detail['arrays']['path']) != detail['arrays']['sha256']:
        raise AssertionError('Spectrum array hash mismatch')
    return detail


class TestRobustnessReport(unittest.TestCase):
    def test_shared_states_pdf_reduction_parent_shapes_and_streamed_arrays(self):
        source, layout, values, conditional, banks = fixture()
        self.assertIn(10, DEFAULT_SPECTRA)
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_report_') as temporary:
            folder = Path(temporary)
            path, output = folder/'input.json', folder/'report.json'
            path.write_text('{"scope":"synthetic wrapper fixture; physical reader is mocked"}\n')
            binning = folder/'main_bins.json'
            binning.write_text(json.dumps(dict(schema='main_common_binning_v1',
                layout_sha256=layout_signature(layout), histograms={'11': [0., 2.]})))
            real_save = np.savez_compressed
            saved = []
            def checked_save(filename, **arrays):
                for name, value in arrays.items():
                    if 'covariance_factor_nominal' in name or 'nonlinear_bias_nominal' in name:
                        self.assertTrue(value.flags.owndata, name)
                saved.append(str(filename))
                real_save(filename, **arrays)
            with mock.patch('robustness_report.read_vectors', return_value=(source, layout, values, conditional)), \
                 mock.patch('robustness_report.np.savez_compressed', side_effect=checked_save):
                with self.assertRaisesRegex(ValueError, 'all eight'):
                    run(path, output, binning, locals=(10, 11))
                self.assertFalse(folder.joinpath('report_spectra').exists())
                run(path, output, binning, locals=(10, 11), full_flavour=False)
            report = json.loads(output.read_text())
            self.assertEqual(len(report['rates']), 10)
            self.assertEqual(len(report['charge_observables']['pairs']), 5)
            self.assertFalse(report['full_flavour'])
            self.assertEqual(report['binning']['retained_bins'], 480)
            rate = report['rates']['top_W__plus']['R04_b25']['input']
            self.assertAlmostEqual(rate['reference_S']['scale']['central'], 1000.)
            self.assertAlmostEqual(rate['reference_Pi_over_S']['scale']['central'], 1.14)
            self.assertAlmostEqual(rate['product_shift_change']['scale']['central'], 45.35)
            self.assertAlmostEqual(rate['product_ratio_double_ratio']['scale']['central'], 1.025)
            self.assertEqual(len(rate['target_Pi']['scale']['points']), 81)
            self.assertNotIn('scan_S', rate)
            start = int(layout['offsets'][banks['R04_b25', 'plus']])
            members = {}
            for variant in ('S', 'Pi'):
                members[variant] = sum(values[group_key(*STATES[0], 'plus', variant, f)][:, start, 82:].mean(axis=0)
                                       for f in FLAVOURS[:2])
            expected_pdf = uncertainty(members['Pi']/members['S'])
            observed_pdf = rate['reference_Pi_over_S']['pdf']
            self.assertGreater(observed_pdf['error_symmetric'], 0.)
            self.assertAlmostEqual(observed_pdf['error_symmetric'], expected_pdf['error_symmetric'])
            migration = report['migrations']['top_W__plus']['R04_b30/fiducial_1b_over_R04_b25']
            self.assertAlmostEqual(migration['reference_S']['scale']['central'], (.92+9.*.9062)/10.)
            leading = spectrum(report, 'top_W__plus', 'R04_b25_h10')
            subleading = spectrum(report, 'top_W__plus', 'R04_b25_h11')
            self.assertEqual(leading['normalization_rate_bins'], ['fiducial_1b'])
            self.assertEqual(subleading['normalization_rate_bins'], ['fiducial_2b'])
            self.assertEqual(len(subleading['normalized']), 1)
            self.assertAlmostEqual(leading['normalized'][0]['reference_S']['scale']['central'], .195)
            self.assertAlmostEqual(subleading['normalized'][0]['reference_S']['scale']['central'], .9)
            self.assertAlmostEqual(leading['measured_range_coverage']['reference_S']['scale']['central'], .9)
            charge = report['charge_observables']
            rate = charge['rates']['top_W']['R04_b25']['input']
            self.assertAlmostEqual(rate['sum']['reference_S']['scale']['central'], 1600.)
            self.assertAlmostEqual(rate['plus_over_minus']['reference_S']['scale']['central'], 5./3.)
            self.assertAlmostEqual(rate['asymmetry']['reference_S']['scale']['central'], .25)
            combined = spectrum(charge, 'top_W', 'R04_b25_h10')
            self.assertAlmostEqual(combined['normalized'][0]['combined']['reference_S']['scale']['central'], .15375)
            self.assertNotAlmostEqual(.15375, .5*(.195+.085))
            covariance = report['shared_nominal_covariance']['joint_all_rates_nominal']
            self.assertEqual(covariance['arrays']['sha256'], digest(covariance['arrays']['path']))
            with np.load(covariance['arrays']['path']) as arrays:
                factor = arrays[covariance['factor_array']]
                delta = LABELS.index('product_shift_change')
                top = next(i for i, row in enumerate(source['comparisons']) if row['name'] == 'top_W' and row['charge'] == 'plus')
                associated = next(i for i, row in enumerate(source['comparisons']) if row['name'] == 'associated_W' and row['charge'] == 'plus')
                total = next(i for i, row in enumerate(source['comparisons']) if row['name'] == 'all_W' and row['charge'] == 'plus')
                expected_variance = sum(np.var(row[:, start, 0], ddof=1)/len(row)
                    for key, row in values.items() if key[:2] == STATES[1] and key[2] == 'plus')
                self.assertAlmostEqual(factor[:, top, delta, 0]@factor[:, associated, delta, 0], -expected_variance, places=9)
                np.testing.assert_allclose(factor[:, total, delta],
                    factor[:, top, delta]+factor[:, associated, delta], atol=2.e-10)
                self.assertEqual(arrays['c0_R04_b25_rates'].shape, (13, 9, 183))
                self.assertFalse(any('_h10_' in name or '_h11_' in name for name in arrays.files))
            with np.load(leading['arrays']['path']) as arrays:
                self.assertEqual(arrays['normalized'].shape, (13, 2, 183))
                self.assertEqual(arrays['weights'].tolist(), WEIGHTS)
            with np.load(combined['arrays']['path']) as arrays:
                self.assertEqual(arrays['normalized'].shape, (3, 9, 2, 183))
            joint = report['shared_nominal_covariance']['joint_R04_b25_h10_normalized_nominal']
            self.assertNotEqual(joint['arrays']['path'], covariance['arrays']['path'])
            self.assertEqual(joint['arrays']['sha256'], digest(joint['arrays']['path']))
            with np.load(joint['arrays']['path']) as arrays:
                self.assertEqual(list(arrays[joint['factor_array']].shape), joint['shape'])
            self.assertGreater(len(saved), 1)
            diagnostic = report['retraining_diagnostics']['onshell__0.0__plus__S__'+''.join(FLAVOURS[0])]
            self.assertGreater(diagnostic['variance_ratio'][0], 1.)
            self.assertEqual(report['source_hashes']['robustness_report.py'], digest(STUDY/'scripts/robustness_report.py'))
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                run(path, output)

    def test_wrong_weights_or_charge_bins_fail_before_output(self):
        source, layout, values, conditional, banks = fixture()
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_guards_') as temporary:
            folder = Path(temporary)
            path, output = folder/'input.json', folder/'report.json'
            path.write_text('{"scope":"synthetic rejection fixture"}\n')
            wrong = dict(layout, weights=WEIGHTS[:82])
            with mock.patch('robustness_report.read_vectors', return_value=(source, wrong, values, conditional)):
                with self.assertRaisesRegex(ValueError, '183-weight'):
                    run(path, output, locals=(10, 11), full_flavour=False)
            wrong = dict(layout, edges=layout['edges'].copy())
            h = banks['R04_b25', 'minus']+9
            wrong['edges'][layout['offsets'][h], 0] += .1
            with mock.patch('robustness_report.read_vectors', return_value=(source, wrong, values, conditional)):
                with self.assertRaisesRegex(ValueError, 'different bin boundaries'):
                    run(path, output, locals=(10, 11), full_flavour=False)
            self.assertFalse(output.exists())
            self.assertFalse(folder.joinpath('report_spectra').exists())


if __name__ == '__main__':
    unittest.main()
