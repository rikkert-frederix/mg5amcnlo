"""Tests of correlated scale operations, HwU inputs, and run preparation."""
import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.various.banner import RunCardNLO


def module(name):
    path = Path(MG5DIR) / 'Template/fNLO/FixedOrderAnalysis' / (name+'.py')
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


scales = module('ttw_product_scales')
setup = module('ttw_product_setup')


class TestTTWProductTools(unittest.TestCase):
    def test_correlated_scale_ratios_and_groups(self):
        points = sorted(scales.POINTS)
        s = {k: k[0]/k[1]+k[2]+2.*k[3] for k in points}
        s['central'] = 4.
        p = {k: 2.*v for k, v in s.items()}
        report = scales.compare_values(s, p)
        self.assertEqual(report['strict']['log_responses']['production']['absolute'], 0.)
        self.assertGreater(report['strict']['log_responses']['decay']['absolute'], 0.)
        self.assertEqual(report['strict']['log_responses']['antitop']['absolute'],
                         2.*report['strict']['log_responses']['top']['absolute'])
        for group, count in [('production7', 7), ('top3', 3), ('antitop3', 3),
                             ('decay3', 3), ('decay9', 9), ('combined21', 21),
                             ('combined63', 63), ('common3', 3), ('all81_diagnostic', 81)]:
            band = report['product_over_strict'][group]
            self.assertEqual(band['npoints'], count)
            self.assertEqual(band['envelope'], [2., 2.])
        shared_s = {k[:3]: v for k,v in s.items() if k != 'central' and k[2] == k[3]}
        shared_s['central'] = s['central']
        shared_p = {k: 2.*v for k,v in shared_s.items()}
        shared = scales.compare_values(shared_s, shared_p)
        self.assertEqual(shared['strict']['combined21'], report['strict']['combined21'])
        self.assertEqual(shared['strict']['all27_diagnostic']['npoints'], 27)
        with self.assertRaisesRegex(ValueError, 'grids differ'):
            scales.compare_values(shared_s, p)
        p[(.5, 2., .5, 2.)] = 999.
        report = scales.compare_values(s, p)
        self.assertEqual(report['product_over_strict']['combined63']['envelope'], [2., 2.])
        self.assertGreater(report['product_over_strict']['all81_diagnostic']['envelope'][1], 100.)
        p[(1., 1., .5, 2.)] = 100.
        report = scales.compare_values(s, p)
        self.assertEqual(report['product_over_strict']['combined21']['envelope'], [2., 2.])
        self.assertGreater(report['product_over_strict']['combined63']['envelope'][1], 2.)
        s[(1., 1., .5, .5)] = 0.
        report = scales.compare_values(s, p)
        self.assertIsNone(report['product_over_strict']['decay3']['envelope'])
        self.assertEqual(report['product_over_strict']['decay3']['nvalid'], 2)

    def test_signed_label_parsing(self):
        for label in ['dyn=3 muR=1 muF=.5 d-6=2 d6=1',
                      'dyn=3 muR=1 muF=.5 d6=1 d-6=2D0']:
            self.assertEqual(scales.scale_point(label), (1., .5, 1., 2.))
        self.assertEqual(scales.scale_point('dyn=3 muR=1 muF=.5 d6=2'), (1., .5, 2.))
        for suffix in ['d6=1 d6=2', 'd6=1 d24=2', 'd-6=1']:
            with self.assertRaisesRegex(ValueError, 'unique'):
                scales.scale_point('dyn=3 muR=1 muF=1 '+suffix)

    @staticmethod
    def fixture(path, multiplier=1., missing_point=False, shared=False):
        points = sorted(scales.SHARED_POINTS if shared else scales.POINTS)
        if missing_point:
            points.pop()
        lines = ['##& xmin & xmax & central value & dy & ' + ' & '.join(
            ('dyn=3 muR=%g muF=%g d6=%g' % point if shared else
             'dyn=3 muR=%g muF=%g d-6=%g d6=%g' % (point[0],point[1],point[3],point[2]))
            for point in points), '']
        for config in ['R04_b25','R03_b25','R05_b25','R04_b30','R04_b40']:
            for charge in ['W+','W-']:
                for i in range(21):
                    title = '%s %s %s' % (config, charge,
                            'rates: stages' if i==0 else '1b test%d' % i)
                    base = [100.,80.,70.,20.,10.,5.,3.,2.,0.] if i==0 else [4.]
                    lines.append('<histogram> %d "%s"' % (len(base), title))
                    for j, value in enumerate(base):
                        value *= multiplier * (1. if charge=='W+' else .5)
                        varied = [value * (k[0]/k[1]+sum(k[2:]))/(len(k)-1) for k in points]
                        lines.append(' '.join('%+.9e' % v for v in
                                     [j+.5,j+1.5,value,abs(value)*.01]+varied))
                    lines.extend(['<\\histogram>', ''])
        path.write_text('\n'.join(lines))

    def test_hwu_normalization_charge_and_missing_points(self):
        with tempfile.TemporaryDirectory(prefix='ttw_scales_') as directory:
            directory = Path(directory)
            s, p, bad = (directory / name for name in ['s.HwU', 'p.HwU', 'bad.HwU'])
            self.fixture(s)
            self.fixture(p, multiplier=2.)
            self.fixture(bad, missing_point=True)
            strict, product = scales.load_sum([s]), scales.load_sum([p])
            report = scales.make_report(strict, product)
            self.assertEqual(report['scale_point_count'], 81)
            self.assertEqual(report['scale_axes'][-2:], ['muR_top', 'muR_antitop'])
            hist = report['histograms']['R04_b25 W+ 1b test1'][0]
            self.assertAlmostEqual(hist['absolute']['product_over_strict']['central'], 2.)
            self.assertAlmostEqual(hist['normalized_to_fiducial']['product_over_strict']['central'], 1.)
            derived = report['derived']['R04_b25 W+']
            self.assertAlmostEqual(derived['acceptance_1b']['strict']['central'], .2)
            self.assertAlmostEqual(derived['charge_ratio_fiducial_1b']['strict']['central'], 2.)
            self.assertAlmostEqual(derived['charge_asymmetry_fiducial_1b']['strict']['central'], 1./3.)
            json.dumps(report, allow_nan=False)
            with self.assertRaisesRegex(ValueError, '27'):
                scales.load_sum([bad])
            with self.assertRaisesRegex(ValueError, 'twice'):
                scales.load_sum([s, s])
            summed = scales.load_sum([s, p])
            self.assertEqual(summed['R04_b25 W+ 1b test1'][0]['values']['central'], 12.)
            self.fixture(bad, shared=True)
            shared = scales.load_sum([bad])
            self.assertEqual(scales.make_report(shared, shared)['scale_point_count'], 27)
            with self.assertRaisesRegex(ValueError, 'grids differ'):
                scales.load_sum([s, bad])

    def test_process_commands(self):
        text = setup.process_commands('minus', ['mu','e','mu'], '/tmp/TTW_test')
        self.assertIn('import model loop_sm-no_b_mass', text)
        self.assertIn('t t~ w- QCD=2 QED=1 [QCD]', text)
        self.assertIn('w+ > mu+ vm', text)
        self.assertIn('w- > mu- vm~', text)
        self.assertNotIn(' -f', text)
        with self.assertRaises(ValueError):
            setup.process_commands('plus', ['e','e','mu'], '/tmp/a; command')

    def test_configuration_records_widths_and_all_scales(self):
        with tempfile.TemporaryDirectory(prefix='ttw_setup_') as directory:
            process = Path(directory)
            cards = process / 'Cards'
            cards.mkdir()
            (process / 'SubProcesses').mkdir()
            (process / 'SubProcesses/decay_chain_parameters.f90').write_text('DECAY_SCALE_GROUPING')
            (process / 'FixedOrderAnalysis').mkdir()
            (process / 'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90').write_text('fixture')
            template = Path(MG5DIR) / 'Template/fNLO/Cards/run_card.dat'
            RunCardNLO().write(str(cards / 'run_card.dat'), template=str(template))
            (cards / 'run_card_default.dat').write_bytes((cards / 'run_card.dat').read_bytes())
            (cards / 'decay_card.dat').write_text('original decay card\n')
            (cards / 'param_card.dat').write_text('''BLOCK MASS
  5 0.0
  6 173.0
  24 80.4
DECAY 6 1.5
DECAY 24 2.05
''')
            args = argparse.Namespace(process_dir=str(process), variant='S', top_width_lo=1.5,
                                      top_width_nlo=1.35, width_source='synthetic test only',
                                      pdf_id=None, production_scale='core-ht-half', ecm=13000.,
                                      seed=42, points=10, grid_points=10, iterations=1,
                                      decay_scales='separate')
            with contextlib.redirect_stdout(io.StringIO()):
                setup.configure(args)
            decay = (cards / 'decay_card.dat').read_text()
            self.assertIn('INDEPENDENT = decay_scale_variation_mode', decay)
            self.assertIn('SIGNED_PDG = decay_scale_grouping', decay)
            self.assertIn('1, 0.5, 2 = decay_scale_factors', decay)
            self.assertIn('AUTO = decay_width_scale_mode(6)', decay)
            run = RunCardNLO(str(cards / 'run_card.dat'))
            self.assertEqual(list(run['rw_rscale']), [1., .5, 2.])
            self.assertEqual(list(run['rw_fscale']), [1., .5, 2.])
            self.assertFalse(run['cut_decays'])
            self.assertEqual(run['ptj'], 0.)
            archive = process / 'study_cards/S_core-ht-half_separate_42'
            manifest = json.loads((archive / 'manifest.json').read_text())
            self.assertEqual(manifest['top_width_reference_scale'], 173.)
            self.assertEqual(manifest['decay_scale_grouping'], 'separate')
            self.assertEqual((archive / 'before/decay_card.dat').read_text(), 'original decay card\n')
            with self.assertRaisesRegex(ValueError, 'already exists'):
                setup.configure(args)
            args.decay_scales = 'shared'
            with contextlib.redirect_stdout(io.StringIO()):
                setup.configure(args)
            self.assertIn('SPECIES = decay_scale_grouping', (cards / 'decay_card.dat').read_text())
            self.assertTrue((process / 'study_cards/S_core-ht-half_shared_42/manifest.json').is_file())
            (process / 'SubProcesses/decay_chain_parameters.f90').write_text('old export')
            with self.assertRaisesRegex(ValueError, 'Re-export'):
                setup.configure(args)


if __name__ == '__main__':
    unittest.main()
