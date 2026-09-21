"""Tests of correlated scale operations, HwU inputs, and run preparation."""
import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

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

    def test_bw_process_commands(self):
        for charge, assoc in [('plus', 'mu+ vm'), ('minus', 'mu- vm~')]:
            for treatment in setup.W_TREATMENTS:
                for corrected in ('both', 't', 'tbar', 'neither'):
                    text = setup.process_commands(charge, ['e', 'mu', 'mu'], '/tmp/TTW_bw',
                                                  corrected, True, True, treatment)
                    self.assertIn('# ttW W treatment: ' + treatment, text)
                    self.assertEqual(text.count('[real=QCD]'),
                                     1 + (corrected in ('both', 't')) +
                                     (corrected in ('both', 'tbar')))
                    if treatment == 'onshell':
                        self.assertIn('t > w+ b QED=1', text)
                    else:
                        self.assertIn('t > b e+ ve QED=2', text)
                        self.assertIn('t~ > b~ mu- vm~ QED=2', text)
                        self.assertNotIn('t > w+', text)
                    if treatment == 'all-bw':
                        self.assertIn('t t~ %s QCD=2 QED=2' % assoc, text)
                        self.assertNotIn('w+ >', text)
                        self.assertNotIn('w- >', text)
                    else:
                        self.assertIn('QCD=2 QED=1', text)
                        self.assertIn(' > ' + assoc, text)
        with self.assertRaisesRegex(ValueError, 'Unknown W'):
            setup.process_commands('plus', ['e', 'e', 'mu'], '/tmp/TTW', w_treatment='bad')

    def test_massive_decay_commands_and_export_guard(self):
        from models.check_param_card import ParamCard
        for treatment in setup.W_TREATMENTS:
            text = setup.process_commands('plus', ['e', 'e', 'mu'], '/tmp/TTW',
                                          w_treatment=treatment, decay_bottom_mass=4.8)
            self.assertIn('import model loop_sm-no_b_mass\nset decay_bottom_mass 4.8', text)
            self.assertIn('define p = g u c d s b u~ c~ d~ s~ b~', text)
        for mass in (-1., float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'bottom mass'):
                setup.process_commands('plus', ['e', 'e', 'mu'], '/tmp/TTW', decay_bottom_mass=mass)
        with tempfile.TemporaryDirectory() as directory:
            process = Path(directory)
            cards = process / 'Cards'
            cards.mkdir()
            path = cards / 'param_card.dat'
            path.write_text('BLOCK MASS\n 5 0.\nBLOCK DECAYMASS\n 5 4.8\n')
            param = ParamCard(str(path))
            with self.assertRaisesRegex(ValueError, 'without a generated'):
                setup.export_decay_bottom_mass(process, param)
            scheme = cards / 'decay_mass_scheme.json'
            scheme.write_text(json.dumps(dict(format=1, production_model='loop_sm-no_b_mass',
                                              decay_model='loop_sm', parameter=['decaymass', 5],
                                              alpha_s_flavours=5)))
            with self.assertRaisesRegex(ValueError, 'mass lookup'):
                setup.export_decay_bottom_mass(process, param)
            source = process / 'Source/MODEL/get_mass_width_fcts.f'
            source.parent.mkdir(parents=True)
            source.write_text('GET_DECAY_MASS_FROM_ID=DC_MDL_MB')
            self.assertEqual(setup.export_decay_bottom_mass(process, param), (4.8, scheme))
            param['decaymass'].get((5,)).value = 0.
            with self.assertRaisesRegex(ValueError, 'requires positive'):
                setup.export_decay_bottom_mass(process, param)

    @staticmethod
    def topology_fixture(process, treatment):
        """Small format-4 fixtures with real node/leaf/core semantics."""
        directory = process / 'SubProcesses/P0_test'
        directory.mkdir(parents=True, exist_ok=True)
        if treatment == 'onshell':
            nodes = [(1, 0, 6, 'NODE 2 LEAF 3'), (2, 1, 24, 'LEAF 1 LEAF 2'),
                     (3, 0, -6, 'NODE 4 LEAF 6'), (4, 3, -24, 'LEAF 4 LEAF 5'),
                     (5, 0, 24, 'LEAF 7 LEAF 8')]
            leaves = [(2, -11), (2, 12), (1, 5), (4, 11), (4, -12), (3, -5),
                      (5, -13), (5, 14)]
        else:
            nodes = [(1, 0, 6, 'LEAF 1 LEAF 2 LEAF 3'),
                     (2, 0, -6, 'LEAF 4 LEAF 5 LEAF 6')]
            leaves = [(1, 5), (1, -11), (1, 12), (2, -5), (2, 11), (2, -12)]
            if treatment == 'top-bw':
                nodes.append((3, 0, 24, 'LEAF 7 LEAF 8'))
                leaves.extend([(3, -13), (3, 14)])
        final = [6, -6, 14, -13] if treatment == 'all-bw' else [6, -6, 24]
        lines = ['FORMAT 4']
        lines += ['NODE %d %d %d 0 0 %d %s' % (i, parent, pdg, len(children.split())//2, children)
                  for i, parent, pdg, children in nodes]
        lines += ['DECAY_LEAF %d %d %d' % (i, parent, pdg)
                  for i, (parent, pdg) in enumerate(leaves, 1)]
        lines += ['CONTEXT 1 BORN 1 %d 10' % (len(final)+2),
                  'CORE_LEG 1 1 2 I', 'CORE_LEG 1 2 -1 I']
        lines += ['CORE_LEG 1 %d %d F' % (i, pdg) for i, pdg in enumerate(final, 3)]
        path = directory / 'decay_chain_info.dat'
        path.write_text('\n'.join(lines) + '\nEND\n')
        path.with_name('decay_internal_widths.json').write_text(json.dumps(
            dict(format=1, pdgs=[] if treatment == 'onshell' else [24])))
        return path

    def test_export_w_treatment_checks_actual_core_and_decays(self):
        with tempfile.TemporaryDirectory(prefix='ttw_topology_') as directory:
            process = Path(directory)
            with self.assertRaisesRegex(ValueError, 'Missing exported'):
                setup.export_w_treatment(process)
            for treatment in setup.W_TREATMENTS:
                path = self.topology_fixture(process, treatment)
                self.assertEqual(setup.export_w_treatment(process), (treatment, [path]))
            original = path.read_text()
            path.write_text(original.replace('CORE_LEG 1 6 -13 F', 'CORE_LEG 1 6 13 F'))
            with self.assertRaisesRegex(ValueError, 'associated'):
                setup.export_w_treatment(process)
            path.write_text(original.replace('DECAY_LEAF 3 1 12', 'DECAY_LEAF 3 1 14'))
            with self.assertRaisesRegex(ValueError, 'leptonic'):
                setup.export_w_treatment(process)

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
            self.topology_fixture(process, 'onshell')
            args = argparse.Namespace(process_dir=str(process), variant='S', top_width_lo=1.5,
                                      top_width_nlo=1.35, width_source='synthetic test only',
                                      pdf_id=None, production_scale='core-w-ht-half', ecm=13000.,
                                      seed=42, points=10, grid_points=10, iterations=1,
                                      decay_scales='separate', w_treatment='onshell',
                                      top_width_w_treatment='onshell')
            with self.assertRaisesRegex(ValueError, 'W-system CORE HT/2'):
                setup.configure(args)
            (process / 'SubProcesses/decay_chain_parameters.f90').write_text(
                'DECAY_SCALE_GROUPING PRODUCTION_SCALE_GROUPING')
            (process / 'SubProcesses/setscales.f90').write_text('w_system_core_ht_half')
            (process / 'SubProcesses/decay_chain_scales.f90').write_text('function w_system_core_ht_half')
            with contextlib.redirect_stdout(io.StringIO()):
                setup.configure(args)
            decay = (cards / 'decay_card.dat').read_text()
            self.assertIn('INDEPENDENT = decay_scale_variation_mode', decay)
            self.assertIn('SIGNED_PDG = decay_scale_grouping', decay)
            self.assertIn('1, 0.5, 2 = decay_scale_factors', decay)
            self.assertIn('AUTO = decay_width_scale_mode(6)', decay)
            self.assertIn('W_SYSTEM = production_scale_grouping', decay)
            run = RunCardNLO(str(cards / 'run_card.dat'))
            self.assertEqual(list(run['rw_rscale']), [1., .5, 2.])
            self.assertEqual(list(run['rw_fscale']), [1., .5, 2.])
            self.assertFalse(run['cut_decays'])
            self.assertEqual(run['ptj'], 0.)
            self.assertEqual(run['maxjetflavor'], 5)
            self.assertEqual(run['req_acc_fo'], -1.)
            archive = process / 'study_cards/S_onshell_core-w-ht-half_separate_42'
            manifest = json.loads((archive / 'manifest.json').read_text())
            self.assertEqual(manifest['top_width_reference_scale'], 173.)
            self.assertEqual(manifest['decay_scale_grouping'], 'separate')
            self.assertEqual(manifest['w_treatment'], 'onshell')
            self.assertEqual(manifest['top_width_w_treatment'], 'onshell')
            self.assertEqual(manifest['production_bottom_mass'], 0.)
            self.assertEqual(manifest['decay_bottom_mass'], 0.)
            self.assertEqual(manifest['top_width_bottom_mass'], 0.)
            self.assertEqual(manifest['alpha_s_flavours'], 5)
            self.assertEqual(len(manifest['topology_hashes']), 1)
            self.assertEqual(manifest['production_scale_grouping'], 'W_SYSTEM')
            self.assertEqual(manifest['production_sampling'], 'flat')
            self.assertIsNone(manifest['production_sampling_parameters'])
            self.assertNotIn('production_phase_space_sampling', decay)
            self.assertEqual(len(manifest['scale_runtime_hashes']), 3)
            self.assertEqual((archive / 'before/decay_card.dat').read_text(), 'original decay card\n')
            with self.assertRaisesRegex(ValueError, 'already exists'):
                setup.configure(args)
            args.decay_scales = 'shared'
            for invalid in (0., 1., -2., float('nan')):
                args.accuracy = invalid
                with self.assertRaisesRegex(ValueError, 'Integration accuracy'):
                    setup.configure(args)
            args.accuracy = .03
            args.job_seconds = 60.
            with contextlib.redirect_stdout(io.StringIO()):
                setup.configure(args)
            run = RunCardNLO(str(cards / 'run_card.dat'))
            self.assertEqual(run['req_acc_fo'], .03)
            self.assertEqual(run['fo_job_target_time'], 60.)
            self.assertIn('SPECIES = decay_scale_grouping', (cards / 'decay_card.dat').read_text())
            self.assertTrue((process / 'study_cards/S_onshell_core-w-ht-half_shared_42/manifest.json').is_file())
            before = (cards / 'decay_card.dat').read_text()
            args.w_treatment = 'top-bw'
            with self.assertRaisesRegex(ValueError, 'Export has W treatment'):
                setup.configure(args)
            self.topology_fixture(process, 'top-bw')
            with self.assertRaisesRegex(ValueError, 'requires top total widths'):
                setup.configure(args)
            self.assertEqual((cards / 'decay_card.dat').read_text(), before)
            self.assertFalse((process / 'study_cards/S_top-bw_core-w-ht-half_shared_42').exists())
            args.top_width_w_treatment = 'bw'
            internal_width_path = process / 'SubProcesses/P0_test/decay_internal_widths.json'
            internal_width_path.unlink()
            with self.assertRaisesRegex(ValueError, 'preserve internal decay W widths'):
                setup.configure(args)
            self.assertEqual((cards / 'decay_card.dat').read_text(), before)
            param_before = (cards / 'param_card.dat').read_bytes()
            for treatment in ('top-bw', 'all-bw'):
                self.topology_fixture(process, treatment)
                args.w_treatment = treatment
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
                    setup.configure(args)
                decay = (cards / 'decay_card.dat').read_text()
                self.assertIn('lo_decay_width(6)', decay)
                self.assertEqual('lo_decay_width(24)' in decay, treatment == 'top-bw')
                self.assertEqual((cards / 'param_card.dat').read_bytes(), param_before)
                run = RunCardNLO(str(cards / 'run_card.dat'))
                self.assertFalse(run['fixed_ren_scale'])
                self.assertFalse(run['fixed_fac_scale'])
                self.assertFalse(run['fixed_qes_scale'])
                self.assertEqual(list(run['dynamical_scale_choice']), [3])
                archive = process / ('study_cards/S_%s_core-w-ht-half_shared_42' % treatment)
                manifest = json.loads((archive / 'manifest.json').read_text())
                self.assertEqual(manifest['top_width_w_treatment'], 'bw')
                self.assertEqual(manifest['w_width'], 2.05)
                if treatment == 'all-bw':
                    self.assertIn('associated lepton', manifest['production_core_objects'])
                    self.assertIn('associated W system', manifest['production_scale_objects'])
                    self.assertNotIn('native dynamic CORE HT/2', err.getvalue())
            # Sampling is an opt-in coordinate change, with independent
            # runtime and card provenance. Invalid/old exports fail before writes.
            before = {name: (cards / name).read_bytes()
                      for name in ('param_card.dat', 'run_card.dat', 'decay_card.dat')}
            args.production_sampling = 'bad'
            with self.assertRaisesRegex(ValueError, 'flat or w-current'):
                setup.configure(args)
            args.production_sampling = 'w-current'
            self.topology_fixture(process, 'top-bw')
            args.w_treatment = 'top-bw'
            with self.assertRaisesRegex(ValueError, 'requires an all-bw'):
                setup.configure(args)
            self.topology_fixture(process, 'all-bw')
            args.w_treatment = 'all-bw'
            with self.assertRaisesRegex(ValueError, 'both core paths'):
                setup.configure(args)
            for name, contents in before.items():
                self.assertEqual((cards / name).read_bytes(), contents)
            for name in ('phase_space_kinematics.f90', 'factorized_block_kinematics.f90',
                         'decay_chain_parameters.f90', 'decay_chain_kinematics.f90',
                         'nlo_decay_kinematics.f90'):
                (process / 'SubProcesses' / name).write_bytes(
                    (Path(MG5DIR) / 'Template/fNLO/SubProcesses' / name).read_bytes())
            with contextlib.redirect_stdout(io.StringIO()):
                setup.configure(args)
            archive = process / 'study_cards/S_all-bw_core-w-ht-half_shared_42_wcurrent'
            manifest = json.loads((archive / 'manifest.json').read_text())
            self.assertEqual(manifest['production_sampling'], 'w-current')
            self.assertEqual(manifest['production_sampling_parameters'],
                             dict(mass_GeV=80.4, width_GeV=2.05))
            self.assertEqual(len(manifest['sampling_runtime_hashes']), 5)
            self.assertIn('W_CURRENT = production_phase_space_sampling',
                          (cards / 'decay_card.dat').read_text())
            for name in ('param_card.dat', 'run_card.dat'):
                self.assertEqual((cards / name).read_bytes(), before[name])
            del args.production_sampling
            for choice in ('core-ht-half', 'fixed'):
                args.production_scale = choice
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
                    setup.configure(args)
                self.assertNotIn('production_scale_grouping', (cards / 'decay_card.dat').read_text())
                run = RunCardNLO(str(cards / 'run_card.dat'))
                self.assertEqual(run['fixed_ren_scale'], choice == 'fixed')
                self.assertEqual('scale definition' in err.getvalue(), choice == 'core-ht-half')
            # A massive export must declare the same mass in its matrix
            # elements and both supplied widths before any card is changed.
            param_path = cards / 'param_card.dat'
            param_path.write_text(param_path.read_text() + 'BLOCK DECAYMASS\n 5 4.8\n')
            scheme_path = cards / 'decay_mass_scheme.json'
            scheme_path.write_text(json.dumps(dict(
                format=1, production_model='loop_sm-no_b_mass', decay_model='loop_sm',
                parameter=['decaymass', 5], alpha_s_flavours=5)))
            mass_source = process / 'Source/MODEL/get_mass_width_fcts.f'
            mass_source.parent.mkdir(parents=True)
            mass_source.write_text('GET_DECAY_MASS_FROM_ID=DC_MDL_MB')
            before = {name: (cards / name).read_bytes()
                      for name in ('param_card.dat', 'run_card.dat', 'decay_card.dat')}
            with self.assertRaisesRegex(ValueError, 'Requested decay bottom mass'):
                setup.configure(args)
            args.decay_bottom_mass = 4.8
            for value in (None, 0., 4.7, float('nan')):
                args.top_width_bottom_mass = value
                with self.assertRaisesRegex(ValueError, 'top-width-bottom-mass'):
                    setup.configure(args)
                for name, contents in before.items():
                    self.assertEqual((cards / name).read_bytes(), contents)
            args.top_width_bottom_mass = 4.8
            with contextlib.redirect_stdout(io.StringIO()):
                setup.configure(args)
            archive = process / 'study_cards/S_all-bw_fixed_shared_42_mb4p8'
            manifest = json.loads((archive / 'manifest.json').read_text())
            self.assertEqual(manifest['production_bottom_mass'], 0.)
            self.assertEqual(manifest['decay_bottom_mass'], 4.8)
            self.assertEqual(manifest['top_width_bottom_mass'], 4.8)
            self.assertEqual(manifest['bottom_mass_scheme'], 'on-shell')
            self.assertEqual(manifest['alpha_s_flavours'], 5)
            self.assertTrue(manifest['decay_mass_scheme_hash'])
            self.assertEqual((archive / scheme_path.name).read_bytes(), scheme_path.read_bytes())
            (process / 'SubProcesses/decay_chain_parameters.f90').write_text('old export')
            with self.assertRaisesRegex(ValueError, 'Re-export'):
                setup.configure(args)

    def test_reference_comparisons_are_point_matched(self):
        for points in (scales.POINTS, scales.SHARED_POINTS):
            rs = {k: k[0]*10. for k in points}
            rs['central'] = 10.
            rp = {k: 1.1*v for k, v in rs.items()}
            s = {k: 2.*v for k, v in rs.items()}
            p = {k: 1.2*v for k, v in s.items()}
            result = scales.compare_reference_values(rs, rp, s, p)
            self.assertAlmostEqual(result['shift_change']['central'], 3.)
            self.assertAlmostEqual(result['relative_shift_change']['central'], .1)
            self.assertAlmostEqual(result['product_over_strict_double_ratio']['central'], 1.2/1.1)
            band = result['product_over_strict_double_ratio']['combined21']['envelope']
            self.assertAlmostEqual(band[0], band[1])
            s['central'] = 0.
            result = scales.compare_reference_values(rs, rp, s, p)
            self.assertEqual(result['shift_change']['central'], 23.)
            self.assertIsNone(result['relative_shift_change']['central'])
            self.assertIsNone(result['product_over_strict_double_ratio']['central'])
            s['central'] = None
            self.assertIsNone(scales.compare_reference_values(rs, rp, s, p)['shift_change']['central'])

    def test_reference_report_includes_shapes_acceptances_and_charge(self):
        with tempfile.TemporaryDirectory(prefix='ttw_reference_') as directory:
            files = [Path(directory) / name for name in ('rs.HwU', 'rp.HwU', 's.HwU', 'p.HwU')]
            for path, multiplier in zip(files, (1., 1.1, 2., 2.4)):
                self.fixture(path, multiplier)
            rs, rp, s, p = (scales.load_sum([path]) for path in files)
            reference, target = scales.make_report(rs, rp), scales.make_report(s, p)
            result = scales.make_reference_report(reference, target)
            hist = result['histograms']['R04_b25 W+ 1b test1'][0]
            self.assertAlmostEqual(hist['absolute']['shift_change']['central'], 1.2)
            self.assertAlmostEqual(hist['normalized_to_fiducial']['relative_shift_change']['central'], 0.)
            derived = result['derived']['R04_b25 W+']
            self.assertAlmostEqual(derived['acceptance_2b']['strict_over_reference']['central'], 1.)
            self.assertAlmostEqual(derived['charge_asymmetry_fiducial_1b']['relative_shift_change']['central'], 0.)
            json.dumps(result, allow_nan=False)
            output = Path(directory) / 'comparison.json'
            argv = ['ttw_product_scales.py', '--strict', str(files[2]), '--product', str(files[3]),
                    '--w-treatment', 'all-bw', '--reference-strict', str(files[0]),
                    '--reference-product', str(files[1]), '--output', str(output)]
            with mock.patch('sys.argv', argv), contextlib.redirect_stdout(io.StringIO()):
                scales.main()
            cli = json.loads(output.read_text())
            self.assertEqual(cli['w_treatment'], 'all-bw')
            self.assertEqual(cli['reference']['w_treatment'], 'onshell')
            self.assertEqual(cli['production_scale'], 'core-w-ht-half')
            self.assertTrue(cli['reference_comparison']['scale_definitions_match'])
            self.assertNotEqual(scales.production_scale_definition('core-ht-half', 'all-bw'),
                                cli['production_scale_definition'])
            self.assertEqual(scales.production_scale_definition('core-ht-half', 'onshell'),
                             cli['production_scale_definition'])
            # The same interface can compare central choices within one W treatment.
            argv[argv.index('all-bw')] = 'onshell'
            argv[-1] = str(Path(directory) / 'central-scales.json')
            argv += ['--reference-production-scale', 'fixed']
            with mock.patch('sys.argv', argv), contextlib.redirect_stdout(io.StringIO()):
                scales.main()
            cli = json.loads(Path(argv[argv.index('--output')+1]).read_text())
            self.assertEqual(cli['reference']['production_scale'], 'fixed')
            self.assertEqual(cli['production_scale'], 'core-w-ht-half')
            self.assertFalse(cli['reference_comparison']['scale_definitions_match'])
            self.assertIn('not just W-width', cli['reference_comparison']['scale_definition_note'])
            argv[argv.index('--reference-production-scale') + 1] = 'core-w-ht-half'
            argv[argv.index('--output') + 1] = str(Path(directory) / 'bottom-masses.json')
            argv += ['--decay-bottom-mass', '4.8', '--reference-decay-bottom-mass', '0']
            with mock.patch('sys.argv', argv), contextlib.redirect_stdout(io.StringIO()):
                scales.main()
            cli = json.loads(Path(argv[argv.index('--output')+1]).read_text())
            self.assertEqual(cli['decay_bottom_mass'], 4.8)
            self.assertEqual(cli['reference']['decay_bottom_mass'], 0.)
            self.assertTrue(cli['reference_comparison']['scale_definitions_match'])
            target['histograms']['R04_b25 W+ 1b test1'][0]['edges'] = [0., 1.]
            with self.assertRaisesRegex(ValueError, 'bin edges differ'):
                scales.make_reference_report(reference, target)

    def test_setup_defaults_to_w_system_scale(self):
        argv = ['ttw_product_setup.py', 'configure', '--process-dir', '/tmp/TTW', '--variant', 'S',
                '--top-width-lo', '1.5', '--top-width-nlo', '1.35',
                '--top-width-w-treatment', 'bw', '--w-treatment', 'all-bw', '--width-source', 'test only']
        with mock.patch('sys.argv', argv), mock.patch.object(setup, 'configure') as configure:
            setup.main()
        self.assertEqual(configure.call_args[0][0].production_scale, 'core-w-ht-half')


if __name__ == '__main__':
    unittest.main()
