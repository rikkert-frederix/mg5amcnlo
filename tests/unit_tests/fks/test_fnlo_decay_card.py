"""Numerical tests of fNLO decay cards, running widths and local scales."""

import math
import os
import shutil
import subprocess
import tempfile

import tests.unit_tests as unittest
from madgraph import MG5DIR
from madgraph.fks import fks_decay


class TestFNLODecayCard(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which('gfortran')
        if not compiler:
            raise unittest.SkipTest('gfortran is required for fNLO runtime tests')
        cls.build = tempfile.TemporaryDirectory(prefix='fnlo_decay_card_')
        cls.addClassCleanup(cls.build.cleanup)
        template = os.path.join(MG5DIR, 'Template', 'fNLO', 'SubProcesses')
        fixtures = os.path.join(MG5DIR, 'tests', 'input_files', 'fks_decay')
        cls.executable = os.path.join(cls.build.name, 'check_card')
        sources = [os.path.join(fixtures, 'decay_card_runtime_stubs.f90')]
        sources += [os.path.join(template, name + '.f90') for name in (
            'decay_chain_parameters', 'factorized_phase_space',
            'dummy_fct', 'decay_chain_scales')]
        sources.append(os.path.join(fixtures, 'decay_card_runtime_driver.f90'))
        result = subprocess.run(
            [compiler, '-O0', '-g', '-fcheck=all', '-ffree-line-length-none',
             '-o', cls.executable] + sources,
            cwd=cls.build.name, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True)
        if result.returncode:
            raise AssertionError(result.stdout)

    def card(self, **options):
        defaults = dict(nlo_width_pdgs={6}, nlo_widths={6: 1.8},
                        decay_scale_variation_mode='INDEPENDENT',
                        decay_scale_factors=(1., .5, 2.))
        defaults.update(options)
        return fks_decay.decay_card_text({6: 2.}, {6: 100.}, **defaults)

    def run_card(self, card, mode='scales', error=None):
        with tempfile.TemporaryDirectory(dir=self.build.name) as directory:
            with open(os.path.join(directory, 'decay_card.dat'), 'w') as stream:
                stream.write(card)
            result = subprocess.run([self.executable, mode], cwd=directory,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True)
        if error:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(error, result.stdout)
            return
        self.assertEqual(result.returncode, 0, result.stdout)
        return {line.split()[0]: line.split()[1:]
                for line in result.stdout.splitlines() if line.split()}

    @staticmethod
    def alpha(mu):
        return .1 / (1. + .1 * math.log(mu / 100.))

    def width(self, mu):
        return 2. - .2 * self.alpha(mu) / self.alpha(100.)

    def test_dynamic_scales_and_widths_follow_each_occurrence(self):
        data = self.run_card(self.card(decay_dynamical_scale_choices={6: 3}))
        self.assertEqual(list(map(float, data['SCALES'])), [30., 40.])
        for actual, expected in zip(data['WIDTHS'], [self.width(30.), self.width(80.)]):
            self.assertAlmostEqual(float(actual), expected, places=13)
        expected = -sum((self.width(mu)-2.)/2. for mu in [60., 80.])
        self.assertAlmostEqual(float(data['COUNTERTERM'][0]), expected, places=13)
        self.assertEqual(float(data['DENOMINATOR'][0]), 1.)
        for actual, mu in zip(data['LOCAL_WIDTHS'], [60., 80.]):
            self.assertAlmostEqual(float(actual), 2./self.width(mu), places=13)
        self.assertAlmostEqual(float(data['COUPLING'][0]),
                               4.*math.pi*self.alpha(60.), places=13)

    def test_fixed_scale_reweighting_matches_direct_scale_change(self):
        varied = self.run_card(self.card())
        direct_card = fks_decay.decay_card_text(
            {6: 2.}, {6: 200.}, nlo_width_pdgs={6},
            nlo_widths={6: self.width(200.)})
        direct = self.run_card(direct_card)
        for key in ['COUNTERTERM', 'DENOMINATOR', 'COUPLING', 'LOCAL_WIDTHS']:
            for first, second in zip(varied[key], direct[key]):
                self.assertAlmostEqual(float(first), float(second), places=13)

    def test_mixed_orders_and_global_lo(self):
        for production, decay, expected in [('LO', 'NLO', ['F', 'T', 'T', 'F']),
                                             ('NLO', 'LO', ['T', 'F', 'F', 'F']),
                                             ('LO', 'LO', ['F', 'F', 'F', 'F'])]:
            data = self.run_card(self.card(production_order=production,
                                           decay_order=decay), 'orders')
            self.assertEqual(data['ORDERS'], expected)
            self.assertAlmostEqual(float(data['COUNTERTERM'][0]),
                                   .2 if decay == 'NLO' else 0.)
        data = self.run_card(self.card(decay_perturbative_orders={6: 'LO'}), 'orders')
        self.assertEqual(data['ORDERS'], ['T', 'F', 'F', 'F'])
        data = self.run_card(self.card(nlo_decay_combination='MULTIPLICATIVE'), 'all_lo')
        self.assertEqual(data['ORDERS'], ['F', 'F', 'F', 'F'])
        self.assertEqual(float(data['COUNTERTERM'][0]), 0.)

    def test_alpha_s_dependent_born_rejects_dynamic_scales(self):
        self.run_card(self.card(decay_dynamical_scale_choices={6: 3}), 'qcd_born',
                      error='event-by-event decay scales require an alpha_s-independent Born')
        self.run_card(self.card(decay_width_scale_modes={6: 'AUTO'}), 'qcd_born',
                      error='AUTO widths require an alpha_s-independent Born')
        self.run_card(self.card(), 'qcd_born', error='no LO width at every factor')

    def test_explicit_widths(self):
        options = dict(lo_width_variations={(6, .5): 1.9, (6, 2.): 2.1},
                       nlo_width_variations={(6, .5): 1.7, (6, 2.): 1.95})
        data = self.run_card(self.card(**options), 'qcd_born')
        self.assertAlmostEqual(float(data['DENOMINATOR'][0]), (2./2.1)**2)
        self.assertAlmostEqual(float(data['COUNTERTERM'][0]), -2.*(1.95-2.1)/2.1)
        # Switching to LO must also work without deleting unused width tables.
        self.run_card(self.card(decay_order='LO', **options), 'orders')

    def test_only_unversioned_assignment_cards_are_supported(self):
        card = self.card()
        self.assertNotIn('= format', card)
        self.assertNotIn('dummy_width', card.lower())
        self.run_card(card)
        for record in ['FORMAT 3', 'FORMAT 4', 'FORMAT 5', 'FORMAT 6',
                       'END', 'PRODUCTION_ORDER LO']:
            with self.subTest(record=record):
                self.run_card(card + '\n' + record,
                              error='expected value = parameter')
        for record, name in [('6 = format', 'FORMAT'),
                             ('0.1 = dummy_width_ratio', 'DUMMY_WIDTH_RATIO')]:
            with self.subTest(record=record):
                self.run_card(card + '\n' + record, error='unknown keyword ' + name)
        self.run_card(card + '\n2.0 = decay_width(6)',
                      error='unexpected index on DECAY_WIDTH')
        self.run_card(card + '\n2.0 = lo_decay_width 6',
                      error='malformed decay-card parameter name')
        self.run_card(card + '\n6 2.0 = lo_decay_width',
                      error='require a PDG index')

    def test_comments_case_and_parameter_order(self):
        card = self.card(decay_dynamical_scale_choices={6: 3})
        records = [line.split('!', 1)[0].strip() for line in card.splitlines()
                   if line.strip() and not line.lstrip().startswith('#')]
        reordered = '\n'.join(line.swapcase() + ' # trailing comment'
                              for line in reversed(records))
        self.assertEqual(self.run_card(card), self.run_card(reordered))
        self.assertEqual(self.run_card(card), self.run_card(
            card.replace('decay_scale_factors', 'Decay_Scale_Factors')))
        self.run_card(card + '\nLO = production_order\n', error='duplicate PRODUCTION_ORDER')
        self.run_card(card + '\n0 = decay_dynamical_scale_choice(999)\n', error='unknown PDG')

    def test_standalone_decay_running_width_normalization(self):
        data = self.run_card(self.card(decay_dynamical_scale_choices={6: 3}), 'standalone')
        self.assertAlmostEqual(float(data['DENOMINATOR'][0]),
                               1.8/self.width(60.), places=13)

    def test_undecayed_scales_do_not_initialize_decay_card(self):
        data = self.run_card('', 'no_decay')
        self.assertEqual(data['NO_DECAY'], ['T', 'T', 'T'])

    def test_widths_and_scales_cover_the_same_species(self):
        self.run_card(self.card() + '\n100. = decay_ren_scale(999)\n',
                      error='a decay scale has no physical width')
        self.run_card(self.card() + '\n1. = lo_decay_width(999)\n',
                      error='a physical width has no renormalisation scale')
