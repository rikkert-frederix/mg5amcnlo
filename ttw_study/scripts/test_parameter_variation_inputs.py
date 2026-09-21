"""Physical-input/card tests; no cross-section integration or fake result."""
import contextlib
import copy
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from campaign import ROOT, STUDY, digest
from parameter_variation_inputs import SCENARIOS, validate, validate_scenario
from parameter_variation_cards import configure
from madgraph.various.banner import RunCardNLO
from models.check_param_card import ParamCard
from tests.unit_tests.fks import test_ttw_product_tools

INPUTS = STUDY/'inputs/parameter_variation_inputs_bounded_v1.json'


class TestParameterVariationInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads(INPUTS.read_text())
        cls.benchmark, cls.preflight = validate(cls.record)

    def test_all_width_scales_and_installed_couplings(self):
        import lhapdf
        lhapdf.setVerbosity(0)
        self.assertEqual(sum(len(s['top_widths']) for s in self.record['scenarios']), 48)
        for scenario in self.record['scenarios']:
            pdf = lhapdf.mkPDF(scenario['PDF']['name'], 0)
            for row in scenario['coupling_checkpoints']:
                self.assertAlmostEqual(pdf.alphasQ(row['Q_GeV']), row['alpha_s'], places=14)
            for row in scenario['top_widths']:
                self.assertAlmostEqual(row['gamma_lo']+row['qcd_coefficient']*pdf.alphasQ(row['muR_GeV']),
                                       row['gamma_nlo_pdf'], places=13)

    def test_inconsistent_mass_width_scale_and_pdf_are_rejected(self):
        scenario = self.record['scenarios'][5]
        for key, value in [('top_yukawa_GeV', 172.5), ('decay_scale_reference_GeV', 172.5),
                           ('fixed_W_width_GeV', 2.1), ('production_mb_GeV', 4.8)]:
            wrong = copy.deepcopy(scenario)
            wrong[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_scenario(wrong, SCENARIOS[5], self.benchmark, self.preflight)
        for change in ('missing_scale', 'wrong_pdf', 'wrong_width'):
            wrong = copy.deepcopy(scenario)
            if change == 'missing_scale': wrong['top_widths'].pop()
            if change == 'wrong_pdf': wrong['PDF']['id'] = 333900
            if change == 'wrong_width': wrong['top_widths'][0]['gamma_nlo_pdf'] += .001
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_scenario(wrong, SCENARIOS[5], self.benchmark, self.preflight)

    @staticmethod
    def fixture(process, mode):
        (process/'Cards').mkdir(parents=True)
        (process/'SubProcesses').mkdir()
        (process/'FixedOrderAnalysis').mkdir()
        shutil.copy2(ROOT/'models/loop_sm/restrict_no_b_mass.dat', process/'Cards/param_card.dat')
        template = ROOT/'Template/fNLO/Cards/run_card.dat'
        RunCardNLO().write(str(process/'Cards/run_card.dat'), template=str(template))
        shutil.copy2(process/'Cards/run_card.dat', process/'Cards/run_card_default.dat')
        (process/'Cards/decay_card.dat').write_text('Initial synthetic fixture\n')
        source = ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90'
        shutil.copy2(source, process/'FixedOrderAnalysis'/source.name)
        for name in ('setscales.f90', 'decay_chain_scales.f90', 'decay_chain_parameters.f90',
                     'phase_space_kinematics.f90', 'factorized_block_kinematics.f90',
                     'decay_chain_kinematics.f90', 'nlo_decay_kinematics.f90'):
            shutil.copy2(ROOT/'Template/fNLO/SubProcesses'/name, process/'SubProcesses'/name)
        test_ttw_product_tools.TestTTWProductTools.topology_fixture(process, mode)

    def test_all_cards_move_only_matched_physical_inputs(self):
        import lhapdf
        lhapdf.setVerbosity(0)
        with tempfile.TemporaryDirectory(prefix='parameter_scan_card_test_') as folder:
            scenarios = [(row, 'onshell') for row in self.record['scenarios']]
            scenarios += [(self.record['scenarios'][5], 'all-bw'),
                          (self.record['scenarios'][3], 'top-bw')]
            for index, (scenario, mode) in enumerate(scenarios):
                with self.subTest(scenario=scenario['name'], mode=mode):
                    process = Path(folder)/str(index)
                    self.fixture(process, mode)
                    variant = 'S' if index % 2 == 0 else 'Pi'
                    with contextlib.redirect_stdout(io.StringIO()):
                        name, archive = configure(process, INPUTS, scenario['name'], variant, 900001+index, mode)
                    manifest = json.loads((archive/'manifest.json').read_text())
                    param = ParamCard(str(archive/'param_card.dat'))
                    run = RunCardNLO(str(archive/'run_card.dat'))
                    self.assertEqual(float(param['mass'].get((6,)).value), scenario['mt_GeV'])
                    self.assertEqual(float(param['yukawa'].get((6,)).value), scenario['mt_GeV'])
                    self.assertEqual(float(param['mass'].get((5,)).value), 0.)
                    self.assertAlmostEqual(float(param['mass'].get((24,)).value), 80.385, places=11)
                    self.assertEqual(run['ebeam1'], scenario['beam_energy_GeV'])
                    self.assertEqual(run['ebeam2'], scenario['beam_energy_GeV'])
                    self.assertEqual(list(run['lhaid']), [scenario['PDF']['id']])
                    self.assertEqual(list(run['reweight_pdf']), [False])
                    self.assertEqual(list(run['rw_rscale']), [1., .5, 2.])
                    self.assertEqual(list(run['rw_fscale']), [1., .5, 2.])
                    self.assertEqual(manifest['top_width_reference_scale'], scenario['mt_GeV'])
                    self.assertEqual(manifest['w_width'], scenario['fixed_W_width_GeV'])
                    self.assertEqual(manifest['production_scale_grouping'], 'W_SYSTEM')
                    self.assertEqual(manifest['parameter_variation_inputs_sha256'], digest(INPUTS))
                    decay = {}
                    for line in (archive/'decay_card.dat').read_text().splitlines():
                        line = line.split('!', 1)[0].strip()
                        if line and not line.startswith('#') and '=' in line:
                            value, key = line.split('=', 1)
                            decay[key.strip()] = value.strip()
                    self.assertEqual(float(decay['decay_ren_scale(6)']), scenario['mt_GeV'])
                    self.assertEqual(decay['decay_scale_grouping'], 'SIGNED_PDG')
                    self.assertEqual(decay['decay_width_scale_mode(6)'], 'AUTO')
                    self.assertEqual(decay['nlo_decay_combination'],
                                     'ADDITIVE' if variant == 'S' else 'MULTIPLICATIVE')
                    g0, gnlo = (float(decay[key+'(6)']) for key in ('lo_decay_width', 'nlo_decay_width'))
                    pdf = lhapdf.mkPDF(scenario['PDF']['name'], 0)
                    for width in scenario['top_widths']:
                        if width['w_treatment'] == ('onshell' if mode == 'onshell' else 'bw'):
                            running = g0+(gnlo-g0)*pdf.alphasQ(width['muR_GeV'])/pdf.alphasQ(scenario['mt_GeV'])
                            self.assertAlmostEqual(running, width['gamma_nlo_pdf'], places=12)
                    self.assertEqual('lo_decay_width(24)' in decay, mode != 'all-bw')
                    if mode != 'all-bw':
                        self.assertEqual(float(decay['lo_decay_width(24)']), scenario['fixed_W_width_GeV'])
                    self.assertFalse((process/'Events'/name).exists())
                    for filename, sha in manifest['hashes'].items():
                        if filename.endswith('.dat'): self.assertEqual(digest(archive/filename), sha)
                    # The archived common setup still retains its own original
                    # PDF-reweight option and hash; only the final scan is used.
                    base = json.loads(Path(manifest['base_configuration']).read_text())
                    self.assertTrue(base['settings']['reweight_pdf'])
                    self.assertNotEqual(base['hashes']['run_card.dat'], manifest['hashes']['run_card.dat'])
                    with self.assertRaises(ValueError):
                        configure(process, INPUTS, scenario['name'], variant, 900001+index, mode)


if __name__ == '__main__':
    unittest.main()
