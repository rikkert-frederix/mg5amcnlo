"""Massive top decays must not change the massless production calculation."""

import copy
import json
import os
import tempfile
from pathlib import Path

import tests.unit_tests as unittest
from madgraph import InvalidCmd
from madgraph.fks import fks_decay_masses, fks_helas_objects
from madgraph.interface.master_interface import MasterCmd
from madgraph.iolibs import export_fks, file_writers
from models import model_reader


class TestFKSDecayMasses(unittest.TestCase):
    def setUp(self):
        self.addCleanup(os.chdir, os.getcwd())
        self.command = MasterCmd()
        self.command.exec_cmd('import model loop_sm-no_b_mass', printcmd=False, precmd=True)

    def test_private_model_preserves_production_and_shares_inputs(self):
        production = self.command._curr_model
        interactions = [(i['id'], dict(i['couplings'])) for i in production['interactions']]
        expressions = {p.name: p.expr for group in production['couplings'].values() for p in group}
        particles = [(p['pdg_code'], p['mass'], copy.deepcopy(p['counterterm']))
                     for p in production['particles']]
        decay = fks_decay_masses.massive_bottom_decay_model(production, 4.8)
        self.assertEqual(production.get_nflav(), 5)
        self.assertEqual(production.get_particle(5)['mass'], 'ZERO')
        self.assertEqual(decay.get_particle(5)['mass'], 'dc_mdl_MB')
        self.assertTrue(decay.get_particle(5)['counterterm'])
        self.assertEqual(interactions, [(i['id'], dict(i['couplings']))
                                        for i in production['interactions']])
        self.assertEqual(particles, [(p['pdg_code'], p['mass'], p['counterterm'])
                                    for p in production['particles']])
        self.assertEqual(expressions, {p.name: p.expr for group in production['couplings'].values()
                                       for p in group if not p.name.startswith('DC_')})
        reader = model_reader.ModelReader(production)
        reader.set_parameters_and_couplings()
        self.assertEqual(reader['parameter_dict']['dc_mdl_MB'], 4.8)
        for name in ('MT', 'MW', 'WW', 'WT'):
            self.assertEqual(reader['parameter_dict']['mdl_' + name],
                             reader['parameter_dict']['dc_mdl_' + name])
        self.assertEqual(reader['coupling_dict']['GC_11'], reader['coupling_dict']['DC_GC_11'])
        self.assertEqual(reader['coupling_dict']['GC_5'], reader['coupling_dict']['DC_GC_5'])
        self.assertIs(fks_decay_masses.massive_bottom_decay_model(production, 4.8), decay)
        with self.assertRaisesRegex(InvalidCmd, 'Re-import'):
            fks_decay_masses.massive_bottom_decay_model(production, 4.7)

    def test_option_guards(self):
        for value in ('-1', 'nan', 'inf', '4 5'):
            with self.assertRaisesRegex(InvalidCmd, 'finite and nonnegative'):
                self.command.exec_cmd('set decay_bottom_mass ' + value, printcmd=False, precmd=True)
        self.command.exec_cmd('set decay_bottom_mass 4.8', printcmd=False, precmd=True)
        with self.assertRaisesRegex(InvalidCmd, 'requires fNLO generation'):
            self.command.exec_cmd('generate u u~ > t t~, t > b w+', printcmd=False, precmd=True)
        with self.assertRaisesRegex(InvalidCmd, 'explicit top decay'):
            self.command.exec_cmd('generate u u~ > t t~ [real=QCD]', printcmd=False, precmd=True)
        self.command.exec_cmd('import model loop_sm', printcmd=False, precmd=True)
        with self.assertRaisesRegex(InvalidCmd, 'loop_sm-no_b_mass'):
            self.command.exec_cmd('generate u u~ > t t~ [real=QCD], t > b w+ [real=QCD]',
                                  printcmd=False, precmd=True)

    def test_bottom_beams_and_decay_fks_models_are_distinct(self):
        self.command.exec_cmd('set decay_bottom_mass 4.8', printcmd=False, precmd=True)
        self.command.exec_cmd(
            'generate b b~ > t t~ [real=QCD], (t > b w+ [real=QCD]), (t~ > b~ w- [real=QCD])',
            printcmd=False, precmd=True)
        bundle = self.command._fks_multi_proc
        production = bundle.production
        for amplitude in production.get_born_amplitudes():
            self.assertEqual(amplitude['process']['model'].get_particle(5)['mass'], 'ZERO')
        for member in bundle.decays:
            for amplitude in member.get_born_amplitudes():
                self.assertEqual(amplitude['process']['model'].get_particle(5)['mass'], 'dc_mdl_MB')
        for real in production.get_real_amplitudes():
            self.assertEqual(real['process']['model'].get_particle(5)['mass'], 'ZERO')
        with self.assertRaisesRegex(InvalidCmd, 'before generating'):
            self.command.exec_cmd('set decay_bottom_mass 0', printcmd=False, precmd=True)

    def test_mass_getters_have_distinct_fortran_paths(self):
        model = self.command._curr_model
        fks_decay_masses.massive_bottom_decay_model(model, 4.8)
        with tempfile.TemporaryDirectory() as directory:
            makeinc = Path(directory) / 'makeinc.inc'
            makeinc.write_text('MODEL = couplings.o\n')
            exporter = object.__new__(export_fks.ProcessExporterFortranFKS)
            path = Path(directory) / 'mass.f'
            with file_writers.FortranWriter(str(path)) as writer:
                exporter.write_get_mass_width_file(writer, str(makeinc), model)
            source = path.read_text().upper()
            production = source.split('DOUBLE PRECISION FUNCTION GET_WIDTH_FROM_ID')[0]
            self.assertNotIn('DC_MDL_MB', production)
            self.assertIn('GET_DECAY_MASS_FROM_ID=DC_MDL_MB', source)
            self.assertIn('IEEE_IS_FINITE(DC_MDL_MB)', source)
            self.assertIn('ABS(ID).EQ.5', source)

    def test_bw_bundle_keeps_massive_density_providers(self):
        self.command.exec_cmd('set decay_bottom_mass 4.8', printcmd=False, precmd=True)
        self.command.exec_cmd(
            'generate u d~ > t t~ mu+ vm QCD=2 QED=2 [real=QCD], '
            '(t > b e+ ve QED=2 [real=QCD]), (t~ > b~ e- ve~ QED=2 [real=QCD])',
            printcmd=False, precmd=True)
        helas = fks_helas_objects.FKSHelasMultiProcess(self.command._fks_multi_proc)
        for matrix_element in helas.get_matrix_elements():
            for node, component in matrix_element.spin_density_plan['components'].items():
                model = component['born']['matrix_element']['processes'][0]['model']
                self.assertEqual(model.get_particle(5)['mass'], 'ZERO' if node == 0 else 'dc_mdl_MB')

    def test_mixed_w_mode_keeps_internal_bw_widths(self):
        self.command.exec_cmd('set decay_bottom_mass 4.8', printcmd=False, precmd=True)
        self.command.exec_cmd(
            'generate u d~ > t t~ w+ QCD=2 QED=1 [real=QCD], '
            '(t > b e+ ve QED=2 [real=QCD]), (t~ > b~ e- ve~ QED=2 [real=QCD]), '
            'w+ > mu+ vm', printcmd=False, precmd=True)
        helas = fks_helas_objects.FKSHelasMultiProcess(self.command._fks_multi_proc)
        with tempfile.TemporaryDirectory() as directory:
            os.chdir(directory)
            for matrix_element in helas.get_matrix_elements():
                export_fks.ProcessExporterFortranFKS.write_spin_density_internal_widths(
                    matrix_element.spin_density_plan)
                self.assertEqual(json.loads(Path('decay_internal_widths.json').read_text()),
                                 {'format': 1, 'pdgs': [24]})
                currents = [wf for wf in matrix_element.born_me.get_all_wavefunctions()
                            if abs(wf.get('pdg_code')) == 24 and wf.get('mothers')]
                self.assertTrue(any(wf.get('decay_node_id') and wf.get('width') == 'ZERO'
                                    for wf in currents))
                self.assertTrue(any(not wf.get('decay_node_id') and wf.get('width') != 'ZERO'
                                    for wf in currents))
