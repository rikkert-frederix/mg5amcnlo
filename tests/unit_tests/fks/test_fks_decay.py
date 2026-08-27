################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of MadGraph5_aMC@NLO.
#
################################################################################

"""Tests for decay chains attached to an undecayed FKS production core."""

from __future__ import absolute_import

import copy
import os
import sys
import tempfile

root_path = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
sys.path.insert(0, os.path.join(root_path, '..', '..'))

import tests.unit_tests as unittest

from madgraph import InvalidCmd, MG5DIR
from madgraph.fks import fks_common
from madgraph.fks import fks_decay
from madgraph.fks import fks_helas_objects
from madgraph.interface.master_interface import MasterCmd


class TestFKSDecayChains(unittest.TestCase):

    @staticmethod
    def generate(process):
        command = MasterCmd()
        command.exec_cmd('import model loop_sm', printcmd=False, precmd=True)
        command.exec_cmd('generate %s' % process,
                         printcmd=False, precmd=True)
        return command

    @staticmethod
    def core_snapshot(fks_multi):
        result = []
        for born in fks_multi['born_processes']:
            result.append({
                'born': [
                    leg.get('id')
                    for leg in born.born_amp['process']['legs']],
                'reals': [[
                    leg.get('id') for leg in real.process['legs']]
                    for real in born.real_amps],
                'infos': [[(
                    info['i'], info['j'], info['ij'],
                    tuple(info['splitting_type']),
                    info['need_color_links'])
                    for info in real.fks_infos]
                    for real in born.real_amps],
                'partners': [copy.deepcopy(real.fks_j_from_i)
                             for real in born.real_amps]})
        return result

    @staticmethod
    def golden(name):
        path = os.path.join(
            MG5DIR, 'tests', 'input_files', 'fks_decay', name)
        with open(path) as stream:
            return stream.read()

    def assert_local_widths(self, matrix_element, expected_nodes):
        expected_widths = {
            6: 'FNLO_DECAY_DUMMY_WIDTH_RATIO()*mdl_MT',
            24: 'FNLO_DECAY_DUMMY_WIDTH_RATIO()*mdl_MW'}
        found_nodes = set()
        for wavefunction in matrix_element.get_all_wavefunctions():
            pdg = abs(wavefunction.get('pdg_code'))
            if pdg not in expected_widths:
                continue
            node_id = wavefunction.get('decay_node_id')
            if node_id:
                found_nodes.add(node_id)
                self.assertEqual(
                    wavefunction.get('width'), expected_widths[pdg])
            else:
                self.assertEqual(wavefunction.get('width'), 'ZERO')
        self.assertEqual(found_nodes, set(expected_nodes))

    def test_fks_skeleton_is_unchanged_and_undecayed(self):
        plain = self.generate('u u~ > t t~ [real=QCD]')
        decayed = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b, w+ > u d~)')

        self.assertEqual(
            self.core_snapshot(plain._fks_multi_proc),
            self.core_snapshot(decayed._fks_multi_proc))
        fks_process = decayed._fks_multi_proc['born_processes'][0]
        self.assertEqual(
            [leg.get('id') for leg in
             fks_process.born_amp['process']['legs']],
            [2, -2, 6, -6])
        for real in fks_process.real_amps:
            self.assertNotIn(24, [leg.get('id')
                                  for leg in real.process['legs']])
            self.assertNotIn(5, [leg.get('id')
                                 for leg in real.process['legs']])

    def test_nested_decay_widths_color_links_and_metadata(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b, w+ > u d~)')
        model = command._curr_model
        top_width = model.get_particle(6).get('width')
        w_width = model.get_particle(24).get('width')

        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)
        self.assertEqual(len(helas['matrix_elements']), 1)
        matrix_element = helas['matrix_elements'][0]

        self.assertEqual(
            [leg.get('id') for leg in
             matrix_element.born_me['processes'][0].get_legs_with_decays()],
            [2, -2, 2, -1, 5, -6])
        self.assertEqual(
            [(node['pdg'], node['qcd_order'], node['carrier_leaf'])
             for node in matrix_element.decay_metadata['nodes']],
            [(6, 0, 3), (24, 0, 0)])

        self.assert_local_widths(matrix_element.born_me, [1, 2])
        for real in matrix_element.real_processes:
            self.assert_local_widths(real.matrix_element, [1, 2])
        self.assertEqual(model.get_particle(6).get('width'), top_width)
        self.assertEqual(model.get_particle(24).get('width'), w_width)

        # The W daughters are visible legs 3 and 4.  They must never carry a
        # production colour link; the top is represented by visible b leg 5.
        visible_links = [tuple(link['link'])
                         for link in matrix_element.color_links]
        self.assertTrue(all(3 not in link and 4 not in link
                            for link in visible_links))
        self.assertIn((5, 5), visible_links)
        self.assertIn({
            'core_first': 3,
            'core_second': 3,
            'visible_first': 5,
            'visible_second': 5,
            'generated_index': 6},
            matrix_element.decay_metadata['color_links'])

        metadata_text = fks_decay.decay_chain_info_text(
            matrix_element.decay_metadata)
        self.assertEqual(
            metadata_text,
            self.golden('decay_chain_info_nested.dat'))
        with tempfile.TemporaryDirectory() as output_dir:
            fks_decay.write_decay_chain_info(
                output_dir, matrix_element.decay_metadata)
            with open(os.path.join(output_dir, 'decay_chain_info.dat')) \
                    as stream:
                self.assertEqual(stream.read(), metadata_text)

        card_text = fks_decay.decay_card_text(
            {6: 1.4915, 24: 2.0476}, {6: 173.0, 24: 80.419})
        self.assertEqual(card_text, (
            '# FNLO_DECAY_CARD\n'
            '# Runtime parameters for fixed-on-shell decay chains.\n'
            '# DECAY_WIDTH entries are physical total widths in GeV.\n'
            '# DECAY_REN_SCALE entries are independent decay scales in GeV.\n'
            'FORMAT 2\n'
            'DUMMY_WIDTH_RATIO 1.0000000000000001e-01\n'
            'PRODUCTION_REN_SCALE_MOMENTA CORE\n'
            'DECAY_WIDTH 6 1.4915000000000000e+00\n'
            'DECAY_REN_SCALE 6 1.7300000000000000e+02\n'
            'DECAY_WIDTH 24 2.0476000000000001e+00\n'
            'DECAY_REN_SCALE 24 8.0418999999999997e+01\n'
            'END\n'))

        decayed_scale_text = fks_decay.decay_card_text(
            {6: 1.4915}, {6: 173.0},
            production_scale_momenta='decayed')
        self.assertIn(
            'PRODUCTION_REN_SCALE_MOMENTA DECAYED\n',
            decayed_scale_text)

    def test_qcd_decay_order_is_stored_per_node(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b, w+ > u d~ g)')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)
        matrix_element = helas['matrix_elements'][0]
        self.assertEqual(
            [(node['pdg'], node['qcd_order'])
             for node in matrix_element.decay_metadata['nodes']],
            [(6, 0), (24, 1)])

    def test_two_decays_and_both_madloop_modes(self):
        command = self.generate(
            'u u~ > t t~ [QCD], '
            '(t > w+ b, w+ > u d~), '
            '(t~ > w- b~, w- > d u~)')

        for optimized in [False, True]:
            helas = fks_helas_objects.FKSHelasMultiProcess(
                command._fks_multi_proc, loop_optimized=optimized)
            matrix_element = helas['matrix_elements'][0]
            self.assertEqual(
                [leg.get('id') for leg in matrix_element.born_me[
                    'processes'][0].get_legs_with_decays()],
                [2, -2, 2, -1, 5, 1, -2, -5])
            self.assertIsNotNone(matrix_element.virt_matrix_element)
            self.assert_local_widths(
                matrix_element.born_me, [1, 2, 3, 4])
            for real in matrix_element.real_processes:
                self.assert_local_widths(
                    real.matrix_element, [1, 2, 3, 4])
            self.assert_local_widths(
                matrix_element.virt_matrix_element, [1, 2, 3, 4])
            self.assertTrue(matrix_element.virt_matrix_element[
                'born_color_basis'])
            self.assertTrue(matrix_element.virt_matrix_element[
                'loop_color_basis'])
            for diagram in matrix_element.virt_matrix_element.get_loop_diagrams():
                for amplitude in diagram.get_loop_amplitudes():
                    self.assertEqual(
                        sum(amplitude.get('pairing')),
                        len(amplitude.get('mothers')))
            if optimized:
                self.assertEqual(
                    fks_decay.decay_chain_info_text(
                        matrix_element.decay_metadata),
                    self.golden('decay_chain_info_two_tops.dat'))

    def test_only_equivalent_concrete_assignments_are_combined(self):
        command = self.generate(
            'u u~ > z z [real=QCD], (z > l+ l-)')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc)
        self.assertEqual(len(helas['matrix_elements']), 2)
        self.assertEqual(
            sum(len(matrix_element.born_me['processes'])
                for matrix_element in helas['matrix_elements']),
            3)
        visible_processes = [[tuple(
            leg.get('id') for leg in process.get_legs_with_decays())
            for process in matrix_element.born_me['processes']]
            for matrix_element in helas['matrix_elements']]
        self.assertEqual(visible_processes, [[
            (2, -2, -11, 11, -11, 11),
            (2, -2, -13, 13, -13, 13)], [
            (2, -2, -11, 11, -13, 13)]])

    def test_mixed_decayed_and_undecayed_processes(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], (t > w+ b)')
        command.exec_cmd(
            'add process d d~ > t t~ [real=QCD]',
            printcmd=False, precmd=True)

        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc)
        self.assertEqual(len(helas['matrix_elements']), 2)
        self.assertEqual(
            [matrix_element.decay_metadata is not None
             for matrix_element in helas['matrix_elements']],
            [True, False])
        self.assertTrue(all(
            matrix_element.born_me['color_basis']
            for matrix_element in helas['matrix_elements']))

    def test_color_splitting_decay_is_rejected(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], (t > t g)')
        with self.assertRaisesRegex(
                fks_common.FKSProcessError,
                'exactly one colour-carrying child'):
            fks_helas_objects.FKSHelasMultiProcess(
                command._fks_multi_proc)

    def test_flat_singlet_subdecay_uses_generated_topology(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], (t > b j j)')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)

        self.assertEqual(len(helas['matrix_elements']), 1)
        matrix_element = helas['matrix_elements'][0]
        visible_processes = [tuple(
            leg.get('id') for leg in process.get_legs_with_decays())
            for process in matrix_element.born_me['processes']]
        self.assertEqual(visible_processes, [
            (2, -2, 5, 2, -1, -6),
            (2, -2, 5, 4, -3, -6)])
        for real in matrix_element.real_processes:
            self.assertEqual(len(real.matrix_element['processes']), 2)

        metadata = matrix_element.decay_metadata
        self.assertEqual(metadata['forced_species'], [6])
        self.assertEqual(
            [(node['pdg'], node['qcd_order'], node['carrier_leaf'])
             for node in metadata['nodes']],
            [(6, 0, 1)])
        self.assertEqual(metadata['leaves'][0]['pdg'], 5)

        # The generated HELAS topology is t -> b W with W -> j j.  The
        # b is visible leg 3 and carries the production top colour; the
        # two W daughters (legs 4 and 5) must not enter colour links.
        visible_links = [tuple(link['link'])
                         for link in matrix_element.color_links]
        self.assertIn((3, 3), visible_links)
        self.assertTrue(all(
            4 not in link and 5 not in link for link in visible_links))

    def test_grouped_decay_flavours_cover_every_real_subprocess(self):
        command = self.generate(
            'g g > t t~ [real=QCD], (t > b j j)')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)

        self.assertEqual(len(helas['matrix_elements']), 1)
        matrix_element = helas['matrix_elements'][0]
        expected_decays = set([(2, -1), (4, -3)])
        for real in matrix_element.real_processes:
            decays_by_core_process = {}
            for process in real.matrix_element['processes']:
                pdgs = tuple(
                    leg.get('id')
                    for leg in process.get_legs_with_decays())
                core_process = pdgs[:3] + pdgs[5:]
                decays_by_core_process.setdefault(
                    core_process, set()).add(pdgs[3:5])
            self.assertTrue(decays_by_core_process)
            self.assertTrue(all(
                decays == expected_decays
                for decays in decays_by_core_process.values()))

    def test_equivalent_qcd_production_flavours_share_subprocesses(self):
        command = self.generate(
            'p p > t t~ [QCD], t > j j b')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)

        # These three matrix elements become the three P* directories: gg,
        # q q~ and q~ q.  Up- and down-type light flavours differ in electric
        # charge, but their pure-QCD production matrix elements are identical.
        self.assertEqual(len(helas['matrix_elements']), 3)
        initial_state_groups = []
        for matrix_element in helas['matrix_elements']:
            initial_state_groups.append(frozenset(
                tuple(leg.get('id') for leg in
                      process.get_legs_with_decays()[:2])
                for process in matrix_element.born_me['processes']))
        self.assertEqual(set(initial_state_groups), set([
            frozenset([(21, 21)]),
            frozenset([(2, -2), (4, -4), (1, -1), (3, -3)]),
            frozenset([(-2, 2), (-4, 4), (-1, 1), (-3, 3)])]))

    def test_inequivalent_production_couplings_remain_separate(self):
        command = self.generate(
            'p p > z [QCD], z > e+ e-')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)

        # Removing electric charge from the metadata-layout key must not
        # override the ordinary matrix-element tags.  The distinct up- and
        # down-type Z couplings still require separate generated code.
        initial_state_groups = []
        for matrix_element in helas['matrix_elements']:
            initial_state_groups.append(frozenset(
                tuple(leg.get('id') for leg in
                      process.get_legs_with_decays()[:2])
                for process in matrix_element.born_me['processes']))
        self.assertEqual(set(initial_state_groups), set([
            frozenset([(2, -2), (4, -4)]),
            frozenset([(1, -1), (3, -3)]),
            frozenset([(-2, 2), (-4, 4)]),
            frozenset([(-1, 1), (-3, 3)])]))

    def test_interface_and_generation_restrictions(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], (t > w+ b)')
        with self.assertRaisesRegex(
                InvalidCmd, 'can only be exported with "output fNLO"'):
            command.check_output([])
        command.check_output(['fNLO'])
        self.assertEqual(command._export_format, 'fNLO')

        process, _ = command.extract_decay_chain_process(
            'u u~ > t t~ [real=QCD], (t > w+ b)')
        options = copy.copy(command.options)
        options['OLP'] = 'GoSam'
        with self.assertRaisesRegex(InvalidCmd, 'native MadLoop'):
            fks_decay.validate_decay_generation(
                process, options, ['QCD'])
        options['OLP'] = 'MadLoop'
        options['low_mem_multicore_nlo_generation'] = True
        with self.assertRaisesRegex(InvalidCmd, 'serial process generation'):
            fks_decay.validate_decay_generation(
                process, options, ['QCD'])
        options['low_mem_multicore_nlo_generation'] = False
        options['complex_mass_scheme'] = True
        with self.assertRaisesRegex(InvalidCmd, 'complex-mass scheme'):
            fks_decay.validate_decay_generation(
                process, options, ['QCD'])
        options['complex_mass_scheme'] = False
        with self.assertRaisesRegex(InvalidCmd, 'QCD corrections only'):
            fks_decay.validate_decay_generation(
                process, options, ['QED'])
        with self.assertRaisesRegex(InvalidCmd, 'EW Sudakov'):
            fks_decay.validate_decay_generation(
                process, options, ['QCD'], ewsudakov=True)

        massless = MasterCmd()
        massless.exec_cmd('import model loop_sm',
                          printcmd=False, precmd=True)
        with self.assertRaisesRegex(InvalidCmd, 'massless decay parent 21'):
            massless.exec_cmd(
                'generate u u~ > z g [real=QCD], (g > u u~)',
                printcmd=False, precmd=True)


if __name__ == '__main__':
    import unittest as unittest_main
    unittest_main.main()
