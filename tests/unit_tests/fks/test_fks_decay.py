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
from madgraph.fks import fks_product
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
            '# LO_DECAY_WIDTH entries are LO physical total widths in GeV.\n'
            '# NLO_DECAY_WIDTH entries are NLO physical total widths in GeV.\n'
            '# Bundled NLO results use the strict O(alpha_s) width expansion.\n'
            '# All NWA denominators use LO widths; NLO-LO enters only linearly.\n'
            '# DECAY_REN_SCALE entries are independent decay scales in GeV.\n'
            'FORMAT 3\n'
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
        nlo_card_text = fks_decay.decay_card_text(
            {6: 1.4915}, {6: 173.0}, nlo_width_pdgs={6},
            nlo_widths={6: 1.3646})
        self.assertIn('FORMAT 4\n', nlo_card_text)
        self.assertIn(
            'LO_DECAY_WIDTH 6 1.4915000000000000e+00\n',
            nlo_card_text)
        self.assertIn(
            'NLO_DECAY_WIDTH 6 1.3646000000000000e+00\n',
            nlo_card_text)
        self.assertNotIn('\nDECAY_WIDTH 6 ', nlo_card_text)
        varied_card_text = fks_decay.decay_card_text(
            {6: 1.4915}, {6: 173.0}, nlo_width_pdgs={6},
            nlo_widths={6: 1.3646},
            decay_scale_variation_mode='independent',
            decay_scale_factors=(1.0, 0.5, 2.0),
            lo_width_variations={
                (6, 0.5): 1.42, (6, 2.0): 1.55},
            nlo_width_variations={
                (6, 0.5): 1.29, (6, 2.0): 1.43})
        self.assertIn('FORMAT 5\n', varied_card_text)
        self.assertIn(
            'DECAY_SCALE_VARIATION_MODE INDEPENDENT\n',
            varied_card_text)
        self.assertIn(
            'DECAY_SCALE_FACTORS 3 1.0000000000000000e+00 '
            '5.0000000000000000e-01 2.0000000000000000e+00\n',
            varied_card_text)
        self.assertIn(
            'LO_DECAY_WIDTH_VARIATION 6 5.0000000000000000e-01 '
            '1.4199999999999999e+00\n', varied_card_text)
        self.assertIn(
            'NLO_DECAY_WIDTH_VARIATION 6 2.0000000000000000e+00 '
            '1.4299999999999999e+00\n', varied_card_text)

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

    def test_nlo_decay_owns_fks_and_is_glued_to_lo_production(self):
        command = self.generate(
            'u u~ > t t~, '
            '(t > w+ b QED^2=2 QCD^2=0 [real=QCD])')
        fks_multi = command._fks_multi_proc

        self.assertTrue(fks_multi.nlo_decay_prototype)
        self.assertEqual(fks_multi.nlo_decay_selector, (6, 1))
        self.assertEqual(fks_multi.nlo_decay_mode, 'real')
        self.assertEqual(len(fks_multi['born_processes']), 1)

        # The FKS family belongs to the decay before HELAS composition.
        fks_process = fks_multi['born_processes'][0]
        self.assertEqual(
            [leg.get('id') for leg in
             fks_process.born_amp['process']['legs']],
            [6, 5, 24])
        self.assertEqual(len(fks_process.real_amps), 1)
        self.assertEqual(
            [leg.get('id') for leg in
             fks_process.real_amps[0].process['legs']],
            [6, 5, 24, 21])
        self.assertEqual(
            [(info['i'], info['j'], info['ij'])
             for info in fks_process.real_amps[0].fks_infos],
            [(4, 2, 2)])

        helas = fks_helas_objects.FKSHelasMultiProcess(
            fks_multi, loop_optimized=False)
        self.assertEqual(len(helas['matrix_elements']), 1)
        matrix_element = helas['matrix_elements'][0]

        self.assertEqual(
            [leg.get('id') for leg in matrix_element.born_me[
                'processes'][0].get_legs_with_decays()],
            [2, -2, 5, 24, -6])
        self.assertEqual(
            [leg.get('id') for leg in matrix_element.real_processes[0].
             matrix_element['processes'][0].get_legs_with_decays()],
            [2, -2, 5, 24, 21, -6])
        self.assertEqual(
            matrix_element.born_me.get_nexternal_ninitial(), (5, 2))
        self.assertEqual(
            matrix_element.real_processes[0].matrix_element.
            get_nexternal_ninitial(), (6, 2))
        for component in [
                matrix_element.born_me,
                matrix_element.real_processes[0].matrix_element]:
            top_wavefunctions = [
                wavefunction for wavefunction in
                component.get_all_wavefunctions()
                if abs(wavefunction.get('pdg_code')) == 6]
            self.assertTrue(top_wavefunctions)
            self.assertEqual(
                set(wavefunction.get('decay_node_id')
                    for wavefunction in top_wavefunctions),
                set([0, 1]))
            for wavefunction in top_wavefunctions:
                expected = (
                    'FNLO_DECAY_DUMMY_WIDTH_RATIO()*mdl_MT'
                    if wavefunction.get('decay_node_id') else 'ZERO')
                self.assertEqual(wavefunction.get('width'), expected)

        metadata = matrix_element.nlo_decay_metadata
        self.assertEqual(metadata['status'], 'INTEGRATION_READY')
        self.assertEqual(metadata['format'], 5)
        self.assertEqual(metadata['production_born_qcd_order'], 4)
        self.assertEqual(metadata['decay_born_qcd_order'], 0)
        self.assertEqual(metadata['parent_pdg'], 6)
        self.assertEqual(len(metadata['fks_maps']), 1)
        fks_mapping = metadata['fks_maps'][0]
        self.assertEqual({
            key: fks_mapping[key]
            for key in ['configuration', 'real_context', 'i', 'j', 'ij']}, {
            'configuration': 1,
            'real_context': 2,
            'i': 4,
            'j': 2,
            'ij': 2})
        self.assertEqual(fks_mapping['targets'], {
            'i': ('LEG', 5),
            'j': ('LEG', 3),
            'ij': ('LEG', 3)})
        self.assertEqual(fks_mapping['partners'], [
            {'local': 1, 'kind': 'NODE', 'target': 1},
            {'local': 2, 'kind': 'LEG', 'target': 3}])
        self.assertEqual(fks_mapping['real_to_born'], {
            1: 1, 2: 2, 3: 3})
        self.assertEqual(metadata['contexts'][0]['production_map'], {
            1: ('LEG', 1), 2: ('LEG', 2),
            3: ('NODE', 1), 4: ('LEG', 5)})
        self.assertEqual(metadata['contexts'][1]['production_map'], {
            1: ('LEG', 1), 2: ('LEG', 2),
            3: ('NODE', 1), 4: ('LEG', 6)})

        projected = matrix_element.get_fks_info_list()[0]
        self.assertEqual(
            [matrix_element.real_processes[0].fks_infos[0][key]
             for key in ['i', 'j', 'ij']],
            [4, 2, 2])
        self.assertEqual(
            [projected['fks_info'][key] for key in ['i', 'j', 'ij']],
            [5, 3, 3])
        self.assertEqual(
            [projected['local_fks_info'][key]
             for key in ['i', 'j', 'ij']],
            [4, 2, 2])
        self.assertEqual(projected['pdgs'], [2, -2, 5, 24, 21, -6])
        self.assertEqual(projected['colors'], [3, -3, 3, 1, 8, -3])
        self.assertEqual(
            projected['massless'],
            [True, True, False, False, True, False])
        self.assertEqual(projected['fks_j_from_i'][5], [3])
        self.assertFalse(projected['ij_massless'])
        self.assertEqual(
            [tuple(link['link']) for link in matrix_element.color_links],
            [(3, 3)])
        self.assertEqual(metadata['color_links'], [
            {'local_first': 1, 'local_second': 1,
             'visible_first': 3, 'visible_second': 3,
             'generated_index': 1},
            {'local_first': 1, 'local_second': 2,
             'visible_first': 3, 'visible_second': 3,
             'generated_index': 1},
            {'local_first': 2, 'local_second': 2,
             'visible_first': 3, 'visible_second': 3,
             'generated_index': 1}])
        info = fks_decay.nlo_decay_info_text(metadata)
        self.assertIn('FORMAT 5\n', info)
        self.assertIn('FORCED_SPECIES 1 6\n', info)
        self.assertIn('TOPOLOGY 1 2 1\n', info)
        self.assertIn('NODE 1 0 6 0 2 2 LEAF 1 LEAF 2\n', info)
        self.assertIn('QCD_ORDERS 4 0\n', info)
        self.assertIn('COUNTS 2 1 1 2\n', info)
        self.assertIn('PRODUCTION_MAP 1 3 NODE 1\n', info)
        self.assertIn('PRODUCTION_MAP 2 4 LEG 6\n', info)
        self.assertIn('LOCAL_MAP 2 1 NODE 1\n', info)
        self.assertIn('LOCAL_MAP 2 4 LEG 5\n', info)
        self.assertIn('FKS_MAP 1 2 4 2 2\n', info)
        self.assertIn('FKS_TARGET 1 I 4 LEG 5\n', info)
        self.assertIn('FKS_TARGET 1 J 2 LEG 3\n', info)
        self.assertIn('FKS_TARGET 1 IJ 2 LEG 3\n', info)
        self.assertIn('FKS_PARTNER 1 1 NODE 1\n', info)
        self.assertIn('FKS_PARTNER 1 2 LEG 3\n', info)
        self.assertIn('REAL_BORN_MAP 1 1 1\n', info)
        self.assertIn('REAL_BORN_MAP 1 2 2\n', info)
        self.assertIn('REAL_BORN_MAP 1 3 3\n', info)
        self.assertIn('COLOR_LINK 1 1 3 3 1\n', info)
        self.assertIn('COLOR_LINK 1 2 3 3 1\n', info)
        self.assertIn('COLOR_LINK 2 2 3 3 1\n', info)

    def test_nlo_decay_groups_multiparticle_production_subprocesses(self):
        command = self.generate(
            'p p > t t~, t > w+ b [QCD]')
        fks_multi = command._fks_multi_proc

        concrete_initials = [
            tuple(amplitude.get('process').get_initial_ids())
            for amplitude in fks_multi.nlo_decay_production_amplitudes]
        self.assertEqual(concrete_initials, [
            (21, 21), (2, -2), (4, -4), (1, -1), (3, -3),
            (-2, 2), (-4, 4), (-1, 1), (-3, 3)])

        helas = fks_helas_objects.FKSHelasMultiProcess(
            fks_multi, loop_optimized=False)
        self.assertEqual(len(helas['matrix_elements']), 3)
        groups = {}
        for matrix_element in helas['matrix_elements']:
            born_initials = frozenset(
                tuple(process.get_initial_ids())
                for process in matrix_element.born_me.get('processes'))
            groups[born_initials] = matrix_element

            real_initials = frozenset(
                tuple(process.get_initial_ids())
                for process in matrix_element.real_processes[0].
                matrix_element.get('processes'))
            virtual_initials = frozenset(
                tuple(process.get_initial_ids())
                for process in matrix_element.virt_matrix_element.
                get('processes'))
            self.assertEqual(real_initials, born_initials)
            self.assertEqual(virtual_initials, born_initials)

            metadata_initial = tuple(
                leg['pdg']
                for leg in matrix_element.nlo_decay_metadata[
                    'production_legs']
                if leg['state'] == 'I')
            self.assertIn(metadata_initial, born_initials)

        expected_groups = set([
            frozenset([(21, 21)]),
            frozenset([(2, -2), (4, -4), (1, -1), (3, -3)]),
            frozenset([(-2, 2), (-4, 4), (-1, 1), (-3, 3)])])
        self.assertEqual(set(groups), expected_groups)
        self.assertEqual(
            groups[frozenset([(21, 21)])].nlo_decay_metadata[
                'virtual_current_count'], 3)
        for initial_states, matrix_element in groups.items():
            if initial_states == frozenset([(21, 21)]):
                continue
            self.assertEqual(
                matrix_element.nlo_decay_metadata['virtual_current_count'], 1)

    def test_nlo_decay_exports_grouped_production_subprocesses(self):
        command = self.generate(
            'p p > t t~, t > w+ b [real=QCD]')

        with tempfile.TemporaryDirectory() as output_dir:
            process_dir = os.path.join(output_dir, 'PROC')
            command.exec_cmd(
                'output fNLO %s' % process_dir,
                printcmd=False, precmd=True)
            subprocess_root = os.path.join(process_dir, 'SubProcesses')
            subprocesses = [
                os.path.join(subprocess_root, name)
                for name in os.listdir(subprocess_root)
                if name.startswith('P') and
                os.path.isdir(os.path.join(subprocess_root, name))]
            self.assertEqual(len(subprocesses), 3)

            grouped_process_counts = []
            for subprocess_dir in subprocesses:
                self.assertTrue(os.path.isfile(os.path.join(
                    subprocess_dir, 'nlo_decay_info.dat')))
                self.assertTrue(os.path.isfile(os.path.join(
                    subprocess_dir, 'matrix_1.f')))
                with open(os.path.join(
                        subprocess_dir, 'born.f')) as stream:
                    process_lines = set(
                        line.strip() for line in stream
                        if line.strip().startswith('C     Process:'))
                grouped_process_counts.append(len(process_lines))
            self.assertEqual(sorted(grouped_process_counts), [1, 4, 4])

    def test_nlo_decay_virtual_is_composed_with_lo_production(self):
        command = self.generate(
            'u u~ > t t~, '
            '(t > w+ b QED^2=2 QCD^2=0 [QCD])')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)
        matrix_element = helas['matrix_elements'][0]
        virtual = matrix_element.virt_matrix_element

        self.assertIsNotNone(virtual)
        self.assertIsNone(matrix_element.nlo_decay_virtual_matrix_element)
        self.assertIn(virtual, helas.get_virt_matrix_elements())
        self.assertTrue(helas['has_loops'])
        self.assertTrue(matrix_element.nlo_decay_metadata['has_virtual'])
        self.assertEqual(
            matrix_element.nlo_decay_metadata['virtual_composition'],
            'CROSSED_PRODUCTION_CURRENT')
        self.assertEqual(
            matrix_element.nlo_decay_metadata['virtual_current_count'], 1)
        self.assertEqual(
            [leg.get('id') for leg in
             virtual['processes'][0].get_legs_with_decays()],
            [2, -2, 5, 24, -6])
        self.assertEqual(virtual.get_nexternal_ninitial(), (5, 2))
        self.assertTrue(virtual.get_loop_diagrams())
        self.assertTrue(virtual.get_born_diagrams())
        self.assertTrue(virtual.get('born_color_basis'))
        self.assertTrue(virtual.get('loop_color_basis'))
        self.assertEqual(
            [entry[0] for entry in virtual.get_split_orders_mapping()[0]],
            [(6, 2)])

        external = sorted(set(
            (wavefunction.get('number_external'),
             wavefunction.get('leg_state'))
            for wavefunction in fks_decay._all_wavefunctions(virtual)
            if (not wavefunction.get('mothers') and
                not wavefunction.get('is_loop'))))
        self.assertEqual(external, [
            (1, False), (2, False), (3, True),
            (4, True), (5, True)])
        initial_flow = set(
            (wavefunction.get('number_external'),
             wavefunction.get('is_part'))
            for wavefunction in fks_decay._all_wavefunctions(virtual)
            if (not wavefunction.get('mothers') and
                not wavefunction.get('is_loop') and
                not wavefunction.get('leg_state')))
        self.assertEqual(initial_flow, set([(1, True), (2, False)]))
        connectors = [
            wavefunction for wavefunction in
            fks_decay._all_wavefunctions(virtual)
            if wavefunction.get('decay_node_id') == 1]
        self.assertTrue(connectors)
        for connector in connectors:
            self.assertTrue(connector.get('mothers'))
            self.assertEqual(
                connector.get('width'),
                'FNLO_DECAY_DUMMY_WIDTH_RATIO()*mdl_MT')

    def test_nlo_decay_virtual_composes_in_optimized_representation(self):
        command = self.generate(
            'u u~ > t t~, '
            '(t > w+ b QED^2=2 QCD^2=0 [QCD])')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=True)
        virtual = helas['matrix_elements'][0].virt_matrix_element

        self.assertTrue(virtual.optimized_output)
        self.assertEqual(virtual.get_nexternal_ninitial(), (5, 2))
        self.assertTrue(virtual.get_loop_diagrams())
        self.assertTrue(virtual.get_born_diagrams())
        self.assertTrue(virtual.get('born_color_basis'))
        self.assertTrue(virtual.get('loop_color_basis'))

    def test_nlo_decay_virtual_composes_multi_diagram_production(self):
        for optimized in [False, True]:
            with self.subTest(optimized=optimized):
                command = self.generate(
                    'g g > t t~, '
                    '(t > w+ b QED^2=2 QCD^2=0 [QCD])')
                matrix_element = fks_helas_objects.FKSHelasMultiProcess(
                    command._fks_multi_proc,
                    loop_optimized=optimized)['matrix_elements'][0]
                virtual = matrix_element.virt_matrix_element

                self.assertEqual(
                    matrix_element.nlo_decay_metadata[
                        'virtual_current_count'], 3)
                self.assertEqual(
                    [leg.get('id') for leg in
                     virtual['processes'][0].get_legs_with_decays()],
                    [21, 21, 5, 24, -6])
                self.assertEqual(virtual.get_nexternal_ninitial(), (5, 2))
                self.assertEqual(len(virtual.get_born_diagrams()), 3)
                self.assertEqual(len(virtual.get_loop_diagrams()), 3)
                external_states = set(
                    (wavefunction.get('number_external'),
                     wavefunction.get('state'),
                     wavefunction.get('leg_state'))
                    for wavefunction in
                    fks_decay._all_wavefunctions(virtual)
                    if (not wavefunction.get('mothers') and
                        not wavefunction.get('is_loop')))
                self.assertIn((1, 'initial', False), external_states)
                self.assertIn((2, 'initial', False), external_states)
                self.assertEqual(
                    set(virtual.get('born_color_basis')),
                    set(matrix_element.born_me.get('color_basis')))
                self.assertEqual(
                    set(virtual.get('loop_color_basis')),
                    set(matrix_element.born_me.get('color_basis')))
                for diagram in virtual.get_loop_diagrams():
                    for amplitude in diagram.get_loop_amplitudes():
                        self.assertTrue(any(
                            mother.get('decay_node_id') == 1
                            for mother in amplitude.get('mothers')))

    def test_nlo_decay_combined_virtual_fortran_is_written(self):
        command = self.generate(
            'g g > t t~, '
            '(t > w+ b QED^2=2 QCD^2=0 [QCD])')
        command.exec_cmd(
            'set loop_optimized_output False',
            printcmd=False, precmd=True)

        with tempfile.TemporaryDirectory() as output_dir:
            process_dir = os.path.join(output_dir, 'PROC')
            command.exec_cmd(
                'output fNLO %s' % process_dir,
                printcmd=False, precmd=True)
            subprocess_root = os.path.join(process_dir, 'SubProcesses')
            subprocesses = [
                os.path.join(subprocess_root, name)
                for name in os.listdir(subprocess_root)
                if name.startswith('P') and
                os.path.isdir(os.path.join(subprocess_root, name))]
            self.assertEqual(len(subprocesses), 1)
            subprocess_dir = subprocesses[0]
            virtuals = [
                os.path.join(subprocess_dir, name)
                for name in os.listdir(subprocess_dir)
                if name.startswith('V') and
                os.path.isdir(os.path.join(subprocess_dir, name))]
            self.assertEqual(len(virtuals), 1)
            self.assertFalse(os.path.exists(os.path.join(
                subprocess_dir, 'NLODecayVirtual')))

            virtual_dir = virtuals[0]
            with open(os.path.join(virtual_dir, 'loop_matrix.f')) as stream:
                loop_source = stream.read()
            with open(os.path.join(virtual_dir, 'born_matrix.f')) as stream:
                loop_born_source = stream.read()
            for source in [loop_source, loop_born_source]:
                self.assertIn('P(0,5)', source)
                self.assertIn(
                    'FNLO_DECAY_DUMMY_WIDTH_RATIO()*MDL_MT', source)
            self.assertIn('CALL FFV2_0', loop_source)
            self.assertIn(
                'DOUBLE PRECISION FNLO_DECAY_DUMMY_WIDTH_RATIO',
                loop_source)
            self.assertIn('CALL VVV1P0_1', loop_source)
            self.assertGreaterEqual(loop_source.count(
                'FNLO_DECAY_DUMMY_WIDTH_RATIO()*MDL_MT'), 3)
            self.assertGreaterEqual(loop_born_source.count(
                'FNLO_DECAY_DUMMY_WIDTH_RATIO()*MDL_MT'), 3)
            self.assertTrue(os.path.isfile(os.path.join(
                virtual_dir, 'global_specs.inc')))

            with open(os.path.join(
                    subprocess_dir, 'nlo_decay_info.dat')) as stream:
                metadata = stream.read()
            self.assertIn(
                'VIRTUAL_COMPOSITION CROSSED_PRODUCTION_CURRENT\n',
                metadata)
            self.assertIn('VIRTUAL_CURRENT_COUNT 3\n', metadata)

    def test_nlo_decay_combined_fortran_matrix_elements_are_written(self):
        command = self.generate(
            'u u~ > t t~, '
            '(t > w+ b QED^2=2 QCD^2=0 [real=QCD])')

        with tempfile.TemporaryDirectory() as output_dir:
            process_dir = os.path.join(output_dir, 'PROC')
            command.exec_cmd(
                'output fNLO %s' % process_dir,
                printcmd=False, precmd=True)
            subprocess_root = os.path.join(process_dir, 'SubProcesses')
            subprocesses = [
                name for name in os.listdir(subprocess_root)
                if name.startswith('P') and
                os.path.isdir(os.path.join(subprocess_root, name))]
            self.assertEqual(len(subprocesses), 1)
            subprocess_dir = os.path.join(
                subprocess_root, subprocesses[0])
            born_path = os.path.join(subprocess_dir, 'born.f')
            real_path = os.path.join(subprocess_dir, 'matrix_1.f')
            fks_info_path = os.path.join(subprocess_dir, 'fks_info.inc')

            with open(born_path) as stream:
                born_source = stream.read()
            with open(real_path) as stream:
                real_source = stream.read()
            with open(fks_info_path) as stream:
                fks_info_source = stream.read()
            flat_fks_info = ' '.join(
                fks_info_source.replace('$', ' ').split()).replace(' ,', ',')
            self.assertIn('SUBROUTINE SBORN(P,ANS_SUMMED)', born_source)
            self.assertIn('SUBROUTINE SMATRIX1(P,ANS_SUMMED)', real_source)
            self.assertIn('P(0,5),MDL_MT', born_source)
            self.assertIn('P(0,6),MDL_MT', real_source)
            self.assertEqual(
                born_source.count(
                    'FNLO_DECAY_DUMMY_WIDTH_RATIO()*MDL_MT'), 1)
            self.assertEqual(
                real_source.count(
                    'FNLO_DECAY_DUMMY_WIDTH_RATIO()*MDL_MT'), 2)
            self.assertIn(
                'DOUBLE PRECISION FNLO_DECAY_DUMMY_WIDTH_RATIO',
                born_source)
            self.assertIn(
                'EXTERNAL FNLO_DECAY_DUMMY_WIDTH_RATIO', real_source)
            self.assertIn('DATA FKS_I_D / 5 /', fks_info_source)
            self.assertIn('DATA FKS_J_D / 3 /', fks_info_source)
            self.assertIn(
                'FKS_J_FROM_I_D(1, 5, JPOS), JPOS = 0, 1) / 1, 3 /',
                flat_fks_info)
            self.assertIn(
                'PARTICLE_TYPE_D(1, IPOS), IPOS=1, NEXTERNAL) '
                '/ 3, -3, 3, 1, 8, -3 /',
                flat_fks_info)
            self.assertIn(
                'PDG_TYPE_D(1, IPOS), IPOS=1, NEXTERNAL) '
                '/ 2, -2, 5, 24, 21, -6 /',
                flat_fks_info)
            self.assertIn('DATA IJ_VALUES /0/', born_source)

            self.assertFalse(os.path.exists(os.path.join(
                subprocess_dir, 'NLO_DECAY_SUBTRACTION_INCOMPLETE')))
            with open(os.path.join(
                    subprocess_dir, 'nlo_decay_info.dat')) as stream:
                metadata = stream.read()
            self.assertIn('FORMAT 5\n', metadata)
            self.assertIn('STATUS INTEGRATION_READY\n', metadata)
            self.assertIn('QCD_ORDERS 4 0\n', metadata)
            self.assertIn('PRODUCTION_MAP 1 3 NODE 1\n', metadata)
            self.assertIn('PRODUCTION_MAP 2 4 LEG 6\n', metadata)
            self.assertIn('FKS_TARGET 1 I 4 LEG 5\n', metadata)
            self.assertIn('FKS_PARTNER 1 1 NODE 1\n', metadata)
            self.assertIn('REAL_BORN_MAP 1 2 2\n', metadata)
            for runtime_source in [
                    'nlo_decay_metadata.f90',
                    'nlo_decay_kinematics.f90']:
                self.assertTrue(os.path.isfile(os.path.join(
                    subprocess_dir, runtime_source)))
            with open(os.path.join(
                    subprocess_dir, 'nlo_decay_kinematics.f90')) as stream:
                kinematics = stream.read().lower()
            flat_kinematics = ' '.join(kinematics.split())
            self.assertIn(
                'corrected-parent rest frame', flat_kinematics)
            self.assertIn(
                'no production momentum participates', flat_kinematics)
            with open(os.path.join(
                    subprocess_dir, 'genps_fks.f90')) as stream:
                phase_space = stream.read().lower()
            self.assertIn(
                'generate_nlo_decay_fks_kinematics', phase_space)
            self.assertTrue(os.path.isfile(os.path.join(
                process_dir, 'Cards', 'decay_card.dat')))
            with open(os.path.join(
                    process_dir, 'Cards', 'decay_card.dat')) as stream:
                decay_card = stream.read()
            self.assertIn('FORMAT 4\n', decay_card)
            self.assertIn('LO_DECAY_WIDTH 6 ', decay_card)
            self.assertIn('NLO_DECAY_WIDTH 6 ', decay_card)
            self.assertNotIn('\nDECAY_WIDTH 6 ', decay_card)
            self.assertTrue(os.path.islink(os.path.join(
                subprocess_dir, 'decay_card.dat')))

    def test_nlo_decay_additional_and_nested_decay_trees(self):
        cases = [{
            'process': (
                'u u~ > t t~, '
                '(t > w+ b QED=1 [real=QCD]), (t~ > w- b~)'),
            'born': [2, -2, 5, 24, -24, -5],
            'real': [2, -2, 5, 24, 21, -24, -5],
            'corrected': 6,
            'parents': {-6: 0, 6: 0},
            'forced_species': [6]}, {
            'process': (
                'u u~ > t t~, '
                '(t > w+ b, '
                'w+ > u d~ QED=1 [real=QCD]), (t~ > w- b~)'),
            'born': [2, -2, 2, -1, 5, -24, -5],
            'real': [2, -2, 2, -1, 21, 5, -24, -5],
            'corrected': 24,
            'parents': {-6: 0, 6: 0, 24: 6},
            'forced_species': [6, 24]}, {
            'process': (
                'u u~ > t t~, '
                '(t > w+ b QED=1 [real=QCD], '
                'w+ > u d~), (t~ > w- b~)'),
            'born': [2, -2, 5, 2, -1, -24, -5],
            'real': [2, -2, 5, 2, -1, 21, -24, -5],
            'corrected': 6,
            'parents': {-6: 0, 6: 0, 24: 6},
            'forced_species': [6, 24]}]

        for case in cases:
            with self.subTest(process=case['process']):
                command = self.generate(case['process'])
                helas = fks_helas_objects.FKSHelasMultiProcess(
                    command._fks_multi_proc, loop_optimized=False)
                self.assertEqual(len(helas['matrix_elements']), 1)
                matrix_element = helas['matrix_elements'][0]
                metadata = matrix_element.nlo_decay_metadata

                self.assertEqual([
                    leg.get('id') for leg in matrix_element.born_me[
                        'processes'][0].get_legs_with_decays()], case['born'])
                self.assertEqual([
                    leg.get('id') for leg in matrix_element.real_processes[0].
                    matrix_element['processes'][0].get_legs_with_decays()],
                    case['real'])
                nodes_by_pdg = dict(
                    (node['pdg'], node) for node in metadata['nodes'])
                self.assertEqual(set(nodes_by_pdg), set(case['parents']))
                actual_parents = {}
                for pdg, node in nodes_by_pdg.items():
                    parent = node['parent']
                    actual_parents[pdg] = (
                        metadata['nodes'][parent - 1]['pdg']
                        if parent else 0)
                self.assertEqual(actual_parents, case['parents'])
                self.assertEqual(
                    metadata['corrected_node'],
                    nodes_by_pdg[case['corrected']]['id'])
                self.assertEqual(
                    metadata['forced_species'], case['forced_species'])
                self.assertEqual(metadata['format'], 5)

                for context in metadata['contexts']:
                    production_nodes = [
                        target for kind, target in
                        context['production_map'].values()
                        if kind == 'NODE']
                    self.assertEqual(
                        set(production_nodes),
                        set(node['id'] for node in metadata['nodes']
                            if node['parent'] == 0))
                    self.assertEqual(
                        context['local_map'][1],
                        ('NODE', metadata['corrected_node']))

    def test_nlo_decay_virtual_supports_nested_corrections(self):
        cases = [(
            'u u~ > t t~, '
            '(t > w+ b, w+ > u d~ QED=1 [QCD]), (t~ > w- b~)',
            24), (
            'u u~ > t t~, '
            '(t > w+ b QED=1 [QCD], w+ > u d~), (t~ > w- b~)',
            6)]

        for process, corrected_pdg in cases:
            with self.subTest(process=process):
                command = self.generate(process)
                matrix_element = fks_helas_objects.FKSHelasMultiProcess(
                    command._fks_multi_proc,
                    loop_optimized=False)['matrix_elements'][0]
                virtual = matrix_element.virt_matrix_element
                metadata = matrix_element.nlo_decay_metadata

                self.assertIsNotNone(virtual)
                self.assertTrue(virtual.get_loop_diagrams())
                self.assertTrue(virtual.get_born_diagrams())
                self.assertEqual(
                    virtual.get_nexternal_ninitial(),
                    matrix_element.born_me.get_nexternal_ninitial())
                visible_legs = sorted(
                    matrix_element.born_me['processes'][0].
                    get_legs_with_decays(),
                    key=lambda leg: leg.get('number'))
                model = matrix_element.born_me['processes'][0].get('model')
                self.assertEqual(
                    virtual.get_external_masses()[:-2],
                    [model.get_particle(leg.get('id')).get('mass')
                     for leg in visible_legs])
                self.assertEqual(
                    metadata['nodes'][metadata['corrected_node'] - 1]['pdg'],
                    corrected_pdg)
                self.assertEqual(
                    metadata['virtual_composition'],
                    'CROSSED_PRODUCTION_CURRENT')
                self.assertGreater(metadata['virtual_current_count'], 0)

    def test_nlo_decay_symmetry_is_local_to_the_corrected_decay(self):
        command = self.generate(
            'u u~ > t t~ g, t > w+ b QED=1 [real=QCD]')
        matrix_element = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc,
            loop_optimized=False)['matrix_elements'][0]
        real = matrix_element.real_processes[0].matrix_element
        metadata = matrix_element.nlo_decay_metadata

        born_pdgs = [
            leg.get('id') for leg in matrix_element.born_me[
                'processes'][0].get_legs_with_decays()]
        real_pdgs = [
            leg.get('id') for leg in
            real['processes'][0].get_legs_with_decays()]
        self.assertEqual(born_pdgs.count(21), 1)
        self.assertEqual(real_pdgs.count(21), 2)
        # The production and decay gluons belong to different factorized
        # subprocesses; inserting the decay must not create a spurious 2!.
        self.assertEqual(
            matrix_element.born_me['identical_particle_factor'], 1)
        self.assertEqual(real['identical_particle_factor'], 1)

        real_context = next(
            context for context in metadata['contexts']
            if context['kind'] == 'REAL')
        local_final_pdgs = [
            leg['pdg'] for leg in real_context['local_legs']
            if leg['state'] == 'F']
        self.assertEqual(local_final_pdgs.count(21), 1)
        emitted_target = metadata['fks_maps'][0]['targets']['i']
        self.assertEqual(emitted_target[0], 'LEG')
        self.assertEqual(real_pdgs[emitted_target[1] - 1], 21)
        production_gluons = [
            target for number, (kind, target) in
            real_context['production_map'].items()
            if (kind == 'LEG' and
                metadata['production_legs'][number - 1]['pdg'] == 21)]
        self.assertEqual(len(production_gluons), 1)
        self.assertNotEqual(production_gluons[0], emitted_target[1])

    def test_full_nlo_bundle_combines_production_and_each_decay(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b QED=1 [real=QCD]), '
            '(t~ > w- b~ QED=1 [real=QCD])')
        self.assertTrue(command._fks_multi_proc.full_nlo_decay_bundle)
        self.assertEqual(len(command._fks_multi_proc.members), 3)

        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)
        self.assertEqual(len(helas['matrix_elements']), 1)
        matrix_element = helas['matrix_elements'][0]
        self.assertTrue(matrix_element.contribution_bundle)
        self.assertEqual(
            [entry['kind'] for entry in
             matrix_element.bundle_contributions],
            ['PRODUCTION', 'NLO_DECAY', 'NLO_DECAY'])
        self.assertEqual(
            [entry['parent_pdg'] for entry in
             matrix_element.bundle_contributions],
            [0, 6, -6])
        self.assertEqual(
            [entry['has_virtual'] for entry in
             matrix_element.bundle_contributions],
            [False, False, False])
        self.assertEqual(len(matrix_element.bundle_nlo_decay_metadata), 2)
        self.assertEqual(
            [leg.get('id') for leg in matrix_element.born_me[
                'processes'][0].get_legs_with_decays()],
            [2, -2, 24, 5, -24, -5])

        contributions = matrix_element.bundle_contributions
        self.assertEqual(contributions[0]['first'], 1)
        for previous, current in zip(contributions, contributions[1:]):
            self.assertEqual(current['first'], previous['last'] + 1)
        infos = matrix_element.get_fks_info_list()
        self.assertEqual(len(infos), contributions[-1]['last'])
        for contribution in contributions:
            owned = [info for info in infos
                     if info['contribution'] == contribution['id']]
            self.assertEqual(
                len(owned), contribution['last'] -
                contribution['first'] + 1)

    def test_simultaneous_tree_current_contraction(self):
        """Contract two decay reals, and a production/decay real pair."""

        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b QED=1 [real=QCD]), '
            '(t~ > w- b~ QED=1 [real=QCD])')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)
        bundled = helas['matrix_elements'][0]
        self.assertEqual(
            bundled.factorized_contraction_mode,
            'HELAS_CURRENT_PRODUCT')
        families = bundled.factorized_decay_current_families
        self.assertEqual(len(families), 2)
        self.assertEqual(
            set(family['selector'] for family in families),
            set([(6, 1), (-6, 1)]))
        for family in families:
            self.assertTrue(family['real_currents'])

        def assert_top_connectors(matrix_element):
            connectors = [
                wavefunction for wavefunction in
                matrix_element.get_all_wavefunctions()
                if (abs(wavefunction.get('pdg_code')) == 6 and
                    wavefunction.get('decay_node_id'))]
            self.assertEqual(
                set(wavefunction.get('decay_node_id')
                    for wavefunction in connectors), set([1, 2]))
            for wavefunction in connectors:
                self.assertEqual(
                    wavefunction.get('width'),
                    'FNLO_DECAY_DUMMY_WIDTH_RATIO()*mdl_MT')

        production_family = bundled.factorized_production_core_family
        self.assertEqual(
            [leg.get('id') for leg in production_family[
                'born_amplitude']['process']['legs']],
            [2, -2, 6, -6])
        self.assertTrue(production_family['real_amplitudes'])
        decay_real_components = [{
            'selector': family['selector'],
            'current': family['real_currents'][0],
            'stage': 'DECAY_%d' % index,
            'state': 'REAL',
            'source_index': 1}
            for index, family in enumerate(families, 1)]
        double_decay_real, context, metadata = \
            fks_decay.compose_simultaneous_tree_matrix_element(
                production_family['born_amplitude'], decay_real_components,
                contraction_id=7)
        double_decay_pdgs = [
            leg.get('id') for leg in
            double_decay_real['processes'][0].get_legs_with_decays()]
        self.assertEqual(
            double_decay_real.get_nexternal_ninitial(), (8, 2))
        self.assertEqual(double_decay_pdgs.count(21), 2)
        self.assertEqual(context['kind'], 'COMPOSITE')
        self.assertEqual(context['source_index'], 7)
        self.assertEqual(
            double_decay_real.fnlo_simultaneous_contraction,
            (('DECAY_1', 'REAL', 1, families[0]['selector']),
             ('DECAY_2', 'REAL', 1, families[1]['selector'])))
        self.assertEqual(
            [component['root_node_id'] for component in
             metadata['simultaneous_components']], [1, 2])
        assert_top_connectors(double_decay_real)
        self.assertTrue(double_decay_real.get('color_basis'))
        self.assertTrue(double_decay_real.get('color_matrix'))

        production_real = next(
            amplitude for amplitude in
            production_family['real_amplitudes']
            if amplitude.get('diagrams'))
        production_decay_components = []
        for index, family in enumerate(families, 1):
            production_decay_components.append({
                'selector': family['selector'],
                'current': (family['real_currents'][0] if index == 1
                            else family['born_current']),
                'stage': 'DECAY_%d' % index,
                'state': 'REAL' if index == 1 else 'BORN',
                'source_index': 1})
        production_decay_real, _, _ = \
            fks_decay.compose_simultaneous_tree_matrix_element(
                production_real, production_decay_components,
                contraction_id=8)
        production_decay_pdgs = [
            leg.get('id') for leg in
            production_decay_real['processes'][0].get_legs_with_decays()]
        self.assertEqual(
            production_decay_real.get_nexternal_ninitial(), (8, 2))
        self.assertEqual(production_decay_pdgs.count(21), 2)
        assert_top_connectors(production_decay_real)

        catalog = bundled.factorized_product_catalog
        self.assertEqual(bundled.factorized_product_mode,
                         'STAGEWISE_NLO_PRODUCT')
        self.assertEqual(
            [(stage.label, len(stage.choices))
             for stage in catalog.stages],
            [('PRODUCTION', 7), ('DECAY_1', 2), ('DECAY_2', 2)])
        self.assertEqual(len(catalog), 28)
        self.assertEqual(catalog.first_order_sector_count, 8)

        triple_reals = [
            sector for sector in catalog.iter_sectors()
            if (sector.states == (fks_product.REAL,) * 3 and
                sector.choices[0].source_index == 1)]
        self.assertEqual(len(triple_reals), 4)
        sector = triple_reals[0]
        self.assertEqual(sector.perturbative_order, 3)
        self.assertEqual(sector.real_order, 3)
        self.assertEqual(sector.counterevent_count, 16)
        counterevents = dict(
            (event.codes, event) for event in sector.iter_counterevents())
        for codes in [('R', 'R', 'R'), ('S', 'R', 'R'),
                      ('R', 'S', 'R'), ('R', 'R', 'S'),
                      ('S', 'S', 'S')]:
            self.assertIn(codes, counterevents)
        self.assertEqual(counterevents[('R', 'R', 'R')].inclusion_sign, 1)
        self.assertEqual(counterevents[('S', 'S', 'S')].inclusion_sign, -1)

        carrier_expectations = {
            ('R', 'R', 'R'): (9, 3),
            ('S', 'R', 'R'): (8, 2),
            ('R', 'S', 'R'): (8, 2),
            ('R', 'R', 'S'): (8, 2),
            ('S', 'S', 'S'): (6, 0)}
        carriers = {}
        for codes, (nexternal, gluons) in carrier_expectations.items():
            carrier = counterevents[codes].build_tree_matrix_element()
            carriers[codes] = carrier
            self.assertEqual(
                carrier.get_nexternal_ninitial(), (nexternal, 2))
            self.assertEqual([
                leg.get('id') for leg in carrier['processes'][0]
                .get_legs_with_decays()].count(21), gluons)

        # Four production FKS configurations share the first real amplitude.
        # Their sectors must therefore share one coherently contracted ME.
        same_tree_sector = triple_reals[1]
        same_tree_rr = next(
            event for event in same_tree_sector.iter_counterevents()
            if event.codes == ('R', 'R', 'R'))
        self.assertIs(
            same_tree_rr.build_tree_matrix_element(),
            carriers[('R', 'R', 'R')])

        phase_space = sector.phase_space(born_dimension=2)
        self.assertEqual(phase_space.dimension, 11)
        mixed_limit = counterevents[('SC', 'S', 'S')]
        event = phase_space.event(
            mixed_limit,
            (10., 11., .2, .3, .4, .5, .6, .7, .8, .9, 1.))
        self.assertEqual(event.born_coordinates, (10., 11.))
        self.assertEqual(
            [(radiation['slot'], radiation['xi'], radiation['y'])
             for radiation in event.radiation_coordinates],
            [('SOFT_COLLINEAR', 0., 1.),
             ('SOFT', 0., .6),
             ('SOFT', 0., .9)])
        self.assertEqual(event.inclusion_sign, 1)
        self.assertEqual(event.matrix_element.get_nexternal_ninitial(),
                         (6, 2))

        product_info = fks_product.product_info_text(catalog)
        self.assertIn('ENUMERATION CARTESIAN_LAZY\n', product_info)
        self.assertIn('COUNTEREVENTS TENSOR_PRODUCT\n', product_info)
        self.assertIn('STAGES 3\n', product_info)
        self.assertIn('SECTORS 28\n', product_info)
        self.assertEqual(product_info.count('\nSTAGE '), 3)
        self.assertEqual(product_info.count('\nCHOICE '), 11)

        layout = catalog.canonical_layout
        self.assertTrue(layout.supports_all_stages)
        self.assertEqual(layout.initial_count, 2)
        self.assertEqual(layout.base_count, 6)
        self.assertEqual(layout.max_count, 9)
        self.assertEqual(layout.emission_slots, {1: 7, 2: 8, 3: 9})
        self.assertEqual(layout.carrier_count, 16)
        born_carrier = layout.carrier((0, 0, 0))
        self.assertEqual(
            born_carrier.matrix_element.
            fnlo_product_selected_squared_orders,
            ((4, 4),))
        self.assertEqual(
            born_carrier.matrix_element.get('processes')[0].get(
                'squared_orders'), {})
        maximum_carrier = layout.carrier((1, 1, 1))
        self.assertEqual(maximum_carrier.local_count, 9)
        self.assertEqual(set(maximum_carrier.local_to_canonical),
                         set(range(1, 10)))

        # A production-level top colour charge is the coherent sum of its
        # coloured decay descendants.  When that decay is real, its gluon is
        # part of the sum as well.  The incoming top of the decay-local
        # process is crossed and therefore carries the opposite sign.
        top_endpoint = layout.color_endpoint_expansion(
            (0, 1, 1), 1, 1, 3, 3)
        self.assertEqual(len(top_endpoint), 2)
        self.assertEqual(set(sign for _, sign in top_endpoint), set([1]))
        decay_parent_endpoint = layout.color_endpoint_expansion(
            (0, 0, 0), 2, 1, 1, 1)
        self.assertEqual(len(decay_parent_endpoint), 1)
        self.assertEqual(decay_parent_endpoint[0][1], -1)

        layout.prepare_color_insertions()
        self.assertGreater(len(
            layout.carrier((0, 0, 0)).color_insertions), 1)
        runtime_plans = fks_product._product_runtime_plans(layout)
        self.assertEqual(max(
            len(term['eikonals']) for plan in runtime_plans
            for term in plan['soft_terms']), 3)
        layout_info = fks_product.product_layout_info_text(layout)
        self.assertIn('STATUS COMPLETE\n', layout_info)
        self.assertIn('BASE_LEGS 6\n', layout_info)
        self.assertIn('MAX_LEGS 9\n', layout_info)
        self.assertEqual(layout_info.count('\nEMISSION_SLOT '), 3)
        self.assertEqual(layout_info.count('\nCARRIER '), 16)

    def test_product_sector_retains_uncorrected_root_decays(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b QED=1 [real=QCD]), (t~ > w- b~)')
        bundled = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc,
            loop_optimized=False)['matrix_elements'][0]
        catalog = bundled.factorized_product_catalog
        self.assertEqual(len(catalog.baseline_decay_currents), 2)
        sector = next(
            candidate for candidate in catalog.iter_sectors()
            if candidate.states == (fks_product.BORN, fks_product.REAL))
        matrix_element = sector.build_tree_matrix_element()
        pdgs = [
            leg.get('id') for leg in
            matrix_element['processes'][0].get_legs_with_decays()]
        self.assertNotIn(6, pdgs)
        self.assertNotIn(-6, pdgs)
        self.assertEqual(pdgs.count(24), 1)
        self.assertEqual(pdgs.count(-24), 1)
        self.assertEqual(pdgs.count(21), 1)

    def test_full_nlo_bundle_supports_nested_corrected_decays(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b QED=1 [real=QCD], '
            'w+ > u d~ QED=1 [real=QCD]), '
            '(t~ > w- b~)')
        matrix_element = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc,
            loop_optimized=False)['matrix_elements'][0]

        self.assertEqual(
            [entry['parent_pdg'] for entry in
             matrix_element.bundle_contributions],
            [0, 6, 24])
        self.assertEqual(
            [metadata['parent_pdg'] for metadata in
             matrix_element.bundle_nlo_decay_metadata],
            [6, 24])
        self.assertEqual(
            [metadata['nodes'][metadata['corrected_node'] - 1]['pdg']
             for metadata in
             matrix_element.bundle_nlo_decay_metadata],
            [6, 24])

    def test_full_nlo_bundle_expands_identical_decay_occurrences(self):
        command = self.generate(
            'u u~ > z z [real=QCD], '
            '(z > b b~ QED=1 [real=QCD])')
        matrix_element = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc,
            loop_optimized=False)['matrix_elements'][0]

        # The production Born retains its 1/2! identical-Z symmetry divisor,
        # while the decay correction is represented once for each physical Z
        # occurrence.  The two terms therefore form the required sum over
        # resonance assignments without an extra combinatorial multiplier.
        self.assertEqual(
            matrix_element.born_me.get('identical_particle_factor'), 2)
        self.assertEqual(
            [entry['kind'] for entry in
             matrix_element.bundle_contributions],
            ['PRODUCTION', 'NLO_DECAY', 'NLO_DECAY'])
        self.assertEqual(
            [entry['parent_pdg'] for entry in
             matrix_element.bundle_contributions],
            [0, 23, 23])
        self.assertEqual(
            [(metadata['corrected_node'],
              metadata['parent_occurrence'])
             for metadata in
             matrix_element.bundle_nlo_decay_metadata],
            [(1, 1), (2, 2)])
        self.assertEqual(
            [metadata['nodes'][metadata['corrected_node'] - 1]['pdg']
             for metadata in
             matrix_element.bundle_nlo_decay_metadata],
            [23, 23])
        self.assertEqual(
            set(real.matrix_element.get('identical_particle_factor')
                for real in matrix_element.real_processes),
            set([2]))

    def test_full_nlo_bundle_distinguishes_identical_decay_modes(self):
        command = self.generate(
            'u u~ > z z [real=QCD], '
            '(z > b b~ QED=1 [real=QCD]), '
            '(z > e+ e- QED=2)')
        matrix_element = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc,
            loop_optimized=False)['matrix_elements'][0]

        # Only the Z assigned to b b~ owns a numerator correction.  Both Z
        # nodes remain in the production topology, so the runtime width
        # expansion can still apply delta Gamma_Z/Gamma_Z to both physical
        # denominators.
        self.assertEqual(
            [entry['kind'] for entry in
             matrix_element.bundle_contributions],
            ['PRODUCTION', 'NLO_DECAY'])
        self.assertEqual(
            matrix_element.bundle_nlo_decay_metadata[0]
            ['parent_occurrence'], 1)
        self.assertEqual(
            [node['pdg'] for node in matrix_element.decay_metadata['nodes']],
            [23, 23])
        born_factor = matrix_element.born_me.get(
            'identical_particle_factor')
        self.assertEqual(
            set(real.matrix_element.get('identical_particle_factor')
                for real in matrix_element.real_processes),
            set([born_factor]))

    def test_full_nlo_bundle_groups_production_subprocesses(self):
        command = self.generate(
            'p p > t t~ [real=QCD], '
            't > w+ b QED=1 [real=QCD]')
        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)
        self.assertEqual(len(helas['matrix_elements']), 3)
        expected_groups = set([
            frozenset([(21, 21)]),
            frozenset([(2, -2), (4, -4), (1, -1), (3, -3)]),
            frozenset([(-2, 2), (-4, 4), (-1, 1), (-3, 3)])])
        actual_groups = set()
        for matrix_element in helas['matrix_elements']:
            actual_groups.add(frozenset(
                tuple(process.get_initial_ids())
                for process in matrix_element.born_me.get('processes')))
            self.assertEqual(
                [entry['kind'] for entry in
                 matrix_element.bundle_contributions],
                ['PRODUCTION', 'NLO_DECAY'])
        self.assertEqual(actual_groups, expected_groups)

    def test_full_nlo_bundle_supports_additional_production_processes(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b QED=1 [real=QCD]), '
            '(t~ > w- b~ QED=1 [real=QCD])')
        command.exec_cmd(
            'add process d d~ > t t~ [real=QCD], '
            '(t > w+ b QED=1 [real=QCD]), '
            '(t~ > w- b~ QED=1 [real=QCD])',
            printcmd=False, precmd=True)

        self.assertTrue(command._fks_multi_proc.full_nlo_decay_bundle)
        self.assertEqual(
            len(command._fks_multi_proc.production['born_processes']), 2)
        self.assertTrue(all(
            len(member.nlo_decay_production_amplitudes) == 2
            for member in command._fks_multi_proc.decays))

        helas = fks_helas_objects.FKSHelasMultiProcess(
            command._fks_multi_proc, loop_optimized=False)
        self.assertEqual(len(helas['matrix_elements']), 1)
        matrix_element = helas['matrix_elements'][0]
        self.assertEqual(set(
            tuple(process.get_initial_ids())
            for process in matrix_element.born_me.get('processes')),
            set([(2, -2), (1, -1)]))
        self.assertEqual(
            [entry['kind'] for entry in
             matrix_element.bundle_contributions],
            ['PRODUCTION', 'NLO_DECAY', 'NLO_DECAY'])

    def test_full_nlo_bundle_rejects_mismatched_added_decay_trees(self):
        command = self.generate(
            'u u~ > t t~ [real=QCD], '
            '(t > w+ b QED=1 [real=QCD], w+ > u d~), '
            '(t~ > w- b~)')
        with self.assertRaisesRegex(
                InvalidCmd, 'same corrected decay definitions'):
            command.exec_cmd(
                'add process d d~ > t t~ [real=QCD], '
                '(t > w+ b QED=1 [real=QCD], w+ > c s~), '
                '(t~ > w- b~)',
                printcmd=False, precmd=True)

    def test_full_nlo_bundle_exports_virtual_dispatchers(self):
        command = self.generate(
            'u u~ > t t~ [QCD], '
            '(t > w+ b QED=1 [QCD]), '
            '(t~ > w- b~ QED=1 [QCD])')
        command.exec_cmd(
            'set loop_optimized_output False',
            printcmd=False, precmd=True)

        with tempfile.TemporaryDirectory() as output_dir:
            process_dir = os.path.join(output_dir, 'PROC')
            command.exec_cmd(
                'output fNLO %s' % process_dir,
                printcmd=False, precmd=True)
            subprocess_root = os.path.join(process_dir, 'SubProcesses')
            subprocesses = [
                os.path.join(subprocess_root, name)
                for name in os.listdir(subprocess_root)
                if name.startswith('P') and
                os.path.isdir(os.path.join(subprocess_root, name))]
            self.assertEqual(len(subprocesses), 1)
            subprocess_dir = subprocesses[0]
            self.assertEqual(sorted(
                name for name in os.listdir(subprocess_dir)
                if name.startswith('VContribution') and
                os.path.isdir(os.path.join(subprocess_dir, name))),
                ['VContribution1', 'VContribution2', 'VContribution3'])
            with open(os.path.join(
                    subprocess_dir,
                    'virtual_contribution_chooser.f')) as stream:
                chooser = ' '.join(
                    stream.read().replace('$', ' ').split())
            for contribution in [1, 2, 3]:
                self.assertIn(
                    'CALL FNLOC%d_SLOOPMATRIX_THRES(P, RAW_ANS' %
                    contribution,
                    chooser)
                self.assertIn(
                    'CALL FNLOC%d_FORCE_STABILITY_CHECK' % contribution,
                    chooser)
            self.assertIn('ANS(0:3,1) = RAW_ANS(0:3,0)', chooser)
            self.assertNotIn('FNLOC1_GETORDPOWFROMINDEX_ML5', chooser)
            with open(os.path.join(
                    subprocess_dir, 'virtual_libraries.inc')) as stream:
                libraries = stream.read()
            self.assertIn(
                'libMadLoop_1.a libMadLoop_2.a libMadLoop_3.a',
                libraries)
            with open(os.path.join(
                    subprocess_dir,
                    'nlo_contribution_info.dat')) as stream:
                metadata = stream.read()
            self.assertIn('FORMAT 3\n', metadata)
            self.assertIn('COUNT 3\n', metadata)
            self.assertIn('VIRTUAL_GRIDS 3\n', metadata)
            self.assertEqual(metadata.count('\nVIRTUAL_GRID '), 3)
            with open(os.path.join(
                    subprocess_dir,
                    'multiplicative_product_info.dat')) as stream:
                product_metadata = stream.read()
            self.assertIn('PRESCRIPTION STAGEWISE_NLO_PRODUCT\n',
                          product_metadata)
            self.assertIn('ENUMERATION CARTESIAN_LAZY\n',
                          product_metadata)
            self.assertIn('COUNTEREVENTS TENSOR_PRODUCT\n',
                          product_metadata)
            self.assertIn('STAGES 3\n', product_metadata)
            self.assertEqual(product_metadata.count(' FINITE '), 3)
            self.assertEqual(product_metadata.count('\nVIRTUAL_ORDER '), 3)
            self.assertTrue(os.path.isfile(os.path.join(
                subprocess_dir, 'multiplicative_product.f90')))
            self.assertTrue(os.path.isfile(os.path.join(
                subprocess_dir, 'multiplicative_product_kinematics.f90')))
            with open(os.path.join(
                    subprocess_dir,
                    'multiplicative_product_layout.dat')) as stream:
                product_layout = stream.read()
            self.assertIn('STATUS COMPLETE\n', product_layout)
            self.assertIn('INITIAL_LEGS 2\n', product_layout)
            self.assertIn('BASE_LEGS 6\n', product_layout)
            self.assertIn('MAX_LEGS 9\n', product_layout)
            self.assertEqual(product_layout.count('\nEMISSION_SLOT '), 3)
            self.assertEqual(product_layout.count('\nCARRIER '), 16)
            for prefix in ['product_carrier_',
                           'product_carrier_amplitudes_',
                           'product_carrier_contraction_']:
                self.assertEqual(len([
                    name for name in os.listdir(subprocess_dir)
                    if name.startswith(prefix) and name.endswith('.f') and
                    name[len(prefix):-2].isdigit()]),
                    16)
            with open(os.path.join(
                    subprocess_dir, 'product_carrier_001.f')) as stream:
                ordinary_carrier = ' '.join(stream.read().split())
            self.assertIn('DATA CHOSEN_SO_CONFIGS/.TRUE./',
                          ordinary_carrier)
            with open(os.path.join(
                    subprocess_dir,
                    'product_carrier_contraction_001.f')) as stream:
                correlated_carrier = ' '.join(stream.read().split())
            self.assertIn(
                'DATA (KEEP_AMP_PAIR_FLAT(I),I=1,1) /.TRUE./',
                correlated_carrier)
            with open(os.path.join(
                    subprocess_dir,
                    'multiplicative_product_generated.f90')) as stream:
                generated_product = stream.read()
            self.assertIn('module multiplicative_product_generated',
                          generated_product)
            self.assertIn('subroutine generated_product_mapper',
                          generated_product)
            self.assertIn('subroutine generated_product_carrier',
                          generated_product)
            self.assertIn('subroutine generated_product_kernel',
                          generated_product)
            self.assertIn('PRODUCT_CARRIER_001_CONTRACT',
                          generated_product)
            self.assertNotIn('\n     $', generated_product)
            with open(os.path.join(
                    subprocess_dir, 'driver_mintFO.f90')) as stream:
                driver_source = stream.read()
            self.assertIn('call initialize_multiplicative_product()',
                          driver_source)
            self.assertTrue(os.path.isfile(os.path.join(
                subprocess_dir, 'nlo_decay_info_2.dat')))
            self.assertTrue(os.path.isfile(os.path.join(
                subprocess_dir, 'nlo_decay_info_3.dat')))
            generated_links = len([
                name for name in os.listdir(subprocess_dir)
                if name.startswith('b_sf_') and name.endswith('.f')])
            for contribution in [2, 3]:
                with open(os.path.join(
                        subprocess_dir,
                        'nlo_decay_info_%d.dat' % contribution)) as stream:
                    counts = next(
                        line for line in stream if line.startswith('COUNTS '))
                self.assertEqual(int(counts.split()[3]), generated_links)

    def test_nlo_decay_generation_restrictions(self):
        command = self.generate(
            'u u~ > t t~, '
            '(t > w+ b QED^2=2 QCD^2=0 [real=QCD])')
        with self.assertRaisesRegex(
                InvalidCmd, 'can only be exported with "output fNLO"'):
            command.check_output([])
        command.check_output(['fNLO'])
        with self.assertRaisesRegex(
                InvalidCmd, 'generate/add-process combinations'):
            command.exec_cmd(
                'add process d d~ > t t~ [real=QCD]',
                printcmd=False, precmd=True)

        explicit_lo = self.generate(
            'u u~ > t t~ [LOonly], '
            '(t > w+ b QED^2=2 QCD^2=0 [real=QCD])')
        self.assertTrue(explicit_lo._fks_multi_proc.nlo_decay_prototype)
        self.assertEqual(
            [leg.get('id') for leg in explicit_lo._fks_multi_proc.
             nlo_decay_production_amplitudes[0]['process']['legs']],
            [2, -2, 6, -6])

        amplitude_order = self.generate(
            'u u~ > t t~, t > w+ b QED=1 [real=QCD]')
        decay_process = amplitude_order._fks_multi_proc[
            'born_processes'][0].born_amp.get('process')
        self.assertEqual(
            decay_process.get('born_sq_orders'), {'QED': 2, 'QCD': 0})
        self.assertEqual(
            decay_process.get('squared_orders'), {'QED': 2, 'QCD': 2})

        with self.assertRaisesRegex(
                InvalidCmd, 'exactly one perturbatively corrected decay'):
            self.generate(
                'u u~ > t t~, '
                '(t > w+ b QED^2=2 QCD^2=0 [real=QCD]), '
                '(t~ > w- b~ QED^2=2 QCD^2=0 [real=QCD])')

if __name__ == '__main__':
    import unittest as unittest_main
    unittest_main.main()
