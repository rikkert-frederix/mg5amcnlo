##############################################################################
#
# Copyright (c) 2010 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################
from __future__ import absolute_import
from cmd import Cmd
""" Basic test of the command interface """

import unittest
import madgraph

import madgraph.interface.master_interface as master_cmd
import madgraph.interface.madgraph_interface  as mg_cmd
import madgraph.interface.extended_cmd as ext_cmd
import madgraph.interface.amcatnlo_interface as mecmd
import madgraph.interface.amcatnlo_run_interface as run_mecmd
import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.export_fks as export_fks
import madgraph.various.misc as misc
import madgraph.various.banner as banner_mod
import os
import tempfile
from io import StringIO
import logging

root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
root_path = os.path.dirname(root_path)
# root_path is ./tests
pjoin = os.path.join

class MGerror(Exception): pass

class TestMadEventCmd(unittest.TestCase):
    """ check if the ValidCmd works correctly """

    def test_fnlo_output_format(self):
        """The fNLO token selects its template instead of becoming a path."""

        interface = object.__new__(mecmd.aMCatNLOInterface)
        interface._fks_multi_proc = object()
        interface._curr_model = {'name': 'sm'}

        with tempfile.TemporaryDirectory() as output_root:
            interface.writing_dir = output_root
            args = ['fNLO']
            interface.check_output(args)

            self.assertEqual(interface._export_format, 'fNLO')
            self.assertEqual(args, [])
            self.assertEqual(interface._export_dir,
                             os.path.realpath(pjoin(
                                 output_root, 'PROC_fNLO_sm_0')))

            custom_path = pjoin(output_root, 'custom')
            args = ['fNLO', custom_path]
            interface.check_output(args)
            self.assertEqual(interface._export_dir,
                             os.path.realpath(custom_path))
            self.assertEqual(args, [])

    def test_fnlo_rejects_non_qcd_corrections(self):
        """The reduced fNLO template accepts QCD corrections only."""

        class FakeFKSMultiProcess(object):

            def __init__(self, ewsudakov=False):
                self.ewsudakov = ewsudakov

            def get(self, name):
                if name == 'ewsudakov':
                    return self.ewsudakov
                return None

        interface = object.__new__(mecmd.aMCatNLOInterface)
        interface._curr_model = {'name': 'sm'}

        with tempfile.TemporaryDirectory() as output_root:
            interface.writing_dir = output_root
            interface._generate_info = 'p p > t t~ [QCD]'
            interface._fks_multi_proc = FakeFKSMultiProcess()
            interface.check_output(['fNLO'])

            interface._generate_info = 'p p > t t~ [QED]'
            interface._fks_multi_proc = FakeFKSMultiProcess()
            with self.assertRaisesRegex(
                    interface.InvalidCmd, 'QCD corrections only'):
                interface.check_output(['fNLO'])

            interface._generate_info = 'p p > t t~ [QCD]'
            interface._fks_multi_proc = FakeFKSMultiProcess(ewsudakov=True)
            with self.assertRaisesRegex(
                    interface.InvalidCmd, 'EW Sudakov corrections'):
                interface.check_output(['fNLO'])

        self.assertEqual(
            export_fks.get_nlo_correction_orders('p p > t t~ [QCD]'),
            ['QCD'])
        self.assertEqual(
            export_fks.get_nlo_correction_orders(
                'p p > t t~ [real=QCD] --no_warning=duplicate'),
            ['QCD'])
        self.assertEqual(
            export_fks.get_nlo_correction_orders('p p > t t~ [LOonly]'),
            ['QCD'])

    def test_fnlo_exporter_rejects_mixed_born_orders(self):
        """The exporter enforces the reduced fNLO physics invariant."""

        class FakeBorn(object):

            def __init__(self, squared_orders, amplitude_orders=None):
                self.squared_orders = squared_orders
                self.amplitude_orders = (amplitude_orders if
                                         amplitude_orders is not None else
                                         [((2, 0), (1,))])

            def get_split_orders_mapping(self):
                return self.squared_orders, self.amplitude_orders

        class FakeRealProcess(object):

            def __init__(self, squared_orders, amplitude_orders=None):
                self.matrix_element = FakeBorn(
                    squared_orders, amplitude_orders)

        class FakeMatrixElement(object):

            def __init__(self, squared_orders=((4, 0),),
                         real_squared_orders=((6, 0),),
                         amplitude_orders=None,
                         real_amplitude_orders=None,
                         virtual_squared_orders=None,
                         splitting_types=('QCD',), extra=False,
                         ewsudakov=False):
                self.born_me = FakeBorn(squared_orders, amplitude_orders)
                self.real_processes = [FakeRealProcess(
                    real_squared_orders, real_amplitude_orders)]
                self.virt_matrix_element = (FakeBorn(
                    virtual_squared_orders) if
                    virtual_squared_orders is not None else None)
                self.extra_cnt_me_list = [object()] if extra else []
                self.ewsudakov = ewsudakov
                self.splitting_types = splitting_types

            def get_fks_info_list(self):
                return [{'fks_info': {
                    'splitting_type': list(self.splitting_types)}}]

        exporter = object.__new__(export_fks.ProcessExporterFortranFKS)
        exporter.opt = {'fks_template': 'fNLO'}
        exporter.validate_fnlo_matrix_element(FakeMatrixElement())

        with self.assertRaisesRegex(Exception, 'one maximally QCD-like Born'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                squared_orders=((4, 0), (2, 2))))
        with self.assertRaisesRegex(Exception, 'one maximally QCD-like Born'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                amplitude_orders=[((2, 0), (1,)), ((0, 2), (2,))]))
        with self.assertRaisesRegex(Exception, 'one QCD-corrected real'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                real_squared_orders=((6, 0), (4, 2))))
        with self.assertRaisesRegex(Exception, 'one QCD-corrected real'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                real_amplitude_orders=[
                    ((3, 0), (1,)), ((1, 2), (2,))]))
        with self.assertRaisesRegex(Exception, 'one QCD-corrected virtual'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                virtual_squared_orders=((6, 0), (4, 2))))
        with self.assertRaisesRegex(Exception, 'QCD FKS splittings only'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                splitting_types=('QED',)))
        with self.assertRaisesRegex(Exception, 'extra counterterms'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                extra=True))
        with self.assertRaisesRegex(Exception, 'EW Sudakov corrections'):
            exporter.validate_fnlo_matrix_element(FakeMatrixElement(
                ewsudakov=True))

        born_template_path = pjoin(
            madgraph.MG5DIR, 'madgraph', 'iolibs', 'template_files',
            'bornmatrix_splitorders_fks.inc')
        with open(born_template_path) as stream:
            born_template = stream.read()
        with open(pjoin(
                madgraph.MG5DIR, 'madgraph', 'iolibs', 'template_files',
                'bornmatrix_qcd_fks.inc')) as stream:
            born_wrapper = stream.read()
        reduced_template = exporter.specialize_fnlo_born_template(
            born_wrapper, born_template).lower()
        self.assertNotIn('split_type_used', reduced_template)
        self.assertNotIn('has_ewsudakov', reduced_template)
        self.assertNotIn('subroutine sborn_onehel', reduced_template)
        self.assertNotIn('keep_order', reduced_template)
        self.assertNotIn('force_ijglu_zero', reduced_template)
        self.assertIn('subroutine sborn_splitorders', reduced_template)

        for wrapper_name, generic_name, suffix_marker in [
                ('born_hel_qcd_fks.inc',
                 'born_hel_splitorders_fks.inc',
                 'SUBROUTINE SBORN_HEL_SPLITORDERS'),
                ('realmatrix_qcd_fks.inc',
                 'realmatrix_splitorders_fks.inc',
                 'SUBROUTINE SMATRIX%(proc_prefix)s_SPLITORDERS'),
                ('b_sf_xxx_qcd_fks.inc',
                 'b_sf_xxx_splitorders_fks.inc',
                 'SUBROUTINE SB_SF_%(ilink)3.3d_SPLITORDERS')]:
            with open(pjoin(
                    madgraph.MG5DIR, 'madgraph', 'iolibs', 'template_files',
                    wrapper_name)) as stream:
                wrapper = stream.read()
            with open(pjoin(
                    madgraph.MG5DIR, 'madgraph', 'iolibs', 'template_files',
                    generic_name)) as stream:
                generic_template = stream.read()
            reduced_template = exporter.compose_fnlo_matrix_template(
                wrapper, generic_template, suffix_marker).lower()
            for obsolete_selector in [
                    'keep_order', 'firsttime', 'born_orders', 'nlo_orders',
                    'split_type', 'ewsudakov']:
                self.assertNotIn(obsolete_selector, reduced_template)

    def test_fnlo_run_card_include_is_reduced(self):
        """Hidden matching and EW fields do not leak into run_card.inc."""

        run_card = banner_mod.RunCardNLO()
        self.assertFalse(run_card['cut_decays'])
        run_mecmd.prepare_fixed_order_only_run_card(run_card)
        include = StringIO()
        run_card.write_include_file(None, include)
        include_text = include.getvalue().lower()

        for parameter in run_mecmd.FNLO_OMITTED_RUN_CARD_PARAMETERS:
            self.assertNotIn(parameter, include_text)
        self.assertIn('cut_decays = .false.', include_text)

        with tempfile.TemporaryDirectory() as output_dir:
            card_path = pjoin(output_dir, 'run_card.dat')
            template = pjoin(
                madgraph.MG5DIR, 'Template', 'fNLO', 'Cards',
                'run_card.dat')
            run_card.write(
                card_path, template=template, python_template=True)
            with open(card_path) as stream:
                rendered_card = stream.read().lower()
            self.assertIn('false = cut_decays', rendered_card)

        run_card.set('cut_decays', True)
        include = StringIO()
        run_card.write_include_file(None, include)
        self.assertIn('cut_decays = .true.', include.getvalue().lower())

    def test_fnlo_decay_cut_opt_out_is_wired(self):
        """Decay leaves are hidden only in the generation-cut event view."""

        subprocess_dir = pjoin(
            madgraph.MG5DIR, 'Template', 'fNLO', 'SubProcesses')
        with open(pjoin(subprocess_dir, 'decay_chain_kinematics.f90')) \
                as stream:
            kinematics = stream.read().lower()
        self.assertIn('call set_decay_cut_mask(context)', kinematics)
        self.assertIn(
            'visible = leaf_visible_leg(context, leaf)', kinematics)
        self.assertIn(
            'event_from_decay(visible) = .true.', kinematics)

        with open(pjoin(subprocess_dir, 'nlo_decay_kinematics.f90')) \
                as stream:
            nlo_kinematics = stream.read().lower()
        self.assertIn(
            'call set_nlo_decay_cut_mask(context)', nlo_kinematics)
        self.assertIn(
            'visible_leg = nlo_decay_leaf_visible(context, leaf)',
            nlo_kinematics)
        self.assertIn(
            'nlo_decay_local_target_kind(context, leg)', nlo_kinematics)
        self.assertIn(
            'event_from_decay(visible_leg) = .true.', nlo_kinematics)

        with open(pjoin(subprocess_dir, 'cuts.f90')) as stream:
            cuts = stream.read().lower()
        mask_call = cuts.index(
            'call apply_decay_cut_mask(pp,istatus,ipdg)')
        cut_call = cuts.index('passcuts = passcuts_user(pp,istatus,ipdg)')
        self.assertLess(mask_call, cut_call)
        self.assertIn('if (cut_decays) return', cuts)
        self.assertIn('if (event_from_decay(i)) then', cuts)
        self.assertIn('istatus(i)=decay_cut_status', cuts)
        self.assertIn('ipdg(i)=decay_cut_pdg', cuts)
        self.assertIn('passcuts_pdgs(p_reco,istatus)', cuts)
        self.assertIn('if (istatus(i).ne.1) cycle', cuts)

        with open(pjoin(subprocess_dir, 'setcuts_bridge.f')) as stream:
            setcuts_bridge = stream.read().lower()
        self.assertIn('from_decay=.false.', setcuts_bridge)
        self.assertIn('if (has_nlo_decay()) then', setcuts_bridge)
        self.assertIn('call set_nlo_decay_tau_min_impl', setcuts_bridge)
        self.assertIn('if (has_decay_chains()) then', setcuts_bridge)
        self.assertIn('call set_decay_tau_min_impl', setcuts_bridge)

        with open(pjoin(subprocess_dir, 'setcuts.f90')) as stream:
            setcuts = stream.read().lower()
        self.assertIn(
            'do leg = nincoming + 1, context_core_count(context)', setcuts)
        self.assertIn('case (decay_node_target)', setcuts)
        self.assertIn('production_mass = production_mass + leg_mass', setcuts)
        self.assertIn(
            'do leg = nincoming + 1, nlo_decay_production_count()',
            setcuts)
        self.assertIn('case (nlo_decay_node_target)', setcuts)

        with open(pjoin(subprocess_dir, 'genps_born.f90')) as stream:
            born_phase_space = stream.read().lower()
        decay_start = born_phase_space.index(
            'subroutine generate_decay_born_phase_space')
        decay_end = born_phase_space.index(
            'end subroutine generate_decay_born_phase_space', decay_start)
        self.assertIn(
            'call set_tau_min()', born_phase_space[decay_start:decay_end])
        nlo_decay_start = born_phase_space.index(
            'subroutine generate_nlo_decay_born_phase_space')
        nlo_decay_end = born_phase_space.index(
            'end subroutine generate_nlo_decay_born_phase_space',
            nlo_decay_start)
        self.assertIn(
            'call set_tau_min()',
            born_phase_space[nlo_decay_start:nlo_decay_end])

        decay_mask_path = pjoin(
            madgraph.MG5DIR, 'Template', 'NLO', 'SubProcesses',
            'decay_cut_mask.inc')
        with open(decay_mask_path) as stream:
            decay_mask = stream.read().lower()
        self.assertIn('logical from_decay(nexternal)', decay_mask)
        self.assertIn('common /to_decay_cut_mask/from_decay', decay_mask)

    def test_fnlo_rejects_non_msbar_pdf_scheme(self):
        """The reduced fNLO template supports MSbar factorization only."""

        run_card = banner_mod.RunCardNLO()
        run_card.set('pdfscheme', 1)
        with self.assertRaisesRegex(
                banner_mod.InvalidRunCard, 'MSbar PDF factorization scheme'):
            run_mecmd.prepare_fixed_order_only_run_card(run_card)

    def test_fnlo_test_me_failures_propagate(self):
        """A crashed test executable cannot be reported as passing."""

        with tempfile.TemporaryDirectory() as me_dir:
            subprocess_dir = pjoin(me_dir, 'SubProcesses', 'P0_test')
            os.makedirs(subprocess_dir)
            with open(pjoin(me_dir, 'test_ME_input.txt'), 'w') as stream:
                stream.write('-2 -2\n1 1\n0\n0\n')

            original_compile = run_mecmd.misc.compile
            original_call = run_mecmd.misc.call
            run_mecmd.misc.compile = lambda *args, **opts: None
            run_mecmd.misc.call = lambda *args, **opts: 7
            try:
                with self.assertRaisesRegex(
                        run_mecmd.aMCatNLOError,
                        'test_ME failed .* exit status 7'):
                    run_mecmd.compile_dir(
                        me_dir, 'P0_test', 'fixed_order',
                        {'reweightonly': True}, ['test_ME'], 'madevent', 0)
            finally:
                run_mecmd.misc.compile = original_compile
                run_mecmd.misc.call = original_call

            log_path = pjoin(subprocess_dir, 'test_ME.log')
            with open(log_path, 'w') as stream:
                stream.write(
                    'Fatal error in NLO-decay phase space: '
                    'insufficient energy\n')
            with self.assertRaises(run_mecmd.aMCatNLOError):
                run_mecmd.aMCatNLOCmd.parse_test_mx_log(None, log_path)

    def test_fnlo_exporter_and_template(self):
        """The fNLO factory route uses the fixed-order-only template."""

        class EmptyFKSMultiProcess(object):

            def get(self, name):
                return []

            def get_virt_amplitudes(self):
                return []

        interface = master_cmd.MasterCmd()
        interface._curr_amps = []
        interface._fks_multi_proc = EmptyFKSMultiProcess()
        interface._curr_model = {'running_elements': []}
        interface._export_dir = pjoin(root_path, 'unused_fnlo_output')
        interface.options['loop_optimized_output'] = False

        exporter = export_v4.ExportV4Factory(
            interface, False, output_type='fnlo', group_subprocesses=False)

        template_dir = pjoin(madgraph.MG5DIR, 'Template', 'fNLO')
        self.assertEqual(exporter.get_fks_template_dir(), template_dir)
        self.assertFalse(os.path.exists(pjoin(template_dir, 'Utilities')))
        self.assertFalse(os.path.exists(pjoin(template_dir, 'MCatNLO')))
        analysis_dir = pjoin(template_dir, 'FixedOrderAnalysis')
        analysis_files = os.listdir(analysis_dir)
        self.assertFalse(any(name.startswith('analysis_root_')
                             for name in analysis_files))
        self.assertFalse(any(name.startswith('analysis_td_')
                             for name in analysis_files))
        for removed_analysis_support in [
                'dbook.f', 'dbook.inc', 'rbook_be8.cc', 'rbook_fe8.f']:
            self.assertNotIn(removed_analysis_support, analysis_files)
        subprocess_dir = pjoin(template_dir, 'SubProcesses')
        runtime_common = pjoin(template_dir, 'Source',
                               'fnlo_runtime_common.f90')
        process_common = pjoin(subprocess_dir, 'fnlo_process_common.f')
        self.assertTrue(os.path.exists(runtime_common))
        self.assertTrue(os.path.exists(process_common))
        for removed_common_include in [
                'pineappl_common.inc', 'q_es.inc',
                'reweight_pineappl.inc', 'timing_variables.inc']:
            self.assertFalse(os.path.lexists(pjoin(
                subprocess_dir, removed_common_include)))

        common_owners = {
            os.path.realpath(runtime_common),
            os.path.realpath(process_common),
        }
        common_sources = set()
        for source_root in [pjoin(template_dir, 'Source'), subprocess_dir,
                            analysis_dir]:
            for directory, _, filenames in os.walk(source_root):
                for filename in filenames:
                    if not filename.lower().endswith(('.f', '.f90')):
                        continue
                    source_path = pjoin(directory, filename)
                    with open(source_path) as stream:
                        if any(line.lstrip().lower().startswith('common')
                               for line in stream):
                            common_sources.add(os.path.realpath(source_path))
        self.assertEqual(common_sources, common_owners)

        compatibility_includes = {
            os.path.realpath(pjoin(template_dir, 'Source', 'run.inc')),
            os.path.realpath(pjoin(template_dir, 'Source', 'alfas.inc')),
            os.path.realpath(pjoin(template_dir, 'Source', 'PDF',
                                   'pdf.inc')),
        }
        common_includes = set()
        for source_root in [pjoin(template_dir, 'Source'), subprocess_dir,
                            analysis_dir]:
            for directory, _, filenames in os.walk(source_root):
                for filename in filenames:
                    if not filename.lower().endswith('.inc'):
                        continue
                    source_path = pjoin(directory, filename)
                    with open(source_path) as stream:
                        if any(line.lstrip().lower().startswith('common')
                               for line in stream):
                            common_includes.add(os.path.realpath(source_path))
        self.assertEqual(common_includes, compatibility_includes)

        self.assertTrue(os.path.exists(pjoin(
            subprocess_dir, 'driver_mintFO.f90')))
        self.assertTrue(os.path.exists(pjoin(
            subprocess_dir, 'driver_mintFO_bridge.f')))
        self.assertTrue(os.path.exists(pjoin(
            subprocess_dir, 'genps_born.f90')))
        for physics_module in [
                'phase_space_kinematics.f90',
                'fks_qcd_splitting.f90',
                'fks_soft_kernels.f90',
                'fks_model_state.f90',
                'fks_singular.f90',
                'fks_contributions.f90',
                'fks_weights.f90',
                'fks_channel_map.f90',
                'fks_diagnostics.f90',
                'fks_random.f90']:
            self.assertTrue(os.path.exists(pjoin(
                subprocess_dir, physics_module)))
        self.assertFalse(os.path.exists(pjoin(
            subprocess_dir, 'kinematic_runtime_state.f90')))
        self.assertFalse(os.path.exists(pjoin(
            subprocess_dir, 'driver_mintMC.f')))
        self.assertFalse(os.path.exists(pjoin(
            subprocess_dir, 'montecarlocounter.f')))
        self.assertFalse(os.path.exists(pjoin(
            subprocess_dir, 'montecarlocounter_alt.f')))
        self.assertFalse(os.path.exists(pjoin(
            subprocess_dir, 'fks_powers.inc')))
        for removed_source in [
                'cluster.f90', 'cluster_bridge.f', 'recmom.f90',
                'cuts_bridge.f', 'genps_fks_bridge.f',
                'iproc_map.f90', 'iproc_map_bridge.f',
                'open_output_files_bridge.f',
                'pineappl_interface.cc',
                'pineappl_interface_dummy.f90',
                'pineappl_interface_dummy_bridge.f',
                'ewsudakov_functions.f90',
                'ewsudakov_functions_dummy.f90', 'sa_ewsudakov.f90',
                'sub_f2py_ewsudakov.f90']:
            self.assertFalse(os.path.exists(pjoin(
                subprocess_dir, removed_source)))

        with open(pjoin(template_dir, 'Cards', 'FKS_params.dat')) as stream:
            fks_parameters = stream.read().lower()
        for removed_order_filter in [
                '#qcd^2==', '#selectedcouplingorders',
                '#vetoedcontributiontypes']:
            self.assertNotIn(removed_order_filter, fks_parameters)

        for source_name in [
                'FKSParams.f90', 'fks_weights.f90',
                'process_dimensions.f90', 'process_dimensions_bridge.f']:
            with open(pjoin(subprocess_dir, source_name)) as stream:
                source = stream.read().lower()
            for removed_order_state in [
                    'qcd_squared_selected', 'selectedcouplingorders',
                    'vetoedcontributiontypes', 'nlo_orders']:
                self.assertNotIn(removed_order_state, source)

        self.assertFalse(os.path.lexists(pjoin(template_dir, 'Source',
                                              'cuts.inc')))
        self.assertFalse(os.path.lexists(pjoin(subprocess_dir, 'cuts.inc')))
        for template_source in ['run.inc', 'run_state.f90']:
            with open(pjoin(template_dir, 'Source', template_source)) as stream:
                source = stream.read().lower()
            for matching_name in [
                    'ickkw', 'xqcut', 'xmtc', 'ktscheme', 'chcluster',
                    'pdfwgt', 'to_cluster', 'to_specxpt']:
                self.assertNotIn(matching_name, source)
            self.assertNotIn('xbk', source)

        with open(pjoin(
                madgraph.MG5DIR, 'madgraph', 'iolibs', 'template_files',
                'parton_lum_n_fnlo.inc')) as stream:
            luminosity_template = stream.read().lower()
        self.assertIn('dlum_%(n_me)d(lum,bjorken_x)',
                      luminosity_template)
        self.assertNotIn('xbk', luminosity_template)

        with open(pjoin(subprocess_dir, 'fks_singular.f90')) as stream:
            fks_singular = stream.read()
        self.assertNotIn('compute_MC_subt_term', fks_singular)
        self.assertNotIn('replace_MC_subt', fks_singular)
        self.assertNotIn('factor_n1body_NLOPS', fks_singular)
        self.assertNotIn('subroutine AP_reduced', fks_singular)
        self.assertNotIn('subroutine eikonal_Ireg', fks_singular)
        self.assertNotIn('subroutine set_cms_stuff', fks_singular)
        self.assertNotIn('split_type', fks_singular.lower())
        self.assertNotIn('extra_cnt', fks_singular.lower())
        self.assertFalse(os.path.exists(pjoin(
            subprocess_dir, 'fks_event_kinematics.f90')))

        with open(pjoin(subprocess_dir, 'genps_born.f90')) as stream:
            genps_born = stream.read().lower()
        with open(pjoin(subprocess_dir, 'genps_fks.f90')) as stream:
            genps_fks = stream.read().lower()
        self.assertIn('subroutine generate_born_phase_space', genps_born)
        self.assertIn('subroutine generate_momenta_born', genps_born)
        self.assertNotIn('subroutine generate_fks_kinematics', genps_born)
        self.assertNotIn('subroutine generate_momenta_massless_final',
                         genps_born)
        self.assertIn('subroutine generate_fks_kinematics', genps_fks)
        self.assertIn('subroutine generate_momenta_massless_final', genps_fks)
        self.assertNotIn('subroutine generate_momenta_born', genps_fks)
        self.assertNotIn('subroutine generate_tau', genps_fks)
        self.assertNotIn('fks_as_is', genps_fks)
        self.assertNotIn('-2:2', genps_fks)
        self.assertIn('event_generation_order', genps_fks)
        self.assertIn('counterevent_loop: do event_position', genps_fks)
        self.assertIn('if (skip_counterevents) exit counterevent_loop',
                      genps_fks)
        self.assertNotIn('skipped_counterevents', genps_fks)
        self.assertNotIn('event_ycm', genps_fks)
        self.assertNotIn('fill_fks_commons', genps_fks)
        self.assertIn('store_fks_event', genps_fks)
        self.assertNotIn('goto 111', genps_fks)
        self.assertNotIn('111 continue', genps_fks)
        self.assertNotIn('use fks_singular_module', genps_born)
        self.assertNotIn('use fks_singular_module', genps_fks)

        with open(pjoin(subprocess_dir, 'fnlo_process_common.f')) as stream:
            process_common = stream.read().lower()
        self.assertIn('real_event=3', process_common)
        self.assertIn('event_momenta(0:3,nexternal,', process_common)
        self.assertIn('event_jacobian(soft_counterevent:real_event)',
                      process_common)
        self.assertIn(
            'ybst_til_tolab(soft_counterevent:real_event)',
            process_common)
        self.assertIn(
            'ybst_til_tocm(soft_counterevent:real_event)',
            process_common)
        self.assertNotIn('event_ycm', process_common)
        self.assertNotIn('parton_cms_stuff', process_common)
        for split_event_state in [
                'p_ev', 'p1_cnt', 'wgt_cnt', 'pswgt_cnt', 'jac_cnt',
                'xi_i_fks_ev', 'xi_i_fks_cnt', 'p_i_fks_ev',
                'p_i_fks_cnt']:
            self.assertNotIn(split_event_state, process_common)
        for removed_mixed_order_state in [
                'split_type', 'iextra_cnt', 'isplitorder_born',
                'isplitorder_cnt', 'amp_split_size_born']:
            self.assertNotIn(removed_mixed_order_state, process_common)

        for counterevent_source in [
                'driver_mintFO.f90', 'fks_Sij.f90',
                'fks_singular.f90', 'test_soft_col_limits.f90']:
            with open(pjoin(subprocess_dir, counterevent_source)) as stream:
                self.assertNotIn('-2:2', stream.read())

        with open(pjoin(subprocess_dir,
                        'test_soft_col_limits.f90')) as stream:
            limit_test = stream.read()
        self.assertNotIn('xmcsubt', limit_test)
        self.assertNotIn('amp_split_mc', limit_test)
        limit_test_lower = limit_test.lower()
        self.assertIn('initial_ebeam = ebeam', limit_test_lower)
        self.assertIn('max(initial_ebeam(1)', limit_test_lower)
        self.assertIn('nlo_decay_minimum_production_mass()',
                      limit_test_lower)
        self.assertNotIn('soft tests skipped', limit_test_lower)
        self.assertIn('nlo_decay_fks_sister_mass', limit_test_lower)

        with open(pjoin(subprocess_dir,
                        'nlo_decay_metadata.f90')) as stream:
            nlo_decay_metadata = stream.read().lower()
        self.assertIn("case ('fks_partner')", nlo_decay_metadata)
        self.assertIn("case ('color_link')", nlo_decay_metadata)
        self.assertIn('nlo_decay_partner_local', nlo_decay_metadata)
        self.assertIn('nlo_decay_map_color_link', nlo_decay_metadata)

        with open(pjoin(subprocess_dir, 'fks_singular.f90')) as stream:
            fks_singular = stream.read().lower()
        self.assertIn('sbornsoft_nlo_decay', fks_singular)
        self.assertIn('evaluate_nlo_decay_fks_sij', fks_singular)
        self.assertIn('nlo_decay_map_color_link', fks_singular)
        self.assertIn('link_multiplier*link_weight*eik', fks_singular)

        with open(pjoin(subprocess_dir,
                        'nlo_decay_kinematics.f90')) as stream:
            nlo_decay_kinematics = stream.read().lower()
        self.assertIn('local_event_cache', nlo_decay_kinematics)
        self.assertIn('get_nlo_decay_event_momenta',
                      nlo_decay_kinematics)

        
    def test_v31_syntax_crash(self):
        """Check that process with ambiguous syntax correctly crashes if the flag is not set correctly
        """
        #cmd = mecmd.aMCatNLOInterface()
        #category = set()
        #valid_command = [c for c in dir(cmd) if c.startswith('do_')]

        interface = master_cmd.MasterCmd()
        interface.no_notification()
        interface.do_import('model loop_qcd_qed_sm')
        
        #run_cmd('import model %s' % model)

        stream = StringIO()
        handler_stream = logging.StreamHandler(stream)
        log = logging.getLogger('cmdprint')
        log.setLevel(logging.CRITICAL)

        for handler in log.handlers: 
            log.removeHandler(handler)
        log.addHandler(handler_stream)


        def check_message(line):
            """return False if not warning is raised, return True is a warning is raised"""
            
            stream.seek(0)
            stream.truncate(0)
            myprocdef = interface.extract_process(line)
            try: 
                interface.proc_validity(myprocdef,'aMCatNLO_all')
            except Exception as error:
                if '1804.10017' in str(error):
                    raise MGerror

            text = stream.getvalue()
            if '1804.10017' in text:
                return True
            else:
                return False

        # force the option to not by bypassed
        interface.options['acknowledged_v3.1_syntax'] = False

        # check case where the code crash
        self.assertRaises(MGerror, check_message, "p p > t t~ QED=1 [QED]")
        self.assertRaises(MGerror,check_message, "p p > t t~ QCD=1 [QED QCD]")

        # check case where the code write a warning (critical level)
        self.assertRaises(MGerror, check_message, "p p > t t~ QED=1 [QCD]")
        self.assertRaises(MGerror, check_message, "p p > t t~ QCD=1 QED=0 [QCD]")
        self.assertRaises(MGerror, check_message, "p p > t t~ QED=98 [QCD]")
        self.assertRaises(MGerror, check_message, "p p > t t~ QED=99 QCD=1 [QCD]")
        # check case where the code does not complain
        self.assertFalse(check_message("p p > t t~ QED=0 [QCD]"))
        self.assertFalse(check_message("p p > t t~ / z QCD=0 [QCD]"))
        self.assertFalse(check_message("p p > t t~ QCD=1 QED^2==2 [QCD]"))
        self.assertFalse(check_message("p p > w+ w- QED=99 [QCD]"))
        self.assertFalse(check_message("p p > w+ w- j j $ h QED<=99 [QCD]"))

        # force the option to not  be by bypassed
        interface.options['acknowledged_v3.1_syntax'] = True 
        # and check that no crash/warning is raised anymore
        self.assertFalse(check_message( "p p > t t~ QED=1 [QED]"))
        self.assertFalse(check_message( "p p > t t~ QCD=1 [QED QCD]"))
        self.assertFalse(check_message("p p > t t~ QED=1 [QCD]"))
        self.assertFalse(check_message("p p > t t~ QCD=1 QED=0 [QCD]"))
        self.assertFalse(check_message("p p > t t~ QED=98 [QCD]"))
        self.assertFalse(check_message("p p > t t~ QED=99 QCD=1 [QCD]"))
