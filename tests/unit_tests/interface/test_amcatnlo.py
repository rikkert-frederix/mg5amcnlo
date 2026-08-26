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

    def test_fnlo_run_card_include_is_reduced(self):
        """Hidden matching and EW fields do not leak into run_card.inc."""

        run_card = banner_mod.RunCardNLO()
        run_mecmd.prepare_fixed_order_only_run_card(run_card)
        include = StringIO()
        run_card.write_include_file(None, include)
        include_text = include.getvalue().lower()

        for parameter in run_mecmd.FNLO_OMITTED_RUN_CARD_PARAMETERS:
            self.assertNotIn(parameter, include_text)

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
                'fks_event_kinematics.f90',
                'fks_qcd_splitting.f90',
                'fks_soft_kernels.f90',
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
        for removed_source in [
                'cluster.f90', 'cluster_bridge.f', 'recmom.f90',
                'ewsudakov_functions.f90',
                'ewsudakov_functions_dummy.f90', 'sa_ewsudakov.f90',
                'sub_f2py_ewsudakov.f90']:
            self.assertFalse(os.path.exists(pjoin(
                subprocess_dir, removed_source)))

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

        with open(pjoin(subprocess_dir, 'fks_singular.f90')) as stream:
            fks_singular = stream.read()
        self.assertNotIn('compute_MC_subt_term', fks_singular)
        self.assertNotIn('replace_MC_subt', fks_singular)
        self.assertNotIn('factor_n1body_NLOPS', fks_singular)
        self.assertNotIn('subroutine AP_reduced', fks_singular)
        self.assertNotIn('subroutine eikonal_Ireg', fks_singular)
        self.assertNotIn('subroutine set_cms_stuff', fks_singular)

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
        self.assertNotIn('use fks_singular_module', genps_born)
        self.assertNotIn('use fks_singular_module', genps_fks)

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
