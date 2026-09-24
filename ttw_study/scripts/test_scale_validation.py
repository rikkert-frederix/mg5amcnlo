import copy
import contextlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

import numpy as np

from campaign import ROOT, STUDY, digest
from load_results import POINTS, assert_matched_physics
from run_scale_validation import (HOOK, assert_technical_match, cases,
                                  compare, install_hook, overlapping_points,
                                  physical_scale_coordinates, run_case)
import test_load_results as load_fixtures


class TestScaleValidation(unittest.TestCase):
    def test_direct_weight_coordinates_are_not_relabelled(self):
        self.assertEqual(len(cases()),6)
        for case in cases():
            pairs=overlapping_points(case['factors'])
            self.assertEqual(pairs[0],(0,1+POINTS.index(tuple(case['factors'])),case['factors']))
            expected=37 if case['name'] in ('asymmetric_t2_at_half','production_r2_f_half') else 82
            self.assertEqual(len(pairs),expected)
            for direct,reference,absolute in pairs[1:]:
                self.assertEqual(tuple(absolute),POINTS[reference-1])
                self.assertEqual(tuple(absolute),tuple(a*b for a,b in zip(POINTS[direct-1],case['factors'])))
        point=physical_scale_coordinates(cases()[0])
        self.assertEqual(point['decay_numerator_scales_GeV'],[345.,86.25])
        self.assertEqual(point['top_width_coefficient_reference_GeV'],172.5)

    def samples(self,case):
        reference=load_fixtures.TestMatchedIdentity().samples()['S']
        reference['report']['manifest'].update(variant='S',w_treatment=case['mode'],decay_bottom_mass=0.)
        reference['parameter_signature']='same verified parameters'
        reference['seeds']=np.arange(10).reshape(5,2)
        direct=copy.deepcopy(reference)
        direct['seeds']+=100
        manifest=direct['report']['manifest']
        manifest.update(technical_test=case,physical_scale_coordinates=physical_scale_coordinates(case),
                        custom_user_hook={'sha256':digest(HOOK)} if case['group']=='asymmetric' else None)
        manifest['settings'].update(iseed=51001,mur_over_ref=case['factors'][0],
                                    muf_over_ref=case['factors'][1],qes_over_ref=case['qes'])
        return direct,reference

    def test_only_declared_scale_changes_enter_an_identity(self):
        for case in cases()[:4]:
            direct,reference=self.samples(case)
            assert_technical_match(direct,reference,case)
            with self.assertRaises(ValueError):
                assert_matched_physics([direct,reference])
            for change in ('pdf','width','hook','seed','coordinates'):
                bad=copy.deepcopy(direct)
                if change=='pdf':
                    bad['report']['manifest']['settings']['lhaid']=14400
                elif change=='width':
                    bad['report']['manifest']['top_width_nlo']=1.4
                elif change=='hook':
                    bad['report']['manifest']['custom_user_hook']={'sha256':'wrong'}
                elif change=='seed':
                    bad['seeds']=reference['seeds']
                else:
                    bad['report']['manifest']['physical_scale_coordinates']['qes_over_ref']=3.
                with self.subTest(case=case['name'],change=change), self.assertRaises(ValueError):
                    assert_technical_match(bad,reference,case)

    def test_disabled_split_policy_matches_an_older_control(self):
        config=next(case for case in cases() if case['name']=='qes_half')
        direct,reference=self.samples(config)
        direct['report']['manifest']['settings'].update(
            fo_split_outlier_threshold=0.,fo_split_outlier_variance_fraction=.95,
            fo_split_outlier_min_splits=8)
        assert_technical_match(direct,reference,config)
        direct['report']['manifest']['settings']['fo_split_outlier_threshold']=10.
        with self.assertRaisesRegex(ValueError,'split-exclusion'):
            assert_technical_match(direct,reference,config)

    def test_hook_installation_is_private_and_never_replaces_a_used_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            process=Path(temporary)
            (process/'SubProcesses').mkdir()
            (process/'Events').mkdir()
            (process/'Events/.keep').touch()
            target=process/'SubProcesses/dummy_fct.f90'
            default=ROOT/'Template/fNLO/SubProcesses/dummy_fct.f90'
            shutil.copy2(default,target)
            old=digest(default)
            record=install_hook(process)
            self.assertEqual(digest(target),record['sha256'])
            self.assertEqual(digest(default),old)
            self.assertEqual(digest(process/'study_cards/asymmetric_user_hook/default_dummy_fct.f90'),old)
            with self.assertRaises(ValueError):
                install_hook(process)
        with tempfile.TemporaryDirectory() as temporary:
            process=Path(temporary)
            (process/'SubProcesses').mkdir()
            (process/'Events/existing_run').mkdir(parents=True)
            shutil.copy2(default,process/'SubProcesses/dummy_fct.f90')
            with self.assertRaises(ValueError):
                install_hook(process)

    def test_card_override_archive_preserves_the_unexecuted_base(self):
        import tests.unit_tests.fks.test_ttw_product_tools as fixtures
        from madgraph.various.banner import RunCardNLO
        source=STUDY/'processes/TTWplus_onshell_eemu_mb0p0_both_r9/study_cards/S_onshell_core-w-ht-half_separate_32001'
        if not source.is_dir():
            self.skipTest('Physical pilot cards are required for this workflow regression')
        with tempfile.TemporaryDirectory() as temporary:
            study=Path(temporary)
            (study/'inputs').mkdir()
            (study/'logs').mkdir()
            shutil.copy2(STUDY/'inputs/benchmark.json',study/'inputs/benchmark.json')
            for i,config in enumerate(cases()):
                process=study/('fake_'+config['name'])
                cards=process/'Cards'
                cards.mkdir(parents=True)
                (process/'Events').mkdir()
                sub=process/'SubProcesses'
                sub.mkdir()
                (process/'FixedOrderAnalysis').mkdir()
                for name in ('run_card.dat','decay_card.dat','param_card.dat','FO_analyse_card.dat'):
                    shutil.copy2(source/name,cards/name)
                shutil.copy2(source/'run_card.dat',cards/'run_card_default.dat')
                for name in ('decay_chain_parameters.f90','setscales.f90','decay_chain_scales.f90','dummy_fct.f90',
                             'phase_space_kinematics.f90','factorized_block_kinematics.f90',
                             'decay_chain_kinematics.f90','nlo_decay_kinematics.f90'):
                    shutil.copy2(ROOT/'Template/fNLO/SubProcesses'/name,sub/name)
                shutil.copy2(ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90',
                             process/'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90')
                fixtures.TestTTWProductTools.topology_fixture(process,config['mode'])
                (study/'inputs'/(process.name+'_generation.json')).write_text(json.dumps({'args':{}}))
                hook=install_hook(process) if config['group']=='asymmetric' else None
                def fake_launch(command,cwd,log,path,record):
                    name=command[command.index('-n')+1]
                    (cwd/'Events'/name).mkdir()
                    (cwd/'Events'/name/'MADatNLO.HwU').write_text('mock launch, not a physical histogram\n')
                    record.update(status='finished',returncode=0)
                with mock.patch('run_scale_validation.STUDY',study), \
                     mock.patch('run_scale_validation.launch',side_effect=fake_launch), \
                     mock.patch('run_scale_validation.audit'), mock.patch('run_scale_validation.harvest'), \
                     contextlib.redirect_stdout(io.StringIO()):
                    result=run_case(process,config,'S',81001+i,.05,{},hook)
                archive=process/'study_cards'/result['run']
                manifest=json.loads((archive/'manifest.json').read_text())
                base=Path(manifest['base_configuration'])
                baseline=json.loads(base.read_text())
                self.assertEqual(digest(base),manifest['base_configuration_sha256'])
                for name in ('param_card.dat','run_card.dat','decay_card.dat','FO_analyse_card.dat'):
                    self.assertEqual(digest(base.parent/name),baseline['hashes'][name])
                    self.assertEqual(digest(archive/name),manifest['hashes'][name])
                self.assertFalse((process/'Events'/baseline['run_name']).exists())
                run=RunCardNLO(str(archive/'run_card.dat'))
                for key,value in (('mur_over_ref',config['factors'][0]),
                                  ('muf_over_ref',config['factors'][1]),('qes_over_ref',config['qes'])):
                    self.assertEqual(run[key],value)
                    self.assertEqual(manifest['settings'][key],value)
                self.assertEqual(list(run['lhaid']),[331700])
                self.assertTrue(run['reweight_scale'][0])
                self.assertTrue(run['reweight_pdf'][0])
                self.assertEqual(run['fixed_ren_scale'],config['production_scale']=='fixed')
                if hook:
                    self.assertEqual(digest(archive/HOOK.name),digest(HOOK))
                    self.assertIn('-1 = decay_dynamical_scale_choice(6)',(archive/'decay_card.dat').read_text())

    def test_joint_comparison_uses_absolute_coordinates_and_parent_rates(self):
        config=cases()[0]
        direct,reference=self.samples(config)
        offsets=np.r_[0,9,np.arange(10,30)]
        base=np.r_[[100.,80.,70.,20.,10.,5.,3.,2.,0.],np.arange(1.,21.)]
        for row,scale,noise in ((direct,config['factors'],[.6,.9,1.,1.1,1.4]),
                                (reference,[1.]*4,[.8,.9,1.,1.1,1.2])):
            coordinates=[np.asarray(scale)]+[np.asarray(p)*scale for p in POINTS]
            weights=np.r_[[sum(p) for p in coordinates],np.full(101,sum(scale))]
            row.update(offsets=offsets,edges=np.asarray([[0.,1.]]*29),
                       titles=['R04_b25 W+ rates: stages']+['spectrum %d'%i for i in range(2,22)],
                       strata=np.zeros(5,dtype=int),
                       contributions=np.asarray(noise)[:,None,None]*base[None,:,None]*weights[None,None,:]/5.)
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary)
            paths=[directory/name for name in ('direct.npz','reference.npz')]
            for path in paths:
                path.write_text('synthetic inputs, loaded from test fixture')
            destination=directory/'comparison.json'
            with mock.patch('run_scale_validation.load_batches',side_effect=[direct,reference]):
                compare(*paths,config,destination)
            report=json.loads(destination.read_text())
            row=report['rates']['fiducial_2b']
            np.testing.assert_allclose(row['difference_pb'],0.,atol=1.e-13)
            self.assertGreater(row['difference_mc_error_pb'][0],0.)
            self.assertEqual(report['reference_points'][0],[1.,1.,2.,.5])
            with np.load(destination.with_suffix('.npz')) as data:
                for local in (3,4,6,11,12,20):
                    np.testing.assert_allclose(data['normalized_%d'%local][2],0.,atol=1.e-13)
                    parent=10. if local in (11,12) else 20.
                    self.assertAlmostEqual(data['normalized_%d'%local][0,0,0],(local-1.)/parent)

    def test_compiled_signed_hook_matches_the_existing_reweighting_runtime(self):
        import tests.unit_tests.fks.test_fnlo_decay_card as fixture
        fixture.TestFNLODecayCard.setUpClass()
        case=fixture.TestFNLODecayCard()
        try:
            template=ROOT/'Template/fNLO/SubProcesses'
            inputs=ROOT/'tests/input_files/fks_decay'
            build=Path(case.build.name)
            executable=build/'check_asymmetric'
            sources=[inputs/'decay_card_runtime_stubs.f90',
                     template/'decay_chain_parameters.f90',template/'factorized_phase_space.f90',HOOK,
                     *[template/(name+'.f90') for name in ('decay_chain_scales','setscales',
                        'fnlo_scale_variations','weight_lines','spin_density_weight_lines')],
                     inputs/'decay_card_runtime_driver.f90']
            result=subprocess.run(['gfortran','-O0','-g','-fcheck=all','-ffree-line-length-none',
                                   '-o',str(executable),*map(str,sources)],cwd=build,
                                  stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            self.assertEqual(result.returncode,0,result.stdout)
            standard=case.executable
            for combination in ('ADDITIVE','MULTIPLICATIVE'):
                case.executable=standard
                options=dict(decay_scale_grouping='SIGNED_PDG',nlo_decay_combination=combination,
                             decay_width_scale_modes={6:'AUTO'})
                reference=case.run_card(case.card(**options),'split')
                case.executable=str(executable)
                direct=case.run_card(case.card(**options,decay_dynamical_scale_choices={6:-1},
                                              decay_scale_variation_mode='NONE'))
                np.testing.assert_allclose(list(map(float,direct['SCALES'])),[200.,50.])
                for key in ('COUNTERTERM','PRODUCT_WIDTH','LOCAL_WIDTHS','DENSITY_MULTIPLIERS',
                            'COUPLING','ANTITOP_COUPLING'):
                    np.testing.assert_allclose(list(map(float,direct[key])),
                                               list(map(float,reference[key])),rtol=2.e-14,atol=1.e-14)
                grid=case.run_card(case.card(**options,decay_dynamical_scale_choices={6:-1}),'grid')
                self.assertEqual(int(grid['GRID_COUNT'][0]),81)
        finally:
            fixture.TestFNLODecayCard.doClassCleanups()


if __name__=='__main__':
    unittest.main()
