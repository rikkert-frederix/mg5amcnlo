import contextlib
import copy
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from campaign import ROOT, STUDY, digest
from run_small_mass_checks import cases, completed_exports, configure, recovery_plan, validate_cases
from small_mass_inputs import validate


class TestSmallMassChecks(unittest.TestCase):
    def inputs(self):
        path=STUDY/'inputs/small_mass_inputs_v1.json'
        return path,json.loads(path.read_text())

    def test_width_matching_rejects_missing_coordinates_and_wrong_coupling(self):
        _,inputs=self.inputs()
        benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
        validate(inputs,benchmark)
        for change in ('missing','duplicate','width','scale','W_width'):
            altered=copy.deepcopy(inputs)
            if change=='missing': altered['top_widths'].pop()
            elif change=='duplicate': altered['top_widths'][-1]=copy.deepcopy(altered['top_widths'][0])
            elif change=='width': altered['top_widths'][0]['gamma_nlo_pdf']*=1.01
            elif change=='scale': altered['top_widths'][0]['muR_GeV']*=2.
            else: altered['fixed_W_width_GeV']*=.5
            with self.assertRaises(ValueError):
                validate(altered,benchmark)

    def test_complete_mass_and_mode_inventory_with_fresh_reference_first(self):
        rows=cases()
        self.assertEqual(len(rows),12)
        self.assertEqual(len({r['seed'] for r in rows}),12)
        for mode in ('onshell','all-bw'):
            selected=[r for r in rows if r['w_treatment']==mode]
            self.assertEqual([r['decay_bottom_mass'] for r in selected],[0.,0.,1.,1.,.1,.1])
            self.assertEqual([r['variant'] for r in selected],['S','Pi']*3)

    def stopped_queue(self):
        selected=cases()
        return dict(status='stopped',cases=selected,current_case=selected[4],
                    jobs=[dict(case=case,process='kept',audit='frozen') for case in selected[:4]])

    def test_recovery_preserves_completed_jobs_and_only_reseeds_failed_case(self):
        previous=self.stopped_queue()
        original=copy.deepcopy(previous)
        plan=recovery_plan(previous,85013)
        self.assertEqual(plan['carried_jobs'],previous['jobs'])
        self.assertEqual(plan['cases'][:4],previous['cases'][:4])
        self.assertEqual(plan['cases'][5:],previous['cases'][5:])
        self.assertEqual(plan['cases'][4],dict(previous['current_case'],seed=85013))
        self.assertEqual(plan['failed_index'],4)
        self.assertEqual(len(plan['cases'][len(plan['carried_jobs']):]),8)
        plan['carried_jobs'][0]['audit']='altered copy'
        self.assertEqual(previous,original)
        validate_cases(plan['cases'])

    def test_recovery_rejects_unverified_order_status_and_reused_seed(self):
        previous=self.stopped_queue()
        for change in ('status','missing','order','wrong_case'):
            broken=copy.deepcopy(previous)
            if change=='status': broken['status']='running'
            elif change=='missing': broken['jobs'].pop()
            elif change=='order': broken['jobs'].reverse()
            else: broken['current_case']=dict(broken['current_case'],seed=99999)
            with self.subTest(change=change),self.assertRaises(ValueError):
                recovery_plan(broken,85013)
        for seed in (85001,85005,85012,0,900000001,True,85013.,None,[]):
            with self.subTest(seed=seed),self.assertRaises(ValueError):
                recovery_plan(previous,seed)

    def test_physical_inventory_cannot_change_when_replacing_a_seed(self):
        for change in ('mass','variant','order','missing','duplicate_seed','unknown_field'):
            broken=cases()
            if change=='mass': broken[4]['decay_bottom_mass']=.2
            elif change=='variant': broken[4]['variant']='Pi'
            elif change=='order': broken.reverse()
            elif change=='missing': broken.pop()
            elif change=='duplicate_seed': broken[4]['seed']=broken[0]['seed']
            else: broken[4]['extra']=True
            with self.subTest(change=change),self.assertRaises(ValueError):
                validate_cases(broken)

    def test_carried_exports_follow_actual_jobs_across_recoveries(self):
        previous=dict(carried_exports={'old':dict(process='old',production_identity='old_proof')},
            exports={'same_mass':dict(process='new',production_identity='new_proof'),
                     'failed':dict(process='unused')})
        jobs=[dict(process='old'),dict(process='old'),dict(process='new')]
        retained=completed_exports(previous,jobs)
        self.assertEqual(set(retained),{'old','new'})
        self.assertEqual([retained[p]['production_identity'] for p in ('old','new')],
                         ['old_proof','new_proof'])
        retained['old']['process']='changed copy'
        self.assertEqual(previous['carried_exports']['old']['process'],'old')
        with self.assertRaises(ValueError):
            completed_exports(previous,[dict(process='missing')])

    def test_actual_cards_use_matching_small_mass_widths_and_5fs_production(self):
        from models.check_param_card import ParamCard
        from madgraph.various.banner import RunCardNLO
        from tests.unit_tests.fks.test_ttw_product_tools import TestTTWProductTools
        path,inputs=self.inputs()
        with tempfile.TemporaryDirectory(prefix='small_mass_cards_') as temporary:
            for case in cases():
                process=Path(temporary)/str(case['seed'])
                mode=case['w_treatment']
                massive=case['decay_bottom_mass']>0
                source=STUDY/'processes'/('TTWplus_%s_eemu_mb%s_both_%s'%(
                    mode,'4p8' if massive else '0p0','width_mass_v4' if massive else 'alignment_controls_v1'))
                cards=process/'Cards'
                cards.mkdir(parents=True)
                (process/'Events').mkdir()
                sub=process/'SubProcesses'
                sub.mkdir()
                (process/'FixedOrderAnalysis').mkdir()
                for name in ('param_card.dat','run_card.dat','run_card_default.dat','decay_card.dat','FO_analyse_card.dat'):
                    shutil.copy2(source/'Cards'/name,cards/name)
                if massive:
                    shutil.copy2(source/'Cards/decay_mass_scheme.json',cards/'decay_mass_scheme.json')
                    model=process/'Source/MODEL'
                    model.mkdir(parents=True)
                    shutil.copy2(source/'Source/MODEL/get_mass_width_fcts.f',model/'get_mass_width_fcts.f')
                    param=ParamCard(str(cards/'param_card.dat'))
                    param['decaymass'].get((5,)).value=case['decay_bottom_mass']
                    param.write(str(cards/'param_card.dat'),precision=16)
                for name in ('decay_chain_parameters.f90','setscales.f90','decay_chain_scales.f90',
                    'factorized_block_kinematics.f90','phase_space_kinematics.f90',
                    'decay_chain_kinematics.f90','nlo_decay_kinematics.f90'):
                    shutil.copy2(ROOT/'Template/fNLO/SubProcesses'/name,sub/name)
                shutil.copy2(ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90',
                    process/'FixedOrderAnalysis/analysis_HwU_pp_ttxw_product.f90')
                TestTTWProductTools.topology_fixture(process,mode)
                with contextlib.redirect_stdout(io.StringIO()):
                    name,archive=configure(process,case,inputs,path,.01,
                        exclude_split_outliers=case['seed']==85001)
                manifest=json.loads((archive/'manifest.json').read_text())
                self.assertEqual(manifest['small_mass_inputs_sha256'],digest(path))
                self.assertEqual(manifest['small_mass_test'],case)
                param=ParamCard(str(archive/'param_card.dat'))
                actual=param['decaymass'].get((5,)).value if massive else 0.
                self.assertEqual(actual,case['decay_bottom_mass'])
                self.assertEqual(param['mass'].get((5,)).value,0.)
                width=next(r for r in inputs['top_widths'] if r['mb_GeV']==actual and
                           r['w_treatment']==('onshell' if mode=='onshell' else 'bw') and r['scale_factor']==1.)
                self.assertEqual(param['decay'].get((6,)).value,width['gamma_nlo_pdf'])
                self.assertEqual(manifest['top_width_lo'],width['gamma_lo'])
                self.assertEqual(manifest['top_width_nlo'],width['gamma_nlo_pdf'])
                run_card=RunCardNLO(str(archive/'run_card.dat'))
                self.assertFalse(run_card['fixed_ren_scale'])
                self.assertEqual(run_card['lhaid'],[331700])
                self.assertEqual(run_card['npoints_fo_grid'],1000)
                self.assertEqual(run_card['fo_split_outlier_threshold'],10. if case['seed']==85001 else 0.)
                if case['seed']==85001:
                    self.assertEqual(manifest['settings']['fo_split_outlier_min_splits'],8)
                self.assertFalse(list((process/'Events').iterdir()))


if __name__=='__main__':
    unittest.main()
