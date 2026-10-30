import argparse
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

import run_precision_campaign as precision


class TestPrecisionCampaign(unittest.TestCase):
    def test_complete_rounds_and_work_allocation(self):
        rows=precision.cases(10,901001)
        self.assertEqual(len(rows),320)
        self.assertEqual(len({r['seed'] for r in rows}),320)
        for round_index in range(10):
            selected=rows[32*round_index:32*(round_index+1)]
            self.assertEqual(len({(r['charge'],tuple(r['flavours']),r['variant']) for r in selected}),32)
            self.assertTrue(all(r['decay_bottom_mass']==0 and r['w_treatment']=='onshell' for r in selected))
        factor,accuracy=precision.allocation(16,10)
        self.assertEqual(factor,8)
        self.assertAlmostEqual(accuracy,.015/np.sqrt(8))
        self.assertAlmostEqual((.015/accuracy)**2*10/5,16)

    def test_invalid_allocations(self):
        for count,seed in ((4,1),(10,0),(10,900000000)):
            with self.assertRaises(ValueError):precision.cases(count,seed)
        for factor in (0,np.nan,np.inf):
            with self.assertRaises(ValueError):precision.allocation(factor,10)

    def test_execution_resumes_completed_integration_and_stops_after_one_comparison(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'state.json';cases=precision.cases(5,902001)
            path.write_text(json.dumps(dict(status='stopped',current_case=cases[0],jobs=[],
                cases=cases,tag='fixture',settings=dict(accuracy=.005))))
            def process(case):
                return root/'processes'/('TTW%s_onshell_%s_mb0p0_both_fixture'%(case['charge'],''.join(case['flavours'])))
            def finish(case):
                p=process(case);name='%s_onshell_core-w-ht-half_separate_%d'%(case['variant'],case['seed'])
                cards=p/'study_cards'/name;cards.mkdir(parents=True)
                (cards/'execution.json').write_text('{"status":"finished"}')
                event=p/'Events'/name;event.mkdir(parents=True)
                (event/'MADatNLO.HwU').write_text('fixture')
            finish(cases[0])
            def run(args):
                finish(next(c for c in cases if c['seed']==args.seed))
            def generate(args):
                process(dict(charge=args.charge,flavours=args.flavours)).mkdir(parents=True)
            with mock.patch.object(precision,'STUDY',root),                 mock.patch.object(precision,'unchanged',side_effect=AssertionError('No source scans')),                 mock.patch.object(precision,'generate',side_effect=generate),                 mock.patch.object(precision,'run',side_effect=run) as launch:
                precision.execute(path,1)
                state=json.loads(path.read_text());self.assertEqual(state['status'],'checkpoint')
                self.assertEqual(len(state['jobs']),32);self.assertIsNone(state['current_case'])
                self.assertEqual(launch.call_count,31)
                with self.assertRaisesRegex(ValueError,'next uncompleted'):precision.execute(path,1)
                precision.execute(path,2)
                state=json.loads(path.read_text());self.assertEqual(len(state['jobs']),64)
                self.assertEqual(launch.call_count,63)
                self.assertTrue(all(c.args[0].accuracy==.005 for c in launch.call_args_list))
                self.assertEqual(len({j['case']['seed'] for j in state['jobs']}),64)

    def test_incomplete_existing_integration_is_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'state.json';cases=precision.cases(5,902001)
            path.write_text(json.dumps(dict(status='stopped',current_case=cases[0],jobs=[],
                cases=cases,tag='fixture',settings=dict(accuracy=.005))))
            cards=root/'processes/TTWplus_onshell_eee_mb0p0_both_fixture/study_cards/S_onshell_core-w-ht-half_separate_902001'
            cards.mkdir(parents=True)
            with mock.patch.object(precision,'STUDY',root),mock.patch.object(precision,'run') as launch:
                with self.assertRaisesRegex(ValueError,'Incomplete existing'):precision.execute(path,1)
                launch.assert_not_called()
            self.assertEqual(json.loads(path.read_text())['status'],'stopped')

    def test_analysis_averages_flavours_preserves_zero_bins_and_writes_plots(self):
        from test_main_report import TestMainReport
        _,layout,vectors,conditional=TestMainReport().fixture(2)
        cases=precision.cases(5,901001)
        jobs=[dict(case=c) for c in cases]
        def collect(job,**kwargs):
            c=job['case'];key=c['charge'],c['variant'],tuple(c['flavours']);i=c['replica']
            variance=np.broadcast_to(conditional[key][i,:,None],vectors[key][i].shape)
            return layout,{(c['seed'],1)},vectors[key][i],variance,dict(case=c)
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);state=root/'state.json';output=root/'plots'
            state.write_text(json.dumps(dict(status='finished',jobs=jobs,cases=cases,tag='fixture',source_hashes={})))
            with mock.patch.object(precision,'unchanged'),mock.patch.object(precision,'collect_job',side_effect=collect),\
                 mock.patch.object(precision,'assert_matched_physics'),mock.patch.object(precision,'assert_same_layout'),\
                 mock.patch('main_binning.apply_plan',side_effect=lambda l,v,c,p:(l,v,c,{})):
                precision.analyse(state,output,1)
            m=json.loads((output/'manifest.json').read_text())
            self.assertEqual(m['completed_rounds'],5)
            self.assertEqual(len(m['assets']),16)
            self.assertFalse(m['publication_certified'])
            expected=sum(vectors['plus','S',f][:,3,0].mean() for f in precision.FLAVOURS)
            self.assertAlmostEqual(m['rates']['plus']['one_b']['S_pb']['value'],expected)
            with self.assertRaisesRegex(ValueError,'overwrite'):
                with mock.patch.object(precision,'unchanged'):precision.analyse(state,output,1)


if __name__=='__main__':unittest.main()
