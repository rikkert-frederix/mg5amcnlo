"""Full companion coverage and release/scheduling tests; no MC is launched."""
import collections
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from campaign import STUDY, digest
import run_robustness_campaign as runner


class TestRobustnessCampaign(unittest.TestCase):
    @staticmethod
    def folder(path):
        path = Path(path)
        (path/'inputs').mkdir()
        (path/'inputs/benchmark.json').write_bytes((STUDY/'inputs/benchmark.json').read_bytes())
        return path

    @staticmethod
    def release(folder, record, jobs=()):
        evidence = folder/'synthetic_evidence.txt'
        if not evidence.exists():
            evidence.write_text('SYNTHETIC DRIVER TEST ONLY; not a scientific release.\n')
        release = dict(status=runner.RELEASED, benchmark_sha256=record['benchmark_sha256'],
            campaign_tag=record['tag'], campaign_plan_sha256=record['plan_sha256'], reference_jobs=list(jobs),
            evidence={name:{str(evidence):digest(evidence)} for name in runner.EVIDENCE})
        path = folder/'synthetic_release.json'
        path.write_text(json.dumps(release))
        return path

    def test_all_five_states_are_available_for_every_requested_contrast(self):
        rows = runner.cases()
        self.assertEqual(len(rows),640)
        self.assertEqual({row['seed'] for row in rows},set(range(300001,300641)))
        self.assertEqual(len({runner.export_key(row) for row in rows}),64)
        groups = collections.defaultdict(list)
        available = {(charge,flavour,'onshell',0.) for charge in ('plus','minus') for flavour in runner.FLAVOURS}
        for row in rows:
            key = runner.export_key(row)
            if row['decay_bottom_mass']:
                self.assertIn(key[:-1]+(0.,),available)
            available.add(key)
            groups[key+(row['variant'],)].append(row['replica'])
            self.assertEqual(row['production_sampling'],'w-current' if row['w_treatment']=='all-bw' else 'flat')
            self.assertEqual(runner.run_name(row).endswith('_mb4p8'),row['decay_bottom_mass']==4.8)
        self.assertEqual(len(groups),128)
        self.assertTrue(all(values==list(range(5)) for values in groups.values()))
        self.assertEqual(len(runner.comparisons()),10)
        for row in runner.comparisons():
            self.assertEqual({tuple(f) for f in row['flavours']},set(runner.FLAVOURS))
            for state in (row['reference'],row['target']):
                for flavour in runner.FLAVOURS:
                    self.assertIn((row['charge'],flavour,state['w_treatment'],state['decay_bottom_mass']),available)

    def test_prepare_requires_valid_rounds_and_keeps_frozen_allocation(self):
        for replicas,seed in ((4,300001),(True,300001),(5,0),(5,900000000),(5,1.5)):
            with self.assertRaises(ValueError): runner.cases(replicas,seed)
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_prepare_') as directory:
            folder = self.folder(directory)
            with mock.patch.object(runner,'STUDY',folder), mock.patch.object(runner,'generate') as generate:
                path = runner.prepare('synthetic_v1')
                record = json.loads(path.read_text())
                self.assertEqual(record['status'],runner.PREPARED)
                self.assertFalse(record['jobs'])
                self.assertFalse(record['exports'])
                self.assertEqual(record['max_cores'],64)
                runner.unchanged(record)
                release = folder/'unreleased.json'
                release.write_text('{"status":"technical validation incomplete"}\n')
                with self.assertRaisesRegex(ValueError,'Require validation'):
                    runner.execute(path,release)
                generate.assert_not_called()
                self.assertEqual(json.loads(path.read_text())['status'],runner.PREPARED)
                with self.assertRaises(ValueError): runner.prepare('synthetic_v1')
                for kind in ('settings','cases','sources'):
                    altered = copy.deepcopy(record)
                    if kind=='settings': altered['settings']['production_scale']='fixed'
                    elif kind=='cases': altered['cases'].pop()
                    else: altered['source_hashes'].pop(next(iter(altered['source_hashes'])))
                    with self.assertRaisesRegex(ValueError,'allocation or source inventory'):
                        runner.unchanged(altered)
                altered = copy.deepcopy(record)
                altered['cases'].pop()
                altered['plan_sha256'] = runner.plan_digest(altered)
                with self.assertRaisesRegex(ValueError,'full required coverage'): runner.unchanged(altered)

    def test_release_requires_complete_baselines_and_rejects_real_pilots(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_release_') as directory:
            folder = self.folder(directory)
            with mock.patch.object(runner,'STUDY',folder):
                record = json.loads(runner.prepare('synthetic_v1').read_text())
                path = self.release(folder,record)
                with self.assertRaisesRegex(ValueError,'all eight flavours and both charges'):
                    runner.release_evidence(path,record)
                release = json.loads(path.read_text())
                release['campaign_plan_sha256']='different initial allocation'
                path.write_text(json.dumps(release))
                with self.assertRaisesRegex(ValueError,'Require validation'): runner.release_evidence(path,record)
                controls = json.loads((STUDY/'inputs/rng_pilot_queue_alignment_controls_v1.json').read_text())
                jobs = [job for job in controls['jobs'] if '_onshell_' in job['run']]
                path = self.release(folder,record,jobs)
                # These actual audited controls are single pilots, never five main retrainings.
                with self.assertRaisesRegex(ValueError,'independently retrained matched main'):
                    runner.release_evidence(path,record)
                path = self.release(folder,record)
                release = json.loads(path.read_text())
                release['evidence'].pop('small_mass_limit')
                path.write_text(json.dumps(release))
                with self.assertRaisesRegex(ValueError,'small_mass_limit'): runner.release_evidence(path,record)

    def test_serial_execution_reuses_exports_and_checks_massless_production_first(self):
        # Reduce only this temporary scheduling fixture to one flavour. The
        # separate coverage test above uses the real eight-flavour inventory.
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_execute_') as directory:
            folder = self.folder(directory)
            flavours = (('e','e','mu'),)
            generated, called, identities = [], [], []
            references = {(charge,flavours[0],'onshell',0.):folder/('MAIN_'+charge) for charge in ('plus','minus')}
            def generate(args):
                case = vars(args)
                process = runner.process_path(case,args.tag)
                process.mkdir(parents=True)
                generated.append(process)
            def support(path,process):
                path.write_text('synthetic support fixture only\n')
            def production(reference,candidate,path):
                self.assertIn(reference,list(references.values())+generated)
                self.assertTrue('mb0p0' in reference.name or reference.name.startswith('MAIN_'))
                self.assertIn('mb4p8',candidate.name)
                identities.append((reference,candidate))
                path.write_text('synthetic production identity fixture only\n')
            def run(args):
                self.assertTrue(args.main_run)
                self.assertIsNone(args.grid_reference)
                self.assertEqual(args.production_scale,'core-w-ht-half')
                self.assertEqual(args.ecm,13000.)
                self.assertEqual(args.accuracy,.015)
                called.append(args)
            with mock.patch.object(runner,'STUDY',folder), mock.patch.object(runner,'FLAVOURS',flavours), \
                 mock.patch.object(runner,'release_evidence',return_value=({}, {(1,1)},dict(references))), \
                 mock.patch.object(runner,'generate',side_effect=generate), \
                 mock.patch.object(runner,'check_support',side_effect=support), \
                 mock.patch.object(runner,'benchmark_param_card'), \
                 mock.patch.object(runner,'check_production',side_effect=production), \
                 mock.patch.object(runner,'run',side_effect=run), mock.patch.object(runner,'audit'), \
                 mock.patch.object(runner,'harvest'), mock.patch.object(runner,'read'), \
                 mock.patch.object(runner,'audit_virtuals'), \
                 mock.patch.object(runner,'verified_pairs',side_effect=lambda job:{(job['case']['seed'],2)}):
                path = runner.prepare('synthetic_execution')
                record = json.loads(path.read_text())
                release = self.release(folder,record)
                runner.execute(path,release)
                final = json.loads(path.read_text())
                self.assertEqual(final['status'],runner.DONE)
                self.assertEqual(len(final['jobs']),80)
                self.assertEqual(len(generated),8)
                self.assertEqual(len(identities),4)
                self.assertEqual(len(called),80)
                self.assertEqual([job['case'] for job in final['jobs']],record['cases'])
                self.assertEqual(final['reference_stage_initializations'],1)
                self.assertEqual(final['companion_stage_initializations'],80)
                with self.assertRaisesRegex(ValueError,'never restart'): runner.execute(path,release)

    def test_stream_overlap_stops_without_marking_failed_job_complete(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_overlap_') as directory:
            folder = self.folder(directory)
            references = {(charge,tuple(flavour),'onshell',0.):folder/('MAIN_'+charge)
                          for charge in ('plus','minus') for flavour in runner.FLAVOURS}
            def support(path,process): path.write_text('synthetic support fixture only\n')
            with mock.patch.object(runner,'STUDY',folder), \
                 mock.patch.object(runner,'release_evidence',return_value=({}, {(1,1)},references)), \
                 mock.patch.object(runner,'generate'), mock.patch.object(runner,'check_support',side_effect=support), \
                 mock.patch.object(runner,'run') as run, mock.patch.object(runner,'audit'), \
                 mock.patch.object(runner,'harvest'), mock.patch.object(runner,'read'), \
                 mock.patch.object(runner,'audit_virtuals'), mock.patch.object(runner,'verified_pairs',return_value={(1,1)}):
                path = runner.prepare('synthetic_overlap')
                release = self.release(folder,json.loads(path.read_text()))
                with self.assertRaisesRegex(ValueError,'streams overlap'): runner.execute(path,release)
                final = json.loads(path.read_text())
                self.assertEqual(final['status'],'stopped')
                self.assertFalse(final['jobs'])
                self.assertIn('streams overlap',final['error'])
                self.assertEqual(run.call_count,1)
                with self.assertRaisesRegex(ValueError,'never restart'): runner.execute(path,release)


if __name__=='__main__':
    unittest.main()
