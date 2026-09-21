"""Parameter selection/release checks; never launch generation or integration."""
import collections
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from campaign import STUDY, digest
import run_parameter_campaign as runner
from parameter_variation_inputs import SCENARIOS


class TestParameterCampaign(unittest.TestCase):
    @staticmethod
    def selection():
        return dict(schema='parameter_scan_selection_v1', reason='Synthetic coverage test, not a physical allocation',
            contexts=[dict(scenario=name, charge=charge, flavours=list(flavours), w_treatment='onshell')
                for name, charge, flavours in itertools.product([row[0] for row in SCENARIOS[1:]],
                    ('plus', 'minus'), itertools.product(('e', 'mu'), repeat=3))])

    def test_full_requested_family_coverage_and_independent_retraining_rounds(self):
        rows = runner.cases(self.selection())
        self.assertEqual(len(rows), 7*2*8*2*5)
        self.assertEqual(len({row['seed'] for row in rows}), len(rows))
        groups = collections.defaultdict(list)
        for row in rows:
            groups[row['scenario'], row['charge'], tuple(row['flavours']), row['variant']].append(row['replica'])
        self.assertEqual(len(groups), 7*2*8*2)
        self.assertTrue(all(replicas == list(range(5)) for replicas in groups.values()))
        for replica in range(5):
            self.assertEqual(len([row for row in rows if row['replica'] == replica]), 7*2*8*2)

    def test_invalid_selection_retraining_or_seed_is_rejected(self):
        selection = self.selection()
        for replicas, seed in ((4, 200001), (5, 0), (5, 900000000)):
            with self.assertRaises(ValueError): runner.cases(selection, replicas, seed)
        selection['contexts'].append(selection['contexts'][0])
        with self.assertRaises(ValueError): runner.cases(selection)
        selection = self.selection()
        selection['contexts'][0]['scenario'] = 'unmatched_widths'
        with self.assertRaises(ValueError): runner.cases(selection)

    def test_prepare_does_not_release_or_launch_and_case_edits_are_detected(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_parameter_selection_') as folder:
            folder = Path(folder)
            (folder/'inputs').mkdir()
            (folder/'scripts').symlink_to(STUDY/'scripts', target_is_directory=True)
            (folder/'inputs/benchmark.json').write_bytes((STUDY/'inputs/benchmark.json').read_bytes())
            selection = self.selection()
            selection['contexts'] = selection['contexts'][:1]
            selection_path = folder/'selection.json'
            selection_path.write_text(json.dumps(selection))
            inputs = STUDY/'inputs/parameter_variation_inputs_bounded_v1.json'
            with mock.patch.object(runner, 'STUDY', folder), mock.patch.object(runner, 'generate') as generate:
                record_path = runner.prepare('synthetic_test', selection_path, inputs)
                record = json.loads(record_path.read_text())
                self.assertEqual(record['status'], runner.PREPARED)
                self.assertEqual(len(record['cases']), 10)
                self.assertFalse(record['jobs'])
                release = folder/'unreleased.json'
                release.write_text(json.dumps(dict(status='main precision still unresolved')))
                with self.assertRaisesRegex(ValueError, 'main S/Pi convergence'):
                    runner.execute(record_path, release)
                generate.assert_not_called()
                self.assertEqual(json.loads(record_path.read_text())['status'], runner.PREPARED)
                with self.assertRaises(ValueError): runner.prepare('synthetic_test', selection_path, inputs)
                record['cases'][0]['variant'] = 'P'
                record_path.write_text(json.dumps(record))
                with self.assertRaisesRegex(ValueError, 'cases differ'):
                    runner.execute(record_path, release)
                generate.assert_not_called()

    def test_release_rejects_insufficient_or_repeated_physical_baseline_references(self):
        controls = json.loads((STUDY/'inputs/rng_pilot_queue_alignment_controls_v1.json').read_text())
        references = [job for job in controls['jobs'] if '_onshell_' in job['run']]
        selection = dict(schema='parameter_scan_selection_v1', reason='Synthetic release-rejection test',
            contexts=[dict(scenario='as117', charge='plus', flavours=['e', 'e', 'mu'], w_treatment='onshell')])
        record = dict(cases=runner.cases(selection), benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
                      selection_sha256='synthetic selection only')
        with tempfile.TemporaryDirectory(prefix='synthetic_parameter_release_') as directory:
            folder = Path(directory)
            evidence = folder/'synthetic_evidence.txt'
            evidence.write_text('SYNTHETIC RELEASE-REJECTION TEST; not a main convergence inspection.\n')
            release = dict(status=runner.RELEASED, benchmark_sha256=record['benchmark_sha256'],
                selection_sha256=record['selection_sha256'], reference_jobs=references,
                evidence={name: {str(evidence): digest(evidence)} for name in
                          ('main_convergence_inspection', 'main_reduction')})
            path = folder/'rejected_release.json'
            path.write_text(json.dumps(release))
            with self.assertRaisesRegex(ValueError, 'at least five'):
                runner.release_evidence(path, record)
            release['reference_jobs'] = references*5
            path.write_text(json.dumps(release))
            with self.assertRaisesRegex(ValueError, 'reuse actual stage streams'):
                runner.release_evidence(path, record)


if __name__ == '__main__':
    unittest.main()
