import json
import math
import tempfile
import unittest
from pathlib import Path

import analyse_tail_pilot as analysis


TRACE = '''POINT 3 1 2 2 0.0 2.0000000000000000D+16
X .25 .75
LINE 1 1 2 4 1.0D+16 1 0 0 .1 .2 100 100
LINE 2 4 2 4 1.0D+00 1 0 0 .1 .2 100 100
LINE 3 5 2 4 -1.0D+16 -1 0 0 .1 .2 100 100
PRODUCTION 4 8.0 1.0D-5 2.0D-5 1.0D-4 .1 1.0D+6
END
SUMMARY 4 0.0 2.0 2.0 1 1 0
MOMENT 2 4 0.0 2.0 2.0
TYPEMOMENT 2 1 4 0.0 2.0 2.0
'''


class TestTailAnalysis(unittest.TestCase):
    def test_precise_sum_exposes_native_cancellation(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'record_breakers.dat'
            path.write_text(TRACE)
            result = analysis.parse_trace(path)
        point = result['points'][0]
        self.assertEqual(point['diagnostics']['sum'], 1.0)
        self.assertEqual(point['diagnostics']['native_vs_fsum_absolute'], 1.0)
        self.assertEqual(point['diagnostics']['production_local_subtraction']['sum'], 1.0)
        self.assertEqual(point['production_coordinates']['xi'], 1e-5)
        self.assertEqual(result['summaries'][0]['points'], 4)
        self.assertEqual(result['component_moments'][1]['type'], 1)

    def test_incomplete_record_and_invalid_coordinates(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'record_breakers.dat'
            path.write_text('POINT 1 1 1 2 0 0\nX .5 .5\n')
            result = analysis.parse_trace(path)
            self.assertTrue(result['incomplete_trailing_point'])
            self.assertEqual(result['points'], [])
            path.write_text('POINT 1 1 1 2 0 0\nX .5\nEND\n')
            with self.assertRaisesRegex(ValueError, 'X length'):
                analysis.parse_trace(path)

    def test_components_and_zero_denominator(self):
        result = analysis.conditioning([1., -1.])
        self.assertIsNone(result['cancellation_ratio'])
        self.assertTrue(result['exact_zero_with_nonzero_terms'])
        cleaned = analysis.clean({'a': float('inf'), 'b': [float('nan')]})
        self.assertEqual(cleaned, {'a': None, 'b': [None]})
        json.dumps(cleaned, allow_nan=False)

    def test_report_replay_and_full_moments(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ('first', 'replay'):
                path = root / name
                path.mkdir()
                path.joinpath('record_breakers.dat').write_text(TRACE)
                path.joinpath('stdout.log').write_text('Maximum number of iterations is: 1\nFinal result: 1D-4 +/- 2D-5\n')
            root.joinpath('replay/manifest.json').write_text(json.dumps({'seed': 123, 'replay_of': 'first'}))
            report = analysis.make_report(root)
        self.assertEqual(report['worker_count'], 2)
        self.assertTrue(report['replays'][0]['saved_records_equal'])
        self.assertTrue(report['replays'][0]['summaries_equal'])
        self.assertEqual(report['workers'][0]['native']['error_pb'], 2e-5)
        self.assertAlmostEqual(report['workers'][0]['all_point_moments']['independent_point_standard_error'],
                               math.sqrt(2. / 12.))

    def test_multiplicative_leaves(self):
        trace = '''POINT 5 1 2 1 2.0 4.0
VOLUME .125
X .5
LEAF 0 3.0
LEAF 1 -1.0
END
'''
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'record_breakers.dat'
            path.write_text(trace)
            result = analysis.parse_trace(path)
        point = result['points'][0]
        self.assertEqual(point['volume'], .125)
        self.assertEqual(point['diagnostics']['sum'], 2.)
        self.assertEqual(point['diagnostics']['cancellation_ratio'], 2.)
        self.assertEqual(point['diagnostics']['decomposition'], 'multiplicative_leaves')

    def test_campaign_metadata_excludes_smoke_and_replay_from_ensembles(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            allocations = [('smoke', 1., 1, False), ('main', 2., 2, False),
                           ('main', 4., 3, False), ('replay', 2., 2, True)]
            executions = []
            for ordinal, (label, value, seed, replay) in enumerate(allocations):
                path = root / str(ordinal)
                path.mkdir()
                path.joinpath('record_breakers.dat').write_text(TRACE)
                path.joinpath('log.txt').write_text(f'Final result: {value} +/- 1.0\n')
                executions.append(dict(label=label, same_stream_replay=replay,
                    jobs=[dict(directory=str(path), seed=seed, variant='S', points=4)]))
            campaign = root / 'campaign.json'
            campaign.write_text(json.dumps(dict(executions=executions)))
            report = analysis.make_report(root, campaign_record=campaign)
        ensemble = report['independent_ensembles'][0]
        self.assertEqual(ensemble['workers'], 2)
        self.assertEqual(ensemble['native_mean_pb'], 3.)
        self.assertEqual(ensemble['independent_scatter_error_pb'], 1.)
        self.assertEqual(ensemble['seeds'], [2, 3])
        self.assertEqual(len(report['replays']), 1)

    def test_training_excluded_and_iterations_separate_ensembles(self):
        rows = []
        for label, iterations in [('S_training_grid', 4), ('warmup', 2), ('main', 1), ('main', 2)]:
            rows.append(dict(metadata={'campaign_job': dict(execution_label=label,
                variant='S', points=100, iterations=iterations, seed=iterations)},
                native=dict(value_pb=1., error_pb=.1)))
        ensembles = analysis.independent_ensembles(rows)
        self.assertEqual([row['role'] for row in rows[:2]], ['training', 'training'])
        self.assertEqual(len(ensembles), 2)
        self.assertEqual({e['requested_iterations'] for e in ensembles}, {1, 2})

    def test_native_virtual_cpu_and_efficiency(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            root.joinpath('log.txt').write_text('accumulated results V  1 = .3D-2 +/- .2D-2 ( 10 % )\n'
                'Final result: 1.0 +/- .1\nTime spent in Total : 20.0\n')
            native = analysis.native_log(root)
        self.assertEqual(native['production_virtual']['error_pb'], .002)
        self.assertEqual(native['cpu_seconds']['Total'], 20.)
        rows = [dict(metadata={'campaign_job': dict(execution_label='main',variant='S',points=100)},
            native=native) for _ in range(2)]
        ensemble = analysis.independent_ensembles(rows)[0]
        self.assertAlmostEqual(ensemble['native_variance_times_cpu_pb2_s'], .2)
        self.assertAlmostEqual(ensemble['production_virtual_variance_times_cpu_pb2_s'], .00008)


if __name__ == '__main__':
    unittest.main()
