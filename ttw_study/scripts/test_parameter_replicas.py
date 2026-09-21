"""Synthetic collection/reduction test with mocked physical-file validators.

The input/card/HwU/stream validators are tested separately on physical
controls. Here temporary complete-run arrays test inventory, grouping,
baseline reuse and the separation of conditional and retraining variances.
"""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from campaign import STUDY, digest, save
from load_results import WEIGHTS
from parameter_statistics import LABELS, estimate
from parameter_variation_results import vectors
from run_parameter_campaign import DONE, cases
import parameter_replicas as collector
from test_parameter_variation_results import synthetic_hwu


class TestParameterReplicas(unittest.TestCase):
    def fixture(self, folder):
        selection = dict(schema='parameter_scan_selection_v1', reason='Synthetic collection test only',
            contexts=[dict(scenario=scenario, charge='plus', flavours=list(flavour), w_treatment='onshell')
                for scenario in ('as117', 'as119') for flavour in (('e', 'e', 'e'), ('e', 'e', 'mu'))])
        selection_path = folder/'selection.json'
        save(selection_path, selection)
        inputs = STUDY/'inputs/parameter_variation_inputs_bounded_v1.json'
        queue = dict(status=DONE, tag='synthetic_scan', cases=cases(selection), jobs=[], replicas=5, seed_start=200001,
            selection=str(selection_path), selection_sha256=digest(selection_path), inputs=str(inputs),
            inputs_sha256=digest(inputs), benchmark_sha256=digest(STUDY/'inputs/benchmark.json'), source_hashes={})
        layout = vectors(synthetic_hwu(1.))
        results, streams, references = {}, {}, []
        spread = np.asarray([.94, .97, 1., 1.03, 1.06])
        reference_cases = [dict(scenario='baseline', charge='plus', flavours=list(flavour), w_treatment='onshell',
                                replica=r, variant=v, seed=300001+i, decay_bottom_mass=0.)
            for i, (r, flavour, v) in enumerate((r, flavour, v) for r in range(5)
                for flavour in (('e', 'e', 'e'), ('e', 'e', 'mu')) for v in ('S', 'Pi'))]
        for index, (origin, case) in enumerate([('reference', row) for row in reference_cases]+
                                             [('scan', row) for row in queue['cases']]):
            scenario, variant = case['scenario'], case['variant']
            amplitude = (1. if case['flavours'][-1] == 'e' else 4.)*spread[case['replica']]
            amplitude *= dict(baseline=1., as117=1.05, as119=.95)[scenario]
            if variant == 'Pi': amplitude *= dict(baseline=1.2, as117=1.224, as119=1.18)[scenario]
            values = amplitude*layout['values']
            if variant == 'S':
                values = values.copy()
                for h, title in enumerate(layout['titles']):
                    if ' rates:' in title:
                        a = layout['offsets'][h]
                        values[a+6] += values[a+7]
                        values[a+7:a+9] = 0.
            if origin == 'reference':
                values = np.concatenate([values, np.repeat(values[:, :1], 101, axis=1)], axis=1)
            weights = WEIGHTS if origin == 'reference' else WEIGHTS[:82]
            args = dict(charge='plus', flavours=case['flavours'], w_treatment='onshell', decay_bottom_mass=0.,
                        corrected='both', tag='synthetic_scan_'+scenario)
            manifest = dict(variant=variant, w_treatment='onshell', decay_bottom_mass=0.,
                settings=dict(iseed=case['seed'], lhaid={'baseline':331700, 'as117':333900, 'as119':334100}[scenario]),
                top_width_lo=1.5, top_width_nlo=1.3, w_width=2.1, top_width_reference_scale=172.5,
                top_width_w_treatment='onshell', top_width_bottom_mass=0., production_scale='core-w-ht-half',
                decay_scale_grouping='separate', hashes={'analysis_source':'synthetic measurement', 'param_card.dat':'synthetic hash'},
                parameter_scenario=dict(name=scenario))
            execution = dict(export_generation=dict(args=args), pilot=False, grid_reference=None,
                             source_hashes={}, parameter_case=case)
            report = dict(variant=variant, manifest=manifest, execution=execution)
            audit_path = folder/('synthetic_%d_audit.json' % index)
            save(audit_path, report)
            result = dict(layout, values=values, weights=weights, report=report, parameter_signature='synthetic_'+scenario)
            results[str(audit_path)] = result
            batch_path = folder/('synthetic_%d_batches.npz' % index)
            batch_weights = np.asarray([.16, .18, .20, .22, .24])
            np.savez_compressed(batch_path, contributions=batch_weights[:, None, None]*values[None],
                strata=np.zeros(5, dtype=int), edges=layout['edges'], offsets=layout['offsets'],
                titles=layout['titles'], weights=weights)
            save(batch_path.with_suffix('.json'), dict(refinement_rng_audit=dict(distinct_initialization_pairs=1)))
            job = dict(case=case, variant=variant, audit=str(audit_path), batches=str(batch_path),
                       run='synthetic_%d' % index, process=str(folder))
            streams[job['run']] = {(1000+index, 2000+index)}
            (references if origin == 'reference' else queue['jobs']).append(job)
        release = dict(reference_jobs=references, data_scope='synthetic fixture only')
        release_path = folder/'synthetic_release.json'
        save(release_path, release)
        queue.update(release=str(release_path), release_sha256=digest(release_path))
        queue_path = folder/'synthetic_queue.json'
        save(queue_path, queue)
        reference_streams = set().union(*(streams[job['run']] for job in references))
        return queue_path, results, streams, release, reference_streams

    def test_complete_collection_preserves_shared_baselines_and_retraining_variance(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_parameter_vectors_') as directory:
            folder = Path(directory)
            queue, results, streams, release, reference_streams = self.fixture(folder)
            output = folder/'synthetic_vectors.json'
            loader = lambda path: copy.deepcopy(results[str(path)])
            with mock.patch.object(collector, 'load_baseline', side_effect=loader), \
                 mock.patch.object(collector, 'load_scan', side_effect=loader), \
                 mock.patch.object(collector, 'release_evidence', return_value=(release, reference_streams.copy())), \
                 mock.patch.object(collector, 'verified_pairs', side_effect=lambda job: streams[job['run']]):
                collector.collect(queue, output)
            report, layout, ensembles, conditional = collector.read(output)
            self.assertEqual(len(report['samples']), 60)
            self.assertEqual(len(ensembles), 12)
            self.assertEqual(sum(key[0] == 'reference' for key in ensembles), 4)
            self.assertEqual(report['reference_stage_initializations'], 20)
            self.assertEqual(report['scan_stage_initializations'], 40)
            self.assertTrue(all(value.shape == (5, 290, 82) for value in ensembles.values()))
            reduced = {key: value[:, [0, 3, 4]] for key, value in ensembles.items()}
            result = estimate(reduced, report['comparisons'], full_flavour=False)
            self.assertAlmostEqual(result['value'][0, LABELS.index('reference_S'), 0, 0], 500.)
            self.assertAlmostEqual(result['value'][0, LABELS.index('reference_Pi_over_S'), 0, 0], 1.2)
            self.assertAlmostEqual(result['value'][0, LABELS.index('product_ratio_change'), 0, 0], .024)
            expected_variance = (100.**2+400.**2)*np.var([.94, .97, 1., 1.03, 1.06], ddof=1)/5
            actual = result['mc_error'][0, LABELS.index('reference_S'), 0, 0]**2
            self.assertAlmostEqual(actual, expected_variance, places=9)
            diagnostic = sum(values[:, 0].sum()/25 for key, values in conditional.items()
                             if key[0] == 'reference' and key[4] == 'S')
            self.assertGreater(diagnostic, 0.)
            self.assertNotAlmostEqual(actual, expected_variance+diagnostic)
            with self.assertRaises(ValueError): collector.collect(queue, output)
            # A complete-looking queue with one dropped run fails before any
            # physical-output loader can bless an incomplete ensemble.
            broken = json.loads(queue.read_text())
            broken['jobs'].pop()
            with self.assertRaisesRegex(ValueError, 'retrainings'): collector.inventory(broken)

    def test_scan_overlap_with_baseline_is_rejected_before_collection_output(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_parameter_overlap_') as directory:
            folder = Path(directory)
            queue, results, streams, release, reference_streams = self.fixture(folder)
            output = folder/'rejected_vectors.json'
            loader = lambda path: copy.deepcopy(results[str(path)])
            with mock.patch.object(collector, 'load_baseline', side_effect=loader), \
                 mock.patch.object(collector, 'load_scan', side_effect=loader), \
                 mock.patch.object(collector, 'release_evidence', return_value=(release, reference_streams.copy())), \
                 mock.patch.object(collector, 'verified_pairs', return_value={next(iter(reference_streams))}):
                with self.assertRaisesRegex(ValueError, 'share actual stage streams'):
                    collector.collect(queue, output)
            self.assertFalse(output.exists())
            self.assertFalse(output.with_suffix('.npz').exists())


if __name__ == '__main__':
    unittest.main()
