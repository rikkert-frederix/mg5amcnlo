"""Complete synthetic collection and actual archived W/mass physics checks.

Only the collection fixture mocks physical-file validators. A separate test
uses the real twenty pilot archives and four production-identity records.
No synthetic vector is placed in the physical result directory.
"""
import contextlib
import copy
import io
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from campaign import STUDY, digest, save
from full_flavour_statistics import FLAVOURS
from load_results import WEIGHTS, load
from robustness_statistics import LABELS, estimate, group_key, select
from run_robustness_campaign import DONE, STATES
from test_w_treatment_pilot_report import pilot_jobs
import robustness_replicas as collector


class TestRobustnessReplicas(unittest.TestCase):
    @staticmethod
    def fixture(folder):
        queue = json.loads((STUDY/'inputs/robustness_campaign_fullflavour_v1.json').read_text())
        queue.update(status=DONE, jobs=[], exports={})
        titles = ['W+ rates: synthetic', 'W- rates: synthetic']
        edges = np.tile(np.column_stack([np.arange(9), np.arange(1, 10)]), (2, 1))
        layout = dict(titles=titles, edges=edges, offsets=np.asarray([0, 9, 18]), weights=WEIGHTS)
        references = [dict(replica=r, charge=c, flavours=list(f), variant=v,
                           w_treatment='onshell', decay_bottom_mass=0., seed=100001+i)
            for i, (r, c, f, v) in enumerate(itertools.product(range(5), ('plus', 'minus'), FLAVOURS, ('S', 'Pi')))]
        results, streams, reference_jobs, processes = {}, {}, [], {}
        source = next(iter(queue['source_hashes']))
        spread = np.asarray([.8, .9, 1., 1.1, 1.2])
        proof = folder/'synthetic_export_proof.json'
        save(proof, dict(scope='Temporary synthetic collector fixture; physical validators are mocked'))
        for index, (origin, case) in enumerate([('reference', row) for row in references]+
                                              [('companion', row) for row in queue['cases']]):
            mode, mass, charge, variant, flavour = (case['w_treatment'], case['decay_bottom_mass'],
                case['charge'], case['variant'], tuple(case['flavours']))
            state = STATES.index((mode, mass))
            context = mode, mass, charge, flavour
            process = processes.setdefault(context, folder/('synthetic_process_%d' % len(processes)))
            run = 'synthetic_%d' % index
            amplitude = (1.+FLAVOURS.index(flavour))*(1.+.1*state)*spread[case['replica']]
            if variant == 'Pi': amplitude *= 1.05+.01*state
            values = np.zeros((18, len(WEIGHTS)))
            first = 0 if charge == 'plus' else 9
            values[first:first+9] = amplitude*np.arange(1., 10.)[:, None]
            args = dict(charge=charge, flavours=list(flavour), w_treatment=mode, decay_bottom_mass=mass,
                        corrected='both', tag=queue['tag'] if origin == 'companion' else 'synthetic_main')
            manifest = dict(variant=variant, w_treatment=mode, decay_bottom_mass=mass,
                settings=dict(iseed=case['seed'], lhaid=331700), production_scale='core-w-ht-half',
                production_sampling='w-current' if mode == 'all-bw' else 'flat', decay_scale_grouping='separate',
                top_width_lo=1.5+state*.01, top_width_nlo=1.3+state*.01, w_width=2.1,
                top_width_reference_scale=172.5, top_width_w_treatment='onshell' if mode == 'onshell' else 'bw',
                top_width_bottom_mass=mass, hashes={'analysis_source':'synthetic', 'param_card.dat':'synthetic'})
            execution = dict(pilot=False, grid_reference=None, export_generation=dict(args=args),
                source_hashes={source: queue['source_hashes'][source]})
            report = dict(path=str(process/'Events'/run/'MADatNLO.HwU'), variant=variant, manifest=manifest, execution=execution)
            audit_path = folder/(run+'_audit.json')
            save(audit_path, report)
            results[str(audit_path)] = dict(layout, report=report, values=values, parameter_signature='synthetic_%d' % state)
            batch_path = folder/(run+'_batches.npz')
            np.savez_compressed(batch_path, **layout, strata=np.zeros(5, dtype=int),
                contributions=np.asarray([.16, .18, .2, .22, .24])[:, None, None]*values[None])
            save(batch_path.with_suffix('.json'), dict(refinement_rng_audit=dict(distinct_initialization_pairs=1)))
            job = dict(case=case, process=str(process), run=run, variant=variant, audit=str(audit_path), batches=str(batch_path))
            streams[run] = {(1000+index, 2000+index)}
            if origin == 'companion':
                exported = dict(process=str(process), local_support=str(proof), local_support_sha256=digest(proof))
                if mass:
                    exported.update(production_identity=str(proof), production_identity_sha256=digest(proof))
                queue['exports'][process.name] = exported
                job['export_checks'] = exported
                queue['jobs'].append(job)
            else:
                reference_jobs.append(job)
        release = dict(reference_jobs=reference_jobs, scope='Synthetic fixture only')
        release_path = folder/'synthetic_release.json'
        save(release_path, release)
        queue.update(release=str(release_path), release_sha256=digest(release_path))
        queue_path = folder/'synthetic_queue.json'
        save(queue_path, queue)
        reference_streams = set().union(*(streams[job['run']] for job in reference_jobs))
        return queue_path, results, streams, release, reference_streams

    def patches(self, fixture):
        queue, results, streams, release, reference_streams = fixture
        stack = contextlib.ExitStack()
        stack.enter_context(mock.patch.object(collector, 'load', side_effect=lambda path: copy.deepcopy(results[str(path)])))
        stack.enter_context(mock.patch.object(collector, 'release_evidence', return_value=(release, reference_streams.copy(), {})))
        stack.enter_context(mock.patch.object(collector, 'verified_pairs', side_effect=lambda job: streams[job['run']]))
        stack.enter_context(mock.patch.object(collector, 'generated_evidence'))
        stack.enter_context(mock.patch.object(collector, 'validate_comparisons', return_value=16))
        stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
        return stack

    def test_complete_eight_hundred_run_collection_and_shared_references(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_vectors_') as directory:
            folder = Path(directory)
            fixture = self.fixture(folder)
            output = folder/'synthetic_vectors.json'
            with self.patches(fixture):
                collector.collect(fixture[0], output)
            record, layout, arrays, conditional = collector.read(output)
            self.assertEqual(len(record['samples']), 800)
            self.assertEqual(len(arrays), 160)
            self.assertEqual(sum(key[:2] == STATES[0] for key in arrays), 32)
            self.assertEqual(record['reference_stage_initializations'], 160)
            self.assertEqual(record['companion_stage_initializations'], 640)
            self.assertEqual(layout['weights'], WEIGHTS)
            self.assertTrue(all(row.shape == (5, 18, 183) for row in arrays.values()))
            comparison = record['comparisons'][:1]
            subset = {key: row[:, :1, :1] for key, row in select(arrays, comparison).items()}
            result = estimate(subset, comparison)
            self.assertAlmostEqual(result['value'][0, 0, 0, 0], 36.)
            self.assertAlmostEqual(result['value'][0, LABELS.index('reference_Pi_over_S'), 0, 0], 1.05)
            expected = sum(i*i for i in range(1, 9))*np.var([.8, .9, 1., 1.1, 1.2], ddof=1)/5
            self.assertAlmostEqual(result['mc_error'][0, 0, 0, 0]**2, expected, places=12)
            self.assertGreater(sum(row.sum() for row in conditional.values()), 0.)
            with self.assertRaisesRegex(ValueError, 'overwrite'):
                collector.collect(fixture[0], output)
            broken = json.loads(fixture[0].read_text())
            broken['jobs'].pop()
            with self.assertRaisesRegex(ValueError, 'retrainings'):
                collector.inventory(broken)

    def test_stream_overlap_stops_before_a_collected_archive(self):
        with tempfile.TemporaryDirectory(prefix='synthetic_robustness_overlap_') as directory:
            folder = Path(directory)
            fixture = self.fixture(folder)
            queue, results, streams, release, reference_streams = fixture
            first = json.loads(queue.read_text())['jobs'][0]
            streams[first['run']] = {next(iter(reference_streams))}
            output = folder/'rejected_vectors.json'
            with self.patches(fixture):
                with self.assertRaisesRegex(ValueError, 'share actual stage streams'):
                    collector.collect(queue, output)
            self.assertFalse(output.exists())
            self.assertFalse(output.with_suffix('.npz').exists())

    def test_actual_twenty_pilot_archives_cover_all_physical_edges(self):
        queue = json.loads((STUDY/'inputs/width_mass_pilot_queue_width_mass_v4.json').read_text())
        jobs = [job for charge in ('plus', 'minus') for job in pilot_jobs(charge)]
        jobs += [job for job in queue['jobs'] if job['case']['decay_bottom_mass'] == 4.8]
        representatives = {}
        for job in jobs:
            result = load(job['audit'])
            collector.generated_evidence(result)
            generation = result['report']['execution']['export_generation']['args']
            key = group_key(generation['w_treatment'], generation['decay_bottom_mass'], generation['charge'],
                            result['report']['variant'], generation['flavours'])
            representatives[key] = {name: result[name] for name in
                ('report', 'parameter_signature', 'titles', 'edges', 'offsets', 'weights')}
        identities = {row['process']: row['production_identity'] for row in queue['exports'].values()
                      if 'production_identity' in row}
        self.assertEqual(len(representatives), 20)
        self.assertEqual(collector.validate_comparisons(representatives, identities, full_flavour=False), 2)
        with self.assertRaisesRegex(ValueError, 'both charges and all eight'):
            collector.validate_comparisons(representatives, identities)
        with self.assertRaisesRegex(ValueError, 'Missing massive'):
            collector.validate_comparisons(representatives, {}, full_flavour=False)
        changed = copy.deepcopy(representatives)
        key = next(key for key in changed if key[:2] == ('all-bw', 4.8))
        changed[key]['report']['manifest']['top_width_nlo'] = 1.
        with self.assertRaises(ValueError):
            collector.validate_comparisons(changed, identities, full_flavour=False)


if __name__ == '__main__':
    unittest.main()
