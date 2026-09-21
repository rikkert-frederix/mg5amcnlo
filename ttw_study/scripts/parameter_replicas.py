#!/usr/bin/env python3
"""Collect complete scan retrainings and their shared physical baselines.

This keeps one vector per independently retrained integration. Conditional
within-run variances are diagnostics only; physics errors come from the
spread of complete run vectors and must not receive that noise twice.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from load_results import assert_matched_physics, assert_same_layout, load as load_baseline
from main_replicas import collapse_batches
from parameter_statistics import group_key
from parameter_variation_results import SCALE_WEIGHTS, load as load_scan
from run_parameter_campaign import DONE, cases, release_evidence
from run_width_mass_pilots import verified_pairs

COLLECTED = 'complete independent parameter vectors collected; comparisons and convergence required'


def inventory(queue):
    if queue.get('status') != DONE or not queue.get('cases'):
        raise ValueError('Require a complete selected parameter campaign')
    if (digest(queue['selection']) != queue['selection_sha256']
            or digest(queue['inputs']) != queue['inputs_sha256']
            or digest(STUDY/'inputs/benchmark.json') != queue['benchmark_sha256']):
        raise ValueError('Changed parameter selection or physical inputs')
    expected = cases(json.loads(Path(queue['selection']).read_text()), queue['replicas'], queue['seed_start'])
    if queue['cases'] != expected or [job['case'] for job in queue['jobs']] != expected:
        raise ValueError('Missing, repeated or reordered parameter retrainings')
    return queue['replicas']


def collect(queue_path, output):
    queue_path, output = Path(queue_path).resolve(), Path(output).resolve()
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite parameter retraining vectors')
    queue = json.loads(queue_path.read_text())
    count = inventory(queue)
    release_path = Path(queue['release'])
    if digest(release_path) != queue['release_sha256']:
        raise ValueError('Parameter-campaign release changed')
    release, streams = release_evidence(release_path, queue)
    reference_count = len(streams)
    groups, samples, physics, comparisons = {}, [], {}, {}
    layout = None
    jobs = [('reference', job) for job in release['reference_jobs']]+[('scan', job) for job in queue['jobs']]
    for origin, job in jobs:
        result = (load_baseline if origin == 'reference' else load_scan)(job['audit'])
        report = result['report']
        manifest, execution = report['manifest'], report['execution']
        generated = execution['export_generation']['args']
        charge, flavour, mode = generated['charge'], tuple(generated['flavours']), manifest['w_treatment']
        variant = job['variant']
        if (report['variant'] != variant or variant not in ('S', 'Pi')
                or generated['w_treatment'] != mode or generated['decay_bottom_mass'] != 0.
                or generated['corrected'] != 'both' or execution.get('grid_reference') is not None):
            raise ValueError('Unmatched process or imported grid in the parameter comparison')
        if origin == 'scan':
            case = job['case']
            scenario = case['scenario']
            expected = dict(charge=case['charge'], flavours=case['flavours'],
                w_treatment=case['w_treatment'], decay_bottom_mass=0., corrected='both',
                tag=queue['tag']+'_'+scenario)
            if (execution['pilot'] or any(generated[key] != value for key, value in expected.items())
                    or case['variant'] != variant or manifest['settings']['iseed'] != case['seed']
                    or execution['parameter_case'] != case or manifest['parameter_scenario']['name'] != scenario
                    or execution['source_hashes'] != queue['source_hashes']):
                raise ValueError('Scan run differs from its prepared context or frozen source')
            pairs = verified_pairs(job)
            if streams & pairs:
                raise ValueError('Scan/reference retrainings share actual stage streams')
            streams.update(pairs)
            stage_count = len(pairs)
            comparison = comparisons.setdefault((scenario, mode, charge), set())
            comparison.add(flavour)
        else:
            scenario = 'baseline'
            # release_evidence already checked all reference archives, including
            # pair uniqueness across references. Do not count them twice.
            metadata = json.loads(Path(job['batches']).with_suffix('.json').read_text())
            stage_count = metadata['refinement_rng_audit']['distinct_initialization_pairs']
        result['values'] = result['values'][:, :82]
        result['weights'] = list(SCALE_WEIGHTS)
        if layout is None:
            layout = {key: result[key] for key in ('titles', 'edges', 'offsets', 'weights')}
        else:
            assert_same_layout([layout, result])
        projected = dict(result, report=copy.deepcopy(report))
        projected['report']['execution']['export_generation']['args'].update(charge='plus', flavours=['e', 'e', 'e'])
        physics_key = origin, scenario, mode
        if physics_key in physics:
            assert_matched_physics([physics[physics_key], projected])
        else:
            # Only the physical signature is used by assert_matched_physics.
            physics[physics_key] = {key: projected[key] for key in ('report', 'parameter_signature')}
        for h, title in enumerate(result['titles']):
            if ('W-' if charge == 'plus' else 'W+') in title:
                a, b = map(int, result['offsets'][h:h+2])
                if np.any(result['values'][a:b]):
                    raise ValueError('A single-charge scan/reference populated the other charge')
        with np.load(job['batches']) as batches:
            if (batches['titles'].tolist() != layout['titles']
                    or not np.array_equal(batches['edges'], layout['edges'])
                    or not np.array_equal(batches['offsets'], layout['offsets'])
                    or batches['weights'][:82].tolist() != SCALE_WEIGHTS):
                raise ValueError('Retraining batch layout changed')
            mean, variance = collapse_batches(batches['contributions'][:, :, :82], batches['strata'], result['values'])
        key = group_key(origin, scenario, mode, charge, variant, flavour)
        entry = groups.setdefault(key, dict(origin=origin, scenario=scenario, w_treatment=mode, charge=charge,
            variant=variant, flavours=list(flavour), vectors=[], conditional_variances=[], samples=[]))
        entry['vectors'].append(mean)
        entry['conditional_variances'].append(variance[:, 0])
        sample = dict(origin=origin, scenario=scenario, charge=charge, flavours=list(flavour),
            w_treatment=mode, variant=variant, seed=manifest['settings']['iseed'],
            replica=job['case']['replica'] if origin == 'scan' else len(entry['vectors'])-1,
            audit=job['audit'], audit_sha256=digest(job['audit']),
            batches=job['batches'], batches_sha256=digest(job['batches']), stage_initializations=stage_count)
        entry['samples'].append(len(samples))
        samples.append(sample)
        print('Collected independent parameter/reference run', len(samples), '/', len(jobs), job['run'], flush=True)
    arrays, definitions = dict(layout), {}
    for key, entry in sorted(groups.items()):
        if len(entry['vectors']) < 5 or (entry['origin'] == 'scan' and
                [samples[i]['replica'] for i in entry['samples']] != list(range(count))):
            raise ValueError('Missing complete retraining ensemble')
        name = '__'.join(key[:-1]+(''.join(key[-1]),))
        arrays[name] = np.asarray(entry.pop('vectors'))
        arrays[name+'__conditional_variances'] = np.asarray(entry.pop('conditional_variances'))
        definitions[name] = entry
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'), **arrays)
    save(output, dict(created_utc=now(), status=COLLECTED, input=str(queue_path), input_sha256=digest(queue_path),
        benchmark_sha256=queue['benchmark_sha256'], arrays_sha256=digest(output.with_suffix('.npz')),
        release=str(release_path), release_sha256=digest(release_path), groups=definitions, samples=samples,
        comparisons=[dict(scenario=s, w_treatment=m, charge=c, flavours=[list(f) for f in sorted(flavours)])
                     for (s, m, c), flavours in sorted(comparisons.items())],
        reference_stage_initializations=reference_count, scan_stage_initializations=len(streams)-reference_count,
        cross_run_rng_overlap=0, script_sha256=digest(Path(__file__)),
        convention='One row per complete independent retraining. Average within each group, then sum '
                    'exact ordered flavour means at matched coordinates. Reuse shared baseline ensembles '
                    'once for covariance across scenarios. Stored nominal conditional variances are '
                    'diagnostics only and must not be added to the between-retraining variance.',
        limitations='Selected flavour coverage is explicit. Only all eight assignments constitute an '
                    'inclusive direct e/mu prediction. Convergence and resolved sensitivity remain unproven.'))
    print(output, flush=True)


def read(path):
    path = Path(path).resolve()
    record = json.loads(path.read_text())
    if (record['status'] != COLLECTED or digest(record['input']) != record['input_sha256']
            or digest(record['release']) != record['release_sha256']
            or digest(path.with_suffix('.npz')) != record['arrays_sha256']
            or digest(STUDY/'inputs/benchmark.json') != record['benchmark_sha256']):
        raise ValueError('Changed or incomplete independent parameter vectors')
    inventory(json.loads(Path(record['input']).read_text()))
    ensembles, conditional = {}, {}
    with np.load(path.with_suffix('.npz')) as arrays:
        layout = {key: arrays[key].copy() for key in ('edges', 'offsets')}
        layout.update(titles=arrays['titles'].tolist(), weights=arrays['weights'].tolist())
        if layout['weights'] != SCALE_WEIGHTS:
            raise ValueError('Parameter collection has the wrong canonical weight coordinates')
        for name, definition in record['groups'].items():
            key = group_key(definition['origin'], definition['scenario'], definition['w_treatment'],
                            definition['charge'], definition['variant'], definition['flavours'])
            if key in ensembles:
                raise ValueError('Duplicate parameter retraining group')
            ensembles[key] = arrays[name].copy()
            conditional[key] = arrays[name+'__conditional_variances'].copy()
            if (ensembles[key].ndim != 3 or ensembles[key].shape[1:] != (len(layout['edges']), 82)
                    or len(ensembles[key]) < 5 or conditional[key].shape != ensembles[key].shape[:2]
                    or not np.isfinite(ensembles[key]).all() or not np.isfinite(conditional[key]).all()):
                raise ValueError('Invalid full-run parameter vectors or conditional diagnostics')
    return record, layout, ensembles, conditional


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    collect(args.queue, args.output)
