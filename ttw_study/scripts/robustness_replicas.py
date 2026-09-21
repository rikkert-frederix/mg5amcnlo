#!/usr/bin/env python3
"""Collect complete W/mass retrainings and shared main on-shell references."""
import argparse
import copy
import itertools
import json
from pathlib import Path

import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from full_flavour_statistics import FLAVOURS
from load_results import WEIGHTS, assert_matched_physics, assert_same_layout, load
from main_replicas import collapse_batches
from mass_effects import validate_physics as validate_mass
from robustness_statistics import group_key, select
from run_robustness_campaign import DONE, STATES, cases, comparisons, release_evidence, unchanged
from run_width_mass_pilots import verified_pairs
from w_treatment_pilot_report import validate_physics as validate_w

COLLECTED = 'complete independent robustness vectors collected; comparisons and convergence required'


def inventory(queue):
    if queue.get('status') != DONE:
        raise ValueError('Require a complete full-flavour robustness campaign')
    unchanged(queue)
    expected = cases(queue['replicas'], queue['seed_start'])
    if queue['cases'] != expected or [job['case'] for job in queue['jobs']] != expected:
        raise ValueError('Missing, repeated or reordered full-flavour robustness retrainings')
    return queue['replicas']


def generated_evidence(result):
    report = result['report']
    manifest = report['manifest']
    process = Path(report['path']).parents[2]
    for kind in ('topology_hashes', 'internal_width_hashes'):
        if not manifest[kind] or any(digest(process/path) != sha for path, sha in manifest[kind].items()):
            raise ValueError('Generated W topology or internal-width evidence changed')
    for relative in manifest['internal_width_hashes']:
        if json.loads((process/relative).read_text())['pdgs'] != ([] if manifest['w_treatment'] == 'onshell' else [24]):
            raise ValueError('Wrong internal W propagator-width prescription')
    for kind in ('scale_runtime_hashes', 'sampling_runtime_hashes'):
        if not manifest[kind] or any(digest(process/'SubProcesses'/path) != sha for path, sha in manifest[kind].items()):
            raise ValueError('Generated scale or sampling source changed')


def validate_comparisons(representatives, identities, full_flavour=True):
    """Validate all physical W/mass edges without changing loaded evidence.

    Representatives carry only metadata/layout; subsequent retrainings are
    required to match them exactly. Existing numeric-card, width, topology,
    scale and production/FKS validators apply to the actual archived files.
    """
    contexts = sorted({(key[2], key[4]) for key in representatives})
    if full_flavour and set(contexts) != set(itertools.product(('plus', 'minus'), FLAVOURS)):
        raise ValueError('Physical validation requires both charges and all eight flavours')
    expected = {group_key(mode, mass, charge, variant, flavour)
                for charge, flavour in contexts for mode, mass in STATES for variant in ('S', 'Pi')}
    if not expected or set(representatives) != expected:
        raise ValueError('Missing physical W/mass state or S/Pi representative')
    for charge, flavour in contexts:
        def rows(states):
            return [representatives[group_key(mode, mass, charge, variant, flavour)]
                    for mode, mass in states for variant in ('S', 'Pi')]
        validate_w(rows(STATES[:3]))
        for mode in ('onshell', 'all-bw'):
            chosen = rows(((mode, 0.), (mode, 4.8)))
            processes = {str(Path(row['report']['path']).parents[2]) for row in chosen[2:]}
            if not processes <= set(identities):
                raise ValueError('Missing massive production-source identity')
            validate_mass(chosen, [Path(identities[process]) for process in sorted(processes)])
    return len(contexts)


def collect(queue_path, output):
    queue_path, output = Path(queue_path).resolve(), Path(output).resolve()
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite robustness retraining vectors')
    queue = json.loads(queue_path.read_text())
    count = inventory(queue)
    release_path = Path(queue['release'])
    if digest(release_path) != queue['release_sha256']:
        raise ValueError('Robustness scientific release changed')
    release, streams, _ = release_evidence(release_path, queue)
    reference_count = len(streams)
    jobs = [('reference', job) for job in release['reference_jobs']]+[('companion', job) for job in queue['jobs']]
    groups, samples, representatives, identities, physical_states = {}, [], {}, {}, {}
    layout = None
    for origin, job in jobs:
        result = load(job['audit'])
        report = result['report']
        manifest, execution = report['manifest'], report['execution']
        generated = execution['export_generation']['args']
        mode, mass, charge, variant, flavour = (manifest['w_treatment'], manifest['decay_bottom_mass'],
            generated['charge'], job['variant'], tuple(generated['flavours']))
        if ((mode, mass) not in STATES or charge not in ('plus', 'minus') or flavour not in FLAVOURS
                or variant not in ('S', 'Pi') or report['variant'] != variant or manifest['variant'] != variant
                or generated['w_treatment'] != mode or generated['decay_bottom_mass'] != mass
                or generated['corrected'] != 'both' or execution['pilot'] or execution.get('grid_reference') is not None
                or manifest['production_scale'] != 'core-w-ht-half' or manifest['decay_scale_grouping'] != 'separate'
                or not execution['source_hashes']
                or any(queue['source_hashes'].get(p) != sha for p, sha in execution['source_hashes'].items())
                or str(Path(report['path']).parents[2]) != job['process']
                or Path(report['path']).parent.name != job['run']):
            raise ValueError('Robustness run differs from its frozen process or independent-retraining convention')
        generated_evidence(result)
        if origin == 'companion':
            case = job['case']
            expected = dict(charge=case['charge'], flavours=case['flavours'], w_treatment=case['w_treatment'],
                decay_bottom_mass=case['decay_bottom_mass'], corrected='both', tag=queue['tag'])
            if (any(generated[key] != value for key, value in expected.items()) or variant != case['variant']
                    or manifest['settings']['iseed'] != case['seed']
                    or manifest['production_sampling'] != case['production_sampling']):
                raise ValueError('Companion differs from its allocated case')
            exported = queue['exports'][Path(job['process']).name]
            if exported != job['export_checks'] or exported['process'] != job['process']:
                raise ValueError('Companion export proof differs from the campaign inventory')
            for kind in ('local_support', 'production_identity'):
                if kind in exported and digest(exported[kind]) != exported[kind+'_sha256']:
                    raise ValueError('Companion phase-space or production identity evidence changed')
            if mass:
                identities[job['process']] = exported['production_identity']
            pairs = verified_pairs(job)
            if streams & pairs:
                raise ValueError('Robustness/reference retrainings share actual stage streams')
            streams.update(pairs)
            stage_count = len(pairs)
        else:
            if (mode, mass) != STATES[0]:
                raise ValueError('Main reference is not massless on-shell')
            # The release validator already checked every reference stream.
            metadata = json.loads(Path(job['batches']).with_suffix('.json').read_text())
            stage_count = metadata['refinement_rng_audit']['distinct_initialization_pairs']
        if result['weights'] != WEIGHTS:
            raise ValueError('Require all 81 canonical scale points and 101 PDF members')
        if layout is None:
            layout = {key: result[key] for key in ('titles', 'edges', 'offsets', 'weights')}
        else:
            assert_same_layout([layout, result])
        light = {key: result[key] for key in ('report', 'parameter_signature', 'titles', 'edges', 'offsets', 'weights')}
        projected = dict(light, report=copy.deepcopy(report))
        projected['report']['execution']['export_generation']['args'].update(charge='plus', flavours=['e', 'e', 'e'])
        state = mode, mass
        if state in physical_states:
            assert_matched_physics([physical_states[state], projected])
        else:
            physical_states[state] = projected
        for h, title in enumerate(result['titles']):
            if ('W-' if charge == 'plus' else 'W+') in title:
                a, b = map(int, result['offsets'][h:h+2])
                if np.any(result['values'][a:b]):
                    raise ValueError('A single-charge robustness sample populated the other charge')
        with np.load(job['batches']) as batches:
            if (batches['titles'].tolist() != layout['titles'] or batches['weights'].tolist() != WEIGHTS
                    or not np.array_equal(batches['edges'], layout['edges'])
                    or not np.array_equal(batches['offsets'], layout['offsets'])):
                raise ValueError('Robustness batch layout differs from the final audit')
            mean, variance = collapse_batches(batches['contributions'], batches['strata'], result['values'])
        key = group_key(mode, mass, charge, variant, flavour)
        representatives.setdefault(key, light)
        entry = groups.setdefault(key, dict(w_treatment=mode, decay_bottom_mass=mass, charge=charge, variant=variant,
            flavours=list(flavour), origin=origin, vectors=[], conditional_variances=[], samples=[]))
        sample = dict(origin=origin, w_treatment=mode, decay_bottom_mass=mass, charge=charge, variant=variant,
            flavours=list(flavour), seed=manifest['settings']['iseed'],
            replica=job['case']['replica'] if origin == 'companion' else len(entry['vectors']),
            audit=job['audit'], audit_sha256=digest(job['audit']), batches=job['batches'],
            batches_sha256=digest(job['batches']), stage_initializations=stage_count)
        entry['vectors'].append(mean)
        entry['conditional_variances'].append(variance[:, 0].copy())
        entry['samples'].append(len(samples))
        samples.append(sample)
        print('Collected independent robustness/reference run', len(samples), '/', len(jobs), job['run'], flush=True)
    validated_contexts = validate_comparisons(representatives, identities)
    arrays, definitions = dict(layout), {}
    for i, (key, entry) in enumerate(sorted(groups.items())):
        if len(entry['vectors']) < 5 or (entry['origin'] == 'companion' and
                [samples[j]['replica'] for j in entry['samples']] != list(range(count))):
            raise ValueError('Missing complete independent robustness retraining ensemble')
        name = 'group_%03d' % i
        arrays[name] = np.asarray(entry.pop('vectors'))
        arrays[name+'__conditional_variances'] = np.asarray(entry.pop('conditional_variances'))
        definitions[name] = entry
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'), **arrays)
    sources = ('robustness_replicas.py', 'robustness_statistics.py', 'main_replicas.py', 'mass_effects.py',
               'w_treatment_pilot_report.py', 'load_results.py', 'replica_statistics.py')
    save(output, dict(created_utc=now(), status=COLLECTED, input=str(queue_path), input_sha256=digest(queue_path),
        benchmark_sha256=queue['benchmark_sha256'], arrays_sha256=digest(output.with_suffix('.npz')),
        release=str(release_path), release_sha256=digest(release_path), groups=definitions, samples=samples,
        comparisons=comparisons(), reference_stage_initializations=reference_count,
        companion_stage_initializations=len(streams)-reference_count, cross_run_rng_overlap=0,
        validated_charge_flavour_contexts=validated_contexts,
        physics_identities={path: digest(path) for path in sorted(set(identities.values()))},
        reduction_source_hashes={str(STUDY/'scripts'/name): digest(STUDY/'scripts'/name) for name in sources},
        convention='One row per complete independently retrained integration. Average each ordered flavour, '
                   'then sum all eight at matched scale/PDF coordinates. Shared physical states occur once. '
                   'Conditional nominal variances diagnose training stability; never add them to retraining variance.',
        limitations='A collected complete campaign is not a convergence certificate. Common bins, rates/shapes, '
                    'charge observables, precision and sensitivity must still be reduced and inspected.'))
    print(output, flush=True)


def read(path):
    path = Path(path).resolve()
    record = json.loads(path.read_text())
    if (record['status'] != COLLECTED or digest(record['input']) != record['input_sha256']
            or digest(record['release']) != record['release_sha256']
            or digest(path.with_suffix('.npz')) != record['arrays_sha256']
            or digest(STUDY/'inputs/benchmark.json') != record['benchmark_sha256']
            or record['comparisons'] != comparisons()):
        raise ValueError('Changed or incomplete independent robustness archive')
    inventory(json.loads(Path(record['input']).read_text()))
    ensembles, conditional = {}, {}
    with np.load(path.with_suffix('.npz')) as arrays:
        layout = {key: arrays[key].copy() for key in ('edges', 'offsets')}
        layout.update(titles=arrays['titles'].tolist(), weights=arrays['weights'].tolist())
        if layout['weights'] != WEIGHTS:
            raise ValueError('Wrong canonical scale/PDF coordinates in robustness archive')
        for name, definition in record['groups'].items():
            key = group_key(*(definition[field] for field in
                ('w_treatment', 'decay_bottom_mass', 'charge', 'variant', 'flavours')))
            if key in ensembles:
                raise ValueError('Duplicate robustness state/flavour group')
            value, variance = arrays[name].copy(), arrays[name+'__conditional_variances'].copy()
            if (value.ndim != 3 or len(value) < 5 or value.shape[1:] != (len(layout['edges']), len(WEIGHTS))
                    or variance.shape != value.shape[:2] or not np.isfinite(value).all()
                    or not np.isfinite(variance).all() or np.any(variance < 0.)):
                raise ValueError('Invalid full-run robustness vectors or conditional diagnostics')
            ensembles[key], conditional[key] = value, variance
    if set(select(ensembles, record['comparisons'])) != set(ensembles):
        raise ValueError('Unreferenced physical state in robustness archive')
    return record, layout, ensembles, conditional


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    collect(args.queue, args.output)
