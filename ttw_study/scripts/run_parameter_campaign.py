#!/usr/bin/env python3
"""Serialized, independently retrained S/Pi parameter-scan integrations.

The selection is explicit and must be released after inspecting main S/Pi
convergence. Preparation is not a physical result or a release decision.
No parameter campaign is selected automatically by this script.
"""
import argparse
import collections
import fcntl
import itertools
import json
from pathlib import Path
import sys

from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, save
from harvest_splits import harvest
from load_results import load as load_baseline
from parameter_variation_cards import configure
from parameter_variation_inputs import SCENARIOS, validate
from parameter_variation_results import audit, read_batches
from run_inclusive import launch
from run_width_mass_pilots import verified_pairs

PREPARED = 'parameter campaign prepared; main S/Pi convergence release required'
DONE = 'parameter integrations finished; matched comparisons and convergence required'
RELEASED = 'parameter integration released after main S/Pi convergence inspection'
VARIANTS = ('S', 'Pi')


def contexts(selection):
    if selection.get('schema') != 'parameter_scan_selection_v1' or not selection.get('reason', '').strip():
        raise ValueError('Require an explicit parameter selection and scientific reason')
    rows = selection['contexts']
    if not rows:
        raise ValueError('Empty parameter selection')
    valid_scenarios = {row[0] for row in SCENARIOS}
    seen = set()
    for row in rows:
        if (set(row) != {'scenario', 'charge', 'flavours', 'w_treatment'}
                or row['scenario'] not in valid_scenarios or row['charge'] not in ('plus', 'minus')
                or len(row['flavours']) != 3 or any(f not in ('e', 'mu') for f in row['flavours'])
                or row['w_treatment'] not in ('onshell', 'top-bw', 'all-bw')):
            raise ValueError('Unsupported parameter-scan context')
        key = row['scenario'], row['charge'], tuple(row['flavours']), row['w_treatment']
        if key in seen:
            raise ValueError('Repeated parameter-scan context')
        seen.add(key)
    return rows


def cases(selection, replicas=5, seed_start=200001):
    rows = contexts(selection)
    if not isinstance(replicas, int) or replicas < 5:
        raise ValueError('Require at least five independent retrainings per selected prediction')
    count = replicas*len(rows)*len(VARIANTS)
    if not isinstance(seed_start, int) or not 0 < seed_start <= 900000000-count:
        raise ValueError('Invalid parameter-campaign seed range')
    return [dict(row, replica=replica, variant=variant, seed=seed_start+i, decay_bottom_mass=0.)
            for i, (replica, row, variant) in enumerate(itertools.product(range(replicas), rows, VARIANTS))]


def prepare(tag, selection_path, inputs_path, replicas=5, seed_start=200001, accuracy=.015):
    if not tag.replace('_', '').isalnum() or not 0. < accuracy < 1.:
        raise ValueError('Invalid parameter-campaign tag or accuracy')
    selection_path, inputs_path = Path(selection_path).resolve(), Path(inputs_path).resolve()
    selection = json.loads(selection_path.read_text())
    validate(json.loads(inputs_path.read_text()))
    selected = cases(selection, replicas, seed_start)
    destination = STUDY/'inputs'/('parameter_campaign_'+tag+'.json')
    if destination.exists():
        raise ValueError('Existing parameter campaign; inspect before further work')
    sources = {name: digest(ROOT/name) for name in RUNTIME_SOURCES}
    for name in ('run_parameter_campaign.py', 'parameter_variation_inputs.py', 'parameter_variation_cards.py',
                 'parameter_variation_results.py', 'stage_batch_reader.py', 'run_inclusive.py',
                 'harvest_splits.py', 'read_splits.py', 'rng_history.py', 'load_results.py',
                 'audit_virtuals.py', 'run_width_mass_pilots.py', 'compare_current_batches.py'):
        sources['ttw_study/scripts/'+name] = digest(STUDY/'scripts'/name)
    record = dict(created_utc=now(), status=PREPARED, tag=tag, cases=selected, jobs=[],
        source_hashes=sources, benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
        inputs=str(inputs_path), inputs_sha256=digest(inputs_path),
        selection=str(selection_path), selection_sha256=digest(selection_path),
        selection_reason=selection['reason'], max_cores=64, accuracy=accuracy, replicas=replicas, seed_start=seed_start,
        selected_scenarios=sorted({row['scenario'] for row in selected}),
        scope='Only explicitly selected S/Pi, charge, ordered flavour and W-treatment contexts. '
              'Member 0 and all 81 signed scales, matched decay coefficients and fresh grid retraining.',
        precision='Total-rate accuracy is an initial allocation. Reassess fiducial and shape errors, '
                  'whole-run covariance and matched Pi/S shifts before claiming parameter sensitivity.',
        remaining='The full study requires all seven bounded variations and actual direct e/mu flavour '
                  'sums for its inclusive predictions. A completed selected queue alone does not establish that scope.')
    save(destination, record)
    return destination


def release_evidence(path, record):
    release = json.loads(Path(path).read_text())
    if (release.get('status') != RELEASED or release.get('benchmark_sha256') != record['benchmark_sha256']
            or release.get('selection_sha256') != record['selection_sha256']):
        raise ValueError('Require inspection of main S/Pi convergence for this exact selection')
    for name in ('main_convergence_inspection', 'main_reduction'):
        evidence = release.get('evidence', {}).get(name)
        if not evidence or any(digest(p) != sha for p, sha in evidence.items()):
            raise ValueError('Missing or changed scientific release evidence: '+name)
    required = {(row['charge'], tuple(row['flavours']), row['w_treatment'], variant)
                for row in record['cases'] for variant in VARIANTS}
    coverage, streams = collections.Counter(), set()
    jobs = release.get('reference_jobs', [])
    if not jobs:
        raise ValueError('Missing audited baseline reference integrations')
    for job in jobs:
        original = load_baseline(job['audit'])['report']
        manifest = original['manifest']
        generated = original['execution']['export_generation']['args']
        key = generated['charge'], tuple(generated['flavours']), manifest['w_treatment'], job['variant']
        if (key not in required or manifest['variant'] != job['variant']
                or manifest['settings']['lhaid'] != 331700
                or manifest['settings']['ebeam1'] != 6500. or manifest['settings']['ebeam2'] != 6500.
                or manifest['decay_bottom_mass'] != 0. or manifest['top_width_reference_scale'] != 172.5
                or manifest['production_scale'] != 'core-w-ht-half' or generated['corrected'] != 'both'):
            raise ValueError('Unmatched main-benchmark reference context')
        pairs = verified_pairs(job)
        if streams & pairs:
            raise ValueError('Reference integrations reuse actual stage streams')
        streams.update(pairs)
        coverage[key] += 1
    if set(coverage) != required or any(n < 5 for n in coverage.values()):
        raise ValueError('Require at least five independent audited baseline S/Pi retrainings for each context')
    return release, streams


def run_case(process, case, record):
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        name, archive = configure(process, record['inputs'], case['scenario'], case['variant'],
                                  case['seed'], case['w_treatment'], record['accuracy'])
        execution = dict(max_cores=64, pilot=False, parameter_case=case,
            source_hashes=record['source_hashes'], grid_reference=None,
            requested_total_rate_accuracy=record['accuracy'],
            export_generation=json.loads((STUDY/'inputs'/(process.name+'_generation.json')).read_text()))
        launch(['taskset', '-c', '0-63', sys.executable, str(process/'bin/calculate_xsect'),
                'NLO', '-f', '-n', name, '--multicore', '--nb_core=64'], process,
               STUDY/'logs'/(process.name+'_'+name+'.log'), archive/'execution.json', execution)
        hwu = process/'Events'/name/'MADatNLO.HwU'
        if not hwu.is_file():
            execution.update(status='failed: no final histogram', missing_output=str(hwu))
            save(archive/'execution.json', execution)
            raise ValueError('Missing final parameter-scan histogram')
        execution['outputs'] = {str(hwu): digest(hwu)}
        save(archive/'execution.json', execution)
        audit_path = STUDY/'results/audits'/(process.name+'_'+name+'.json')
        audit(hwu, audit_path)
        workers = STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
        harvest(audit_path, workers)
    batches = STUDY/'results'/(process.name+'_'+name+'_batches.npz')
    read_batches(workers.with_suffix('.json'), batches)
    virtuals = STUDY/'results'/(process.name+'_'+name+'_virtuals.json')
    audit_virtuals(workers.with_suffix('.json'), virtuals)
    return dict(case=case, variant=case['variant'], process=str(process), run=name,
        audit=str(audit_path), workers=str(workers), batches=str(batches),
        analytic_virtual_checks=str(virtuals), finished_utc=now())


def execute(path, release_path):
    path, release_path = Path(path).resolve(), Path(release_path).resolve()
    record = json.loads(path.read_text())
    if record['status'] != PREPARED or record['jobs']:
        raise ValueError('Inspect existing execution; never restart a parameter campaign implicitly')
    def unchanged():
        if (digest(record['inputs']) != record['inputs_sha256']
                or digest(record['selection']) != record['selection_sha256']
                or digest(STUDY/'inputs/benchmark.json') != record['benchmark_sha256']
                or any(digest(ROOT/p) != sha for p, sha in record['source_hashes'].items())):
            raise ValueError('Frozen parameter-campaign inputs or source changed')
        expected = cases(json.loads(Path(record['selection']).read_text()), record['replicas'], record['seed_start'])
        if record['cases'] != expected:
            raise ValueError('Parameter cases differ from the frozen explicit selection')
    unchanged()
    _, streams = release_evidence(release_path, record)
    reference_count = len(streams)
    record.update(status='running selected parameter integrations; convergence not certified',
        release=str(release_path), release_sha256=digest(release_path),
        reference_stage_pairs=reference_count)
    save(path, record)
    processes = set()
    try:
        for case in record['cases']:
            unchanged()
            export_tag = record['tag']+'_'+case['scenario']
            process = STUDY/'processes'/('TTW%s_%s_%s_mb0p0_both_%s' % (
                case['charge'], case['w_treatment'], ''.join(case['flavours']), export_tag))
            record['current_case'] = case
            save(path, record)
            if process not in processes:
                generate(argparse.Namespace(charge=case['charge'], flavours=case['flavours'],
                    w_treatment=case['w_treatment'], decay_bottom_mass=0., corrected='both', tag=export_tag))
                processes.add(process)
            job = run_case(process, case, record)
            pairs = verified_pairs(job)
            if streams & pairs:
                raise ValueError('Scan streams overlap another scan or its baseline reference')
            streams.update(pairs)
            record['jobs'].append(job)
            record['distinct_scan_stage_pairs'] = len(streams)-reference_count
            save(path, record)
        record.update(status=DONE, current_case=None, finished_utc=now())
    except BaseException as error:
        record.update(status='stopped', error=repr(error), stopped_utc=now())
        save(path, record)
        raise
    save(path, record)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--tag', required=True)
    p.add_argument('--selection', type=Path, required=True)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--replicas', type=int, default=5)
    p.add_argument('--seed-start', type=int, default=200001)
    p.add_argument('--accuracy', type=float, default=.015)
    p = sub.add_parser('execute')
    p.add_argument('record', type=Path)
    p.add_argument('--release', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare':
        print(prepare(args.tag, args.selection, args.inputs, args.replicas, args.seed_start, args.accuracy))
    else:
        execute(args.record, args.release)
