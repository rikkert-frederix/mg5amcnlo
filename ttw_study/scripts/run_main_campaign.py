#!/usr/bin/env python3
"""Full direct-flavour on-shell campaign with independent grid retraining.

Preparation writes an immutable case/source record. Execution requires a
separate evidence record releasing the main integrations after inspection
of the technical/reference results. This is an internal scientific gate,
not a request for additional user authorization.
"""
import argparse
import fcntl
import itertools
import json
from pathlib import Path

from audit_results import audit
from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, run, save
from harvest_splits import harvest
from read_splits import read
from run_width_mass_pilots import verified_pairs

VARIANTS = ('S', 'Pi', 'LO', 'P', 'D', 'PiD')
FLAVOURS = tuple(itertools.product(('e', 'mu'), repeat=3))


def cases(replicas=5, seed_start=100001):
    if not isinstance(replicas, int) or replicas < 5:
        raise ValueError('At least five independent retrainings per flavour/prescription are required')
    count = replicas*2*8*6
    if not isinstance(seed_start, int) or not 0 < seed_start <= 900000000-count:
        raise ValueError('Invalid main seed range')
    # Each completed round contains every physical flavour once for all six
    # prescriptions and both charges. Repeated rounds are averaged, not summed.
    return [dict(replica=replica, charge=charge, flavours=list(flavours), variant=variant,
                 w_treatment='onshell', decay_bottom_mass=0., seed=seed_start+i)
            for i, (replica, charge, flavours, variant) in enumerate(itertools.product(
                range(replicas), ('plus', 'minus'), FLAVOURS, VARIANTS))]


def prepare(tag, replicas, seed_start, accuracy, grid_points):
    if not tag.replace('_', '').isalnum() or not 0 < accuracy < 1 or grid_points < 200:
        raise ValueError('Invalid main campaign settings')
    path = STUDY/'inputs'/('main_campaign_'+tag+'.json')
    if path.exists():
        raise ValueError('Existing campaign record; inspect before further work')
    sources = {p: digest(ROOT/p) for p in RUNTIME_SOURCES}
    for name in ('run_main_campaign.py', 'audit_results.py', 'audit_virtuals.py',
                 'harvest_splits.py', 'read_splits.py', 'rng_history.py', 'run_width_mass_pilots.py'):
        sources['ttw_study/scripts/'+name] = digest(STUDY/'scripts'/name)
    record = dict(created_utc=now(), status='prepared; technical/reference release required',
        tag=tag, cases=cases(replicas, seed_start), jobs=[], source_hashes=sources,
        benchmark_sha256=digest(STUDY/'inputs/benchmark.json'), max_cores=64,
        settings=dict(accuracy=accuracy, points=1000, grid_points=grid_points, iterations=3,
                      job_seconds=1., production_scale='core-w-ht-half', ecm=13000.),
        scope='All eight ordered direct e/mu assignments, both charges, six on-shell prescriptions. '
              'Five or more independently retrained runs per prediction; no flavour rescaling or grid imports.',
        precision='Accuracy is an initial allocation on the total rate, not fiducial/shape certification. '
                  'Measure between-retraining variance and freeze common bins after full flavour coverage. '
                  'Require sub-percent fiducial errors and 1--2% retained populated shapes, or explicit sensitivity.',
        remaining_scope='BW/massive companions, parameter/PDF-family/alpha-s/mass/13.6-TeV scans, '
                        'coefficient attribution where resolved, final figures and paper remain required.')
    save(path, record)
    print(path)


def release_evidence(path):
    report = json.loads(path.read_text())
    if report['status'] != 'main integration released after technical and reference inspection':
        raise ValueError('Main integrations have not been released by the scientific audit')
    if report['benchmark_sha256'] != digest(STUDY/'inputs/benchmark.json'):
        raise ValueError('Release refers to a different physical benchmark')
    for name in ('ttw_reference', 'ttbar_reference', 'direct_scales', 'inclusive_normalization'):
        evidence = report['evidence'][name]
        if not evidence or any(digest(p) != sha for p, sha in evidence.items()):
            raise ValueError('Missing or changed release evidence: '+name)
    return report


def execute(path, release):
    record = json.loads(path.read_text())
    if record['status'] != 'prepared; technical/reference release required' or record['jobs']:
        raise ValueError('Inspect existing execution; this entry point never restarts an existing run')
    release_evidence(release)
    def unchanged():
        if (any(digest(ROOT/p) != sha for p, sha in record['source_hashes'].items())
                or digest(STUDY/'inputs/benchmark.json') != record['benchmark_sha256']):
            raise ValueError('Frozen main source or physical benchmark changed')
    unchanged()
    record.update(status='running main integrations; convergence not certified',
                  release=str(release.resolve()), release_sha256=digest(release))
    save(path, record)
    processes, streams = {}, set()
    try:
        for case in record['cases']:
            unchanged()
            key = case['charge'], tuple(case['flavours'])
            process = STUDY/'processes'/('TTW%s_onshell_%s_mb0p0_both_%s' % (
                case['charge'], ''.join(case['flavours']), record['tag']))
            if key not in processes:
                generate(argparse.Namespace(charge=case['charge'], flavours=case['flavours'],
                    w_treatment='onshell', decay_bottom_mass=0., corrected='both', tag=record['tag']))
                processes[key] = process
            name = '%s_onshell_core-w-ht-half_separate_%d' % (case['variant'], case['seed'])
            record['current_case'] = case
            save(path, record)
            run(argparse.Namespace(process_dir=str(process), variant=case['variant'],
                seed=case['seed'], main_run=True, grid_reference=None, **record['settings']))
            hwu = process/'Events'/name/'MADatNLO.HwU'
            audit_path = STUDY/'results/audits'/(process.name+'_'+name+'.json')
            workers = STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
            with (STUDY/'mg5.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                audit(hwu)
                harvest(audit_path, workers)
            batches = STUDY/'results'/(process.name+'_'+name+'_batches.npz')
            read(workers.with_suffix('.json'), batches)
            job = dict(case=case, variant=case['variant'], process=str(process), run=name,
                       audit=str(audit_path), workers=str(workers), batches=str(batches))
            pairs = verified_pairs(job)
            if pairs & streams:
                raise ValueError('Actual main-stage random streams overlap earlier main runs')
            streams.update(pairs)
            if case['variant'] in ('D', 'S', 'PiD', 'Pi'):
                virtuals = STUDY/'results'/(process.name+'_'+name+'_virtuals.json')
                audit_virtuals(workers.with_suffix('.json'), virtuals)
                job['analytic_virtual_checks'] = str(virtuals)
            job['finished_utc'] = now()
            record['jobs'].append(job)
            record['distinct_main_stage_pairs'] = len(streams)
            save(path, record)
        record.update(status='main integrations finished; flavour sums and convergence analysis required',
                      current_case=None, finished_utc=now())
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
    p.add_argument('--replicas', type=int, default=5)
    p.add_argument('--seed-start', type=int, default=100001)
    p.add_argument('--accuracy', type=float, default=.015)
    p.add_argument('--grid-points', type=int, default=1000)
    p = sub.add_parser('execute')
    p.add_argument('record', type=Path)
    p.add_argument('--release', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare(args.tag, args.replicas, args.seed_start, args.accuracy, args.grid_points)
    else:
        execute(args.record, args.release)
