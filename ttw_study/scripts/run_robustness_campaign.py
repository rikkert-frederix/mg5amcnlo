#!/usr/bin/env python3
"""Full-flavour W-treatment and decay-bottom-mass S/Pi companions.

The massless on-shell reference comes from independently retrained main
runs. Four additional states provide the three W-treatment contrasts and
the two required 0/4.8 GeV mass comparisons. Preparation never launches
MG5; execution requires the corresponding scientific validation evidence.
"""
import argparse
import collections
import copy
import fcntl
import hashlib
import itertools
import json
from pathlib import Path

from audit_results import audit
from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card, digest, generate, now, run, save
from check_phase_space_support import run as check_support
from full_flavour_statistics import FLAVOURS
from harvest_splits import harvest
from load_results import assert_matched_physics, assert_same_layout, load
from production_invariance import run as check_production
from read_splits import read
from run_width_mass_pilots import run_name, verified_pairs

PREPARED = 'full-flavour robustness campaign prepared; validation and main-convergence release required'
RELEASED = 'full-flavour robustness integrations released after validation and main S/Pi convergence inspection'
DONE = 'full-flavour robustness integrations finished; matched reductions and convergence required'
VARIANTS = ('S', 'Pi')
STATES = (('onshell', 0.), ('top-bw', 0.), ('all-bw', 0.), ('onshell', 4.8), ('all-bw', 4.8))
EVIDENCE = ('main_convergence_inspection', 'main_reduction', 'w_treatment_validation',
            'inclusive_normalization', 'massive_decay_validation', 'narrow_w_limit', 'small_mass_limit')


def cases(replicas=5, seed_start=300001):
    if type(replicas) is not int or replicas < 5:
        raise ValueError('Require at least five independently retrained full-flavour rounds')
    count = replicas*2*len(FLAVOURS)*len(STATES[1:])*len(VARIANTS)
    if type(seed_start) is not int or not 0 < seed_start <= 900000000-count:
        raise ValueError('Invalid robustness seed range')
    return [dict(replica=replica, charge=charge, flavours=list(flavours), w_treatment=mode,
                 decay_bottom_mass=mass, variant=variant, seed=seed_start+i,
                 production_sampling='w-current' if mode == 'all-bw' else 'flat')
            for i, (replica, charge, flavours, (mode, mass), variant) in enumerate(itertools.product(
                range(replicas), ('plus', 'minus'), FLAVOURS, STATES[1:], VARIANTS))]


def comparisons():
    pairs = (('top_W', STATES[1], STATES[0]), ('associated_W', STATES[2], STATES[1]),
             ('all_W', STATES[2], STATES[0]), ('bottom_mass_onshell', STATES[3], STATES[0]),
             ('bottom_mass_all_bw', STATES[4], STATES[2]))
    def state(row):
        return dict(w_treatment=row[0], decay_bottom_mass=row[1])
    return [dict(name=name, charge=charge, target=state(target), reference=state(reference),
                 flavours=[list(row) for row in FLAVOURS])
            for charge in ('plus', 'minus') for name, target, reference in pairs]


def export_key(case):
    return case['charge'], tuple(case['flavours']), case['w_treatment'], case['decay_bottom_mass']


def process_path(case, tag):
    return STUDY/'processes'/('TTW%s_%s_%s_mb%s_both_%s' % (
        case['charge'], case['w_treatment'], ''.join(case['flavours']),
        str(case['decay_bottom_mass']).replace('.', 'p'), tag))


def source_hashes():
    paths = set(RUNTIME_SOURCES)
    paths.update('ttw_study/scripts/'+name for name in (
        'run_robustness_campaign.py', 'run_width_mass_pilots.py', 'run_sampler_pilots.py',
        'audit_results.py', 'audit_virtuals.py', 'harvest_splits.py', 'read_splits.py',
        'rng_history.py', 'compare_current_batches.py', 'load_results.py', 'full_flavour_statistics.py',
        'replica_statistics.py', 'production_invariance.py', 'inclusive_report.py',
        'bw_branching.py', 'check_phase_space_support.py'))
    paths.update(('madgraph/fks/fks_decay_masses.py', 'models/model_reader.py',
        'Template/fNLO/Source/kin_functions.f90',
        'tests/input_files/fks_decay/phase_space_test_dimensions.f90',
        'tests/input_files/fks_decay/bw_support_checks.f90'))
    return {name: digest(ROOT/name) for name in sorted(paths)}


def plan_digest(record):
    keys = ('tag','replicas','seed_start','cases','comparisons','settings','source_hashes','benchmark_sha256','max_cores')
    return hashlib.sha256(json.dumps({key: record[key] for key in keys}, sort_keys=True,
                                    allow_nan=False).encode()).hexdigest()


def prepare(tag, replicas=5, seed_start=300001, accuracy=.015, grid_points=1000):
    if (not tag.replace('_', '').isalnum() or not 0. < accuracy < 1.
            or type(grid_points) is not int or grid_points < 200):
        raise ValueError('Invalid robustness tag or initial accuracy/grid allocation')
    selected = cases(replicas, seed_start)
    destination = STUDY/'inputs'/('robustness_campaign_'+tag+'.json')
    if destination.exists():
        raise ValueError('Existing robustness campaign; inspect before further work')
    record = dict(created_utc=now(), status=PREPARED, tag=tag, replicas=replicas, seed_start=seed_start,
        cases=selected, jobs=[], exports={}, comparisons=comparisons(), source_hashes=source_hashes(),
        benchmark_sha256=digest(STUDY/'inputs/benchmark.json'), max_cores=64,
        settings=dict(accuracy=accuracy, points=1000, grid_points=grid_points, iterations=3,
                      job_seconds=1., production_scale='core-w-ht-half', ecm=13000.),
        reference_scope='At least five independent main massless on-shell S/Pi runs per ordered flavour '
                        'and charge. References are shared, never counted as new integrations.',
        scope='All eight ordered direct e/mu assignments in both charges, S/Pi in massless top-bw/all-bw '
              'and in onshell/all-bw with 4.8 GeV decay-only bottom mass. Production remains massless 5FS. '
              'Common dynamic W-system scale, matched top widths and all 81 scales plus 101 PDF members.',
        precision='Total-rate accuracy is an initial allocation. Measure full-flavour fiducial/shape '
                  'precision and mixed-term sensitivity using independent whole-run covariance; '
                  'inspect common bins and retraining convergence before publication.',
        remaining='Conditional top-bw massive or mb=4.6/5.0 sensitivity, selected PDF/alpha-s/top-mass/energy '
                  'checks, resolved-shift attribution and final paper conclusions remain separate tasks.')
    record['plan_sha256'] = plan_digest(record)
    save(destination, record)
    return destination


def unchanged(record):
    if record.get('plan_sha256') != plan_digest(record):
        raise ValueError('Frozen robustness allocation or source inventory changed')
    if (digest(STUDY/'inputs/benchmark.json') != record['benchmark_sha256']
            or any(digest(ROOT/name) != expected for name, expected in record['source_hashes'].items())):
        raise ValueError('Frozen robustness source or benchmark changed')
    if (record['cases'] != cases(record['replicas'], record['seed_start'])
            or record['comparisons'] != comparisons()):
        raise ValueError('Robustness cases or comparisons differ from full required coverage')


def release_evidence(path, record):
    release = json.loads(Path(path).read_text())
    if (release.get('status') != RELEASED or release.get('benchmark_sha256') != record['benchmark_sha256']
            or release.get('campaign_tag') != record['tag']
            or release.get('campaign_plan_sha256') != record['plan_sha256']):
        raise ValueError('Require validation and main S/Pi convergence inspection for this campaign')
    for name in EVIDENCE:
        evidence = release.get('evidence', {}).get(name)
        if not evidence or any(digest(p) != sha for p, sha in evidence.items()):
            raise ValueError('Missing or changed robustness release evidence: '+name)
    required = set(itertools.product(('plus', 'minus'), FLAVOURS, VARIANTS))
    coverage, streams, processes = collections.Counter(), set(), {}
    reference, layout = None, None
    jobs = release.get('reference_jobs', [])
    for job in jobs:
        result = load(job['audit'])
        report = result['report']
        manifest, execution = report['manifest'], report['execution']
        generation = execution['export_generation']['args']
        key = generation['charge'], tuple(generation['flavours']), job['variant']
        if (key not in required or execution['pilot'] or execution.get('grid_reference') is not None
                or manifest['variant'] != job['variant'] or manifest['w_treatment'] != 'onshell'
                or manifest['decay_bottom_mass'] != 0. or generation['corrected'] != 'both'
                or generation['decay_bottom_mass'] != 0. or generation['w_treatment'] != 'onshell'
                or manifest['settings']['lhaid'] != 331700
                or manifest['settings']['ebeam1'] != 6500. or manifest['settings']['ebeam2'] != 6500.
                or manifest['top_width_reference_scale'] != 172.5
                or manifest['production_scale'] != 'core-w-ht-half'
                or manifest['decay_scale_grouping'] != 'separate'):
            raise ValueError('Require independently retrained matched main on-shell S/Pi references')
        if (not execution['source_hashes'] or
                any(record['source_hashes'].get(p) != sha for p, sha in execution['source_hashes'].items())):
            raise ValueError('Main reference uses a different frozen runtime source')
        if (str(Path(report['path']).parents[2]) != job['process']
                or Path(report['path']).parent.name != job['run']):
            raise ValueError('Reference job label differs from its actual audited output')
        projected = dict(result, report=copy.deepcopy(report))
        projected['report']['execution']['export_generation']['args'].update(charge='plus', flavours=['e','e','e'])
        if reference is None:
            reference = projected
            layout = {name: result[name] for name in ('titles','edges','offsets','weights')}
        else:
            assert_matched_physics([reference, projected])
            assert_same_layout([layout, result])
        pairs = verified_pairs(job)
        if streams & pairs:
            raise ValueError('Baseline main references reuse actual stage streams')
        streams.update(pairs)
        coverage[key] += 1
        processes.setdefault(key[:2]+('onshell', 0.), Path(job['process']))
    if set(coverage) != required or any(count < 5 for count in coverage.values()):
        raise ValueError('Require five or more independent main S/Pi references for all eight flavours and both charges')
    return release, streams, processes


def execute(path, release_path):
    path, release_path = Path(path).resolve(), Path(release_path).resolve()
    record = json.loads(path.read_text())
    if record['status'] != PREPARED or record['jobs'] or record['exports']:
        raise ValueError('Inspect existing execution; never restart a robustness campaign implicitly')
    unchanged(record)
    _, streams, references = release_evidence(release_path, record)
    reference_count = len(streams)
    record.update(status='running full-flavour robustness integrations; convergence not certified',
                  release=str(release_path), release_sha256=digest(release_path))
    save(path, record)
    try:
        for case in record['cases']:
            unchanged(record)
            key = export_key(case)
            process = process_path(case, record['tag'])
            name = run_name(case)
            if key not in references:
                record.update(current_case=case, current_process=str(process))
                save(path, record)
                generate(argparse.Namespace(charge=case['charge'], flavours=case['flavours'],
                    w_treatment=case['w_treatment'], decay_bottom_mass=case['decay_bottom_mass'],
                    corrected='both', tag=record['tag']))
                support = STUDY/'inputs'/(process.name+'_phase_space_support.json')
                check_support(support, process)
                exported = dict(process=str(process), local_support=str(support), local_support_sha256=digest(support))
                if case['decay_bottom_mass']:
                    baseline = references[key[:-1]+(0.,)]
                    identity = STUDY/'inputs'/(process.name+'_production_identity.json')
                    with (STUDY/'mg5.lock').open('a') as lock:
                        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        benchmark_param_card(process, json.loads((STUDY/'inputs/benchmark.json').read_text()),
                                             case['w_treatment'])
                        check_production(baseline, process, identity)
                    exported.update(production_identity=str(identity), production_identity_sha256=digest(identity))
                references[key] = process
                record['exports'][process.name] = exported
                save(path, record)
            if ((process/'Events'/name).exists() or (process/'study_cards'/name).exists()
                    or (STUDY/'logs'/(process.name+'_'+name+'.log')).exists()):
                raise ValueError('Existing robustness run artifacts; preserve them and inspect')
            unchanged(record)
            record.update(current_case=case, current_process=str(process))
            save(path, record)
            run(argparse.Namespace(process_dir=str(process), variant=case['variant'], seed=case['seed'],
                production_sampling=case['production_sampling'], main_run=True, grid_reference=None, **record['settings']))
            audit_path = STUDY/'results/audits'/(process.name+'_'+name+'.json')
            workers = STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
            with (STUDY/'mg5.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                audit(process/'Events'/name/'MADatNLO.HwU')
                harvest(audit_path, workers)
            batches = STUDY/'results'/(process.name+'_'+name+'_batches.npz')
            read(workers.with_suffix('.json'), batches)
            virtuals = STUDY/'results'/(process.name+'_'+name+'_virtuals.json')
            audit_virtuals(workers.with_suffix('.json'), virtuals)
            job = dict(case=case, variant=case['variant'], process=str(process), run=name,
                audit=str(audit_path), workers=str(workers), batches=str(batches),
                analytic_virtual_checks=str(virtuals), export_checks=record['exports'][process.name], finished_utc=now())
            pairs = verified_pairs(job)
            if pairs & streams:
                raise ValueError('Robustness training/refinement streams overlap main references or previous companions')
            streams.update(pairs)
            record['jobs'].append(job)
            record.update(reference_stage_initializations=reference_count,
                          companion_stage_initializations=len(streams)-reference_count)
            save(path, record)
        record.update(status=DONE, current_case=None, current_process=None, finished_utc=now())
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
    p.add_argument('--seed-start', type=int, default=300001)
    p.add_argument('--accuracy', type=float, default=.015)
    p.add_argument('--grid-points', type=int, default=1000)
    p = sub.add_parser('execute')
    p.add_argument('record', type=Path)
    p.add_argument('--release', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'prepare':
        print(prepare(args.tag, args.replicas, args.seed_start, args.accuracy, args.grid_points))
    else:
        execute(args.record, args.release)
