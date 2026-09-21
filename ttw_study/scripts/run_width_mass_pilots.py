#!/usr/bin/env python3
"""Fresh, serialized width/mass S/Pi pilots after the RNG-audited controls.

One explicit e,e,mu assignment per charge. These are technical companions,
not the required full-flavour main result or a certification of precision.
"""
import argparse
import fcntl
import json
from pathlib import Path
import time

from audit_results import audit
from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card, digest, generate, now, run, save
from card_writer_transition import validate as validate_writer_transition
from check_decay_scale_cancellation import run as check_decay_scales
from check_phase_space_support import run as check_support
from compare_current_batches import all_stage_pairs, batch_record
from harvest_splits import harvest
from joint_comparison import run as compare_variants
from production_invariance import run as check_production
from read_splits import read
from run_rng_pilots import DONE as PREVIOUS_DONE
from run_sampler_pilots import run_name as massless_run_name, run_settings

DONE='width/mass technical pilots finished; main statistics and continuity remain required'


def cases(seed_start=60001):
    if not isinstance(seed_start,int) or not 0 < seed_start <= 900000000-16:
        raise ValueError('Invalid companion seed range')
    exports=(('plus','onshell',4.8),('plus','all-bw',4.8),('plus','top-bw',0.),
             ('minus','onshell',0.),('minus','top-bw',0.),('minus','all-bw',0.),
             ('minus','onshell',4.8),('minus','all-bw',4.8))
    return [dict(charge=charge,w_treatment=mode,decay_bottom_mass=mass,variant=variant,
                 production_sampling='w-current' if mode=='all-bw' else 'flat',seed=seed_start+2*i+j)
            for i,(charge,mode,mass) in enumerate(exports) for j,variant in enumerate(('S','Pi'))]


def export_key(case):
    return '%s_%s_mb%s' % (case['charge'],case['w_treatment'],str(case['decay_bottom_mass']).replace('.','p'))


def process_path(case,tag):
    return STUDY/'processes'/('TTW%s_%s_eemu_mb%s_both_%s' % (
        case['charge'],case['w_treatment'],str(case['decay_bottom_mass']).replace('.','p'),tag))


def run_name(case):
    name=massless_run_name(case)
    return name+('_mb%s' % ('%.12g' % case['decay_bottom_mass']).replace('.','p')
                 if case['decay_bottom_mass'] else '')


def verified_pairs(job):
    batches=Path(job['batches'])
    metadata,record=batch_record(batches)
    if (digest(batches)!=metadata['arrays_sha256'] or digest(metadata['input'])!=metadata['input_sha256']
            or digest(record['audit'])!=record['audit_sha256']
            or not metadata['status'].startswith('equal-count conditional')):
        raise ValueError('Changed or unsuccessful pilot batch evidence')
    pairs=all_stage_pairs(Path(job['process']),record)
    if len(pairs)!=metadata['refinement_rng_audit']['distinct_initialization_pairs']:
        raise ValueError('Stage stream inventory differs from audited batch evidence')
    return pairs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',type=Path,required=True)
    parser.add_argument('--accuracy',type=float,default=.03)
    parser.add_argument('--seed-start',type=int,default=60001)
    parser.add_argument('--writer-transition',type=Path)
    args=parser.parse_args()
    selected_cases=cases(args.seed_start)
    if not args.tag.replace('_','').isalnum() or not 0.<args.accuracy<1.:
        raise ValueError('Invalid width/mass pilot tag/accuracy')
    destination=STUDY/'inputs'/('width_mass_pilot_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing width/mass queue; do not overwrite it')
    extra=[Path(__file__).resolve(),STUDY/'inputs/benchmark.json',
           ROOT/'madgraph/fks/fks_decay_masses.py',ROOT/'madgraph/iolibs/export_fks.py',
           ROOT/'models/model_reader.py',
           *[STUDY/'scripts'/name for name in (
               'run_rng_pilots.py','run_sampler_pilots.py','production_invariance.py','card_writer_transition.py',
               'audit_results.py','audit_virtuals.py','load_results.py','harvest_splits.py',
               'read_splits.py','rng_history.py','check_phase_space_support.py',
               'check_decay_scale_cancellation.py','inclusive_report.py','bw_branching.py',
               'joint_comparison.py','joint_report.py','replica_statistics.py',
               'pilot_report.py','compare_current_batches.py')],
           STUDY/'inputs/bw_branching_normalization.json',
           *[ROOT/'tests/input_files/fks_decay'/name for name in (
               'phase_space_test_dimensions.f90','bw_support_checks.f90','current_proposal_checks.f90')],
           ROOT/'Template/fNLO/Source/kin_functions.f90']
    if args.writer_transition:
        extra.append(args.writer_transition.resolve())
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in extra})
    record=dict(created_utc=now(),status='waiting for stage-audited controls',source_hashes=sources,
                after_queue=str(args.after_queue.resolve()),max_cores=64,pilot=True,
                flavours=['e','e','mu'],cases=selected_cases,jobs=[],exports={},comparisons={},
                limitations='One ordered flavour assignment. No main precision, small-mass or narrow-W '
                            'limit, full-flavour sum, literature match or final physics conclusion.')
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen width/mass pilot source changed; inspect before further launches')
    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                raise ValueError('Stage-audited predecessor stopped')
            if previous['status']==PREVIOUS_DONE:
                break
            time.sleep(15)
        unchanged()
        if len(previous['jobs'])!=4:
            raise ValueError('Incomplete predecessor controls')
        if args.writer_transition:
            transition=json.loads(args.writer_transition.read_text())
            validate_writer_transition(previous,transition)
            record['writer_transition']=dict(path=str(args.writer_transition.resolve()),
                                             sha256=digest(args.writer_transition),source_transition=transition['source_transition'])
        elif any(digest(ROOT/p)!=sha for p,sha in previous['source_hashes'].items()):
            raise ValueError('Changed predecessor source requires the specific audited writer transition')
        record['predecessor_sha256']=digest(args.after_queue)
        seen=set()
        references={}
        for job in previous['jobs']:
            pairs=verified_pairs(job)
            if seen & pairs:
                raise ValueError('Cross-control training/refinement stream collision')
            seen.update(pairs)
            references[('plus',job['case']['w_treatment'])]=Path(job['process'])
        for case in selected_cases:
            unchanged()
            key=export_key(case)
            process=process_path(case,args.tag)
            name=run_name(case)
            record.update(status='preparing fresh width/mass pilot',
                          current_job=dict(process=str(process),run=name,case=case))
            save(destination,record)
            if key not in record['exports']:
                generate(argparse.Namespace(charge=case['charge'],w_treatment=case['w_treatment'],
                    flavours=['e','e','mu'],corrected='both',decay_bottom_mass=case['decay_bottom_mass'],tag=args.tag))
                support=STUDY/'inputs'/('%s_%s_phase_space_support.json' % (args.tag,key))
                check_support(support,process)
                exported=dict(process=str(process),local_support=str(support))
                if case['decay_bottom_mass']:
                    reference=references[(case['charge'],case['w_treatment'])]
                    identity=STUDY/'inputs'/('%s_%s_production_identity.json' % (args.tag,key))
                    with (STUDY/'mg5.lock').open('a') as lock:
                        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
                        benchmark_param_card(process,json.loads((STUDY/'inputs/benchmark.json').read_text()),case['w_treatment'])
                        check_production(reference,process,identity)
                    exported['production_identity']=str(identity)
                else:
                    references[(case['charge'],case['w_treatment'])]=process
                record['exports'][key]=exported
            unchanged()
            record['status']='running fresh width/mass pilot'
            save(destination,record)
            run(run_settings(process,case,args.accuracy))
            with (STUDY/'mg5.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
                audit(process/'Events'/name/'MADatNLO.HwU')
                audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
                workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
                harvest(audit_path,workers)
            virtuals=STUDY/'results'/(args.tag+'_'+name+'_virtual_checks.json')
            audit_virtuals(workers.with_suffix('.json'),virtuals)
            batches=STUDY/'results'/(args.tag+'_'+name+'_batches.npz')
            read(workers.with_suffix('.json'),batches)
            job=dict(case=case,variant=case['variant'],process=str(process),run=name,audit=str(audit_path),
                     workers=str(workers),batches=str(batches),analytic_virtual_checks=str(virtuals),finished_utc=now())
            pairs=verified_pairs(job)
            if seen & pairs:
                raise ValueError('Cross-pilot training/refinement stream collision')
            seen.update(pairs)
            job['cross_pilot_rng_check']='No initialization pair overlap with preceding controls/pilots'
            decay_scales=STUDY/'results'/(args.tag+'_'+name+'_inclusive_decay_scales.json')
            check_decay_scales(batches,decay_scales)
            job['inclusive_decay_scales']=str(decay_scales)
            record['jobs'].append(job)
            record['distinct_cross_pilot_initialization_pairs']=len(seen)
            save(destination,record)
            if case['variant']=='Pi':
                strict=next(j for j in record['jobs'] if j['process']==str(process) and j['variant']=='S')
                comparison=STUDY/'results'/('%s_%s_joint_comparison.json' % (args.tag,key))
                compare_variants(Path(strict['batches']),batches,comparison)
                record['comparisons'][key]=str(comparison)
                save(destination,record)
        record.update(status=DONE,finished_utc=now(),current_job=None)
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
