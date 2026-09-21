#!/usr/bin/env python3
"""Fresh S/Pi kernel pilots with stage-audited RNG after the independent LO check.

These replace the unexecuted/cross-round-unverified sampler controls without
overwriting them or transferring their grids. They are not main statistics.
"""
import argparse
import fcntl
import json
import math
from pathlib import Path
import time

from audit_results import audit
from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, run, save
from check_phase_space_support import run as check_support
from harvest_splits import harvest
from read_splits import read
from run_sampler_pilots import run_name, run_settings

REFERENCE_DONE='independent LO current pilot finished; inspect agreement and convergence'
DONE='stage-audited S/Pi pilots finished; inspect physics and convergence'


def cases(seed_start=57001):
    if not isinstance(seed_start,int) or not 0 < seed_start <= 900000000-4:
        raise ValueError('Invalid control seed range')
    return [dict(w_treatment=mode,variant=variant,production_sampling=sampling,seed=seed_start+i)
            for i,(mode,variant,sampling) in enumerate((
                ('all-bw','S','w-current'),('all-bw','Pi','w-current'),
                ('onshell','S','flat'),('onshell','Pi','flat')))]


def reference_gate(comparison):
    if (comparison['status']!='independent LO current normalization pilot; inspect convergence'
            or comparison['associated_branching_factor_applied']
            or len(comparison['decayed_pb'])!=183 or len(comparison['expected_pb'])!=183):
        raise ValueError('Require the matched complete independent LO current comparison')
    expected=comparison['expected_pb'][0]
    error=comparison['nominal_combined_mc_error_pb']
    pull=comparison['nominal_pull']
    if (not all(math.isfinite(v) for v in (expected,error,pull)) or expected<=0. or error<=0.
            or abs(pull)>3. or error/expected>.01):
        raise ValueError('Independent LO nominal normalization is unresolved; inspect before further pilots')
    return dict(nominal_pull=pull,nominal_relative_combined_error=error/expected,
                limitations='Nominal pilot compatibility only, not full convergence or varied-weight coverage')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',type=Path,required=True)
    parser.add_argument('--accuracy',type=float,default=.03)
    parser.add_argument('--seed-start',type=int,default=57001)
    parser.add_argument('--after-alignment-pilot',type=Path)
    args=parser.parse_args()
    selected_cases=cases(args.seed_start)
    if not args.tag.replace('_','').isalnum() or not 0.<args.accuracy<1.:
        raise ValueError('Invalid stage-audited pilot tag/accuracy')
    destination=STUDY/'inputs'/('rng_pilot_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing stage-audited queue')
    extra=[Path(__file__).resolve(),STUDY/'inputs/benchmark.json',
           *[STUDY/'scripts'/name for name in (
               'run_sampler_pilots.py','audit_results.py','audit_virtuals.py','load_results.py',
               'harvest_splits.py','read_splits.py','rng_history.py','check_phase_space_support.py')],
           *[ROOT/'tests/input_files/fks_decay'/name for name in (
               'phase_space_test_dimensions.f90','bw_support_checks.f90','current_proposal_checks.f90')],
           ROOT/'Template/fNLO/Source/kin_functions.f90']
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in extra})
    record=dict(created_utc=now(),status='waiting for independent LO current reference',
                source_hashes=sources,after_queue=str(args.after_queue.resolve()),max_cores=64,
                cases=selected_cases,jobs=[],pilot=True,charge='plus',flavours=['e','e','mu'],
                limitations='No full flavour sum, main precision, narrow-W or massive-bottom completion')
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen stage-audited pilot source changed')
    try:
        if args.after_alignment_pilot:
            record.update(status='waiting for full Born-alignment pilot audit',
                          after_alignment_pilot=str(args.after_alignment_pilot.resolve()))
            save(destination,record)
            while True:
                pilot=json.loads(args.after_alignment_pilot.read_text())
                if pilot['status']=='stopped':
                    raise ValueError('Born-alignment integration failed; inspect before controls')
                if pilot['status'].startswith('completed technical pilot;'):
                    break
                time.sleep(15)
            if any(digest(ROOT/p)!=sha for p,sha in pilot['source_hashes'].items()):
                raise ValueError('Alignment pilot does not validate the current runtime source')
            batches=Path(pilot['batches'])
            batch_metadata=json.loads(batches.with_suffix('.json').read_text())
            if (digest(batches)!=batch_metadata['arrays_sha256'] or
                    not batch_metadata['status'].startswith('equal-count conditional') or
                    not batch_metadata.get('refinement_rng_audit')):
                raise ValueError('Alignment pilot lacks successful unchanged worker/RNG audit')
            record['alignment_pilot_gate']={
                'path':str(args.after_alignment_pilot.resolve()),
                'sha256':digest(args.after_alignment_pilot),'batches_sha256':digest(batches)}
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                raise ValueError('Independent LO predecessor stopped; inspect before further pilots')
            if previous['status']==REFERENCE_DONE:
                break
            time.sleep(15)
        unchanged()
        report_path=Path(previous['comparison'])
        report=json.loads(report_path.read_text())
        if any(digest(p)!=sha for p,sha in report['inputs'].items()):
            raise ValueError('Independent LO comparison input changed')
        record.update(predecessor_sha256=digest(args.after_queue),reference_comparison=str(report_path),
                      reference_comparison_sha256=digest(report_path),reference_gate=reference_gate(report))
        exports={}
        for case in selected_cases:
            unchanged()
            mode=case['w_treatment']
            process=STUDY/'processes'/('TTWplus_%s_eemu_mb0p0_both_%s' % (mode,args.tag))
            name=run_name(case)
            record.update(status='preparing fresh stage-audited pilot',
                          current_job=dict(process=str(process),run=name,case=case))
            save(destination,record)
            if mode not in exports:
                generate(argparse.Namespace(charge='plus',w_treatment=mode,flavours=['e','e','mu'],
                                            corrected='both',decay_bottom_mass=0.,tag=args.tag))
                support=STUDY/'inputs'/('%s_%s_phase_space_support.json' % (args.tag,mode))
                check_support(support,process)
                exports[mode]=dict(process=str(process),local_support=str(support))
                record['exports']=exports
            unchanged()
            record['status']='running fresh stage-audited pilot'
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
            record['jobs'].append(dict(case=case,variant=case['variant'],process=str(process),run=name,
                                       audit=str(audit_path),workers=str(workers),batches=str(batches),
                                       analytic_virtual_checks=str(virtuals),finished_utc=now()))
            save(destination,record)
        record.update(status=DONE,finished_utc=now(),current_job=None)
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
