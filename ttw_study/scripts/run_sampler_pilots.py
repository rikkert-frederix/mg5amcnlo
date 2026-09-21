#!/usr/bin/env python3
"""Fresh physical checks of the corrected kernel and opt-in current proposal.

Six W+ e,e,mu pilots, never main-campaign precision samples. Accept either
the planned archived-control/source-guard transition or an explicitly supplied
snapshot of a cancelled, nonconvergent control. Never relabel that failed
run as completed. No grids or estimates are transferred to these exports.
"""
import argparse
import fcntl
import json
from pathlib import Path
import time

from audit_results import audit
from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, run, save
from check_phase_space_support import run as check_support
from harvest_splits import harvest

DONE = 'sampler pilots finished; physical agreement and convergence require inspection'


def cases():
    return [dict(w_treatment=mode, variant=variant, production_sampling=sampling,
                 seed=53001+i)
            for i,(mode,variant,sampling) in enumerate((
                ('all-bw','LO','w-current'), ('all-bw','LO','flat'),
                ('all-bw','S','w-current'), ('all-bw','Pi','w-current'),
                ('onshell','S','flat'), ('onshell','Pi','flat')))]


def run_name(case):
    name='%s_%s_core-w-ht-half_separate_%d' % (
        case['variant'],case['w_treatment'],case['seed'])
    return name+('_wcurrent' if case['production_sampling']=='w-current' else '')


def run_settings(process,case,accuracy):
    # The linear map is a deliberately bounded diagnostic after the previous
    # adaptive attempt escalated by 100x. Its measured error remains binding;
    # reaching this fixed count does not establish map agreement.
    fixed=(case['w_treatment']=='all-bw' and case['variant']=='LO'
           and case['production_sampling']=='flat')
    return argparse.Namespace(process_dir=str(process),variant=case['variant'],seed=case['seed'],
                              production_sampling=case['production_sampling'],
                              accuracy=-1. if fixed else accuracy,
                              points=262144 if fixed else 1000,grid_points=200,
                              iterations=1 if fixed else 3,job_seconds=1.,grid_reference=None,
                              production_scale='core-w-ht-half',ecm=13000.,main_run=False)


def drained_control(previous):
    """Only accept the planned source-guard stop, not an arbitrary failed run."""
    if previous['status'] != 'stopped':
        return None
    expected=repr(RuntimeError('Source changed after queue preparation; inspect and start a fresh queue'))
    if previous.get('error') != expected:
        raise ValueError('Predecessor stopped for a reason other than the planned source guard')
    target=previous.get('current_job',{})
    jobs=[job for job in previous['jobs'] if job.get('process')==target.get('process')
          and job.get('run')==target.get('run')]
    if len(jobs)!=1 or jobs[0].get('case')!=dict(
            charge='plus',w_treatment='all-bw',decay_bottom_mass=0.) or jobs[0].get('variant')!='S':
        raise ValueError('The old all-bw S control must finish and be archived before transition')
    return jobs[0]


def cancelled_control(previous, diagnostic):
    if previous['status'] != 'stopped':
        return None
    target=previous.get('current_job',{})
    if (diagnostic.get('status')!='frozen nonconverged adaptive diagnostic; excluded from physics'
            or target.get('process')!=diagnostic.get('process')
            or target.get('run')!=diagnostic.get('run')
            or target.get('case')!=dict(charge='plus',w_treatment='all-bw',decay_bottom_mass=0.)
            or not previous.get('error','').startswith("RuntimeError('Integration failed or lacks HwU output;")):
        raise ValueError('Failed-control snapshot does not match the stopped integration')
    if any(j.get('process')==target['process'] and j.get('run')==target['run']
           for j in previous['jobs']):
        raise ValueError('A failed-control transition must not relabel a completed result')
    return dict(process=target['process'],run=target['run'],
                diagnostic_archive=diagnostic['archive'],excluded_from_physics=True)


def sampler_source_guard(previous):
    """Inspect a planned drain of this pilot queue, never an integration failure."""
    if previous['status']!='stopped':
        return None
    if previous.get('error')!=repr(ValueError(
            'Frozen sampler source/input changed; use a fresh inspected queue')):
        raise ValueError('Sampler predecessor did not stop at its source guard')
    current=previous.get('current_job') or {}
    matches=[j for j in previous['jobs'] if j.get('process')==current.get('process')
             and j.get('run')==current.get('run')]
    if (len(matches)!=1 or matches[0]!=previous['jobs'][-1]
            or matches[0]['case']!=current.get('case')):
        raise ValueError('Sampler source guard must follow its last fully archived run')
    return matches[0]


def verify_archived_control(control):
    audit_path=Path(control['audit'])
    workers=Path(control['workers'])
    report=json.loads(audit_path.read_text())
    record=json.loads(workers.with_suffix('.json').read_text())
    if (report['execution']['status']!='finished' or digest(report['path'])!=report['output_sha256']
            or digest(workers)!=record['archive_sha256'] or digest(audit_path)!=record['audit_sha256']
            or Path(report['path']).parent.name!=control['run']
            or str(Path(report['path']).parents[2])!=control['process']):
        raise ValueError('Source-guard transition requires an unchanged completed output/archive')
    return dict(audit_sha256=digest(audit_path),workers_record_sha256=digest(workers.with_suffix('.json')))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--failed-control',type=Path,
                        help='explicit immutable snapshot of the cancelled nonconvergent predecessor')
    parser.add_argument('--accuracy',type=float,default=.03)
    args=parser.parse_args()
    if not args.tag.replace('_','').isalnum() or not 0.<args.accuracy<1.:
        raise ValueError('Invalid sampler pilot tag/accuracy')
    destination=STUDY/'inputs'/('sampler_pilot_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing sampler queue; inspect its state before further work')
    extra=[Path(__file__).resolve(),STUDY/'inputs/benchmark.json',
           *[STUDY/'scripts'/name for name in ('audit_results.py','audit_virtuals.py',
             'check_phase_space_support.py','harvest_splits.py','rng_history.py')],
           *[ROOT/'tests/input_files/fks_decay'/name for name in (
             'phase_space_test_dimensions.f90','bw_support_checks.f90','current_proposal_checks.f90')],
           ROOT/'Template/fNLO/Source/kin_functions.f90']
    diagnostic=None
    if args.failed_control:
        args.failed_control=args.failed_control.resolve()
        diagnostic=json.loads(args.failed_control.read_text())
        extra.append(args.failed_control)
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in extra})
    record=dict(created_utc=now(),status='waiting for archived old-source control',
                after_queue=str(args.after_queue.resolve()),source_hashes=sources,
                max_cores=64,pilot=True,charge='plus',flavours=['e','e','mu'],
                cases=cases(),jobs=[],current_job=None,
                failed_control_snapshot=str(args.failed_control) if args.failed_control else None,
                linear_LO_control='fresh training, 262144 points per stratum, one fixed refinement; '
                                  'convergence and map agreement not guaranteed by this count',
                limitations='Independent proposal/kernel pilots only. No full flavour sum, '
                            'narrow-W limit, main precision or agreement certified by queue completion.')
    save(destination,record)

    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen sampler source/input changed; use a fresh inspected queue')

    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            control=(cancelled_control(previous,diagnostic) if diagnostic else drained_control(previous))
            if control is not None:
                break
            time.sleep(15)
        unchanged()
        # Read-only evidence for the intentional transition; never modify the
        # old queue, its process, cards, source fingerprints or raw vectors.
        if diagnostic:
            failed_execution=Path(control['process'])/'study_cards'/control['run']/'execution.json'
            failed=json.loads(failed_execution.read_text())
            if (failed['status']!='failed' or failed['outputs']
                    or failed['child_pid']!=diagnostic['calculator_pid']
                    or (Path('/proc')/str(failed['child_pid'])).exists()
                    or digest(diagnostic['archive'])!=diagnostic['archive_sha256']):
                raise ValueError('Require a stopped, failed execution and unchanged diagnostic archive')
            record.update(cancelled_old_source_control=control,
                          failed_execution_sha256=digest(failed_execution),
                          failed_snapshot_sha256=digest(args.failed_control))
        else:
            audit_record=json.loads(Path(control['audit']).read_text())
            workers=Path(control['workers']).with_suffix('.json')
            worker_record=json.loads(workers.read_text())
            if (audit_record['execution']['status']!='finished'
                    or digest(audit_record['path'])!=audit_record['output_sha256']
                    or digest(control['workers'])!=worker_record['archive_sha256']
                    or digest(control['audit'])!=worker_record['audit_sha256']):
                raise ValueError('Old-source control did not preserve a completed consistent archive')
            record.update(old_source_control=control,old_control_audit_sha256=digest(control['audit']))
        record['predecessor_sha256']=digest(args.after_queue)
        save(destination,record)
        exports={}
        for case in cases():
            unchanged()
            mode=case['w_treatment']
            process=STUDY/'processes'/('TTWplus_%s_eemu_mb0p0_both_%s' % (mode,args.tag))
            name=run_name(case)
            record.update(status='preparing fresh sampler pilot',
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
            record['status']='running fresh sampler pilot'
            save(destination,record)
            run(run_settings(process,case,args.accuracy))
            with (STUDY/'mg5.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
                audit(process/'Events'/name/'MADatNLO.HwU')
                audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
                workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
                harvest(audit_path,workers)
            job=dict(case=case,variant=case['variant'],process=str(process),run=name,
                     audit=str(audit_path),workers=str(workers),finished_utc=now())
            if case['variant'] in ('S','Pi'):
                virtuals=STUDY/'results'/(args.tag+'_'+name+'_virtual_checks.json')
                audit_virtuals(workers.with_suffix('.json'),virtuals)
                job['analytic_virtual_checks']=str(virtuals)
            record['jobs'].append(job)
            save(destination,record)
        record.update(status=DONE,finished_utc=now(),current_job=None)
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
