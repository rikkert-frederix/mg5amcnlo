#!/usr/bin/env python3
"""Serialized stable-production and W/bottom-mass validation pilots.

These are not main-campaign precision samples. Every export is fresh, uses
the frozen benchmark, and is audited before proceeding. Raw worker vectors
are preserved before the next prescription reconfigures that export.
"""
import argparse
import fcntl
import json
from pathlib import Path
import time

from audit_results import audit
from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, run, save
from harvest_splits import harvest
from run_inclusive import run as run_stable


def cases(charges):
    return [dict(charge=charge,w_treatment=mode,decay_bottom_mass=mass)
            for charge in charges
            for mode,mass in (('top-bw',0.),('all-bw',0.),('onshell',4.8),('all-bw',4.8))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--charges',nargs='+',choices=('plus','minus'),default=['plus','minus'])
    parser.add_argument('--accuracy',type=float,default=.05)
    args = parser.parse_args()
    if not args.tag.replace('_','').isalnum() or len(set(args.charges)) != len(args.charges):
        raise ValueError('Invalid tag or repeated charge')
    destination = STUDY/'inputs'/('validation_pilot_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Queue already exists; inspect its record before resuming')
    sources = {name:digest(ROOT/name) for name in RUNTIME_SOURCES}
    sources.update({str(path.relative_to(ROOT)):digest(path) for path in (
        Path(__file__).resolve(),STUDY/'scripts/run_inclusive.py',STUDY/'scripts/harvest_splits.py')})
    record = dict(started_utc=now(),status='waiting for complete on-shell pilot queue',
                  after_queue=str(args.after_queue.resolve()),source_hashes=sources,
                  max_cores=64,jobs=[],current_job=None,pilot=True)
    save(destination,record)

    def unchanged():
        if any(digest(ROOT/name) != checksum for name,checksum in sources.items()):
            raise RuntimeError('Source changed after queue preparation; inspect and start a fresh queue')

    try:
        while True:
            predecessor = json.loads(args.after_queue.read_text())
            if predecessor['status'] == 'stopped':
                raise RuntimeError('Preceding queue stopped; inspect it before starting validation')
            if predecessor['status'] == 'pilot queue finished; main precision not certified':
                break
            time.sleep(15)
        record['after_queue_sha256'] = digest(args.after_queue)
        references = {}
        for job in predecessor['jobs']:
            result = json.loads(Path(job['audit']).read_text())
            if result['variant'] == 'S':
                charge = result['execution']['export_generation']['args']['charge']
                references[charge] = Path(job['process'])/'study_cards'/job['run']/'manifest.json'
        # Physically verify the fixed-count saved-grid controller before it
        # supplies a future main-campaign replica. Reuse only the same audited
        # r5 training stage, never any of its histogram estimates.
        unchanged()
        control_tag=args.tag+'_fixed_count'
        process=STUDY/'processes'/('TTWplus_onshell_eemu_mb0p0_both_'+control_tag)
        control_seed=39001
        control_name='S_onshell_core-w-ht-half_separate_%d' % control_seed
        record.update(status='checking physical fixed-count grid resume',
                      current_job=dict(process=str(process),run=control_name))
        save(destination,record)
        generate(argparse.Namespace(charge='plus',w_treatment='onshell',flavours=['e','e','mu'],
                                    corrected='both',decay_bottom_mass=0.,tag=control_tag))
        grid_reference=STUDY/'processes/TTWplus_onshell_eemu_mb0p0_both_r5/study_cards/S_onshell_core-w-ht-half_separate_32001/manifest.json'
        run(argparse.Namespace(process_dir=str(process),variant='S',seed=control_seed,
                               accuracy=-1.,points=65536,grid_points=200,job_seconds=1.,
                               grid_reference=grid_reference,iterations=1,
                               production_scale='core-w-ht-half',ecm=13000.,main_run=False))
        with (STUDY/'mg5.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
            audit(process/'Events'/control_name/'MADatNLO.HwU')
            audit_path=STUDY/'results/audits'/(process.name+'_'+control_name+'.json')
            batch_path=STUDY/'results/workers'/(process.name+'_'+control_name+'.tar.gz')
            harvest(audit_path,batch_path)
        batches=json.loads(batch_path.with_suffix('.json').read_text())
        counts={}
        for worker in batches['workers']:
            job=worker['job']
            if job['niters_done'] != 1 or job['accuracy'] != 0.:
                raise ValueError('Fixed-count resume retained an adaptive target or changed iterations')
            key=(job['p_dir'],job['channel'])
            counts[key]=counts.get(key,0)+job['npoints_done']
        if len(counts) != 2 or set(counts.values()) != {65536}:
            raise ValueError('Physical resumed sample did not use the requested independent counts')
        record['jobs'].append(dict(kind='physical fixed-count resume passed',process=str(process),
                                   run=control_name,audit=str(audit_path),workers=str(batch_path),
                                   points_per_stratum=65536,finished_utc=now()))
        save(destination,record)
        for i,charge in enumerate(args.charges):
            unchanged()
            record.update(status='running stable-production references',current_job=dict(charge=charge))
            save(destination,record)
            stable_tag = args.tag+'_normalization'
            run_stable(argparse.Namespace(charge=charge,tag=stable_tag,
                                           reference_manifest=str(references[charge]),
                                           accuracy=.01,seed=41001+100*i))
            process = STUDY/'processes'/('TTW%s_stable_%s' % (charge,stable_tag))
            record['jobs'].append(dict(kind='stable LO/P',charge=charge,process=str(process),finished_utc=now()))
            save(destination,record)
        for i,case in enumerate(cases(args.charges)):
            unchanged()
            generation = argparse.Namespace(**case,flavours=['e','e','mu'],corrected='both',tag=args.tag)
            process = STUDY/'processes'/('TTW%s_%s_eemu_mb%s_both_%s' % (
                case['charge'],case['w_treatment'],str(case['decay_bottom_mass']).replace('.','p'),args.tag))
            record.update(status='generating robustness pilot',current_job=dict(process=str(process),case=case))
            save(destination,record)
            generate(generation)
            for j,variant in enumerate(('S','Pi')):
                unchanged()
                seed = 43001+100*i+10*j
                name = '%s_%s_core-w-ht-half_separate_%d' % (variant,case['w_treatment'],seed)
                if case['decay_bottom_mass']:
                    name += '_mb%s' % ('%.12g' % case['decay_bottom_mass']).replace('.','p')
                record.update(status='running robustness pilot',current_job=dict(process=str(process),run=name,case=case))
                save(destination,record)
                run(argparse.Namespace(process_dir=str(process),variant=variant,seed=seed,
                                       accuracy=args.accuracy,points=1000,grid_points=200,
                                       job_seconds=1.,grid_reference=None,iterations=3,
                                       production_scale='core-w-ht-half',ecm=13000.,main_run=False))
                # Prevent another study launcher from reconfiguring the process
                # while its still-present split vectors are being preserved.
                with (STUDY/'mg5.lock').open('a') as lock:
                    fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
                    hwu = process/'Events'/name/'MADatNLO.HwU'
                    audit(hwu)
                    audit_path = STUDY/'results/audits'/(process.name+'_'+name+'.json')
                    batch_path = STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
                    harvest(audit_path,batch_path)
                record['jobs'].append(dict(kind='decayed S/Pi',case=case,variant=variant,
                                           process=str(process),run=name,audit=str(audit_path),
                                           workers=str(batch_path),finished_utc=now()))
                save(destination,record)
        record.update(status='validation pilots finished; scientific validation still required',
                      current_job=None,finished_utc=now())
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__ == '__main__':
    main()
