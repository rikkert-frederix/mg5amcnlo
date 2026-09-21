#!/usr/bin/env python3
"""Serialized absolute LO/NLO stable-W and full-current normalization checks.

Both charges use native stable-top production with all 101 PDF members.
Native stable-W runs use HT/2, equal to the decayed W-system scale;
native full-current runs and their decayed companions use mt+MW/2.
"""
import argparse
import copy
import fcntl
import json
from pathlib import Path
import shutil
import sys
import time

from audit_results import audit as audit_decayed
from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card, digest, generate, now, run, save
from compare_current_batches import batch_record
from harvest_splits import harvest
from native_inclusive import audit as audit_native, native_vectors
from native_current_batches import read as read_native
from read_splits import read as read_decayed
from rng_history import history_paths, validate_history
from run_inclusive import launch
from run_scale_validation import DONE as PREDECESSOR_DONE
from run_width_mass_pilots import verified_pairs

DONE='absolute normalization integrations finished; inspect joint comparisons and convergence'


def cases(seed_start=82001):
    if not isinstance(seed_start,int) or not 0<seed_start<=900000000-14:
        raise ValueError('Invalid normalization seed range')
    rows=[dict(kind=kind,charge=charge,variant=variant)
          for kind in ('stable','current') for charge in ('plus','minus') for variant in ('LO','P')]
    # The completed technical predecessor provides plus-current fixed S/Pi.
    rows += [dict(kind='decayed',charge=charge,variant=variant)
             for charge in ('plus','minus')
             for variant in (('LO','P') if charge=='plus' else ('LO','P','S','Pi'))]
    return [dict(row,seed=seed_start+i) for i,row in enumerate(rows)]


def stopped_prefix(previous):
    """Validate and return the completed prefix of a stopped queue."""
    if previous.get('status')!='stopped':
        raise ValueError('Recovery requires a stopped normalization queue')
    selected=previous.get('cases')
    failed=previous.get('current_case')
    completed=previous.get('jobs')
    if not isinstance(selected,list) or not isinstance(completed,list) or not isinstance(failed,dict):
        raise ValueError('Stopped normalization queue lacks a recoverable case record')
    if len({case.get('seed') for case in selected})!=len(selected):
        raise ValueError('Stopped normalization queue has duplicate seeds')
    try:
        index=selected.index(failed)
    except ValueError as error:
        raise ValueError('Stopped case is absent from its selected queue') from error
    if (len(completed)!=index
            or [job.get('case') for job in completed]!=selected[:index]):
        raise ValueError('Completed normalization jobs are not an ordered prefix')
    return selected,failed,completed,index


def recovery_plan(previous,replacement_seed):
    """Return a fresh logical queue after one stopped, fully archived case.

    Earlier completed jobs are retained only after the caller revalidates their
    batch evidence.  The failed case gets a new seed and every later case is
    run in a new export tag, so no existing output is restarted or overwritten.
    """
    selected,failed,completed,index=stopped_prefix(previous)
    if (not isinstance(replacement_seed,int) or isinstance(replacement_seed,bool)
            or not 0<replacement_seed<=900000000
            or replacement_seed in {case['seed'] for case in selected}):
        raise ValueError('Recovery seed must be new and in the valid RANMAR range')
    repaired=copy.deepcopy(selected)
    replacement=copy.deepcopy(failed)
    replacement['seed']=replacement_seed
    repaired[index]=replacement
    return dict(cases=repaired,carried_jobs=copy.deepcopy(completed),failed_case=copy.deepcopy(failed),
                replacement_case=replacement,failed_index=index)


def adoption_plan(previous,adopted_job):
    """Continue after a reader-only stop using a separately re-audited case."""
    selected,failed,completed,index=stopped_prefix(previous)
    if not isinstance(adopted_job,dict) or adopted_job.get('case')!=failed:
        raise ValueError('Adopted normalization job does not match the stopped case')
    return dict(cases=copy.deepcopy(selected),carried_jobs=copy.deepcopy(completed)+[copy.deepcopy(adopted_job)],
                failed_case=copy.deepcopy(failed),adopted_case=copy.deepcopy(failed),
                adopted_job=copy.deepcopy(adopted_job),failed_index=index)


def adopted_native_job(case,batches):
    """Build one queue job from a completed, independently re-audited native run."""
    batches=Path(batches).resolve()
    metadata,record=batch_record(batches)
    audit_path=Path(record['audit']).resolve()
    audit=json.loads(audit_path.read_text())
    manifest=audit['manifest']
    actual=dict(kind=manifest['native_kind'],charge=manifest['charge'],variant=manifest['variant'],
                seed=record['seed'])
    if actual!=case:
        raise ValueError('Adopted native batch evidence has different physical case/seed')
    process=Path(audit['path']).resolve().parents[2]
    if process.parent!=STUDY/'processes':
        raise ValueError('Adopted native batch evidence is outside the study exports')
    run=manifest['run_name']
    execution=process/'study_cards'/run/'execution.json'
    result=dict(case=copy.deepcopy(case),variant=case['variant'],process=str(process),run=run,
                audit=str(audit_path),execution=str(execution),
                workers=str(Path(metadata['input']).with_suffix('.gz')),batches=str(batches))
    execution_record=json.loads(execution.read_text())
    if execution_record.get('status')!='finished' or execution_record.get('returncode')!=0:
        raise ValueError('Adopted native execution is not a successful completed run')
    result['finished_utc']=execution_record['finished_utc']
    return result


def native_commands(process,kind,charge):
    if kind not in ('stable','current') or charge not in ('plus','minus'):
        raise ValueError('Unknown native normalization process')
    final=('w+' if charge=='plus' else 'w-') if kind=='stable' else (
        'mu+ vm' if charge=='plus' else 'mu- vm~')
    return '\n'.join([
        'set auto_update 0','set automatic_html_opening False','set notification_center False',
        'set nb_core 64','set run_mode 2','set low_mem_multicore_nlo_generation False',
        'set lhapdf '+str(ROOT/'LHAPDF/bin/lhapdf-config'),
        'set fastjet '+str(STUDY/'local/fastjet/bin/fastjet-config'),
        'import model loop_sm-no_b_mass','define p = g u c d s b u~ c~ d~ s~ b~',
        'generate p p > t t~ %s QCD=2 QED=%d [QCD]' % (final,1 if kind=='stable' else 2),
        'output fNLO '+str(process),''])


def native_settings(reference,case,accuracy):
    if reference['w_treatment']!='onshell' or reference['production_scale']!='core-w-ht-half':
        raise ValueError('Require the frozen on-shell W-system reference')
    result=dict(reference['settings'])
    fixed=case['kind']=='current'
    result.update(iseed=case['seed'],req_acc_fo=accuracy,npoints_fo=1000,niters_fo=3,
                  npoints_fo_grid=1000,niters_fo_grid=1,fo_job_target_time=1.,
                  fo_job_min_splits=5,
                  fixed_ren_scale=fixed,fixed_fac_scale=fixed,fixed_qes_scale=fixed,
                  mur_ref_fixed=212.6925,muf_ref_fixed=212.6925,qes_ref_fixed=212.6925)
    return result


def generate_native(process,case,sources):
    from models.check_param_card import ParamCard
    if process.exists():
        raise ValueError('Existing native normalization export; never restart implicitly')
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        command_path=STUDY/'inputs'/(process.name+'.mg5')
        command_path.write_text(native_commands(process,case['kind'],case['charge']))
        generation=dict(source_hashes=sources,max_cores=64,charge=case['charge'],
                        native_kind=case['kind'],factorized_decay_chains=False,
                        command_card_sha256=digest(command_path))
        launch(['taskset','-c','0-63',sys.executable,str(ROOT/'bin/mg5_aMC'),str(command_path)],ROOT,
               STUDY/'logs'/(process.name+'_generation.log'),
               STUDY/'inputs'/(process.name+'_generation.json'),generation)
        if not (process/'Cards/run_card.dat').exists() or (process/'Cards/decay_card.dat').exists():
            raise ValueError('Require an undecayed native export')
        benchmark_param_card(process,benchmark,'all-bw' if case['kind']=='current' else 'onshell')
        card=ParamCard(str(process/'Cards/param_card.dat'))
        card['decay'].get((6,)).value=0.
        card['decay'].get((24,)).value=(benchmark['W_width']['gamma_nlo_pdf']
                                      if case['kind']=='current' else 0.)
        card.write(str(process/'Cards/param_card.dat'),precision=16)
        (process/'Cards/FO_analyse_card.dat').write_text(
            'FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE=analysis_HwU_template.o\n'
            'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n')


def native_run(process,case,reference_path,accuracy,sources):
    from madgraph.various.banner import RunCardNLO
    from madgraph.various.histograms import HwUList
    reference=json.loads(reference_path.read_text())
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    kind,variant=case['kind'],case['variant']
    name='%s_%s_%s_%d' % (variant,kind,'fixed' if kind=='current' else 'core_ht_half',case['seed'])
    archive=process/'study_cards'/name
    archive.mkdir(parents=True)
    options=native_settings(reference,case,accuracy)
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        card=RunCardNLO(str(process/'Cards/run_card.dat'))
        for key,value in options.items():
            card.set(key,value,user=True,raiseerror=True)
        card.write(str(process/'Cards/run_card.dat'),template=str(process/'Cards/run_card_default.dat'))
        for filename in ('run_card.dat','param_card.dat','FO_analyse_card.dat'):
            shutil.copy2(process/'Cards'/filename,archive/filename)
        manifest=dict(variant=variant,native_kind=kind,charge=case['charge'],run_name=name,
            settings=options,source_hashes=sources,benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
            reference_manifest=str(reference_path),reference_manifest_sha256=digest(reference_path),
            files={p.name:digest(p) for p in archive.iterdir()},max_cores=64,
            scale_definition=('fixed muR=muF=QES=212.6925 GeV' if kind=='current' else
                              'stable t,tbar,W and production radiation mT sum / 2'),
            production_scale='fixed' if kind=='current' else 'core-w-ht-half',
            phase_space='native diagram-based genps_born; no factorized decay/core map',
            genps_born_sha256=digest(process/'SubProcesses/genps_born.f90'),
            no_branching_factor_applied=True,associated_flavour='mu' if kind=='current' else None,
            external_amplitude_widths_zero=[6] if kind=='current' else [6,24],
            internal_W_width_GeV=benchmark['W_width']['gamma_nlo_pdf'] if kind=='current' else None)
        save(archive/'manifest.json',manifest)
        execution=dict(variant=variant,charge=case['charge'],pilot=True,max_cores=64,source_hashes=sources)
        launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
                'LO' if variant=='LO' else 'NLO','-f','-n',name,'--multicore','--nb_core=64'],process,
               STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
        hwu=process/'Events'/name/'MADatNLO.HwU'
        if not hwu.is_file():
            execution.update(status='failed: no final histogram',missing_output=str(hwu))
            save(archive/'execution.json',execution)
            raise ValueError('Native calculator returned without its final histogram')
        vectors=native_vectors(HwUList(str(hwu),raw_labels=True),variant)
        history={str(p.relative_to(process)):p.read_bytes() for p in history_paths(process,name)}
        execution.update(output=str(hwu),output_sha256=digest(hwu),
            canonical_weights_pb=vectors['values'][0,vectors['indices']].tolist(),
            nominal_mc_error_pb=float(vectors['values'][0,vectors['columns'].index('dy')]),
            branching_e=benchmark['W_width']['branching_e'],
            rng_history_audit=validate_history(Path(execution['log']).read_text(),case['seed'],history),
            rng_history_hashes={name:digest(process/name) for name in history})
        save(archive/'execution.json',execution)
        audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
        workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
        audit_native(archive/'execution.json',audit_path,workers)
    batches=STUDY/'results'/(process.name+'_'+name+'_batches.npz')
    read_native(workers.with_suffix('.json'),batches,variant=variant,kind=kind)
    return dict(case=case,variant=variant,process=str(process),run=name,audit=str(audit_path),
                execution=str(archive/'execution.json'),workers=str(workers),batches=str(batches))


def decayed_run(process,case,accuracy):
    name='%s_all-bw_fixed_separate_%d_wcurrent' % (case['variant'],case['seed'])
    run(argparse.Namespace(process_dir=str(process),variant=case['variant'],seed=case['seed'],
        production_sampling='w-current',accuracy=accuracy,points=1000,grid_points=1000,
        iterations=3,job_seconds=1.,grid_reference=None,production_scale='fixed',ecm=13000.,main_run=False))
    audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
    workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        audit_decayed(process/'Events'/name/'MADatNLO.HwU')
        harvest(audit_path,workers)
    batches=STUDY/'results'/(process.name+'_'+name+'_batches.npz')
    read_decayed(workers.with_suffix('.json'),batches)
    job=dict(case=case,variant=case['variant'],process=str(process),run=name,audit=str(audit_path),
             workers=str(workers),batches=str(batches))
    if case['variant'] in ('S','Pi'):
        virtuals=STUDY/'results'/(process.name+'_'+name+'_virtual_checks.json')
        audit_virtuals(workers.with_suffix('.json'),virtuals)
        job['analytic_virtual_checks']=str(virtuals)
    return job


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--reference-manifest',required=True,type=Path)
    parser.add_argument('--native-accuracy',type=float,default=.005)
    parser.add_argument('--decayed-accuracy',type=float,default=.01)
    parser.add_argument('--seed-start',type=int,default=82001)
    parser.add_argument('--resume-stopped-queue',type=Path)
    parser.add_argument('--replacement-seed',type=int)
    parser.add_argument('--adopt-failed-batches',type=Path)
    args=parser.parse_args()
    if (not args.tag.replace('_','').isalnum()
            or not all(0<a<1 for a in (args.native_accuracy,args.decayed_accuracy))
            or (args.resume_stopped_queue is None and
                (args.replacement_seed is not None or args.adopt_failed_batches is not None))
            or (args.resume_stopped_queue is not None and
                ((args.replacement_seed is None)==(args.adopt_failed_batches is None)))):
        raise ValueError('Invalid normalization tag/accuracy')
    args.reference_manifest=args.reference_manifest.resolve()
    args.after_queue=args.after_queue.resolve()
    recovery=None
    if args.resume_stopped_queue:
        args.resume_stopped_queue=args.resume_stopped_queue.resolve()
        previous=json.loads(args.resume_stopped_queue.read_text())
        if args.adopt_failed_batches:
            args.adopt_failed_batches=args.adopt_failed_batches.resolve()
            recovery=adoption_plan(previous,adopted_native_job(previous['current_case'],args.adopt_failed_batches))
        else:
            recovery=recovery_plan(previous,args.replacement_seed)
        if (Path(previous.get('after_queue','')).resolve()!=args.after_queue
                or digest(args.after_queue)!=previous.get('predecessor_sha256')):
            raise ValueError('Stopped normalization queue has a different technical predecessor')
        reference_key=str(args.reference_manifest.relative_to(ROOT))
        if previous.get('source_hashes',{}).get(reference_key)!=digest(args.reference_manifest):
            raise ValueError('Stopped normalization queue has a different frozen reference manifest')
        selected=recovery['cases']
    else:
        selected=cases(args.seed_start)
    destination=STUDY/'inputs'/('normalization_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing normalization queue; inspect, never overwrite')
    extra=[Path(__file__).resolve(),args.reference_manifest,STUDY/'inputs/benchmark.json',
        ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_template.f90',
        ROOT/'Template/fNLO/SubProcesses/genps_born.f90',
        *[STUDY/'scripts'/name for name in (
            'run_inclusive.py','native_inclusive.py','native_current_batches.py','audit_native_current.py',
            'audit_results.py','audit_virtuals.py','load_results.py','harvest_splits.py','read_splits.py',
            'rng_history.py','run_width_mass_pilots.py','compare_current_batches.py','run_scale_validation.py')]]
    if recovery:
        extra.append(args.resume_stopped_queue)
    if args.adopt_failed_batches:
        extra.append(args.adopt_failed_batches)
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in extra})
    record=dict(created_utc=now(),status='waiting for direct scale checks',source_hashes=sources,
        max_cores=64,after_queue=str(args.after_queue),cases=selected,
        jobs=(recovery['carried_jobs'] if recovery else []),
        native_accuracy=args.native_accuracy,decayed_accuracy=args.decayed_accuracy,
        scope='Both charges, native LO/P stable W and full associated current. '
              'Fresh fixed-scale all-BW LO/P companions and W- S/Pi; W+ fixed S/Pi from predecessor. '
              'Stable reference also tests completed on-shell/top-BW S/Pi controls. '
              'All weights and actual stage streams preserved; no final convergence claim.')
    if recovery:
        record['recovery']=dict(stopped_queue=str(args.resume_stopped_queue),
            stopped_queue_sha256=digest(args.resume_stopped_queue),failed_case=recovery['failed_case'],
            failed_index=recovery['failed_index'],carried_job_count=len(recovery['carried_jobs']))
        if 'replacement_case' in recovery:
            record['recovery'].update(replacement_case=recovery['replacement_case'],
                convention='The stopped case remains archived but is excluded because it has fewer than five '
                           'independent batches per stratum. A new seed and fresh export are required; '
                           'previous successful jobs are revalidated before reuse.')
        else:
            record['recovery'].update(adopted_case=recovery['adopted_case'],
                adopted_batches=str(args.adopt_failed_batches),
                adopted_batches_sha256=digest(args.adopt_failed_batches),
                convention='The stopped native calculation completed; a repaired reader independently '
                           're-audited its immutable raw archive. The validated case is adopted without '
                           'rerunning it, and all prior successful jobs are revalidated before reuse.')
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen normalization source/input changed')
    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                raise ValueError('Technical predecessor stopped; inspect before further MC')
            if previous['status']==PREDECESSOR_DONE:
                break
            time.sleep(15)
        unchanged()
        record.update(status=('running recovered absolute normalization integrations' if recovery else
                              'running absolute normalization integrations'),
                      predecessor_sha256=digest(args.after_queue))
        save(destination,record)
        processes,streams=set(),set()
        for job in record['jobs']:
            pairs=verified_pairs(job)
            if streams & pairs:
                raise ValueError('Carried normalization stage streams overlap')
            streams.update(pairs)
        for case in selected[len(record['jobs']):]:
            unchanged()
            kind,charge=case['kind'],case['charge']
            process=STUDY/'processes'/('TTW%s_%s_%s' % (charge,
                'all-bw_eemu_mb0p0_both' if kind=='decayed' else kind+'_native',args.tag))
            record['current_case']=case
            save(destination,record)
            if process not in processes:
                if kind=='decayed':
                    generate(argparse.Namespace(charge=charge,w_treatment='all-bw',flavours=['e','e','mu'],
                        corrected='both',decay_bottom_mass=0.,tag=args.tag))
                else:
                    generate_native(process,case,sources)
                processes.add(process)
            job=(decayed_run(process,case,args.decayed_accuracy) if kind=='decayed' else
                 native_run(process,case,args.reference_manifest,args.native_accuracy,sources))
            pairs=verified_pairs(job)
            if streams & pairs:
                raise ValueError('Actual normalization stage streams overlap earlier runs')
            streams.update(pairs)
            job['finished_utc']=now()
            record['jobs'].append(job)
            record['distinct_stage_pairs']=len(streams)
            save(destination,record)
        record.update(status=DONE,finished_utc=now(),current_case=None)
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
