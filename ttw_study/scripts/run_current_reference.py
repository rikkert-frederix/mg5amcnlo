#!/usr/bin/env python3
"""Independent LO production-current normalization at one common fixed scale.

Generate stable tops plus the direct associated current without decay chains:
this uses native diagram-based genps_born, not the factorized current helper.
Compare it times the two independently checked top leptonic branching
densities to a separate all-bw decayed LO run. No associated-W branching
factor is added. This is a normalization pilot, not a narrow-W/full-flavour test.
"""
import argparse
import fcntl
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

from audit_results import audit
from campaign import (ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card,
                      digest, generate, now, run, save)
from harvest_splits import harvest
from inclusive_report import input_bin, production_parameters, validate_bw_normalization
from load_results import canonical_columns, load
from pilot_report import clean
from run_inclusive import launch
from run_sampler_pilots import DONE, sampler_source_guard, verify_archived_control
from rng_history import history_paths, validate_history


def commands(process):
    return '\n'.join([
        'set auto_update 0','set automatic_html_opening False','set notification_center False',
        'set nb_core 64','set run_mode 2','set low_mem_multicore_nlo_generation False',
        'set lhapdf '+str(ROOT/'LHAPDF/bin/lhapdf-config'),
        'set fastjet '+str(STUDY/'local/fastjet/bin/fastjet-config'),
        'import model loop_sm-no_b_mass','define p = g u c d s b u~ c~ d~ s~ b~',
        'generate p p > t t~ mu+ vm QCD=2 QED=2 [QCD]',
        'output fNLO '+str(process),''])


def settings(reference,accuracy,seed):
    result=dict(reference['settings'])
    result.update(iseed=seed,req_acc_fo=accuracy,npoints_fo=1000,niters_fo=3,
                  npoints_fo_grid=200,niters_fo_grid=1,fo_job_target_time=1.,
                  fixed_ren_scale=True,fixed_fac_scale=True,fixed_qes_scale=True,
                  mur_ref_fixed=212.6925,muf_ref_fixed=212.6925,qes_ref_fixed=212.6925)
    return result


def native_run(tag,reference_path,accuracy,sources):
    from madgraph.various.banner import RunCardNLO
    from madgraph.various.histograms import HwUList
    from models.check_param_card import ParamCard
    reference=json.loads(reference_path.read_text())
    if (reference['variant']!='LO' or reference['w_treatment']!='all-bw'
            or reference['decay_bottom_mass']!=0.):
        raise ValueError('Require the frozen massless-decay all-bw LO reference configuration')
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    process=STUDY/'processes'/('TTWplus_current_stable_'+tag)
    if process.exists():
        raise ValueError('Existing independent-current export')
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        command_path=STUDY/'inputs'/(process.name+'.mg5')
        command_path.write_text(commands(process))
        generation=dict(source_hashes=sources,max_cores=64,charge='plus',associated_flavour='mu',
                        factorized_decay_chains=False,command_card_sha256=digest(command_path))
        launch(['taskset','-c','0-63',sys.executable,str(ROOT/'bin/mg5_aMC'),str(command_path)],ROOT,
               STUDY/'logs'/(process.name+'_generation.log'),
               STUDY/'inputs'/(process.name+'_generation.json'),generation)
        if not (process/'Cards/run_card.dat').exists() or (process/'Cards/decay_card.dat').exists():
            raise ValueError('Independent current must be an undecayed native export')
        benchmark_param_card(process,benchmark,'all-bw')
        param=ParamCard(str(process/'Cards/param_card.dat'))
        param['decay'].get((6,)).value=0.  # External on-shell top connectors only.
        if param['decay'].get((24,)).value!=benchmark['W_width']['gamma_nlo_pdf']:
            raise ValueError('Internal associated-current W width must remain finite and matched')
        param.write(str(process/'Cards/param_card.dat'),precision=16)
        (process/'Cards/FO_analyse_card.dat').write_text(
            'FO_ANALYSIS_FORMAT=HwU\nFO_ANALYSE=analysis_HwU_template.o\n'
            'FO_EXTRALIBS=\nFO_EXTRAPATHS=\nFO_INCLUDEPATHS=\n')
        seed=55001
        name='LO_current_fixed_%d' % seed
        archive=process/'study_cards'/name
        archive.mkdir(parents=True)
        options=settings(reference,accuracy,seed)
        card=RunCardNLO(str(process/'Cards/run_card.dat'))
        for key,value in options.items():
            card.set(key,value,user=True,raiseerror=True)
        card.write(str(process/'Cards/run_card.dat'),template=str(process/'Cards/run_card_default.dat'))
        for filename in ('run_card.dat','param_card.dat','FO_analyse_card.dat'):
            shutil.copy2(process/'Cards'/filename,archive/filename)
        manifest=dict(variant='LO',run_name=name,settings=options,source_hashes=sources,
                      benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
                      reference_manifest=str(reference_path),reference_sha256=digest(reference_path),
                      files={p.name:digest(p) for p in archive.iterdir()},max_cores=64,
                      scale_definition='fixed production muR=muF=QES=mt+MW/2=212.6925 GeV',
                      phase_space='native diagram-based genps_born; no factorized decay/core map',
                      genps_born_sha256=digest(process/'SubProcesses/genps_born.f90'),
                      external_amplitude_widths_zero=[6],internal_W_width_GeV=benchmark['W_width']['gamma_nlo_pdf'],
                      no_branching_factor_applied=True)
        save(archive/'manifest.json',manifest)
        execution=dict(variant='LO',charge='plus',pilot=True,max_cores=64,source_hashes=sources)
        launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
                'LO','-f','-n',name,'--multicore','--nb_core=64'],process,
               STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
        hwu=process/'Events'/name/'MADatNLO.HwU'
        if not hwu.exists():
            execution['status']='failed: no final HwU'
            save(archive/'execution.json',execution)
            raise ValueError('No independent current output')
        totals=[h for h in HwUList(str(hwu),raw_labels=True) if h.title.strip()=='total rate']
        if len(totals)!=1 or totals[0].bins[0].boundaries!=(.5,1.5):
            raise ValueError('Unexpected inclusive current histogram')
        row=totals[0].bins[0]
        columns=list(row.wgts)
        weights=[row.wgts[columns[i]] for i in canonical_columns(columns,'LO')]
        execution.update(output=str(hwu),output_sha256=digest(hwu),canonical_weights_pb=weights,
                         nominal_mc_error_pb=row.wgts['dy'])
        history={str(p.relative_to(process)):p.read_bytes() for p in history_paths(process,name)}
        execution['rng_history_audit']=validate_history(
            Path(execution['log']).read_text(),seed,history)
        execution['rng_history_hashes']={name:digest(process/name) for name in history}
        save(archive/'execution.json',execution)
    return archive/'execution.json'


def compare(native_path,audit_path,output):
    if output.exists():
        raise ValueError('Refuse to overwrite an independent-current comparison')
    native=json.loads(native_path.read_text())
    manifest=json.loads((native_path.parent/'manifest.json').read_text())
    if (native['status']!='finished' or digest(native['output'])!=native['output_sha256']
            or any(digest(native_path.parent/name)!=sha for name,sha in manifest['files'].items())):
        raise ValueError('Independent current archive changed')
    decayed=load(audit_path)
    dm=decayed['report']['manifest']
    generation=decayed['report']['execution']['export_generation']['args']
    if (generation['charge']!='plus' or generation['flavours']!=['e','e','mu']
            or dm['variant']!='LO' or dm['w_treatment']!='all-bw'
            or dm['production_scale']!='fixed' or dm['decay_bottom_mass']!=0.):
        raise ValueError('Wrong decayed process/order/scale for the independent LO current')
    expected_settings=dict(manifest['settings'],iseed=dm['settings']['iseed'])
    if dm['settings']!=expected_settings or dm['settings']['iseed']==manifest['settings']['iseed']:
        raise ValueError('Independent current and decayed run settings do not match')
    archive=Path(decayed['report']['path']).parents[2]/'study_cards'/dm['run_name']
    if production_parameters(archive/'param_card.dat')!=production_parameters(native_path.parent/'param_card.dat'):
        raise ValueError('Independent production parameters differ')
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    if digest(STUDY/'inputs/benchmark.json')!=manifest['benchmark_sha256']:
        raise ValueError('Independent-current physical benchmark changed')
    branching=benchmark['W_width']['branching_e']
    proof=validate_bw_normalization(dm,branching)
    if dm['w_width']!=manifest['internal_W_width_GeV'] or not manifest['no_branching_factor_applied']:
        raise ValueError('Wrong physical current width/normalization')
    index=input_bin(decayed)
    actual=decayed['values'][index]
    expected=np.asarray(native['canonical_weights_pb'])*branching**2
    error=np.hypot(decayed['errors'][index],native['nominal_mc_error_pb']*branching**2)
    report=dict(created_utc=now(),status='independent LO current normalization pilot; inspect convergence',
                inputs={str(p):digest(p) for p in (native_path,audit_path)},
                branching_factor=branching**2,associated_branching_factor_applied=False,
                top_current_normalization=proof,decayed_pb=actual,expected_pb=expected,
                residuals_pb=actual-expected,nominal_combined_mc_error_pb=error,
                nominal_pull=(actual[0]-expected[0])/error,
                scope='One ordered assignment, full-current LO at matched fixed production scales. '
                      'Native production map versus factorized production/decay maps. '
                      'No narrow-W limit, NLO normalization, fiducial comparison or varied-weight MC errors certified.')
    save(output,clean(report))
    print('Independent LO current residual:',report['nominal_pull'],'combined nominal MC errors',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--reference-manifest',required=True,type=Path)
    parser.add_argument('--accuracy',type=float,default=.005)
    parser.add_argument('--after-sampler-source-guard',action='store_true',
                        help='explicitly accept only a completed archived run followed by its source guard')
    args=parser.parse_args()
    if not args.tag.replace('_','').isalnum() or not 0.<args.accuracy<1.:
        raise ValueError('Invalid current reference tag or accuracy')
    destination=STUDY/'inputs'/('current_reference_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing current-reference queue')
    args.reference_manifest=args.reference_manifest.resolve()
    extra=[Path(__file__).resolve(),args.reference_manifest,STUDY/'inputs/benchmark.json',
           STUDY/'inputs/bw_branching_normalization.json',
           *[STUDY/'scripts'/name for name in ('run_inclusive.py','inclusive_report.py','load_results.py',
             'audit_results.py','harvest_splits.py','bw_branching.py','run_sampler_pilots.py','rng_history.py')],
           ROOT/'Template/fNLO/FixedOrderAnalysis/analysis_HwU_template.f90',
           ROOT/'Template/fNLO/SubProcesses/genps_born.f90']
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in extra})
    record=dict(created_utc=now(),status='waiting for sampler pilots',source_hashes=sources,
                max_cores=64,after_queue=str(args.after_queue.resolve()),pilot=True,jobs=[],
                after_sampler_source_guard=args.after_sampler_source_guard)
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen independent-current input/source changed')
    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                if not args.after_sampler_source_guard:
                    raise ValueError('Sampler queue stopped; inspect before the independent reference')
                control=sampler_source_guard(previous)
                record['source_guard_control']=control
                record['source_guard_evidence']=verify_archived_control(control)
                break
            if previous['status']==DONE:
                break
            time.sleep(15)
        unchanged()
        record.update(status='running native production current',predecessor_sha256=digest(args.after_queue))
        save(destination,record)
        native=native_run(args.tag,args.reference_manifest,args.accuracy,sources)
        record['jobs'].append(dict(kind='native LO current',execution=str(native)))
        save(destination,record)
        unchanged()
        process=STUDY/'processes'/('TTWplus_all-bw_eemu_mb0p0_both_'+args.tag)
        generate(argparse.Namespace(charge='plus',w_treatment='all-bw',flavours=['e','e','mu'],
                                    corrected='both',decay_bottom_mass=0.,tag=args.tag))
        name='LO_all-bw_fixed_separate_55002_wcurrent'
        record.update(status='running matched decayed LO current',current_job=dict(process=str(process),run=name))
        save(destination,record)
        run(argparse.Namespace(process_dir=str(process),variant='LO',seed=55002,
                               production_sampling='w-current',accuracy=args.accuracy,points=1000,
                               grid_points=200,iterations=3,job_seconds=1.,grid_reference=None,
                               production_scale='fixed',ecm=13000.,main_run=False))
        with (STUDY/'mg5.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
            audit(process/'Events'/name/'MADatNLO.HwU')
            audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
            workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
            harvest(audit_path,workers)
        record['jobs'].append(dict(kind='factorized LO current',audit=str(audit_path),workers=str(workers)))
        report=STUDY/'results'/(args.tag+'_inclusive_current_normalization.json')
        compare(native,audit_path,report)
        record.update(status='independent LO current pilot finished; inspect agreement and convergence',
                      comparison=str(report),finished_utc=now(),current_job=None)
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
