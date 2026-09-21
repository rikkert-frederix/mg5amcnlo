#!/usr/bin/env python3
"""Fresh generated mb -> 0 S/Pi checks with unchanged five-flavour production.

One W+ e,e,mu assignment, masses 0, 1 and 0.1 GeV, on-shell and all-BW.
This initial allocation tests continuity; it does not replace the massive
benchmark, full-flavour predictions or independent-retraining convergence.
"""
import argparse
import copy
import fcntl
import json
from pathlib import Path
import shutil
import sys
import time

from audit_results import audit
from audit_virtuals import audit as audit_virtuals
from campaign import ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card, digest, generate, now, save, setup_module
from check_phase_space_support import run as check_support
from harvest_splits import harvest
from production_invariance import run as check_production
from read_splits import read
from run_inclusive import launch
from run_narrow_w_refinements import DONE as PREVIOUS_DONE
from run_width_mass_pilots import verified_pairs
from small_mass_inputs import validate as validate_widths

DONE='generated small-mass pilots finished; inspect continuity and statistical sensitivity'


def cases(seed_start=85001):
    if not isinstance(seed_start,int) or not 0<seed_start<=900000000-12:
        raise ValueError('Invalid small-mass seed range')
    return [dict(charge='plus',flavours=['e','e','mu'],w_treatment=mode,
        decay_bottom_mass=mass,variant=variant,seed=seed_start+i,
        production_sampling='w-current' if mode=='all-bw' else 'flat')
        for i,(mode,mass,variant) in enumerate(
            (mode,mass,variant) for mode in ('onshell','all-bw')
            for mass in (0.,1.,.1) for variant in ('S','Pi'))]


def local_benchmark(width_inputs):
    base=json.loads((STUDY/'inputs/benchmark.json').read_text())
    if width_inputs['benchmark_sha256']!=digest(STUDY/'inputs/benchmark.json'):
        raise ValueError('Small-mass widths refer to a different physical benchmark')
    validate_widths(width_inputs,base)
    # This in-memory extension supplies the explicitly matched widths to
    # the standard parameter writer. The frozen benchmark file is not edited.
    result=copy.deepcopy(base)
    result['top_widths']=copy.deepcopy(width_inputs['top_widths'])
    return result


def configure(process,case,inputs,inputs_path,accuracy,exclude_split_outliers=False):
    benchmark=local_benchmark(inputs)
    mb,width=benchmark_param_card(process,benchmark,case['w_treatment'])
    if mb!=case['decay_bottom_mass']:
        raise ValueError('Actual decay model mass differs from the continuity case')
    setup_module().configure(argparse.Namespace(process_dir=str(process),variant=case['variant'],
        top_width_lo=width['gamma_lo'],top_width_nlo=width['gamma_nlo_pdf'],
        w_treatment=case['w_treatment'],top_width_w_treatment=width['w_treatment'],
        decay_bottom_mass=mb,top_width_bottom_mass=mb,
        width_source=str(inputs_path)+'; sha256='+digest(inputs_path),pdf_id=331700,
        production_scale='core-w-ht-half',production_sampling=case['production_sampling'],
        decay_scales='separate',ecm=13000.,seed=case['seed'],points=1000,grid_points=1000,
        iterations=3,accuracy=accuracy,job_seconds=1.))
    base_name='%s_%s_core-w-ht-half_separate_%d'%(case['variant'],case['w_treatment'],case['seed'])
    if case['production_sampling']=='w-current':
        base_name+='_wcurrent'
    if mb:
        base_name+='_mb'+('%.12g'%mb).replace('.','p')
    base=process/'study_cards'/base_name
    manifest=json.loads((base/'manifest.json').read_text())
    if exclude_split_outliers:
        from madgraph.various.banner import RunCardNLO
        settings=dict(fo_split_outlier_threshold=10.,fo_split_outlier_variance_fraction=.95,
                      fo_split_outlier_min_splits=8)
        card=RunCardNLO(str(process/'Cards/run_card.dat'))
        for key,value in settings.items():
            card.set(key,value,user=True,raiseerror=True)
        card.write(str(process/'Cards/run_card.dat'),template=str(process/'Cards/run_card_default.dat'))
        manifest['settings'].update(settings)
        manifest['split_outlier_policy']=dict(settings=settings,
            limitation='Opt-in data-dependent exclusion; retained-sample errors omit selection bias.')
    name=base_name+'__small_mass'
    archive=process/'study_cards'/name
    if archive.exists() or (process/'Events'/name).exists():
        raise ValueError('Existing generated small-mass run')
    archive.mkdir()
    for filename in ('param_card.dat','run_card.dat','decay_card.dat','FO_analyse_card.dat'):
        shutil.copy2(process/'Cards'/filename,archive/filename)
        manifest['hashes'][filename]=digest(archive/filename)
    manifest.update(run_name=name,small_mass_test=case,small_mass_inputs=str(inputs_path),
        small_mass_inputs_sha256=digest(inputs_path),generated_event_weight_rescaling=1.,
        base_configuration=str(base/'manifest.json'),base_configuration_sha256=digest(base/'manifest.json'),
        status='configured generated small-mass control; no integration performed')
    save(archive/'manifest.json',manifest)
    return name,archive


def run_case(process,case,inputs,inputs_path,accuracy,sources,exclude_split_outliers=False):
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        name,archive=configure(process,case,inputs,inputs_path,accuracy,exclude_split_outliers)
        execution=dict(max_cores=64,pilot=True,small_mass_test=case,source_hashes=sources,
            export_generation=json.loads((STUDY/'inputs'/(process.name+'_generation.json')).read_text()),
            requested_total_rate_accuracy=accuracy,grid_reference=None)
        launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
            'NLO','-f','-n',name,'--multicore','--nb_core=64'],process,
            STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
        hwu=process/'Events'/name/'MADatNLO.HwU'
        if not hwu.is_file():
            execution.update(status='failed: no final histogram',missing_output=str(hwu))
            save(archive/'execution.json',execution)
            raise ValueError('Missing final small-mass histogram; inspect worker logs')
        execution['outputs']={str(hwu):digest(hwu)}
        save(archive/'execution.json',execution)
        audit(hwu)
        audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
        workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
        harvest(audit_path,workers)
    batches=STUDY/'results'/(process.name+'_'+name+'_batches.npz')
    read(workers.with_suffix('.json'),batches)
    virtuals=STUDY/'results'/(process.name+'_'+name+'_virtuals.json')
    audit_virtuals(workers.with_suffix('.json'),virtuals)
    return dict(case=case,variant=case['variant'],process=str(process),run=name,audit=str(audit_path),
        workers=str(workers),batches=str(batches),analytic_virtual_checks=str(virtuals),finished_utc=now())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--width-inputs',required=True,type=Path)
    parser.add_argument('--accuracy',type=float,default=.01)
    parser.add_argument('--seed-start',type=int,default=85001)
    parser.add_argument('--exclude-split-outliers',action='store_true',
        help='Archive and exclude at most one extreme equal-count split per stratum; filtered errors omit selection bias')
    args=parser.parse_args()
    if not args.tag.replace('_','').isalnum() or not 0<args.accuracy<1:
        raise ValueError('Invalid generated continuity settings')
    args.width_inputs=args.width_inputs.resolve()
    inputs=json.loads(args.width_inputs.read_text())
    benchmark=local_benchmark(inputs)
    selected=cases(args.seed_start)
    destination=STUDY/'inputs'/('small_mass_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing small-mass queue; inspect before further action')
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    files=[Path(__file__).resolve(),args.width_inputs,STUDY/'inputs/benchmark.json',
        *[STUDY/'scripts'/name for name in ('small_mass_inputs.py','run_inclusive.py','run_narrow_w_refinements.py',
            'check_phase_space_support.py','production_invariance.py','audit_results.py','audit_virtuals.py',
            'harvest_splits.py','read_splits.py','load_results.py','rng_history.py','run_width_mass_pilots.py',
            'run_narrow_w_pilots.py','narrow_w_inputs.py','run_normalization_checks.py')],
        *[ROOT/'tests/input_files/fks_decay'/name for name in (
            'phase_space_test_dimensions.f90','bw_support_checks.f90')],
        ROOT/'Template/fNLO/Source/kin_functions.f90']
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in files})
    record=dict(created_utc=now(),status='waiting for narrow-W reference retrainings',max_cores=64,
        after_queue=str(args.after_queue.resolve()),width_inputs=str(args.width_inputs),
        source_hashes=sources,cases=selected,jobs=[],exports={},accuracy=args.accuracy,
        scope='Generated IR-safe small-mass S/Pi continuity at mb=0,1,0.1 GeV, W+ e,e,mu, on-shell/all-BW. '
              'Dynamic W-system production scales and matched five-flavour top widths. '
              'Pilot allocation only; remeasure fiducial/shape sensitivity and retraining convergence.')
    record['exclude_split_outliers']=args.exclude_split_outliers
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen small-mass input/source changed')
    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                raise ValueError('Narrow-W reference predecessor stopped')
            if previous['status']==PREVIOUS_DONE:
                break
            time.sleep(15)
        unchanged()
        record.update(status='running generated small-mass controls',predecessor_sha256=digest(args.after_queue))
        references,streams={},set()
        for case in selected:
            unchanged()
            key=case['w_treatment']+'_mb'+str(case['decay_bottom_mass']).replace('.','p')
            process=STUDY/'processes'/('TTWplus_%s_eemu_mb%s_both_%s'%(
                case['w_treatment'],str(case['decay_bottom_mass']).replace('.','p'),args.tag))
            record['current_case']=case
            save(destination,record)
            if key not in record['exports']:
                generate(argparse.Namespace(charge='plus',flavours=['e','e','mu'],corrected='both',
                    w_treatment=case['w_treatment'],decay_bottom_mass=case['decay_bottom_mass'],tag=args.tag))
                support=STUDY/'inputs'/(args.tag+'_'+key+'_phase_space_support.json')
                check_support(support,process)
                exported=dict(process=str(process),local_support=str(support))
                if case['decay_bottom_mass']:
                    identity=STUDY/'inputs'/(args.tag+'_'+key+'_production_identity.json')
                    with (STUDY/'mg5.lock').open('a') as lock:
                        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
                        benchmark_param_card(process,benchmark,case['w_treatment'])
                        check_production(references[case['w_treatment']],process,identity)
                    exported['production_identity']=str(identity)
                else:
                    references[case['w_treatment']]=process
                record['exports'][key]=exported
            job=run_case(process,case,inputs,args.width_inputs,args.accuracy,sources,args.exclude_split_outliers)
            pairs=verified_pairs(job)
            if streams & pairs:
                raise ValueError('Small-mass controls share actual stage streams')
            streams.update(pairs)
            record['jobs'].append(job)
            record['distinct_stage_pairs']=len(streams)
            save(destination,record)
        record.update(status=DONE,current_case=None,finished_utc=now())
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
