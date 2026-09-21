#!/usr/bin/env python3
"""Generated common-scale narrow-W controls, serialized after massive pilots.

Retain full-weight raw batches. Queue completion tests neither the full flavour
sum nor main precision; narrow-limit rates/shapes still need a joint analysis.
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
from narrow_w_inputs import COMMON_SCALE, FACTORS, rescaling
from read_splits import read
from run_inclusive import launch
from run_width_mass_pilots import DONE as PREVIOUS_DONE, verified_pairs

DONE='generated narrow-W pilots finished; inspect rescaled limits, shapes and convergence'


def cases():
    modes=[('onshell',1.),*[(mode,factor) for factor in FACTORS for mode in ('top-bw','all-bw')]]
    return [dict(variant=variant,w_treatment=mode,width_factor=factor,decay_bottom_mass=0.,
                 seed=65001+7*i+j,production_sampling='w-current' if mode=='all-bw' else 'flat')
            for i,variant in enumerate(('LO','S','Pi')) for j,(mode,factor) in enumerate(modes)]


def width_for(case,inputs):
    mode='onshell' if case['w_treatment']=='onshell' else 'bw'
    matching=[row for row in inputs['rows'] if row['width_mode']==mode
              and row['width_factor']==case['width_factor']
              and row['decay_bottom_mass_GeV']==case['decay_bottom_mass']]
    if len(matching)!=1 or matching[0]['mu_reference_GeV']!=COMMON_SCALE:
        raise ValueError('Missing or ambiguous matched common-scale width')
    return matching[0]


def decay_card(case,width):
    from madgraph.fks.fks_decay import decay_card_text
    widths,scales={6:width['gamma_lo']},{6:COMMON_SCALE}
    if case['w_treatment']!='all-bw':
        widths[24],scales[24]=width['w_width_GeV'],80.385
    sampling=case['production_sampling']=='w-current'
    order='LO' if case['variant']=='LO' else 'NLO'
    return decay_card_text(widths,scales,nlo_width_pdgs={6},nlo_widths={6:width['gamma_nlo']},
        production_order=order,decay_order=order,
        nlo_decay_combination='MULTIPLICATIVE' if case['variant']=='Pi' else 'ADDITIVE',
        decay_scale_variation_mode='INDEPENDENT',decay_scale_factors=(1.,.5,2.),
        decay_scale_grouping='SIGNED_PDG',decay_width_scale_modes={6:'AUTO'},
        production_scale_momenta='CORE',production_scale_grouping='NONE',
        production_phase_space_sampling='W_CURRENT' if sampling else 'FLAT',
        production_sampling_mass=80.385 if sampling else None,
        production_sampling_width=width['w_width_GeV'] if sampling else None)


def configure(process,case,inputs,inputs_path,accuracy,grid_points=200):
    from models.check_param_card import ParamCard
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    width=width_for(case,inputs)
    benchmark_param_card(process,benchmark,case['w_treatment'])
    param=ParamCard(str(process/'Cards/param_card.dat'))
    param['decay'].get((24,)).value=width['w_width_GeV']
    param['decay'].get((6,)).value=width['gamma_nlo']
    param.write(str(process/'Cards/param_card.dat'),precision=16)
    setup_module().configure(argparse.Namespace(process_dir=str(process),variant=case['variant'],
        top_width_lo=width['gamma_lo'],top_width_nlo=width['gamma_nlo'],w_treatment=case['w_treatment'],
        top_width_w_treatment=width['width_mode'],decay_bottom_mass=case['decay_bottom_mass'],
        top_width_bottom_mass=case['decay_bottom_mass'],width_source=str(inputs_path)+'; sha256='+digest(inputs_path),
        pdf_id=331700,production_scale='fixed',production_sampling=case['production_sampling'],
        decay_scales='separate',ecm=13000.,seed=case['seed'],points=1000,grid_points=grid_points,
        iterations=3,accuracy=accuracy,job_seconds=1.))
    base_name='%s_%s_fixed_separate_%d' % (case['variant'],case['w_treatment'],case['seed'])
    if case['production_sampling']=='w-current':
        base_name+='_wcurrent'
    base=process/'study_cards'/base_name
    manifest=copy.deepcopy(json.loads((base/'manifest.json').read_text()))
    name=base_name+'__narrow_'+('%g' % case['width_factor']).replace('.','p')+'_common'
    archive=process/'study_cards'/name
    if archive.exists() or (process/'Events'/name).exists():
        raise ValueError('Existing generated narrow-W control')
    archive.mkdir()
    # The standard setup snapshot remains an unlaunched base. This new
    # explicit archive changes both decay numerator and width-reference scales.
    (process/'Cards/decay_card.dat').write_text(decay_card(case,width))
    for filename in ('param_card.dat','decay_card.dat','run_card.dat','FO_analyse_card.dat'):
        shutil.copy2(process/'Cards'/filename,archive/filename)
        manifest['hashes'][filename]=digest(archive/filename)
    manifest.update(run_name=name,top_width_reference_scale=COMMON_SCALE,narrow_W_test=case,
        narrow_W_inputs=str(inputs_path.resolve()),narrow_W_inputs_sha256=digest(inputs_path),
        physical_scale_coordinates=dict(common_fixed_scale_GeV=COMMON_SCALE,
                                        weights_relative_to_common_scale=True,top_width_mode='AUTO'),
        comparison_cross_section_rescaling=rescaling(case['width_factor']),
        generated_event_weight_rescaling=1.,base_configuration=str(base/'manifest.json'),
        base_configuration_sha256=digest(base/'manifest.json'),
        status='configured common-scale narrow-W pilot; no integration performed')
    save(archive/'manifest.json',manifest)
    return name,archive


def run_case(process,case,inputs,inputs_path,accuracy,sources,grid_points=200):
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        name,archive=configure(process,case,inputs,inputs_path,accuracy,grid_points)
        execution=dict(max_cores=64,pilot=True,narrow_W_test=case,source_hashes=sources,
            export_generation=json.loads((STUDY/'inputs'/(process.name+'_generation.json')).read_text()),
            requested_total_rate_accuracy=accuracy,grid_reference=None)
        launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
                'NLO','-f','-n',name,'--multicore','--nb_core=64'],process,
               STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
        hwu=process/'Events'/name/'MADatNLO.HwU'
        if not hwu.is_file():
            raise ValueError('No completed narrow-W histogram')
        execution['outputs']={str(hwu):digest(hwu)}
        save(archive/'execution.json',execution)
        audit(hwu)
        audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
        workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
        harvest(audit_path,workers)
    return dict(process=str(process),run=name,case=case,variant=case['variant'],
                audit=str(audit_path),workers=str(workers))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--width-inputs',required=True,type=Path)
    parser.add_argument('--lo-accuracy',type=float,default=.005)
    parser.add_argument('--nlo-accuracy',type=float,default=.03)
    args=parser.parse_args()
    if not args.tag.replace('_','').isalnum() or any(not 0.<v<1. for v in (args.lo_accuracy,args.nlo_accuracy)):
        raise ValueError('Invalid narrow-W pilot tag/accuracy')
    destination=STUDY/'inputs'/('narrow_w_pilot_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Refuse to overwrite a narrow-W queue')
    inputs=json.loads(args.width_inputs.read_text())
    if (inputs['benchmark_sha256']!=digest(STUDY/'inputs/benchmark.json')
            or inputs['script_sha256']!=digest(STUDY/'scripts/narrow_w_inputs.py')
            or inputs['common_scale_GeV']!=COMMON_SCALE or inputs['width_factors']!=list(FACTORS)):
        raise ValueError('Narrow-W input convention/evidence changed')
    files=[Path(__file__).resolve(),args.width_inputs.resolve(),STUDY/'inputs/benchmark.json',
           *[STUDY/'scripts'/name for name in ('narrow_w_inputs.py','run_inclusive.py',
             'run_width_mass_pilots.py','check_phase_space_support.py','audit_results.py',
             'audit_virtuals.py','harvest_splits.py','read_splits.py','load_results.py',
             'rng_history.py','compare_current_batches.py')],
           *[ROOT/'tests/input_files/fks_decay'/name for name in (
             'phase_space_test_dimensions.f90','bw_support_checks.f90')],
           ROOT/'Template/fNLO/Source/kin_functions.f90']
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in files})
    record=dict(created_utc=now(),status='waiting for massive and W-scheme companions',max_cores=64,
                source_hashes=sources,after_queue=str(args.after_queue.resolve()),width_inputs=str(args.width_inputs.resolve()),
                cases=cases(),jobs=[],exports={},charge='plus',flavours=['e','e','mu'],pilot=True,
                limitations='No full flavour sum or main precision. Rescaled generated rate/shape narrowing '
                            'must be inspected; deterministic width continuity alone does not establish it.')
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen narrow-W pilot source changed')
    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                raise ValueError('Preceding massive companion queue stopped')
            if previous['status']==PREVIOUS_DONE:
                break
            time.sleep(15)
        unchanged()
        if len(previous['jobs'])!=16 or any(digest(ROOT/p)!=sha for p,sha in previous['source_hashes'].items()):
            raise ValueError('Incomplete or changed massive/W-scheme predecessor')
        record['predecessor_sha256']=digest(args.after_queue)
        controls_path=Path(previous['after_queue'])
        if digest(controls_path)!=previous['predecessor_sha256']:
            raise ValueError('Original stage-audited control record changed')
        controls=json.loads(controls_path.read_text())
        seen=set()
        for job in controls['jobs']+previous['jobs']:
            pairs=verified_pairs(job)
            if seen & pairs:
                raise ValueError('Overlapping predecessor training/refinement streams')
            seen.update(pairs)
        for case in cases():
            unchanged()
            mode=case['w_treatment']
            process=STUDY/'processes'/('TTWplus_%s_eemu_mb0p0_both_%s' % (mode,args.tag))
            record.update(status='preparing generated narrow-W pilot',current_case=case)
            save(destination,record)
            if mode not in record['exports']:
                generate(argparse.Namespace(charge='plus',w_treatment=mode,flavours=['e','e','mu'],
                                             corrected='both',decay_bottom_mass=0.,tag=args.tag))
                support=STUDY/'inputs'/('%s_%s_phase_space_support.json' % (args.tag,mode))
                check_support(support,process)
                record['exports'][mode]=dict(process=str(process),support=str(support))
            unchanged()
            record['status']='running generated narrow-W pilot'
            save(destination,record)
            job=run_case(process,case,inputs,args.width_inputs,
                         args.lo_accuracy if case['variant']=='LO' else args.nlo_accuracy,sources)
            if case['variant']!='LO':
                virtuals=STUDY/'results'/(args.tag+'_'+job['run']+'_virtual_checks.json')
                audit_virtuals(Path(job['workers']).with_suffix('.json'),virtuals)
                job['analytic_virtual_checks']=str(virtuals)
            batches=STUDY/'results'/(args.tag+'_'+job['run']+'_batches.npz')
            read(Path(job['workers']).with_suffix('.json'),batches)
            job['batches']=str(batches)
            pairs=verified_pairs(job)
            if seen & pairs:
                raise ValueError('Narrow-W control overlaps another pilot stream')
            seen.update(pairs)
            job['finished_utc']=now()
            record['jobs'].append(job)
            save(destination,record)
        record.update(status=DONE,finished_utc=now(),current_case=None)
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
