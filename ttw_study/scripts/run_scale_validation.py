#!/usr/bin/env python3
"""Serialized direct scale/QES pilots and central-production-scale diagnostics.

All use fresh technical exports, the frozen main physical inputs and all
81 scale/101 PDF weights. A user hook changes only the signed top numerator
scales in the dedicated asymmetric export; no main/default source is edited.
Every modified configuration gets a new archive, distinct from setup's
unchanged base snapshot. Scientific agreement is never inferred from exit 0.
"""
import argparse
import copy
import fcntl
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

from audit_results import audit
from audit_virtuals import audit as audit_virtuals
from campaign import (ROOT, STUDY, RUNTIME_SOURCES, benchmark_param_card,
                      digest, generate, now, save, setup_module)
from harvest_splits import harvest
from joint_comparison import paired_variants, read as load_batches
from load_results import POINTS, assert_same_layout
from pilot_report import RATES, clean
from read_splits import read as read_batches
from replica_statistics import ratio
from run_inclusive import launch
from run_width_mass_pilots import verified_pairs

HOOK = STUDY/'assets/asymmetric_decay_scale.f90'
DONE = 'scale validation pilots finished; inspect identities and convergence'


def reference_controls(path):
    """Select completed controls with complete actual refinement seed history."""
    from run_rng_pilots import DONE as CONTROLS_DONE
    queue = json.loads(path.read_text())
    if queue['status'] != CONTROLS_DONE or len(queue['jobs']) != 4:
        raise ValueError('Require four completed RNG-corrected on-shell/all-BW controls')
    references, streams = {}, set()
    for job in queue['jobs']:
        audit_record = json.loads(Path(job['audit']).read_text())
        manifest = audit_record['manifest']
        args = audit_record['execution']['export_generation']['args']
        key = manifest['w_treatment'], job['variant']
        if (key in references or key[0] not in ('onshell', 'all-bw') or key[1] not in ('S', 'Pi')
                or args['charge'] != 'plus' or args['flavours'] != ['e', 'e', 'mu']
                or args['decay_bottom_mass'] != 0. or args['corrected'] != 'both'
                or manifest['production_scale'] != 'core-w-ht-half'):
            raise ValueError('Unexpected technical reference physics')
        pairs = verified_pairs(job)
        if pairs & streams:
            raise ValueError('Technical reference controls reuse actual random streams')
        streams.update(pairs)
        references[key] = Path(job['batches'])
    return references, streams


def cases():
    return [dict(name='asymmetric_t2_at_half', mode='onshell', group='asymmetric',
                 production_scale='core-w-ht-half', factors=[1.,1.,2.,.5], qes=1.),
            *[dict(name='qes_'+label, mode='onshell', group='qes',
                   production_scale='core-w-ht-half', factors=[1.]*4, qes=value)
              for label,value in (('half',.5),('two',2.))],
            dict(name='production_r2_f_half', mode='all-bw', group='bw',
                 production_scale='core-w-ht-half', factors=[2.,.5,1.,1.], qes=1.),
            *[dict(name='central_'+label, mode='all-bw', group='bw',
                   production_scale=scale, factors=[1.]*4, qes=1.)
              for label,scale in (('native','core-ht-half'),('fixed','fixed'))]]


def physical_scale_coordinates(case):
    return dict(production_scale=case['production_scale'],
                numerator_factors=case['factors'], qes_over_ref=case['qes'],
                top_width_coefficient_reference_GeV=172.5,
                decay_numerator_scales_GeV=[172.5*f for f in case['factors'][2:]],
                weights_relative_to_these_numerator_scales=True)


def decay_card(variant, benchmark, mode, asymmetric=False):
    from madgraph.fks.fks_decay import decay_card_text
    widths = next(row for row in benchmark['top_widths'] if row['mb_GeV']==0.
                  and row['scale_factor']==1.
                  and row['w_treatment']==('onshell' if mode=='onshell' else 'bw'))
    lo, references = {6:widths['gamma_lo']}, {6:172.5}
    if mode != 'all-bw':
        lo[24],references[24] = benchmark['W_width']['gamma_nlo_pdf'],80.385
    return decay_card_text(
        lo,references,nlo_width_pdgs={6},nlo_widths={6:widths['gamma_nlo_pdf']},
        production_order='NLO',decay_order='NLO',
        nlo_decay_combination='MULTIPLICATIVE' if variant=='Pi' else 'ADDITIVE',
        decay_scale_variation_mode='INDEPENDENT',decay_scale_factors=(1.,.5,2.),
        decay_scale_grouping='SIGNED_PDG',decay_width_scale_modes={6:'AUTO'},
        decay_dynamical_scale_choices={6:-1} if asymmetric else {},
        production_scale_momenta='CORE',production_scale_grouping='W_SYSTEM')


def install_hook(process):
    """Configure the documented user hook in a new, never-compiled export."""
    target = process/'SubProcesses/dummy_fct.f90'
    default = ROOT/'Template/fNLO/SubProcesses/dummy_fct.f90'
    if target.is_symlink() or digest(target)!=digest(default):
        raise ValueError('Fresh hook must be a private copy of the default')
    event_entries=[p for p in (process/'Events').iterdir()
                   if not (p.name=='.keep' and p.is_file() and not p.is_symlink() and p.stat().st_size==0)]
    if any(process.glob('SubProcesses/P*/dummy_fct.o')) or event_entries:
        raise ValueError('Do not install a hook into a used/compiled export')
    archive = process/'study_cards/asymmetric_user_hook'
    archive.mkdir(parents=True)
    shutil.copy2(target,archive/'default_dummy_fct.f90')
    shutil.copy2(HOOK,target)
    shutil.copy2(HOOK,archive/'asymmetric_decay_scale.f90')
    record = dict(source=str(HOOK),sha256=digest(HOOK),
                  default_sha256=digest(default),installed_before_compilation=True,
                  interface='fixed_user_decay_scale, choice -1 for abs(PDG)=6',
                  changed_generator_or_main_default=False)
    save(archive/'installation.json',record)
    return record


def overlapping_points(factors):
    """Relative direct-run coordinates -> absolute reference coordinates."""
    pairs = [(0,1+POINTS.index(tuple(factors)),list(factors))]
    for index,point in enumerate(POINTS,1):
        absolute = tuple(a*b for a,b in zip(point,factors))
        if absolute in POINTS:
            pairs.append((index,1+POINTS.index(absolute),list(absolute)))
    return pairs


def assert_technical_match(direct, reference, case):
    assert_same_layout([direct,reference])
    a,b = [row['report']['manifest'] for row in (direct,reference)]
    for key in ('variant','w_treatment','decay_bottom_mass','top_width_lo',
                'top_width_nlo','w_width','top_width_reference_scale',
                'top_width_w_treatment','top_width_bottom_mass','decay_scale_grouping'):
        if a[key]!=b[key]:
            raise ValueError('Technical comparison physical mismatch: '+key)
    for row in (direct,reference):
        args = row['report']['execution']['export_generation']['args']
        if args['charge']!='plus' or args['flavours']!=['e','e','mu'] or args['corrected']!='both':
            raise ValueError('Technical comparison requires the same ordered W+ e,e,mu process')
    if direct['parameter_signature']!=reference['parameter_signature']:
        raise ValueError('Technical comparison parameter mismatch')
    if a['hashes']['analysis_source']!=b['hashes']['analysis_source']:
        raise ValueError('Technical comparison analysis mismatch')
    ignored = {'iseed','req_acc_fo','npoints_fo','niters_fo','npoints_fo_grid',
               'niters_fo_grid','fo_job_target_time','mur_over_ref','muf_over_ref','qes_over_ref'}
    physics = lambda m: {k:v for k,v in m['settings'].items() if k not in ignored}
    if physics(a)!=physics(b) or a['production_scale']!=b['production_scale']:
        raise ValueError('Unexpected run-setting/central-scale change in identity test')
    for key,value in zip(('mur_over_ref','muf_over_ref'),case['factors'][:2]):
        if a['settings'][key]!=value or b['settings'].get(key,1.)!=1.:
            raise ValueError('Unmatched production reference factor')
    if a['settings']['qes_over_ref']!=case['qes'] or b['settings'].get('qes_over_ref',1.)!=1.:
        raise ValueError('Unmatched QES reference factor')
    if a['physical_scale_coordinates']!=physical_scale_coordinates(case) or b.get('technical_test'):
        raise ValueError('Wrong technical scale convention or nonstandard reference')
    hook=a.get('custom_user_hook')
    if case['group']=='asymmetric':
        if not hook or hook['sha256']!=digest(HOOK):
            raise ValueError('Missing verified asymmetric user hook')
    elif hook:
        raise ValueError('Unexpected user hook in a production/QES test')
    if b.get('custom_user_hook') or b.get('physical_scale_coordinates'):
        raise ValueError('Reference must use the unmodified central scale convention')
    if set(map(tuple,direct['seeds'])) & set(map(tuple,reference['seeds'])):
        raise ValueError('Random-stream collision between independently sampled tests')


def compare(direct_path, reference_path, case, destination):
    direct,reference = load_batches(direct_path),load_batches(reference_path)
    assert_technical_match(direct,reference,case)
    pairs = overlapping_points(case['factors'])
    di,ri = ([p[k] for p in pairs] for k in (0,1))
    report = dict(created_utc=now(),case=case,
                  status='direct/reweighted identity evaluated; agreement not automatically certified',
                  inputs={str(p):digest(p) for p in (direct_path,reference_path)},
                  reference_points=[p[2] for p in pairs],
                  first_point='direct nominal versus the requested reference reweight',
                  convention='Independent runs, conditional split-batch MC errors; '
                             'correlated per-bin pulls are not a global chi-squared test.',
                  rates={},spectra={})
    def estimate(a,b,normalized=False):
        def transform(x,y):
            if normalized:
                x,y = ratio(x[:-1],x[-1]),ratio(y[:-1],y[-1])
            return np.asarray([x,y,x-y,ratio(x,y)])
        return paired_variants(a,b,direct['strata'],reference['strata'],transform)
    h = next(i for i,title in enumerate(direct['titles']) if title.startswith('R04_b25')
             and 'W+' in title and ' rates:' in title)
    start = int(direct['offsets'][h])
    rate_inputs = [row['contributions'][:,start:start+9,:][:,:,indices]
                   for row,indices in ((direct,di),(reference,ri))]
    rates = estimate(*rate_inputs)
    for i,label in enumerate(RATES):
        report['rates'][label] = dict(direct_pb=rates['value'][0,i],
                                     reference_pb=rates['value'][1,i],
                                     difference_pb=rates['value'][2,i],
                                     difference_mc_error_pb=rates['mc_error'][2,i],
                                     difference_pull=ratio(rates['value'][2,i],rates['mc_error'][2,i]))
    arrays = dict(rate_values=rates['value'],rate_mc_errors=rates['mc_error'])
    for local in (3,4,6,11,12,20):
        lo,hi = map(int,direct['offsets'][h+local-1:h+local+1])
        parent = start+(4 if local in (11,12) else 3)
        inputs = [np.concatenate([row['contributions'][:,lo:hi,:][:,:,indices],
                                  row['contributions'][:,parent:parent+1,:][:,:,indices]],axis=1)
                  for row,indices in ((direct,di),(reference,ri))]
        result = estimate(*inputs,normalized=True)
        pull = ratio(result['value'][2,:,0],result['mc_error'][2,:,0])
        finite = np.isfinite(pull)
        report['spectra'][direct['titles'][h+local-1]] = dict(
            finite_nominal_pulls=int(finite.sum()),
            maximum_absolute_nominal_pull=float(np.abs(pull[finite]).max()) if finite.any() else None,
            nominal_bins_above_three=int(np.sum(np.abs(pull[finite])>3.)),
            normalized_to='fiducial_2b' if local in (11,12) else 'fiducial_1b')
        arrays['normalized_%d'%local] = result['value']
        arrays['normalized_mc_errors_%d'%local] = result['mc_error']
        arrays['edges_%d'%local] = direct['edges'][lo:hi]
    if destination.exists() or destination.with_suffix('.npz').exists():
        raise ValueError('Existing technical comparison')
    np.savez_compressed(destination.with_suffix('.npz'),**arrays)
    report['arrays_sha256'] = digest(destination.with_suffix('.npz'))
    save(destination,clean(report))


def run_case(process, case, variant, seed, accuracy, sources, hook):
    from madgraph.various.banner import RunCardNLO
    benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
    setup = setup_module()
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        mb,width = benchmark_param_card(process,benchmark,case['mode'])
        setup.configure(argparse.Namespace(
            process_dir=str(process),variant=variant,top_width_lo=width['gamma_lo'],
            top_width_nlo=width['gamma_nlo_pdf'],w_treatment=case['mode'],
            top_width_w_treatment=width['w_treatment'],decay_bottom_mass=mb,
            top_width_bottom_mass=mb,width_source=str(STUDY/'inputs/benchmark.json')+
            '; sha256='+digest(STUDY/'inputs/benchmark.json'),pdf_id=331700,
            production_scale=case['production_scale'],decay_scales='separate',
            production_sampling='w-current' if case['mode']=='all-bw' else 'flat',
            ecm=13000.,seed=seed,points=1000,grid_points=200,iterations=3,
            accuracy=accuracy,job_seconds=1.))
        base_name='%s_%s_%s_separate_%d'%(variant,case['mode'],case['production_scale'],seed)
        if case['mode']=='all-bw':
            base_name += '_wcurrent'
        base=process/'study_cards'/base_name
        manifest=copy.deepcopy(json.loads((base/'manifest.json').read_text()))
        name=base_name+'__'+case['name']
        archive=process/'study_cards'/name
        if archive.exists() or (process/'Events'/name).exists():
            raise ValueError('Existing technical run')
        archive.mkdir()
        settings=manifest['settings']
        settings.update(mur_over_ref=case['factors'][0],muf_over_ref=case['factors'][1],
                        qes_over_ref=case['qes'])
        card=RunCardNLO(str(process/'Cards/run_card.dat'))
        for key,value in settings.items():
            card.set(key,value,user=True,raiseerror=True)
        card.write(str(process/'Cards/run_card.dat'),template=str(process/'Cards/run_card_default.dat'))
        if hook:
            if digest(process/'SubProcesses/dummy_fct.f90')!=hook['sha256']:
                raise ValueError('Dedicated asymmetric hook changed')
            (process/'Cards/decay_card.dat').write_text(decay_card(variant,benchmark,case['mode'],True))
            shutil.copy2(HOOK,archive/HOOK.name)
        for filename in ('param_card.dat','decay_card.dat','run_card.dat','FO_analyse_card.dat'):
            shutil.copy2(process/'Cards'/filename,archive/filename)
            manifest['hashes'][filename]=digest(archive/filename)
        manifest.update(run_name=name,technical_test=case,physical_scale_coordinates=physical_scale_coordinates(case),
                        custom_user_hook=hook,base_configuration=str(base/'manifest.json'),
                        base_configuration_sha256=digest(base/'manifest.json'),
                        status='configured technical pilot; no integration performed')
        save(archive/'manifest.json',manifest)
        execution=dict(max_cores=64,pilot=True,technical_test=case,source_hashes=sources,
                       export_generation=json.loads((STUDY/'inputs'/(process.name+'_generation.json')).read_text()),
                       custom_user_hook=hook,requested_total_rate_accuracy=accuracy)
        launch(['taskset','-c','0-63',sys.executable,str(process/'bin/calculate_xsect'),
                'NLO','-f','-n',name,'--multicore','--nb_core=64'],process,
               STUDY/'logs'/(process.name+'_'+name+'.log'),archive/'execution.json',execution)
        hwu=process/'Events'/name/'MADatNLO.HwU'
        if not hwu.is_file():
            raise ValueError('No final technical histogram')
        execution['outputs']={str(hwu):digest(hwu)}
        save(archive/'execution.json',execution)
        audit(hwu)
        audit_path=STUDY/'results/audits'/(process.name+'_'+name+'.json')
        workers=STUDY/'results/workers'/(process.name+'_'+name+'.tar.gz')
        harvest(audit_path,workers)
    return dict(process=str(process),run=name,audit=str(audit_path),workers=str(workers))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',required=True)
    parser.add_argument('--after-queue',required=True,type=Path)
    parser.add_argument('--validation-queue',required=True,type=Path)
    parser.add_argument('--accuracy',type=float,default=.03)
    args=parser.parse_args()
    if not args.tag.replace('_','').isalnum() or not 0.<args.accuracy<1.:
        raise ValueError('Invalid technical tag/accuracy')
    destination=STUDY/'inputs'/('scale_validation_queue_'+args.tag+'.json')
    if destination.exists():
        raise ValueError('Existing technical queue; inspect before resuming')
    extra=[Path(__file__).resolve(),HOOK,STUDY/'inputs/benchmark.json',
           *[STUDY/'scripts'/name for name in ('run_inclusive.py','audit_results.py',
             'harvest_splits.py','read_splits.py','load_results.py','joint_comparison.py',
             'replica_statistics.py','audit_virtuals.py','run_width_mass_pilots.py')],
           *[ROOT/'Template/fNLO/SubProcesses'/name for name in ('dummy_fct.f90',
             'setscales.f90','decay_chain_parameters.f90','decay_chain_scales.f90',
             'fks_weights.f90','spin_density_weight_lines.f90')]]
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({str(p.relative_to(ROOT)):digest(p) for p in extra})
    record=dict(created_utc=now(),status='waiting for ttbar reference pilots',max_cores=64,
                source_hashes=sources,after_queue=str(args.after_queue.resolve()),
                validation_queue=str(args.validation_queue.resolve()),jobs=[],cases=cases(),
                pilot=True,charge='plus',flavours=['e','e','mu'])
    save(destination,record)
    def unchanged():
        if any(digest(ROOT/p)!=sha for p,sha in sources.items()):
            raise ValueError('Frozen technical source/input changed')
    try:
        while True:
            previous=json.loads(args.after_queue.read_text())
            if previous['status']=='stopped':
                raise ValueError('Preceding reference queue stopped')
            if previous['status']=='ttbar reference pilot finished; inspect comparison and convergence':
                break
            time.sleep(15)
        unchanged()
        reference_batches, streams = reference_controls(args.validation_queue)
        record['reference_queues_sha256']={str(p):digest(p) for p in (args.after_queue,args.validation_queue)}
        exports={}
        for i,case in enumerate(cases()):
            unchanged()
            if case['group'] not in exports:
                tag=args.tag+'_'+case['group']
                process=STUDY/'processes'/('TTWplus_%s_eemu_mb0p0_both_%s'%(case['mode'],tag))
                record.update(status='generating technical export',current_case=case)
                save(destination,record)
                generate(argparse.Namespace(charge='plus',flavours=['e','e','mu'],
                         w_treatment=case['mode'],decay_bottom_mass=0.,corrected='both',tag=tag))
                hook=None
                if case['group']=='asymmetric':
                    with (STUDY/'mg5.lock').open('a') as lock:
                        fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
                        hook=install_hook(process)
                exports[case['group']]=(process,hook)
            process,hook=exports[case['group']]
            for j,variant in enumerate(('S','Pi')):
                unchanged()
                record.update(status='running technical pilot',current_case=case,current_variant=variant)
                save(destination,record)
                job=run_case(process,case,variant,71001+100*i+j,args.accuracy,sources,hook)
                job['variant'] = variant
                record['jobs'].append(job)
                save(destination,record)
                batches=STUDY/'results'/('%s_%s_%s_batches.npz'%(args.tag,case['name'],variant))
                read_batches(Path(job['workers']).with_suffix('.json'),batches)
                job['batches'] = str(batches)
                pairs = verified_pairs(job)
                if pairs & streams:
                    raise ValueError('Direct-scale run shares actual streams with another control')
                streams.update(pairs)
                virtuals=STUDY/'results'/('%s_%s_%s_virtuals.json'%(args.tag,case['name'],variant))
                audit_virtuals(Path(job['workers']).with_suffix('.json'),virtuals)
                job['analytic_virtual_checks'] = str(virtuals)
                if not case['name'].startswith('central_'):
                    target=STUDY/'results'/('%s_%s_%s_comparison.json'%(args.tag,case['name'],variant))
                    compare(batches,reference_batches[case['mode'],variant],case,target)
                    job.update(batches=str(batches),comparison=str(target))
                save(destination,record)
        record.update(status=DONE,finished_utc=now(),current_case=None,current_variant=None)
    except BaseException as error:
        record.update(status='stopped',error=repr(error),stopped_utc=now())
        save(destination,record)
        raise
    save(destination,record)


if __name__=='__main__':
    main()
