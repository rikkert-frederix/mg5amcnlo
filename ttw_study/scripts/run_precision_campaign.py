#!/usr/bin/env python3
"""Prepare, execute and plot a fresh full-flavour high-statistics S/Pi campaign.

Run after sourcing ttw_study/env.sh. Preparation never starts MG5. Execution
is serialized on 64 cores and can stop at complete flavour-round boundaries.
Execution reads no previous archives and retains final histograms only.
The allocation uses adaptive native integrations with independent retraining;
precision and the known inclusive-normalization flag must be remeasured.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import fcntl
from functools import partial
import gzip
import itertools
import json
import math
from pathlib import Path
import shutil

import numpy as np

from campaign import ROOT, STUDY, RUNTIME_SOURCES, digest, generate, now, run, save
from full_flavour_statistics import FLAVOURS, estimate
from joint_report import normalization_bins
from load_results import assert_matched_physics, assert_same_layout
from pilot_report import clean
from replica_statistics import ratio

BASELINE=STUDY/'inputs/main_campaign_main_onshell_v3.json'
REPORT=STUDY/'results/main_report_main_onshell_v3_coarse_v1.json'
BINNING=STUDY/'inputs/main_binning_main_onshell_v3_coarse_v1.json'
VARIANTS=('S','Pi')
PRIORITY=(3,6,11,12)


def cases(replicas,seed_start):
    if replicas<5 or seed_start<1 or seed_start+32*replicas>900000000:
        raise ValueError('Require at least five retrainings and a valid fresh seed range')
    return [dict(replica=r,charge=c,flavours=list(f),variant=v,w_treatment='onshell',
        decay_bottom_mass=0.,seed=seed_start+i) for i,(r,c,f,v) in enumerate(
            itertools.product(range(replicas),('plus','minus'),FLAVOURS,VARIANTS))]


def allocation(work_factor,replicas,baseline_replicas=5,baseline_accuracy=.015):
    if not math.isfinite(work_factor) or work_factor<1 or replicas<5:
        raise ValueError('Work factor must be >=1 and retrainings >=5')
    per_run=work_factor*baseline_replicas/replicas
    return per_run,baseline_accuracy/math.sqrt(per_run)


def forecast(baseline,work_factor,replicas):
    groups={};inputs={};workers_bytes=0
    for job in baseline['jobs']:
        case=job['case']
        if case['variant'] not in VARIANTS:continue
        path=Path(job['audit']);a=json.loads(path.read_text())
        inputs[str(path)]=digest(path)
        key=(case['charge'],tuple(case['flavours']),case['variant'])
        groups.setdefault(key,[]).append(a['execution']['elapsed_s'])
        workers_bytes+=Path(job['workers']).stat().st_size+Path(job['batches']).stat().st_size
    if len(groups)!=32 or any(len(x)!=5 for x in groups.values()):
        raise ValueError('Timing calibration requires all 160 completed baseline S/Pi runs')
    per_run,accuracy=allocation(work_factor,replicas)
    mean=sum(sum(t) for t in groups.values())*work_factor
    median=sum(float(np.median(t)) for t in groups.values())*per_run*replicas
    return dict(calibration_audit_hashes=inputs,baseline_native_integration_hours=mean/work_factor/3600,
        work_factor=work_factor,per_run_work_factor=per_run,accuracy=accuracy,
        mean_based_native_days=mean/86400,median_based_native_days=median/86400,
        planning_days_with_20_percent_overhead=mean/86400*1.2,
        rough_retained_worker_and_batch_GiB=workers_bytes*work_factor/1024**3,
        limitations='Heuristic T~1/error^2 scaling of actual 64-core wall times, not a completion '
        'deadline or a confidence interval. Rare S runs dominate the mean. Export, training '
        'and archive costs do not scale identically; storage estimate also excludes generated '
        'process files and reductions. Grid adaptation, cancellations and tails can change scaling.')


def shape_forecast(report_path,factor):
    report=json.loads(report_path.read_text());rows=[]
    for title,entry in report['spectra'].items():
        if not title.startswith('R04_b25') or not any(s in title for s in
                ('HTjets','Dphi same-sign','subleading bjet pT','mlb minimax')):continue
        path=Path(entry['path'])
        if digest(path)!=entry['sha256']:raise ValueError('Changed precision calibration spectrum')
        data=json.loads(gzip.decompress(path.read_bytes()));errors=[]
        for row in data['normalized']:
            s,p=[row[v]['scale']['central'] for v in VARIANTS]
            r=row['Pi_over_S']['scale']['central'];e=row['Pi_over_S']['nominal_retraining_mc_error']
            if s is not None and p is not None and min(s,p)>=.005 and r and e is not None:
                errors.append(e/abs(r))
        if errors:rows.append(dict(title=title,populated_bins=len(errors),
            baseline_p90_relative_ratio_error=float(np.quantile(errors,.9)),
            forecast_p90_relative_ratio_error=float(np.quantile(errors,.9))/math.sqrt(factor)))
    return rows


def prepare(args):
    if not args.tag.replace('_','').isalnum():raise ValueError('Invalid tag')
    if args.grid_points<200:raise ValueError('Grid training requires at least 200 points')
    path=STUDY/'inputs'/('precision_campaign_'+args.tag+'.json')
    if path.exists():raise ValueError('Refuse to overwrite campaign')
    baseline=json.loads(BASELINE.read_text())
    if not baseline['status'].startswith('main integrations finished'):raise ValueError('Baseline incomplete')
    inventory=cases(args.replicas,args.seed_start)
    # Check reserved as well as executed seeds in all existing preparation records.
    used=set()
    def seeds(x):
        if isinstance(x,dict):
            for k,v in x.items():
                if k in ('seed','iseed') and isinstance(v,int):used.add(v)
                else:seeds(v)
        elif isinstance(x,list):
            for item in x:seeds(item)
    for p in (STUDY/'inputs').glob('*.json'):seeds(json.loads(p.read_text()))
    if used&{c['seed'] for c in inventory}:raise ValueError('Seed range already used or reserved')
    if not math.isfinite(args.accuracy) or args.accuracy<=0:
        raise ValueError('Accuracy must be positive and finite')
    work_factor=(.015/args.accuracy)**2*args.replicas/5
    timing=forecast(baseline,work_factor,args.replicas)
    timing['accuracy']=args.accuracy
    extra=('run_precision_campaign.py','main_replicas.py','audit_results.py','audit_virtuals.py',
        'harvest_splits.py','read_splits.py','rng_history.py','run_width_mass_pilots.py',
        'full_flavour_statistics.py','replica_statistics.py','main_binning.py','main_figures.py')
    sources={p:digest(ROOT/p) for p in RUNTIME_SOURCES}
    sources.update({'ttw_study/scripts/'+n:digest(STUDY/'scripts'/n) for n in extra})
    references={str(p):digest(p) for p in (BASELINE,REPORT,BINNING,STUDY/'inputs/benchmark.json')}
    save(path,dict(schema='precision_sp_campaign_v1',created_utc=now(),status='prepared',tag=args.tag,
        cases=inventory,jobs=[],current_case=None,source_hashes=sources,reference_hashes=references,
        benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),max_cores=64,
        settings=dict(accuracy=timing['accuracy'],points=1000,grid_points=args.grid_points,
            iterations=3,job_seconds=1.,production_scale='core-w-ht-half',ecm=13000.),
        timing_forecast=timing,shape_forecast=shape_forecast(REPORT,work_factor),
        precision='Exact req_acc_FO=%.6g; %d independent full-flavour run numbers. '
            'Expected %.3gx smaller MC errors if variance scales normally; no promise for every tail. '%
            (args.accuracy,args.replicas,math.sqrt(work_factor))+
            'All signed scales and PDF members retained. No old grids imported or flavour rescaling.',
        known_followup='Inclusive S-P residual and adaptive-stopping coverage remain open. '
            'This campaign measures additional retraining spread; independent fixed-count '
            'normalization controls remain necessary before a publication-precision claim.'))
    print(path)
    print(json.dumps({k:v for k,v in timing.items() if k!='calibration_audit_hashes'},indent=2))


def unchanged(record):
    for p,h in record['source_hashes'].items():
        if digest(ROOT/p)!=h:raise ValueError('Frozen source changed: '+p)
    for p,h in record['reference_hashes'].items():
        if digest(p)!=h:raise ValueError('Frozen reference changed: '+p)


def collect_job(job,**kwargs):
    """Read final histograms only; errors come from independent full runs."""
    from madgraph.various.histograms import HwUList
    from load_results import canonical_columns, WEIGHTS
    histograms=HwUList(job['histogram'],raw_labels=True)
    columns=list(histograms[0].bins[0].wgts)
    selected=canonical_columns(columns,job['variant'])
    values=[];edges=[];offsets=[0];titles=[]
    for hist in histograms:
        data=np.asarray([[row.wgts[columns[i]] for i in selected] for row in hist.bins])
        values.append(data);edges.extend(row.boundaries for row in hist.bins)
        offsets.append(offsets[-1]+len(hist.bins));titles.append(hist.title)
    values=np.concatenate(values)
    layout=dict(titles=titles,edges=np.asarray(edges),offsets=np.asarray(offsets),weights=WEIGHTS)
    return layout,set(),values,np.zeros_like(values),dict(case=job['case'],histogram=job['histogram'])


def execute(path,run_number):
    # The controller lock also makes resuming a stopped invocation safe.
    initial=json.loads(path.read_text())
    with (STUDY/('precision_'+initial['tag']+'.lock')).open('a') as controller_lock:
        fcntl.flock(controller_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        record=json.loads(path.read_text())
        done=len(record['jobs'])
        if [j['case'] for j in record['jobs']]!=record['cases'][:done]:
            raise ValueError('Completed jobs are not an ordered prefix')
        expected=done//32+1
        if run_number!=expected or not 1<=run_number<=len(record['cases'])//32:
            raise ValueError('Run number must be the next uncompleted run: %d'%expected)
        stop=32*run_number
        record.update(status='running',started_utc=record.get('started_utc',now()),
            execution_mode='lean: final histograms only; no automatic archive/RNG/virtual audits')
        save(path,record)
        try:
            for i,case in enumerate(record['cases'][done:stop],start=done):
                process=STUDY/'processes'/('TTW%s_onshell_%s_mb0p0_both_%s'%(
                    case['charge'],''.join(case['flavours']),record['tag']))
                if not process.exists():
                    print('Generating',process.name,flush=True)
                    generate(argparse.Namespace(charge=case['charge'],flavours=case['flavours'],
                        w_treatment='onshell',decay_bottom_mass=0.,corrected='both',tag=record['tag']))
                name='%s_onshell_core-w-ht-half_separate_%d'%(case['variant'],case['seed'])
                cards=process/'study_cards'/name
                histogram=process/'Events'/name/'MADatNLO.HwU'
                execution=cards/'execution.json'
                complete=(execution.exists() and histogram.is_file() and
                    json.loads(execution.read_text()).get('status')=='finished')
                record['current_case']=case;save(path,record)
                if complete:
                    print('Reusing completed integration',name,flush=True)
                else:
                    if cards.exists() or histogram.parent.exists():
                        raise ValueError('Incomplete existing integration; inspect before retrying: '+name)
                    print('Integration %d/%d: %s %s %s seed=%d accuracy=%g'%(i+1,stop,
                        case['charge'],''.join(case['flavours']),case['variant'],case['seed'],
                        record['settings']['accuracy']),flush=True)
                    run(argparse.Namespace(process_dir=str(process),variant=case['variant'],seed=case['seed'],
                        main_run=True,grid_reference=None,**record['settings']))
                    if not histogram.is_file():raise ValueError('Integration produced no final histogram')
                job=dict(case=case,variant=case['variant'],process=str(process),run=name,
                    histogram=str(histogram),execution=str(execution),finished_utc=now())
                record['jobs'].append(job);record['current_case']=None;save(path,record)
            record.update(status='finished' if stop==len(record['cases']) else 'checkpoint',updated_utc=now())
            record.pop('error',None)
        except BaseException as error:
            record.update(status='stopped',error=repr(error),stopped_utc=now())
            save(path,record);raise
        save(path,record)


def analyse(path,output,workers):
    from main_binning import apply_plan
    from main_figures import XLABELS
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    record=json.loads(path.read_text())
    n=len(record['jobs'])
    if record['status'] not in ('checkpoint','finished') or n%32 or n<160:
        raise ValueError('Need at least five complete rounds at a clean checkpoint')
    if output.exists():raise ValueError('Refuse to overwrite plot/reduction directory')
    if not 1<=workers<=64:raise ValueError('Require 1--64 archive readers')
    output.mkdir(parents=True)
    with (STUDY/'mg5.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        collector=partial(collect_job,tag=record['tag'],source_hashes=record['source_hashes'])
        pool=ProcessPoolExecutor(max_workers=workers) if workers>1 else None
        vectors={};conditional={};streams=set();reference=None;provenance=[]
        try:
            results=pool.map(collector,record['jobs']) if pool else map(collector,record['jobs'])
            for projected,pairs,mean,variance,row in results:
                if reference is not None:
                    assert_same_layout([reference,projected])
                else:reference=projected
                if pairs & streams:raise ValueError('Collected runs share actual streams')
                streams.update(pairs);case=row['case'];key=case['charge'],case['variant'],tuple(case['flavours'])
                vectors.setdefault(key,[]).append(mean);conditional.setdefault(key,[]).append(variance[:,0])
                provenance.append(row)
        finally:
            if pool:pool.shutdown(wait=True,cancel_futures=True)
        layout={k:reference[k] for k in ('titles','edges','offsets','weights')}
        vectors={k:np.asarray(v) for k,v in vectors.items()}
        conditional={k:np.asarray(v) for k,v in conditional.items()}
        np.savez_compressed(output/'complete_vectors.npz',**layout,
            **{'__'.join((k[0],k[1],''.join(k[2]))):v for k,v in vectors.items()})
        layout,vectors,_,binning=apply_plan(layout,vectors,conditional,json.loads(BINNING.read_text()))
        summaries={};assets={};arrays={};rate_results={}
        for charge in ('plus','minus'):
            marker='W+' if charge=='plus' else 'W-'
            h=next(i for i,t in enumerate(layout['titles']) if t.startswith('R04_b25 '+marker+' rates:'))
            start=int(layout['offsets'][h]);selected={k:v[:,:,0] for k,v in vectors.items() if k[0]==charge}
            rate_results[charge]={}
            for offset,label in ((0,'input'),(3,'one_b'),(4,'two_b')):
                def rates(s):
                    a,b=s[charge,'S'][start+offset],s[charge,'Pi'][start+offset]
                    return np.array([a,b,b-a,ratio(b,a)])
                r=estimate(selected,charges=(charge,),variants=VARIANTS,transform=rates)
                rate_results[charge][label]={name:dict(value=float(r['value'][i]),mc_error=float(r['mc_error'][i]))
                    for i,name in enumerate(('S_pb','Pi_pb','difference_pb','ratio'))}
            for local in PRIORITY:
                hist=h+local-1;a,b=map(int,layout['offsets'][hist:hist+2]);edges=layout['edges'][a:b]
                parent=[start+i for i in normalization_bins(local)]
                def transform(s):
                    x,y=s[charge,'S'],s[charge,'Pi'];sx=x[a:b]/x[parent].sum();sy=y[a:b]/y[parent].sum()
                    return np.stack((sx,sy,ratio(sy,sx),y[a:b]-x[a:b]))
                r=estimate(selected,charges=(charge,),variants=VARIANTS,transform=transform)
                arrays[charge+'__'+str(local)+'__values']=r['value'];arrays[charge+'__'+str(local)+'__errors']=r['mc_error']
                centers=edges.mean(axis=1);widths=np.diff(edges,axis=1).ravel()
                fig,axes=plt.subplots(3,1,sharex=True,figsize=(6.5,7),layout='constrained')
                for i,label in enumerate(VARIANTS):axes[0].errorbar(centers,r['value'][i]/widths,
                    yerr=r['mc_error'][i]/widths,fmt='o-',markersize=3,label=label)
                axes[1].errorbar(centers,r['value'][2],yerr=r['mc_error'][2],fmt='o',markersize=3)
                axes[1].axhline(1,color='grey',lw=.7);axes[0].legend();axes[0].set_ylabel('Normalized bin density')
                axes[1].set_ylabel('Normalized Pi / S')
                axes[2].errorbar(centers,1000*r['value'][3]/widths,
                    yerr=1000*r['mc_error'][3]/widths,fmt='o',markersize=3)
                axes[2].axhline(0,color='grey',lw=.7)
                axes[2].set_ylabel('Pi - S [fb/'+('rad' if local==6 else 'GeV')+']')
                axes[2].set_xlabel(XLABELS[local])
                fig.suptitle(marker+'; mb=0; all Ws on shell; full direct e/mu sum')
                fig.supxlabel('Independent-retraining MC errors; draft common bins; precision assessment pending',fontsize=8)
                for extension in ('pdf','png'):
                    target=output/('shape_'+charge+'_'+str(local)+'.'+extension);fig.savefig(target);assets[str(target)]=digest(target)
                plt.close(fig)
                summaries[layout['titles'][hist]]=dict(relative_ratio_errors=(r['mc_error'][2]/np.abs(r['value'][2])).tolist())
        np.savez_compressed(output/'shape_reductions.npz',**arrays)
        save(output/'manifest.json',clean(dict(created_utc=now(),campaign=str(path.resolve()),campaign_sha256=digest(path),
            completed_rounds=n//32,planned_rounds=len(record['cases'])//32,provenance=provenance,
            original_vectors_sha256=digest(output/'complete_vectors.npz'),
            reductions_sha256=digest(output/'shape_reductions.npz'),rates=rate_results,shapes=summaries,
            assets=assets,binning=binning,publication_certified=False,
            convention='Average independent retrainings per flavour, sum all eight flavours, then normalize. '
                'Independent S/Pi are not paired. Bin-parent covariance retained. MC bars only; scale/PDF weights '
                'remain in complete_vectors.npz. Failed precision targets require merging or further statistics.')))
    print(output/'manifest.json')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('prepare');p.add_argument('--tag',required=True);p.add_argument('--replicas',type=int,default=10)
    p.add_argument('--accuracy',type=float,default=.005);p.add_argument('--seed-start',type=int,default=900001)
    p.add_argument('--grid-points',type=int,default=1000)
    p=sub.add_parser('execute');p.add_argument('record',type=Path);p.add_argument('--run-number',type=int,required=True,help='Next full-flavour run number (1-based); executes one 32-integration comparison')
    p=sub.add_parser('analyse');p.add_argument('record',type=Path);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--workers',type=int,default=8)
    args=parser.parse_args()
    if args.action=='prepare':prepare(args)
    elif args.action=='execute':execute(args.record,args.run_number)
    else:analyse(args.record,args.output_dir,args.workers)
