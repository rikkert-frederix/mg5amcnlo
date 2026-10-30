#!/usr/bin/env python3
"""Inspect complete main ensembles without discarding or rescaling any run.

Independent linear contrasts, inclusive identities and paired decay-scale
cancellations complement the nonlinear full-flavour report. The few-run
variance estimates and multiple correlated comparisons are diagnostics,
not normal-distribution significances or a publication release.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from full_flavour_statistics import FLAVOURS, VARIANTS, estimate
from joint_report import normalization_bins
from load_results import POINTS, WEIGHTS
from main_observables import LABELS, normalize_spectrum, prescription_observables
from main_report import feature_vectors, read_main, shape_precision
from pilot_report import RATES, clean
from replica_statistics import ratio

INCLUSIVE_CONTRASTS=(('D_minus_LO',{'D':1.,'LO':-1.}),
    ('PiD_minus_LO',{'PiD':1.,'LO':-1.}),('S_minus_P',{'S':1.,'P':-1.}),
    ('Pi_minus_P',{'Pi':1.,'P':-1.}),('Pi_minus_S',{'Pi':1.,'S':-1.}))


def linear_summary(ensembles, coefficients=None):
    """Independent mean sum and Welch--Satterthwaite variance diagnostic."""
    rows=[np.asarray(row,dtype=float) for row in ensembles]
    coefficients=np.ones(len(rows)) if coefficients is None else np.asarray(coefficients,dtype=float)
    if (not rows or coefficients.shape!=(len(rows),) or not np.isfinite(coefficients).all()
            or any(row.ndim<2 or len(row)<5 or not np.isfinite(row).all() for row in rows)
            or len({row.shape[1:] for row in rows})!=1):
        raise ValueError('Require finite matched independent ensembles with at least five runs')
    means=[c*row.mean(axis=0) for c,row in zip(coefficients,rows)]
    variances=[c*c*np.var(row,axis=0,ddof=1)/len(row) for c,row in zip(coefficients,rows)]
    value=sum(means)
    variance=sum(variances)
    effective_dof=ratio(variance**2,sum(v*v/(len(row)-1) for v,row in zip(variances,rows)))
    return dict(value=value,mc_error=np.sqrt(variance),effective_dof=effective_dof,
                max_ensemble_variance_fraction=ratio(np.max(variances,axis=0),variance))


def contrast(vectors,charge,coefficients,selection):
    rows,weights=[],[]
    for variant,coefficient in coefficients.items():
        for flavour in FLAVOURS:
            rows.append(vectors[charge,variant,flavour][:,selection])
            weights.append(coefficient)
    return linear_summary(rows,weights)


def diagnostic(result,labels):
    value,error=result['value'],result['mc_error']
    pull=ratio(value,error)
    finite=np.isfinite(pull)
    maximum=int(np.argmax(np.where(finite,np.abs(pull),-1.))) if finite.any() else None
    return dict(tested_nonzero_error_entries=int(finite.sum()),
        entries_above_three_mc_errors=int(np.sum(finite&(np.abs(pull)>3.))),
        maximum_absolute_pull=float(np.abs(pull.flat[maximum])) if maximum is not None else None,
        maximum_label=labels[maximum] if maximum is not None else None,
        value=value,mc_error=error,pull=pull,effective_dof=result['effective_dof'],
        interpretation='Ratios to estimated MC errors; correlated comparisons and finite retraining '
                       'counts preclude a global significance interpretation.')


def run(path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite convergence evidence')
    source,layout,vectors,conditional=read_main(path)
    report=dict(created_utc=now(),status='main convergence inspected; numerical flags require scientific assessment',
        input=str(path.resolve()),input_sha256=digest(path),input_arrays_sha256=source['arrays_sha256'],
        benchmark_sha256=source['benchmark_sha256'],runs=len(source['samples']),
        distinct_stage_initializations=source['distinct_stage_initializations'],
        cross_run_rng_overlap=source['cross_run_rng_overlap'],rates={},inclusive={},inclusive_charge_sum={},
        inclusive_decay_scales={},strict_identity={},retraining_groups={},original_shape_precision={},
        no_runs_discarded=True,publication_certified=False,
        limitations='Five independent runs per ordered flavour determine the errors; they do not '
            'establish rare-tail coverage or absence of adaptive-stopping bias. Effective degrees '
            'of freedom assume independent finite-variance approximately Gaussian run estimates. '
            'Neither conditional variances nor a fitted variance inflation are added to run scatter. '
            'A local residual above three estimated errors is an inspection flag, not a global test.')
    arrays={}
    input_starts={}
    nominal_vectors={key:value[:,:,:1] for key,value in vectors.items()}
    for charge in ('plus','minus'):
        marker='W+' if charge=='plus' else 'W-'
        rate_h=[h for h,t in enumerate(layout['titles']) if marker in t and ' rates:' in t]
        nominal=next(h for h in rate_h if layout['titles'][h].startswith('R04_b25'))
        start=int(layout['offsets'][nominal])
        input_starts[charge]=start
        for h in rate_h:
            a,b=map(int,layout['offsets'][h:h+2])
            bank={}
            for variant in VARIANTS:
                result=contrast(vectors,charge,{variant:1.},slice(a,b))
                bank[variant]={label:dict(value_pb=float(result['value'][j,0]),
                    mc_error_pb=float(result['mc_error'][j,0]),
                    relative_mc_error=float(ratio(result['mc_error'][j,0],result['value'][j,0])),
                    effective_dof=float(result['effective_dof'][j,0]),
                    max_flavour_variance_fraction=float(result['max_ensemble_variance_fraction'][j,0]))
                    for j,label in enumerate(RATES)}
            report['rates'][layout['titles'][h]]=bank
        for name,coefficients in INCLUSIVE_CONTRASTS:
            result=contrast(vectors,charge,coefficients,start)
            report['inclusive'][charge+'__'+name]=diagnostic(result,WEIGHTS)
        for variant in ('S','Pi','D','PiD'):
            points=[point for point in POINTS if point[2:]!=(1.,1.)]
            indices=[1+POINTS.index(point) for point in points]
            references=[1+POINTS.index(point[:2]+(1.,1.)) for point in points]
            differences=[vectors[charge,variant,f][:,start,indices]-
                         vectors[charge,variant,f][:,start,references] for f in FLAVOURS]
            report['inclusive_decay_scales'][charge+'__'+variant]=diagnostic(
                linear_summary(differences),[WEIGHTS[index] for index in indices])
        active=np.concatenate([np.arange(layout['offsets'][h],layout['offsets'][h+1])
            for h,t in enumerate(layout['titles']) if marker in t])
        result=contrast(vectors,charge,{'S':1.,'P':-1.,'D':-1.,'LO':1.},active)
        for name in ('value','mc_error','effective_dof'):
            arrays['strict_identity_'+charge+'_'+name]=result[name]
        arrays['strict_identity_'+charge+'_original_bin_indices']=active
        pull=ratio(result['value'],result['mc_error'])
        finite=np.isfinite(pull)
        nominal_pull=pull[:,0]
        valid=np.isfinite(nominal_pull)
        maximum=int(np.argmax(np.where(valid,np.abs(nominal_pull),-1.)))
        original=int(active[maximum]); h=int(np.searchsorted(layout['offsets'],original,side='right')-1)
        report['strict_identity'][charge]=dict(nominal_nonzero_error_bins=int(valid.sum()),
            nominal_bins_above_three_mc_errors=int(np.sum(valid&(np.abs(nominal_pull)>3.))),
            maximum_nominal_absolute_pull=float(abs(nominal_pull[maximum])),
            maximum_nominal_title=layout['titles'][h],maximum_nominal_edges=layout['edges'][original],
            all_weights_nonzero_error_entries=int(finite.sum()),
            all_weights_entries_above_three_mc_errors=int(np.sum(finite&(np.abs(pull)>3.))),
            maximum_all_weights_absolute_pull=float(np.nanmax(np.abs(pull))))
        for variant in VARIANTS:
            selected=[start+i for i in (0,3,4)]
            for flavour in FLAVOURS:
                key=charge,variant,flavour
                rows=vectors[key][:,selected,0]; n=len(rows)
                between=np.var(rows,axis=0,ddof=1)/n
                within=conditional[key][:,selected].sum(axis=0)/n**2
                centered=rows-rows.mean(axis=0)
                report['retraining_groups']['__'.join((charge,variant,''.join(flavour)))]=dict(
                    rate_labels=['input','fiducial_1b','fiducial_2b'],values_pb=rows,
                    between_variance_of_mean=between,conditional_variance_of_mean=within,
                    between_over_conditional_variance=ratio(between,within),
                    largest_run_fraction_of_centered_sum_squares=ratio(
                        np.max(centered**2,axis=0),np.sum(centered**2,axis=0)),
                    leave_one_run_out_mean_pb=(rows.sum(axis=0)-rows)/(n-1))
        for local in (3,4,5,6,7,8,11,12,16,20):
            h=nominal+local-1
            a,b=map(int,layout['offsets'][h:h+2])
            parents=[start+i for i in normalization_bins(local)]
            selected=feature_vectors(nominal_vectors,charge,a,b,parents)
            normalized=estimate(selected,charges=(charge,),transform=lambda sums:
                prescription_observables(sums,charge,normalize_spectrum))
            report['original_shape_precision'][layout['titles'][h]]=shape_precision(normalized)
            arrays['original_shape_'+str(h)+'_edges']=layout['edges'][a:b]
            arrays['original_shape_'+str(h)+'_values']=normalized['value'][:,:,0]
            arrays['original_shape_'+str(h)+'_mc_errors']=normalized['mc_error'][:,:,0]
        print('Inspected main convergence:',charge,flush=True)
    for name,coefficients in INCLUSIVE_CONTRASTS:
        rows,weights=[],[]
        for charge in ('plus','minus'):
            for variant,coefficient in coefficients.items():
                for flavour in FLAVOURS:
                    rows.append(vectors[charge,variant,flavour][:,input_starts[charge]])
                    weights.append(coefficient)
        report['inclusive_charge_sum'][name]=diagnostic(linear_summary(rows,weights),WEIGHTS)
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),original_shape_labels=LABELS,**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),source_hashes={
        name:digest(STUDY/'scripts'/name) for name in ('main_convergence.py','main_report.py',
        'main_replicas.py','full_flavour_statistics.py','replica_statistics.py')})
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('replicas',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.replicas,args.output)
