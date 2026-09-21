#!/usr/bin/env python3
"""Narrow-W comparisons with a retrained on-shell mean and shared covariance.

Keep the original on-shell estimate as an independent convergence control.
Only the new reference uncertainty is measured between retrainings; the
original on-shell and BW controls still have conditional batch uncertainties.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from load_results import assert_same_layout, load
from narrow_w_report import SPECTRA, selected_vectors, validate_job
from narrow_w_statistics import estimate_retrained_reference
from pilot_report import RATES, clean
from run_narrow_w_pilots import DONE as PILOTS_DONE, cases as pilot_cases
from run_narrow_w_refinements import DONE as REFERENCE_DONE, cases as reference_cases
from run_width_mass_pilots import verified_pairs


def run(pilot_path,reference_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a retrained narrow-W comparison')
    pilots=json.loads(pilot_path.read_text())
    reference=json.loads(reference_path.read_text())
    if pilots['status']!=PILOTS_DONE or [j['case'] for j in pilots['jobs']]!=pilot_cases():
        raise ValueError('Require the complete original narrow-W pilots')
    if reference['status']!=REFERENCE_DONE or not reference['cases']:
        raise ValueError('Require completed independent narrow-W reference retrainings')
    expected=reference_cases(reference['cases'][0]['seed'],len(reference['cases']))
    if reference['cases']!=expected or [j['case'] for j in reference['jobs']]!=expected:
        raise ValueError('Reference retraining inventory differs from the prepared cases')
    inputs_path=Path(pilots['width_inputs'])
    if digest(inputs_path)!=digest(reference['width_inputs']):
        raise ValueError('Narrow-W width conventions differ')
    inputs=json.loads(inputs_path.read_text())
    if inputs['benchmark_sha256']!=digest(STUDY/'inputs/benchmark.json'):
        raise ValueError('Narrow-W benchmark changed')
    old_jobs=[j for j in pilots['jobs'] if j['variant']=='S']
    if len(old_jobs)!=7:
        raise ValueError('Require the original S reference and six BW controls')
    selected,refs,provenance,counts,streams=[],[],[],[],set()
    signature,layout=None,None
    for index,job in enumerate(reference['jobs']+old_jobs):
        result=load(job['audit'])
        current=validate_job(job,result,inputs,inputs_path)
        if signature is not None and current!=signature:
            raise ValueError('Reference and BW samples have different physical inputs')
        if layout is not None:
            assert_same_layout([layout,result])
        signature=current
        layout={k:result[k] for k in ('titles','offsets','edges')}
        pairs=verified_pairs(job)
        if streams & pairs:
            raise ValueError('Reference/BW runs share actual random initialization pairs')
        streams.update(pairs)
        counts.append(len(pairs))
        row=selected_vectors(job,result)
        if index<len(reference['jobs']):
            refs.append({kind:{key:value.sum(axis=0) for key,value in row[kind].items()}
                         for kind in ('rates','spectra')})
        else:
            selected.append(row)
        provenance.append(dict(case=job['case'],audit=job['audit'],audit_sha256=digest(job['audit']),
            batches=job['batches'],batches_sha256=digest(job['batches'])))
        print('Matched',job['run'],flush=True)
    strata=[row['strata'] for row in selected]
    factors=[1.]+[job['case']['width_factor'] for job in old_jobs]
    report=dict(created_utc=now(),variant='S',
        status='retrained-reference narrow-W comparison evaluated; inspect precision and convergence',
        inputs={str(p.resolve()):digest(p) for p in (pilot_path,reference_path,inputs_path)},
        reference_replica_count=len(refs),samples=provenance,
        comparison_samples=[dict(w_treatment='onshell',width_factor=1.,estimator='mean of fresh retrainings')]+
                           [job['case'] for job in old_jobs],
        stage_stream_counts=counts,distinct_initialization_pairs=len(streams),cross_sample_rng_overlap=0,
        rates={},spectra={},
        convention='Average complete on-shell reference cross-section vectors before ratios or normalization. '
            'Delete one whole reference run or one batch from one independent control stratum at a time. '
            'Reference spread includes reference sampling noise; never add it twice. '
            'Multiply each BW cross section by epsilon cubed once; share the reference covariance across controls. '
            'All nominal/81 scale/101 PDF coordinates are retained; nominal covariance factors are archived.',
        limitations='Only the fresh on-shell reference has between-retraining uncertainty. '
            'The old on-shell and six BW controls remain single learned-grid conditional estimates. '
            'Five reference runs do not by themselves establish tail coverage or convergence of the BW samples. '
            'Coarse diagnostic bins and one ordered W+ assignment; no full-flavour precision claim.')
    arrays={}
    def retain(prefix,result):
        arrays[prefix]=result['value']
        arrays[prefix+'_mc_errors']=result['mc_error']
        arrays[prefix+'_reference_retraining_mc_errors']=result['mc_error_reference_retraining']
        arrays[prefix+'_conditional_samples_mc_errors']=result['mc_error_conditional_samples']
        arrays[prefix+'_covariance_factor_nominal']=result['covariance_factor'][...,0]
        arrays[prefix+'_nonlinear_bias_estimate']=result['nonlinear_bias_estimate']
    for index,title in enumerate(selected[0]['rates']):
        result=estimate_retrained_reference(np.asarray([r['rates'][title] for r in refs]),
            [r['rates'][title] for r in selected],strata,factors)
        retain('rates_%d'%index,result)
        report['rates'][title]={label:dict(
            rescaled_pb=result['value'][0,:,i,0],rescaled_mc_error_pb=result['mc_error'][0,:,i,0],
            difference_pb=result['value'][1,:,i,0],difference_mc_error_pb=result['mc_error'][1,:,i,0],
            ratio=result['value'][2,:,i,0],ratio_mc_error=result['mc_error'][2,:,i,0],
            ratio_reference_retraining_error=result['mc_error_reference_retraining'][2,:,i,0],
            ratio_conditional_samples_error=result['mc_error_conditional_samples'][2,:,i,0])
            for i,label in enumerate(RATES)}
        del result
    for local in SPECTRA:
        reference_vectors=np.asarray([r['spectra'][local] for r in refs])
        vectors=[r['spectra'][local] for r in selected]
        for normalized in (False,True):
            result=estimate_retrained_reference(
                reference_vectors if normalized else reference_vectors[:,:-1],
                vectors if normalized else [v[:,:-1] for v in vectors],strata,factors,normalized)
            prefix=('normalized_' if normalized else 'absolute_')+str(local)
            retain(prefix,result)
            if normalized:
                pulls=np.divide(result['value'][1,:,:,0],result['mc_error'][1,:,:,0],
                    out=np.full_like(result['value'][1,:,:,0],np.nan),where=result['mc_error'][1,:,:,0]>0.)
                report['spectra'][str(local)]=dict(title=layout['titles'][local-1],
                    parent_region='2b' if local in (11,12) else '1b',nominal_normalized_residuals_over_error=pulls,
                    maximum_absolute_nominal_residual_by_sample=[
                        float(np.nanmax(np.abs(p))) if np.isfinite(p).any() else None for p in pulls])
            del result
        arrays['edges_'+str(local)]=selected[0]['edges'][local]
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),script_sha256=digest(Path(__file__)),
                  statistics_sha256=digest(STUDY/'scripts/narrow_w_statistics.py'))
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pilots',type=Path)
    parser.add_argument('--reference-queue',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.pilots,args.reference_queue,args.output)
