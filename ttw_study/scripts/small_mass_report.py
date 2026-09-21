#!/usr/bin/env python3
"""Generated small-mass S/Pi rates and shapes with a shared massless control."""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from load_results import load
from mass_effects import LABELS, contrasts, validate_physics
from narrow_w_report import SPECTRA, selected_vectors
from pilot_report import RATES, clean
from read_splits import split_ensembles
from replica_statistics import independent_jackknife_many, ratio
from run_small_mass_checks import DONE, cases, local_benchmark
from run_width_mass_pilots import verified_pairs


def estimate(vectors,strata,normalized=False):
    """Jointly compare masses 1/0.1 to 0, including common-control covariance."""
    if len(vectors)!=6 or len(strata)!=6:
        raise ValueError('Require independent S/Pi samples at masses 0,1,0.1 in that order')
    groups=[split_ensembles(v,s) for v,s in zip(vectors,strata)]
    boundaries=np.cumsum([0]+[len(group) for group in groups])
    def transform(*means):
        rows=[sum(means[a:b]) for a,b in zip(boundaries[:-1],boundaries[1:])]
        if normalized:
            rows=[ratio(row[:-1],row[-1]) for row in rows]
        return np.asarray([contrasts(*rows[:2],*rows[2:4]),contrasts(*rows[:2],*rows[4:6])])
    return independent_jackknife_many([row for group in groups for row in group],transform)


def run(queue_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a generated small-mass comparison')
    queue=json.loads(queue_path.read_text())
    if queue['status']!=DONE or not queue['cases']:
        raise ValueError('Require the complete generated small-mass queue')
    expected=cases(queue['cases'][0]['seed'])
    if queue['cases']!=expected or [j['case'] for j in queue['jobs']]!=expected:
        raise ValueError('Small-mass case inventory differs from the complete allocation')
    inputs_path=Path(queue['width_inputs'])
    local_benchmark(json.loads(inputs_path.read_text()))
    report=dict(created_utc=now(),status='generated small-mass comparisons evaluated; inspect continuity sensitivity',
        inputs={str(p.resolve()):digest(p) for p in (queue_path,inputs_path)},
        nonzero_decay_masses_GeV=[1.,.1],labels=LABELS,modes={},
        convention='Six independent S/Pi sample/beam-stratum ensembles per W treatment. '
            'Keep the same massless samples across the two finite-mass comparisons and propagate that covariance. '
            'Form ratios from summed vectors and normalize bins by their matching fiducial rate before weight reduction. '
            'All nominal/81 scale/101 PDF coordinates are retained, with nominal covariance factors.',
        limitations='One ordered W+ e,e,mu assignment and conditional learned-grid errors. '
            'Compatibility at finite small masses bounds sensitivity; it does not prove the mathematical massless limit, '
            'independent-retraining coverage, rare-tail convergence, or full-flavour precision. '
            'Coarse spectra are diagnostic bins, not the frozen publication choice.')
    arrays={}
    streams=set()
    def retain(prefix,result):
        arrays[prefix]=result['value']
        arrays[prefix+'_mc_errors']=result['mc_error']
        arrays[prefix+'_covariance_factor_nominal']=result['covariance_factor'][...,0]
        arrays[prefix+'_nonlinear_bias_estimate']=result['nonlinear_bias_estimate']
    for mode in ('onshell','all-bw'):
        jobs=[j for j in queue['jobs'] if j['case']['w_treatment']==mode]
        results=[load(j['audit']) for j in jobs]
        identities=[]
        for mass,index in ((1.,2),(.1,4)):
            proof=Path(queue['exports'][mode+'_mb'+str(mass).replace('.','p')]['production_identity'])
            validate_physics(results[:2]+results[index:index+2],[proof],width_inputs=inputs_path)
            identities.append(proof)
        selected,source_rows,counts=[],[],[]
        for job,result in zip(jobs,results):
            manifest=result['report']['manifest']
            generation=result['report']['execution']['export_generation']['args']
            expected_generation=dict(charge='plus',flavours=['e','e','mu'],corrected='both',
                w_treatment=mode,decay_bottom_mass=job['case']['decay_bottom_mass'])
            if (manifest['small_mass_test']!=job['case'] or manifest['production_scale']!='core-w-ht-half' or
                    manifest['decay_scale_grouping']!='separate' or manifest['generated_event_weight_rescaling']!=1. or
                    any(generation[k]!=v for k,v in expected_generation.items())):
                raise ValueError('Generated small-mass convention differs from its allocation')
            pairs=verified_pairs(job)
            if pairs & streams:
                raise ValueError('Small-mass comparisons reuse actual training/refinement streams')
            streams.update(pairs)
            counts.append(len(pairs))
            selected.append(selected_vectors(job,result))
            source_rows.append(dict(case=job['case'],audit=job['audit'],audit_sha256=digest(job['audit']),
                batches=job['batches'],batches_sha256=digest(job['batches'])))
        strata=[row['strata'] for row in selected]
        mode_report=dict(samples=source_rows,stage_stream_counts=counts,rates={},spectra={},
            production_identity_inputs={str(p):digest(p) for p in identities})
        for index,title in enumerate(selected[0]['rates']):
            result=estimate([r['rates'][title] for r in selected],strata)
            retain('%s_rates_%d'%(mode,index),result)
            mode_report['rates'][title]={label:{name:dict(value=result['value'][:,k,i,0],
                mc_error=result['mc_error'][:,k,i,0]) for k,name in enumerate(LABELS)}
                for i,label in enumerate(RATES)}
            del result
        for local in SPECTRA:
            vectors=[r['spectra'][local] for r in selected]
            for normalized in (False,True):
                result=estimate(vectors if normalized else [v[:,:-1] for v in vectors],strata,normalized)
                prefix=mode+('_normalized_' if normalized else '_absolute_')+str(local)
                retain(prefix,result)
                if normalized:
                    residuals=result['value'][:,:2,:,0]-1.
                    errors=result['mc_error'][:,:2,:,0]
                    pulls=np.divide(residuals,errors,out=np.full_like(residuals,np.nan),where=errors>0.)
                    mode_report['spectra'][str(local)]=dict(title=results[0]['titles'][local-1],
                        parent_region='2b' if local in (11,12) else '1b',
                        nominal_mass_ratios_residuals_over_error=pulls,
                        maximum_absolute_nominal_residual_by_mass_and_variant=[
                            [float(np.nanmax(np.abs(p))) if np.isfinite(p).any() else None for p in row]
                            for row in pulls])
                del result
            arrays[mode+'_edges_'+str(local)]=selected[0]['edges'][local]
        report['modes'][mode]=mode_report
        print('Evaluated small-mass controls:',mode,flush=True)
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),script_sha256=digest(Path(__file__)),
                  physics_gate_sha256=digest(STUDY/'scripts/mass_effects.py'),cross_sample_rng_overlap=0,
                  distinct_initialization_pairs=len(streams))
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.queue,args.output)
