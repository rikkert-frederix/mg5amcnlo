#!/usr/bin/env python3
"""Pilot-based allocation diagnostics, never a substitute for main predictions.

Apply T ~ T0*(e0/e)^2 to measured conditional MC errors. The effective
flavour-count range is an allocation scenario, not a physics flavour rescaling.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY,digest,now,save
from joint_report import estimate
from load_results import load
from narrow_w_report import selected_vectors
from pilot_report import RATES,clean
from replica_statistics import ratio
from run_width_mass_pilots import verified_pairs


def forecast(measured_relative_errors,measured_total_error,target_total_error,replicas,pilot_seconds):
    if (not 0<measured_total_error<1 or not 0<target_total_error<1 or replicas<5 or
            not np.isfinite(pilot_seconds) or pilot_seconds<=0):
        raise ValueError('Invalid measured-error or allocation inputs')
    factor=target_total_error/measured_total_error
    per_run=np.asarray(measured_relative_errors)*factor
    return dict(projected_per_run_relative_error=per_run,
        projected_five_or_more_run_relative_error=per_run/np.sqrt(replicas),
        projected_full_sum_relative_error_effective_flavours_1_to_8=np.stack(
            [per_run/np.sqrt(replicas),per_run/np.sqrt(8*replicas)]),
        naive_seconds_per_run=pilot_seconds/factor**2,
        effective_flavours_needed_for_one_percent_fiducial=(per_run/(.01*np.sqrt(replicas)))**2)


def run(main_path,output):
    if output.exists():raise ValueError('Refuse to overwrite an allocation forecast')
    main=json.loads(main_path.read_text())
    if main['status']!='prepared; technical/reference release required' or main['jobs']:
        raise ValueError('This forecast applies to the unexecuted main allocation')
    replicas=max(case['replica'] for case in main['cases'])+1
    paths=[STUDY/'inputs/rng_pilot_queue_alignment_controls_v1.json',
           STUDY/'inputs/width_mass_pilot_queue_width_mass_v4.json']
    jobs=[job for path in paths for job in json.loads(path.read_text())['jobs']]
    rows,inputs,streams=[],{str(p.resolve()):digest(p) for p in [main_path]+paths},set()
    for job in jobs:
        result=load(job['audit'])
        execution=result['report']['execution']
        generation=execution['export_generation']['args']
        if (generation['w_treatment']!='onshell' or generation['decay_bottom_mass']!=0. or
                result['report']['variant'] not in ('S','Pi')):
            continue
        pairs=verified_pairs(job)
        if streams & pairs:raise ValueError('Repeated stage stream in allocation pilots')
        streams.update(pairs)
        inputs[str(Path(job['audit']).resolve())]=digest(job['audit'])
        inputs[str(Path(job['batches']).resolve())]=digest(job['batches'])
        # This helper's spectrum selection is for W+. For W- use the same
        # actual local IDs/parent bins in its own charge block.
        with np.load(job['batches']) as data:
            values=data['contributions'][:,:,0:1]
            strata=data['strata'].copy()
        charge=generation['charge']; marker='W+' if charge=='plus' else 'W-'
        h=next(i for i,title in enumerate(result['titles']) if title.startswith('R04_b25')
               and marker in title and ' rates:' in title)
        start=int(result['offsets'][h])
        rates=estimate(values[:,start:start+9],strata)
        relative=ratio(rates['mc_error'][:,0],rates['value'][:,0])
        measured_total=float(relative[0])
        prediction=forecast(relative,measured_total,main['settings']['accuracy'],replicas,execution['elapsed_s'])
        spectra={}
        for local in (3,4,6,11,12,20):
            a,b=map(int,result['offsets'][h+local-1:h+local+1])
            parent=start+(4 if local in (11,12) else 3)
            features=np.concatenate([values[:,a:b],values[:,parent:parent+1]],axis=1)
            shape=estimate(features,strata,lambda total:ratio(total[:-1],total[-1]))
            rel=ratio(shape['mc_error'][:,0],shape['value'][:,0])
            valid=np.isfinite(rel)&(rel>0)
            projected=rel[valid]*main['settings']['accuracy']/measured_total/np.sqrt(replicas)
            spectra[result['titles'][h+local-1]]=dict(
                positive_nonzero_error_bins=int(valid.sum()),
                measured_relative_error_quantiles=np.quantile(rel[valid],[.1,.5,.9]) if valid.any() else None,
                projected_combined_error_quantiles_effective_flavours_1=(np.quantile(projected,[.1,.5,.9]) if valid.any() else None),
                projected_combined_error_quantiles_effective_flavours_8=(np.quantile(projected/np.sqrt(8),[.1,.5,.9]) if valid.any() else None),
                projected_bins_below_two_percent_effective_flavours_1=int(np.sum(projected<.02)),
                projected_bins_below_two_percent_effective_flavours_8=int(np.sum(projected/np.sqrt(8)<.02)))
        manifest=result['report']['manifest']
        rows.append(dict(charge=charge,variant=result['report']['variant'],ordered_flavours=generation['flavours'],
            pilot_run=job['run'],pilot_elapsed_seconds=execution['elapsed_s'],
            pilot_requested_total_accuracy=manifest['settings']['req_acc_fo'],
            measured_total_relative_error=measured_total,rate_labels=RATES,
            measured_rate_relative_errors=relative,forecast=prediction,spectra=spectra,
            pilot_grid_points=manifest['settings']['npoints_fo_grid'],
            main_grid_points=main['settings']['grid_points'],stage_initializations=len(pairs)))
        print(charge,result['report']['variant'],'measured total error',measured_total,flush=True)
    if {(row['charge'],row['variant']) for row in rows}!={(c,v) for c in ('plus','minus') for v in ('S','Pi')} or len(rows)!=4:
        raise ValueError('Require four valid matched S/Pi allocation pilots')
    save(output,clean(dict(created_utc=now(),status='pilot-based main allocation forecast; remeasure actual full-flavour convergence',
        inputs=inputs,script_sha256=digest(Path(__file__)),main_total_accuracy=main['settings']['accuracy'],
        independent_retrainings=replicas,pilots=rows,
        convention='Use measured pilot errors with inverse-square cost scaling. Under comparable per-flavour '
            'relative variance, the effective number of positive flavour contributions is (sum sigma)^2/sum sigma^2 '
            'and lies between 1 and 8. Both endpoints are allocation scenarios; no full-flavour rate is inferred.',
        limitations='Pilot errors are conditional on their trained grids. Runtime/source and 200-to-1000-point '
            'training changes, overhead, rare tails and adaptive stopping invalidate precise cost extrapolation. '
            'Between-retraining variation and the actual flavour weights remain unmeasured. The 1.5% total-rate '
            'request alone does not certify any fiducial/shape target. Rebin or increase statistics based on main data; '
            'LO/P/D/PiD costs must be remeasured, not inferred from these S/Pi pilots.')))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('main',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.main,args.output)
