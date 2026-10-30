#!/usr/bin/env python3
"""Project nominal-PDF stage terms from the main signed decay-scale weights.

Production coordinates and PDFs stay fixed during each projection. Undo the
known top-width denominators before fitting the bilinear product numerator.
This is an algebraic diagnostic of stored weights, not a new subset export,
a replacement for independent controls, or a convergence certification.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from full_flavour_statistics import FLAVOURS
from load_results import POINTS
from main_convergence import linear_summary
from main_report import read_main
from pilot_report import RATES, clean

PRODUCTION=tuple(sorted({p[:2] for p in POINTS}))
CENTRAL=PRODUCTION.index((1.,1.))
TERMS=('000','100','010','001','110','101','011','111',
       'A0','A1','A2','A3','strict_reconstructed',
       'production_decay_second_order','decay_second_order','total_second_order',
       'product_reconstructed','delta_reconstructed','higher_remainder')


def width_parameters(benchmark):
    rows={r['scale_factor']:r for r in benchmark['top_widths']
          if r['mb_GeV']==0. and r['w_treatment']=='onshell'}
    if set(rows)!={.5,1.,2.}:
        raise ValueError('Require the three matched massless on-shell widths')
    ratios={f:r['alpha_s_pdf']/rows[1.]['alpha_s_pdf'] for f,r in rows.items()}
    g=rows[1.]['gamma_nlo_pdf']/rows[1.]['gamma_lo']-1.
    for f,r in rows.items():
        if not np.isclose(r['gamma_nlo_pdf']/r['gamma_lo']-1.,g*ratios[f],rtol=1.e-12,atol=1.e-14):
            raise ValueError('Widths do not have the required single coupling power')
    return ratios,g


def project(values,ratios,g,product):
    """Return [...,production point,coefficient], checking all 81 weights.

    The four product coefficients multiply 1,r_t,r_tbar,r_t*r_tbar after
    multiplying by (1+g*r_t)(1+g*r_tbar). Strict coefficients omit the last
    monomial and retain the expanded width counterterms. The outer two
    scale factors determine the projection; the five remaining points
    independently check its form at every bin and in every complete run.
    """
    x=np.asarray(values,dtype=float)
    if x.ndim<2 or x.shape[-1]<82 or not np.isfinite(x).all():
        raise ValueError('Require finite vectors with the complete signed scale grid')
    if set(ratios)!={.5,1.,2.} or ratios[1.]!=1. or not all(np.isfinite(list(ratios.values()))):
        raise ValueError('Require finite relative couplings at the three scales')
    lo,hi=ratios[2.],ratios[.5]
    delta=hi-lo
    if abs(delta)<1.e-3 or any(1.+g*r<=0. for r in ratios.values()):
        raise ValueError('Degenerate coupling projection or nonpositive width')
    result=[];maximum=0.; scaled_maximum=0.
    # HwU text and accumulated signed contributions have finite precision.
    # A per-bin floor also covers cancellations in nearly zero bins.
    magnitude=np.max(np.abs(x[...,1:82]),axis=-1)
    tolerance=2.e-5*magnitude+2.e-10
    for prod in PRODUCTION:
        def value(t,a):
            y=x[...,1+POINTS.index(prod+(t,a))]
            return y*(1.+g*ratios[t])*(1.+g*ratios[a]) if product else y
        ll,hl,lh,hh=value(2.,2.),value(.5,2.),value(2.,.5),value(.5,.5)
        cross=(hh-hl-lh+ll)/delta**2 if product else np.zeros_like(ll)
        top=(hl-ll)/delta-lo*cross
        antitop=(lh-ll)/delta-lo*cross
        base=ll-lo*(top+antitop)-lo*lo*cross
        coefficients=np.stack((base,top,antitop,cross) if product else (base,top,antitop),axis=-1)
        for t in (.5,1.,2.):
            for a in (.5,1.,2.):
                predicted=base+ratios[t]*top+ratios[a]*antitop+ratios[t]*ratios[a]*cross
                if product:predicted=predicted/((1.+g*ratios[t])*(1.+g*ratios[a]))
                observed=x[...,1+POINTS.index(prod+(t,a))]
                residual=np.abs(predicted-observed)
                maximum=max(maximum,float(residual.max()))
                scaled_maximum=max(scaled_maximum,float(np.max(residual/tolerance)))
                if np.any(residual>tolerance):
                    raise ValueError('Stored decay-scale weights do not close the claimed polynomial')
        result.append(coefficients)
    return np.stack(result,axis=-2),dict(maximum_absolute_closure_pb=maximum,
        maximum_fraction_of_rounding_tolerance=scaled_maximum,
        rounding_tolerance='2e-5 * maximum absolute scale weight in that run/bin + 2e-10 pb',
        all_81_scale_points_and_original_bins_checked=True)


def stage_contribution(coefficients,variant,g):
    """Linear contribution of one independent Pi/PiD ensemble to every term."""
    if variant not in ('Pi','PiD'):
        raise ValueError('Stage reconstruction requires product ensembles')
    a,b,c,d=np.moveaxis(np.asarray(coefficients),-1,0)
    zero=np.zeros_like(a)
    stages=(zero,a,zero,zero,b,c,zero,d) if variant=='Pi' else (a,-a,b,c,-b,-c,d,-d)
    a0=stages[0]
    a1=stages[1]+stages[2]+stages[3]
    a2=stages[4]+stages[5]+stages[6]
    a3=stages[7]
    strict=a0+a1-2.*g*a0
    mixed=stages[4]+stages[5]-2.*g*stages[1]
    decay=stages[6]-2.*g*(stages[2]+stages[3])+3.*g*g*a0
    second=a2-2.*g*a1+3.*g*g*a0
    product=(a0+a1+a2+a3)/(1.+g)**2
    shift=product-strict
    return np.stack((*stages,a0,a1,a2,a3,strict,mixed,decay,second,product,shift,shift-second),axis=-1)


def compact(summary,labels):
    return {label:dict(value_pb=float(summary['value'][CENTRAL,j]),
        mc_error_pb=float(summary['mc_error'][CENTRAL,j]),
        effective_dof=float(summary['effective_dof'][CENTRAL,j])) for j,label in enumerate(labels)}


def run(path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite stage-projection evidence')
    source,layout,vectors,_=read_main(path)
    ratios,g=width_parameters(json.loads((STUDY/'inputs/benchmark.json').read_text()))
    banks=[h for h,t in enumerate(layout['titles']) if ' rates:' in t]
    selected=np.concatenate([np.arange(layout['offsets'][h],layout['offsets'][h+1]) for h in banks])
    coefficients={}; checks={}; arrays=dict(rate_original_bin_indices=selected,
        production_coordinates=PRODUCTION,term_labels=TERMS)
    for key,x in vectors.items():
        charge,variant,flavour=key
        if variant not in ('S','D','Pi','PiD'):continue
        projected,closure=project(x,ratios,g,variant in ('Pi','PiD'))
        coefficients[key]=projected[:,selected].copy()
        name='__'.join((charge,variant,''.join(flavour)))
        arrays['coefficients__'+name]=coefficients[key]
        checks[name]=closure
    report=dict(created_utc=now(),status='stored-weight stage projection; independent subset validation pending',
        input=str(path.resolve()),input_sha256=digest(path),input_arrays_sha256=source['arrays_sha256'],
        benchmark_sha256=source['benchmark_sha256'],coupling_ratios=ratios,width_correction=g,
        production_coordinates=PRODUCTION,term_labels=TERMS,closure_checks=checks,rates={},
        strict_decay_diagnostics={},new_integrations=0,publication_certified=False,
        conventions='Each stage label lists production/top/antitop numerator order in {0,1}; '
            'all numerator terms use the common LO-width denominator. Production scales and '
            'member-0 PDFs stay fixed during each decay-scale projection. Independent groups '
            'are never paired. All terms are formed within a complete run before their covariance '
            'is estimated. NPZ stores complete projected run coefficients for shared covariance.',
        limitations='This reconstructs the scale-weight polynomial; it is not an independent '
            'eight-subset export or a new integration. It neither resolves adaptive-stopping '
            'coverage nor replaces direct S/P/D/LO controls. Only member-0 coefficients can '
            'be projected because PDF replicas lack the joint decay-scale grid. No missing '
            'same-stage second-order correction is supplied. Differential attribution is pending.')
    for charge in ('plus','minus'):
        product_rows=[stage_contribution(coefficients[charge,v,f],v,g)
                      for v in ('Pi','PiD') for f in FLAVOURS]
        summary=linear_summary(product_rows)
        arrays[charge+'__term_values']=summary['value']
        arrays[charge+'__term_mc_errors']=summary['mc_error']
        # Compare reconstructed prescriptions to independent direct controls.
        closure_rows=[]
        for v in ('Pi','PiD'):
            for f in FLAVOURS:
                q=coefficients[charge,v,f];a,b,c,_=np.moveaxis(q,-1,0);z=np.zeros_like(a)
                closure_rows.append(np.stack((a if v=='Pi' else z,a if v=='PiD' else z,
                    a+b+c-2*g*a if v=='PiD' else z,a if v=='Pi' else b+c-2*g*a),axis=-1))
        for i,v in enumerate(('P','LO','D','S')):
            for f in FLAVOURS:
                x=vectors[charge,v,f][:,selected][:,:,[1+POINTS.index(p+(1.,1.)) for p in PRODUCTION]]
                row=np.zeros(x.shape+(4,));row[...,i]=-x;closure_rows.append(row)
        closure_summary=linear_summary(closure_rows)
        for h in banks:
            if ('W+' if charge=='plus' else 'W-') not in layout['titles'][h]:continue
            start=int(np.flatnonzero(selected==layout['offsets'][h])[0])
            report['rates'][layout['titles'][h]]={label:dict(
                terms=compact({k:v[start+j] for k,v in summary.items()},TERMS),
                independent_control_residuals=compact({k:v[start+j] for k,v in closure_summary.items()},
                    ('P','LO','D','S'))) for j,label in enumerate(RATES)}
        for variant,reference in (('S','P'),('D','LO')):
            rows=[]
            for f in FLAVOURS:
                q=coefficients[charge,variant,f];base,top,anti=np.moveaxis(q,-1,0)
                rows.append(np.stack((base,top,anti,top+anti,base,base+top+anti),axis=-1))
            for f in FLAVOURS:
                x=vectors[charge,reference,f][:,selected][:,:,[1+POINTS.index(p+(1.,1.)) for p in PRODUCTION]]
                row=np.zeros(x.shape+(6,));row[...,4]=-x;row[...,5]=-x;rows.append(row)
            reduced=linear_summary(rows)
            h=next(h for h in banks if layout['titles'][h].startswith('R04_b25 '+('W+' if charge=='plus' else 'W-')))
            start=int(np.flatnonzero(selected==layout['offsets'][h])[0])
            report['strict_decay_diagnostics'][charge+'__'+variant]=compact(
                {k:v[start] for k,v in reduced.items()},('baseline','top','antitop','total_decay',
                    'baseline_minus_independent_control','total_minus_independent_control'))
        print('Projected stored main stage weights:',charge,flush=True)
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),source_hashes={
        name:digest(STUDY/'scripts'/name) for name in ('main_stage_projection.py','main_convergence.py',
            'main_report.py','main_replicas.py','full_flavour_statistics.py','replica_statistics.py')})
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('replicas',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    run(args.replicas,args.output)
