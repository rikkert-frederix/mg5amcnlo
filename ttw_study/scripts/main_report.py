#!/usr/bin/env python3
"""Full-flavour main rates/shapes with independent-retraining covariance.

Original bins remain available in the input archive. An optional common
plan merges complete run vectors before reduction. Per-spectrum compressed
JSON contains the pointwise scale/PDF reductions; bin choices and final
precision still require inspection of the measured convergence.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from covariance_storage import nominal_copy

from campaign import ROOT, STUDY, digest, now, save
from full_flavour_statistics import FLAVOURS,VARIANTS,estimate,inventory as ensemble_inventory
from joint_report import derived_rates,normalization_bins
from load_results import WEIGHTS
from main_observables import (ABSOLUTE_LABELS,CHARGE_LABELS,CHARGE_SHAPE_LABELS,LABELS,charge_observables,
    normalize_spectrum,prescription_observables)
from main_replicas import inventory as campaign_inventory
from pilot_report import RATES,clean,describe
from replica_statistics import ratio

DERIVED=('acceptance_1b','acceptance_2b','two_b_given_one_b','veto_0extra',
         'fraction_1extra','fraction_2extra','fraction_3plus_extra')


def read_main(path):
    record=json.loads(path.read_text())
    if record['status']!='complete independent main vectors collected; flavour reduction and convergence required':
        raise ValueError('Require collected complete independent main vectors')
    if (digest(path.with_suffix('.npz'))!=record['arrays_sha256'] or
            digest(record['input'])!=record['input_sha256'] or
            digest(STUDY/'inputs/benchmark.json')!=record['benchmark_sha256']):
        raise ValueError('Main replica archive, inventory or benchmark changed')
    count=campaign_inventory(json.loads(Path(record['input']).read_text()))
    vectors,conditional={},{}
    with np.load(path.with_suffix('.npz')) as data:
        layout={k:data[k].copy() for k in ('edges','offsets')}
        layout.update(titles=list(data['titles']),weights=list(data['weights']))
        if layout['weights']!=WEIGHTS or int(layout['offsets'][-1])!=len(layout['edges']):
            raise ValueError('Unexpected main weight or histogram layout')
        for name,group in record['groups'].items():
            key=group['charge'],group['variant'],tuple(group['flavours'])
            if key in vectors or group['replicas']!=list(range(count)):
                raise ValueError('Duplicated or incomplete retraining group')
            value=data[name].copy()
            expected=(count,len(layout['edges']),len(WEIGHTS))
            if value.shape!=expected:
                raise ValueError('Main group has incomplete bins, weights or retrainings')
            vectors[key]=value
            # Retain only nominal within-run diagnostics in memory. The
            # source archive keeps their complete scale/PDF coordinates.
            conditional[key]=data[name+'__conditional_variances'][:,:,0].copy()
    ensemble_inventory(vectors,('plus','minus'),VARIANTS)
    return record,layout,vectors,conditional


def feature_vectors(vectors,charge,start,stop,parent=None):
    result={}
    for key,values in vectors.items():
        if key[0]!=charge:
            continue
        row=values[:,start:stop]
        if parent is not None:
            row=np.concatenate([row,values[:,parent].sum(axis=1)[:,None,:]],axis=1)
        result[key]=row
    return result


def summaries(result,labels,bin_index=None):
    output={}
    for index,label in enumerate(labels):
        value=result['value'][index] if bin_index is None else result['value'][index,bin_index]
        error=result['mc_error'][index] if bin_index is None else result['mc_error'][index,bin_index]
        bias=result['nonlinear_bias_estimate'][index] if bin_index is None else result['nonlinear_bias_estimate'][index,bin_index]
        row=describe(value)
        row.update(nominal_retraining_mc_error=error[0],nominal_jackknife_bias_estimate=bias[0])
        if label in VARIANTS:
            row['nominal_relative_mc_error']=ratio(error[0],value[0])
        if label in ('Pi_minus_S','PiD_minus_D','production_related_product_shift'):
            row['nominal_shift_over_mc_error']=ratio(value[0],error[0])
            row['shift_exceeds_three_mc_errors']=bool(np.isfinite(error[0]) and error[0]>0. and abs(value[0])>3*error[0])
        if label=='S_minus_P_minus_D_plus_LO':
            row['nominal_identity_residual_over_mc_error']=ratio(value[0],error[0])
        output[label]=row
    return output


def shape_precision(result):
    output={}
    for index,variant in enumerate(VARIANTS):
        value,error=result['value'][index,:,0],result['mc_error'][index,:,0]
        valid=np.isfinite(value)&np.isfinite(error)&(value>0.)&(error>0.)
        relative=error[valid]/value[valid]
        output[variant]=dict(positive_nonzero_error_bins=int(valid.sum()),
            bins_below_two_percent=int(np.sum(relative<.02)),
            bins_with_undefined_error=int(np.sum(~np.isfinite(error))),
            bins_with_nonpositive_central=int(np.sum(value<=0.)),
            relative_error_quantiles=np.quantile(relative,[.1,.5,.9]) if relative.size else None)
    return output


def run(path,output,binning_path=None):
    spectra_dir=output.parent/(output.stem+'_spectra')
    if output.exists() or output.with_suffix('.npz').exists() or spectra_dir.exists():
        raise ValueError('Refuse to overwrite main reductions or partial spectrum archives')
    source,layout,vectors,conditional=read_main(path)
    binning=None
    if binning_path is not None:
        from main_binning import apply_plan
        content=binning_path.read_bytes()
        plan=json.loads(content)
        layout,vectors,conditional,binning=apply_plan(layout,vectors,conditional,plan)
        binning.update(path=str(binning_path.resolve()),sha256=hashlib.sha256(content).hexdigest(),plan=plan)
    spectra_dir.mkdir(parents=True)
    nbin=len(layout['edges'])
    arrays={name:np.full((len(labels),nbin,len(WEIGHTS)),np.nan) for name,labels in
            (('absolute',ABSOLUTE_LABELS),('absolute_mc_errors',ABSOLUTE_LABELS),
             ('normalized',LABELS),('normalized_mc_errors',LABELS))}
    arrays.update(layout,absolute_labels=ABSOLUTE_LABELS,normalized_labels=LABELS,
                  charge_labels=CHARGE_LABELS,charge_shape_labels=CHARGE_SHAPE_LABELS)
    report=dict(created_utc=now(),status='full-flavour main observables reduced; inspect precision and convergence',
        input=str(path.resolve()),input_sha256=digest(path),input_arrays_sha256=source['arrays_sha256'],
        benchmark_sha256=source['benchmark_sha256'],rates={},derived_rates={},spectra={},
        migrations={},charge_rates={},charge_spectra={},retraining_diagnostics={},binning=binning,
        convention='Average independently retrained complete vectors within each ordered flavour, '
            'then sum eight disjoint assignments per charge/prescription. Delete from one flavour ensemble '
            'at a time; never pair independent runs by replica index. Form all ratios and normalized shapes '
            'at identical scale/PDF coordinates before reduction. Nominal covariance factors are archived; '
            'the source run vectors retain the full cross-weight covariance. No extra within-run variance is added.',
        limitations='Numerical precision flags do not prove coverage, rare-tail convergence or absence of '
            'adaptive-stopping bias. Full original ranges are retained, with no folded continuous overflows; '
            'freeze common publication bins only after inspection. W-treatment, massive-decay, parameter '
            'variations and conditional coefficient attribution remain separate required results.')
    def keep(prefix,result):
        arrays[prefix]=result['value']
        arrays[prefix+'_mc_errors']=result['mc_error']
        arrays[prefix+'_covariance_factor_nominal']=nominal_copy(result['covariance_factor'])
        arrays[prefix+'_nonlinear_bias_nominal']=nominal_copy(result['nonlinear_bias_estimate'])
    def write_spectrum(name,detail):
        target=spectra_dir/(name+'.json.gz')
        target.write_bytes(gzip.compress(json.dumps(clean(detail),sort_keys=True,allow_nan=False).encode(),mtime=0))
        return dict(path=str(target.resolve()),sha256=digest(target))
    rates=[(h,title,'plus' if 'W+' in title else 'minus') for h,title in enumerate(layout['titles'])
           if ' rates:' in title and ('W+' in title or 'W-' in title)]
    if len(rates)!=10 or any(sum(charge==c for _,_,charge in rates)!=5 for c in ('plus','minus')):
        raise ValueError('Require five complete cut configurations for each charge')
    baseline={charge:next((h,title) for h,title,c in rates if c==charge and title.startswith('R04_b25'))
              for charge in ('plus','minus')}
    for h,title,charge in rates:
        start,stop=map(int,layout['offsets'][h:h+2])
        if stop-start!=9:
            raise ValueError('Incomplete main rate histogram')
        selected=feature_vectors(vectors,charge,start,stop)
        absolute=estimate(selected,charges=(charge,),transform=lambda sums:prescription_observables(sums,charge))
        arrays['absolute'][:,start:stop]=absolute['value']
        arrays['absolute_mc_errors'][:,start:stop]=absolute['mc_error']
        keep('rates_%d'%h,absolute)
        report['rates'][title]={label:summaries(absolute,ABSOLUTE_LABELS,i) for i,label in enumerate(RATES)}
        derived=estimate(selected,charges=(charge,),transform=lambda sums:prescription_observables(sums,charge,derived_rates))
        keep('derived_rates_%d'%h,derived)
        report['derived_rates'][title]={label:summaries(derived,LABELS,i) for i,label in enumerate(DERIVED)}
        del absolute,derived
        for local in range(2,22):
            index=h+local-1
            a,b=map(int,layout['offsets'][index:index+2])
            parents=[start+i for i in normalization_bins(local)]
            selected=feature_vectors(vectors,charge,a,b,parents)
            absolute=estimate({k:v[:,:-1] for k,v in selected.items()},charges=(charge,),
                transform=lambda sums:prescription_observables(sums,charge))
            normalized=estimate(selected,charges=(charge,),
                transform=lambda sums:prescription_observables(sums,charge,normalize_spectrum))
            coverage=estimate(selected,charges=(charge,),transform=lambda sums:
                prescription_observables(sums,charge,lambda total:normalize_spectrum(total).sum(axis=0)))
            for prefix,result in (('absolute',absolute),('normalized',normalized)):
                arrays[prefix][:,a:b]=result['value']
                arrays[prefix+'_mc_errors'][:,a:b]=result['mc_error']
                arrays[prefix+'_covariance_factor_nominal_%d'%index]=nominal_copy(result['covariance_factor'])
            precision=shape_precision(normalized)
            detail=dict(title=layout['titles'][index],edges=layout['edges'][a:b],
                normalization_rate_bins=[RATES[i] for i in normalization_bins(local)],
                units='Absolute entries are cross sections per bin in pb; normalized entries are bin fractions.',
                absolute=[summaries(absolute,ABSOLUTE_LABELS,i) for i in range(b-a)],
                normalized=[summaries(normalized,LABELS,i) for i in range(b-a)],
                measured_range_coverage=summaries(coverage,LABELS),precision=precision,
                interpretation='The strict S-P-D+LO identity is provided only for absolute cross sections. '
                    'Normalized production-related shifts are differences of full prescriptions, not a homogeneous coefficient.')
            report['spectra'][layout['titles'][index]]=dict(write_spectrum('histogram_%d'%index,detail),precision=precision)
            del absolute,normalized,coverage,detail
        base_start=int(layout['offsets'][baseline[charge][0]])
        selected={key:value[:,[start+3,start+4,base_start+3,base_start+4]]
                  for key,value in vectors.items() if key[0]==charge}
        migration=estimate(selected,charges=(charge,),transform=lambda sums:
            prescription_observables(sums,charge,lambda row:ratio(row[:2],row[2:])))
        keep('migration_%d'%h,migration)
        report['migrations'][title]={label:summaries(migration,LABELS,i)
            for i,label in enumerate(('one_b_ratio_to_R04_b25','two_b_ratio_to_R04_b25'))}
        print('Reduced all main bins:',title,flush=True)
    for h,title,charge in rates:
        if charge!='plus':continue
        minus_title=title.replace('W+','W-')
        mh=layout['titles'].index(minus_title)
        starts={'plus':int(layout['offsets'][h]),'minus':int(layout['offsets'][mh])}
        selected={key:value[:,starts[key[0]]:starts[key[0]]+9] for key,value in vectors.items()}
        result=estimate(selected,transform=charge_observables)
        keep('charge_rates_%d'%h,result)
        report['charge_rates'][title.replace('W+','W+/W-')]={label:{variant:{name:dict(
            describe(result['value'][i,v,j]),nominal_retraining_mc_error=result['mc_error'][i,v,j,0])
            for i,name in enumerate(CHARGE_LABELS)} for v,variant in enumerate(VARIANTS)}
            for j,label in enumerate(RATES)}
        del result
        for local in range(2,22):
            indices={'plus':h+local-1,'minus':mh+local-1}
            selected={}
            for c in ('plus','minus'):
                a,b=map(int,layout['offsets'][indices[c]:indices[c]+2])
                parent=[starts[c]+i for i in normalization_bins(local)]
                selected.update(feature_vectors(vectors,c,a,b,parent))
            pa,pb=map(int,layout['offsets'][indices['plus']:indices['plus']+2])
            ma,mb=map(int,layout['offsets'][indices['minus']:indices['minus']+2])
            if not np.array_equal(layout['edges'][pa:pb],layout['edges'][ma:mb]):
                raise ValueError('Charge-comparison bin boundaries differ')
            detail=dict(title=layout['titles'][indices['plus']].replace('W+','W+/W-'),
                edges=layout['edges'][pa:pb],normalization_rate_bins=[RATES[i] for i in normalization_bins(local)],
                interpretation='The combined shape normalizes the summed cross section. Shape ratios/asymmetries '
                    'compare separately normalized charges and are distinct from absolute differential charge ratios.')
            for normalized in (False,True):
                result=estimate(selected if normalized else {key:row[:,:-1] for key,row in selected.items()},
                    transform=(lambda sums:charge_observables(sums,normalize_spectrum)) if normalized else charge_observables)
                labels=CHARGE_SHAPE_LABELS if normalized else CHARGE_LABELS
                prefix=('charge_normalized_' if normalized else 'charge_absolute_')+str(indices['plus'])
                keep(prefix,result)
                detail['normalized' if normalized else 'absolute']=[{variant:{name:dict(
                    describe(result['value'][i,v,j]),nominal_retraining_mc_error=result['mc_error'][i,v,j,0])
                    for i,name in enumerate(labels)} for v,variant in enumerate(VARIANTS)} for j in range(pb-pa)]
                del result
            report['charge_spectra'][detail['title']]=write_spectrum('charge_histogram_%d'%indices['plus'],detail)
            del detail
    for key,value in vectors.items():
        charge,variant,flavour=key
        start=int(layout['offsets'][baseline[charge][0]])
        rows=[start+i for i in (0,3,4)]
        n=len(value)
        between=np.var(value[:,rows,0],axis=0,ddof=1)/n
        within=conditional[key][:,rows].sum(axis=0)/n**2
        report['retraining_diagnostics']['__'.join((charge,variant,''.join(flavour)))]=dict(
            retrainings=n,rate_labels=['input','fiducial_1b','fiducial_2b'],
            between_retraining_variance_of_mean=between,
            mean_conditional_variance_of_combination=within,variance_ratio=ratio(between,within),
            interpretation='Diagnostic comparison with few retrainings; no automatic rescaling or variance addition.')
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),
        source_hashes={name:digest(STUDY/'scripts'/name) for name in ('main_report.py','main_observables.py',
            'main_replicas.py','full_flavour_statistics.py','replica_statistics.py','pilot_report.py','pdf_statistics.py',
            'main_binning.py','rebinning.py','covariance_storage.py')})
    report['source_hashes']['Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py']=digest(
        ROOT/'Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py')
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('replicas',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--binning',type=Path,help='Declared common bin-plan JSON, applied to whole-run vectors')
    args=parser.parse_args()
    run(args.replicas,args.output,args.binning)
