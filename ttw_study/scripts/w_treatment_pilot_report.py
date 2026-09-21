#!/usr/bin/env python3
"""Joint three-treatment W-width diagnostics from actual S/Pi pilot batches.

These are single ordered-flavour, conditional-on-grid comparisons. They do
not replace the full-flavour independently retrained robustness campaign.
"""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from compare_current_batches import all_stage_pairs, batch_record
from covariance_storage import nominal_copy
from inclusive_report import production_parameters
from joint_comparison import read
from joint_report import derived_rates, normalization_bins
from load_results import assert_matched_physics, assert_same_layout
from parameter_statistics import LABELS as PARAMETER_LABELS, contrast
from pilot_report import RATES, clean, describe
from read_splits import split_ensembles
from replica_statistics import independent_jackknife_many, ratio

MODES = ('onshell','top-bw','all-bw')
COMPARISONS = ('top_W','associated_W','all_W')
PAIRS = ((0,1),(1,2),(0,2))
LABELS = tuple(label.replace('scan','target') for label in PARAMETER_LABELS)
DERIVED = ('acceptance_1b','acceptance_2b','two_b_given_one_b','veto_0extra',
           'fraction_1extra','fraction_2extra','fraction_3plus_extra')
SPECTRA = (4,10,11,12,20)


def validate_physics(results):
    from models.check_param_card import ParamCard
    if len(results)!=6 or [row['report']['variant'] for row in results]!=['S','Pi']*3:
        raise ValueError('Require onshell, top-bw, all-bw S/Pi samples in that order')
    assert_same_layout(results)
    benchmark = json.loads((STUDY/'inputs/benchmark.json').read_text())
    inputs = benchmark['inputs']
    nominal = results[0]['report']['manifest']
    projected = []
    for i,result in enumerate(results):
        mode = MODES[i//2]
        manifest = result['report']['manifest']
        generation = result['report']['execution']['export_generation']['args']
        treatment = 'onshell' if mode=='onshell' else 'bw'
        width = next(row for row in benchmark['top_widths'] if row['mb_GeV']==0.
                     and row['w_treatment']==treatment and row['scale_factor']==1.)
        expected = dict(w_treatment=mode,decay_bottom_mass=0.,production_bottom_mass=0.,
            tops_on_shell=True,alpha_s_flavours=5,top_width_bottom_mass=0.,
            top_width_reference_scale=inputs['mt'],top_width_w_treatment=treatment,
            top_width_lo=width['gamma_lo'],top_width_nlo=width['gamma_nlo_pdf'],
            w_width=inputs['Gamma_W_fixed'],decay_scale_grouping='separate',
            production_scale='core-w-ht-half',production_scale_grouping='W_SYSTEM')
        if (any(manifest[key]!=value for key,value in expected.items())
                or generation['w_treatment']!=mode or generation['decay_bottom_mass']!=0.
                or generation['corrected']!='both' or manifest.get('custom_user_hook')
                or manifest.get('physical_scale_coordinates') or manifest.get('technical_test')):
            raise ValueError('W-treatment comparison has unmatched benchmark widths, masses or scales')
        settings = manifest['settings']
        if (any(settings[key] for key in ('fixed_ren_scale','fixed_fac_scale','fixed_qes_scale'))
                or any(settings[key]!=1. for key in ('mur_over_ref','muf_over_ref'))
                or settings.get('qes_over_ref',1.)!=1. or settings['dynamical_scale_choice']!=3
                or settings['lhaid']!=331700 or settings['ebeam1']!=6500. or settings['ebeam2']!=6500.):
            raise ValueError('Require common dynamic W-system scales and the benchmark PDF/energy')
        process = Path(result['report']['path']).parents[2]
        archive = process/'study_cards'/Path(result['report']['path']).parent.name
        card = ParamCard(str(archive/'param_card.dat'))
        alpha = next(row['alpha_s'] for row in benchmark['coupling_checkpoints'] if row['Q_GeV']==inputs['mz'])
        parameters = {('mass',6):inputs['mt'],('yukawa',6):inputs['mt'],('mass',24):inputs['mw'],
            ('mass',23):inputs['mz'],('mass',5):0.,('sminputs',1):inputs['alpha_Gmu_inverse'],
            ('sminputs',2):inputs['GF'],('sminputs',3):alpha,('decay',24):inputs['Gamma_W_fixed']}
        if 'decaymass' in card or any(not np.isclose(float(card[block].get((pdg,)).value),value,rtol=1.e-12,atol=0.)
                                    for (block,pdg),value in parameters.items()):
            raise ValueError('Archived numeric parameters differ from the W-treatment benchmark')
        for kind in ('topology_hashes','internal_width_hashes'):
            if not manifest[kind] or any(digest(process/path)!=sha for path,sha in manifest[kind].items()):
                raise ValueError('W-treatment topology/internal-width evidence changed')
        for relative in manifest['internal_width_hashes']:
            pdgs = json.loads((process/relative).read_text())['pdgs']
            if pdgs != ([] if mode=='onshell' else [24]):
                raise ValueError('Wrong W propagator-width prescription')
        for kind in ('scale_runtime_hashes','sampling_runtime_hashes'):
            if not manifest[kind] or any(digest(process/'SubProcesses'/path)!=sha for path,sha in manifest[kind].items()):
                raise ValueError('Generated scale/sampling sources changed')
        # Project only the W-treatment and its matched top-width change.
        # Numeric weak inputs and the actual fixed W width were checked above.
        row = dict(result,report=copy.deepcopy(result['report']),
                   parameter_signature=production_parameters(archive/'param_card.dat'))
        row['report']['execution']['export_generation']['args']['w_treatment']='onshell'
        for field in ('top_width_lo','top_width_nlo','top_width_w_treatment'):
            row['report']['manifest'][field]=nominal[field]
        projected.append(row)
    assert_matched_physics(projected)
    generation = results[0]['report']['execution']['export_generation']['args']
    return dict(charge=generation['charge'],flavours=generation['flavours'])


def estimate(vectors,strata,transform=None):
    if len(vectors)!=6 or len(strata)!=6:
        raise ValueError('Require six independent W-treatment/variant batch ensembles')
    groups = [split_ensembles(values,labels) for values,labels in zip(vectors,strata)]
    boundaries = np.cumsum([0]+[len(group) for group in groups])
    def combined(*means):
        rows = [sum(means[a:b]) for a,b in zip(boundaries[:-1],boundaries[1:])]
        if transform is not None:
            rows = [transform(row) for row in rows]
        return np.asarray([contrast(rows[2*reference],rows[2*reference+1],rows[2*target],rows[2*target+1])
                           for reference,target in PAIRS])
    return independent_jackknife_many([values for group in groups for values in group],combined)


def summaries(result,index):
    output = {}
    for c,name in enumerate(COMPARISONS):
        output[name] = {}
        for i,label in enumerate(LABELS):
            value,error,bias = (result[key][c,i,index] for key in ('value','mc_error','nonlinear_bias_estimate'))
            row = dict(describe(value),nominal_conditional_mc_error=clean(error[0]),
                       nominal_jackknife_bias_estimate=clean(bias[0]))
            if label in ('product_shift_change','product_ratio_change'):
                row.update(change_over_conditional_mc_error=clean(ratio(value[0],error[0])),
                    exceeds_three_conditional_mc_errors=bool(np.isfinite(error[0]) and error[0]>0. and abs(value[0])>3.*error[0]))
            output[name][label] = row
    return output


def run(paths,output,locals=SPECTRA):
    output = Path(output).resolve()
    paths = [Path(path).resolve() for path in paths]
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite W-treatment pilot reductions')
    if len(set(paths))!=6 or not locals or len(set(locals))!=len(locals) or any(i not in range(2,22) for i in locals):
        raise ValueError('Require six distinct batches and distinct spectrum IDs 2--21')
    results = []
    for path in paths:
        results.append(read(path))
        print('Loaded W-treatment batches:',len(results),'/ 6',flush=True)
    context = validate_physics(results)
    streams,counts = set(),[]
    for path,result in zip(paths,results):
        metadata,record = batch_record(path)
        pairs = all_stage_pairs(Path(result['report']['path']).parents[2],record)
        if len(pairs)!=metadata['refinement_rng_audit']['distinct_initialization_pairs'] or streams & pairs:
            raise ValueError('W-treatment training/refinement stream inventory overlaps or changed')
        streams.update(pairs); counts.append(len(pairs))
    strata = [row['strata'] for row in results]
    layout = {key:results[0][key] for key in ('titles','edges','offsets','weights')}
    arrays = dict(layout,comparison_labels=COMPARISONS,contrast_labels=LABELS)
    report = dict(created_utc=now(),status='conditional three-W-treatment pilot comparison; full-flavour retraining required',
        **context,inputs={str(path):digest(path) for path in paths},benchmark_sha256=digest(STUDY/'inputs/benchmark.json'),
        stage_stream_counts=counts,cross_sample_rng_overlap=0,rates={},derived_rates={},spectra={},
        comparison_definitions={name:dict(reference=MODES[r],target=MODES[t]) for name,(r,t) in zip(COMPARISONS,PAIRS)},
        convention='Three W comparisons share the six independent S/Pi samples. Joint batch deletion '
                    'retains shared-mode covariance, bin/parent correlations and matching scale/PDF coordinates. '
                    'All top widths match their W treatment; no extra BW branching factor is applied.',
        limitations='Single ordered flavour and conditional-on-grid MC errors. No main full-flavour result, '
                    'independent-retraining coverage, global significance or final W/mass robustness conclusion. '
                    'Native bins are retained; no continuous overflow folding or publication bin choice.')
    def keep(prefix,result):
        arrays[prefix]=result['value']; arrays[prefix+'_mc_errors']=result['mc_error']
        arrays[prefix+'_covariance_factor_nominal']=nominal_copy(result['covariance_factor'])
        arrays[prefix+'_nonlinear_bias_nominal']=nominal_copy(result['nonlinear_bias_estimate'])
    banks = []
    for h,title in enumerate(layout['titles']):
        if ' rates:' not in title or ('W+' if context['charge']=='plus' else 'W-') not in title:
            continue
        first = int(layout['offsets'][h]); config=title.split()[0]
        banks.append((h,first,config))
        vectors = [row['contributions'][:,first:first+9] for row in results]
        result = estimate(vectors,strata)
        keep(config+'_rates',result)
        report['rates'][config]={name:summaries(result,i) for i,name in enumerate(RATES)}
        result = estimate(vectors,strata,derived_rates)
        keep(config+'_derived',result)
        report['derived_rates'][config]={name:summaries(result,i) for i,name in enumerate(DERIVED)}
        print('Reduced W-treatment rates:',config,flush=True)
    if len(banks)!=5:
        raise ValueError('Require all five active-charge rate banks')
    vectors=[np.concatenate([row['contributions'][:,first:first+9,:1] for h,first,config in banks],axis=1)
             for row in results]
    result=estimate(vectors,strata)
    keep('joint_all_rates_nominal',result)
    report['joint_rate_feature_labels']=[config+'/'+name for h,first,config in banks for name in RATES]
    del result,vectors
    h,first,config=next(row for row in banks if row[2]=='R04_b25')
    for local in locals:
        a,b=map(int,layout['offsets'][h+local-1:h+local+1])
        parents=[first+i for i in normalization_bins(local)]
        vectors=[np.concatenate([row['contributions'][:,a:b],
                                 row['contributions'][:,parents].sum(axis=1)[:,None,:]],axis=1) for row in results]
        absolute=estimate(vectors,strata,lambda row:row[:-1])
        normalized=estimate(vectors,strata,lambda row:ratio(row[:-1],row[-1]))
        coverage=estimate(vectors,strata,lambda row:ratio(row[:-1],row[-1]).sum(axis=0)[None])
        prefix=config+'_h'+str(local)
        keep(prefix+'_absolute',absolute); keep(prefix+'_normalized',normalized)
        report['spectra'][prefix]=dict(title=layout['titles'][h+local-1],edges=layout['edges'][a:b],
            normalization_rate_bins=[RATES[i] for i in normalization_bins(local)],
            measured_range_coverage=summaries(coverage,0),
            nominal_normalized_product_ratio_changes={name:dict(value=normalized['value'][c,LABELS.index('product_ratio_change'),:,0],
                mc_error=normalized['mc_error'][c,LABELS.index('product_ratio_change'),:,0]) for c,name in enumerate(COMPARISONS)},
            pointwise_values_and_errors=prefix+'_absolute / '+prefix+'_normalized',
            units='Absolute bin integrals in pb; normalized quantities are bin fractions. Divide by bin width only for densities.')
        del absolute,normalized,coverage,vectors
        print('Reduced W-treatment spectrum:',prefix,flush=True)
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),
        source_hashes={name:digest(STUDY/'scripts'/name) for name in ('w_treatment_pilot_report.py','parameter_statistics.py',
            'joint_comparison.py','joint_report.py','load_results.py','compare_current_batches.py','read_splits.py',
            'rng_history.py','replica_statistics.py','covariance_storage.py','pilot_report.py','inclusive_report.py')})
    report['source_hashes']['Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py']=digest(ROOT/'Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py')
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples',nargs=6,type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--histograms',nargs='+',type=int,default=SPECTRA)
    args=parser.parse_args()
    run(args.samples,args.output,args.histograms)
