#!/usr/bin/env python3
"""All-BW S/Pi central-production-scale diagnostics with shared controls."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from campaign import ROOT, STUDY, digest, now, save
from joint_report import derived_rates
from load_results import assert_matched_physics, assert_same_layout, load
from mass_effects import contrasts
from narrow_w_report import SPECTRA, selected_vectors
from pilot_report import RATES, clean, describe
from read_splits import split_ensembles
from replica_statistics import independent_jackknife_many, ratio
from run_scale_validation import DONE, cases, physical_scale_coordinates
from run_width_mass_pilots import verified_pairs

SAMPLES=('S_grouped','Pi_grouped','S_native','Pi_native','S_fixed','Pi_fixed')
CONTRASTS=('S_diagnostic_over_grouped','Pi_diagnostic_over_grouped',
           'Pi_minus_S_grouped','Pi_minus_S_diagnostic','change_in_Pi_minus_S',
           'Pi_over_S_grouped','Pi_over_S_diagnostic','change_in_Pi_over_S','double_ratio')
LABELS=SAMPLES+tuple(mode+'__'+label for mode in ('native','fixed') for label in CONTRASTS)
DERIVED=('acceptance_1b','acceptance_2b','two_b_given_one_b','veto_0extra',
         'fraction_1extra','fraction_2extra','fraction_3plus_extra')


def observables(rows):
    if len(rows)!=6:
        raise ValueError('Require grouped, native and fixed S/Pi samples in that order')
    return np.concatenate([np.asarray(rows),contrasts(*rows[:2],*rows[2:4]),
                           contrasts(*rows[:2],*rows[4:6])])


def estimate(vectors,strata,transform=None):
    if len(vectors)!=6 or len(strata)!=6:
        raise ValueError('Require six independent sample/stratum ensembles')
    groups=[split_ensembles(vector,labels) for vector,labels in zip(vectors,strata)]
    boundaries=np.cumsum([0]+[len(group) for group in groups])
    def combined(*means):
        rows=[sum(means[a:b]) for a,b in zip(boundaries[:-1],boundaries[1:])]
        return observables(rows if transform is None else [transform(row) for row in rows])
    return independent_jackknife_many([row for group in groups for row in group],combined)


def validate_physics(results):
    if len(results)!=6 or [r['report']['variant'] for r in results]!=['S','Pi']*3:
        raise ValueError('Require grouped, native and fixed S/Pi samples in that order')
    assert_same_layout(results)
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    width=next(row for row in benchmark['top_widths'] if row['mb_GeV']==0.
               and row['scale_factor']==1. and row['w_treatment']=='bw')
    reference=results[0]['report']['manifest']
    projected=[]
    choices=[None,None,cases()[-2],cases()[-2],cases()[-1],cases()[-1]]
    for result,choice in zip(results,choices):
        manifest=result['report']['manifest']
        generation=result['report']['execution']['export_generation']['args']
        scale=choice['production_scale'] if choice else 'core-w-ht-half'
        if any(generation[key]!=value for key,value in dict(charge='plus',flavours=['e','e','mu'],
                w_treatment='all-bw',decay_bottom_mass=0.,corrected='both').items()):
            raise ValueError('Require the same ordered massless-decay all-BW W+ sample')
        expected=dict(production_scale=scale,production_scale_grouping='W_SYSTEM' if choice is None else 'NONE',
            w_treatment='all-bw',decay_bottom_mass=0.,decay_scale_grouping='separate',
            top_width_reference_scale=172.5,top_width_w_treatment='bw',top_width_bottom_mass=0.,
            top_width_lo=width['gamma_lo'],top_width_nlo=width['gamma_nlo_pdf'],
            w_width=benchmark['W_width']['gamma_nlo_pdf'])
        if any(manifest[key]!=value for key,value in expected.items()) or manifest.get('custom_user_hook'):
            raise ValueError('Central-scale diagnostic has unmatched widths, scale grouping or decay scales')
        settings=manifest['settings']
        if any(settings[key]!=(scale=='fixed') for key in ('fixed_ren_scale','fixed_fac_scale','fixed_qes_scale')):
            raise ValueError('Declared central scale differs from the actual fixed/dynamic settings')
        if (any(settings[key]!=212.6925 for key in ('mur_ref_fixed','muf_ref_fixed','qes_ref_fixed'))
                or settings['dynamical_scale_choice']!=3 or settings.get('qes_over_ref',1.)!=1.
                or settings['mur_over_ref']!=1. or settings['muf_over_ref']!=1.):
            raise ValueError('Wrong production reference scales or a technical QES/scale variation')
        if choice is None:
            if manifest.get('technical_test') or manifest.get('physical_scale_coordinates'):
                raise ValueError('Grouped reference must use the unchanged central configuration')
        elif (manifest.get('technical_test')!=choice or
              manifest.get('physical_scale_coordinates')!=physical_scale_coordinates(choice)):
            raise ValueError('Diagnostic scale coordinates do not match its declared case')
        # Only the explicitly checked production-scale definition is projected
        # out of the general physics comparison. All other inputs must match.
        row=dict(result,report=copy.deepcopy(result['report']))
        target=row['report']['manifest']
        target['production_scale']=reference['production_scale']
        target.pop('physical_scale_coordinates',None)
        target['settings'].setdefault('qes_over_ref',1.)
        for key in ('fixed_ren_scale','fixed_fac_scale','fixed_qes_scale'):
            target['settings'][key]=reference['settings'][key]
        projected.append(row)
    assert_matched_physics(projected)


def summaries(result,bin_index):
    return {label:dict(describe(result['value'][i,bin_index]),
                nominal_conditional_mc_error=result['mc_error'][i,bin_index,0],
                nominal_nonlinear_bias_estimate=result['nonlinear_bias_estimate'][i,bin_index,0])
            for i,label in enumerate(LABELS)}


def run(queue_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a central-scale comparison')
    queue=json.loads(queue_path.read_text())
    if queue['status']!=DONE or len(queue['jobs'])!=12 or queue['cases']!=cases():
        raise ValueError('Require the complete twelve-case technical queue')
    if any(digest(ROOT/path)!=sha for path,sha in queue['source_hashes'].items()):
        raise ValueError('Frozen technical source or physical benchmark changed')
    controls_path=Path(queue['validation_queue'])
    if digest(controls_path)!=queue['reference_queues_sha256'][str(controls_path)]:
        raise ValueError('Central-scale reference queue changed')
    controls=json.loads(controls_path.read_text())
    indexed={}
    for job in controls['jobs']:
        audit=json.loads(Path(job['audit']).read_text())
        if audit['manifest']['w_treatment']=='all-bw':
            if ('grouped',job['variant']) in indexed:
                raise ValueError('Duplicated grouped-scale reference')
            indexed['grouped',job['variant']]=job
    for job in queue['jobs']:
        audit=json.loads(Path(job['audit']).read_text())
        case=audit['manifest']['technical_test']
        if case['name'] in ('central_native','central_fixed'):
            key=case['name'].removeprefix('central_'),job['variant']
            if key in indexed:
                raise ValueError('Duplicated central-scale sample')
            indexed[key]=job
    expected=[(mode,variant) for mode in ('grouped','native','fixed') for variant in ('S','Pi')]
    if set(indexed)!=set(expected):
        raise ValueError('Missing a grouped/native/fixed central-scale sample')
    jobs=[indexed[key] for key in expected]
    results=[load(job['audit']) for job in jobs]
    validate_physics(results)
    inputs={str(p.resolve()):digest(p) for p in (queue_path,controls_path)}
    selected,counts,streams=[],[],set()
    for job,result in zip(jobs,results):
        pairs=verified_pairs(job)
        if streams & pairs:
            raise ValueError('Central-scale diagnostic samples share actual stage streams')
        streams.update(pairs)
        counts.append(len(pairs))
        for field in ('audit','batches','analytic_virtual_checks'):
            inputs[str(Path(job[field]).resolve())]=digest(job[field])
        virtual=json.loads(Path(job['analytic_virtual_checks']).read_text())
        metadata=json.loads(Path(job['batches']).with_suffix('.json').read_text())
        worker=json.loads(Path(metadata['input']).read_text())
        if (virtual['status']!='all archived top/antitop analytic virtual checks passed' or
                virtual['maximum_relative_difference']>1.e-8 or
                virtual['tolerance']!=1.e-8 or virtual['checked_points']<=0 or
                Path(virtual['input']).resolve()!=Path(metadata['input']).resolve() or
                virtual['input_sha256']!=metadata['input_sha256'] or
                Path(job['audit']).resolve()!=Path(worker['audit']).resolve()):
            raise ValueError('Central-scale analytic virtual checks did not pass')
        selected.append(selected_vectors(job,result))
    strata=[row['strata'] for row in selected]
    arrays=dict(labels=LABELS,sample_labels=SAMPLES)
    report=dict(created_utc=now(),status='central-production-scale diagnostics evaluated; inspect physical robustness',
        inputs=inputs,samples=SAMPLES,labels=LABELS,stage_stream_counts=counts,
        distinct_initialization_pairs=len(streams),cross_sample_rng_overlap=0,
        rates={},derived_rates={},spectra={},
        production_scales=dict(grouped='CORE W-system HT/2',native='CORE lepton-resolved HT/2',
                               fixed_GeV=212.6925,top_and_antitop_reference_GeV=172.5),
        convention='Six independent samples retain the same grouped S/Pi controls across native and fixed '
            'diagnostics, including their shared covariance. Compare corresponding relative multiplier/PDF '
            'coordinates before reducing envelopes. Normalize every spectrum by its matched fiducial rate.',
        limitations='Different central production-scale definitions imply different physical predictions, '
            'not an identity or reweighting cancellation. These are one-assignment all-BW pilots with '
            'conditional learned-grid covariance. They do not replace full-flavour, retraining, or '
            'publication-bin convergence. Coarse diagnostic bins remain separate from publication choices.')
    def retain(prefix,result):
        arrays[prefix]=result['value']
        arrays[prefix+'_mc_errors']=result['mc_error']
        arrays[prefix+'_covariance_factor_nominal']=result['covariance_factor'][...,0]
        arrays[prefix+'_nonlinear_bias_estimate']=result['nonlinear_bias_estimate']
    for index,title in enumerate(selected[0]['rates']):
        vectors=[row['rates'][title] for row in selected]
        result=estimate(vectors,strata)
        retain('rates_%d'%index,result)
        report['rates'][title]={label:summaries(result,i) for i,label in enumerate(RATES)}
        result=estimate(vectors,strata,derived_rates)
        retain('derived_%d'%index,result)
        report['derived_rates'][title]={label:summaries(result,i) for i,label in enumerate(DERIVED)}
        del result
    for local in SPECTRA:
        vectors=[row['spectra'][local] for row in selected]
        entry=dict(title=results[0]['titles'][local-1],edges=selected[0]['edges'][local],
                   parent_region='2b' if local in (11,12) else '1b')
        for normalized in (False,True):
            transform=(lambda row:ratio(row[:-1],row[-1])) if normalized else (lambda row:row[:-1])
            result=estimate(vectors,strata,transform)
            label='normalized' if normalized else 'absolute'
            retain(label+'_'+str(local),result)
            entry[label]=[summaries(result,i) for i in range(result['value'].shape[1])]
            del result
        report['spectra'][str(local)]=entry
        arrays['edges_'+str(local)]=selected[0]['edges'][local]
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),
        source_hashes={str(p):digest(p) for p in [Path(__file__).resolve(),
            *[STUDY/'scripts'/name for name in ('mass_effects.py','narrow_w_report.py','load_results.py',
                'read_splits.py','replica_statistics.py','pilot_report.py','joint_report.py',
                'run_scale_validation.py','run_width_mass_pilots.py','rebinning.py')]]})
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.queue,args.output)
