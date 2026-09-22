#!/usr/bin/env python3
"""Four independent S0,Pi0,Sm,Pim ensembles; form mass contrasts per weight."""
import argparse
import copy
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from campaign import STUDY, digest, now, save
from compare_current_batches import all_stage_pairs, batch_record
from inclusive_report import production_parameters
from joint_comparison import read
from joint_report import derived_rates, normalization_bins
from load_results import assert_matched_physics, assert_same_layout
from pilot_report import RATES, clean, describe
from production_invariance import compare_subprocess, model_functions
from read_splits import split_ensembles
from replica_statistics import independent_jackknife_many, ratio

LABELS=('S_mass_over_zero','Pi_mass_over_zero','Pi_minus_S_zero',
        'Pi_minus_S_mass','mass_change_in_Pi_minus_S','Pi_over_S_zero',
        'Pi_over_S_mass','mass_change_in_Pi_over_S','double_ratio')


def contrasts(s0,p0,sm,pm):
    r0,rm=ratio(p0,s0),ratio(pm,sm)
    return np.asarray([ratio(sm,s0),ratio(pm,p0),p0-s0,pm-sm,
                       (pm-sm)-(p0-s0),r0,rm,rm-r0,ratio(rm,r0)])


def estimate(vectors,strata,transform=contrasts):
    if len(vectors)!=4 or len(strata)!=4:
        raise ValueError('Exactly four independent predictions are required')
    groups=[split_ensembles(v,s) for v,s in zip(vectors,strata)]
    boundaries=np.cumsum([0]+[len(group) for group in groups])
    def combined(*means):
        return transform(*[sum(means[a:b]) for a,b in zip(boundaries[:-1],boundaries[1:])])
    return independent_jackknife_many([row for group in groups for row in group],combined)


def validate_physics(results,identities,*,width_inputs=None):
    from models.check_param_card import ParamCard
    if [r['report']['variant'] for r in results]!=['S','Pi','S','Pi']:
        raise ValueError('Require S0,Pi0,Sm,Pim in that order')
    assert_same_layout(results)
    assert_matched_physics(results[:2])
    assert_matched_physics(results[2:])
    manifests=[r['report']['manifest'] for r in results]
    masses=[m['decay_bottom_mass'] for m in manifests]
    if masses[:2]!=[0.,0.] or masses[2]<=0. or masses[2]!=masses[3]:
        raise ValueError('Require matched zero/positive decay masses')
    benchmark=json.loads((STUDY/'inputs/benchmark.json').read_text())
    width_rows=benchmark['top_widths']
    if width_inputs is not None:
        from small_mass_inputs import validate
        supplied=json.loads(width_inputs.read_text())
        if supplied['benchmark_sha256']!=digest(STUDY/'inputs/benchmark.json'):
            raise ValueError('Custom mass widths use a different benchmark')
        validate(supplied,benchmark)
        if any(m.get('small_mass_inputs_sha256')!=digest(width_inputs) for m in manifests):
            raise ValueError('Custom mass widths were not recorded by every generated sample')
        width_rows=supplied['top_widths']
    projected=[]
    processes=[]
    for result,manifest,mass in zip(results,manifests,masses):
        process=Path(result['report']['path']).parents[2]
        processes.append(process)
        archive=process/'study_cards'/Path(result['report']['path']).parent.name
        card=ParamCard(str(archive/'param_card.dat'))
        actual_mass=card['decaymass'].get((5,)).value if 'decaymass' in card else 0.
        if actual_mass!=mass or card['mass'].get((5,)).value!=0.:
            raise ValueError('Archived production/decay masses differ from the requested scheme')
        if card['decay'].get((24,)).value!=manifest['w_width']:
            raise ValueError('Archived internal W width differs from its manifest')
        mode=manifest['top_width_w_treatment']
        matches=[r for r in width_rows if r['mb_GeV']==mass and
                 r['w_treatment']==mode and r['scale_factor']==1.]
        if len(matches)!=1 or manifest['top_width_reference_scale']!=172.5:
            raise ValueError('Mass comparison requires the matched benchmark width convention')
        for field,expected in [('top_width_lo',matches[0]['gamma_lo']),
                               ('top_width_nlo',matches[0]['gamma_nlo_pdf'])]:
            if not np.isclose(manifest[field],expected,rtol=1.e-12,atol=0.):
                raise ValueError('Unmatched massive/zero top width')
        # Compare all other settings without mutating the loaded evidence.
        projected_result=dict(result,report=copy.deepcopy(result['report']),
                              parameter_signature=production_parameters(archive/'param_card.dat'))
        m=projected_result['report']['manifest']
        for field in ('top_width_lo','top_width_nlo','top_width_bottom_mass'):
            m[field]=manifests[0][field]
        projected_result['report']['execution']['export_generation']['args']['decay_bottom_mass']=0.
        projected.append(projected_result)
    # Gamma_W was verified against each archived card above; compare it
    # through the manifest along with the other non-mass settings.
    assert_matched_physics(projected)
    covered=set()
    for path in identities:
        proof=json.loads(path.read_text())
        reference,candidate=Path(proof['reference']),Path(proof['candidate'])
        if (proof['status'] not in ('production source and FKS-region identity passed',
                'production amplitude and FKS-region identity passed with kinematic-helper transition') or
                reference not in processes[:2] or candidate not in processes[2:]):
            raise ValueError('Production identity does not cover these samples')
        if digest(candidate/'Cards/decay_mass_scheme.json')!=proof['scheme_sha256']:
            raise ValueError('Mass scheme changed since production audit')
        for process in (reference,candidate):
            production_getter,_=model_functions(process/'Source/MODEL/get_mass_width_fcts.f')
            if hashlib.sha256(production_getter.encode()).hexdigest()!=proof['production_getter_sha256']:
                raise ValueError('Production mass/width getter changed since its identity audit')
        for row in proof['subprocesses'].values():
            current=compare_subprocess(Path(row['reference_directory']),Path(row['candidate_directory']),
                allow_kinematic_update=bool(row.get('kinematic_matrix_transition')))
            if any(current[key]!=row[key] for key in current):
                raise ValueError('Production source/FKS evidence changed')
        covered.add(candidate)
    if covered!=set(processes[2:]):
        raise ValueError('Missing unchanged-production evidence for a massive sample')
    return masses[2]


def run(paths,identities,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Never overwrite a mass comparison')
    results=[read(p) for p in paths]
    mass=validate_physics(results,identities)
    streams=set()
    stream_counts=[]
    for path,result in zip(paths,results):
        metadata,record=batch_record(path)
        pairs=all_stage_pairs(Path(result['report']['path']).parents[2],record)
        if pairs & streams:
            raise ValueError('Cross-sample training/refinement RNG overlap')
        streams.update(pairs)
        stream_counts.append(len(pairs))
    strata=[r['strata'] for r in results]
    shape=(len(LABELS),)+results[0]['values'].shape
    arrays={name:np.full(shape,np.nan) for name in
            ('absolute','absolute_mc_errors','normalized','normalized_mc_errors')}
    report=dict(created_utc=now(),status='conditional four-sample mass pilot; convergence required',
        inputs={str(p.resolve()):digest(p) for p in paths+identities},decay_mass_GeV=mass,
        stage_stream_counts=stream_counts,cross_sample_rng_overlap=0,rates={},derived_rates={},spectra={},
        convention='Four independent predictions and independent beam strata. Ratios of sums, '
            'not sums of ratios; paired numerator/denominator bins within each worker. '
            'All 81 scales and 101 PDFs transformed before uncertainty reductions.',
        limitations='One ordered flavour pilot, conditional on learned grids; not five retrainings, '
            'tail coverage, small-mass continuity, full-flavour normalization or publication precision.')
    def summaries(row):
        return {label:dict(describe(row['value'][i]),central_joint_mc_error=row['mc_error'][i,0])
                for i,label in enumerate(LABELS)}
    for h,title in enumerate(results[0]['titles']):
        if ' rates:' not in title:
            continue
        first=int(results[0]['offsets'][h])
        if not any(np.any(r['values'][first:first+9,0]) for r in results):
            continue
        vectors=[r['contributions'][:,first:first+9,:] for r in results]
        rates=estimate(vectors,strata)
        for name,key in [('absolute','value'),('absolute_mc_errors','mc_error')]:
            arrays[name][:,first:first+9]=rates[key]
        report['rates'][title]={label:summaries({key:rates[key][:,i] for key in ('value','mc_error')})
                                for i,label in enumerate(RATES)}
        derived=estimate(vectors,strata,lambda *rows:contrasts(*map(derived_rates,rows)))
        report['derived_rates'][title]={label:summaries({key:derived[key][:,i] for key in ('value','mc_error')})
            for i,label in enumerate(('acceptance_1b','acceptance_2b','two_b_given_one_b',
                                     'veto_0extra','fraction_1extra','fraction_2extra','fraction_3plus_extra'))}
        for local_id in range(2,22):
            index=h+local_id-1
            start,stop=map(int,results[0]['offsets'][index:index+2])
            vectors=[]
            for r in results:
                denominator=r['contributions'][:,[first+i for i in normalization_bins(local_id)],:].sum(axis=1)
                vectors.append(np.concatenate([r['contributions'][:,start:stop,:],denominator[:,None,:]],axis=1))
            absolute=estimate(vectors,strata,lambda *rows:contrasts(*[row[:-1] for row in rows]))
            normalized=estimate(vectors,strata,lambda *rows:contrasts(*[ratio(row[:-1],row[-1]) for row in rows]))
            for prefix,row in [('absolute',absolute),('normalized',normalized)]:
                arrays[prefix][:,start:stop]=row['value']
                arrays[prefix+'_mc_errors'][:,start:stop]=row['mc_error']
            report['spectra'][results[0]['titles'][index]]=dict(
                normalization_rate_bins=[RATES[i] for i in normalization_bins(local_id)],
                undefined_double_ratio_mc_errors=int(np.sum(~np.isfinite(normalized['mc_error'][-1,:,0]))))
        print('Mass comparison: '+title,flush=True)
    np.savez_compressed(output.with_suffix('.npz'),**arrays,labels=LABELS,
        edges=results[0]['edges'],offsets=results[0]['offsets'],titles=results[0]['titles'])
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),script_sha256=digest(Path(__file__)))
    save(output,clean(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples',nargs=4,type=Path,required=True,metavar=('S0','PI0','SM','PIM'))
    parser.add_argument('--production-identity',action='append',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--after-pilot',type=Path)
    args=parser.parse_args()
    if args.after_pilot:
        while True:
            pilot=json.loads(args.after_pilot.read_text())
            if pilot['status']=='stopped':
                raise ValueError('Preceding massive pilot failed; no mass comparison')
            if pilot['status'].startswith('completed technical pilot;'):
                if Path(pilot['batches']).resolve()!=args.samples[3].resolve():
                    raise ValueError('Pilot gate and massive product sample differ')
                break
            time.sleep(15)
    run(args.samples,args.production_identity,args.output)
