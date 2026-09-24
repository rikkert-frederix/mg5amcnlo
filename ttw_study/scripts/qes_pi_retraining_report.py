#!/usr/bin/env python3
"""Compare independent Pi QES retrainings and the earlier conditional pair."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from campaign import STUDY, digest, now, save
from joint_comparison import paired_variants, read as load_batches
from load_results import assert_same_layout
from pilot_report import RATES, clean
from production_invariance import compare_subprocess
from replica_statistics import ratio
from run_qes_pi_retraining import DONE
from run_scale_validation import assert_technical_match
from run_width_mass_pilots import verified_pairs

SPECTRA=(3,4,6,11,12,20)


def pair_values(half,two):
    return np.asarray([half,two,half-two,ratio(half,two)])


def conditional_contrast(new,old):
    difference=new['difference_pb']-old['difference_pb']
    error=math.hypot(new['conditional_mc_error_pb'],old['mc_error_pb'])
    return dict(old=old,new=new,between_pair_difference_pb=difference,
                conditional_mc_error_pb=error,
                between_pair_conditional_pull=difference/error if error>0 else None)


def run(queue_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a QES retraining comparison')
    queue=json.loads(queue_path.read_text())
    jobs=queue['jobs']
    if (queue['status']!=DONE or len(jobs)!=2 or
            [(job['case']['name'],job['case']['qes'],job['variant'],job['seed']) for job in jobs]!=
            [('qes_half',.5,'Pi',93001),('qes_two',2.,'Pi',93002)] or
            jobs[0]['process']!=jobs[1]['process'] or
            digest(queue['source_identity'])!=queue['source_identity_sha256']):
        raise ValueError('Incomplete or changed independent Pi QES pair')
    identity=json.loads(Path(queue['source_identity']).read_text())
    for row in identity['subprocesses'].values():
        current=compare_subprocess(Path(row['reference_directory']),Path(row['candidate_directory']),
                                   allow_kinematic_update=True)
        if any(current[key]!=row[key] for key in current):
            raise ValueError('QES production source transition changed')
    streams=set()
    for job in jobs:
        pairs=verified_pairs(job)
        if streams & pairs:
            raise ValueError('QES retrainings share actual stage streams')
        streams.update(pairs)
        virtual=json.loads(Path(job['analytic_virtual_checks']).read_text())
        if (virtual['status']!='all archived top/antitop analytic virtual checks passed' or
                virtual['maximum_relative_difference']>virtual['tolerance'] or
                virtual['tolerance']!=1.e-8 or digest(virtual['input'])!=virtual['input_sha256']):
            raise ValueError('Incomplete QES virtual validation')
    first,second=[load_batches(Path(job['batches'])) for job in jobs]
    assert_same_layout([first,second])
    if (first['contributions'].shape[-1]!=183 or
            second['contributions'].shape[-1]!=183):
        raise ValueError('QES pair lacks nominal, scale or PDF weights')
    reference=load_batches(Path(queue['reference_batches']))
    for result,job in zip((first,second),jobs):
        assert_technical_match(result,reference,job['case'])
    old_path=STUDY/'results/technical_validation_tech_v4_qes_two_vs_half_Pi_v1.json'
    old_influence_path=STUDY/'results/technical_validation_tech_v4_Pi_rate_influence_inspection_v1.json'
    old=json.loads(old_path.read_text())
    old_influence=json.loads(old_influence_path.read_text())
    if (old['cross_run_stage_collisions']!=0 or
            old_influence['status']!='Pi QES raw-rate influence inspected; independent retraining required before a QES-cancellation conclusion'):
        raise ValueError('Original Pi QES evidence is incomplete')
    inputs={str(p.resolve()):digest(p) for p in
            [queue_path,Path(queue['source_identity']),Path(queue['reference_batches']),
             old_path,old_influence_path]+[Path(job[field]) for job in jobs
              for field in ('batches','audit','analytic_virtual_checks','comparison')]}
    report=dict(created_utc=now(),
        status='independent Pi QES pair evaluated; inspect cancellation and retraining consistency',
        inputs=inputs,source_transition=queue['source_identity'],
        new_stage_stream_counts=[len(verified_pairs(job)) for job in jobs],
        new_pair_stage_overlap=0,old_pair_stage_overlap=old['cross_run_stage_collisions'],
        convention='Use matching physical scale and PDF weight coordinates. Each direct run has a fresh grid; '
            'delete one conditional batch within one independent run at a time. The common central reference '
            'cancels in the direct half-minus-two contrast.',
        rates={},spectra={},original_pair_comparison={})
    arrays={}
    h=next(i for i,title in enumerate(first['titles']) if title.startswith('R04_b25')
           and 'W+' in title and ' rates:' in title)
    start=int(first['offsets'][h])
    rates=paired_variants(first['contributions'][:,start:start+9,:],
        second['contributions'][:,start:start+9,:],first['strata'],second['strata'],pair_values)
    arrays['rate_values']=rates['value']
    arrays['rate_mc_errors']=rates['mc_error']
    for i,name in enumerate(RATES):
        value,error=rates['value'][:,i,0],rates['mc_error'][:,i,0]
        report['rates'][name]=dict(half_pb=value[0],two_pb=value[1],
            difference_pb=value[2],conditional_mc_error_pb=error[2],
            conditional_pull=value[2]/error[2] if error[2]>0 else None,
            half_over_two=value[3],ratio_conditional_mc_error=error[3])
    rejected=[row['contributions'][:,start:start+1,:]-
              row['contributions'][:,start+1:start+2,:] for row in (first,second)]
    rejection=paired_variants(*rejected,first['strata'],second['strata'],pair_values)
    arrays['rejected_values']=rejection['value']
    arrays['rejected_mc_errors']=rejection['mc_error']
    value,error=rejection['value'][:,0,0],rejection['mc_error'][:,0,0]
    report['rates']['rejected_by_lepton_selection']=dict(
        half_pb=value[0],two_pb=value[1],difference_pb=value[2],
        conditional_mc_error_pb=error[2],
        conditional_pull=value[2]/error[2] if error[2]>0 else None)
    for local in SPECTRA:
        index=h+local-1
        lo,hi=map(int,first['offsets'][index:index+2])
        parent=start+(4 if local in (11,12) else 3)
        vectors=[np.concatenate([row['contributions'][:,lo:hi,:],
                                 row['contributions'][:,parent:parent+1,:]],axis=1)
                 for row in (first,second)]
        normalized=paired_variants(*vectors,first['strata'],second['strata'],
            lambda half,two:pair_values(ratio(half[:-1],half[-1]),ratio(two[:-1],two[-1])))
        pulls=ratio(normalized['value'][2,:,0],normalized['mc_error'][2,:,0])
        finite=np.isfinite(pulls)
        report['spectra'][first['titles'][index]]=dict(
            nominal_populated_bins=int(finite.sum()),
            maximum_absolute_nominal_pull=float(np.max(np.abs(pulls[finite]))) if finite.any() else None,
            nominal_bins_above_three=int(np.sum(np.abs(pulls[finite])>3)))
        arrays['normalized_%d'%local]=normalized['value']
        arrays['normalized_mc_errors_%d'%local]=normalized['mc_error']
        arrays['edges_%d'%local]=first['edges'][lo:hi]
    for name,original in (
            ('input',old_influence['observations']['input_nominal']),
            ('leptons',old_influence['observations']['accepted_leptons_nominal']),
            ('rejected_by_lepton_selection',old_influence['observations']['rejected_by_lepton_selection_nominal'])):
        current=report['rates'][name]
        report['original_pair_comparison'][name]=conditional_contrast(current,
            dict(difference_pb=original['difference_pb'],mc_error_pb=original['mc_error_pb']))
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),
                  script_sha256=digest(Path(__file__)),
                  limitations='Conditional trained-grid errors do not prove coverage of between-retraining variance. '
                    'The old and new pairs share physical inputs but use separately audited generic loop-helper sources. '
                    'No global significance is inferred from correlated scale points or spectrum bins.')
    save(output,clean(report))
    print(output,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('queue',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.queue,args.output)
