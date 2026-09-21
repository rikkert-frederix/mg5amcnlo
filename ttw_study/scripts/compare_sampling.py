#!/usr/bin/env python3
"""Compare same-physics, same-prescription sampling proposals without certification."""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import digest, now, save
from joint_comparison import paired_variants, read
from load_results import assert_matched_physics, assert_same_layout
from pilot_report import RATES, clean
from replica_statistics import ratio


def compare(first_path,second_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a sampling comparison')
    first,second=read(first_path),read(second_path)
    assert_same_layout([first,second])
    assert_matched_physics([first,second])
    manifests=[row['report']['manifest'] for row in (first,second)]
    if manifests[0]['variant']!=manifests[1]['variant']:
        raise ValueError('Sampling controls must use the same perturbative prescription')
    if manifests[0].get('sampling_runtime_hashes')!=manifests[1].get('sampling_runtime_hashes'):
        raise ValueError('Sampling-map control must use the same corrected source kernel')
    if set(map(tuple,first['seeds'])) & set(map(tuple,second['seeds'])):
        raise ValueError('Independent comparison has repeated random streams')
    report=dict(created_utc=now(),status='sampling control comparison; agreement and convergence not certified',
                inputs={str(p.resolve()):digest(p) for p in (first_path,second_path)},
                controls=[dict(run=m['run_name'],variant=m['variant'],
                               sampling=m.get('production_sampling','flat'),
                               proposal=m.get('production_sampling_parameters')) for m in manifests],
                convention='Independent variant/stratum batch means; conditional on trained grids. '
                           'Correlated scale/PDF weights remain paired within each run. '
                           'Large-tail nonconvergence can invalidate a Gaussian interpretation of pulls.',
                rates={},acceptances={},batch_concentration={})
    def estimate(a,b,transform):
        return paired_variants(a,b,first['strata'],second['strata'],transform)
    def values(a,b):
        return np.asarray([a,b,b-a,ratio(b,a)])
    arrays={}
    for h,title in enumerate(first['titles']):
        if ' rates:' not in title:
            continue
        start=int(first['offsets'][h])
        if not np.any(first['values'][start:start+9,0]):
            continue
        a,b=[row['contributions'][:,start:start+9,:] for row in (first,second)]
        rate=estimate(a,b,values)
        report['rates'][title]={label:dict(
            first_pb=rate['value'][0,i],second_pb=rate['value'][1,i],
            second_minus_first_pb=rate['value'][2,i],
            second_minus_first_mc_error_pb=rate['mc_error'][2,i],
            nominal_difference_pull=ratio(rate['value'][2,i,0],rate['mc_error'][2,i,0]),
            second_over_first=rate['value'][3,i],second_over_first_mc_error=rate['mc_error'][3,i])
            for i,label in enumerate(RATES)}
        accept=estimate(a,b,lambda x,y:values(ratio(x[[3,4]],x[0]),ratio(y[[3,4]],y[0])))
        report['acceptances'][title]=dict(values=accept['value'],mc_errors=accept['mc_error'],
                                         labels=['first','second','second_minus_first','second_over_first'])
        arrays['rates_'+str(h)]=rate['value']
        arrays['rate_errors_'+str(h)]=rate['mc_error']
        arrays['acceptances_'+str(h)]=accept['value']
        arrays['acceptance_errors_'+str(h)]=accept['mc_error']
    for label,row in zip(('first','second'),(first,second)):
        h=next(i for i,title in enumerate(row['titles']) if ' rates:' in title
               and row['values'][row['offsets'][i],0]!=0.)
        samples=row['contributions'][:,int(row['offsets'][h]),0]
        strata=[]
        for s in sorted(set(row['strata'])):
            x=samples[row['strata']==s]
            strata.append(dict(stratum=int(s),batches=len(x),sum_pb=sum(x),
                               minimum_contribution_pb=min(x),maximum_contribution_pb=max(x),
                               largest_absolute_contribution_fraction=max(abs(x))/sum(abs(x))))
        report['batch_concentration'][label]=strata
    np.savez_compressed(output.with_suffix('.npz'),**arrays)
    report.update(arrays_sha256=digest(output.with_suffix('.npz')),script_sha256=digest(Path(__file__)))
    save(output,clean(report))
    for title,rows in report['rates'].items():
        if title.startswith('R04_b25'):
            print(title,flush=True)
            for key in ('input','fiducial_1b','fiducial_2b'):
                r=rows[key]
                print(key,'difference',r['second_minus_first_pb'][0],
                      '+/-',r['second_minus_first_mc_error_pb'][0],
                      'conditional pull',r['nominal_difference_pull'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',required=True,type=Path)
    parser.add_argument('--second',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    compare(args.first,args.second,args.output)
