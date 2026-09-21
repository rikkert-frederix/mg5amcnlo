#!/usr/bin/env python3
"""Conditional joint-batch S/Pi comparisons, preserving bin/rate correlations."""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import digest,now,save
from joint_report import derived_rates,normalization_bins
from load_results import assert_matched_physics,assert_same_layout,load
from pilot_report import RATES,clean,describe
from read_splits import split_ensembles
from replica_statistics import independent_jackknife_many,ratio


def read(path):
    metadata = json.loads(path.with_suffix('.json').read_text())
    if digest(path) != metadata['arrays_sha256'] or not metadata['status'].startswith('equal-count conditional'):
        raise ValueError('Require an unchanged successful batch audit')
    source = Path(metadata['input'])
    if digest(source) != metadata['input_sha256']:
        raise ValueError('Batch source record changed')
    record = json.loads(source.read_text())
    if digest(record['audit']) != record['audit_sha256']:
        raise ValueError('Combined-output audit changed')
    result = load(record['audit'])
    with np.load(path) as data:
        result.update(contributions=data['contributions'].copy(),strata=data['strata'].copy(),
                      seeds=data['random_seeds'].copy())
    return result


def paired_variants(first,second,first_strata,second_strata,transform):
    """Independent variant/stratum ensembles, never arbitrary cross-run pairing."""
    a = split_ensembles(first,first_strata)
    b = split_ensembles(second,second_strata)
    count = len(a)
    return independent_jackknife_many(a+b,lambda *means: transform(sum(means[:count]),sum(means[count:])))


def comparison_values(strict,product):
    return np.asarray([strict,product,product-strict,ratio(product,strict)])


def run(strict_path,product_path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a joint comparison')
    strict,product = read(strict_path),read(product_path)
    if strict['report']['variant'] != 'S' or product['report']['variant'] != 'Pi':
        raise ValueError('Provide S first and Pi second')
    assert_same_layout([strict,product])
    assert_matched_physics([strict,product])
    if set(map(tuple,strict['seeds'])) & set(map(tuple,product['seeds'])):
        raise ValueError('Random stream collision between nominally independent variants')
    report = dict(created_utc=now(),status='conditional joint-batch S/Pi pilot; full-run convergence required',
                  inputs={str(p.resolve()):digest(p) for p in (strict_path,product_path)},
                  convention='Ratios use the sums of stratum means, not means of per-worker ratios. '
                             'Delete from only one independent variant/stratum ensemble at a time. '
                             'Envelopes/PDF uncertainties are reduced after each derived observable.',
                  rates={},derived_rates={},spectra={})
    shape = strict['values'].shape
    arrays = {key:np.full((4,)+shape,np.nan) for key in ('absolute','absolute_mc_errors','normalized','normalized_mc_errors')}
    labels = ('S','Pi','Pi_minus_S','Pi_over_S')
    derived_labels = ('acceptance_1b','acceptance_2b','two_b_given_one_b',
                      'veto_0extra','fraction_1extra','fraction_2extra','fraction_3plus_extra')
    def estimate(a,b,transform):
        return paired_variants(a,b,strict['strata'],product['strata'],transform)
    def summarize(row,dimensionful):
        result = {}
        for i,label in enumerate(labels):
            result[label] = describe(row['value'][i])
            key = 'central_joint_mc_error_pb' if dimensionful and i<3 else 'central_joint_mc_error'
            result[label][key] = row['mc_error'][i,0]
            if i == 2:
                result[label]['signal_over_joint_mc_error'] = ratio(row['value'][i,0],row['mc_error'][i,0])
        return result
    for h,title in enumerate(strict['titles']):
        if ' rates:' not in title:
            continue
        first = int(strict['offsets'][h])
        if not np.any(product['values'][first:first+9,0]):
            continue
        a = strict['contributions'][:,first:first+9,:]
        b = product['contributions'][:,first:first+9,:]
        rates = estimate(a,b,comparison_values)
        arrays['absolute'][:,first:first+9] = rates['value']
        arrays['absolute_mc_errors'][:,first:first+9] = rates['mc_error']
        report['rates'][title] = {label:summarize({key:row[:,i] for key,row in rates.items()
                                                  if key in ('value','mc_error')},True)
                                 for i,label in enumerate(RATES)}
        derived = estimate(a,b,lambda s,p:comparison_values(derived_rates(s),derived_rates(p)))
        report['derived_rates'][title] = {label:summarize({key:row[:,i] for key,row in derived.items()
                                                          if key in ('value','mc_error')},False)
                                         for i,label in enumerate(derived_labels)}
        for local_id in range(2,22):
            index = h+local_id-1
            start,stop = map(int,strict['offsets'][index:index+2])
            inputs = []
            for result in (strict,product):
                denominator = result['contributions'][:,[first+i for i in normalization_bins(local_id)],:].sum(axis=1)
                inputs.append(np.concatenate([result['contributions'][:,start:stop,:],denominator[:,None,:]],axis=1))
            absolute = estimate(*inputs,lambda s,p:comparison_values(s[:-1],p[:-1]))
            normalized = estimate(*inputs,lambda s,p:comparison_values(ratio(s[:-1],s[-1]),ratio(p[:-1],p[-1])))
            arrays['absolute'][:,start:stop] = absolute['value']
            arrays['absolute_mc_errors'][:,start:stop] = absolute['mc_error']
            arrays['normalized'][:,start:stop] = normalized['value']
            arrays['normalized_mc_errors'][:,start:stop] = normalized['mc_error']
            selected = (normalized['value'][0,:,0]>0.) & (normalized['value'][1,:,0]>0.)
            errors = normalized['mc_error'][3,selected,0]
            report['spectra'][strict['titles'][index]] = dict(
                positive_in_both_bins=int(selected.sum()),
                undefined_ratio_mc_errors=int(np.sum(~np.isfinite(normalized['mc_error'][3,:,0]))),
                normalized_ratio_mc_error_quantiles=(np.quantile(errors,[.1,.5,.9]) if len(errors) else None),
                normalization_rate_bins=[RATES[i] for i in normalization_bins(local_id)])
        print('Compared joint vectors for '+title,flush=True)
    np.savez_compressed(output.with_suffix('.npz'),**arrays,labels=labels,
                        edges=strict['edges'],offsets=strict['offsets'],titles=strict['titles'])
    report['arrays_sha256'] = digest(output.with_suffix('.npz'))
    report['script_sha256'] = digest(Path(__file__))
    save(output,clean(report))
    print(output,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strict',required=True,type=Path)
    parser.add_argument('--product',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    run(args.strict,args.product,args.output)
