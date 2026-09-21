#!/usr/bin/env python3
"""Pilot covariance for rates, acceptances and rate-normalized spectra.

Uses validated equal-count split batches conditional on trained grids.
Continuous spectra are bin-integrated, with no folded overflow; the correct
rates histogram supplies each normalization. This never certifies sparse
tails or replaces subsequent independent-run convergence checks.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from campaign import digest, now, save
from pilot_report import RATES, clean, describe
from read_splits import split_ensembles
from replica_statistics import independent_jackknife_many, ratio


def normalization_bins(local_id):
    if local_id in range(2,11) or local_id in (19,20,21):
        return [3]
    if local_id in (11,12,13,14):
        return [4]
    if local_id in (15,18):
        return [6,7,8]
    if local_id == 16:
        return [7,8]
    if local_id == 17:
        return [5]
    raise ValueError('No spectrum normalization for this local histogram ID')


def derived_rates(total):
    return np.asarray([ratio(total[3],total[0]),ratio(total[4],total[0]),
                       ratio(total[4],total[3]),*[ratio(total[i],total[4]) for i in range(5,9)]])


def estimate(contributions,strata,transform=lambda total: total):
    ensembles = split_ensembles(contributions,strata)
    return independent_jackknife_many(ensembles,lambda *means: transform(sum(means)))


def run(path,output):
    if output.exists() or output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a joint pilot report')
    metadata = json.loads(path.with_suffix('.json').read_text())
    if digest(path) != metadata['arrays_sha256'] or not metadata['status'].startswith('equal-count conditional'):
        raise ValueError('Require an unmodified successful batch audit')
    with np.load(path) as data:
        values = data['contributions']
        strata = data['strata']
        old_errors = data['individual_errors']
        titles,offsets,edges = list(data['titles']),data['offsets'],data['edges']
    nominal = values.sum(axis=0)
    diagonal_hwu = np.sqrt(np.sum(old_errors**2,axis=0))
    report = dict(created_utc=now(),input=str(path.resolve()),input_sha256=digest(path),
                  status='conditional split-batch pilot; full-run convergence not certified',
                  convention='Delete from one subprocess/channel stratum at a time; '
                             'all bins and scale/PDF weights within that batch remain paired. '
                             'Independent stratum variances add. Envelopes/PDF reduction follow ratios. '
                             'Spectra are bin integrals normalized by the matching rate, not by a truncated sum.',
                  rates={},derived_rates={},spectra={},migrations={})
    names = ('acceptance_1b','acceptance_2b','two_b_given_one_b',
             'veto_0extra','fraction_1extra','fraction_2extra','fraction_3plus_extra')
    absolute_error = np.full(nominal.shape,np.nan)
    normalized = np.full(nominal.shape,np.nan)
    normalized_error = np.full(nominal.shape,np.nan)
    rate_starts = []
    for h,title in enumerate(titles):
        if ' rates:' not in title:
            continue
        first = int(offsets[h])
        if not np.any(nominal[first:first+9,0]):
            continue
        rate_starts.append(first)
        rates = estimate(values[:,first:first+9,:],strata)
        absolute_error[first:first+9] = rates['mc_error']
        derived = estimate(values[:,first:first+9,:],strata,derived_rates)
        report['rates'][title] = {}
        for j,label in enumerate(RATES):
            row = describe(rates['value'][j])
            row.update(central_joint_mc_error_pb=rates['mc_error'][j,0],
                       central_hwu_mc_error_pb=diagonal_hwu[first+j],
                       joint_to_hwu_error_ratio=ratio(rates['mc_error'][j,0],diagonal_hwu[first+j]))
            report['rates'][title][label] = row
        report['derived_rates'][title] = {}
        for j,label in enumerate(names):
            row = describe(derived['value'][j])
            row.update(central_joint_mc_error=derived['mc_error'][j,0],
                       central_jackknife_bias_estimate=derived['nonlinear_bias_estimate'][j,0])
            report['derived_rates'][title][label] = row
        for local_id in range(2,22):
            index = h+local_id-1
            start,stop = map(int,offsets[index:index+2])
            denominator = values[:,[first+i for i in normalization_bins(local_id)],:].sum(axis=1)
            combined = np.concatenate([values[:,start:stop,:],denominator[:,None,:]],axis=1)
            absolute = estimate(combined,strata)
            shape = estimate(combined,strata,lambda total: ratio(total[:-1],total[-1]))
            coverage = estimate(combined,strata,lambda total: ratio(total[:-1].sum(axis=0),total[-1]))
            absolute_error[start:stop] = absolute['mc_error'][:-1]
            normalized[start:stop] = shape['value']
            normalized_error[start:stop] = shape['mc_error']
            selected = (shape['value'][:,0] > 0.) & (shape['mc_error'][:,0] > 0.)
            relative = shape['mc_error'][selected,0]/shape['value'][selected,0]
            report['spectra'][titles[index]] = dict(
                normalization_rate_bins=[RATES[i] for i in normalization_bins(local_id)],
                central_coverage=coverage['value'][0],central_coverage_mc_error=coverage['mc_error'][0],
                positive_populated_bins=int(selected.sum()),
                normalized_bins_below_2percent=int(np.sum(relative < .02)),
                normalized_relative_error_quantiles=(np.quantile(relative,[.1,.5,.9]) if relative.size else None))
        print('Reduced covariance for '+title,flush=True)
    if len(rate_starts) != 5:
        raise ValueError('Expected all five cut configurations for one active charge')
    baseline = rate_starts[0]
    for start in rate_starts[1:]:
        title = next(titles[h] for h in range(len(titles)) if offsets[h] == start)
        contrast = values[:,[start+3,start+4,baseline+3,baseline+4],:]
        migration = estimate(contrast,strata,lambda total: ratio(total[:2],total[2:]))
        report['migrations'][title] = {
            label:dict(distribution=describe(migration['value'][j]),
                       central_joint_mc_error=migration['mc_error'][j,0])
            for j,label in enumerate(('one_b_ratio_to_R04_b25','two_b_ratio_to_R04_b25'))}
    output.parent.mkdir(parents=True,exist_ok=True)
    arrays = output.with_suffix('.npz')
    np.savez_compressed(arrays,values=nominal,absolute_joint_errors=absolute_error,
                        normalized=normalized,normalized_joint_errors=normalized_error,
                        edges=edges,offsets=offsets,titles=titles)
    report['arrays_sha256'] = digest(arrays)
    save(output,clean(report))
    print(output,flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batches',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    run(args.batches,args.output)
