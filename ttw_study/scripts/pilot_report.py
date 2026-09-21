#!/usr/bin/env python3
"""Quantify pilot precision and strict/product algebra without certifying physics.

Inputs are one independently seeded audit per variant of a single ordered
flavour/charge sample. This is not an integration-replica combiner.
"""
import argparse
import importlib.util
from pathlib import Path

import numpy as np

from campaign import ROOT, digest, now, save
from load_results import POINTS, assert_matched_physics, assert_same_layout, load, strict_identity
from pdf_statistics import uncertainty
from replica_statistics import ratio

RATES = ('input','leptons','mass_requirements','fiducial_1b','fiducial_2b',
         '2b_0extra','2b_1extra','2b_2extra','2b_3plus_extra')
spec = importlib.util.spec_from_file_location(
    'provided_scale_comparison',ROOT/'Template/fNLO/FixedOrderAnalysis/ttw_product_scales.py')
scales = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scales)


def clean(value):
    if isinstance(value,np.ndarray):
        return clean(value.tolist())
    if isinstance(value,np.generic):
        return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):
        return None
    if isinstance(value,dict):
        return {key:clean(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):
        return [clean(item) for item in value]
    return value


def describe(values):
    if values.shape != (183,):
        raise ValueError('Expected canonical nominal, 81-scale and 101-PDF weights')
    scale_values = {'central':clean(values[0])}
    scale_values.update({point:clean(values[i+1]) for i,point in enumerate(POINTS)})
    return dict(scale=scales.describe(scale_values),pdf=uncertainty(values[82:]))


def independent_linear(results,index,coefficients):
    """A point-matched linear contrast, with independent-variant MC errors."""
    selected = [results[name] for name in coefficients]
    seeds = [row['report']['manifest']['settings']['iseed'] for row in selected]
    if len(set(seeds)) != len(seeds):
        raise ValueError('Independent contrast requires distinct variant seeds')
    value = sum(coefficient*results[name]['values'][index]
                for name,coefficient in coefficients.items())
    error = np.sqrt(sum((coefficient*results[name]['errors'][index])**2
                       for name,coefficient in coefficients.items()))
    result = describe(value)
    result.update(central_mc_error_pb=error,
                  central_signal_to_mc_error=ratio(value[0],error),
                  shift_resolved_at_3sigma=bool(error > 0. and abs(value[0]) > 3*error))
    return result


def make_report(paths):
    results = {}
    for path in paths:
        result = load(path)
        variant = result['report']['variant']
        if variant in results:
            raise ValueError('Multiple estimates of a variant require a replica combiner, not a sum')
        results[variant] = result
    assert_same_layout(list(results.values()))
    assert_matched_physics(list(results.values()))
    reference = next(iter(results.values()))
    report = dict(created_utc=now(),status='pilot diagnostic; no publication precision certification',
                  inputs={str(Path(path).resolve()):digest(path) for path in paths},
                  precision={},rates={},strict_identity=None,
                  convention='All ratios and differences precede scale envelopes and PDF reduction. '
                             'Difference MC errors assume independently seeded variants. '
                             'No ratio, acceptance or normalized-shape MC error is inferred from individual HwU bins.')
    for variant,result in results.items():
        central,error = result['values'][:,0],result['errors']
        selected = (central > 0.) & (error > 0.)
        relative = error[selected]/central[selected]
        report['precision'][variant] = dict(
            elapsed_s=result['report']['execution']['elapsed_s'],
            positive_populated_bins=int(selected.sum()),
            relative_mc_error_quantiles=(dict(zip(('q10','median','q90'),np.quantile(relative,[.1,.5,.9])))
                                         if relative.size else None),
            bins_below_2percent=int((relative<.02).sum()),
            bins_below_1percent=int((relative<.01).sum()),
            nonzero_bins_with_zero_error=int(((central != 0.) & (error == 0.)).sum()))
    if {'LO','P','D','S'} <= results.keys():
        identity = strict_identity({name:results[name] for name in ('LO','P','D','S')})
        finite = np.isfinite(identity['central_pull'])
        pulls = identity['central_pull'][finite]
        report['strict_identity'] = dict(
            note=identity['note'],finite_pull_bins=int(finite.sum()),
            max_absolute_pull=float(np.max(np.abs(pulls))) if pulls.size else None,
            bins_above_3sigma=int((np.abs(pulls)>3.).sum()),
            bins_above_5sigma=int((np.abs(pulls)>5.).sum()),
            nonzero_difference_with_zero_error=int(((identity['central_mc_error'] == 0.) &
                                                    (identity['difference'][:,0] != 0.)).sum()))
    for h,title in enumerate(reference['titles']):
        if ' rates:' not in title:
            continue
        start = reference['offsets'][h]
        if not any(np.any(result['values'][start:start+9,0]) for result in results.values()):
            continue
        rates = report['rates'][title] = {}
        for j,label in enumerate(RATES):
            i = start+j
            row = rates[label] = {}
            for variant,result in results.items():
                row[variant] = describe(result['values'][i])
                row[variant]['central_mc_error_pb'] = result['errors'][i]
                row[variant]['relative_mc_error'] = ratio(result['errors'][i],abs(result['values'][i,0]))
                row[variant]['cost_multiplier_to_0p5percent'] = (
                    row[variant]['relative_mc_error']/.005)**2
            if {'S','Pi'} <= results.keys():
                s,p = results['S'],results['Pi']
                row['Pi_minus_S'] = independent_linear(results,i,{'Pi':1.,'S':-1.})
                row['Pi_over_S'] = describe(ratio(p['values'][i],s['values'][i]))
            if {'S','P'} <= results.keys():
                row['S_over_P'] = describe(ratio(results['S']['values'][i],results['P']['values'][i]))
            if {'PiD','D'} <= results.keys():
                row['PiD_minus_D'] = independent_linear(results,i,{'PiD':1.,'D':-1.})
            if {'Pi','S','PiD','D'} <= results.keys():
                row['production_related_product_shift'] = independent_linear(
                    results,i,{'Pi':1.,'S':-1.,'PiD':-1.,'D':1.})
                row['production_related_product_shift']['definition'] = (
                    '(Pi-S)-(PiD-D); not a pure homogeneous coefficient in alpha_s')
            if report['strict_identity'] is not None:
                row['S_minus_P_minus_D_plus_LO'] = dict(
                    central_pb=identity['difference'][i,0],
                    central_mc_error_pb=identity['central_mc_error'][i],
                    central_pull=identity['central_pull'][i])
    return clean(report)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audits',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Refuse to overwrite an archived pilot report')
    save(args.output,make_report(args.audits))
    print(args.output)
