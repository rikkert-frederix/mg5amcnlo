#!/usr/bin/env python3
"""Compare strict/product HwU files with matched 81- or shared 27-point grids.

Input files within one prediction must be disjoint charge/flavour samples,
not repeated estimates of the same integral. Output contents are bin integrals.
"""
import argparse
import itertools
import json
import math
from pathlib import Path
import re
import sys

FACTORS = (0.5, 1., 2.)
SHARED_POINTS = set(itertools.product(FACTORS, repeat=3))
POINTS = set(itertools.product(FACTORS, repeat=4))
SEVEN = {(r, f) for r, f in itertools.product(FACTORS, repeat=2) if .5 <= r/f <= 2.}
SHARED_GROUPS = {
    'production7': {(r, f, 1.) for r, f in SEVEN},
    'decay3': {(1., 1., d) for d in FACTORS},
    'combined21': {(r, f, d) for r, f in SEVEN for d in FACTORS},
    'common3': {(d, d, d) for d in FACTORS},
    'all27_diagnostic': SHARED_POINTS,
}
# Coordinates are always (production muR, production muF, top, antitop),
# regardless of the order of d6/d-6 in the runtime's signed-PDG labels.
GROUPS = {
    'production7': {(r, f, 1., 1.) for r, f in SEVEN},
    'top3': {(1., 1., d, 1.) for d in FACTORS},
    'antitop3': {(1., 1., 1., d) for d in FACTORS},
    'decay3': {(1., 1., d, d) for d in FACTORS},
    'decay9': {(1., 1., t, a) for t, a in itertools.product(FACTORS, repeat=2)},
    'combined21': {(r, f, d, d) for r, f in SEVEN for d in FACTORS},
    'combined63': {(r, f, t, a) for r, f in SEVEN
                   for t, a in itertools.product(FACTORS, repeat=2)},
    'common3': {(d, d, d, d) for d in FACTORS},
    'all81_diagnostic': POINTS,
}
NUMBER = r'([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)'
SCALE = re.compile(r'muR\s*=\s*'+NUMBER+r'\s+muF\s*=\s*'+NUMBER, re.I)
DECAY_SCALE = re.compile(r'(?:^|\s)d(-?\d+)\s*=\s*'+NUMBER, re.I)


def scale_point(label):
    match = SCALE.search(label)
    if match:
        def number(x):
            return float(x.replace('D', 'E').replace('d', 'e'))
        axes = DECAY_SCALE.findall(label)
        if not axes:
            return None
        decay = {int(pdg): number(value) for pdg, value in axes}
        if len(decay) != len(axes) or set(decay) not in ({6}, {6, -6}):
            raise ValueError('Expected unique d6 or d6/d-6 decay scale axes: ' + label)
        return tuple(map(number, match.groups())) + tuple(decay[pdg] for pdg in
                     ([6, -6] if -6 in decay else [6]))
    return None


def grid(values):
    points = set(values)-{'central'}
    if points == POINTS:
        return POINTS, GROUPS
    if points == SHARED_POINTS:
        return SHARED_POINTS, SHARED_GROUPS
    raise ValueError('Expected all 27 shared or 81 signed-PDG production/decay scale points')


def ratio(a, b):
    """Undefined/negative-denominator ratios remain explicitly absent."""
    return a/b if a is not None and b is not None and b > 0. else None


def describe(values):
    """Envelope only after computing the derived observable at each point."""
    points, groups = grid(values)
    result = {'central': values['central'], 'points': {}}
    for point in sorted(points):
        result['points'][','.join('%g' % x for x in point)] = values[point]
    for name, group in groups.items():
        selected = [values[point] for point in group]
        valid = [v for v in selected if v is not None and math.isfinite(v)]
        result[name] = dict(npoints=len(selected), nvalid=len(valid),
                            envelope=[min(valid), max(valid)] if len(valid)==len(selected) else None)
    result['log_responses'] = {}
    responses = [('production',(.5,.5,1.),(2.,2.,1.)),
                 ('decay',(1.,1.,.5),(1.,1.,2.))]
    if points == POINTS:
        responses = [('production',(.5,.5,1.,1.),(2.,2.,1.,1.)),
                     ('decay',(1.,1.,.5,.5),(1.,1.,2.,2.)),
                     ('top',(1.,1.,.5,1.),(1.,1.,2.,1.)),
                     ('antitop',(1.,1.,1.,.5),(1.,1.,1.,2.))]
    for name, low, high in responses:
        a,b = values[low],values[high]
        slope = (b-a)/math.log(4.) if a is not None and b is not None else None
        result['log_responses'][name] = dict(absolute=slope, relative=ratio(slope, values['central']))
    return result


def compare_values(strict, product):
    points, _ = grid(strict)
    if grid(product)[0] != points:
        raise ValueError('Strict/product scale grids differ; use matching decay-scale grouping')
    keys = ['central'] + sorted(points)
    delta = {k: product[k]-strict[k] if product[k] is not None and strict[k] is not None else None
             for k in keys}
    return dict(strict=describe(strict), product=describe(product), difference=describe(delta),
                product_over_strict=describe({k: ratio(product[k], strict[k]) for k in keys}))


def load_sum(paths):
    for parent in Path(__file__).resolve().parents:
        if (parent / 'madgraph/various/histograms.py').is_file():
            sys.path.insert(0, str(parent))
            break
    from madgraph.various.histograms import HwUList
    result = {}
    expected_points = None
    resolved = [str(Path(p).resolve()) for p in paths]
    if len(resolved)!=len(set(resolved)):
        raise ValueError('An input file was provided twice')
    for filename in resolved:
        current = {}
        for hist in HwUList(filename, raw_labels=True):
            if not re.match(r'^R0[345]_b(?:25|30|40) W[+-] ', hist.title):
                continue
            if hist.title in current:
                raise ValueError('Duplicate histogram: ' + hist.title)
            rows = []
            for b in hist.bins:
                weights = dict(b.wgts)
                values = {'central': weights['central value']}
                for label, value in weights.items():
                    point = scale_point(label)
                    if point is not None:
                        if point in values:
                            raise ValueError('Duplicate scale point (use one central-scale choice per file)')
                        values[point] = value
                points, _ = grid(values)
                if expected_points is None:
                    expected_points = points
                if points != expected_points:
                    raise ValueError('Scale grids differ between input samples or bins')
                if not all(math.isfinite(v) for v in values.values()):
                    raise ValueError('Non-finite weight in ' + hist.title)
                rows.append(dict(edges=list(b.boundaries), values=values,
                                 variance=weights['dy']**2))
            current[hist.title] = rows
        if len(current)!=210:
            raise ValueError('Expected 210 product-study histograms in ' + filename)
        if not result:
            result = current
            continue
        if current.keys()!=result.keys():
            raise ValueError('Histogram titles differ between disjoint input samples')
        for title, rows in current.items():
            if len(rows)!=len(result[title]):
                raise ValueError('Different bin counts: ' + title)
            for target, row in zip(result[title], rows):
                if target['edges']!=row['edges']:
                    raise ValueError('Different bin edges: ' + title)
                for key, value in row['values'].items():
                    target['values'][key] += value
                target['variance'] += row['variance']
    return result


def make_report(strict, product):
    if strict.keys()!=product.keys():
        raise ValueError('Strict/product histogram sets differ')
    rates = {title[:10]: (rows, product[title]) for title, rows in strict.items() if ' rates:' in title}
    output = dict(
        conventions='bin integrals; correlated scale points; 2b conditional spectra normalized to inclusive 2b rate',
        uncertainty_note='MC errors for differences assume independent strict/product integrations. '
        'No MC errors for ratios or normalized shapes without bin/rate covariance. '
        'Common seeds alone do not establish covariance. Zero/nonpositive denominators are null.',
        histograms={}, derived={})
    first_values = next(iter(strict.values()))[0]['values']
    points, _ = grid(first_values)
    output['scale_axes'] = ['muR_production', 'muF_production'] + (
        ['muR_top', 'muR_antitop'] if points == POINTS else ['muR_decay_shared'])
    output['scale_point_count'] = len(points)
    for title, srows in strict.items():
        prows = product[title]
        if len(srows)!=len(prows):
            raise ValueError('Strict/product bin counts differ')
        records = []
        for s, p in zip(srows, prows):
            if s['edges']!=p['edges']:
                raise ValueError('Strict/product bin edges differ')
            record = dict(edges=s['edges'], absolute=compare_values(s['values'], p['values']),
                          mc_error_strict=math.sqrt(s['variance']),
                          mc_error_product=math.sqrt(p['variance']),
                          mc_error_difference_independent=math.sqrt(s['variance']+p['variance']))
            tier = 3 if ' 1b ' in title else 4 if ' 2b ' in title else None
            if tier is not None:
                sr, pr = rates[title[:10]]
                sn = {k: ratio(v, sr[tier]['values'][k]) for k, v in s['values'].items()}
                pn = {k: ratio(v, pr[tier]['values'][k]) for k, v in p['values'].items()}
                record['normalized_to_fiducial'] = compare_values(sn, pn)
            records.append(record)
        output['histograms'][title] = records
    for prefix, (sr, pr) in rates.items():
        derived = {}
        for name, num, den in [('acceptance_1b',3,0), ('acceptance_2b',4,0),
                               ('two_b_given_one_b',4,3), ('zero_extra_fraction_2b',5,4)]:
            def fractions(rows):
                return {k: ratio(v,rows[den]['values'][k]) for k,v in rows[num]['values'].items()}
            derived[name] = compare_values(fractions(sr), fractions(pr))
        output['derived'][prefix] = derived
        if prefix.endswith('W+'):
            sm, pm = rates[prefix[:-1]+'-']
            for name, stage in [('input',0), ('fiducial_1b',3), ('fiducial_2b',4)]:
                def charges(plus, minus, asymmetry=False):
                    return {k: ratio(v-minus[stage]['values'][k], v+minus[stage]['values'][k])
                            if asymmetry else ratio(v,minus[stage]['values'][k])
                            for k,v in plus[stage]['values'].items()}
                derived['charge_ratio_'+name] = compare_values(charges(sr,sm), charges(pr,pm))
                derived['charge_asymmetry_'+name] = compare_values(charges(sr,sm,True), charges(pr,pm,True))
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--strict', nargs='+', required=True)
    parser.add_argument('--product', nargs='+', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        report = make_report(load_sum(args.strict), load_sum(args.product))
        report['inputs'] = dict(strict=args.strict, product=args.product)
        # Exclusive creation prevents silently replacing a comparison.
        with open(args.output, 'x') as stream:
            json.dump(report, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print('Wrote '+args.output)
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, str(error)+'\n')


if __name__ == '__main__':
    main()
