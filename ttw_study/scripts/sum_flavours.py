#!/usr/bin/env python3
"""Sum exactly eight ordered direct e/mu assignments, never rescale a subset.

This accepts one independent estimate per ordered assignment. Repeated
estimates of the same assignment need a replica combiner, not another term
in this sum. Derived-observable MC covariance requires the joint vectors.
"""
import argparse
import itertools
from pathlib import Path

import numpy as np

from campaign import digest, now, save
from load_results import assert_matched_physics, assert_same_layout, load

FLAVOURS = set(itertools.product(('e','mu'),repeat=3))


def combine(results):
    if len(results) != 8:
        raise ValueError('Require exactly all eight ordered direct e/mu assignments')
    flavours, seeds, comparable, variants = [], [], [], set()
    for result in results:
        report=result['report']
        generation=report['execution']['export_generation']['args']
        flavours.append(tuple(generation['flavours']))
        seeds.append(report['manifest']['settings']['iseed'])
        variants.add(report['variant'])
        # Remove ONLY the explicitly summed flavour labels from comparison.
        # The original audit object is never modified. Charge, widths, masses,
        # cuts, PDF, scale definitions and corrected stages must still match.
        generation=dict(generation,flavours=['e','e','e'])
        exported=dict(report['execution']['export_generation'],args=generation)
        execution=dict(report['execution'],export_generation=exported)
        comparable.append(dict(result,report=dict(report,execution=execution)))
    if set(flavours) != FLAVOURS:
        raise ValueError('Missing, duplicated or non-direct flavour assignment')
    if len(set(seeds)) != 8:
        raise ValueError('Independent-error sum requires distinct integration seeds')
    if len(variants) != 1:
        raise ValueError('Cannot sum different perturbative prescriptions')
    assert_same_layout(results)
    assert_matched_physics(comparable)
    reference=results[0]
    generation=reference['report']['execution']['export_generation']['args']
    return dict(values=sum(result['values'] for result in results),
                errors=np.sqrt(sum(result['errors']**2 for result in results)),
                edges=reference['edges'],offsets=reference['offsets'],titles=reference['titles'],
                weights=reference.get('weights'),variant=next(iter(variants)),
                charge=generation['charge'],w_treatment=generation['w_treatment'],
                decay_bottom_mass=generation['decay_bottom_mass'],
                flavours=sorted(flavours),seeds=seeds,
                status='complete direct-flavour algebraic sum; MC convergence not certified',
                mc_convention='Central-bin errors added in quadrature for independently seeded '
                              'flavour estimates. This alone supplies no acceptance/shape covariance.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('audits',nargs='+',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.npz').exists():
        raise ValueError('Refuse to overwrite a flavour-sum archive')
    result=combine([load(path) for path in args.audits])
    arrays={key:result.pop(key) for key in ('values','errors','edges','offsets','titles','weights')}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output.with_suffix('.npz'),**arrays)
    result.update(created_utc=now(),kind='disjoint_flavour_sum',
                  inputs={str(path.resolve()):digest(path) for path in args.audits},
                  arrays_sha256=digest(args.output.with_suffix('.npz')),
                  script_sha256=digest(Path(__file__)))
    save(args.output,result)
    print(args.output)


if __name__ == '__main__':
    main()
