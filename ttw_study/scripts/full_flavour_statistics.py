"""Independent retraining ensembles for complete direct e/mu predictions.

Average repeated estimates of each ordered assignment, then sum its eight
disjoint flavour means. Delete from one ensemble at a time: equal replica
indices in independent flavour/charge/prescription samples are not pairing.
"""
import itertools

import numpy as np

from replica_statistics import independent_jackknife_many

FLAVOURS=tuple(itertools.product(('e','mu'),repeat=3))
VARIANTS=('S','Pi','LO','P','D','PiD')


def inventory(ensembles,charges,variants):
    charges,variants=tuple(charges),tuple(variants)
    if (not charges or not variants or len(charges)!=len(set(charges)) or
            len(variants)!=len(set(variants)) or not set(charges)<= {'plus','minus'} or
            not set(variants)<=set(VARIANTS)):
        raise ValueError('Require distinct supported charges and prescriptions')
    order=list(itertools.product(charges,variants,FLAVOURS))
    if set(ensembles)!=set(order):
        raise ValueError('Require exactly eight ordered e/mu ensembles per selected charge/prescription')
    arrays=[np.asarray(ensembles[key],dtype=float) for key in order]
    if (any(a.ndim<2 or len(a)<5 or not np.isfinite(a).all() for a in arrays) or
            len({a.shape[1:] for a in arrays})!=1):
        raise ValueError('Require at least five finite, matched full-run vectors in every flavour ensemble')
    return order,arrays


def estimate(ensembles,charges=('plus','minus'),variants=VARIANTS,transform=None):
    """Transform a dictionary {(charge,variant): complete flavour sum}.

    Retain the final feature axis as the full common scale/PDF grid when
    reducing physical data. The callback forms ratios/shapes/charge algebra
    before any uncertainty-envelope or PDF-member reduction.
    """
    order,arrays=inventory(ensembles,charges,variants)
    predictions=list(dict.fromkeys(key[:2] for key in order))
    indices={key:[i for i,row in enumerate(order) if row[:2]==key] for key in predictions}
    def combined(*means):
        totals={key:sum(means[i] for i in selected) for key,selected in indices.items()}
        return (np.asarray([totals[key] for key in predictions]) if transform is None else transform(totals))
    return independent_jackknife_many(arrays,combined)
