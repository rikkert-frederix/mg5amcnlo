"""Full-flavour W/mass contrasts from independent complete retrainings.

A physical state has one S/Pi ensemble per charge and ordered flavour.
Shared states enter covariance once, including when a state is a target
in one comparison and the reference in another.
"""
import itertools

import numpy as np

from full_flavour_statistics import FLAVOURS
from parameter_charge_statistics import (CHARGE_LABELS, CHARGE_SHAPE_LABELS,
    CHARGE_DERIVED_LABELS, LABELS as PARAMETER_CHARGE_LABELS, charge_prediction,
    contrast as charge_contrast)
from parameter_statistics import LABELS as PARAMETER_LABELS, contrast
from replica_statistics import independent_jackknife_many
from run_robustness_campaign import STATES, comparisons as campaign_comparisons

LABELS = tuple(label.replace('scan', 'target') for label in PARAMETER_LABELS)
CHARGE_CONTRAST_LABELS = tuple(label.replace('scan', 'target') for label in PARAMETER_CHARGE_LABELS)


def state_key(state):
    key = state['w_treatment'], state['decay_bottom_mass']
    if key not in STATES:
        raise ValueError('Unsupported W-treatment/decay-mass state')
    return key


def group_key(mode, mass, charge, variant, flavours):
    return mode, float(mass), charge, variant, tuple(flavours)


def comparison_keys(comparison, full_flavour=True):
    flavours = list(map(tuple, comparison['flavours']))
    allowed = {row['name']: row for row in campaign_comparisons() if row['charge'] == 'plus'}
    name, charge = comparison['name'], comparison['charge']
    if (name not in allowed or charge not in ('plus', 'minus') or not flavours
            or len(flavours) != len(set(flavours)) or not set(flavours) <= set(FLAVOURS)):
        raise ValueError('Invalid robustness comparison or ordered flavour selection')
    if full_flavour and set(flavours) != set(FLAVOURS):
        raise ValueError('A full direct e/mu prediction requires all eight ordered flavours')
    for origin in ('reference', 'target'):
        if state_key(comparison[origin]) != state_key(allowed[name][origin]):
            raise ValueError('Comparison name does not match the specified physical states')
    return {(origin, variant): [group_key(*state_key(comparison[origin]), charge, variant, flavour)
                for flavour in flavours]
            for origin, variant in itertools.product(('reference', 'target'), ('S', 'Pi'))}


def select(ensembles, comparisons, full_flavour=True):
    required = {key for row in comparisons for keys in comparison_keys(row, full_flavour).values()
                for key in keys}
    if not required <= set(ensembles):
        raise ValueError('Missing selected robustness ensemble')
    return {key: ensembles[key] for key in sorted(required)}


def inputs(ensembles, mappings):
    required = {key for mapping in mappings for keys in mapping.values() for key in keys}
    if not required or set(ensembles) != required:
        raise ValueError('Missing or unexpected robustness retraining ensemble')
    order = sorted(required)
    arrays = [np.asarray(ensembles[key], dtype=float) for key in order]
    if (any(row.ndim < 2 or len(row) < 5 or not np.isfinite(row).all() for row in arrays)
            or len({row.shape[1:] for row in arrays}) != 1):
        raise ValueError('Require at least five finite matched full-run vectors per ensemble')
    indices = {key: i for i, key in enumerate(order)}
    selections = [{label: [indices[key] for key in keys] for label, keys in mapping.items()}
                  for mapping in mappings]
    return order, arrays, selections


def estimate(ensembles, comparisons, transform=None, full_flavour=True):
    """Return [comparison, contrast, ...features] with shared-state covariance.

    Sum flavour means before normalization or ratios. Keep matched scale
    and PDF coordinates as the last feature axis until after this step.
    """
    names = [(row['name'], row['charge']) for row in comparisons]
    if not names or len(names) != len(set(names)):
        raise ValueError('Empty or duplicate robustness comparison')
    mappings = [comparison_keys(row, full_flavour) for row in comparisons]
    order, arrays, selections = inputs(ensembles, mappings)
    def combined(*means):
        output = []
        for selection in selections:
            sums = {key: sum(means[i] for i in chosen) for key, chosen in selection.items()}
            output.append(contrast(sums['reference', 'S'], sums['reference', 'Pi'],
                                   sums['target', 'S'], sums['target', 'Pi'], transform))
        return np.asarray(output)
    result = independent_jackknife_many(arrays, combined)
    result['ensemble_order'] = order
    return result


def paired_comparisons(pair):
    return [dict(pair, charge=charge) for charge in ('plus', 'minus')]


def pair_comparisons(comparisons, full_flavour=True):
    grouped = {}
    for row in comparisons:
        comparison_keys(row, full_flavour)
        charges = grouped.setdefault(row['name'], {})
        if row['charge'] in charges:
            raise ValueError('Duplicate robustness charge comparison')
        charges[row['charge']] = row
    pairs, unpaired = [], []
    for name, charges in grouped.items():
        example = next(iter(charges.values()))
        pair = {key: example[key] for key in ('name', 'target', 'reference')}
        if set(charges) != {'plus', 'minus'}:
            unpaired.append(dict(pair, available_charges=sorted(charges)))
            continue
        plus, minus = [set(map(tuple, charges[charge]['flavours'])) for charge in ('plus', 'minus')]
        if plus != minus:
            raise ValueError('Charge comparison requires identical ordered flavour selections')
        pairs.append(dict(pair, flavours=[list(flavour) for flavour in sorted(plus)]))
    return pairs, unpaired


def estimate_charges(ensembles, pairs, transform=None, full_flavour=True):
    """Return [pair, charge quantity, contrast, ...features].

    Charge sums use raw cross sections before shape or acceptance
    normalization. Absolute shifts remain defined at zero asymmetry.
    """
    if not pairs or len({row['name'] for row in pairs}) != len(pairs):
        raise ValueError('Empty or duplicate paired robustness comparison')
    mappings = [{(row['charge'], origin, variant): keys
                 for row in paired_comparisons(pair)
                 for (origin, variant), keys in comparison_keys(row, full_flavour).items()}
                for pair in pairs]
    order, arrays, selections = inputs(ensembles, mappings)
    def combined(*means):
        output = []
        for selection in selections:
            sums = {key: sum(means[i] for i in chosen) for key, chosen in selection.items()}
            predictions = [charge_prediction(sums['plus', origin, variant],
                                            sums['minus', origin, variant], transform)
                           for origin, variant in (('reference', 'S'), ('reference', 'Pi'),
                                                   ('target', 'S'), ('target', 'Pi'))]
            output.append(charge_contrast(*predictions))
        return np.asarray(output)
    result = independent_jackknife_many(arrays, combined)
    result['ensemble_order'] = order
    return result
