"""Charge observables of matched parameter scans with independent run errors.

Charge sums are formed before any acceptance or shape normalization. Ratios
and asymmetries use matching scale coordinates. A reference used by several
scenarios is one shared ensemble, while different charges remain independent.
"""
import numpy as np

from parameter_statistics import comparison_keys
from replica_statistics import independent_jackknife_many, ratio


LABELS = ('reference_S', 'reference_Pi', 'scan_S', 'scan_Pi',
          'reference_Pi_minus_S', 'scan_Pi_minus_S',
          'S_scan_minus_reference', 'Pi_scan_minus_reference', 'product_shift_change')
CHARGE_LABELS = ('sum', 'plus_over_minus', 'asymmetry')
CHARGE_SHAPE_LABELS = ('combined', 'plus_shape_over_minus_shape', 'shape_asymmetry')
CHARGE_DERIVED_LABELS = ('combined', 'plus_over_minus', 'asymmetry')


def paired_comparisons(pair):
    return [dict(pair, charge=charge) for charge in ('plus', 'minus')]


def pair_comparisons(comparisons, full_flavour=True):
    """Find complete charge pairs and record missing counterparts explicitly."""
    grouped = {}
    for comparison in comparisons:
        comparison_keys(comparison, full_flavour)
        key = comparison['scenario'], comparison['w_treatment']
        charges = grouped.setdefault(key, {})
        if comparison['charge'] in charges:
            raise ValueError('Duplicate parameter charge comparison')
        charges[comparison['charge']] = comparison
    pairs, unpaired = [], []
    for (scenario, mode), charges in grouped.items():
        row = dict(scenario=scenario, w_treatment=mode)
        if set(charges) != {'plus', 'minus'}:
            unpaired.append(dict(row, available_charges=sorted(charges),
                reason='Charge comparison requires both independently integrated charges'))
            continue
        plus, minus = (set(map(tuple, charges[charge]['flavours'])) for charge in ('plus', 'minus'))
        if plus != minus:
            raise ValueError('Charge comparisons require the same ordered flavour selection')
        pairs.append(dict(row, flavours=[list(flavour) for flavour in sorted(plus)]))
    return pairs, unpaired


def charge_prediction(plus, minus, transform=None):
    """Return combined prediction, charge ratio and charge asymmetry."""
    combined = plus+minus
    if transform is not None:
        combined, plus, minus = (transform(row) for row in (combined, plus, minus))
    return np.asarray([combined, ratio(plus, minus), ratio(plus-minus, plus+minus)])


def contrast(reference_s, reference_pi, scan_s, scan_pi):
    """Absolute changes also apply when a charge asymmetry crosses zero."""
    delta0, delta1 = reference_pi-reference_s, scan_pi-scan_s
    return np.stack([reference_s, reference_pi, scan_s, scan_pi, delta0, delta1,
                     scan_s-reference_s, scan_pi-reference_pi, delta1-delta0], axis=1)


def estimate(ensembles, pairs, transform=None, full_flavour=True):
    """Return [pair, charge quantity, contrast, ...features] and covariance."""
    if not pairs:
        raise ValueError('Empty paired parameter comparison')
    identifiers = [(pair['scenario'], pair['w_treatment']) for pair in pairs]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError('Duplicate paired parameter comparison')
    mappings = [{comparison['charge']: comparison_keys(comparison, full_flavour)
                 for comparison in paired_comparisons(pair)} for pair in pairs]
    required = {key for mapping in mappings for charge in mapping.values()
                for keys in charge.values() for key in keys}
    if set(ensembles) != required:
        raise ValueError('Missing or unexpected paired reference/scan ensemble')
    order = sorted(required)
    arrays = [np.asarray(ensembles[key], dtype=float) for key in order]
    if (any(row.ndim < 2 or len(row) < 5 or not np.isfinite(row).all() for row in arrays)
            or len({row.shape[1:] for row in arrays}) != 1):
        raise ValueError('Require at least five finite matched full-run vectors per ensemble')
    indices = {key: i for i, key in enumerate(order)}
    selections = [{(charge, origin, variant): [indices[key] for key in keys]
                   for charge, mapping in mappings_for_pair.items()
                   for (origin, variant), keys in mapping.items()}
                  for mappings_for_pair in mappings]
    def combined(*means):
        output = []
        for selection in selections:
            sums = {key: sum(means[i] for i in selected) for key, selected in selection.items()}
            predictions = [charge_prediction(sums['plus', origin, variant],
                                              sums['minus', origin, variant], transform)
                           for origin, variant in (('reference', 'S'), ('reference', 'Pi'),
                                                   ('scan', 'S'), ('scan', 'Pi'))]
            output.append(contrast(*predictions))
        return np.asarray(output)
    result = independent_jackknife_many(arrays, combined)
    result['ensemble_order'] = order
    return result
