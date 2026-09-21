"""Point-matched parameter contrasts with independent retraining ensembles.

Shared baseline groups are passed to the estimator once, preserving their
covariance across different parameter scenarios. Neither scale coordinates
nor equal replica indices in separate ensembles create artificial pairing.
"""
import itertools

import numpy as np

from full_flavour_statistics import FLAVOURS
from replica_statistics import independent_jackknife_many, ratio

LABELS = ('reference_S', 'reference_Pi', 'scan_S', 'scan_Pi',
          'reference_Pi_minus_S', 'scan_Pi_minus_S', 'reference_Pi_over_S', 'scan_Pi_over_S',
          'S_scan_over_reference', 'Pi_scan_over_reference',
          'product_shift_change', 'product_ratio_change', 'product_ratio_double_ratio')


def group_key(origin, scenario, mode, charge, variant, flavours):
    return origin, scenario, mode, charge, variant, tuple(flavours)


def comparison_keys(comparison, full_flavour=True):
    scenario, mode, charge = (comparison[key] for key in ('scenario', 'w_treatment', 'charge'))
    flavours = [tuple(row) for row in comparison['flavours']]
    if (not flavours or len(set(flavours)) != len(flavours) or not set(flavours) <= set(FLAVOURS)
            or charge not in ('plus', 'minus') or mode not in ('onshell', 'top-bw', 'all-bw')):
        raise ValueError('Invalid selected parameter-comparison context')
    if full_flavour and set(flavours) != set(FLAVOURS):
        raise ValueError('A full direct e/mu prediction requires all eight ordered flavours')
    return {(origin, variant): [group_key(origin, 'baseline' if origin == 'reference' else scenario,
                mode, charge, variant, flavour) for flavour in flavours]
            for origin, variant in itertools.product(('reference', 'scan'), ('S', 'Pi'))}


def contrast(reference_s, reference_pi, scan_s, scan_pi, transform=None):
    if transform is not None:
        reference_s, reference_pi, scan_s, scan_pi = [transform(row)
            for row in (reference_s, reference_pi, scan_s, scan_pi)]
    delta0, delta1 = reference_pi-reference_s, scan_pi-scan_s
    r0, r1 = ratio(reference_pi, reference_s), ratio(scan_pi, scan_s)
    return np.asarray([reference_s, reference_pi, scan_s, scan_pi, delta0, delta1, r0, r1,
                       ratio(scan_s, reference_s), ratio(scan_pi, reference_pi),
                       delta1-delta0, r1-r0, ratio(r1, r0)])


def estimate(ensembles, comparisons, transform=None, full_flavour=True):
    """Return [comparison, contrast, ...features] values/errors/covariance.

    Ratios and normalization act on flavour sums of ensemble means. For large
    spectra, reduce each scenario's scale weights separately and evaluate a
    joint nominal slice to retain cross-scenario shared-baseline covariance.
    """
    if not comparisons:
        raise ValueError('Empty parameter comparison')
    labels = [(row['scenario'], row['w_treatment'], row['charge']) for row in comparisons]
    if len(labels) != len(set(labels)):
        raise ValueError('Duplicate parameter comparison')
    mappings = [comparison_keys(row, full_flavour) for row in comparisons]
    required = {key for mapping in mappings for keys in mapping.values() for key in keys}
    if set(ensembles) != required:
        raise ValueError('Missing or unexpected reference/scan retraining ensemble')
    order = sorted(required)
    arrays = [np.asarray(ensembles[key], dtype=float) for key in order]
    if (any(row.ndim < 2 or len(row) < 5 or not np.isfinite(row).all() for row in arrays)
            or len({row.shape[1:] for row in arrays}) != 1):
        raise ValueError('Require at least five finite matched full-run vectors per ensemble')
    indices = {key: i for i, key in enumerate(order)}
    selections = [{name: [indices[key] for key in keys] for name, keys in mapping.items()}
                  for mapping in mappings]
    def combined(*means):
        outputs = []
        for selection in selections:
            sums = {key: sum(means[i] for i in selected) for key, selected in selection.items()}
            outputs.append(contrast(sums['reference', 'S'], sums['reference', 'Pi'],
                                    sums['scan', 'S'], sums['scan', 'Pi'], transform))
        return np.asarray(outputs)
    result = independent_jackknife_many(arrays, combined)
    result['ensemble_order'] = order
    return result


def select(ensembles, comparisons, full_flavour=True):
    """Select exactly the shared groups needed by a comparison subset."""
    required = {key for comparison in comparisons for values in comparison_keys(comparison, full_flavour).values()
                for key in values}
    if not required <= set(ensembles):
        raise ValueError('Missing selected parameter-comparison ensemble')
    return {key: ensembles[key] for key in sorted(required)}
