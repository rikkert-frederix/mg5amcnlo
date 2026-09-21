"""Exact contiguous rebinning of bin-integrated joint weight vectors."""
import numpy as np


def spans(original, requested):
    original,requested = np.asarray(original,dtype=float),np.asarray(requested,dtype=float)
    if original.ndim != 2 or original.shape[1] != 2 or len(original) == 0:
        raise ValueError('Require lower/upper boundaries of every original bin')
    if not np.isfinite(original).all() or not np.isfinite(requested).all():
        raise ValueError('Nonfinite bin boundary')
    if not np.array_equal(original[1:,0],original[:-1,1]) or not (original[:,1] > original[:,0]).all():
        raise ValueError('Original bins must be ordered and contiguous')
    if requested.ndim != 1 or len(requested) < 2 or not (np.diff(requested) > 0.).all():
        raise ValueError('Requested boundaries must be strictly increasing')
    existing = np.r_[original[:,0],original[-1,1]]
    if requested[0] != existing[0] or requested[-1] != existing[-1]:
        raise ValueError('Rebinning must retain the full original range, including empty bins')
    positions = []
    for edge in requested:
        matches = np.flatnonzero(existing == edge)
        if len(matches) != 1:
            raise ValueError('Requested boundary splits an original bin')
        positions.append(int(matches[0]))
    return list(zip(positions[:-1],positions[1:]))


def rebin(values,original,requested,axis=-2):
    """Sum raw joint vectors BEFORE nonlinear transformations or MC reduction.

    Never quadrature-add per-bin errors here: signed counterevents can correlate
    different bins. Convert the final merged integrals to densities separately.
    """
    values = np.asarray(values)
    if values.shape[axis] != len(original):
        raise ValueError('Bin axis does not match boundaries')
    moved = np.moveaxis(values,axis,0)
    result = np.asarray([moved[first:last].sum(axis=0) for first,last in spans(original,requested)])
    return np.moveaxis(result,0,axis)
