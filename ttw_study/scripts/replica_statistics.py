"""Replica-level MC uncertainty with correlations retained within each run.

Inputs are independent, equally weighted integration replicas of the SAME
prediction. Sum disjoint flavours at matched scale/PDF points; for their MC
errors use independent_jackknife_many, not arbitrary cross-flavour run pairing.
Never treat PDF members, scale weights, or histogram bins as MC replicas.
"""
import numpy as np


def _replicas(values):
    values = np.asarray(values, dtype=float)
    if values.ndim < 2 or len(values) < 5:
        raise ValueError('At least five independent vector-valued integration replicas are required')
    if not np.isfinite(values).all():
        raise ValueError('Non-finite integration replica')
    return values


def covariance_factor(values):
    """Flatten F's trailing feature axes: F.T @ F is covariance of the mean."""
    values = _replicas(values)
    return (values-values.mean(axis=0))/np.sqrt(len(values)*(len(values)-1))


def ratio(numerator, denominator):
    """No ratio is defined for a zero/negative denominator."""
    numerator, denominator = np.broadcast_arrays(numerator, denominator)
    result = np.full(numerator.shape, np.nan)
    np.divide(numerator, denominator, out=result, where=denominator > 0.)
    return result


def jackknife(values, transform):
    """Evaluate a nonlinear observable of the mean, with leave-one-run errors.

    The reported central value is transform(mean), not mean(transform).
    Whole-run deletion preserves bin/normalization/scale/PDF correlations.
    A nonfinite leave-one-out result explicitly invalidates that uncertainty.
    """
    values = _replicas(values)
    count = len(values)
    total = values.sum(axis=0)
    central = np.asarray(transform(total/count), dtype=float)
    deleted = np.asarray([transform((total-row)/(count-1)) for row in values])
    valid = np.isfinite(central) & np.isfinite(deleted).all(axis=0)
    factors = np.sqrt((count-1)/count)*(deleted-deleted.mean(axis=0))
    error = np.sqrt(np.sum(factors*factors, axis=0))
    bias_estimate = (count-1)*(deleted.mean(axis=0)-central)
    return dict(value=central, mc_error=np.where(valid,error,np.nan),
                covariance_factor=factors, valid=valid,
                nonlinear_bias_estimate=np.where(valid,bias_estimate,np.nan))


def independent_jackknife(first, second, transform):
    """Nonlinear comparison of two explicitly independent replica ensembles.

    Delete from one ensemble at a time. Arbitrarily pairing independent runs
    would introduce a spurious estimated cross-ensemble covariance.
    """
    return independent_jackknife_many([first,second],transform)


def independent_jackknife_many(ensembles, transform):
    """Independent flavour/charge/variant ensembles, including four-run shifts.

    Delete a whole run from ONE ensemble at a time. This retains correlations
    among that run's bins and weights without manufacturing cross-ensemble
    correlations. The central transform receives all ensemble means as args.
    """
    ensembles = [_replicas(values) for values in ensembles]
    if not ensembles:
        raise ValueError('At least one integration ensemble is required')
    means = [values.mean(axis=0) for values in ensembles]
    partial = []
    for i,values in enumerate(ensembles):
        partial.append(jackknife(values,lambda mean: transform(*(means[:i]+[mean]+means[i+1:]))))
    valid = np.logical_and.reduce([row['valid'] for row in partial])
    variance = sum(row['mc_error']**2 for row in partial)
    return dict(value=np.asarray(transform(*means)),valid=valid,
                mc_error=np.where(valid,np.sqrt(variance),np.nan),
                covariance_factor=np.concatenate([row['covariance_factor'] for row in partial]),
                nonlinear_bias_estimate=sum(row['nonlinear_bias_estimate'] for row in partial))
