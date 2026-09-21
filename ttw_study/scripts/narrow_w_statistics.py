"""Joint narrow-W comparisons with a shared independently sampled reference.

This is a statistics kernel, not validation of generated physics inputs.
The caller must audit common scales, matched top widths, unchanged weak
couplings, epsilon*Gamma_W, full support and distinct actual RNG streams.
"""
import numpy as np
from read_splits import split_ensembles
from replica_statistics import independent_jackknife_many, ratio


def width_factors(factors, count):
    factors=np.asarray(factors,dtype=float)
    if (factors.shape!=(count,) or count<2 or not np.isfinite(factors).all()
            or np.any(factors<=0.) or np.any(factors>1.) or factors[0]!=1.):
        raise ValueError('Require one physical-width reference and positive narrowing factors')
    return factors


def comparisons(means,factors,normalized=False):
    """Return scaled predictions, differences and ratios to sample zero.

    For normalized spectra, each input ends with its correct parent rate.
    Epsilon cubed applies to numerator and denominator alike, never as an
    additional event weight or a second branching factor.
    """
    factors=width_factors(factors,len(means))
    if len({np.asarray(row).shape for row in means})!=1:
        raise ValueError('All samples must have the same bins and weights')
    scaled=[np.asarray(row)*factor**3 for row,factor in zip(means,factors)]
    if normalized:
        if any(row.ndim<2 or len(row)<2 for row in scaled):
            raise ValueError('Normalized spectra require bins followed by a parent rate')
        scaled=[ratio(row[:-1],row[-1]) for row in scaled]
    reference=scaled[0]
    return np.asarray([scaled,[row-reference for row in scaled],
                       [ratio(row,reference) for row in scaled]])


def estimate(vectors,strata,factors,normalized=False):
    """Delete from one sample/beam stratum at a time; share reference errors.

    The covariance factors retain correlations between width points induced
    by their common reference, as well as bin/parent and weight correlations.
    These remain conditional learned-grid errors, not retraining coverage.
    """
    factors=width_factors(factors,len(vectors))
    if len(strata)!=len(vectors):
        raise ValueError('Every sample needs its own beam-stratum labels')
    groups=[split_ensembles(v,s) for v,s in zip(vectors,strata)]
    boundaries=np.cumsum([0]+[len(group) for group in groups])
    def transform(*means):
        samples=[sum(means[a:b]) for a,b in zip(boundaries[:-1],boundaries[1:])]
        return comparisons(samples,factors,normalized)
    return independent_jackknife_many([row for group in groups for row in group],transform)


def estimate_retrained_reference(reference_replicas,vectors,strata,factors,normalized=False):
    """Average whole retrained references; keep BW errors conditional on grids.

    Each reference entry is the complete cross-section vector from one
    independent integration. Its between-run spread already includes that
    integration's sampling noise, so within-reference batch variances must
    not be added again. The remaining inputs are independent conditional
    split batches, just as in estimate(). Report both uncertainty components.
    """
    factors=width_factors(factors,len(vectors)+1)
    if len(strata)!=len(vectors):
        raise ValueError('Every conditional sample needs its own beam strata')
    groups=[split_ensembles(v,s) for v,s in zip(vectors,strata)]
    boundaries=np.cumsum([1]+[len(group) for group in groups])
    def transform(*means):
        samples=[means[0]]+[sum(means[a:b]) for a,b in zip(boundaries[:-1],boundaries[1:])]
        return comparisons(samples,factors,normalized)
    result=independent_jackknife_many(
        [reference_replicas]+[row for group in groups for row in group],transform)
    count=len(reference_replicas)
    for label,selected in [('reference_retraining',slice(None,count)),('conditional_samples',slice(count,None))]:
        factors=result['covariance_factor'][selected]
        result['mc_error_'+label]=np.where(result['valid'],np.sqrt(np.sum(factors*factors,axis=0)),np.nan)
    return result
