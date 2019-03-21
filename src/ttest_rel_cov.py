# adapted from
# ./site-packages/scipy/stats/stats.py
from collections import namedtuple
import numpy as np
from numpy import ma
from scipy.stats import distributions
# from scipy.stats import _chk2_asarray, _contains_nan, Ttest_relResult
# from scipy.stats import Ttest_relResult

Ttest_relResult = namedtuple('Ttest_relResult', ('statistic', 'pvalue'))

def _chk2_asarray(a, b, axis):
    if axis is None:
        a = np.ravel(a)
        b = np.ravel(b)
        outaxis = 0
    else:
        a = np.asarray(a)
        b = np.asarray(b)
        outaxis = axis

    if a.ndim == 0:
        a = np.atleast_1d(a)
    if b.ndim == 0:
        b = np.atleast_1d(b)

    return a, b, outaxis


def _contains_nan(a, nan_policy='propagate'):
    policies = ['propagate', 'raise', 'omit']
    if nan_policy not in policies:
        raise ValueError("nan_policy must be one of {%s}" %
                         ', '.join("'%s'" % s for s in policies))
    try:
        # Calling np.sum to avoid creating a huge array into memory
        # e.g. np.isnan(a).any()
        with np.errstate(invalid='ignore'):
            contains_nan = np.isnan(np.sum(a))
    except TypeError:
        # If the check cannot be properly performed we fallback to omiting
        # nan values and raising a warning. This can happen when attempting to
        # sum things that are not numbers (e.g. as in the function `mode`).
        contains_nan = False
        nan_policy = 'omit'
        warnings.warn("The input array could not be properly checked for nan "
                      "values. nan values will be ignored.", RuntimeWarning)

    if contains_nan and nan_policy == 'raise':
        raise ValueError("The input contains nan values")

    return (contains_nan, nan_policy)


def _ttest_finish(df, t):
    """Common code between all 3 t-test functions."""
    prob = distributions.t.sf(np.abs(t), df) * 2  # use np.abs to get upper tail
    if t.ndim == 0:
        t = t[()]

    return t, prob


def ttest_rel_cov(a, b, axis=0, nan_policy='propagate', ncov=0):
    """
    scipy.stats.ttest_rel but with nbatch
    """

    a, b, axis = _chk2_asarray(a, b, axis)

    cna, npa = _contains_nan(a, nan_policy)
    cnb, npb = _contains_nan(b, nan_policy)
    contains_nan = cna or cnb
    if npa == 'omit' or npb == 'omit':
        nan_policy = 'omit'

    if contains_nan and nan_policy == 'omit':
        a = ma.masked_invalid(a)
        b = ma.masked_invalid(b)
        m = ma.mask_or(ma.getmask(a), ma.getmask(b))
        aa = ma.array(a, mask=m, copy=True)
        bb = ma.array(b, mask=m, copy=True)
        return mstats_basic.ttest_rel(aa, bb, axis)

    if a.shape[axis] != b.shape[axis]:
        raise ValueError('unequal length arrays')

    if a.size == 0 or b.size == 0:
        return np.nan, np.nan

    n = a.shape[axis]
    df = float(n - 1 - ncov)

    d = (a - b).astype(np.float64)
    v = np.var(d, axis, ddof=1)
    dm = np.mean(d, axis)
    denom = np.sqrt(v / float(n))

    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.divide(dm, denom)
    t, prob = _ttest_finish(df, t)

    return Ttest_relResult(t, prob)
