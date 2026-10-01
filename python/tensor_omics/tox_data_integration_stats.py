r"""tox_data_integration_stats

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test

A permutation test estimating an empirical p-value for each study's weighted global JSD
against the consensus pmf. Under the null hypothesis that a study is exchangeable with the
pool, every reference point's pooled consensus histogram counts are repeatedly resampled
without replacement (GSL's multivariate hypergeometric distribution) into each study's own
draw size, the pmfs, their equal-weight consensus and each study's JSD/weighted global JSD
against it are recomputed from that resample, and the observed value is compared against the
resulting null distribution (`p_values`). From the same resamples, a second p-value
(`p_values_observed_consensus`) compares each resampled study against the fixed observed
consensus pmf instead, the null this test used before; it is kept for comparison, and
`p_values` is the primary result.

Generalized to K studies; the pooled resampling pool (`tmp_mean_pmf_counts`) is a work copy
reset every permutation, so the caller's own `mean_pmf_counts` is left untouched.

Python binding, generated from tox_data_integration_stats. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.gjct_permutation_test_c.restype = None
_lib.gjct_permutation_test_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GJCT_PERMUTATION_TEST_ARGUMENTS = ("n_permutations", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf", "mean_pmf_included_n_reps", "included_n_reps", "global_jsd_observed", "p_values", "p_values_observed_consensus", "random_seed", "ierr",)
#: For a derived argument, the one the caller passed it in
_GJCT_PERMUTATION_TEST_ARGUMENT_SOURCES = (None, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", None, None, None, None, None, None, None, None, None,)

def gjct_permutation_test(
        n_permutations,
        mean_pmf_counts,
        mean_pmf,
        mean_pmf_included_n_reps,
        included_n_reps,
        global_jsd_observed,
        random_seed=42,
):
    r"""Estimate how likely each study's observed weighted global JSD is to occur by chance

    Ported from 125-stabilize-jscomp's `gjct_permutation_test_helper`, generalized from its
    hardcoded 2-study shuffle to a K-study GSL `random_multiv_hypergeom` resample from the
    pooled `mean_pmf_counts` (without replacement, per reference point, per study). Impure:
    draws random numbers, so it carries `ierr` for the one genuine runtime failure a validated
    caller cannot foresee -- `create_rng` failing to allocate the GSL generator -- folded in as
    `ERR_ALLOC_FAIL`; a later `random_multiv_hypergeom` failure (which validated,
    internally-consistent inputs should never trigger) is likewise folded in, first failure
    only, without stopping the resampling already in flight, matching
    :func:`tensor_omics.bootstrap_histogram`'s own
    precedent.

    The consensus pmf is recomputed from the permuted pmfs in every permutation, as the
    equal-weight mean of all studies' pmfs, exactly as
    :func:`tensor_omics.create_mean_pmf` builds
    the observed one. So each permuted study's JSD is taken against the consensus of the
    permuted studies, and the null statistic is computed exactly like the observed one. This is
    required because the equal-weight mixture is not invariant under the permutation when the
    studies differ in size: the permuted studies follow the pooled, replicate-weighted
    distribution, whose equal-weight mean moves with every draw, so a fixed observed consensus
    would shift the whole null distribution.

    A second p-value, `p_values_observed_consensus`, is computed from the same permutations (no
    further random draws) with the null this test used before: each permuted study's JSD
    against the fixed observed consensus `mean_pmf`. It is kept to compare the two nulls on data
    with known batch effects; `p_values` stays the primary result.

    Both p-values apply the `(1+count)/(n+1)` Laplace add-one correction --
    `p_values(i) = anint(count(i)+1)/(n_permutations+1)` -- so a study whose observed JSD is
    never reached by any permutation gets `p = 1/(n_permutations+1)`, never `p = 0.0` exactly.

    Parameters
    ----------
    n_permutations : int
        Number of permutations to perform
        The minimum valid value is `0`.
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_bins, n_points,), column-major (order='F')
        Absolute counts of a residual per bin for the consensus pmf -- the pool each
        permutation resamples from without replacement, per reference point
        The minimum valid value is `0`.
    mean_pmf : np.ndarray[np.float64] of shape (n_bins, n_points,), column-major (order='F')
        The consensus pmf built from all studies' pmfs; each permuted study is compared
        against it for `p_values_observed_consensus` only
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates (included ones) per reference point, per study -- how
        many elements are drawn (without replacement) from the pooled pool per study
        The minimum valid value is `0`.
    global_jsd_observed : np.ndarray[np.float64] of shape (n_studies,)
        Observed weighted global JSD of each study against the consensus pmf
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.

    Returns
    -------
    dict
        with keys:

        p_values : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study: the Laplace-corrected fraction of permutations whose
            resampled global JSD reached or exceeded the observed value -- see the correction
            note above
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values_observed_consensus : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study from the same permutations as `p_values`, but with each
            permuted study's JSD taken against the fixed observed consensus `mean_pmf` instead of
            the consensus of the permuted studies -- the null this test used before it recomputed
            the consensus per permutation. Kept for comparison only; `p_values` is the primary
            result. The same Laplace correction applies
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_stats::gjct_permutation_test`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        mean_pmf_counts = np.asfortranarray(mean_pmf_counts, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_pmf_counts' must be an array of np.int32: {error}") from None
    if mean_pmf_counts.ndim != 2:
        raise ValueError(f"'mean_pmf_counts' must have 2 dimensions, but has {mean_pmf_counts.ndim}")
    try:
        mean_pmf = np.asfortranarray(mean_pmf, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_pmf' must be an array of np.float64: {error}") from None
    if mean_pmf.ndim != 2:
        raise ValueError(f"'mean_pmf' must have 2 dimensions, but has {mean_pmf.ndim}")
    try:
        mean_pmf_included_n_reps = np.ascontiguousarray(mean_pmf_included_n_reps, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_pmf_included_n_reps' must be an array of np.int32: {error}") from None
    if mean_pmf_included_n_reps.ndim != 1:
        raise ValueError(f"'mean_pmf_included_n_reps' must have 1 dimension, but has {mean_pmf_included_n_reps.ndim}")
    try:
        included_n_reps = np.asfortranarray(included_n_reps, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'included_n_reps' must be an array of np.int32: {error}") from None
    if included_n_reps.ndim != 2:
        raise ValueError(f"'included_n_reps' must have 2 dimensions, but has {included_n_reps.ndim}")
    try:
        global_jsd_observed = np.ascontiguousarray(global_jsd_observed, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'global_jsd_observed' must be an array of np.float64: {error}") from None
    if global_jsd_observed.ndim != 1:
        raise ValueError(f"'global_jsd_observed' must have 1 dimension, but has {global_jsd_observed.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bins = mean_pmf_counts.shape[0]
    n_points = mean_pmf_counts.shape[1]
    n_studies = included_n_reps.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if mean_pmf.shape[0] != n_bins:
        raise ValueError(f"'mean_pmf' has {mean_pmf.shape[0]} along axis 0, but "
            f"'mean_pmf_counts' implies n_bins == {n_bins}"
        )
    if mean_pmf.shape[1] != n_points:
        raise ValueError(f"'mean_pmf' has {mean_pmf.shape[1]} along axis 1, but "
            f"'mean_pmf_counts' implies n_points == {n_points}"
        )
    if mean_pmf_included_n_reps.shape[0] != n_points:
        raise ValueError(f"'mean_pmf_included_n_reps' has {mean_pmf_included_n_reps.shape[0]} along axis 0, but "
            f"'mean_pmf_counts' implies n_points == {n_points}"
        )
    if included_n_reps.shape[0] != n_points:
        raise ValueError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but "
            f"'mean_pmf_counts' implies n_points == {n_points}"
        )
    if global_jsd_observed.shape[0] != n_studies:
        raise ValueError(f"'global_jsd_observed' has {global_jsd_observed.shape[0]} along axis 0, but "
            f"'included_n_reps' implies n_studies == {n_studies}"
        )

    # outputs and work arrays, which the caller never sees
    p_values = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values_observed_consensus = np.empty((n_studies,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.gjct_permutation_test_c(
        ctypes.byref(ctypes.c_int(n_permutations)),
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_studies)),
        mean_pmf_counts,
        mean_pmf,
        mean_pmf_included_n_reps,
        included_n_reps,
        global_jsd_observed,
        p_values,
        p_values_observed_consensus,
        ctypes.byref(ctypes.c_int(random_seed)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GJCT_PERMUTATION_TEST_ARGUMENTS, _GJCT_PERMUTATION_TEST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    p_values.flags.writeable = False
    p_values_observed_consensus.flags.writeable = False

    return {
        "p_values": p_values,
        "p_values_observed_consensus": p_values_observed_consensus,
    }
