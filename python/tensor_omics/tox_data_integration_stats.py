r"""tox_data_integration_stats

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test

The permutation test of the JS-Comp-Test:
:func:`tensor_omics.gjct_permutation_test` estimates an
empirical p-value for each of K studies' weighted global JSD against the consensus pmf, by
resampling every reference point's pooled histogram counts into the studies without
replacement. `p_values` is the primary result; `p_values_observed_consensus` uses the older
null (a fixed observed consensus) and is kept for comparison.

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

    In every permutation, each reference point's pooled consensus counts (`mean_pmf_counts`) are
    drawn without replacement (multivariate hypergeometric) into the studies, each study taking
    as many residuals as it observed (`included_n_reps`). The permuted pmfs, their equal-weight
    consensus (built exactly as
    :func:`tensor_omics.create_mean_pmf` builds the
    observed one) and each study's weighted global JSD against that consensus are recomputed, and
    compared with the observed `global_jsd_observed`. `p_values_observed_consensus` comes from the
    same draws, taking each permuted study's JSD against the fixed observed `mean_pmf` instead.

    Both p-values are `(count + 1)/(n_permutations + 1)`, so never exactly 0. With
    `n_permutations = 0` no permutation runs and both are NaN ("not computed").

    Errors (`ierr`): `ERR_INVALID_INPUT`, before any draw, when at some reference point the
    studies' counts do not add up to the pool's total (`mean_pmf_included_n_reps`) or that total
    is not the sum of the pool's bin counts; `ERR_ALLOC_FAIL` when the random number generator
    cannot be created. The p-values are then NaN. A failed draw later on is recorded (first
    failure only) without stopping the remaining permutations.

    Parameters
    ----------
    n_permutations : int
        Number of permutations; `0` skips the test (both p-value outputs are then NaN)
        The minimum valid value is `0`.
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
        Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
        The minimum valid value is `0`.
    mean_pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
        Observed consensus pmf; used only for `p_values_observed_consensus`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates per reference point, per study: the size of each study's
        draw. Every study needs at least one at every reference point
        The minimum valid value is `1`.
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
            Empirical p-value per study, `(count + 1)/(n_permutations + 1)` over the permutations
            whose global JSD reached the observed one; NaN when `n_permutations` is 0
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values_observed_consensus : np.ndarray[np.float64] of shape (n_studies,), read-only
            Like `p_values`, from the same permutations, but against the fixed observed consensus
            `mean_pmf`; kept for comparison only. NaN when `n_permutations` is 0
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
    n_bins = mean_pmf_counts.shape[1]
    n_points = mean_pmf_counts.shape[0]
    n_studies = included_n_reps.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if mean_pmf.shape[1] != n_bins:
        raise ValueError(f"'mean_pmf' has {mean_pmf.shape[1]} along axis 1, but "
            f"'mean_pmf_counts' implies n_bins == {n_bins}"
        )
    if mean_pmf.shape[0] != n_points:
        raise ValueError(f"'mean_pmf' has {mean_pmf.shape[0]} along axis 0, but "
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
