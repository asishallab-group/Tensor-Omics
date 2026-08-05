"""Python binding to Generated from the kernel; do not edit -- regenerate instead.

Generated from tox_data_integration_stats. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.gjct_permutation_test_expert_c.restype = None
_lib.gjct_permutation_test_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GJCT_PERMUTATION_TEST_EXPERT_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "tmp_residuals_S1", "tmp_residuals_S2", "tmp_pool", "tmp_pmf_S1", "tmp_pmf_S2", "tmp_counts", "tmp_included_n_reps_S1", "tmp_included_n_reps_S2", "tmp_js_divergences", "tmp_weights", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2", "ierr",)

_lib.gjct_permutation_test_c.restype = None
_lib.gjct_permutation_test_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GJCT_PERMUTATION_TEST_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2", "ierr",)

def gjct_permutation_test_expert(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        global_jsd_observed,
        n_bins,
        shared_residual_range,
        n_permutations,
        random_seed=None,
        neighbor_mask_S1=None,
        neighbor_mask_S2=None,
):
    r"""Estimate how likely the observed divergence is to occur by chance

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    global_jsd_observed : float
        Observed global JSD value for both studies
    n_bins : int
        Number of equally sized histogram bins used for the studies
    shared_residual_range : float
        Computed residual range for both studies
        The minimum valid value is `0.0`.
    n_permutations : int
        Number of permutations to perform
    random_seed : int, optional
        Seed to use for shuffling
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        jsd_null : np.ndarray[np.float64] of shape (n_permutations,)
            Vector of global divergence values obtained under the null hypothesis
        p_value : float
            Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_stats::gjct_permutation_test`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals_S1 = np.asfortranarray(neighborhood_residuals_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S1' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}")
    try:
        neighborhood_residuals_S2 = np.asfortranarray(neighborhood_residuals_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S2' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}")
    if neighbor_mask_S1 is not None:
        try:
            neighbor_mask_S1 = np.asfortranarray(neighbor_mask_S1, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S1' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S1.ndim != 2:
            raise ValueError(f"'neighbor_mask_S1' must have 2 dimensions, but has {neighbor_mask_S1.ndim}")
    if neighbor_mask_S2 is not None:
        try:
            neighbor_mask_S2 = np.asfortranarray(neighbor_mask_S2, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S2' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S2.ndim != 2:
            raise ValueError(f"'neighbor_mask_S2' must have 2 dimensions, but has {neighbor_mask_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    jsd_null = np.empty((n_permutations,), dtype=np.float64, order='C')
    p_value = ctypes.c_double(0)
    tmp_residuals_S1 = np.empty((n_reps_S1, n_neighbors, n_points,), dtype=np.float64, order='F')
    tmp_residuals_S2 = np.empty((n_reps_S2, n_neighbors, n_points,), dtype=np.float64, order='F')
    tmp_pool = np.empty((n_reps_S1 + n_reps_S2, n_neighbors,), dtype=np.float64, order='F')
    tmp_pmf_S1 = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    tmp_pmf_S2 = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    tmp_counts = np.empty((n_points, n_bins,), dtype=np.int32, order='F')
    tmp_included_n_reps_S1 = np.empty((n_points,), dtype=np.int32, order='C')
    tmp_included_n_reps_S2 = np.empty((n_points,), dtype=np.int32, order='C')
    tmp_js_divergences = np.empty((n_points,), dtype=np.float64, order='C')
    tmp_weights = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.gjct_permutation_test_expert_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(global_jsd_observed)),
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(ctypes.c_int(n_permutations)),
        jsd_null,
        ctypes.byref(p_value),
        tmp_residuals_S1,
        tmp_residuals_S2,
        tmp_pool,
        tmp_pmf_S1,
        tmp_pmf_S2,
        tmp_counts,
        tmp_included_n_reps_S1,
        tmp_included_n_reps_S2,
        tmp_js_divergences,
        tmp_weights,
        None if random_seed is None else ctypes.byref(ctypes.c_int(random_seed)),
        neighbor_mask_S1,
        neighbor_mask_S2,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GJCT_PERMUTATION_TEST_EXPERT_ARGUMENTS)

    return {
        "jsd_null": jsd_null,
        "p_value": p_value.value,
    }

def gjct_permutation_test(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        global_jsd_observed,
        n_bins,
        shared_residual_range,
        n_permutations,
        random_seed=None,
        neighbor_mask_S1=None,
        neighbor_mask_S2=None,
):
    r"""Estimate how likely the observed divergence is to occur by chance

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    global_jsd_observed : float
        Observed global JSD value for both studies
    n_bins : int
        Number of equally sized histogram bins used for the studies
    shared_residual_range : float
        Computed residual range for both studies
        The minimum valid value is `0.0`.
    n_permutations : int
        Number of permutations to perform
    random_seed : int, optional
        Seed to use for shuffling
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        jsd_null : np.ndarray[np.float64] of shape (n_permutations,)
            Vector of global divergence values obtained under the null hypothesis
        p_value : float
            Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_stats::gjct_permutation_test_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals_S1 = np.asfortranarray(neighborhood_residuals_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S1' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}")
    try:
        neighborhood_residuals_S2 = np.asfortranarray(neighborhood_residuals_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S2' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}")
    if neighbor_mask_S1 is not None:
        try:
            neighbor_mask_S1 = np.asfortranarray(neighbor_mask_S1, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S1' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S1.ndim != 2:
            raise ValueError(f"'neighbor_mask_S1' must have 2 dimensions, but has {neighbor_mask_S1.ndim}")
    if neighbor_mask_S2 is not None:
        try:
            neighbor_mask_S2 = np.asfortranarray(neighbor_mask_S2, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S2' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S2.ndim != 2:
            raise ValueError(f"'neighbor_mask_S2' must have 2 dimensions, but has {neighbor_mask_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    jsd_null = np.empty((n_permutations,), dtype=np.float64, order='C')
    p_value = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.gjct_permutation_test_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(global_jsd_observed)),
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(ctypes.c_int(n_permutations)),
        jsd_null,
        ctypes.byref(p_value),
        None if random_seed is None else ctypes.byref(ctypes.c_int(random_seed)),
        neighbor_mask_S1,
        neighbor_mask_S2,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GJCT_PERMUTATION_TEST_ARGUMENTS)

    return {
        "jsd_null": jsd_null,
        "p_value": p_value.value,
    }
