"""tox_data_integration_jsd

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_data_integration_jsd. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.determine_shared_residual_range_expert_c.restype = None
_lib.determine_shared_residual_range_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS = ("abs_residual_pool", "abs_residual_pool_perm", "pool_size", "shared_residual_range", "residual_range_quantile", "ierr",)

_lib.determine_shared_residual_range_c.restype = None
_lib.determine_shared_residual_range_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENTS = ("abs_residual_pool", "pool_size", "shared_residual_range", "residual_range_quantile", "ierr",)

_lib.determine_study_shared_residual_range_expert_c.restype = None
_lib.determine_study_shared_residual_range_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "tmp_abs_residual_pool", "tmp_abs_residual_pool_perm", "shared_residual_range", "residual_range_quantile", "ierr",)

_lib.determine_study_shared_residual_range_c.restype = None
_lib.determine_study_shared_residual_range_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "shared_residual_range", "residual_range_quantile", "ierr",)

_lib.build_residual_histograms_c.restype = None
_lib.build_residual_histograms_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_RESIDUAL_HISTOGRAMS_ARGUMENTS = ("neighborhood_residuals", "n_reps", "n_neighbors", "n_points", "shared_residual_range", "n_bins", "counts", "pmf", "included_n_reps", "neighbor_mask", "ierr",)

_lib.compute_divergence_per_reference_point_c.restype = None
_lib.compute_divergence_per_reference_point_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENTS = ("pmf_S1", "pmf_S2", "n_points", "n_bins", "js_divergences", "ierr",)

_lib.compute_weighted_global_divergence_c.restype = None
_lib.compute_weighted_global_divergence_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENTS = ("js_divergences", "n_points", "included_n_reps_S1", "included_n_reps_S2", "global_js_divergence", "weights", "ierr",)

def determine_shared_residual_range_expert(
        abs_residual_pool,
        abs_residual_pool_perm,
        residual_range_quantile=95.0,
):
    r"""Compute the shared residual range [-R, R] from a pooled set of absolute residuals

    Parameters
    ----------
    abs_residual_pool : np.ndarray[np.float64] of shape (pool_size,)
        The absolute residual values of the concatenated S1,S2 residuals
        NaN is permitted for this value.
    abs_residual_pool_perm : np.ndarray[np.int32] of shape (pool_size,)
        The permutation vector that sorts `abs_residual_pool`
        The minimum valid value is `1`.
        The maximum valid value is `pool_size`.
    residual_range_quantile : float, optional, default 95.0
        Quantile for determining the residual range
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `95.0`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::determine_shared_residual_range`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        abs_residual_pool = np.ascontiguousarray(abs_residual_pool, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'abs_residual_pool' must be an array of np.float64: {error}") from None
    if abs_residual_pool.ndim != 1:
        raise ValueError(f"'abs_residual_pool' must have 1 dimension, but has {abs_residual_pool.ndim}")
    try:
        abs_residual_pool_perm = np.ascontiguousarray(abs_residual_pool_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'abs_residual_pool_perm' must be an array of np.int32: {error}") from None
    if abs_residual_pool_perm.ndim != 1:
        raise ValueError(f"'abs_residual_pool_perm' must have 1 dimension, but has {abs_residual_pool_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    pool_size = abs_residual_pool.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if abs_residual_pool_perm.shape[0] != pool_size:
        raise ValueError(f"'abs_residual_pool_perm' has {abs_residual_pool_perm.shape[0]} along axis 0, but "
            f"'abs_residual_pool' implies pool_size == {pool_size}"
        )

    # outputs and work arrays, which the caller never sees
    shared_residual_range = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.determine_shared_residual_range_expert_c(
        abs_residual_pool,
        abs_residual_pool_perm,
        ctypes.byref(ctypes.c_int(pool_size)),
        ctypes.byref(shared_residual_range),
        ctypes.byref(ctypes.c_double(residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS)

    return shared_residual_range.value

def determine_shared_residual_range(
        abs_residual_pool,
        residual_range_quantile=95.0,
):
    r"""Compute the shared residual range [-R, R] from a pooled set of absolute residuals

    Parameters
    ----------
    abs_residual_pool : np.ndarray[np.float64] of shape (pool_size,)
        The absolute residual values of the concatenated S1,S2 residuals
        NaN is permitted for this value.
    residual_range_quantile : float, optional, default 95.0
        Quantile for determining the residual range
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `95.0`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::determine_shared_residual_range_alloc`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        abs_residual_pool = np.ascontiguousarray(abs_residual_pool, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'abs_residual_pool' must be an array of np.float64: {error}") from None
    if abs_residual_pool.ndim != 1:
        raise ValueError(f"'abs_residual_pool' must have 1 dimension, but has {abs_residual_pool.ndim}")

    # what the inputs already say, rather than asking for it again
    pool_size = abs_residual_pool.shape[0]

    # outputs and work arrays, which the caller never sees
    shared_residual_range = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.determine_shared_residual_range_c(
        abs_residual_pool,
        ctypes.byref(ctypes.c_int(pool_size)),
        ctypes.byref(shared_residual_range),
        ctypes.byref(ctypes.c_double(residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENTS)

    return shared_residual_range.value

def determine_study_shared_residual_range_expert(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        residual_range_quantile=95.0,
):
    r"""Compute the shared residual range [-R, R] from the neighborhood residuals of two studies

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    residual_range_quantile : float, optional, default 95.0
        Quantile for determining the residual range
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `95.0`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::determine_study_shared_residual_range`, whose argument names are
    the ones an error message reports.
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
    tmp_abs_residual_pool = np.empty(((n_reps_S1 + n_reps_S2)*n_neighbors*n_points,), dtype=np.float64, order='C')
    tmp_abs_residual_pool_perm = np.empty(((n_reps_S1 + n_reps_S2)*n_neighbors*n_points,), dtype=np.int32, order='C')
    shared_residual_range = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.determine_study_shared_residual_range_expert_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        tmp_abs_residual_pool,
        tmp_abs_residual_pool_perm,
        ctypes.byref(shared_residual_range),
        ctypes.byref(ctypes.c_double(residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS)

    return shared_residual_range.value

def determine_study_shared_residual_range(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        residual_range_quantile=95.0,
):
    r"""Compute the shared residual range [-R, R] from the neighborhood residuals of two studies

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    residual_range_quantile : float, optional, default 95.0
        Quantile for determining the residual range
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `95.0`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::determine_study_shared_residual_range_alloc`, whose argument names are
    the ones an error message reports.
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
    shared_residual_range = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.determine_study_shared_residual_range_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(shared_residual_range),
        ctypes.byref(ctypes.c_double(residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_ARGUMENTS)

    return shared_residual_range.value

def build_residual_histograms(
        neighborhood_residuals,
        shared_residual_range,
        n_bins,
        neighbor_mask=None,
):
    r"""Summarize the neighborhood residuals in absolute histogram counts and probability mass functions

    Parameters
    ----------
    neighborhood_residuals : np.ndarray[np.float64] of shape (n_reps, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for a study, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    shared_residual_range : float
        Computed residual range (R)
        The minimum valid value is `0.0`.
    n_bins : int
        Number of equally sized histogram bins in range [-R,R]
    neighbor_mask : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
            Absolute counts of a residual per bin
        pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
            `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        included_n_reps : np.ndarray[np.int32] of shape (n_points,)
            Stores the count of non-NaN replicates (included ones)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::build_residual_histograms`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals = np.asfortranarray(neighborhood_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals' must be an array of np.float64: {error}") from None
    if neighborhood_residuals.ndim != 3:
        raise ValueError(f"'neighborhood_residuals' must have 3 dimensions, but has {neighborhood_residuals.ndim}")
    if neighbor_mask is not None:
        try:
            neighbor_mask = np.asfortranarray(neighbor_mask, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask' must be an array of np.bool_: {error}") from None
        if neighbor_mask.ndim != 2:
            raise ValueError(f"'neighbor_mask' must have 2 dimensions, but has {neighbor_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps = neighborhood_residuals.shape[0]
    n_neighbors = neighborhood_residuals.shape[1]
    n_points = neighborhood_residuals.shape[2]

    # outputs and work arrays, which the caller never sees
    counts = np.empty((n_points, n_bins,), dtype=np.int32, order='F')
    pmf = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.build_residual_histograms_c(
        neighborhood_residuals,
        ctypes.byref(ctypes.c_int(n_reps)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(ctypes.c_int(n_bins)),
        counts,
        pmf,
        included_n_reps,
        neighbor_mask,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BUILD_RESIDUAL_HISTOGRAMS_ARGUMENTS)

    return {
        "counts": counts,
        "pmf": pmf,
        "included_n_reps": included_n_reps,
    }

def compute_divergence_per_reference_point(
        pmf_S1,
        pmf_S2,
):
    r"""Compute the Jensen-Shannon divergence per reference point from two histograms

    Parameters
    ----------
    pmf_S1 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
        Computed normalized histogram counts for study 1
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    pmf_S2 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
        Computed normalized histogram counts for study 2
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    js_divergences : np.ndarray[np.float64] of shape (n_points,)
        Jensen-Shannon divergence per reference point

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::compute_divergence_per_reference_point`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pmf_S1 = np.asfortranarray(pmf_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pmf_S1' must be an array of np.float64: {error}") from None
    if pmf_S1.ndim != 2:
        raise ValueError(f"'pmf_S1' must have 2 dimensions, but has {pmf_S1.ndim}")
    try:
        pmf_S2 = np.asfortranarray(pmf_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pmf_S2' must be an array of np.float64: {error}") from None
    if pmf_S2.ndim != 2:
        raise ValueError(f"'pmf_S2' must have 2 dimensions, but has {pmf_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = pmf_S1.shape[0]
    n_bins = pmf_S1.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if pmf_S2.shape[0] != n_points:
        raise ValueError(f"'pmf_S2' has {pmf_S2.shape[0]} along axis 0, but "
            f"'pmf_S1' implies n_points == {n_points}"
        )
    if pmf_S2.shape[1] != n_bins:
        raise ValueError(f"'pmf_S2' has {pmf_S2.shape[1]} along axis 1, but "
            f"'pmf_S1' implies n_bins == {n_bins}"
        )

    # outputs and work arrays, which the caller never sees
    js_divergences = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_divergence_per_reference_point_c(
        pmf_S1,
        pmf_S2,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_bins)),
        js_divergences,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENTS)

    return js_divergences

def compute_weighted_global_divergence(
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
):
    r"""Compute the global weighted Jensen-Shannon divergence from the per-neighbor divergences

    Parameters
    ----------
    js_divergences : np.ndarray[np.float64] of shape (n_points,)
        Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        The minimum valid value is `0.0`.
    included_n_reps_S1 : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN residuals (included ones) in study 1
        The minimum valid value is `0`.
    included_n_reps_S2 : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN residuals (included ones) in study 2
        The minimum valid value is `0`.

    Returns
    -------
    dict
        with keys:

        global_js_divergence : float
            Weighted global Jensen-Shannon divergence
        weights : np.ndarray[np.float64] of shape (n_points,)
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::compute_weighted_global_divergence`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        js_divergences = np.ascontiguousarray(js_divergences, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'js_divergences' must be an array of np.float64: {error}") from None
    if js_divergences.ndim != 1:
        raise ValueError(f"'js_divergences' must have 1 dimension, but has {js_divergences.ndim}")
    try:
        included_n_reps_S1 = np.ascontiguousarray(included_n_reps_S1, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'included_n_reps_S1' must be an array of np.int32: {error}") from None
    if included_n_reps_S1.ndim != 1:
        raise ValueError(f"'included_n_reps_S1' must have 1 dimension, but has {included_n_reps_S1.ndim}")
    try:
        included_n_reps_S2 = np.ascontiguousarray(included_n_reps_S2, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'included_n_reps_S2' must be an array of np.int32: {error}") from None
    if included_n_reps_S2.ndim != 1:
        raise ValueError(f"'included_n_reps_S2' must have 1 dimension, but has {included_n_reps_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = js_divergences.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if included_n_reps_S1.shape[0] != n_points:
        raise ValueError(f"'included_n_reps_S1' has {included_n_reps_S1.shape[0]} along axis 0, but "
            f"'js_divergences' implies n_points == {n_points}"
        )
    if included_n_reps_S2.shape[0] != n_points:
        raise ValueError(f"'included_n_reps_S2' has {included_n_reps_S2.shape[0]} along axis 0, but "
            f"'js_divergences' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    global_js_divergence = ctypes.c_double(0)
    weights = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_weighted_global_divergence_c(
        js_divergences,
        ctypes.byref(ctypes.c_int(n_points)),
        included_n_reps_S1,
        included_n_reps_S2,
        ctypes.byref(global_js_divergence),
        weights,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENTS)

    return {
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
    }
