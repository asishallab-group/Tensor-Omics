r"""tox_data_integration_jsd

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation

The step that turns neighborhood residuals -- as produced by the preprocessing module -- into
a JSD value: residuals are binned into histograms per study, the two distributions are
compared per reference point, and the per-point divergences are weighted into one global
figure for the pair of studies.

Python binding, generated from tox_data_integration_jsd. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.determine_shared_residual_range_c.restype = None
_lib.determine_shared_residual_range_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENTS = ("abs_residual_pool", "pool_size", "shared_residual_range", "residual_range_quantile_level", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENT_SOURCES = (None, "abs_residual_pool", None, None, None,)

_lib.determine_shared_residual_range_expert_c.restype = None
_lib.determine_shared_residual_range_expert_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS = ("abs_residual_pool", "abs_residual_pool_perm", "pool_size", "shared_residual_range", "residual_range_quantile_level", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENT_SOURCES = (None, None, "abs_residual_pool", None, None, None,)

_lib.determine_study_shared_residual_range_c.restype = None
_lib.determine_study_shared_residual_range_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "shared_residual_range", "residual_range_quantile_level", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_ARGUMENT_SOURCES = (None, None, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", None, None, None,)

_lib.build_residual_histograms_c.restype = None
_lib.build_residual_histograms_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_RESIDUAL_HISTOGRAMS_ARGUMENTS = ("neighborhood_residuals", "n_reps", "n_neighbors", "n_points", "shared_residual_range", "n_bins", "counts", "pmf", "included_n_reps", "neighbor_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_BUILD_RESIDUAL_HISTOGRAMS_ARGUMENT_SOURCES = (None, "neighborhood_residuals", "neighborhood_residuals", "neighborhood_residuals", None, "counts", None, None, None, None, None,)

_lib.compute_divergence_per_reference_point_c.restype = None
_lib.compute_divergence_per_reference_point_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENTS = ("pmf_S1", "pmf_S2", "n_points", "n_bins", "js_divergences", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENT_SOURCES = (None, None, "pmf_S1", "pmf_S1", None, None,)

_lib.compute_weighted_global_divergence_c.restype = None
_lib.compute_weighted_global_divergence_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENTS = ("js_divergences", "n_points", "included_n_reps_S1", "included_n_reps_S2", "global_js_divergence", "weights", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENT_SOURCES = (None, "js_divergences", None, None, None, None, None,)

def determine_shared_residual_range(
        abs_residual_pool,
        residual_range_quantile_level=0.95,
):
    r"""Compute the shared residual range [-R, R] from a pooled set of absolute residuals

    This takes the pool already built; `determine_study_shared_residual_range` builds it from
    the neighborhood residuals of two studies first, if that is what is at hand.

    Parameters
    ----------
    abs_residual_pool : np.ndarray[np.float64] of shape (pool_size,)
        The absolute residual values of the concatenated S1,S2 residuals
        NaN is permitted for this value.
    residual_range_quantile_level : float, optional, default 0.95
        Quantile level in [0,1] for determining the residual range
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::determine_shared_residual_range`, whose argument names are
    the ones an error message reports.

    This entry point seeds `abs_residual_pool_perm` and sorts it by `abs_residual_pool`.
    Call `determine_shared_residual_range_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        abs_residual_pool = _np.ascontiguousarray(abs_residual_pool, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'abs_residual_pool' must be an array of np.float64: {_error}", argument="abs_residual_pool") from None
    if abs_residual_pool.ndim != 1:
        raise _ToxArgumentError(f"'abs_residual_pool' must have 1 dimension, but has {abs_residual_pool.ndim}", argument="abs_residual_pool")

    # what the inputs already say, rather than asking for it again
    pool_size = abs_residual_pool.shape[0]

    # outputs and work arrays, which the caller never sees
    shared_residual_range = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.determine_shared_residual_range_c(
        abs_residual_pool,
        _ctypes.byref(_ctypes.c_int(pool_size)),
        _ctypes.byref(shared_residual_range),
        _ctypes.byref(_ctypes.c_double(residual_range_quantile_level)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENTS, _DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENT_SOURCES)

    return shared_residual_range.value

def determine_shared_residual_range_expert(
        abs_residual_pool,
        abs_residual_pool_perm,
        residual_range_quantile_level=0.95,
):
    r"""Compute the shared residual range [-R, R] from a pooled set of absolute residuals

    This takes the pool already built; `determine_study_shared_residual_range` builds it from
    the neighborhood residuals of two studies first, if that is what is at hand.

    Parameters
    ----------
    abs_residual_pool : np.ndarray[np.float64] of shape (pool_size,)
        The absolute residual values of the concatenated S1,S2 residuals
        NaN is permitted for this value.
    abs_residual_pool_perm : np.ndarray[np.int32] of shape (pool_size,)
        The permutation vector that sorts `abs_residual_pool`
        The minimum valid value is `1`.
        The maximum valid value is `pool_size`.
    residual_range_quantile_level : float, optional, default 0.95
        Quantile level in [0,1] for determining the residual range
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::determine_shared_residual_range_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `abs_residual_pool_perm` yourself.
    `determine_shared_residual_range` seeds `abs_residual_pool_perm` and sorts it by `abs_residual_pool`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        abs_residual_pool = _np.ascontiguousarray(abs_residual_pool, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'abs_residual_pool' must be an array of np.float64: {_error}", argument="abs_residual_pool") from None
    if abs_residual_pool.ndim != 1:
        raise _ToxArgumentError(f"'abs_residual_pool' must have 1 dimension, but has {abs_residual_pool.ndim}", argument="abs_residual_pool")
    try:
        abs_residual_pool_perm = _np.ascontiguousarray(abs_residual_pool_perm, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'abs_residual_pool_perm' must be an array of np.int32: {_error}", argument="abs_residual_pool_perm") from None
    if abs_residual_pool_perm.ndim != 1:
        raise _ToxArgumentError(f"'abs_residual_pool_perm' must have 1 dimension, but has {abs_residual_pool_perm.ndim}", argument="abs_residual_pool_perm")

    # what the inputs already say, rather than asking for it again
    pool_size = abs_residual_pool.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if abs_residual_pool_perm.shape[0] != pool_size:
        raise _ToxArgumentError(f"'abs_residual_pool_perm' has {abs_residual_pool_perm.shape[0]} along axis 0, but 'abs_residual_pool' implies pool_size == {pool_size}; 'abs_residual_pool_perm' has shape (pool_size)", argument="abs_residual_pool_perm")

    # outputs and work arrays, which the caller never sees
    shared_residual_range = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.determine_shared_residual_range_expert_c(
        abs_residual_pool,
        abs_residual_pool_perm,
        _ctypes.byref(_ctypes.c_int(pool_size)),
        _ctypes.byref(shared_residual_range),
        _ctypes.byref(_ctypes.c_double(residual_range_quantile_level)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS, _DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENT_SOURCES)

    return shared_residual_range.value

def determine_study_shared_residual_range(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        residual_range_quantile_level=0.95,
):
    r"""Compute the shared residual range [-R, R] from the neighborhood residuals of two studies

    Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
    as `determine_shared_residual_range` does.

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    residual_range_quantile_level : float, optional, default 0.95
        Quantile level in [0,1] for determining the residual range
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::determine_study_shared_residual_range`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals_S1 = _np.asfortranarray(neighborhood_residuals_S1, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_residuals_S1' must be an array of np.float64: {_error}", argument="neighborhood_residuals_S1") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise _ToxArgumentError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}", argument="neighborhood_residuals_S1")
    try:
        neighborhood_residuals_S2 = _np.asfortranarray(neighborhood_residuals_S2, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' must be an array of np.float64: {_error}", argument="neighborhood_residuals_S2") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}", argument="neighborhood_residuals_S2")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighborhood_residuals_S2' has shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighborhood_residuals_S2' has shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")

    # outputs and work arrays, which the caller never sees
    shared_residual_range = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.determine_study_shared_residual_range_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        _ctypes.byref(_ctypes.c_int(n_reps_S1)),
        _ctypes.byref(_ctypes.c_int(n_reps_S2)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(shared_residual_range),
        _ctypes.byref(_ctypes.c_double(residual_range_quantile_level)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_ARGUMENTS, _DETERMINE_STUDY_SHARED_RESIDUAL_RANGE_ARGUMENT_SOURCES)

    return shared_residual_range.value

def build_residual_histograms(
        neighborhood_residuals,
        shared_residual_range,
        n_bins,
        neighbor_mask=None,
):
    r"""Summarize the neighborhood residuals in absolute histogram counts and probability mass functions

    The probability mass function `pmf(residual, bin)` is actually a matrix.

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

        counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F'), read-only
            Absolute counts of a residual per bin
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
            `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps : np.ndarray[np.int32] of shape (n_points,), read-only
            Stores the count of non-NaN replicates (included ones)
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::build_residual_histograms`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals = _np.asfortranarray(neighborhood_residuals, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_residuals' must be an array of np.float64: {_error}", argument="neighborhood_residuals") from None
    if neighborhood_residuals.ndim != 3:
        raise _ToxArgumentError(f"'neighborhood_residuals' must have 3 dimensions, but has {neighborhood_residuals.ndim}", argument="neighborhood_residuals")
    if neighbor_mask is not None:
        try:
            neighbor_mask = _np.asfortranarray(neighbor_mask, dtype=_np.bool_)
        except (TypeError, ValueError) as _error:
            raise _ToxArgumentError(f"'neighbor_mask' must be an array of np.bool_: {_error}", argument="neighbor_mask") from None
        if neighbor_mask.ndim != 2:
            raise _ToxArgumentError(f"'neighbor_mask' must have 2 dimensions, but has {neighbor_mask.ndim}", argument="neighbor_mask")

    # what the inputs already say, rather than asking for it again
    n_reps = neighborhood_residuals.shape[0]
    n_neighbors = neighborhood_residuals.shape[1]
    n_points = neighborhood_residuals.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if neighbor_mask is not None and neighbor_mask.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask' has {neighbor_mask.shape[0]} along axis 0, but 'neighborhood_residuals' implies n_neighbors == {n_neighbors}; 'neighbor_mask' has shape (n_neighbors, n_points)", argument="neighbor_mask")
    if neighbor_mask is not None and neighbor_mask.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask' has {neighbor_mask.shape[1]} along axis 1, but 'neighborhood_residuals' implies n_points == {n_points}; 'neighbor_mask' has shape (n_neighbors, n_points)", argument="neighbor_mask")

    # outputs and work arrays, which the caller never sees
    counts = _np.empty((n_points, n_bins,), dtype=_np.int32, order='F')
    pmf = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    included_n_reps = _np.empty((n_points,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.build_residual_histograms_c(
        neighborhood_residuals,
        _ctypes.byref(_ctypes.c_int(n_reps)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_double(shared_residual_range)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        counts,
        pmf,
        included_n_reps,
        neighbor_mask,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BUILD_RESIDUAL_HISTOGRAMS_ARGUMENTS, _BUILD_RESIDUAL_HISTOGRAMS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    counts.flags.writeable = False
    pmf.flags.writeable = False
    included_n_reps.flags.writeable = False

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

    Takes the probabilities `pmf` produced by `build_residual_histograms`.

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
    js_divergences : np.ndarray[np.float64] of shape (n_points,), read-only
        Jensen-Shannon divergence per reference point
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::compute_divergence_per_reference_point`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pmf_S1 = _np.asfortranarray(pmf_S1, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'pmf_S1' must be an array of np.float64: {_error}", argument="pmf_S1") from None
    if pmf_S1.ndim != 2:
        raise _ToxArgumentError(f"'pmf_S1' must have 2 dimensions, but has {pmf_S1.ndim}", argument="pmf_S1")
    try:
        pmf_S2 = _np.asfortranarray(pmf_S2, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'pmf_S2' must be an array of np.float64: {_error}", argument="pmf_S2") from None
    if pmf_S2.ndim != 2:
        raise _ToxArgumentError(f"'pmf_S2' must have 2 dimensions, but has {pmf_S2.ndim}", argument="pmf_S2")

    # what the inputs already say, rather than asking for it again
    n_points = pmf_S1.shape[0]
    n_bins = pmf_S1.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if pmf_S2.shape[0] != n_points:
        raise _ToxArgumentError(f"'pmf_S2' has {pmf_S2.shape[0]} along axis 0, but 'pmf_S1' implies n_points == {n_points}; 'pmf_S2' has shape (n_points, n_bins)", argument="pmf_S2")
    if pmf_S2.shape[1] != n_bins:
        raise _ToxArgumentError(f"'pmf_S2' has {pmf_S2.shape[1]} along axis 1, but 'pmf_S1' implies n_bins == {n_bins}; 'pmf_S2' has shape (n_points, n_bins)", argument="pmf_S2")

    # outputs and work arrays, which the caller never sees
    js_divergences = _np.empty((n_points,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_divergence_per_reference_point_c(
        pmf_S1,
        pmf_S2,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        js_divergences,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENTS, _COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    js_divergences.flags.writeable = False

    return js_divergences

def compute_weighted_global_divergence(
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
):
    r"""Compute the global weighted Jensen-Shannon divergence from the per-neighbor divergences

    Takes the divergences produced by `compute_divergence_per_reference_point`.

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
        weights : np.ndarray[np.float64] of shape (n_points,), read-only
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_jsd::compute_weighted_global_divergence`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        js_divergences = _np.ascontiguousarray(js_divergences, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'js_divergences' must be an array of np.float64: {_error}", argument="js_divergences") from None
    if js_divergences.ndim != 1:
        raise _ToxArgumentError(f"'js_divergences' must have 1 dimension, but has {js_divergences.ndim}", argument="js_divergences")
    try:
        included_n_reps_S1 = _np.ascontiguousarray(included_n_reps_S1, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'included_n_reps_S1' must be an array of np.int32: {_error}", argument="included_n_reps_S1") from None
    if included_n_reps_S1.ndim != 1:
        raise _ToxArgumentError(f"'included_n_reps_S1' must have 1 dimension, but has {included_n_reps_S1.ndim}", argument="included_n_reps_S1")
    try:
        included_n_reps_S2 = _np.ascontiguousarray(included_n_reps_S2, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'included_n_reps_S2' must be an array of np.int32: {_error}", argument="included_n_reps_S2") from None
    if included_n_reps_S2.ndim != 1:
        raise _ToxArgumentError(f"'included_n_reps_S2' must have 1 dimension, but has {included_n_reps_S2.ndim}", argument="included_n_reps_S2")

    # what the inputs already say, rather than asking for it again
    n_points = js_divergences.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if included_n_reps_S1.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps_S1' has {included_n_reps_S1.shape[0]} along axis 0, but 'js_divergences' implies n_points == {n_points}; 'included_n_reps_S1' has shape (n_points)", argument="included_n_reps_S1")
    if included_n_reps_S2.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps_S2' has {included_n_reps_S2.shape[0]} along axis 0, but 'js_divergences' implies n_points == {n_points}; 'included_n_reps_S2' has shape (n_points)", argument="included_n_reps_S2")

    # outputs and work arrays, which the caller never sees
    global_js_divergence = _ctypes.c_double(0)
    weights = _np.empty((n_points,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_weighted_global_divergence_c(
        js_divergences,
        _ctypes.byref(_ctypes.c_int(n_points)),
        included_n_reps_S1,
        included_n_reps_S2,
        _ctypes.byref(global_js_divergence),
        weights,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENTS, _COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    weights.flags.writeable = False

    return {
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
    }
