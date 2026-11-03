r"""tox_data_integration_jsd

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation

The step that turns neighborhood residuals -- as produced by the preprocessing module -- into
a JSD value: residuals are binned into histograms per study, the two distributions are
compared per reference point, and the per-point divergences are weighted into one global
figure for the pair of studies. For K studies, each study is compared against the consensus pmf
of all of them, built by
:func:`tensor_omics.create_mean_pmf` (with its pooled
counts) or :func:`tensor_omics.create_mean_pmf_only`.

Python binding, generated from tox_data_integration_jsd. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
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
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_RESIDUAL_HISTOGRAMS_ARGUMENTS = ("neighborhood_residuals", "n_reps", "n_neighbors", "n_points", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins", "n_bins_per_point", "counts", "pmf", "included_n_reps", "neighbor_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_BUILD_RESIDUAL_HISTOGRAMS_ARGUMENT_SOURCES = (None, "neighborhood_residuals", "neighborhood_residuals", "neighborhood_residuals", None, None, "counts", None, None, None, None, None, None,)

_lib.calc_pmf_c.restype = None
_lib.calc_pmf_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_PMF_ARGUMENTS = ("counts", "included_n_reps", "n_points", "n_bins", "pmf", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_PMF_ARGUMENT_SOURCES = (None, None, "counts", "counts", None, None,)

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

_lib.jct_compute_jsd_pipeline_c.restype = None
_lib.jct_compute_jsd_pipeline_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_JCT_COMPUTE_JSD_PIPELINE_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "global_js_divergence", "weights", "pmf_S1", "pmf_S2", "neighbor_mask_S1", "neighbor_mask_S2", "ierr",)
#: For a derived argument, the one the caller passed it in
_JCT_COMPUTE_JSD_PIPELINE_ARGUMENT_SOURCES = (None, None, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", "pmf_S1", None, None, None, None, None, None, None, None, None, None, None,)

_lib.create_mean_pmf_c.restype = None
_lib.create_mean_pmf_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CREATE_MEAN_PMF_ARGUMENTS = ("pmfs", "counts", "n_bins", "n_points", "n_studies", "included_n_reps", "mean_pmf", "mean_pmf_included_n_reps", "mean_pmf_counts", "ierr",)
#: For a derived argument, the one the caller passed it in
_CREATE_MEAN_PMF_ARGUMENT_SOURCES = (None, None, "pmfs", "pmfs", "pmfs", None, None, None, None, None,)

_lib.create_mean_pmf_only_c.restype = None
_lib.create_mean_pmf_only_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CREATE_MEAN_PMF_ONLY_ARGUMENTS = ("pmfs", "n_bins", "n_points", "n_studies", "mean_pmf", "ierr",)
#: For a derived argument, the one the caller passed it in
_CREATE_MEAN_PMF_ONLY_ARGUMENT_SOURCES = (None, "pmfs", "pmfs", "pmfs", None, None,)

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
    # check every input by the input contract, converting only when C needs it
    abs_residual_pool = _contract.reals(abs_residual_pool, "abs_residual_pool", ndim=1)
    residual_range_quantile_level = _contract.real(residual_range_quantile_level, "residual_range_quantile_level")

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
    # check every input by the input contract, converting only when C needs it
    abs_residual_pool = _contract.reals(abs_residual_pool, "abs_residual_pool", ndim=1)
    abs_residual_pool_perm = _contract.integers(abs_residual_pool_perm, "abs_residual_pool_perm", ndim=1)
    residual_range_quantile_level = _contract.real(residual_range_quantile_level, "residual_range_quantile_level")

    # what the inputs already say, rather than asking for it again
    pool_size = abs_residual_pool.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if abs_residual_pool_perm.shape[0] != pool_size:
        raise _ToxArgumentError(f"'abs_residual_pool_perm' has {abs_residual_pool_perm.shape[0]} along axis 0, but 'abs_residual_pool' implies pool_size == {pool_size}; 'abs_residual_pool_perm' is declared with shape (pool_size)", argument="abs_residual_pool_perm")

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
    # check every input by the input contract, converting only when C needs it
    neighborhood_residuals_S1 = _contract.reals(neighborhood_residuals_S1, "neighborhood_residuals_S1", ndim=3, order="F")
    neighborhood_residuals_S2 = _contract.reals(neighborhood_residuals_S2, "neighborhood_residuals_S2", ndim=3, order="F")
    residual_range_quantile_level = _contract.real(residual_range_quantile_level, "residual_range_quantile_level")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighborhood_residuals_S2' is declared with shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighborhood_residuals_S2' is declared with shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")

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
        shared_residual_range_low,
        shared_residual_range_high,
        max_n_bins,
        n_bins_per_point,
        neighbor_mask=None,
):
    r"""Summarize the neighborhood residuals in absolute histogram counts and probability mass functions

    The probability mass function `pmf(residual, bin)` is actually a matrix.

    Parameters
    ----------
    neighborhood_residuals : np.ndarray[np.float64] of shape (n_reps, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for a study, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    shared_residual_range_low : np.ndarray[np.float64] of shape (n_points,)
        Lower bound of the histogram range (R_low) for reference point i_point -- e.g. from
        :func:`tensor_omics.determine_bin_count_occupancy`'s
        own `shared_residual_range_low` output
    shared_residual_range_high : np.ndarray[np.float64] of shape (n_points,)
        Upper bound of the histogram range (R_high) for reference point i_point -- e.g. from
        :func:`tensor_omics.determine_bin_count_occupancy`'s
        own `shared_residual_range_high` output
    max_n_bins : int
        Widest histogram bin count used by any reference point in this call -- the array
        extent `counts`/`pmf` are declared with. A reference point whose own
        `n_bins_per_point` is smaller has its remaining columns zero-padded.
    n_bins_per_point : np.ndarray[np.int32] of shape (n_points,)
        Number of equally sized histogram bins in range [shared_residual_range_low(i_point),
        shared_residual_range_high(i_point)] to use for this reference point
        The minimum valid value is `1`.
        The maximum valid value is `max_n_bins`.
    neighbor_mask : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        counts : np.ndarray[np.int32] of shape (n_points, max_n_bins,), column-major (order='F'), read-only
            Absolute counts of a residual per bin, per reference point. Zero-padded beyond
            `n_bins_per_point(i_point)` for each reference point `i_point`
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmf : np.ndarray[np.float64] of shape (n_points, max_n_bins,), column-major (order='F'), read-only
            `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`.
            Zero-padded beyond `n_bins_per_point(i_point)` for each reference point `i_point`
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
    # check every input by the input contract, converting only when C needs it
    neighborhood_residuals = _contract.reals(neighborhood_residuals, "neighborhood_residuals", ndim=3, order="F")
    shared_residual_range_low = _contract.reals(shared_residual_range_low, "shared_residual_range_low", ndim=1)
    shared_residual_range_high = _contract.reals(shared_residual_range_high, "shared_residual_range_high", ndim=1)
    max_n_bins = _contract.integer(max_n_bins, "max_n_bins")
    n_bins_per_point = _contract.integers(n_bins_per_point, "n_bins_per_point", ndim=1)
    if neighbor_mask is not None:
        neighbor_mask = _contract.logicals(neighbor_mask, "neighbor_mask", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_reps = neighborhood_residuals.shape[0]
    n_neighbors = neighborhood_residuals.shape[1]
    n_points = neighborhood_residuals.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if shared_residual_range_low.shape[0] != n_points:
        raise _ToxArgumentError(f"'shared_residual_range_low' has {shared_residual_range_low.shape[0]} along axis 0, but 'neighborhood_residuals' implies n_points == {n_points}; 'shared_residual_range_low' is declared with shape (n_points)", argument="shared_residual_range_low")
    if shared_residual_range_high.shape[0] != n_points:
        raise _ToxArgumentError(f"'shared_residual_range_high' has {shared_residual_range_high.shape[0]} along axis 0, but 'neighborhood_residuals' implies n_points == {n_points}; 'shared_residual_range_high' is declared with shape (n_points)", argument="shared_residual_range_high")
    if n_bins_per_point.shape[0] != n_points:
        raise _ToxArgumentError(f"'n_bins_per_point' has {n_bins_per_point.shape[0]} along axis 0, but 'neighborhood_residuals' implies n_points == {n_points}; 'n_bins_per_point' is declared with shape (n_points)", argument="n_bins_per_point")
    if neighbor_mask is not None and neighbor_mask.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask' has {neighbor_mask.shape[0]} along axis 0, but 'neighborhood_residuals' implies n_neighbors == {n_neighbors}; 'neighbor_mask' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask")
    if neighbor_mask is not None and neighbor_mask.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask' has {neighbor_mask.shape[1]} along axis 1, but 'neighborhood_residuals' implies n_points == {n_points}; 'neighbor_mask' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask")

    # outputs and work arrays, which the caller never sees
    counts = _np.empty((n_points, max_n_bins,), dtype=_np.int32, order='F')
    pmf = _np.empty((n_points, max_n_bins,), dtype=_np.float64, order='F')
    included_n_reps = _np.empty((n_points,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.build_residual_histograms_c(
        neighborhood_residuals,
        _ctypes.byref(_ctypes.c_int(n_reps)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        shared_residual_range_low,
        shared_residual_range_high,
        _ctypes.byref(_ctypes.c_int(max_n_bins)),
        n_bins_per_point,
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

def calc_pmf(
        counts,
        included_n_reps,
):
    r"""Normalize histogram counts into a probability mass function

    The counts-to-pmf half of `build_residual_histograms_impl`, factored out so a caller that
    already holds a study's `counts` -- untouched by anything that perturbed scratch copies
    downstream of it -- can re-derive `pmf` without re-binning residuals.

    Parameters
    ----------
    counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
        Absolute counts of a residual per bin
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates (included ones) per reference point
        The minimum valid value is `0`.

    Returns
    -------
    pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
        `counts` normalized to `0 <= pmf(:, i) <= 1` and `sum(pmf(:, i)) == 1`
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
    Generated from the Fortran procedure `tox_data_integration_jsd::calc_pmf`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    counts = _contract.integers(counts, "counts", ndim=2, order="F")
    included_n_reps = _contract.integers(included_n_reps, "included_n_reps", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_points = counts.shape[0]
    n_bins = counts.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but 'counts' implies n_points == {n_points}; 'included_n_reps' is declared with shape (n_points)", argument="included_n_reps")

    # outputs and work arrays, which the caller never sees
    pmf = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.calc_pmf_c(
        counts,
        included_n_reps,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        pmf,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_PMF_ARGUMENTS, _CALC_PMF_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    pmf.flags.writeable = False

    return pmf

def compute_divergence_per_reference_point(
        pmf_S1,
        pmf_S2,
):
    r"""Compute the Jensen-Shannon divergence per reference point from two histograms

    Takes the probabilities `pmf` produced by `build_residual_histograms`. The natural JSD
    computed from the KL divergences lies in `[0, ln 2]`; the final step below rescales it by
    `LOG_2` onto `[0, 1]`.

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
        Jensen-Shannon divergence per reference point, rescaled by `LOG_2` onto `[0, 1]`
        (rather than its natural `[0, ln 2]` range)
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
    # check every input by the input contract, converting only when C needs it
    pmf_S1 = _contract.reals(pmf_S1, "pmf_S1", ndim=2, order="F")
    pmf_S2 = _contract.reals(pmf_S2, "pmf_S2", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_points = pmf_S1.shape[0]
    n_bins = pmf_S1.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if pmf_S2.shape[0] != n_points:
        raise _ToxArgumentError(f"'pmf_S2' has {pmf_S2.shape[0]} along axis 0, but 'pmf_S1' implies n_points == {n_points}; 'pmf_S2' is declared with shape (n_points, n_bins)", argument="pmf_S2")
    if pmf_S2.shape[1] != n_bins:
        raise _ToxArgumentError(f"'pmf_S2' has {pmf_S2.shape[1]} along axis 1, but 'pmf_S1' implies n_bins == {n_bins}; 'pmf_S2' is declared with shape (n_points, n_bins)", argument="pmf_S2")

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
    # check every input by the input contract, converting only when C needs it
    js_divergences = _contract.reals(js_divergences, "js_divergences", ndim=1)
    included_n_reps_S1 = _contract.integers(included_n_reps_S1, "included_n_reps_S1", ndim=1)
    included_n_reps_S2 = _contract.integers(included_n_reps_S2, "included_n_reps_S2", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_points = js_divergences.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if included_n_reps_S1.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps_S1' has {included_n_reps_S1.shape[0]} along axis 0, but 'js_divergences' implies n_points == {n_points}; 'included_n_reps_S1' is declared with shape (n_points)", argument="included_n_reps_S1")
    if included_n_reps_S2.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps_S2' has {included_n_reps_S2.shape[0]} along axis 0, but 'js_divergences' implies n_points == {n_points}; 'included_n_reps_S2' is declared with shape (n_points)", argument="included_n_reps_S2")

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

def jct_compute_jsd_pipeline(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        n_bins,
        shared_residual_range,
        neighbor_mask_S1=None,
        neighbor_mask_S2=None,
):
    r"""Compute the JSD of two studies from their neighborhood residuals in one call

    Runs :func:`tensor_omics.build_residual_histograms` for
    each study, :func:`tensor_omics.compute_divergence_per_reference_point`
    and :func:`tensor_omics.compute_weighted_global_divergence`,
    with one bin count `n_bins` and one symmetric range `[-shared_residual_range, shared_residual_range]`
    for every reference point.

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    n_bins : int
        Number of equally sized histogram bins used for the studies
    shared_residual_range : float
        Computed residual range for both studies
        The minimum valid value is `0.0`.
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        js_divergences : np.ndarray[np.float64] of shape (n_points,), read-only
            Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps_S1 : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN residuals (included ones) in study 1
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps_S2 : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN residuals (included ones) in study 2
            A result is a value; call `.copy()` to obtain a modifiable array.
        global_js_divergence : float
            Weighted global Jensen-Shannon divergence
        weights : np.ndarray[np.float64] of shape (n_points,), read-only
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmf_S1 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
            Normalized histogram counts for study 1
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmf_S2 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
            Normalized histogram counts for study 2
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
    Generated from the Fortran procedure `tox_data_integration_jsd::jct_compute_jsd_pipeline`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    neighborhood_residuals_S1 = _contract.reals(neighborhood_residuals_S1, "neighborhood_residuals_S1", ndim=3, order="F")
    neighborhood_residuals_S2 = _contract.reals(neighborhood_residuals_S2, "neighborhood_residuals_S2", ndim=3, order="F")
    n_bins = _contract.integer(n_bins, "n_bins")
    shared_residual_range = _contract.real(shared_residual_range, "shared_residual_range")
    if neighbor_mask_S1 is not None:
        neighbor_mask_S1 = _contract.logicals(neighbor_mask_S1, "neighbor_mask_S1", ndim=2, order="F")
    if neighbor_mask_S2 is not None:
        neighbor_mask_S2 = _contract.logicals(neighbor_mask_S2, "neighbor_mask_S2", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighborhood_residuals_S2' is declared with shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighborhood_residuals_S2' is declared with shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")
    if neighbor_mask_S1 is not None and neighbor_mask_S1.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighbor_mask_S1' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S1")
    if neighbor_mask_S1 is not None and neighbor_mask_S1.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighbor_mask_S1' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S1")
    if neighbor_mask_S2 is not None and neighbor_mask_S2.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighbor_mask_S2' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S2")
    if neighbor_mask_S2 is not None and neighbor_mask_S2.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighbor_mask_S2' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S2")

    # outputs and work arrays, which the caller never sees
    js_divergences = _np.empty((n_points,), dtype=_np.float64, order='C')
    included_n_reps_S1 = _np.empty((n_points,), dtype=_np.int32, order='C')
    included_n_reps_S2 = _np.empty((n_points,), dtype=_np.int32, order='C')
    global_js_divergence = _ctypes.c_double(0)
    weights = _np.empty((n_points,), dtype=_np.float64, order='C')
    pmf_S1 = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    pmf_S2 = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.jct_compute_jsd_pipeline_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        _ctypes.byref(_ctypes.c_int(n_reps_S1)),
        _ctypes.byref(_ctypes.c_int(n_reps_S2)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_double(shared_residual_range)),
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
        _ctypes.byref(global_js_divergence),
        weights,
        pmf_S1,
        pmf_S2,
        neighbor_mask_S1,
        neighbor_mask_S2,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _JCT_COMPUTE_JSD_PIPELINE_ARGUMENTS, _JCT_COMPUTE_JSD_PIPELINE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    js_divergences.flags.writeable = False
    included_n_reps_S1.flags.writeable = False
    included_n_reps_S2.flags.writeable = False
    weights.flags.writeable = False
    pmf_S1.flags.writeable = False
    pmf_S2.flags.writeable = False

    return {
        "js_divergences": js_divergences,
        "included_n_reps_S1": included_n_reps_S1,
        "included_n_reps_S2": included_n_reps_S2,
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
        "pmf_S1": pmf_S1,
        "pmf_S2": pmf_S2,
    }

def create_mean_pmf(
        pmfs,
        counts,
        included_n_reps,
):
    r"""Build the consensus pmf and its histogram counts from all studies' pmfs

    `mean_pmf` is the equal-weight mean of all studies' pmfs, `mean_pmf_counts` and
    `mean_pmf_included_n_reps` the sums of their counts. Known limitation: every study is
    compared against a consensus that includes itself, not a leave-one-out background.

    Parameters
    ----------
    pmfs : np.ndarray[np.float64] of shape (n_points, n_bins, n_studies,), column-major (order='F')
        Per-study probabilities of each bin per reference point, from
        :func:`tensor_omics.build_residual_histograms`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    counts : np.ndarray[np.int32] of shape (n_points, n_bins, n_studies,), column-major (order='F')
        Absolute counts of a residual per bin for `pmfs`
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates per reference point, per study; every study needs at
        least one at every reference point
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        mean_pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
            The consensus pmf, `sum(pmfs, dim=3) / n_studies`
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN replicates per reference point, summed over the studies
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F'), read-only
            Residual counts per bin of the consensus, `sum(counts, dim=3)`
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
    Generated from the Fortran procedure `tox_data_integration_jsd::create_mean_pmf`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    pmfs = _contract.reals(pmfs, "pmfs", ndim=3, order="F")
    counts = _contract.integers(counts, "counts", ndim=3, order="F")
    included_n_reps = _contract.integers(included_n_reps, "included_n_reps", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_bins = pmfs.shape[1]
    n_points = pmfs.shape[0]
    n_studies = pmfs.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if counts.shape[0] != n_points:
        raise _ToxArgumentError(f"'counts' has {counts.shape[0]} along axis 0, but 'pmfs' implies n_points == {n_points}; 'counts' is declared with shape (n_points, n_bins, n_studies)", argument="counts")
    if counts.shape[1] != n_bins:
        raise _ToxArgumentError(f"'counts' has {counts.shape[1]} along axis 1, but 'pmfs' implies n_bins == {n_bins}; 'counts' is declared with shape (n_points, n_bins, n_studies)", argument="counts")
    if counts.shape[2] != n_studies:
        raise _ToxArgumentError(f"'counts' has {counts.shape[2]} along axis 2, but 'pmfs' implies n_studies == {n_studies}; 'counts' is declared with shape (n_points, n_bins, n_studies)", argument="counts")
    if included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but 'pmfs' implies n_points == {n_points}; 'included_n_reps' is declared with shape (n_points, n_studies)", argument="included_n_reps")
    if included_n_reps.shape[1] != n_studies:
        raise _ToxArgumentError(f"'included_n_reps' has {included_n_reps.shape[1]} along axis 1, but 'pmfs' implies n_studies == {n_studies}; 'included_n_reps' is declared with shape (n_points, n_studies)", argument="included_n_reps")

    # outputs and work arrays, which the caller never sees
    mean_pmf = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    mean_pmf_included_n_reps = _np.empty((n_points,), dtype=_np.int32, order='C')
    mean_pmf_counts = _np.empty((n_points, n_bins,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.create_mean_pmf_c(
        pmfs,
        counts,
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_studies)),
        included_n_reps,
        mean_pmf,
        mean_pmf_included_n_reps,
        mean_pmf_counts,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CREATE_MEAN_PMF_ARGUMENTS, _CREATE_MEAN_PMF_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    mean_pmf.flags.writeable = False
    mean_pmf_included_n_reps.flags.writeable = False
    mean_pmf_counts.flags.writeable = False

    return {
        "mean_pmf": mean_pmf,
        "mean_pmf_included_n_reps": mean_pmf_included_n_reps,
        "mean_pmf_counts": mean_pmf_counts,
    }

def create_mean_pmf_only(
        pmfs,
):
    r"""Build only the consensus pmf from all studies' pmfs, without its histogram counts

    As :func:`tensor_omics.create_mean_pmf`, for
    callers that need no counts (the bootstrap and the permutation test). Same known limitation: the consensus includes
    the study compared against it.

    Parameters
    ----------
    pmfs : np.ndarray[np.float64] of shape (n_points, n_bins, n_studies,), column-major (order='F')
        Per-study probabilities of each bin per reference point, from
        :func:`tensor_omics.build_residual_histograms`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    mean_pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
        The consensus pmf, `sum(pmfs, dim=3) / n_studies`
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
    Generated from the Fortran procedure `tox_data_integration_jsd::create_mean_pmf_only`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    pmfs = _contract.reals(pmfs, "pmfs", ndim=3, order="F")

    # what the inputs already say, rather than asking for it again
    n_bins = pmfs.shape[1]
    n_points = pmfs.shape[0]
    n_studies = pmfs.shape[2]

    # outputs and work arrays, which the caller never sees
    mean_pmf = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.create_mean_pmf_only_c(
        pmfs,
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_studies)),
        mean_pmf,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CREATE_MEAN_PMF_ONLY_ARGUMENTS, _CREATE_MEAN_PMF_ONLY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    mean_pmf.flags.writeable = False

    return mean_pmf
