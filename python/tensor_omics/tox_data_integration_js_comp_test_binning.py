r"""tox_data_integration_js_comp_test_binning

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Bin Counts

How many histogram bins each reference point of the JS-Comp-Test gets, chosen from the pooled
residuals of all studies in its neighborhood. The occupancy search
:func:`tensor_omics.determine_bin_count_occupancy`
(the default) takes the largest bin count whose every bin holds at least a minimum number of
residuals; Stone's cross-validation
:func:`tensor_omics.determine_bin_count_stone`
takes the bin count with the smallest estimated risk instead (`bin_count_mode`).
:func:`tensor_omics.determine_bin_count_occupancy_exhaustive`
is the brute-force reference of the occupancy search, and
:func:`tensor_omics.estimate_bin_count`
the Sturges/Freedman-Diaconis rule-of-thumb estimate. Every bin count is at most
``MAX_N_BINS``.

Python binding, generated from tox_data_integration_js_comp_test_binning. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.estimate_bin_count_c.restype = None
_lib.estimate_bin_count_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ESTIMATE_BIN_COUNT_ARGUMENTS = ("residuals", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "shared_residual_range", "n_bins", "sturges_bins", "fd_bins", "ierr",)
#: For a derived argument, the one the caller passed it in
_ESTIMATE_BIN_COUNT_ARGUMENT_SOURCES = (None, "residuals", None, None, None, None, None, None, None,)

_lib.estimate_bin_count_expert_c.restype = None
_lib.estimate_bin_count_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ESTIMATE_BIN_COUNT_EXPERT_ARGUMENTS = ("residuals", "residuals_perm", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "shared_residual_range", "n_bins", "sturges_bins", "fd_bins", "ierr",)
#: For a derived argument, the one the caller passed it in
_ESTIMATE_BIN_COUNT_EXPERT_ARGUMENT_SOURCES = (None, None, "residuals", None, None, None, None, None, None, None,)

_lib.determine_bin_count_occupancy_c.restype = None
_lib.determine_bin_count_occupancy_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_BIN_COUNT_OCCUPANCY_ARGUMENTS = ("pooled_residuals", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "occupancy_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "m_min", "m_max", "min_residuals_per_bin", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_BIN_COUNT_OCCUPANCY_ARGUMENT_SOURCES = (None, "pooled_residuals", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.determine_bin_count_occupancy_expert_c.restype = None
_lib.determine_bin_count_occupancy_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_BIN_COUNT_OCCUPANCY_EXPERT_ARGUMENTS = ("pooled_residuals", "pooled_residuals_perm", "n_residuals", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "occupancy_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "tmp_bin_counts", "m_min", "m_max", "min_residuals_per_bin", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_BIN_COUNT_OCCUPANCY_EXPERT_ARGUMENT_SOURCES = (None, None, "pooled_residuals", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.determine_bin_count_occupancy_exhaustive_c.restype = None
_lib.determine_bin_count_occupancy_exhaustive_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_ARGUMENTS = ("pooled_residuals", "n_residuals", "n_pooled_residuals", "shared_residual_range_low", "shared_residual_range_high", "selected_n_bins", "occupancy_failed", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "m_min", "m_max", "min_residuals_per_bin", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_ARGUMENT_SOURCES = (None, "pooled_residuals", None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.determine_bin_count_occupancy_exhaustive_expert_c.restype = None
_lib.determine_bin_count_occupancy_exhaustive_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_EXPERT_ARGUMENTS = ("pooled_residuals", "pooled_residuals_perm", "n_residuals", "n_pooled_residuals", "shared_residual_range_low", "shared_residual_range_high", "selected_n_bins", "occupancy_failed", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "m_min", "m_max", "min_residuals_per_bin", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_EXPERT_ARGUMENT_SOURCES = (None, None, "pooled_residuals", None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.determine_bin_count_stone_c.restype = None
_lib.determine_bin_count_stone_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_BIN_COUNT_STONE_ARGUMENTS = ("pooled_residuals", "n_residuals", "residual_to_study", "n_studies", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "stone_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "bin_count_mode", "m_min", "m_max", "min_residuals_per_bin", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_BIN_COUNT_STONE_ARGUMENT_SOURCES = (None, "pooled_residuals", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.determine_bin_count_stone_expert_c.restype = None
_lib.determine_bin_count_stone_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_BIN_COUNT_STONE_EXPERT_ARGUMENTS = ("pooled_residuals", "pooled_residuals_perm", "n_residuals", "residual_to_study", "n_studies", "max_n_reps_all_studies", "n_neighbors", "selected_n_bins", "stone_failed", "shared_residual_range_low", "shared_residual_range_high", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "tmp_bin_counts", "tmp_residual_to_study_perm", "tmp_study_block_sizes", "bin_count_mode", "m_min", "m_max", "min_residuals_per_bin", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_DETERMINE_BIN_COUNT_STONE_EXPERT_ARGUMENT_SOURCES = (None, None, "pooled_residuals", None, "tmp_study_block_sizes", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def estimate_bin_count(
        residuals,
        max_n_reps_all_studies,
        n_neighbors,
        shared_residual_range,
):
    r"""Estimate the histogram bin count for one (n_points, n_neighbors) candidate

    Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    (`R = shared_residual_range`), each clamped to
    ``MAX_N_BINS`` and returned as
    `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
    `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
    from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
    bin count.

    Parameters
    ----------
    residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals across all studies, reference points and neighbors
        NaN is permitted for this value.
    max_n_reps_all_studies : int
        Maximum number of replicates across all studies
        The minimum valid value is `1`.
    n_neighbors : int
        Neighborhood size of the candidate under test
        The minimum valid value is `1`.
    shared_residual_range : float
        Computed residual range (R)
        The minimum valid value is `0.0`.

    Returns
    -------
    dict
        with keys:

        n_bins : int
            Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)
        sturges_bins : int
            Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS
        fd_bins : int
            Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
            to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::estimate_bin_count`, whose argument names are
    the ones an error message reports.

    This entry point seeds `residuals_perm` and sorts it by `residuals`.
    Call `estimate_bin_count_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        residuals = np.ascontiguousarray(residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residuals' must be an array of np.float64: {error}") from None
    if residuals.ndim != 1:
        raise ValueError(f"'residuals' must have 1 dimension, but has {residuals.ndim}")

    # what the inputs already say, rather than asking for it again
    n_residuals = residuals.shape[0]

    # outputs and work arrays, which the caller never sees
    n_bins = ctypes.c_int(0)
    sturges_bins = ctypes.c_int(0)
    fd_bins = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.estimate_bin_count_c(
        residuals,
        ctypes.byref(ctypes.c_int(n_residuals)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(n_bins),
        ctypes.byref(sturges_bins),
        ctypes.byref(fd_bins),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ESTIMATE_BIN_COUNT_ARGUMENTS, _ESTIMATE_BIN_COUNT_ARGUMENT_SOURCES)

    return {
        "n_bins": n_bins.value,
        "sturges_bins": sturges_bins.value,
        "fd_bins": fd_bins.value,
    }

def estimate_bin_count_expert(
        residuals,
        residuals_perm,
        max_n_reps_all_studies,
        n_neighbors,
        shared_residual_range,
):
    r"""Estimate the histogram bin count for one (n_points, n_neighbors) candidate

    Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    (`R = shared_residual_range`), each clamped to
    ``MAX_N_BINS`` and returned as
    `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
    `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
    from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
    bin count.

    Parameters
    ----------
    residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals across all studies, reference points and neighbors
        NaN is permitted for this value.
    residuals_perm : np.ndarray[np.int32] of shape (n_residuals,)
        Sorting permutation for `residuals`, ascending, NaN last
        The minimum valid value is `1`.
        The maximum valid value is `n_residuals`.
    max_n_reps_all_studies : int
        Maximum number of replicates across all studies
        The minimum valid value is `1`.
    n_neighbors : int
        Neighborhood size of the candidate under test
        The minimum valid value is `1`.
    shared_residual_range : float
        Computed residual range (R)
        The minimum valid value is `0.0`.

    Returns
    -------
    dict
        with keys:

        n_bins : int
            Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)
        sturges_bins : int
            Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS
        fd_bins : int
            Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
            to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::estimate_bin_count_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `residuals_perm` yourself.
    `estimate_bin_count` seeds `residuals_perm` and sorts it by `residuals`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        residuals = np.ascontiguousarray(residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residuals' must be an array of np.float64: {error}") from None
    if residuals.ndim != 1:
        raise ValueError(f"'residuals' must have 1 dimension, but has {residuals.ndim}")
    try:
        residuals_perm = np.ascontiguousarray(residuals_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residuals_perm' must be an array of np.int32: {error}") from None
    if residuals_perm.ndim != 1:
        raise ValueError(f"'residuals_perm' must have 1 dimension, but has {residuals_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_residuals = residuals.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if residuals_perm.shape[0] != n_residuals:
        raise ValueError(f"'residuals_perm' has {residuals_perm.shape[0]} along axis 0, but "
            f"'residuals' implies n_residuals == {n_residuals}"
        )

    # outputs and work arrays, which the caller never sees
    n_bins = ctypes.c_int(0)
    sturges_bins = ctypes.c_int(0)
    fd_bins = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.estimate_bin_count_expert_c(
        residuals,
        residuals_perm,
        ctypes.byref(ctypes.c_int(n_residuals)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(n_bins),
        ctypes.byref(sturges_bins),
        ctypes.byref(fd_bins),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ESTIMATE_BIN_COUNT_EXPERT_ARGUMENTS, _ESTIMATE_BIN_COUNT_EXPERT_ARGUMENT_SOURCES)

    return {
        "n_bins": n_bins.value,
        "sturges_bins": sturges_bins.value,
        "fd_bins": fd_bins.value,
    }

def determine_bin_count_occupancy(
        pooled_residuals,
        max_n_reps_all_studies,
        n_neighbors,
        m_min=3,
        m_max=120,
        min_residuals_per_bin=10,
        gamma_occupancy=1.25,
        lower_residual_range_quantile=0.05,
        upper_residual_range_quantile=0.95,
):
    r"""Determine one neighborhood's occupancy-constrained histogram bin count (Issue #187)

    For one neighborhood's pooled residuals (all its neighbors, all studies), finds the largest
    bin count `M` in `[m_min, m_max]` whose equal-width histogram over the neighborhood's own range
    `[shared_residual_range_low, shared_residual_range_high]` (the
    `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of these residuals)
    has every bin at or above `min_residuals_per_bin`. A geometric search (factor
    `gamma_occupancy`) brackets the first failing `M`, then every integer below it is tested,
    because occupancy is not monotonic in `M`.

    `occupancy_failed = True` (even `m_min` bins fail, or every residual is NaN) means the caller
    should reject the neighborhood; `selected_n_bins` is then `m_min` only so that it is defined.
    `sturges_bins`/`fd_bins` come from
    :func:`tensor_omics.estimate_bin_count` and are
    diagnostics only.

    Parameters
    ----------
    pooled_residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals for one neighborhood, across all its neighbors and all studies
        NaN is permitted for this value.
    max_n_reps_all_studies : int
        Maximum number of replicates across all studies
        The minimum valid value is `1`.
    n_neighbors : int
        Neighborhood size of the candidate under test
        The minimum valid value is `1`.
    m_min : int, optional, default 3
        Smallest candidate bin count the search will ever test (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
        smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a candidate bin count to
        be admissible (n_min)
        The minimum valid value is `0`.
        The default value is `10`.
    gamma_occupancy : float, optional, default 1.25
        Geometric growth factor for the coarse search stage; must exceed 1 or the search
        never advances
        The minimum valid value is `above(1.0)`.
        The default value is `1.25`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for this neighborhood's own lower residual-range bound
        (shared_residual_range_low)
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for this neighborhood's own upper residual-range bound
        (shared_residual_range_high)
        The minimum valid value is `lower_residual_range_quantile`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    dict
        with keys:

        selected_n_bins : int
            The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
            every bin at or above min_residuals_per_bin; m_min when occupancy_failed
        occupancy_failed : bool
            `True` iff even m_min bins fail the occupancy criterion (or every pooled residual is
            NaN); the caller should then reject this neighborhood
        shared_residual_range_low : float
            This neighborhood's lower residual-range bound (R_low): the
            lower_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
        shared_residual_range_high : float
            This neighborhood's upper residual-range bound (R_high): the
            upper_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
        n_pooled_residuals : int
            Count of non-NaN pooled residuals (N_j)
        min_bin_occupancy : int
            Minimum bin count at selected_n_bins; 0 when occupancy_failed
        mean_bin_occupancy : float
            Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            occupancy_failed
        max_bin_occupancy : int
            Maximum bin count at selected_n_bins; 0 when occupancy_failed
        sturges_bins : int
            Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        fd_bins : int
            Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy`, whose argument names are
    the ones an error message reports.

    This entry point seeds `pooled_residuals_perm` and sorts it by `pooled_residuals` and computes `tmp_bin_counts` for you.
    Call `determine_bin_count_occupancy_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_residuals = np.ascontiguousarray(pooled_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals' must be an array of np.float64: {error}") from None
    if pooled_residuals.ndim != 1:
        raise ValueError(f"'pooled_residuals' must have 1 dimension, but has {pooled_residuals.ndim}")

    # what the inputs already say, rather than asking for it again
    n_residuals = pooled_residuals.shape[0]

    # outputs and work arrays, which the caller never sees
    selected_n_bins = ctypes.c_int(0)
    occupancy_failed = ctypes.c_bool(0)
    shared_residual_range_low = ctypes.c_double(0)
    shared_residual_range_high = ctypes.c_double(0)
    n_pooled_residuals = ctypes.c_int(0)
    min_bin_occupancy = ctypes.c_int(0)
    mean_bin_occupancy = ctypes.c_double(0)
    max_bin_occupancy = ctypes.c_int(0)
    sturges_bins = ctypes.c_int(0)
    fd_bins = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.determine_bin_count_occupancy_c(
        pooled_residuals,
        ctypes.byref(ctypes.c_int(n_residuals)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(selected_n_bins),
        ctypes.byref(occupancy_failed),
        ctypes.byref(shared_residual_range_low),
        ctypes.byref(shared_residual_range_high),
        ctypes.byref(n_pooled_residuals),
        ctypes.byref(min_bin_occupancy),
        ctypes.byref(mean_bin_occupancy),
        ctypes.byref(max_bin_occupancy),
        ctypes.byref(sturges_bins),
        ctypes.byref(fd_bins),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_double(gamma_occupancy)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_BIN_COUNT_OCCUPANCY_ARGUMENTS, _DETERMINE_BIN_COUNT_OCCUPANCY_ARGUMENT_SOURCES)

    return {
        "selected_n_bins": selected_n_bins.value,
        "occupancy_failed": occupancy_failed.value,
        "shared_residual_range_low": shared_residual_range_low.value,
        "shared_residual_range_high": shared_residual_range_high.value,
        "n_pooled_residuals": n_pooled_residuals.value,
        "min_bin_occupancy": min_bin_occupancy.value,
        "mean_bin_occupancy": mean_bin_occupancy.value,
        "max_bin_occupancy": max_bin_occupancy.value,
        "sturges_bins": sturges_bins.value,
        "fd_bins": fd_bins.value,
    }

def determine_bin_count_occupancy_expert(
        pooled_residuals,
        pooled_residuals_perm,
        max_n_reps_all_studies,
        n_neighbors,
        m_min=3,
        m_max=120,
        min_residuals_per_bin=10,
        gamma_occupancy=1.25,
        lower_residual_range_quantile=0.05,
        upper_residual_range_quantile=0.95,
):
    r"""Determine one neighborhood's occupancy-constrained histogram bin count (Issue #187)

    For one neighborhood's pooled residuals (all its neighbors, all studies), finds the largest
    bin count `M` in `[m_min, m_max]` whose equal-width histogram over the neighborhood's own range
    `[shared_residual_range_low, shared_residual_range_high]` (the
    `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of these residuals)
    has every bin at or above `min_residuals_per_bin`. A geometric search (factor
    `gamma_occupancy`) brackets the first failing `M`, then every integer below it is tested,
    because occupancy is not monotonic in `M`.

    `occupancy_failed = True` (even `m_min` bins fail, or every residual is NaN) means the caller
    should reject the neighborhood; `selected_n_bins` is then `m_min` only so that it is defined.
    `sturges_bins`/`fd_bins` come from
    :func:`tensor_omics.estimate_bin_count` and are
    diagnostics only.

    Parameters
    ----------
    pooled_residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals for one neighborhood, across all its neighbors and all studies
        NaN is permitted for this value.
    pooled_residuals_perm : np.ndarray[np.int32] of shape (n_residuals,)
        Sorting permutation for `pooled_residuals`, ascending, NaN last
        The minimum valid value is `1`.
        The maximum valid value is `n_residuals`.
    max_n_reps_all_studies : int
        Maximum number of replicates across all studies
        The minimum valid value is `1`.
    n_neighbors : int
        Neighborhood size of the candidate under test
        The minimum valid value is `1`.
    m_min : int, optional, default 3
        Smallest candidate bin count the search will ever test (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
        smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a candidate bin count to
        be admissible (n_min)
        The minimum valid value is `0`.
        The default value is `10`.
    gamma_occupancy : float, optional, default 1.25
        Geometric growth factor for the coarse search stage; must exceed 1 or the search
        never advances
        The minimum valid value is `above(1.0)`.
        The default value is `1.25`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for this neighborhood's own lower residual-range bound
        (shared_residual_range_low)
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for this neighborhood's own upper residual-range bound
        (shared_residual_range_high)
        The minimum valid value is `lower_residual_range_quantile`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    dict
        with keys:

        selected_n_bins : int
            The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
            every bin at or above min_residuals_per_bin; m_min when occupancy_failed
        occupancy_failed : bool
            `True` iff even m_min bins fail the occupancy criterion (or every pooled residual is
            NaN); the caller should then reject this neighborhood
        shared_residual_range_low : float
            This neighborhood's lower residual-range bound (R_low): the
            lower_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
        shared_residual_range_high : float
            This neighborhood's upper residual-range bound (R_high): the
            upper_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
        n_pooled_residuals : int
            Count of non-NaN pooled residuals (N_j)
        min_bin_occupancy : int
            Minimum bin count at selected_n_bins; 0 when occupancy_failed
        mean_bin_occupancy : float
            Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            occupancy_failed
        max_bin_occupancy : int
            Maximum bin count at selected_n_bins; 0 when occupancy_failed
        sturges_bins : int
            Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        fd_bins : int
            Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `pooled_residuals_perm` and `tmp_bin_counts` yourself.
    `determine_bin_count_occupancy` seeds `pooled_residuals_perm` and sorts it by `pooled_residuals` and computes `tmp_bin_counts` for you.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_residuals = np.ascontiguousarray(pooled_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals' must be an array of np.float64: {error}") from None
    if pooled_residuals.ndim != 1:
        raise ValueError(f"'pooled_residuals' must have 1 dimension, but has {pooled_residuals.ndim}")
    try:
        pooled_residuals_perm = np.ascontiguousarray(pooled_residuals_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals_perm' must be an array of np.int32: {error}") from None
    if pooled_residuals_perm.ndim != 1:
        raise ValueError(f"'pooled_residuals_perm' must have 1 dimension, but has {pooled_residuals_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_residuals = pooled_residuals.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if pooled_residuals_perm.shape[0] != n_residuals:
        raise ValueError(f"'pooled_residuals_perm' has {pooled_residuals_perm.shape[0]} along axis 0, but "
            f"'pooled_residuals' implies n_residuals == {n_residuals}"
        )

    # outputs and work arrays, which the caller never sees
    selected_n_bins = ctypes.c_int(0)
    occupancy_failed = ctypes.c_bool(0)
    shared_residual_range_low = ctypes.c_double(0)
    shared_residual_range_high = ctypes.c_double(0)
    n_pooled_residuals = ctypes.c_int(0)
    min_bin_occupancy = ctypes.c_int(0)
    mean_bin_occupancy = ctypes.c_double(0)
    max_bin_occupancy = ctypes.c_int(0)
    sturges_bins = ctypes.c_int(0)
    fd_bins = ctypes.c_int(0)
    tmp_bin_counts = np.empty((256,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.determine_bin_count_occupancy_expert_c(
        pooled_residuals,
        pooled_residuals_perm,
        ctypes.byref(ctypes.c_int(n_residuals)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(selected_n_bins),
        ctypes.byref(occupancy_failed),
        ctypes.byref(shared_residual_range_low),
        ctypes.byref(shared_residual_range_high),
        ctypes.byref(n_pooled_residuals),
        ctypes.byref(min_bin_occupancy),
        ctypes.byref(mean_bin_occupancy),
        ctypes.byref(max_bin_occupancy),
        ctypes.byref(sturges_bins),
        ctypes.byref(fd_bins),
        tmp_bin_counts,
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_double(gamma_occupancy)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_BIN_COUNT_OCCUPANCY_EXPERT_ARGUMENTS, _DETERMINE_BIN_COUNT_OCCUPANCY_EXPERT_ARGUMENT_SOURCES)

    return {
        "selected_n_bins": selected_n_bins.value,
        "occupancy_failed": occupancy_failed.value,
        "shared_residual_range_low": shared_residual_range_low.value,
        "shared_residual_range_high": shared_residual_range_high.value,
        "n_pooled_residuals": n_pooled_residuals.value,
        "min_bin_occupancy": min_bin_occupancy.value,
        "mean_bin_occupancy": mean_bin_occupancy.value,
        "max_bin_occupancy": max_bin_occupancy.value,
        "sturges_bins": sturges_bins.value,
        "fd_bins": fd_bins.value,
    }

def determine_bin_count_occupancy_exhaustive(
        pooled_residuals,
        n_pooled_residuals,
        shared_residual_range_low,
        shared_residual_range_high,
        m_min=3,
        m_max=120,
        min_residuals_per_bin=10,
):
    r"""Exhaustive brute-force reference implementation of Issue #187's occupancy search

    Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    pooled histogram satisfies the occupancy criterion. A reference for validating
    :func:`tensor_omics.determine_bin_count_occupancy`,
    whose faster search could in principle miss a larger admissible `M`. Takes the range and
    `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.

    Parameters
    ----------
    pooled_residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals for one neighborhood, across all its neighbors and all studies
        NaN is permitted for this value.
    n_pooled_residuals : int
        Count of non-NaN pooled residuals (N_j), from a prior
        :func:`tensor_omics.determine_bin_count_occupancy`
        call
        The minimum valid value is `0`.
        The maximum valid value is `n_residuals`.
    shared_residual_range_low : float
        Lower bound of the histogram range (R_low), from a prior
        :func:`tensor_omics.determine_bin_count_occupancy`
        call
    shared_residual_range_high : float
        Upper bound of the histogram range (R_high), from the same prior call as
        `shared_residual_range_low`
        The minimum valid value is `shared_residual_range_low`.
    m_min : int, optional, default 3
        Smallest candidate bin count tested (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count tested (M_max); raised to `m_min` if smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a candidate bin count to
        be admissible (n_min)
        The minimum valid value is `0`.
        The default value is `10`.

    Returns
    -------
    dict
        with keys:

        selected_n_bins : int
            The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
            above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed
        occupancy_failed : bool
            `True` iff no candidate bin count in [m_min, m_max] satisfies the occupancy criterion
        min_bin_occupancy : int
            Minimum bin count at selected_n_bins; 0 when occupancy_failed
        mean_bin_occupancy : float
            Mean bin count at selected_n_bins; 0 when occupancy_failed
        max_bin_occupancy : int
            Maximum bin count at selected_n_bins; 0 when occupancy_failed

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy_exhaustive`, whose argument names are
    the ones an error message reports.

    This entry point seeds `pooled_residuals_perm` and sorts it by `pooled_residuals`.
    Call `determine_bin_count_occupancy_exhaustive_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_residuals = np.ascontiguousarray(pooled_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals' must be an array of np.float64: {error}") from None
    if pooled_residuals.ndim != 1:
        raise ValueError(f"'pooled_residuals' must have 1 dimension, but has {pooled_residuals.ndim}")

    # what the inputs already say, rather than asking for it again
    n_residuals = pooled_residuals.shape[0]

    # outputs and work arrays, which the caller never sees
    selected_n_bins = ctypes.c_int(0)
    occupancy_failed = ctypes.c_bool(0)
    min_bin_occupancy = ctypes.c_int(0)
    mean_bin_occupancy = ctypes.c_double(0)
    max_bin_occupancy = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.determine_bin_count_occupancy_exhaustive_c(
        pooled_residuals,
        ctypes.byref(ctypes.c_int(n_residuals)),
        ctypes.byref(ctypes.c_int(n_pooled_residuals)),
        ctypes.byref(ctypes.c_double(shared_residual_range_low)),
        ctypes.byref(ctypes.c_double(shared_residual_range_high)),
        ctypes.byref(selected_n_bins),
        ctypes.byref(occupancy_failed),
        ctypes.byref(min_bin_occupancy),
        ctypes.byref(mean_bin_occupancy),
        ctypes.byref(max_bin_occupancy),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_ARGUMENTS, _DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_ARGUMENT_SOURCES)

    return {
        "selected_n_bins": selected_n_bins.value,
        "occupancy_failed": occupancy_failed.value,
        "min_bin_occupancy": min_bin_occupancy.value,
        "mean_bin_occupancy": mean_bin_occupancy.value,
        "max_bin_occupancy": max_bin_occupancy.value,
    }

def determine_bin_count_occupancy_exhaustive_expert(
        pooled_residuals,
        pooled_residuals_perm,
        n_pooled_residuals,
        shared_residual_range_low,
        shared_residual_range_high,
        m_min=3,
        m_max=120,
        min_residuals_per_bin=10,
):
    r"""Exhaustive brute-force reference implementation of Issue #187's occupancy search

    Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    pooled histogram satisfies the occupancy criterion. A reference for validating
    :func:`tensor_omics.determine_bin_count_occupancy`,
    whose faster search could in principle miss a larger admissible `M`. Takes the range and
    `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.

    Parameters
    ----------
    pooled_residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals for one neighborhood, across all its neighbors and all studies
        NaN is permitted for this value.
    pooled_residuals_perm : np.ndarray[np.int32] of shape (n_residuals,)
        Sorting permutation for `pooled_residuals`, ascending, NaN last
        The minimum valid value is `1`.
        The maximum valid value is `n_residuals`.
    n_pooled_residuals : int
        Count of non-NaN pooled residuals (N_j), from a prior
        :func:`tensor_omics.determine_bin_count_occupancy`
        call
        The minimum valid value is `0`.
        The maximum valid value is `n_residuals`.
    shared_residual_range_low : float
        Lower bound of the histogram range (R_low), from a prior
        :func:`tensor_omics.determine_bin_count_occupancy`
        call
    shared_residual_range_high : float
        Upper bound of the histogram range (R_high), from the same prior call as
        `shared_residual_range_low`
        The minimum valid value is `shared_residual_range_low`.
    m_min : int, optional, default 3
        Smallest candidate bin count tested (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count tested (M_max); raised to `m_min` if smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a candidate bin count to
        be admissible (n_min)
        The minimum valid value is `0`.
        The default value is `10`.

    Returns
    -------
    dict
        with keys:

        selected_n_bins : int
            The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
            above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed
        occupancy_failed : bool
            `True` iff no candidate bin count in [m_min, m_max] satisfies the occupancy criterion
        min_bin_occupancy : int
            Minimum bin count at selected_n_bins; 0 when occupancy_failed
        mean_bin_occupancy : float
            Mean bin count at selected_n_bins; 0 when occupancy_failed
        max_bin_occupancy : int
            Maximum bin count at selected_n_bins; 0 when occupancy_failed

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::determine_bin_count_occupancy_exhaustive_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `pooled_residuals_perm` yourself.
    `determine_bin_count_occupancy_exhaustive` seeds `pooled_residuals_perm` and sorts it by `pooled_residuals`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_residuals = np.ascontiguousarray(pooled_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals' must be an array of np.float64: {error}") from None
    if pooled_residuals.ndim != 1:
        raise ValueError(f"'pooled_residuals' must have 1 dimension, but has {pooled_residuals.ndim}")
    try:
        pooled_residuals_perm = np.ascontiguousarray(pooled_residuals_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals_perm' must be an array of np.int32: {error}") from None
    if pooled_residuals_perm.ndim != 1:
        raise ValueError(f"'pooled_residuals_perm' must have 1 dimension, but has {pooled_residuals_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_residuals = pooled_residuals.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if pooled_residuals_perm.shape[0] != n_residuals:
        raise ValueError(f"'pooled_residuals_perm' has {pooled_residuals_perm.shape[0]} along axis 0, but "
            f"'pooled_residuals' implies n_residuals == {n_residuals}"
        )

    # outputs and work arrays, which the caller never sees
    selected_n_bins = ctypes.c_int(0)
    occupancy_failed = ctypes.c_bool(0)
    min_bin_occupancy = ctypes.c_int(0)
    mean_bin_occupancy = ctypes.c_double(0)
    max_bin_occupancy = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.determine_bin_count_occupancy_exhaustive_expert_c(
        pooled_residuals,
        pooled_residuals_perm,
        ctypes.byref(ctypes.c_int(n_residuals)),
        ctypes.byref(ctypes.c_int(n_pooled_residuals)),
        ctypes.byref(ctypes.c_double(shared_residual_range_low)),
        ctypes.byref(ctypes.c_double(shared_residual_range_high)),
        ctypes.byref(selected_n_bins),
        ctypes.byref(occupancy_failed),
        ctypes.byref(min_bin_occupancy),
        ctypes.byref(mean_bin_occupancy),
        ctypes.byref(max_bin_occupancy),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_EXPERT_ARGUMENTS, _DETERMINE_BIN_COUNT_OCCUPANCY_EXHAUSTIVE_EXPERT_ARGUMENT_SOURCES)

    return {
        "selected_n_bins": selected_n_bins.value,
        "occupancy_failed": occupancy_failed.value,
        "min_bin_occupancy": min_bin_occupancy.value,
        "mean_bin_occupancy": mean_bin_occupancy.value,
        "max_bin_occupancy": max_bin_occupancy.value,
    }

def determine_bin_count_stone(
        pooled_residuals,
        residual_to_study,
        n_studies,
        max_n_reps_all_studies,
        n_neighbors,
        bin_count_mode='bin_count_stone_equal',
        m_min=3,
        m_max=120,
        min_residuals_per_bin=10,
        lower_residual_range_quantile=0.05,
        upper_residual_range_quantile=0.95,
):
    r"""Choose one neighborhood's histogram bin count by Stone's cross-validation rule

    Chooses the number of equal-width histogram bins for one neighborhood's pooled residuals by
    Stone's leave-one-out cross-validation (Duke STA 293.3 lecture notes, ch. 6, Theorem 2): for
    `n` values over `m` bins of width `h`, the bin count minimising
    `J(m) = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2` (`p_j` the share of values in bin `j`) is
    chosen. An alternative to
    :func:`tensor_omics.determine_bin_count_occupancy`,
    with the same range, floor and diagnostics.

    The range is the `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of
    all non-NaN pooled residuals, so `h = (R_high - R_low)/m`. The score counts only the residuals
    **inside** the range (`n` is their number). Every `m` in `[m_min, m_max]` is tested, ties
    going to the smaller `m`. A single score (pooled mode, or one study) is bit-exact for up to
    about 200,000 residuals, so ties are exact; a sum over studies may differ in the last bit.

    `residual_to_study` says which study each residual comes from, and `bin_count_mode` how the
    studies are combined into one shared bin count:
    ``'bin_count_stone_equal'`` (the
    default) minimises the sum of every study's own score `J_k(m)`,
    ``'bin_count_stone_size_weighted'``
    the sum of `n_k * J_k(m)` with `n_k` the study's inside count, and
    ``'bin_count_stone_pooled'`` the
    score `J(m)` of all studies' residuals together. The two per-study modes leave out a study
    with fewer than 2 residuals inside the range.

    `min_residuals_per_bin` stays a hard floor, with the occupancy search's definition: a bin
    count is admissible only if the pooled histogram with the outside residuals clamped into the
    edge bins has every bin at or above it (`0` turns the floor off). The minimum is taken over
    the admissible bin counts only.

    `stone_failed = True` means no bin count could be chosen: every pooled residual is NaN, the
    range has width 0 or less, no study has 2 residuals inside the range (pooled mode: fewer than
    2 inside in total), or no bin count passes the floor. `selected_n_bins` is then `m_min`, only
    so a caller ignoring the flag has a value to build with, and the occupancy diagnostics are 0.
    Otherwise the diagnostics describe the pooled clamped histogram at `selected_n_bins`.

    Parameters
    ----------
    pooled_residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals for one neighborhood, across all its neighbors and all studies
        NaN is permitted for this value.
    residual_to_study : np.ndarray[np.int32] of shape (n_residuals,)
        Study each pooled residual comes from
        The minimum valid value is `1`.
        The maximum valid value is `n_studies`.
    n_studies : int
        Number of studies
        The minimum valid value is `1`.
    max_n_reps_all_studies : int
        Maximum number of replicates across all studies
        The minimum valid value is `1`.
    n_neighbors : int
        Neighborhood size of the candidate under test
        The minimum valid value is `1`.
    bin_count_mode : str, one of 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_stone_equal'
        How the studies' scores are combined into one bin count

        The default value is `'bin_count_stone_equal'`.
    m_min : int, optional, default 3
        Smallest candidate bin count tested (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count tested (M_max); raised to `m_min` if smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin of the pooled clamped histogram must
        reach for a candidate bin count to be admissible (0 turns the floor off)
        The minimum valid value is `0`.
        The default value is `10`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for this neighborhood's own lower residual-range bound
        (shared_residual_range_low)
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for this neighborhood's own upper residual-range bound
        (shared_residual_range_high)
        The minimum valid value is `lower_residual_range_quantile`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    dict
        with keys:

        selected_n_bins : int
            The chosen bin count: the admissible bin count in [m_min, m_max] with the smallest
            Stone score; m_min when stone_failed
        stone_failed : bool
            `True` iff no bin count could be chosen (see above); the caller should then reject
            this neighborhood rather than build a histogram from selected_n_bins
        shared_residual_range_low : float
            This neighborhood's lower residual-range bound (R_low); 0.0 if every residual is NaN
        shared_residual_range_high : float
            This neighborhood's upper residual-range bound (R_high); 0.0 if every residual is NaN
        n_pooled_residuals : int
            Count of non-NaN pooled residuals (N_j)
        min_bin_occupancy : int
            Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        mean_bin_occupancy : float
            Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            stone_failed
        max_bin_occupancy : int
            Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        sturges_bins : int
            Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        fd_bins : int
            Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::determine_bin_count_stone`, whose argument names are
    the ones an error message reports.

    This entry point seeds `pooled_residuals_perm` and sorts it by `pooled_residuals` and computes `tmp_bin_counts`, `tmp_residual_to_study_perm` and `tmp_study_block_sizes` for you.
    Call `determine_bin_count_stone_expert` to do that yourself.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_residuals = np.ascontiguousarray(pooled_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals' must be an array of np.float64: {error}") from None
    if pooled_residuals.ndim != 1:
        raise ValueError(f"'pooled_residuals' must have 1 dimension, but has {pooled_residuals.ndim}")
    try:
        residual_to_study = np.ascontiguousarray(residual_to_study, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residual_to_study' must be an array of np.int32: {error}") from None
    if residual_to_study.ndim != 1:
        raise ValueError(f"'residual_to_study' must have 1 dimension, but has {residual_to_study.ndim}")
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

    # what the inputs already say, rather than asking for it again
    n_residuals = pooled_residuals.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if residual_to_study.shape[0] != n_residuals:
        raise ValueError(f"'residual_to_study' has {residual_to_study.shape[0]} along axis 0, but "
            f"'pooled_residuals' implies n_residuals == {n_residuals}"
        )

    # outputs and work arrays, which the caller never sees
    selected_n_bins = ctypes.c_int(0)
    stone_failed = ctypes.c_bool(0)
    shared_residual_range_low = ctypes.c_double(0)
    shared_residual_range_high = ctypes.c_double(0)
    n_pooled_residuals = ctypes.c_int(0)
    min_bin_occupancy = ctypes.c_int(0)
    mean_bin_occupancy = ctypes.c_double(0)
    max_bin_occupancy = ctypes.c_int(0)
    sturges_bins = ctypes.c_int(0)
    fd_bins = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.determine_bin_count_stone_c(
        pooled_residuals,
        ctypes.byref(ctypes.c_int(n_residuals)),
        residual_to_study,
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(selected_n_bins),
        ctypes.byref(stone_failed),
        ctypes.byref(shared_residual_range_low),
        ctypes.byref(shared_residual_range_high),
        ctypes.byref(n_pooled_residuals),
        ctypes.byref(min_bin_occupancy),
        ctypes.byref(mean_bin_occupancy),
        ctypes.byref(max_bin_occupancy),
        ctypes.byref(sturges_bins),
        ctypes.byref(fd_bins),
        bin_count_mode,
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_BIN_COUNT_STONE_ARGUMENTS, _DETERMINE_BIN_COUNT_STONE_ARGUMENT_SOURCES)

    return {
        "selected_n_bins": selected_n_bins.value,
        "stone_failed": stone_failed.value,
        "shared_residual_range_low": shared_residual_range_low.value,
        "shared_residual_range_high": shared_residual_range_high.value,
        "n_pooled_residuals": n_pooled_residuals.value,
        "min_bin_occupancy": min_bin_occupancy.value,
        "mean_bin_occupancy": mean_bin_occupancy.value,
        "max_bin_occupancy": max_bin_occupancy.value,
        "sturges_bins": sturges_bins.value,
        "fd_bins": fd_bins.value,
    }

def determine_bin_count_stone_expert(
        pooled_residuals,
        pooled_residuals_perm,
        residual_to_study,
        n_studies,
        max_n_reps_all_studies,
        n_neighbors,
        bin_count_mode='bin_count_stone_equal',
        m_min=3,
        m_max=120,
        min_residuals_per_bin=10,
        lower_residual_range_quantile=0.05,
        upper_residual_range_quantile=0.95,
):
    r"""Choose one neighborhood's histogram bin count by Stone's cross-validation rule

    Chooses the number of equal-width histogram bins for one neighborhood's pooled residuals by
    Stone's leave-one-out cross-validation (Duke STA 293.3 lecture notes, ch. 6, Theorem 2): for
    `n` values over `m` bins of width `h`, the bin count minimising
    `J(m) = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2` (`p_j` the share of values in bin `j`) is
    chosen. An alternative to
    :func:`tensor_omics.determine_bin_count_occupancy`,
    with the same range, floor and diagnostics.

    The range is the `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of
    all non-NaN pooled residuals, so `h = (R_high - R_low)/m`. The score counts only the residuals
    **inside** the range (`n` is their number). Every `m` in `[m_min, m_max]` is tested, ties
    going to the smaller `m`. A single score (pooled mode, or one study) is bit-exact for up to
    about 200,000 residuals, so ties are exact; a sum over studies may differ in the last bit.

    `residual_to_study` says which study each residual comes from, and `bin_count_mode` how the
    studies are combined into one shared bin count:
    ``'bin_count_stone_equal'`` (the
    default) minimises the sum of every study's own score `J_k(m)`,
    ``'bin_count_stone_size_weighted'``
    the sum of `n_k * J_k(m)` with `n_k` the study's inside count, and
    ``'bin_count_stone_pooled'`` the
    score `J(m)` of all studies' residuals together. The two per-study modes leave out a study
    with fewer than 2 residuals inside the range.

    `min_residuals_per_bin` stays a hard floor, with the occupancy search's definition: a bin
    count is admissible only if the pooled histogram with the outside residuals clamped into the
    edge bins has every bin at or above it (`0` turns the floor off). The minimum is taken over
    the admissible bin counts only.

    `stone_failed = True` means no bin count could be chosen: every pooled residual is NaN, the
    range has width 0 or less, no study has 2 residuals inside the range (pooled mode: fewer than
    2 inside in total), or no bin count passes the floor. `selected_n_bins` is then `m_min`, only
    so a caller ignoring the flag has a value to build with, and the occupancy diagnostics are 0.
    Otherwise the diagnostics describe the pooled clamped histogram at `selected_n_bins`.

    Parameters
    ----------
    pooled_residuals : np.ndarray[np.float64] of shape (n_residuals,)
        Pooled signed residuals for one neighborhood, across all its neighbors and all studies
        NaN is permitted for this value.
    pooled_residuals_perm : np.ndarray[np.int32] of shape (n_residuals,)
        Sorting permutation for `pooled_residuals`, ascending, NaN last
        The minimum valid value is `1`.
        The maximum valid value is `n_residuals`.
    residual_to_study : np.ndarray[np.int32] of shape (n_residuals,)
        Study each pooled residual comes from
        The minimum valid value is `1`.
        The maximum valid value is `n_studies`.
    n_studies : int
        Number of studies
        The minimum valid value is `1`.
    max_n_reps_all_studies : int
        Maximum number of replicates across all studies
        The minimum valid value is `1`.
    n_neighbors : int
        Neighborhood size of the candidate under test
        The minimum valid value is `1`.
    bin_count_mode : str, one of 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_stone_equal'
        How the studies' scores are combined into one bin count

        The default value is `'bin_count_stone_equal'`.
    m_min : int, optional, default 3
        Smallest candidate bin count tested (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count tested (M_max); raised to `m_min` if smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin of the pooled clamped histogram must
        reach for a candidate bin count to be admissible (0 turns the floor off)
        The minimum valid value is `0`.
        The default value is `10`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for this neighborhood's own lower residual-range bound
        (shared_residual_range_low)
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for this neighborhood's own upper residual-range bound
        (shared_residual_range_high)
        The minimum valid value is `lower_residual_range_quantile`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    dict
        with keys:

        selected_n_bins : int
            The chosen bin count: the admissible bin count in [m_min, m_max] with the smallest
            Stone score; m_min when stone_failed
        stone_failed : bool
            `True` iff no bin count could be chosen (see above); the caller should then reject
            this neighborhood rather than build a histogram from selected_n_bins
        shared_residual_range_low : float
            This neighborhood's lower residual-range bound (R_low); 0.0 if every residual is NaN
        shared_residual_range_high : float
            This neighborhood's upper residual-range bound (R_high); 0.0 if every residual is NaN
        n_pooled_residuals : int
            Count of non-NaN pooled residuals (N_j)
        min_bin_occupancy : int
            Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        mean_bin_occupancy : float
            Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            stone_failed
        max_bin_occupancy : int
            Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        sturges_bins : int
            Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        fd_bins : int
            Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_binning::determine_bin_count_stone_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `pooled_residuals_perm`, `tmp_bin_counts`, `tmp_residual_to_study_perm` and `tmp_study_block_sizes` yourself.
    `determine_bin_count_stone` seeds `pooled_residuals_perm` and sorts it by `pooled_residuals` and computes `tmp_bin_counts`, `tmp_residual_to_study_perm` and `tmp_study_block_sizes` for you.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_residuals = np.ascontiguousarray(pooled_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals' must be an array of np.float64: {error}") from None
    if pooled_residuals.ndim != 1:
        raise ValueError(f"'pooled_residuals' must have 1 dimension, but has {pooled_residuals.ndim}")
    try:
        pooled_residuals_perm = np.ascontiguousarray(pooled_residuals_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_residuals_perm' must be an array of np.int32: {error}") from None
    if pooled_residuals_perm.ndim != 1:
        raise ValueError(f"'pooled_residuals_perm' must have 1 dimension, but has {pooled_residuals_perm.ndim}")
    try:
        residual_to_study = np.ascontiguousarray(residual_to_study, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residual_to_study' must be an array of np.int32: {error}") from None
    if residual_to_study.ndim != 1:
        raise ValueError(f"'residual_to_study' must have 1 dimension, but has {residual_to_study.ndim}")
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

    # what the inputs already say, rather than asking for it again
    n_residuals = pooled_residuals.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if pooled_residuals_perm.shape[0] != n_residuals:
        raise ValueError(f"'pooled_residuals_perm' has {pooled_residuals_perm.shape[0]} along axis 0, but "
            f"'pooled_residuals' implies n_residuals == {n_residuals}"
        )
    if residual_to_study.shape[0] != n_residuals:
        raise ValueError(f"'residual_to_study' has {residual_to_study.shape[0]} along axis 0, but "
            f"'pooled_residuals' implies n_residuals == {n_residuals}"
        )

    # outputs and work arrays, which the caller never sees
    selected_n_bins = ctypes.c_int(0)
    stone_failed = ctypes.c_bool(0)
    shared_residual_range_low = ctypes.c_double(0)
    shared_residual_range_high = ctypes.c_double(0)
    n_pooled_residuals = ctypes.c_int(0)
    min_bin_occupancy = ctypes.c_int(0)
    mean_bin_occupancy = ctypes.c_double(0)
    max_bin_occupancy = ctypes.c_int(0)
    sturges_bins = ctypes.c_int(0)
    fd_bins = ctypes.c_int(0)
    tmp_bin_counts = np.empty((256,), dtype=np.int32, order='C')
    tmp_residual_to_study_perm = np.empty((n_residuals,), dtype=np.int32, order='C')
    tmp_study_block_sizes = np.empty((n_studies,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.determine_bin_count_stone_expert_c(
        pooled_residuals,
        pooled_residuals_perm,
        ctypes.byref(ctypes.c_int(n_residuals)),
        residual_to_study,
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(selected_n_bins),
        ctypes.byref(stone_failed),
        ctypes.byref(shared_residual_range_low),
        ctypes.byref(shared_residual_range_high),
        ctypes.byref(n_pooled_residuals),
        ctypes.byref(min_bin_occupancy),
        ctypes.byref(mean_bin_occupancy),
        ctypes.byref(max_bin_occupancy),
        ctypes.byref(sturges_bins),
        ctypes.byref(fd_bins),
        tmp_bin_counts,
        tmp_residual_to_study_perm,
        tmp_study_block_sizes,
        bin_count_mode,
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DETERMINE_BIN_COUNT_STONE_EXPERT_ARGUMENTS, _DETERMINE_BIN_COUNT_STONE_EXPERT_ARGUMENT_SOURCES)

    return {
        "selected_n_bins": selected_n_bins.value,
        "stone_failed": stone_failed.value,
        "shared_residual_range_low": shared_residual_range_low.value,
        "shared_residual_range_high": shared_residual_range_high.value,
        "n_pooled_residuals": n_pooled_residuals.value,
        "min_bin_occupancy": min_bin_occupancy.value,
        "mean_bin_occupancy": mean_bin_occupancy.value,
        "max_bin_occupancy": max_bin_occupancy.value,
        "sturges_bins": sturges_bins.value,
        "fd_bins": fd_bins.value,
    }
