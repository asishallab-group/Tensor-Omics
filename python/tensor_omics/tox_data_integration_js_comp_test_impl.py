r"""tox_data_integration_js_comp_test_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search

The data-driven `(n_points, n_neighbors)` parameter-stabilization search this pipeline runs
before the JSD-Comp-Test proper (Issue #126): a GAMMA-decay candidate grid
(:func:`tensor_omics.generate_js_comp_test_candidates`
generates candidate `(n_points, n_neighbors)` pairs only -- each candidate's real
per-neighborhood histogram bin count is decided later, per reference point, by
:func:`tensor_omics.determine_bin_count_occupancy`),
two admissibility gates a candidate must pass before it is bootstrapped
(:func:`tensor_omics.check_neighborhood_overlaps`,
:func:`tensor_omics.check_mean_pmf_min_counts`),
and the plateau check that decides when the search has converged
(:func:`tensor_omics.check_plateau_condition`).
:func:`tensor_omics.determine_bin_count_occupancy_exhaustive`
is a brute-force reference implementation of the same per-point bin-count search, exhaustively
testing every candidate `M` instead of the fast geometric-search-then-refinement the production
routine uses, for validating that the fast search's own result is correct.
`calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
allocates its own. Once a candidate has passed both gates, its bootstrap confidence interval
is resampled from the pooled consensus histogram by
:func:`tensor_omics.bootstrap_histogram` (heap
size recommended by
:func:`tensor_omics.calc_js_comp_test_n_top_k_jsds`).
Later stages of the port add the K-study permutation test and the top-level orchestrator that
wire these building blocks together.

Python binding, generated from tox_data_integration_js_comp_test_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_js_comp_test_candidate_bounds_c.restype = None
_lib.calc_js_comp_test_candidate_bounds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "max_n_points_candidate", "max_n_neighbors_candidate",)

_lib.gather_pooled_neighborhood_residuals_c.restype = None
_lib.gather_pooled_neighborhood_residuals_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENTS = ("residuals", "max_n_reps_all_studies", "max_n_genes_all_studies", "n_neighbors", "n_studies", "neighborhood_indices_point", "pooled_residuals", "ierr",)
#: For a derived argument, the one the caller passed it in
_GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENT_SOURCES = (None, "residuals", "residuals", "neighborhood_indices_point", "residuals", None, None, None,)

_lib.calc_js_comp_test_n_top_k_jsds_c.restype = None
_lib.calc_js_comp_test_n_top_k_jsds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_N_TOP_K_JSDS_ARGUMENTS = ("n_bootstraps", "two_sided_bootstrapping_significance_level", "n_top_k",)

def calc_js_comp_test_candidate_bounds(
        max_n_genes_all_studies,
):
    r"""Recommend upper bounds for the js-comp-test candidate-grid work arrays

    Closed-form from the GAMMA-decay candidate-grid formula
    :func:`tensor_omics.generate_js_comp_test_candidates`
    uses: `max_n_points_candidate` is exactly the grid's first (largest) `n_points` candidate.
    `max_n_neighbors_candidate` is a safe, not necessarily tight, upper bound on the grid's
    largest `n_neighbors` candidate -- reached with the smallest `KX_FACTORS` entry at the
    smallest `n_points_high` the grid loop ever uses, which by construction never drops below
    `n_points_low`.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        max_n_points_candidate : int
            Exact upper bound on the `n_points` candidate the grid ever produces
        max_n_neighbors_candidate : int
            Safe upper bound on the `n_neighbors` candidate the grid ever produces

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_impl::calc_js_comp_test_candidate_bounds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = ctypes.c_int(0)
    max_n_neighbors_candidate = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_js_comp_test_candidate_bounds_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(max_n_points_candidate),
        ctypes.byref(max_n_neighbors_candidate),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "max_n_neighbors_candidate": max_n_neighbors_candidate.value,
    }

def gather_pooled_neighborhood_residuals(
        residuals,
        neighborhood_indices_point,
):
    r"""Pool one reference point's residuals across every neighbor and every study

    Given one reference point's own per-study neighbor gene indices (one column of a larger
    `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)`, as produced by
    :func:`tensor_omics.construct_neighborhoods_ranged`),
    gathers that point's residual values from every neighbor gene, across every study, into one
    flat pooled array. This is the exact same pooling
    :func:`tensor_omics.run_js_comp_test` and
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    perform internally, per reference point, before handing the result to
    :func:`tensor_omics.determine_bin_count_occupancy`'s
    own occupancy search -- published so a caller can reconstruct that exact same input directly
    on real data and feed it to
    :func:`tensor_omics.determine_bin_count_occupancy_exhaustive`
    (or to `determine_bin_count_occupancy` itself), to check whether the fast search and the
    exhaustive reference ever actually disagree in practice, not just on a synthetic fixture.

    Parameters
    ----------
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Matrix of signed residuals per study, NaN explicitly allowed for missing values;
        unvalidated on this export path, so any value passes through as-is
    neighborhood_indices_point : np.ndarray[np.int32] of shape (n_neighbors, n_studies,), column-major (order='F')
        Gene indices of one reference point's neighborhood, per study -- one column of a
        larger neighborhood_indices_all_studies(n_neighbors, n_points, n_studies), as sliced
        by the caller

    Returns
    -------
    pooled_residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies*n_neighbors*n_studies,), read-only
        The pooled residual values for this reference point, across every neighbor and every
        study, laid out exactly as a (max_n_reps_all_studies, n_neighbors, n_studies) array
        would be
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_impl::gather_pooled_neighborhood_residuals`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        residuals = np.asfortranarray(residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residuals' must be an array of np.float64: {error}") from None
    if residuals.ndim != 3:
        raise ValueError(f"'residuals' must have 3 dimensions, but has {residuals.ndim}")
    try:
        neighborhood_indices_point = np.asfortranarray(neighborhood_indices_point, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_indices_point' must be an array of np.int32: {error}") from None
    if neighborhood_indices_point.ndim != 2:
        raise ValueError(f"'neighborhood_indices_point' must have 2 dimensions, but has {neighborhood_indices_point.ndim}")

    # what the inputs already say, rather than asking for it again
    max_n_reps_all_studies = residuals.shape[0]
    max_n_genes_all_studies = residuals.shape[1]
    n_neighbors = neighborhood_indices_point.shape[0]
    n_studies = residuals.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_indices_point.shape[1] != n_studies:
        raise ValueError(f"'neighborhood_indices_point' has {neighborhood_indices_point.shape[1]} along axis 1, but "
            f"'residuals' implies n_studies == {n_studies}"
        )

    # outputs and work arrays, which the caller never sees
    pooled_residuals = np.empty((max_n_reps_all_studies*n_neighbors*n_studies,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.gather_pooled_neighborhood_residuals_c(
        residuals,
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_studies)),
        neighborhood_indices_point,
        pooled_residuals,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENTS, _GATHER_POOLED_NEIGHBORHOOD_RESIDUALS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    pooled_residuals.flags.writeable = False

    return pooled_residuals

def calc_js_comp_test_n_top_k_jsds(
        n_bootstraps,
        two_sided_bootstrapping_significance_level=2.5,
):
    r"""Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval

    Ported from 125-stabilize-jscomp's inline `n_bootstrapping_top_k_jsds` computation in
    `determine_js_comp_test_n_points_n_neighbors_alloc`: sizes the top-k/bottom-k heaps
    :func:`tensor_omics.bootstrap_histogram` uses
    to track a two-sided bootstrap confidence interval's endpoints. E.g. for a 95% two-sided
    interval (2.5% reserved at each end, the default) with `n_bootstraps=1000`,
    `n_top_k = max(1, floor(0.025*1000)) = 25`.

    Parameters
    ----------
    n_bootstraps : int
        Number of bootstrap resamples that will be performed
        The minimum valid value is `1`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Two-sided significance level, as a percentage reserved at each end of the bootstrap
        distribution (e.g. 2.5 for a 95% two-sided interval)
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.

    Returns
    -------
    n_top_k : int
        Number of elements to keep at each end of the bootstrap distribution, at least 1

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_impl::calc_js_comp_test_n_top_k_jsds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    n_top_k = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_js_comp_test_n_top_k_jsds_c(
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        ctypes.byref(ctypes.c_double(two_sided_bootstrapping_significance_level)),
        ctypes.byref(n_top_k),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_JS_COMP_TEST_N_TOP_K_JSDS_ARGUMENTS)

    return n_top_k.value
