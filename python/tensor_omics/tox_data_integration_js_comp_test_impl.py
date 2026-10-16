r"""tox_data_integration_js_comp_test_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search

The global JS-Comp-Test of K studies, with a data-driven choice of its neighborhood parameters.

Fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each):
:func:`tensor_omics.run_js_comp_test_parameter_search`
searches the `(n_points, n_neighbors)` candidates of
:func:`tensor_omics.generate_js_comp_test_candidates`
for a bootstrap plateau, and
:func:`tensor_omics.run_js_comp_test` runs the
final test (histograms, JSDs, weighting, permutation test) for the chosen pair.

Adaptive neighborhoods, grown per reference point while the residual dispersion stays stable:
:func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`
searches the `(k_start, k_step, k_max)` candidates of
:func:`tensor_omics.generate_adaptive_js_comp_test_candidates`,
:func:`tensor_omics.construct_adaptive_neighborhoods`
builds the neighborhoods for the chosen knobs, and
:func:`tensor_omics.run_js_comp_test_adaptive` runs
the final test on them.

Each reference point's histogram bin count is chosen by the occupancy search
:func:`tensor_omics.determine_bin_count_occupancy`
(default) or by Stone's cross-validation
:func:`tensor_omics.determine_bin_count_stone`
(`bin_count_mode`). The building blocks are published too: the admissibility gates
:func:`tensor_omics.check_neighborhood_overlaps` and
:func:`tensor_omics.check_mean_pmf_min_counts`, the
plateau checks
:func:`tensor_omics.check_plateau_condition` and
:func:`tensor_omics.check_effect_size_plateau_condition`,
the bootstrap
:func:`tensor_omics.bootstrap_histogram`, and the
brute-force reference
:func:`tensor_omics.determine_bin_count_occupancy_exhaustive`.
`calc_js_comp_test_candidate_bounds`, `calc_adaptive_js_comp_test_bounds` and
`calc_js_comp_test_n_top_k_jsds` size the work arrays for a caller that allocates its own.

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
_CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "max_n_points_candidate", "max_n_neighbors_candidate", "ierr",)

_lib.calc_adaptive_js_comp_test_bounds_c.restype = None
_lib.calc_adaptive_js_comp_test_bounds_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_ADAPTIVE_JS_COMP_TEST_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "n_studies", "max_n_points_candidate", "ierr",)

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

    For the candidate grid of
    :func:`tensor_omics.generate_js_comp_test_candidates`:
    `max_n_points_candidate` is exactly its largest `n_points` candidate, and
    `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`.

    Rejects a non-positive `max_n_genes_all_studies` with invalid input; both outputs are then
    left undefined.

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

def calc_adaptive_js_comp_test_bounds(
        max_n_genes_all_studies,
        n_studies,
):
    r"""Recommend a reference-point capacity for the adaptive js-comp-test candidate sequence

    Sizes the per-point work arrays of the adaptive parameter search over the candidates of
    :func:`tensor_omics.generate_adaptive_js_comp_test_candidates`.
    A candidate's reference-point count is only known after its neighborhoods are grown, so this
    is a **practical capacity, not a proven bound**:
    `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
    ``ADAPTIVE_POINT_CAPACITY_FACTOR``,
    `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the smallest
    candidate `k_start`. The search rejects a candidate that exceeds it with a capacity status.

    Rejects, with invalid input naming neither argument when the bound is joint, a
    `max_n_genes_all_studies` or `n_studies` below 1, a product `N` beyond the 32-bit range, or a
    first candidate `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))` above
    `huge(1)/4 = 536870911` (so that its `k_max = 4*k_start_1` fits); the output is then
    left undefined.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    n_studies : int
        Number of studies
        The minimum valid value is `1`.

    Returns
    -------
    max_n_points_candidate : int
        Practical upper bound on the number of reference points any adaptive candidate is
        allowed to emerge with

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_impl::calc_adaptive_js_comp_test_bounds`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_adaptive_js_comp_test_bounds_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(max_n_points_candidate),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_ADAPTIVE_JS_COMP_TEST_BOUNDS_ARGUMENTS)

    return max_n_points_candidate.value

def gather_pooled_neighborhood_residuals(
        residuals,
        neighborhood_indices_point,
):
    r"""Pool one reference point's residuals across every neighbor and every study

    Gathers one reference point's residuals from every neighbor gene of every study into one flat
    array, given that point's per-study neighbor gene indices (one column of the
    `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)` of
    :func:`tensor_omics.construct_neighborhoods_ranged`).
    Study `s` fills pool block `(s - 1)*max_n_reps_all_studies*n_neighbors + 1` to
    `s*max_n_reps_all_studies*n_neighbors`. This is the pooling
    :func:`tensor_omics.run_js_comp_test` and
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    feed to the bin-count search, so a caller can pass the same input to
    :func:`tensor_omics.determine_bin_count_occupancy_exhaustive`
    on real data.

    Parameters
    ----------
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Matrix of signed residuals per study, NaN allowed for missing values; not validated
    neighborhood_indices_point : np.ndarray[np.int32] of shape (n_neighbors, n_studies,), column-major (order='F')
        Gene indices of one reference point's neighborhood, per study; each in
        [1, max_n_genes_all_studies]

    Returns
    -------
    pooled_residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies*n_neighbors*n_studies,), read-only
        This reference point's pooled residuals, laid out as a
        (max_n_reps_all_studies, n_neighbors, n_studies) array
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

    Sizes the top-k/bottom-k heaps
    :func:`tensor_omics.bootstrap_histogram` uses
    for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
    the default 2.5 and `n_bootstraps=1000`.

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
