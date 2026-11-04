r"""tox_data_integration_js_comp_test_fixed_k_search

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Fixed-k Parameter Search

The data-driven choice of the JS-Comp-Test's fixed-k neighborhoods, `n_points` reference
points of `n_neighbors` genes each.
:func:`tensor_omics.run_js_comp_test_parameter_search`
walks the `(n_points, n_neighbors)` candidate grid of
:func:`tensor_omics.generate_js_comp_test_candidates`
from finest to coarsest and stops at the first admissible candidate whose bootstrapped JSDs
have reached a plateau, by the criteria of
:func:`tensor_omics.check_plateau_condition` and
:func:`tensor_omics.check_effect_size_plateau_condition`.
By default the grid holds the large neighborhoods
(``'neighborhood_size_large'``);
`neighborhood_size_mode` selects the small ones instead
(``'neighborhood_size_small'``).
`grid_points_factor` and `grid_min_points` set the grid's first `n_points` and its floor.
`calc_js_comp_test_candidate_bounds` sizes the work arrays for a caller that allocates its own.

The recommended way to choose the neighborhoods is the adaptive search,
:func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`,
whose neighborhoods follow the local density of the gene means; this fixed-k search stays
available.

Python binding, generated from tox_data_integration_js_comp_test_fixed_k_search. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.calc_js_comp_test_candidate_bounds_c.restype = None
_lib.calc_js_comp_test_candidate_bounds_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS = ("max_n_genes_all_studies", "max_n_points_candidate", "max_n_neighbors_candidate", "grid_points_factor", "grid_min_points", "ierr",)

_lib.generate_js_comp_test_candidates_c.restype = None
_lib.generate_js_comp_test_candidates_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GENERATE_JS_COMP_TEST_CANDIDATES_ARGUMENTS = ("max_n_genes_all_studies", "candidates_n_points_n_neighbors", "n_candidates", "neighborhood_size_mode", "grid_points_factor", "grid_min_points", "ierr",)

_lib.run_js_comp_test_parameter_search_c.restype = None
_lib.run_js_comp_test_parameter_search_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "max_n_neighbors_candidate", "n_points", "n_neighbors", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_pair_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_ci_overlap", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "neighborhood_size_mode", "delta_abs_tolerance", "divergence_mode", "grid_points_factor", "grid_min_points", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", None, None, None, None, "n_bins_per_point", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def calc_js_comp_test_candidate_bounds(
        max_n_genes_all_studies,
        grid_points_factor=4.0,
        grid_min_points=200,
):
    r"""Recommend upper bounds for the js-comp-test candidate-grid work arrays

    For the candidate grid of
    :func:`tensor_omics.generate_js_comp_test_candidates`
    with the same `grid_points_factor` and `grid_min_points`: `max_n_points_candidate` is exactly
    its largest `n_points` candidate, and `max_n_neighbors_candidate` a safe (not necessarily
    tight) bound on its largest `n_neighbors`, in every `neighborhood_size_mode` (the bound uses
    the large-neighborhood factor at the grid's floor). Pass the same two options as the
    parameter search that is sized with the result.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    grid_points_factor : float, optional, default 4.0
        The grid's `grid_points_factor`, as passed to the parameter search
        The minimum valid value is `above(0.0)`.
        The default value is `4.0`.
    grid_min_points : int, optional, default 200
        The grid's `grid_min_points`, as passed to the parameter search
        The minimum valid value is `1`.
        The maximum valid value is `MAX_POINTS`.
        The default value is `200`.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_fixed_k_search::calc_js_comp_test_candidate_bounds`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    max_n_genes_all_studies = _contract.integer(max_n_genes_all_studies, "max_n_genes_all_studies")
    grid_points_factor = _contract.real(grid_points_factor, "grid_points_factor")
    grid_min_points = _contract.integer(grid_min_points, "grid_min_points")

    # outputs and work arrays, which the caller never sees
    max_n_points_candidate = _ctypes.c_int(0)
    max_n_neighbors_candidate = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.calc_js_comp_test_candidate_bounds_c(
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(max_n_points_candidate),
        _ctypes.byref(max_n_neighbors_candidate),
        _ctypes.byref(_ctypes.c_double(grid_points_factor)),
        _ctypes.byref(_ctypes.c_int(grid_min_points)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_JS_COMP_TEST_CANDIDATE_BOUNDS_ARGUMENTS)

    return {
        "max_n_points_candidate": max_n_points_candidate.value,
        "max_n_neighbors_candidate": max_n_neighbors_candidate.value,
    }

def generate_js_comp_test_candidates(
        max_n_genes_all_studies,
        neighborhood_size_mode='neighborhood_size_large',
        grid_points_factor=4.0,
        grid_min_points=200,
):
    r"""Generate the GAMMA-decay (n_points, n_neighbors) candidate grid

    Starting from `n_points_high = clamp(ceiling(grid_points_factor*sqrt(G)), grid_min_points, MAX_POINTS)`
    (`G` = `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
    `n_points_low = max(grid_min_points, ceiling(N_POINTS_LOW_FRACTION*n_points_high))`, and pairs each
    distinct `n_points` with `n_neighbors = min(G, max(1, floor(G/(f*n_points_high))))` for the
    factor `f` of `KX_FACTORS` that `neighborhood_size_mode` selects: by default `f = 0.25` (the
    large neighborhoods), so `n_neighbors` grows smoothly from one candidate to the next. The cap
    at `G` only binds for a `grid_min_points` below 4 (below 2 for the small neighborhoods), where
    `G/(f*n_points_high)` can exceed the number of genes. Small
    neighborhoods are not the default: a higher JSCompTest value from small neighborhoods only
    shows that they are too small.

    When the next decayed value would fall below `n_points_low`, the floor itself is added as the
    last, coarsest candidate (with `n_neighbors` from `n_points_low`), unless the coarsest
    candidate already is the floor or the grid already holds MAX_CANDIDATE_PAIRS candidates.
    Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
    candidate. Bin counts are not part of a candidate; they are chosen per reference point later.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.
    neighborhood_size_mode : str, one of 'neighborhood_size_large' | 'neighborhood_size_small', optional, default 'neighborhood_size_large'
        Which neighborhood size the grid pairs with each `n_points`

        The default value is `'neighborhood_size_large'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    grid_points_factor : float, optional, default 4.0
        Factor in the grid's first reference-point count `ceiling(grid_points_factor * sqrt(G))`
        The minimum valid value is `above(0.0)`.
        The default value is `4.0`.
    grid_min_points : int, optional, default 200
        The grid's floor: the lower clamp of its first `n_points`, the least `n_points` it
        keeps, and its last, coarsest candidate
        The minimum valid value is `1`.
        The maximum valid value is `MAX_POINTS`.
        The default value is `200`.

    Returns
    -------
    candidates_n_points_n_neighbors : np.ndarray[np.int32] of shape (2, 16,), column-major (order='F'), read-only
        Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
        The first `n_candidates` elements will hold the results.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_fixed_k_search::generate_js_comp_test_candidates`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    max_n_genes_all_studies = _contract.integer(max_n_genes_all_studies, "max_n_genes_all_studies")
    neighborhood_size_mode = _np.array([_contract.mode(neighborhood_size_mode, "neighborhood_size_mode", ("neighborhood_size_large", "neighborhood_size_small",)).ljust(23)], dtype="S23")
    grid_points_factor = _contract.real(grid_points_factor, "grid_points_factor")
    grid_min_points = _contract.integer(grid_min_points, "grid_min_points")

    # outputs and work arrays, which the caller never sees
    candidates_n_points_n_neighbors = _np.empty((2, 16,), dtype=_np.int32, order='F')
    n_candidates = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.generate_js_comp_test_candidates_c(
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        candidates_n_points_n_neighbors,
        _ctypes.byref(n_candidates),
        neighborhood_size_mode,
        _ctypes.byref(_ctypes.c_double(grid_points_factor)),
        _ctypes.byref(_ctypes.c_int(grid_min_points)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GENERATE_JS_COMP_TEST_CANDIDATES_ARGUMENTS)

    # a result is a value: modify a copy, not this
    candidates_n_points_n_neighbors.flags.writeable = False

    return candidates_n_points_n_neighbors[..., :n_candidates.value]

def run_js_comp_test_parameter_search(
        gene_means,
        residuals,
        n_bootstraps,
        join_method,
        max_n_points_candidate,
        max_n_neighbors_candidate,
        min_residuals_per_bin=10,
        min_neighbor_overlap=0.1,
        succeeding_ci_overlap=0.9,
        plateau_mode='plateau_effect_size',
        delta_median_threshold=0.05,
        delta_max_threshold=0.1,
        delta_epsilon=1e-10,
        delta_min_consecutive_transitions=2,
        m_min=3,
        m_max=120,
        gamma_occupancy=1.25,
        lower_residual_range_quantile=0.05,
        upper_residual_range_quantile=0.95,
        two_sided_bootstrapping_significance_level=2.5,
        random_seed=42,
        bin_count_mode='bin_count_stone_equal',
        neighborhood_size_mode='neighborhood_size_large',
        delta_abs_tolerance=0.0001,
        divergence_mode='divergence_pairwise_k2',
        grid_points_factor=4.0,
        grid_min_points=200,
):
    r"""Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting

    Walks the candidates of
    :func:`tensor_omics.generate_js_comp_test_candidates`
    from finest to coarsest resolution, by default the large neighborhoods
    (`neighborhood_size_mode`; `grid_points_factor` and `grid_min_points` shape the grid). For
    each candidate:

    1. Reference points from the pooled gene means
    (:func:`tensor_omics.pool_means`) and every
    study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    (:func:`tensor_omics.check_neighborhood_overlaps`).
    2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    residuals, then every study's histograms.
    3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    `min_residuals_per_bin`
    (:func:`tensor_omics.check_mean_pmf_min_counts`).
    A candidate that passes both gates is *admissible*.
    4. The observed global JSD per study, its bootstrap confidence interval
    (:func:`tensor_omics.bootstrap_histogram`),
    and both plateau criteria: CI overlap with the best candidate so far
    (:func:`tensor_omics.check_plateau_condition`)
    and relative effect-size stability
    (:func:`tensor_omics.check_effect_size_plateau_condition`).
    Both are always computed and traced; `plateau_mode` decides which one stops the search
    (by default the effect-size criterion, see `plateau_mode`). A stop by the effect-size
    criterion returns the candidate that triggered it.

    The returned candidate:

    1. A plateau was found, or the grid has a single candidate and it is admissible: the best
    candidate, `plateau_established = True`.
    2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
    median confidence-interval width across studies, `plateau_established = False`.
    3. No candidate was admissible: the first (finest) candidate, with
    `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
    ranges `0.0`, and `plateau_established = False`.

    Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
    for its meaningful part. Errors (`ierr`): `ERR_INVALID_INPUT` before any work when
    `max_n_points_candidate` or `max_n_neighbors_candidate` is below what
    `calc_js_comp_test_candidate_bounds` returns for the same `grid_points_factor` and
    `grid_min_points` (`n_admissible_evaluated` is then 0, every other output undefined); a
    failure to create the bootstrap's random number generator stops the search, the outputs are
    then undefined.

    Parameters
    ----------
    gene_means : np.ndarray[np.float64] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Per-gene mean expression values for all studies
        NaN is permitted for this value.
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Matrix of signed residuals per study
        NaN is permitted for this value.
    n_bootstraps : int
        Number of bootstraps to perform for a candidate pair
        The minimum valid value is `1`.
    join_method : str, one of 'join_min' | 'join_max' | 'join_median'
        How the studies' confidence-interval overlaps are joined into the plateau condition
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    max_n_points_candidate : int
        Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
        at least what `calc_js_comp_test_candidate_bounds` returns for the same
        `grid_points_factor` and `grid_min_points`, else invalid input (checked by the search,
        since the bound depends on those two optional arguments)
        It is recommended to compute this argument from the `max_n_points_candidate` output produced by :func:`tensor_omics.calc_js_comp_test_candidate_bounds`.
        The minimum valid value is `1`.
    max_n_neighbors_candidate : int
        Upper bound on the grid's `n_neighbors` candidates; at least what
        `calc_js_comp_test_candidate_bounds` returns for the same `grid_points_factor` and
        `grid_min_points`, else invalid input (checked by the search, like
        `max_n_points_candidate`)
        It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by :func:`tensor_omics.calc_js_comp_test_candidate_bounds`.
        The minimum valid value is `1`.
    min_residuals_per_bin : int, optional, default 10
        Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
        one threshold for both, so a point the search marks as failed also fails the gate
        The minimum valid value is `0`.
        The default value is `10`.
    min_neighbor_overlap : float, optional, default 0.1
        Minimum fractional overlap of two consecutive neighborhoods (gate 1)
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.1`.
    succeeding_ci_overlap : float, optional, default 0.9
        Minimum fractional overlap of a candidate's confidence interval with the best one so
        far, per `join_method`. `0` makes the first admissible candidate plateau under the
        CI-overlap criterion (for baseline sweeps)
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.
    plateau_mode : str, one of 'plateau_ci_overlap' | 'plateau_effect_size' | 'plateau_both', optional, default 'plateau_effect_size'
        Which plateau criterion decides when the search stops. The default is the change of
        the observed JSD values from one admissible candidate to the next: it measures
        directly how stable the observed value is. The CI-overlap criterion stays available;
        its interval is a robustness interval of each candidate's observed value (each
        study's histograms drawn again from its own counts)

        The default value is `'plateau_effect_size'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    delta_median_threshold : float, optional, default 0.05
        Upper bound on the median relative JSD change for a transition to qualify
        The minimum valid value is `above(0.0)`.
        The default value is `0.05`.
    delta_max_threshold : float, optional, default 0.1
        Upper bound on the largest relative JSD change for a transition to qualify
        The minimum valid value is `above(0.0)`.
        The default value is `0.10`.
    delta_epsilon : float, optional, default 1e-10
        Floor of the denominators of the relative JSD change and relative CI width
        The minimum valid value is `above(0.0)`.
        The default value is `1.0e-10`.
    delta_min_consecutive_transitions : int, optional, default 2
        Number of consecutive qualifying transitions for an effect-size plateau
        The minimum valid value is `1`.
        The default value is `2`.
    m_min : int, optional, default 3
        Smallest bin count the bin-count search tests (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    gamma_occupancy : float, optional, default 1.25
        Geometric growth factor of the occupancy search; unused in the Stone modes
        The minimum valid value is `above(1.0)`.
        The default value is `1.25`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] of each reference point's lower residual-range bound
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] of each reference point's upper residual-range bound
        The minimum valid value is `lower_residual_range_quantile`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Two-sided significance level of the bootstrap interval, in percent per end
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_stone_equal'
        How each reference point's histogram bin count is chosen, for every candidate; the
        Stone modes test every bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_stone_equal'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    neighborhood_size_mode : str, one of 'neighborhood_size_large' | 'neighborhood_size_small', optional, default 'neighborhood_size_large'
        Which neighborhood size the candidate grid pairs with each `n_points`. Small
        neighborhoods are not the default: if they give higher JSCompTest values, that only
        shows they are too small

        The default value is `'neighborhood_size_large'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    delta_abs_tolerance : float, optional, default 0.0001
        Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
        The minimum valid value is `0.0`.
        The default value is `0.0001`.
    divergence_mode : str, one of 'divergence_symmetric' | 'divergence_pairwise_k2', optional, default 'divergence_pairwise_k2'
        How each study's per-point divergence is computed, for the observed JSDs and every
        bootstrap draw

        The default value is `'divergence_pairwise_k2'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    grid_points_factor : float, optional, default 4.0
        Factor in the candidate grid's first reference-point count
        `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
        The minimum valid value is `above(0.0)`.
        The default value is `4.0`.
    grid_min_points : int, optional, default 200
        The candidate grid's floor: the lower clamp of its first `n_points`, the least
        `n_points` it keeps, and its last, coarsest candidate
        The minimum valid value is `1`.
        The maximum valid value is `MAX_POINTS`.
        The default value is `200`.

    Returns
    -------
    dict
        with keys:

        n_points : int
            The finally chosen candidate's `n_points`
        n_neighbors : int
            The finally chosen candidate's `n_neighbors`
        n_bins_per_point : np.ndarray[np.int32] of shape (max_n_points_candidate,), read-only
            The chosen candidate's bin count per reference point; the first `n_points` entries
            are meaningful, the rest are 0
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The chosen candidate's lower residual-range bound (R_low) per reference point; the
            first `n_points` entries are meaningful, `0.0` when no candidate was admissible; the
            rest are `0.0`
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
            A result is a value; call `.copy()` to obtain a modifiable array.
        best_candidate_pair_confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), read-only
            The chosen candidate's bootstrapped JSD confidence interval (a robustness interval of
            its observed JSD); `-1.0` when no candidate was admissible
            A result is a value; call `.copy()` to obtain a modifiable array.
        plateau_established : bool
            `True` when a plateau was found (by the criterion `plateau_mode` selects) or the
            only candidate was admissible; `False` otherwise, also when the smallest-uncertainty
            fallback returns a real candidate and interval
        trace_n_points : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `n_points`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_n_neighbors : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `n_neighbors`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_global_js_divergence : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_lower : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            (`L_{i,t}`)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_upper : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            (`U_{i,t}`)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_width : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            L_{i,t}`)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_width_relative : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study relative confidence-interval width
            (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_overlap : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            interval with the best candidate's so far, relative to its own width (the values
            `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study relative JSD change from the previous admissible
            candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate; 0 for a study
            whose absolute change is at most `delta_abs_tolerance`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta_median : np.ndarray[np.float64] of shape (16,), read-only
            Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta_max : np.ndarray[np.float64] of shape (16,), read-only
            Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_selected_n_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            `1:trace_n_points(t)` of column `t` are meaningful (the rest are 0), so a Python/R
            caller slices `[:trace_n_points[t], t]`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_occupancy_failed : np.ndarray[np.bool_] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point flag: `True` iff the bin-count search
            failed and the point used `m_min` bins. Jagged like trace_selected_n_bins (`False`
            beyond)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_n_pooled_residuals : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
            like trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_min_bin_occupancy : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_mean_bin_occupancy : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_max_bin_occupancy : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_sturges_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_fd_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
            Jagged like trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
            Jagged like trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_shared_residual_range_high : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
            Jagged like trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_fixed_k_search::run_js_comp_test_parameter_search`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    gene_means = _contract.reals(gene_means, "gene_means", ndim=2, order="F")
    residuals = _contract.reals(residuals, "residuals", ndim=3, order="F")
    n_bootstraps = _contract.integer(n_bootstraps, "n_bootstraps")
    join_method = _np.array([_contract.mode(join_method, "join_method", ("join_min", "join_max", "join_median",)).ljust(11)], dtype="S11")
    max_n_points_candidate = _contract.integer(max_n_points_candidate, "max_n_points_candidate")
    max_n_neighbors_candidate = _contract.integer(max_n_neighbors_candidate, "max_n_neighbors_candidate")
    min_residuals_per_bin = _contract.integer(min_residuals_per_bin, "min_residuals_per_bin")
    min_neighbor_overlap = _contract.real(min_neighbor_overlap, "min_neighbor_overlap")
    succeeding_ci_overlap = _contract.real(succeeding_ci_overlap, "succeeding_ci_overlap")
    plateau_mode = _np.array([_contract.mode(plateau_mode, "plateau_mode", ("plateau_ci_overlap", "plateau_effect_size", "plateau_both",)).ljust(19)], dtype="S19")
    delta_median_threshold = _contract.real(delta_median_threshold, "delta_median_threshold")
    delta_max_threshold = _contract.real(delta_max_threshold, "delta_max_threshold")
    delta_epsilon = _contract.real(delta_epsilon, "delta_epsilon")
    delta_min_consecutive_transitions = _contract.integer(delta_min_consecutive_transitions, "delta_min_consecutive_transitions")
    m_min = _contract.integer(m_min, "m_min")
    m_max = _contract.integer(m_max, "m_max")
    gamma_occupancy = _contract.real(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile = _contract.real(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile = _contract.real(upper_residual_range_quantile, "upper_residual_range_quantile")
    two_sided_bootstrapping_significance_level = _contract.real(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed = _contract.integer(random_seed, "random_seed")
    bin_count_mode = _np.array([_contract.mode(bin_count_mode, "bin_count_mode", ("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled",)).ljust(29)], dtype="S29")
    neighborhood_size_mode = _np.array([_contract.mode(neighborhood_size_mode, "neighborhood_size_mode", ("neighborhood_size_large", "neighborhood_size_small",)).ljust(23)], dtype="S23")
    delta_abs_tolerance = _contract.real(delta_abs_tolerance, "delta_abs_tolerance")
    divergence_mode = _np.array([_contract.mode(divergence_mode, "divergence_mode", ("divergence_symmetric", "divergence_pairwise_k2",)).ljust(22)], dtype="S22")
    grid_points_factor = _contract.real(grid_points_factor, "grid_points_factor")
    grid_min_points = _contract.integer(grid_min_points, "grid_min_points")

    # what the inputs already say, rather than asking for it again
    n_studies = gene_means.shape[1]
    max_n_genes_all_studies = gene_means.shape[0]
    max_n_reps_all_studies = residuals.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if residuals.shape[1] != max_n_genes_all_studies:
        raise _ToxArgumentError(f"'residuals' has {residuals.shape[1]} along axis 1, but 'gene_means' implies max_n_genes_all_studies == {max_n_genes_all_studies}; 'residuals' is declared with shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies)", argument="residuals")
    if residuals.shape[2] != n_studies:
        raise _ToxArgumentError(f"'residuals' has {residuals.shape[2]} along axis 2, but 'gene_means' implies n_studies == {n_studies}; 'residuals' is declared with shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies)", argument="residuals")

    # outputs and work arrays, which the caller never sees
    n_points = _ctypes.c_int(0)
    n_neighbors = _ctypes.c_int(0)
    n_bins_per_point = _np.empty((max_n_points_candidate,), dtype=_np.int32, order='C')
    shared_residual_range_low = _np.empty((max_n_points_candidate,), dtype=_np.float64, order='C')
    shared_residual_range_high = _np.empty((max_n_points_candidate,), dtype=_np.float64, order='C')
    best_candidate_pair_confidence_interval = _np.empty((2, n_studies,), dtype=_np.float64, order='F')
    plateau_established = _ctypes.c_bool(0)
    n_admissible_evaluated = _ctypes.c_int(0)
    trace_n_points = _np.empty((16,), dtype=_np.int32, order='C')
    trace_n_neighbors = _np.empty((16,), dtype=_np.int32, order='C')
    trace_global_js_divergence = _np.empty((n_studies, 16,), dtype=_np.float64, order='F')
    trace_ci_lower = _np.empty((n_studies, 16,), dtype=_np.float64, order='F')
    trace_ci_upper = _np.empty((n_studies, 16,), dtype=_np.float64, order='F')
    trace_ci_width = _np.empty((n_studies, 16,), dtype=_np.float64, order='F')
    trace_ci_width_relative = _np.empty((n_studies, 16,), dtype=_np.float64, order='F')
    trace_ci_overlap = _np.empty((n_studies, 16,), dtype=_np.float64, order='F')
    trace_delta = _np.empty((n_studies, 16,), dtype=_np.float64, order='F')
    trace_delta_median = _np.empty((16,), dtype=_np.float64, order='C')
    trace_delta_max = _np.empty((16,), dtype=_np.float64, order='C')
    trace_selected_n_bins = _np.empty((max_n_points_candidate, 16,), dtype=_np.int32, order='F')
    trace_occupancy_failed = _np.empty((max_n_points_candidate, 16,), dtype=_np.bool_, order='F')
    trace_n_pooled_residuals = _np.empty((max_n_points_candidate, 16,), dtype=_np.int32, order='F')
    trace_min_bin_occupancy = _np.empty((max_n_points_candidate, 16,), dtype=_np.int32, order='F')
    trace_mean_bin_occupancy = _np.empty((max_n_points_candidate, 16,), dtype=_np.float64, order='F')
    trace_max_bin_occupancy = _np.empty((max_n_points_candidate, 16,), dtype=_np.int32, order='F')
    trace_sturges_bins = _np.empty((max_n_points_candidate, 16,), dtype=_np.int32, order='F')
    trace_fd_bins = _np.empty((max_n_points_candidate, 16,), dtype=_np.int32, order='F')
    trace_shared_residual_range_low = _np.empty((max_n_points_candidate, 16,), dtype=_np.float64, order='F')
    trace_shared_residual_range_high = _np.empty((max_n_points_candidate, 16,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.run_js_comp_test_parameter_search_c(
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_reps_all_studies)),
        gene_means,
        residuals,
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        join_method,
        _ctypes.byref(_ctypes.c_int(max_n_points_candidate)),
        _ctypes.byref(_ctypes.c_int(max_n_neighbors_candidate)),
        _ctypes.byref(n_points),
        _ctypes.byref(n_neighbors),
        n_bins_per_point,
        shared_residual_range_low,
        shared_residual_range_high,
        best_candidate_pair_confidence_interval,
        _ctypes.byref(plateau_established),
        _ctypes.byref(n_admissible_evaluated),
        trace_n_points,
        trace_n_neighbors,
        trace_global_js_divergence,
        trace_ci_lower,
        trace_ci_upper,
        trace_ci_width,
        trace_ci_width_relative,
        trace_ci_overlap,
        trace_delta,
        trace_delta_median,
        trace_delta_max,
        trace_selected_n_bins,
        trace_occupancy_failed,
        trace_n_pooled_residuals,
        trace_min_bin_occupancy,
        trace_mean_bin_occupancy,
        trace_max_bin_occupancy,
        trace_sturges_bins,
        trace_fd_bins,
        trace_shared_residual_range_low,
        trace_shared_residual_range_high,
        _ctypes.byref(_ctypes.c_int(min_residuals_per_bin)),
        _ctypes.byref(_ctypes.c_double(min_neighbor_overlap)),
        _ctypes.byref(_ctypes.c_double(succeeding_ci_overlap)),
        plateau_mode,
        _ctypes.byref(_ctypes.c_double(delta_median_threshold)),
        _ctypes.byref(_ctypes.c_double(delta_max_threshold)),
        _ctypes.byref(_ctypes.c_double(delta_epsilon)),
        _ctypes.byref(_ctypes.c_int(delta_min_consecutive_transitions)),
        _ctypes.byref(_ctypes.c_int(m_min)),
        _ctypes.byref(_ctypes.c_int(m_max)),
        _ctypes.byref(_ctypes.c_double(gamma_occupancy)),
        _ctypes.byref(_ctypes.c_double(lower_residual_range_quantile)),
        _ctypes.byref(_ctypes.c_double(upper_residual_range_quantile)),
        _ctypes.byref(_ctypes.c_double(two_sided_bootstrapping_significance_level)),
        _ctypes.byref(_ctypes.c_int(random_seed)),
        bin_count_mode,
        neighborhood_size_mode,
        _ctypes.byref(_ctypes.c_double(delta_abs_tolerance)),
        divergence_mode,
        _ctypes.byref(_ctypes.c_double(grid_points_factor)),
        _ctypes.byref(_ctypes.c_int(grid_min_points)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENTS, _RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    n_bins_per_point.flags.writeable = False
    shared_residual_range_low.flags.writeable = False
    shared_residual_range_high.flags.writeable = False
    best_candidate_pair_confidence_interval.flags.writeable = False
    trace_n_points.flags.writeable = False
    trace_n_neighbors.flags.writeable = False
    trace_global_js_divergence.flags.writeable = False
    trace_ci_lower.flags.writeable = False
    trace_ci_upper.flags.writeable = False
    trace_ci_width.flags.writeable = False
    trace_ci_width_relative.flags.writeable = False
    trace_ci_overlap.flags.writeable = False
    trace_delta.flags.writeable = False
    trace_delta_median.flags.writeable = False
    trace_delta_max.flags.writeable = False
    trace_selected_n_bins.flags.writeable = False
    trace_occupancy_failed.flags.writeable = False
    trace_n_pooled_residuals.flags.writeable = False
    trace_min_bin_occupancy.flags.writeable = False
    trace_mean_bin_occupancy.flags.writeable = False
    trace_max_bin_occupancy.flags.writeable = False
    trace_sturges_bins.flags.writeable = False
    trace_fd_bins.flags.writeable = False
    trace_shared_residual_range_low.flags.writeable = False
    trace_shared_residual_range_high.flags.writeable = False

    return {
        "n_points": n_points.value,
        "n_neighbors": n_neighbors.value,
        "n_bins_per_point": n_bins_per_point,
        "shared_residual_range_low": shared_residual_range_low,
        "shared_residual_range_high": shared_residual_range_high,
        "best_candidate_pair_confidence_interval": best_candidate_pair_confidence_interval,
        "plateau_established": plateau_established.value,
        "trace_n_points": trace_n_points[..., :n_admissible_evaluated.value],
        "trace_n_neighbors": trace_n_neighbors[..., :n_admissible_evaluated.value],
        "trace_global_js_divergence": trace_global_js_divergence[..., :n_admissible_evaluated.value],
        "trace_ci_lower": trace_ci_lower[..., :n_admissible_evaluated.value],
        "trace_ci_upper": trace_ci_upper[..., :n_admissible_evaluated.value],
        "trace_ci_width": trace_ci_width[..., :n_admissible_evaluated.value],
        "trace_ci_width_relative": trace_ci_width_relative[..., :n_admissible_evaluated.value],
        "trace_ci_overlap": trace_ci_overlap[..., :n_admissible_evaluated.value],
        "trace_delta": trace_delta[..., :n_admissible_evaluated.value],
        "trace_delta_median": trace_delta_median[..., :n_admissible_evaluated.value],
        "trace_delta_max": trace_delta_max[..., :n_admissible_evaluated.value],
        "trace_selected_n_bins": trace_selected_n_bins[..., :n_admissible_evaluated.value],
        "trace_occupancy_failed": trace_occupancy_failed[..., :n_admissible_evaluated.value],
        "trace_n_pooled_residuals": trace_n_pooled_residuals[..., :n_admissible_evaluated.value],
        "trace_min_bin_occupancy": trace_min_bin_occupancy[..., :n_admissible_evaluated.value],
        "trace_mean_bin_occupancy": trace_mean_bin_occupancy[..., :n_admissible_evaluated.value],
        "trace_max_bin_occupancy": trace_max_bin_occupancy[..., :n_admissible_evaluated.value],
        "trace_sturges_bins": trace_sturges_bins[..., :n_admissible_evaluated.value],
        "trace_fd_bins": trace_fd_bins[..., :n_admissible_evaluated.value],
        "trace_shared_residual_range_low": trace_shared_residual_range_low[..., :n_admissible_evaluated.value],
        "trace_shared_residual_range_high": trace_shared_residual_range_high[..., :n_admissible_evaluated.value],
    }
