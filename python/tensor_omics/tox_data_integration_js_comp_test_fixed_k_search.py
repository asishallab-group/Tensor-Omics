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
`calc_js_comp_test_candidate_bounds` sizes the work arrays for a caller that allocates its own.

Python binding, generated from tox_data_integration_js_comp_test_fixed_k_search. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.generate_js_comp_test_candidates_c.restype = None
_lib.generate_js_comp_test_candidates_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GENERATE_JS_COMP_TEST_CANDIDATES_ARGUMENTS = ("max_n_genes_all_studies", "candidates_n_points_n_neighbors", "n_candidates", "ierr",)

_lib.run_js_comp_test_parameter_search_c.restype = None
_lib.run_js_comp_test_parameter_search_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "max_n_neighbors_candidate", "n_points", "n_neighbors", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_pair_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_ci_overlap", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", None, None, None, None, "n_bins_per_point", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def generate_js_comp_test_candidates(
        max_n_genes_all_studies,
):
    r"""Generate the GAMMA-decay (n_points, n_neighbors) candidate grid

    Starting from `n_points_high = clamp(ceiling(4*sqrt(G)), MIN_POINTS, MAX_POINTS)` (`G` =
    `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
    `n_points_low = max(MIN_POINTS, ceiling(0.2*n_points_high))`, and pairs each distinct
    `n_points` with up to `size(KX_FACTORS)` distinct `n_neighbors = max(1, floor(G/(f*n_points_high)))`.
    Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
    candidate. Bin counts are not part of a candidate; they are chosen per reference point later.

    Parameters
    ----------
    max_n_genes_all_studies : int
        Maximum number of genes across all studies
        The minimum valid value is `1`.

    Returns
    -------
    candidates_n_points_n_neighbors : np.ndarray[np.int32] of shape (2, 16,), column-major (order='F'), read-only
        Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
        The first `n_candidates` elements will hold the results.
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_fixed_k_search::generate_js_comp_test_candidates`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    candidates_n_points_n_neighbors = np.empty((2, 16,), dtype=np.int32, order='F')
    n_candidates = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.generate_js_comp_test_candidates_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        candidates_n_points_n_neighbors,
        ctypes.byref(n_candidates),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GENERATE_JS_COMP_TEST_CANDIDATES_ARGUMENTS)

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
        plateau_mode='plateau_ci_overlap',
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
        bin_count_mode='bin_count_occupancy',
):
    r"""Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting

    Walks the candidates of
    :func:`tensor_omics.generate_js_comp_test_candidates`
    from finest to coarsest resolution. For each candidate:

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
    Both are always computed and traced; `plateau_mode` decides which one stops the search. A
    stop by the effect-size criterion returns the candidate that triggered it.

    The returned candidate:

    1. A plateau was found, or the grid has a single candidate and it is admissible: the best
    candidate, `plateau_established = True`.
    2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
    median confidence-interval width across studies, `plateau_established = False`.
    3. No candidate was admissible: the first (finest) candidate, with
    `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
    ranges `0.0`, and `plateau_established = False`.

    Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
    for its meaningful part. A failure to create the bootstrap's random number generator is
    returned in `ierr` and stops the search; the outputs are then undefined.

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

    max_n_points_candidate : int
        Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
        at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
        It is recommended to compute this argument from the `max_n_points_candidate` output produced by :func:`tensor_omics.calc_js_comp_test_candidate_bounds`.
        The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
    max_n_neighbors_candidate : int
        Upper bound on the grid's `n_neighbors` candidates; at least what
        `calc_js_comp_test_candidate_bounds` returns, else invalid input
        It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by :func:`tensor_omics.calc_js_comp_test_candidate_bounds`.
        The minimum valid value is `max(1, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2*real(max(MIN_POINTS, min(ceiling(4.0*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
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
    plateau_mode : str, one of 'plateau_ci_overlap' | 'plateau_effect_size' | 'plateau_both', optional, default 'plateau_ci_overlap'
        Which plateau criterion decides when the search stops

        The default value is `'plateau_ci_overlap'`.
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
        Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
        with invalid input reported on upper_residual_range_quantile.
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
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_occupancy'
        How each reference point's histogram bin count is chosen, for every candidate; the
        Stone modes test every bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_occupancy'`.

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
            The chosen candidate's bootstrapped JSD confidence interval; `-1.0` when no candidate
            was admissible
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
            candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_fixed_k_search::run_js_comp_test_parameter_search`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_means = np.asfortranarray(gene_means, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_means' must be an array of np.float64: {error}") from None
    if gene_means.ndim != 2:
        raise ValueError(f"'gene_means' must have 2 dimensions, but has {gene_means.ndim}")
    try:
        residuals = np.asfortranarray(residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residuals' must be an array of np.float64: {error}") from None
    if residuals.ndim != 3:
        raise ValueError(f"'residuals' must have 3 dimensions, but has {residuals.ndim}")
    join_method = np.array([str(join_method).lower().encode().ljust(11)], dtype="S11")
    plateau_mode = np.array([str(plateau_mode).lower().encode().ljust(19)], dtype="S19")
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

    # what the inputs already say, rather than asking for it again
    n_studies = gene_means.shape[1]
    max_n_genes_all_studies = gene_means.shape[0]
    max_n_reps_all_studies = residuals.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if residuals.shape[2] != n_studies:
        raise ValueError(f"'residuals' has {residuals.shape[2]} along axis 2, but "
            f"'gene_means' implies n_studies == {n_studies}"
        )
    if residuals.shape[1] != max_n_genes_all_studies:
        raise ValueError(f"'residuals' has {residuals.shape[1]} along axis 1, but "
            f"'gene_means' implies max_n_genes_all_studies == {max_n_genes_all_studies}"
        )

    # outputs and work arrays, which the caller never sees
    n_points = ctypes.c_int(0)
    n_neighbors = ctypes.c_int(0)
    n_bins_per_point = np.empty((max_n_points_candidate,), dtype=np.int32, order='C')
    shared_residual_range_low = np.empty((max_n_points_candidate,), dtype=np.float64, order='C')
    shared_residual_range_high = np.empty((max_n_points_candidate,), dtype=np.float64, order='C')
    best_candidate_pair_confidence_interval = np.empty((2, n_studies,), dtype=np.float64, order='F')
    plateau_established = ctypes.c_bool(0)
    n_admissible_evaluated = ctypes.c_int(0)
    trace_n_points = np.empty((16,), dtype=np.int32, order='C')
    trace_n_neighbors = np.empty((16,), dtype=np.int32, order='C')
    trace_global_js_divergence = np.empty((n_studies, 16,), dtype=np.float64, order='F')
    trace_ci_lower = np.empty((n_studies, 16,), dtype=np.float64, order='F')
    trace_ci_upper = np.empty((n_studies, 16,), dtype=np.float64, order='F')
    trace_ci_width = np.empty((n_studies, 16,), dtype=np.float64, order='F')
    trace_ci_width_relative = np.empty((n_studies, 16,), dtype=np.float64, order='F')
    trace_ci_overlap = np.empty((n_studies, 16,), dtype=np.float64, order='F')
    trace_delta = np.empty((n_studies, 16,), dtype=np.float64, order='F')
    trace_delta_median = np.empty((16,), dtype=np.float64, order='C')
    trace_delta_max = np.empty((16,), dtype=np.float64, order='C')
    trace_selected_n_bins = np.empty((max_n_points_candidate, 16,), dtype=np.int32, order='F')
    trace_occupancy_failed = np.empty((max_n_points_candidate, 16,), dtype=np.bool_, order='F')
    trace_n_pooled_residuals = np.empty((max_n_points_candidate, 16,), dtype=np.int32, order='F')
    trace_min_bin_occupancy = np.empty((max_n_points_candidate, 16,), dtype=np.int32, order='F')
    trace_mean_bin_occupancy = np.empty((max_n_points_candidate, 16,), dtype=np.float64, order='F')
    trace_max_bin_occupancy = np.empty((max_n_points_candidate, 16,), dtype=np.int32, order='F')
    trace_sturges_bins = np.empty((max_n_points_candidate, 16,), dtype=np.int32, order='F')
    trace_fd_bins = np.empty((max_n_points_candidate, 16,), dtype=np.int32, order='F')
    trace_shared_residual_range_low = np.empty((max_n_points_candidate, 16,), dtype=np.float64, order='F')
    trace_shared_residual_range_high = np.empty((max_n_points_candidate, 16,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_parameter_search_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        gene_means,
        residuals,
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        join_method,
        ctypes.byref(ctypes.c_int(max_n_points_candidate)),
        ctypes.byref(ctypes.c_int(max_n_neighbors_candidate)),
        ctypes.byref(n_points),
        ctypes.byref(n_neighbors),
        n_bins_per_point,
        shared_residual_range_low,
        shared_residual_range_high,
        best_candidate_pair_confidence_interval,
        ctypes.byref(plateau_established),
        ctypes.byref(n_admissible_evaluated),
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
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_double(min_neighbor_overlap)),
        ctypes.byref(ctypes.c_double(succeeding_ci_overlap)),
        plateau_mode,
        ctypes.byref(ctypes.c_double(delta_median_threshold)),
        ctypes.byref(ctypes.c_double(delta_max_threshold)),
        ctypes.byref(ctypes.c_double(delta_epsilon)),
        ctypes.byref(ctypes.c_int(delta_min_consecutive_transitions)),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_double(gamma_occupancy)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(two_sided_bootstrapping_significance_level)),
        ctypes.byref(ctypes.c_int(random_seed)),
        bin_count_mode,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENTS, _RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENT_SOURCES)

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
