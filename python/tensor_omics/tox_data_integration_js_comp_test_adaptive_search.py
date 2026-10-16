r"""tox_data_integration_js_comp_test_adaptive_search

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Adaptive Parameter Search

The data-driven choice of the JS-Comp-Test's adaptive neighborhoods, each grown from `k_start`
pooled entries in rounds of `k_step` up to `k_max` while the residual dispersion stays stable.
:func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`
walks the `(k_start, k_step, k_max)` candidates of
:func:`tensor_omics.generate_adaptive_js_comp_test_candidates`
and stops at the first admissible candidate whose bootstrapped JSDs have reached a plateau, by
the criteria of
:func:`tensor_omics.check_plateau_condition` and
:func:`tensor_omics.check_effect_size_plateau_condition`;
every candidate it tries is logged with an `ADAPTIVE_CANDIDATE_*` status.
:func:`tensor_omics.construct_adaptive_neighborhoods`
builds the neighborhoods for the chosen knobs. `calc_adaptive_js_comp_test_bounds` sizes the
work arrays for a caller that allocates its own.

Python binding, generated from tox_data_integration_js_comp_test_adaptive_search. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.generate_adaptive_js_comp_test_candidates_c.restype = None
_lib.generate_adaptive_js_comp_test_candidates_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GENERATE_ADAPTIVE_JS_COMP_TEST_CANDIDATES_ARGUMENTS = ("max_n_genes_all_studies", "n_studies", "candidates_k_start_k_step_k_max", "n_candidates", "ierr",)

_lib.run_js_comp_test_adaptive_parameter_search_c.restype = None
_lib.run_js_comp_test_adaptive_parameter_search_c.argtypes = (
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
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
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
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
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
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_ADAPTIVE_PARAMETER_SEARCH_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "k_start", "k_step", "k_max", "n_points", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_k_start", "trace_k_step", "trace_k_max", "trace_n_points", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_ci_overlap", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "n_candidates_tried", "candidate_k_start", "candidate_n_points", "candidate_status", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_PARAMETER_SEARCH_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", None, None, None, None, "n_bins_per_point", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def generate_adaptive_js_comp_test_candidates(
        max_n_genes_all_studies,
        n_studies,
):
    r"""Generate the ascending adaptive (k_start, k_step, k_max) candidate sequence

    The adaptive counterpart of
    :func:`tensor_omics.generate_js_comp_test_candidates`.
    An adaptive candidate fixes only its growth knobs (each neighborhood grows from `k_start`
    pooled entries in rounds of `k_step` up to `k_max`); its number of reference points emerges.
    `k_start` shrinks from one candidate to the next, so the candidates **ascend** in the number
    of reference points they are expected to produce.

    With `N = max_n_genes_all_studies * n_studies` (the padded pool size) and
    `floor_k = max(10, 2*n_studies)`:

    - `k_start_1 = max(floor_k, ceiling(0.02 * N))`
    - `k_start_t = floor(k_start_1 * 0.8**(t-1))`
    - stop once `k_start_t < floor_k`, or once 16 candidates are accepted
    - `k_step_t = ceiling(0.25 * k_start_t)`
    - `k_max_t = ceiling(4 * k_start_t)`

    The constants are ``ADAPTIVE_START_FRACTION`` (0.02),
    ``ADAPTIVE_GAMMA`` (0.8),
    ``ADAPTIVE_K_STEP_FRACTION`` (0.25),
    ``ADAPTIVE_K_MAX_FACTOR`` (4),
    ``ADAPTIVE_K_START_MIN_ABS`` (10) and
    ``MAX_CANDIDATE_PAIRS`` (16).
    The values strictly descend without deduplication. At least one candidate is always
    returned, even for data smaller than `k_start_1`, whose construction then fails with a
    too-few-means status in the search.

    Not checked here: the values fit 32-bit integers only while
    `k_start_1 <= floor(huge(1)/ADAPTIVE_K_MAX_FACTOR)` (``ADAPTIVE_K_MAX_FACTOR``);
    :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`
    rejects inputs beyond that.

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
    candidates_k_start_k_step_k_max : np.ndarray[np.int32] of shape (3, 16,), column-major (order='F'), read-only
        Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending
        The first `n_candidates` elements will hold the results.
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_adaptive_search::generate_adaptive_js_comp_test_candidates`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    candidates_k_start_k_step_k_max = np.empty((3, 16,), dtype=np.int32, order='F')
    n_candidates = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.generate_adaptive_js_comp_test_candidates_c(
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(n_studies)),
        candidates_k_start_k_step_k_max,
        ctypes.byref(n_candidates),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GENERATE_ADAPTIVE_JS_COMP_TEST_CANDIDATES_ARGUMENTS)

    # a result is a value: modify a copy, not this
    candidates_k_start_k_step_k_max.flags.writeable = False

    return candidates_k_start_k_step_k_max[..., :n_candidates.value]

def run_js_comp_test_adaptive_parameter_search(
        gene_means,
        residuals,
        n_bootstraps,
        join_method,
        max_n_points_candidate,
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
        tau=0.1,
        mad_distance_factor=1.0,
        max_pooled_residuals=0,
        min_study_neighbors=1,
        bin_count_mode='bin_count_occupancy',
):
    r"""Search the ascending adaptive (k_start, k_step, k_max) candidate sequence for a stable JSD parameter setting (Issue #217)

    The adaptive-neighborhood counterpart of
    :func:`tensor_omics.run_js_comp_test_parameter_search`.
    Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
    :func:`tensor_omics.generate_adaptive_js_comp_test_candidates`
    in order. Their `k_start` shrinks, so the number of reference points that emerges ascends.
    The plateau criteria, `plateau_mode` and the fallbacks are exactly the fixed-k search's.

    **Per candidate**, in order:

    1. Grow its neighborhoods exactly as
    :func:`tensor_omics.construct_adaptive_neighborhoods`
    does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
    `min_study_neighbors`.
    2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    reference points emerged.
    3. Gate 1:
    :func:`tensor_omics.check_neighborhood_overlaps`
    on the pooled ranges.
    4. The pipeline of
    :func:`tensor_omics.run_js_comp_test_adaptive`
    up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
    (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
    count
    (:func:`tensor_omics.check_mean_pmf_min_counts`).
    5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
    first plateau `plateau_mode` selects.

    **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
    `candidate_n_points` and `candidate_status`, until the search stops. The status is the first
    failing check, in the order above:
    ``ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED``
    (too few means or residuals),
    ``ADAPTIVE_CANDIDATE_EMPTY_STUDY``
    (a study neighborhood below `min_study_neighbors`, or a study without a non-NaN residual at
    some point),
    ``ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED``,
    ``ADAPTIVE_CANDIDATE_OVERLAP_FAILED``,
    ``ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED``,
    or else ``ADAPTIVE_CANDIDATE_EVALUATED``
    (admissible, with a trace column).

    **Result**, one of three cases, as in the fixed-k search:

    1. A plateau was found, or the sequence had a single candidate and it passed both gates:
    that candidate's knobs, `plateau_established = True`.
    2. No plateau but at least one admissible candidate: the admissible candidate with the
    smallest median confidence-interval width, `plateau_established = False`.
    3. Nothing admissible, including a single candidate that was rejected: candidate 1's knobs,
    confidence interval `-1.0`, `plateau_established = False`.

    `n_points` is the selected candidate's logged point count; in case 3 the bin counts are
    `m_min` and the ranges `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
    `construct_adaptive_neighborhoods` with the returned knobs and the same construction
    options, then `run_js_comp_test_adaptive` with the same bin-count options
    (`bin_count_mode` included): its `global_js_divergence` equals the selected candidate's
    `trace_global_js_divergence` column.

    **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
    place of `trace_n_neighbors`; the per-point traces are jagged in the same way.

    **Memory.** The per-point arrays are sized by `max_n_points_candidate` (from
    :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`);
    the construction buffers by `N = max_n_genes_all_studies * n_studies` (the largest,
    `tmp_n_neighbors_per_point`, `n_studies * N` integers); the pooling buffers hold
    `max_n_reps_all_studies * N` values each, the size of `residuals`.

    A failure to create the bootstrap's random number generator is returned in `ierr` and
    stops the search; the outputs are then undefined.

    Parameters
    ----------
    gene_means : np.ndarray[np.float64] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Mean expression of every gene in every study, NaN for a missing gene
        NaN is permitted for this value.
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Signed residuals of every replicate of every gene in every study, NaN for a missing value
        NaN is permitted for this value.
    n_bootstraps : int
        Number of bootstraps to perform for a candidate
        The minimum valid value is `1`.
    join_method : str, one of 'join_min' | 'join_max' | 'join_median'
        How the studies' confidence-interval overlaps are joined into the plateau condition

    max_n_points_candidate : int
        Most reference points a candidate may emerge with; a candidate with more is rejected
        with the capacity status. The per-point outputs and traces below are sized by it, so
        the caller has to know it up front
        It is recommended to compute this argument from the `max_n_points_candidate` output produced by :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`.
        The minimum valid value is `1`.
    min_residuals_per_bin : int, optional, default 10
        Minimum count every bin must reach, both in the bin-count search (the floor in every
        `bin_count_mode`) and in the second admissibility gate (the consensus pmf)
        The minimum valid value is `0`.
        The default value is `10`.
    min_neighbor_overlap : float, optional, default 0.1
        Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
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
    tau : float, optional, default 0.1
        Largest relative dispersion increase a growth round may cause and still be committed
        The minimum valid value is `0.0`.
        The default value is `0.1`.
    mad_distance_factor : float, optional, default 1.0
        Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
        The minimum valid value is `0.0`.
        The default value is `1.0`.
    max_pooled_residuals : int, optional, default 0
        Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
        The minimum valid value is `0`.
        The default value is `0`.
    min_study_neighbors : int, optional, default 1
        Fewest entries of each study every neighborhood must have for a candidate to be
        admissible
        The minimum valid value is `1`.
        The default value is `1`.
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

        k_start : int
            The finally chosen candidate's `k_start`
        k_step : int
            The finally chosen candidate's `k_step`
        k_max : int
            The finally chosen candidate's `k_max`
        n_points : int
            Number of reference points the finally chosen candidate's construction produced
        n_bins_per_point : np.ndarray[np.int32] of shape (max_n_points_candidate,), read-only
            The chosen candidate's bin count per reference point; the first `n_points` entries
            are meaningful. Without an admissible candidate, the first
            `min(n_points, max_n_points_candidate)` entries are `m_min`. The rest are 0
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The chosen candidate's lower residual-range bound (R_low) per reference point; the
            first `n_points` entries are meaningful, `0.0` without an admissible candidate; the
            rest are `0.0`
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
            A result is a value; call `.copy()` to obtain a modifiable array.
        best_candidate_confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), read-only
            The bootstrapped JSD confidence interval of the finally chosen candidate;
            `-1.0` throughout where no admissible candidate backs the result
            A result is a value; call `.copy()` to obtain a modifiable array.
        plateau_established : bool
            `True` when a plateau was found (by the criterion `plateau_mode` selects) or the
            only candidate was admissible; `False` otherwise
        trace_k_start : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `k_start`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_k_step : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `k_step`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_k_max : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `k_max`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_n_points : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate number of reference points
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_global_js_divergence : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_lower : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_upper : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_width : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study confidence-interval width
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_width_relative : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study relative confidence-interval width
            (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_ci_overlap : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            interval with the best candidate's so far, relative to its own width; `-1.0` for the
            first admissible candidate
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study relative JSD change from the previous admissible
            candidate; `-1.0` for the first
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
            Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            jagged as trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_min_bin_occupancy : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_mean_bin_occupancy : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_max_bin_occupancy : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_sturges_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
            trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_fd_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
            as trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            jagged as trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_shared_residual_range_high : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            jagged as trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        candidate_k_start : np.ndarray[np.int32] of shape (16,), read-only
            `k_start` of every candidate tried, in sequence order
            The first `n_candidates_tried` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        candidate_n_points : np.ndarray[np.int32] of shape (16,), read-only
            Number of reference points each tried candidate's construction produced (0 for too
            few means)
            The first `n_candidates_tried` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        candidate_status : np.ndarray[np.int32] of shape (16,), read-only
            What happened to every tried candidate: `0` evaluated,
            `1` construction failed,
            `2` a study neighborhood below `min_study_neighbors`
            or a study without a non-NaN residual at some point,
            `3` capacity exceeded,
            `4` overlap gate failed,
            `5` consensus-pmf count gate failed
            The first `n_candidates_tried` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_adaptive_search::run_js_comp_test_adaptive_parameter_search`, whose argument names are
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
    k_start = ctypes.c_int(0)
    k_step = ctypes.c_int(0)
    k_max = ctypes.c_int(0)
    n_points = ctypes.c_int(0)
    n_bins_per_point = np.empty((max_n_points_candidate,), dtype=np.int32, order='C')
    shared_residual_range_low = np.empty((max_n_points_candidate,), dtype=np.float64, order='C')
    shared_residual_range_high = np.empty((max_n_points_candidate,), dtype=np.float64, order='C')
    best_candidate_confidence_interval = np.empty((2, n_studies,), dtype=np.float64, order='F')
    plateau_established = ctypes.c_bool(0)
    n_admissible_evaluated = ctypes.c_int(0)
    trace_k_start = np.empty((16,), dtype=np.int32, order='C')
    trace_k_step = np.empty((16,), dtype=np.int32, order='C')
    trace_k_max = np.empty((16,), dtype=np.int32, order='C')
    trace_n_points = np.empty((16,), dtype=np.int32, order='C')
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
    n_candidates_tried = ctypes.c_int(0)
    candidate_k_start = np.empty((16,), dtype=np.int32, order='C')
    candidate_n_points = np.empty((16,), dtype=np.int32, order='C')
    candidate_status = np.empty((16,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_adaptive_parameter_search_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        gene_means,
        residuals,
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        join_method,
        ctypes.byref(ctypes.c_int(max_n_points_candidate)),
        ctypes.byref(k_start),
        ctypes.byref(k_step),
        ctypes.byref(k_max),
        ctypes.byref(n_points),
        n_bins_per_point,
        shared_residual_range_low,
        shared_residual_range_high,
        best_candidate_confidence_interval,
        ctypes.byref(plateau_established),
        ctypes.byref(n_admissible_evaluated),
        trace_k_start,
        trace_k_step,
        trace_k_max,
        trace_n_points,
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
        ctypes.byref(n_candidates_tried),
        candidate_k_start,
        candidate_n_points,
        candidate_status,
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
        ctypes.byref(ctypes.c_double(tau)),
        ctypes.byref(ctypes.c_double(mad_distance_factor)),
        ctypes.byref(ctypes.c_int(max_pooled_residuals)),
        ctypes.byref(ctypes.c_int(min_study_neighbors)),
        bin_count_mode,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_ADAPTIVE_PARAMETER_SEARCH_ARGUMENTS, _RUN_JS_COMP_TEST_ADAPTIVE_PARAMETER_SEARCH_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    n_bins_per_point.flags.writeable = False
    shared_residual_range_low.flags.writeable = False
    shared_residual_range_high.flags.writeable = False
    best_candidate_confidence_interval.flags.writeable = False
    trace_k_start.flags.writeable = False
    trace_k_step.flags.writeable = False
    trace_k_max.flags.writeable = False
    trace_n_points.flags.writeable = False
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
    candidate_k_start.flags.writeable = False
    candidate_n_points.flags.writeable = False
    candidate_status.flags.writeable = False

    return {
        "k_start": k_start.value,
        "k_step": k_step.value,
        "k_max": k_max.value,
        "n_points": n_points.value,
        "n_bins_per_point": n_bins_per_point,
        "shared_residual_range_low": shared_residual_range_low,
        "shared_residual_range_high": shared_residual_range_high,
        "best_candidate_confidence_interval": best_candidate_confidence_interval,
        "plateau_established": plateau_established.value,
        "trace_k_start": trace_k_start[..., :n_admissible_evaluated.value],
        "trace_k_step": trace_k_step[..., :n_admissible_evaluated.value],
        "trace_k_max": trace_k_max[..., :n_admissible_evaluated.value],
        "trace_n_points": trace_n_points[..., :n_admissible_evaluated.value],
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
        "candidate_k_start": candidate_k_start[..., :n_candidates_tried.value],
        "candidate_n_points": candidate_n_points[..., :n_candidates_tried.value],
        "candidate_status": candidate_status[..., :n_candidates_tried.value],
    }
