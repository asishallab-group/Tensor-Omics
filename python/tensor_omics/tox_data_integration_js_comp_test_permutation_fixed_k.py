r"""tox_data_integration_js_comp_test_permutation_fixed_k

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test, k-NN

Asks whether the assignment of samples to studies (batches) influences the JS-Comp-Test values
at all. Null hypothesis: it does not. The test shuffles the study labels of the samples
(without replacement), reruns the whole JS-Comp-Test on every shuffled data set, and compares
each study's observed global JSD with the values the shuffled data sets give. Study `i`'s p-value is
`(count + 1)/(n_valid + 1)`, where `count` is how many of the `n_valid` valid shuffled data sets give
study `i` a global JSD at least as large as the observed one. The +1 counts the observed data itself.

**Input layout.** All studies arrive in one array `expression(max_n_reps_all_studies,
max_n_genes_all_studies, n_studies)`: gene `g` is the same gene in every study, a gene missing in
a study is NaN there, and a replicate slot `(r, s)` whose values are all NaN is padding (a study
with fewer samples than `max_n_reps_all_studies`). Every other slot is a sample. Samples are
numbered in slot order, replicate fastest, then study.

**Entry points.** One assignment of samples to slots, one k-NN JS-Comp-Test (the fixed
`(n_points, n_neighbors)` of `run_js_comp_test`):
:func:`tensor_omics.run_js_comp_test_on_sample_assignment`.
The whole test, `n_permutations` random assignments plus the observed one:
:func:`tensor_omics.run_js_comp_test_permutation`.
The adaptive counterparts are in ``tox_data_integration_js_comp_test_permutation_adaptive_impl``, the shared helpers in ``tox_data_integration_js_comp_test_permutation_common_impl``.
`calc_js_comp_test_permutation_bounds` and `calc_adaptive_js_comp_test_permutation_bounds`
(``tox_data_integration_js_comp_test_permutation_bounds_impl``) size the work arrays for a caller that allocates its own.

The adaptive routines are the recommended way, as for the JS-Comp-Test itself; the k-NN ones
stay available.

Python binding, generated from tox_data_integration_js_comp_test_permutation_fixed_k. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.run_js_comp_test_on_sample_assignment_c.restype = None
_lib.run_js_comp_test_on_sample_assignment_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
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
_RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "n_points", "n_neighbors", "global_js_divergence", "used_n_points", "used_n_neighbors", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "neighborhood_size_mode", "delta_abs_tolerance", "divergence_mode", "grid_points_factor", "grid_min_points", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, "sample_assignment", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.run_js_comp_test_permutation_c.restype = None
_lib.run_js_comp_test_permutation_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
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
_RUN_JS_COMP_TEST_PERMUTATION_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_points", "n_neighbors", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "neighborhood_size_mode", "delta_abs_tolerance", "divergence_mode", "grid_points_factor", "grid_min_points", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_PERMUTATION_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, None, None, "jsd_null", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def run_js_comp_test_on_sample_assignment(
        expression,
        sample_assignment,
        n_points,
        n_neighbors,
        rerun_parameter_search=False,
        n_bootstraps=1000,
        join_method='join_min',
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
    r"""Run the k-NN JS-Comp-Test on one assignment of samples to studies

    One building block of the study-label Permutation-Test: relabels the samples by `sample_assignment` and runs
    the whole JS-Comp-Test on the relabelled data set. `expression` holds every study (input layout as in the module
    description). `sample_assignment` is a permutation of `1..n_samples`: destination sample `k` receives all genes of source
    sample `sample_assignment(k)`, so every study keeps its slots, sample count and padding. The identity is the observed
    data, and its result is exactly that of a caller doing the steps below on the observed data.

    Steps: per study, the gene means and residuals of the relabelled data
    (:func:`tensor_omics.compute_gene_means`, :func:`tensor_omics.compute_residuals`); the reference points from
    the pooled means (:func:`tensor_omics.pool_means`); then :func:`tensor_omics.run_js_comp_test` with the
    bin-count options and `divergence_mode` as given.

    With `rerun_parameter_search = True` the relabelled data first run
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    with the search options below, and the run uses the `(n_points, n_neighbors)` it returns (also its fallback when no
    candidate was admissible). `n_points` and `n_neighbors` then only size the work arrays, and `used_n_points` and
    `used_n_neighbors` report what the run used. Every call runs the whole JS-Comp-Test, with the rerun option the
    whole parameter search first.

    Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
    number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. In the Fortran expert tier, with the rerun option, also when a capacity is below the
    search's own bound (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`).
    Errors of the search and the run are returned with their code, and `global_js_divergence` is then NaN.
    Known limitation: a gene missing in a whole study travels with the samples, so after relabelling that study may have
    values for it.

    Parameters
    ----------
    expression : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Expression of every sample slot of every gene in every study (see above); NaN for a
        missing value, a missing gene or a padding slot
        NaN is permitted for this value.
    sample_assignment : np.ndarray[np.int32] of shape (n_samples,)
        A permutation of `1..n_samples`: destination sample `k` receives source sample
        `sample_assignment(k)`; the identity is the observed data
        The minimum valid value is `1`.
        The maximum valid value is `n_samples`.
    n_points : int
        Number of reference points of the run (sizes the work arrays only, with the rerun option)
        The minimum valid value is `1`.
    n_neighbors : int
        Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
        rerun option)
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies`.
    rerun_parameter_search : bool, optional, default False
        `True`: run the fixed-k parameter search on the relabelled data first and use its
        `(n_points, n_neighbors)`; `False`: use `n_points` and `n_neighbors`
        The default value is `False`.
    n_bootstraps : int, optional, default 1000
        Number of bootstraps per candidate of the rerun parameter search
        The minimum valid value is `1`.
        The default value is `1000`.
    join_method : str, one of 'join_min' | 'join_max' | 'join_median', optional, default 'join_min'
        How the rerun search joins the studies' confidence-interval overlaps

        The default value is `'join_min'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach in the bin-count search (and,
        in the rerun search, in every consensus bin)
        The minimum valid value is `0`.
        The default value is `10`.
    min_neighbor_overlap : float, optional, default 0.1
        Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.1`.
    succeeding_ci_overlap : float, optional, default 0.9
        Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.
    plateau_mode : str, one of 'plateau_ci_overlap' | 'plateau_effect_size' | 'plateau_both', optional, default 'plateau_effect_size'
        Rerun search: which plateau criterion stops the search; by default, as in the
        searches, the change of the observed JSD values from one admissible candidate to the next

        The default value is `'plateau_effect_size'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    delta_median_threshold : float, optional, default 0.05
        Rerun search: upper bound on the median relative JSD change for a transition to qualify
        The minimum valid value is `above(0.0)`.
        The default value is `0.05`.
    delta_max_threshold : float, optional, default 0.1
        Rerun search: upper bound on the largest relative JSD change for a transition to qualify
        The minimum valid value is `above(0.0)`.
        The default value is `0.10`.
    delta_epsilon : float, optional, default 1e-10
        Rerun search: floor of the denominators of the relative JSD change and CI width
        The minimum valid value is `above(0.0)`.
        The default value is `1.0e-10`.
    delta_min_consecutive_transitions : int, optional, default 2
        Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
        Rerun search: two-sided significance level of the bootstrap interval, in percent per end
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    random_seed : int, optional, default 42
        Rerun search: seed of the bootstrap's random number generator
        The default value is `42`.
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_stone_equal'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_stone_equal'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    neighborhood_size_mode : str, one of 'neighborhood_size_large' | 'neighborhood_size_small', optional, default 'neighborhood_size_large'
        Rerun search: which neighborhood size the k-NN candidate grid pairs with each
        `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
        default: a higher JSCompTest value from small neighborhoods only shows that they are
        too small

        The default value is `'neighborhood_size_large'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    delta_abs_tolerance : float, optional, default 0.0001
        Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
        The minimum valid value is `0.0`.
        The default value is `0.0001`.
    divergence_mode : str, one of 'divergence_symmetric' | 'divergence_pairwise_k2', optional, default 'divergence_pairwise_k2'
        How each study's per-point divergence is computed, in the run (observed and every
        shuffled data set) and in the rerun search

        The default value is `'divergence_pairwise_k2'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    grid_points_factor : float, optional, default 4.0
        Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
        grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
        The minimum valid value is `above(0.0)`.
        The default value is `4.0`.
    grid_min_points : int, optional, default 200
        Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
        floor, its last and coarsest `n_points` candidate
        The minimum valid value is `1`.
        The maximum valid value is `MAX_POINTS`.
        The default value is `200`.

    Returns
    -------
    dict
        with keys:

        global_js_divergence : np.ndarray[np.float64] of shape (n_studies,), read-only
            Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
            two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
            error
            A result is a value; call `.copy()` to obtain a modifiable array.
        used_n_points : int
            Number of reference points the run used (`n_points`, or the search's choice)
        used_n_neighbors : int
            Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_fixed_k::run_js_comp_test_on_sample_assignment`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression = _contract.reals(expression, "expression", ndim=3, order="F")
    sample_assignment = _contract.integers(sample_assignment, "sample_assignment", ndim=1)
    n_points = _contract.integer(n_points, "n_points")
    n_neighbors = _contract.integer(n_neighbors, "n_neighbors")
    rerun_parameter_search = _contract.logical(rerun_parameter_search, "rerun_parameter_search")
    n_bootstraps = _contract.integer(n_bootstraps, "n_bootstraps")
    join_method = _np.array([_contract.mode(join_method, "join_method", ("join_min", "join_max", "join_median",)).ljust(11)], dtype="S11")
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
    n_studies = expression.shape[2]
    max_n_genes_all_studies = expression.shape[1]
    max_n_reps_all_studies = expression.shape[0]
    n_samples = sample_assignment.shape[0]

    # outputs and work arrays, which the caller never sees
    global_js_divergence = _np.empty((n_studies,), dtype=_np.float64, order='C')
    used_n_points = _ctypes.c_int(0)
    used_n_neighbors = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.run_js_comp_test_on_sample_assignment_c(
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_reps_all_studies)),
        expression,
        _ctypes.byref(_ctypes.c_int(n_samples)),
        sample_assignment,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        global_js_divergence,
        _ctypes.byref(used_n_points),
        _ctypes.byref(used_n_neighbors),
        _ctypes.byref(_ctypes.c_bool(rerun_parameter_search)),
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        join_method,
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

    _check_err_code(ierr.value, _RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENTS, _RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    global_js_divergence.flags.writeable = False

    return {
        "global_js_divergence": global_js_divergence,
        "used_n_points": used_n_points.value,
        "used_n_neighbors": used_n_neighbors.value,
    }

def run_js_comp_test_permutation(
        expression,
        n_points,
        n_neighbors,
        n_permutations,
        rerun_parameter_search=False,
        n_bootstraps=1000,
        join_method='join_min',
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
    r"""Permutation-Test of the k-NN JS-Comp-Test that permutes the samples' study labels

    Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value (null
    hypothesis: it does not), with the data, parameters and options of
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`.

    1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
    `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
    shuffled data. Row `i` of `jsd_null` holds its values.
    3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.

    With `rerun_parameter_search = True` every data set, the observed one included, first runs the fixed-k parameter search. `random_seed` then
    also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
    draws. The routine is sequential: a caller who wants parallel permutations calls
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`
    itself, once per assignment.

    Failures: a permutation that fails with invalid input (a safeguard: the relabelling keeps every study's sample count, so the k-NN test does not normally trigger it) gets a NaN row in `jsd_null`, is
    counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
    `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
    returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
    is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own bound
    (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`) fails the same way with the rerun option.
    Known limitation: a gene missing in a whole study travels with the samples.

    Parameters
    ----------
    expression : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Expression of every sample slot of every gene in every study (see above); NaN for a
        missing value, a missing gene or a padding slot
        NaN is permitted for this value.
    n_points : int
        Number of reference points of the run (sizes the work arrays only, with the rerun option)
        The minimum valid value is `1`.
    n_neighbors : int
        Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
        rerun option)
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies`.
    n_permutations : int
        Number of random assignments (rows of `jsd_null`)
        The minimum valid value is `1`.
    rerun_parameter_search : bool, optional, default False
        `True`: run the fixed-k parameter search on the relabelled data first and use its
        `(n_points, n_neighbors)`; `False`: use `n_points` and `n_neighbors`
        The default value is `False`.
    n_bootstraps : int, optional, default 1000
        Number of bootstraps per candidate of the rerun parameter search
        The minimum valid value is `1`.
        The default value is `1000`.
    join_method : str, one of 'join_min' | 'join_max' | 'join_median', optional, default 'join_min'
        How the rerun search joins the studies' confidence-interval overlaps

        The default value is `'join_min'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach in the bin-count search (and,
        in the rerun search, in every consensus bin)
        The minimum valid value is `0`.
        The default value is `10`.
    min_neighbor_overlap : float, optional, default 0.1
        Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.1`.
    succeeding_ci_overlap : float, optional, default 0.9
        Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.
    plateau_mode : str, one of 'plateau_ci_overlap' | 'plateau_effect_size' | 'plateau_both', optional, default 'plateau_effect_size'
        Rerun search: which plateau criterion stops the search; by default, as in the
        searches, the change of the observed JSD values from one admissible candidate to the next

        The default value is `'plateau_effect_size'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    delta_median_threshold : float, optional, default 0.05
        Rerun search: upper bound on the median relative JSD change for a transition to qualify
        The minimum valid value is `above(0.0)`.
        The default value is `0.05`.
    delta_max_threshold : float, optional, default 0.1
        Rerun search: upper bound on the largest relative JSD change for a transition to qualify
        The minimum valid value is `above(0.0)`.
        The default value is `0.10`.
    delta_epsilon : float, optional, default 1e-10
        Rerun search: floor of the denominators of the relative JSD change and CI width
        The minimum valid value is `above(0.0)`.
        The default value is `1.0e-10`.
    delta_min_consecutive_transitions : int, optional, default 2
        Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
        Rerun search: two-sided significance level of the bootstrap interval, in percent per end
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    random_seed : int, optional, default 42
        Seed of the random number generator of the label draws, and of the bootstrap of
        every rerun parameter search (a separate generator)
        The default value is `42`.
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_stone_equal'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_stone_equal'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    neighborhood_size_mode : str, one of 'neighborhood_size_large' | 'neighborhood_size_small', optional, default 'neighborhood_size_large'
        Rerun search: which neighborhood size the k-NN candidate grid pairs with each
        `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
        default: a higher JSCompTest value from small neighborhoods only shows that they are
        too small

        The default value is `'neighborhood_size_large'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    delta_abs_tolerance : float, optional, default 0.0001
        Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
        The minimum valid value is `0.0`.
        The default value is `0.0001`.
    divergence_mode : str, one of 'divergence_symmetric' | 'divergence_pairwise_k2', optional, default 'divergence_pairwise_k2'
        How each study's per-point divergence is computed, in the run (observed and every
        shuffled data set) and in the rerun search

        The default value is `'divergence_pairwise_k2'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    grid_points_factor : float, optional, default 4.0
        Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
        grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
        The minimum valid value is `above(0.0)`.
        The default value is `4.0`.
    grid_min_points : int, optional, default 200
        Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
        floor, its last and coarsest `n_points` candidate
        The minimum valid value is `1`.
        The maximum valid value is `MAX_POINTS`.
        The default value is `200`.

    Returns
    -------
    dict
        with keys:

        global_js_divergence_observed : np.ndarray[np.float64] of shape (n_studies,), read-only
            Weighted global JSD of each study on the observed data (the identity assignment)
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            whose value reached the observed one; NaN when none is valid or on error
            A result is a value; call `.copy()` to obtain a modifiable array.
        jsd_null : np.ndarray[np.float64] of shape (n_permutations, n_studies,), column-major (order='F'), read-only
            Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_failed_permutations : int
            Number of permutations whose JS-Comp-Test failed with invalid input

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_fixed_k::run_js_comp_test_permutation`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression = _contract.reals(expression, "expression", ndim=3, order="F")
    n_points = _contract.integer(n_points, "n_points")
    n_neighbors = _contract.integer(n_neighbors, "n_neighbors")
    n_permutations = _contract.integer(n_permutations, "n_permutations")
    rerun_parameter_search = _contract.logical(rerun_parameter_search, "rerun_parameter_search")
    n_bootstraps = _contract.integer(n_bootstraps, "n_bootstraps")
    join_method = _np.array([_contract.mode(join_method, "join_method", ("join_min", "join_max", "join_median",)).ljust(11)], dtype="S11")
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
    n_studies = expression.shape[2]
    max_n_genes_all_studies = expression.shape[1]
    max_n_reps_all_studies = expression.shape[0]

    # outputs and work arrays, which the caller never sees
    global_js_divergence_observed = _np.empty((n_studies,), dtype=_np.float64, order='C')
    p_values = _np.empty((n_studies,), dtype=_np.float64, order='C')
    jsd_null = _np.empty((n_permutations, n_studies,), dtype=_np.float64, order='F')
    n_failed_permutations = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.run_js_comp_test_permutation_c(
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_reps_all_studies)),
        expression,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_permutations)),
        global_js_divergence_observed,
        p_values,
        jsd_null,
        _ctypes.byref(n_failed_permutations),
        _ctypes.byref(_ctypes.c_bool(rerun_parameter_search)),
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        join_method,
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

    _check_err_code(ierr.value, _RUN_JS_COMP_TEST_PERMUTATION_ARGUMENTS, _RUN_JS_COMP_TEST_PERMUTATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    global_js_divergence_observed.flags.writeable = False
    p_values.flags.writeable = False
    jsd_null.flags.writeable = False

    return {
        "global_js_divergence_observed": global_js_divergence_observed,
        "p_values": p_values,
        "jsd_null": jsd_null,
        "n_failed_permutations": n_failed_permutations.value,
    }
