r"""tox_data_integration_js_comp_test_permutation_adaptive

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test, adaptive

The adaptive counterpart of ``tox_data_integration_js_comp_test_permutation_fixed_k_impl``, with the same input layout, assignment
convention and p-values: one assignment of samples to slots, one adaptive JS-Comp-Test
(:func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`), and the whole test
(:func:`tensor_omics.run_js_comp_test_adaptive_permutation`). The shared helpers are in
``tox_data_integration_js_comp_test_permutation_common_impl``.

Python binding, generated from tox_data_integration_js_comp_test_permutation_adaptive. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.run_js_comp_test_adaptive_on_sample_assignment_c.restype = None
_lib.run_js_comp_test_adaptive_on_sample_assignment_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "k_start", "k_step", "k_max", "global_js_divergence", "used_k_start", "used_k_step", "used_k_max", "n_points", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "delta_abs_tolerance", "divergence_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, "sample_assignment", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.run_js_comp_test_adaptive_permutation_c.restype = None
_lib.run_js_comp_test_adaptive_permutation_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
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
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_ADAPTIVE_PERMUTATION_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "k_start", "k_step", "k_max", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "delta_abs_tolerance", "divergence_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_PERMUTATION_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, None, None, None, "jsd_null", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def run_js_comp_test_adaptive_on_sample_assignment(
        expression,
        sample_assignment,
        k_start,
        k_step,
        k_max,
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
        tau=0.1,
        mad_distance_factor=1.0,
        max_pooled_residuals=0,
        min_study_neighbors=1,
        bin_count_mode='bin_count_stone_equal',
        delta_abs_tolerance=0.0001,
        divergence_mode='divergence_pairwise_k2',
):
    r"""Run the adaptive JS-Comp-Test on one assignment of samples to studies

    The adaptive counterpart of
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`,
    with the same input layout and assignment convention (the identity is the observed data).

    Steps: per study, the gene means and residuals of the relabelled data; the adaptive neighborhoods for the
    growth knobs `k_start`, `k_step`, `k_max` and the construction options
    (:func:`tensor_omics.construct_adaptive_neighborhoods`); then
    :func:`tensor_omics.run_js_comp_test_adaptive` on its ranges, with `divergence_mode` as given. `n_points` is the number
    of neighborhoods the construction built.

    With `rerun_parameter_search = True` the relabelled data first run
    :func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`
    (with its own capacity, so it rejects the same candidates as a direct call), and the construction and run use the knobs
    it returns (also its fallback when no candidate was admissible). `k_start`, `k_step` and `k_max` are then not used, except
    `k_start` to size the work arrays, and `used_k_start`, `used_k_step`, `used_k_max` report what was used. Every call runs the
    whole JS-Comp-Test, with the rerun option the whole parameter search first.

    Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
    number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. Also when the construction does not succeed: its status is not ok (too few means, a
    neighborhood with too few residuals or with fewer than `min_study_neighbors` entries of some study), it builds no neighborhood, or more
    neighborhoods than the work arrays hold. In the Fortran expert tier, with the rerun option, also when the capacity is below
    the search's own.
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
    k_start : int
        Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
        The minimum valid value is `1`.
    k_step : int
        Pooled entries staged per adaptive growth round
        The minimum valid value is `1`.
    k_max : int
        Largest number of pooled entries in one neighborhood
        The minimum valid value is `k_start`.
    rerun_parameter_search : bool, optional, default False
        `True`: run the adaptive parameter search on the relabelled data first and use its
        knobs; `False`: use `k_start`, `k_step` and `k_max`
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
        Fewest entries of each study every neighborhood must have (construction status; in
        the rerun search, candidate admissibility)
        The minimum valid value is `1`.
        The default value is `1`.
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_stone_equal'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_stone_equal'`.
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

    Returns
    -------
    dict
        with keys:

        global_js_divergence : np.ndarray[np.float64] of shape (n_studies,), read-only
            Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
            two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
            error
            A result is a value; call `.copy()` to obtain a modifiable array.
        used_k_start : int
            `k_start` the construction used (`k_start`, or the search's choice)
        used_k_step : int
            `k_step` the construction used (`k_step`, or the search's choice)
        used_k_max : int
            `k_max` the construction used (`k_max`, or the search's choice)
        n_points : int
            Number of neighborhoods the construction built (0 when it built none)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_adaptive::run_js_comp_test_adaptive_on_sample_assignment`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression = _contract.reals(expression, "expression", ndim=3, order="F")
    sample_assignment = _contract.integers(sample_assignment, "sample_assignment", ndim=1)
    k_start = _contract.integer(k_start, "k_start")
    k_step = _contract.integer(k_step, "k_step")
    k_max = _contract.integer(k_max, "k_max")
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
    tau = _contract.real(tau, "tau")
    mad_distance_factor = _contract.real(mad_distance_factor, "mad_distance_factor")
    max_pooled_residuals = _contract.integer(max_pooled_residuals, "max_pooled_residuals")
    min_study_neighbors = _contract.integer(min_study_neighbors, "min_study_neighbors")
    bin_count_mode = _np.array([_contract.mode(bin_count_mode, "bin_count_mode", ("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled",)).ljust(29)], dtype="S29")
    delta_abs_tolerance = _contract.real(delta_abs_tolerance, "delta_abs_tolerance")
    divergence_mode = _np.array([_contract.mode(divergence_mode, "divergence_mode", ("divergence_symmetric", "divergence_pairwise_k2",)).ljust(22)], dtype="S22")

    # what the inputs already say, rather than asking for it again
    n_studies = expression.shape[2]
    max_n_genes_all_studies = expression.shape[1]
    max_n_reps_all_studies = expression.shape[0]
    n_samples = sample_assignment.shape[0]

    # outputs and work arrays, which the caller never sees
    global_js_divergence = _np.empty((n_studies,), dtype=_np.float64, order='C')
    used_k_start = _ctypes.c_int(0)
    used_k_step = _ctypes.c_int(0)
    used_k_max = _ctypes.c_int(0)
    n_points = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.run_js_comp_test_adaptive_on_sample_assignment_c(
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_reps_all_studies)),
        expression,
        _ctypes.byref(_ctypes.c_int(n_samples)),
        sample_assignment,
        _ctypes.byref(_ctypes.c_int(k_start)),
        _ctypes.byref(_ctypes.c_int(k_step)),
        _ctypes.byref(_ctypes.c_int(k_max)),
        global_js_divergence,
        _ctypes.byref(used_k_start),
        _ctypes.byref(used_k_step),
        _ctypes.byref(used_k_max),
        _ctypes.byref(n_points),
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
        _ctypes.byref(_ctypes.c_double(tau)),
        _ctypes.byref(_ctypes.c_double(mad_distance_factor)),
        _ctypes.byref(_ctypes.c_int(max_pooled_residuals)),
        _ctypes.byref(_ctypes.c_int(min_study_neighbors)),
        bin_count_mode,
        _ctypes.byref(_ctypes.c_double(delta_abs_tolerance)),
        divergence_mode,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENTS, _RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    global_js_divergence.flags.writeable = False

    return {
        "global_js_divergence": global_js_divergence,
        "used_k_start": used_k_start.value,
        "used_k_step": used_k_step.value,
        "used_k_max": used_k_max.value,
        "n_points": n_points.value,
    }

def run_js_comp_test_adaptive_permutation(
        expression,
        k_start,
        k_step,
        k_max,
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
        tau=0.1,
        mad_distance_factor=1.0,
        max_pooled_residuals=0,
        min_study_neighbors=1,
        bin_count_mode='bin_count_stone_equal',
        delta_abs_tolerance=0.0001,
        divergence_mode='divergence_pairwise_k2',
):
    r"""Permutation-Test of the adaptive JS-Comp-Test that permutes the samples' study labels

    The adaptive counterpart of
    :func:`tensor_omics.run_js_comp_test_permutation`:
    the same test, each data set run by
    :func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
    with the same data, knobs and options.

    1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
    `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
    shuffled data. Row `i` of `jsd_null` holds its values.
    3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.

    With `rerun_parameter_search = True` every data set, the observed one included, first runs the adaptive parameter search. `random_seed` then
    also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
    draws. The routine is sequential: a caller who wants parallel permutations calls
    :func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
    itself, once per assignment.

    Failures: a permutation that fails with invalid input (a study left without residuals at some reference point, or a construction that does not
    succeed: a status other than ok, no neighborhood, or more neighborhoods than the work arrays hold) gets a NaN row in `jsd_null`, is
    counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
    `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
    returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
    is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own fails the same way with the rerun option.
    Known limitation: a gene missing in a whole study travels with the samples.

    Parameters
    ----------
    expression : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Expression of every sample slot of every gene in every study (see above); NaN for a
        missing value, a missing gene or a padding slot
        NaN is permitted for this value.
    k_start : int
        Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
        The minimum valid value is `1`.
    k_step : int
        Pooled entries staged per adaptive growth round
        The minimum valid value is `1`.
    k_max : int
        Largest number of pooled entries in one neighborhood
        The minimum valid value is `k_start`.
    n_permutations : int
        Number of random assignments (rows of `jsd_null`)
        The minimum valid value is `1`.
    rerun_parameter_search : bool, optional, default False
        `True`: run the adaptive parameter search on the relabelled data first and use its
        knobs; `False`: use `k_start`, `k_step` and `k_max`
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
        Fewest entries of each study every neighborhood must have (construction status; in
        the rerun search, candidate admissibility)
        The minimum valid value is `1`.
        The default value is `1`.
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_stone_equal'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_stone_equal'`.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation_adaptive::run_js_comp_test_adaptive_permutation`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression = _contract.reals(expression, "expression", ndim=3, order="F")
    k_start = _contract.integer(k_start, "k_start")
    k_step = _contract.integer(k_step, "k_step")
    k_max = _contract.integer(k_max, "k_max")
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
    tau = _contract.real(tau, "tau")
    mad_distance_factor = _contract.real(mad_distance_factor, "mad_distance_factor")
    max_pooled_residuals = _contract.integer(max_pooled_residuals, "max_pooled_residuals")
    min_study_neighbors = _contract.integer(min_study_neighbors, "min_study_neighbors")
    bin_count_mode = _np.array([_contract.mode(bin_count_mode, "bin_count_mode", ("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled",)).ljust(29)], dtype="S29")
    delta_abs_tolerance = _contract.real(delta_abs_tolerance, "delta_abs_tolerance")
    divergence_mode = _np.array([_contract.mode(divergence_mode, "divergence_mode", ("divergence_symmetric", "divergence_pairwise_k2",)).ljust(22)], dtype="S22")

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

    _lib.run_js_comp_test_adaptive_permutation_c(
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_genes_all_studies)),
        _ctypes.byref(_ctypes.c_int(max_n_reps_all_studies)),
        expression,
        _ctypes.byref(_ctypes.c_int(k_start)),
        _ctypes.byref(_ctypes.c_int(k_step)),
        _ctypes.byref(_ctypes.c_int(k_max)),
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
        _ctypes.byref(_ctypes.c_double(tau)),
        _ctypes.byref(_ctypes.c_double(mad_distance_factor)),
        _ctypes.byref(_ctypes.c_int(max_pooled_residuals)),
        _ctypes.byref(_ctypes.c_int(min_study_neighbors)),
        bin_count_mode,
        _ctypes.byref(_ctypes.c_double(delta_abs_tolerance)),
        divergence_mode,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _RUN_JS_COMP_TEST_ADAPTIVE_PERMUTATION_ARGUMENTS, _RUN_JS_COMP_TEST_ADAPTIVE_PERMUTATION_ARGUMENT_SOURCES)

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
