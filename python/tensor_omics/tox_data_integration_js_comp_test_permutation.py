r"""tox_data_integration_js_comp_test_permutation

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test

Asks whether the assignment of samples to studies (batches) influences the JS-Comp-Test values
at all. Null hypothesis: it does not. The test shuffles the study labels of the samples
(without replacement), reruns the whole JS-Comp-Test on every shuffled data set, and compares
each study's observed global JSD with the values the shuffled data sets give.

**Input layout.** All studies arrive in one array `expression(max_n_reps_all_studies,
max_n_genes_all_studies, n_studies)`: gene `g` is the same gene in every study, a gene missing in
a study is NaN there, and a replicate slot `(r, s)` whose values are all NaN is padding (a study
with fewer samples than `max_n_reps_all_studies`). Every other slot is a sample. Samples are
numbered in slot order, replicate fastest, then study.

**Entry points.** One assignment of samples to slots, one JS-Comp-Test:
:func:`tensor_omics.run_js_comp_test_on_sample_assignment`
(k-NN neighborhoods, the fixed `(n_points, n_neighbors)` of `run_js_comp_test`) and
:func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
(adaptive neighborhoods, `run_js_comp_test_adaptive`). The whole test, `n_permutations` random
assignments plus the observed one:
:func:`tensor_omics.run_js_comp_test_label_permutation_test`
and
:func:`tensor_omics.run_js_comp_test_adaptive_label_permutation_test`.
`calc_js_comp_test_permutation_bounds` and `calc_adaptive_js_comp_test_permutation_bounds` size
the work arrays for a caller that allocates its own.

**What a p-value means.** Study `i`'s p-value is `(count + 1)/(n_valid + 1)`, where `count` is
how many of the `n_valid` valid shuffled data sets give study `i` a global JSD at least as large
as the observed one (the +1 counts the observed data itself, so it is never 0). A small p-value
(e.g. 0.01) says: hardly any random relabelling makes study `i` look as different from the
others as the real labels do, so its samples differ systematically (a batch effect). A large
one says the real labels are nothing special.

This test is not the residual-level permutation test of `run_js_comp_test`
(:func:`tensor_omics.gjct_permutation_test`), which keeps
the data and the histograms fixed and only reshuffles the pooled histogram counts.

Python binding, generated from tox_data_integration_js_comp_test_permutation. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.run_js_comp_test_on_sample_assignment_c.restype = None
_lib.run_js_comp_test_on_sample_assignment_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
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
_RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "n_points", "n_neighbors", "global_js_divergence", "used_n_points", "used_n_neighbors", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, "sample_assignment", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.run_js_comp_test_adaptive_on_sample_assignment_c.restype = None
_lib.run_js_comp_test_adaptive_on_sample_assignment_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
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
_RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "k_start", "k_step", "k_max", "global_js_divergence", "used_k_start", "used_k_step", "used_k_max", "n_points", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, "sample_assignment", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.run_js_comp_test_label_permutation_test_c.restype = None
_lib.run_js_comp_test_label_permutation_test_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
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
_RUN_JS_COMP_TEST_LABEL_PERMUTATION_TEST_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_points", "n_neighbors", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_LABEL_PERMUTATION_TEST_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, None, None, "jsd_null", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.run_js_comp_test_adaptive_label_permutation_test_c.restype = None
_lib.run_js_comp_test_adaptive_label_permutation_test_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
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
_RUN_JS_COMP_TEST_ADAPTIVE_LABEL_PERMUTATION_TEST_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "k_start", "k_step", "k_max", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_LABEL_PERMUTATION_TEST_ARGUMENT_SOURCES = ("expression", "expression", "expression", None, None, None, None, "jsd_null", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

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
        bin_count_mode='bin_count_occupancy',
):
    r"""Run the k-NN JS-Comp-Test on one assignment of samples to studies

    One building block of the study-label Permutation-Test: relabels the samples by
    `sample_assignment` and runs the whole JS-Comp-Test on the relabelled data set.

    **Input layout.** `expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)`
    holds every study: gene `g` is the same gene in every study, a gene missing in a study is NaN
    there, and a replicate slot whose values are all NaN is padding (a study with fewer samples).
    Every other slot is a sample, numbered in slot order: replicate fastest, then study (sample 1
    is the first non-padding slot of study 1).

    **Assignment.** `sample_assignment` is a permutation of `1..n_samples`: destination sample
    `k` (the `k`-th non-padding slot) receives all genes of source sample `sample_assignment(k)`.
    So every study keeps its slots, its sample count and its padding; only which samples sit in
    them changes. The identity `1..n_samples` is the observed data.

    **Steps.** Per study, the gene means and residuals of the relabelled data
    (:func:`tensor_omics.compute_gene_means`,
    :func:`tensor_omics.compute_residuals`); the
    reference points from the pooled means
    (:func:`tensor_omics.pool_means` on their
    ascending, NaN-last order); then
    :func:`tensor_omics.run_js_comp_test` with
    `n_permutations = 0` (its own residual-level permutation test does not run) and the bin-count
    options as given. On the identity assignment the result is exactly that of a caller doing
    these steps on the observed data with `run_js_comp_test`.

    **Rerunning the parameter search.** With `rerun_parameter_search = True`, the relabelled
    data set first runs
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    with the search options below, and the run then uses the `(n_points, n_neighbors)` it returns
    (also its fallback when no candidate was admissible); `n_points` and `n_neighbors` are then
    not used for the run, only to size the work arrays. `used_n_points` and `used_n_neighbors`
    report what the run used. Run time: every call runs the whole JS-Comp-Test; with the rerun
    option it runs the whole parameter search first.

    **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
    the samples, so after relabelling that study may have values for it.

    Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
    no sample at all, when `n_samples` is not the number of samples in `expression`, or when
    `sample_assignment` is not a permutation of `1..n_samples`. (Fortran expert tier only, where
    the caller supplies the work-array capacities: with the rerun option it also fails when a
    capacity is below the search's own bound, `calc_js_comp_test_candidate_bounds`.) Errors of
    the search and the run are returned with their code (for example `ERR_INVALID_INPUT` when some
    study has no residual at some reference point of the relabelled data). `global_js_divergence`
    is then NaN.

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
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_occupancy'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_occupancy'`.

    Returns
    -------
    dict
        with keys:

        global_js_divergence : np.ndarray[np.float64] of shape (n_studies,), read-only
            Weighted global JSD of each study against the consensus, for this assignment; NaN on error
            A result is a value; call `.copy()` to obtain a modifiable array.
        used_n_points : int
            Number of reference points the run used (`n_points`, or the search's choice)
        used_n_neighbors : int
            Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation::run_js_comp_test_on_sample_assignment`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression = np.asfortranarray(expression, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression' must be an array of np.float64: {error}") from None
    if expression.ndim != 3:
        raise ValueError(f"'expression' must have 3 dimensions, but has {expression.ndim}")
    try:
        sample_assignment = np.ascontiguousarray(sample_assignment, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sample_assignment' must be an array of np.int32: {error}") from None
    if sample_assignment.ndim != 1:
        raise ValueError(f"'sample_assignment' must have 1 dimension, but has {sample_assignment.ndim}")
    join_method = np.array([str(join_method).lower().encode().ljust(11)], dtype="S11")
    plateau_mode = np.array([str(plateau_mode).lower().encode().ljust(19)], dtype="S19")
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

    # what the inputs already say, rather than asking for it again
    n_studies = expression.shape[2]
    max_n_genes_all_studies = expression.shape[1]
    max_n_reps_all_studies = expression.shape[0]
    n_samples = sample_assignment.shape[0]

    # outputs and work arrays, which the caller never sees
    global_js_divergence = np.empty((n_studies,), dtype=np.float64, order='C')
    used_n_points = ctypes.c_int(0)
    used_n_neighbors = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_on_sample_assignment_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        expression,
        ctypes.byref(ctypes.c_int(n_samples)),
        sample_assignment,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        global_js_divergence,
        ctypes.byref(used_n_points),
        ctypes.byref(used_n_neighbors),
        ctypes.byref(ctypes.c_bool(rerun_parameter_search)),
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        join_method,
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

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENTS, _RUN_JS_COMP_TEST_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    global_js_divergence.flags.writeable = False

    return {
        "global_js_divergence": global_js_divergence,
        "used_n_points": used_n_points.value,
        "used_n_neighbors": used_n_neighbors.value,
    }

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
        bin_count_mode='bin_count_occupancy',
):
    r"""Run the adaptive JS-Comp-Test on one assignment of samples to studies

    The adaptive counterpart of
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`,
    with the same input layout and assignment convention: `expression` holds every study (NaN for
    a missing gene, all-NaN replicate slots are padding, samples numbered replicate fastest, then
    study), and destination sample `k` receives all genes of source sample `sample_assignment(k)`;
    the identity is the observed data.

    **Steps.** Per study, the gene means and residuals of the relabelled data; the adaptive
    neighborhoods for the growth knobs `k_start`, `k_step`, `k_max` and the construction options
    (:func:`tensor_omics.construct_adaptive_neighborhoods`);
    then :func:`tensor_omics.run_js_comp_test_adaptive`
    on its ranges with `n_permutations = 0`. `n_points` is the number of neighborhoods the
    construction built. On the identity assignment the result is exactly that of a caller doing
    these steps on the observed data.

    **Rerunning the parameter search.** With `rerun_parameter_search = True`, the relabelled
    data first runs
    :func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`
    (with its own capacity, `calc_adaptive_js_comp_test_bounds`, so it rejects the same candidates
    as a direct call), and the construction and run then use the knobs it returns (also its
    fallback when no candidate was admissible); `k_start`, `k_step` and `k_max` are then not used,
    except `k_start` to size the work arrays. `used_k_start`, `used_k_step` and `used_k_max` report
    what was used. Run time: every call runs the whole JS-Comp-Test; with the rerun option it runs
    the whole parameter search first.

    **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
    the samples, so after relabelling that study may have values for it.

    Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
    no sample at all, when `n_samples` is not the number of samples in `expression`, when
    `sample_assignment` is not a permutation of `1..n_samples`, or when the construction builds no
    neighborhood or more neighborhoods than the work arrays hold. (Fortran expert tier only, where
    the caller supplies the work-array capacity: with the rerun option it also fails when the
    capacity is below the search's own.) Errors of the search and the run are returned with their
    code. `global_js_divergence` is then NaN.

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
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_occupancy'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_occupancy'`.

    Returns
    -------
    dict
        with keys:

        global_js_divergence : np.ndarray[np.float64] of shape (n_studies,), read-only
            Weighted global JSD of each study against the consensus, for this assignment; NaN on error
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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation::run_js_comp_test_adaptive_on_sample_assignment`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression = np.asfortranarray(expression, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression' must be an array of np.float64: {error}") from None
    if expression.ndim != 3:
        raise ValueError(f"'expression' must have 3 dimensions, but has {expression.ndim}")
    try:
        sample_assignment = np.ascontiguousarray(sample_assignment, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sample_assignment' must be an array of np.int32: {error}") from None
    if sample_assignment.ndim != 1:
        raise ValueError(f"'sample_assignment' must have 1 dimension, but has {sample_assignment.ndim}")
    join_method = np.array([str(join_method).lower().encode().ljust(11)], dtype="S11")
    plateau_mode = np.array([str(plateau_mode).lower().encode().ljust(19)], dtype="S19")
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

    # what the inputs already say, rather than asking for it again
    n_studies = expression.shape[2]
    max_n_genes_all_studies = expression.shape[1]
    max_n_reps_all_studies = expression.shape[0]
    n_samples = sample_assignment.shape[0]

    # outputs and work arrays, which the caller never sees
    global_js_divergence = np.empty((n_studies,), dtype=np.float64, order='C')
    used_k_start = ctypes.c_int(0)
    used_k_step = ctypes.c_int(0)
    used_k_max = ctypes.c_int(0)
    n_points = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_adaptive_on_sample_assignment_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        expression,
        ctypes.byref(ctypes.c_int(n_samples)),
        sample_assignment,
        ctypes.byref(ctypes.c_int(k_start)),
        ctypes.byref(ctypes.c_int(k_step)),
        ctypes.byref(ctypes.c_int(k_max)),
        global_js_divergence,
        ctypes.byref(used_k_start),
        ctypes.byref(used_k_step),
        ctypes.byref(used_k_max),
        ctypes.byref(n_points),
        ctypes.byref(ctypes.c_bool(rerun_parameter_search)),
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        join_method,
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

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENTS, _RUN_JS_COMP_TEST_ADAPTIVE_ON_SAMPLE_ASSIGNMENT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    global_js_divergence.flags.writeable = False

    return {
        "global_js_divergence": global_js_divergence,
        "used_k_start": used_k_start.value,
        "used_k_step": used_k_step.value,
        "used_k_max": used_k_max.value,
        "n_points": n_points.value,
    }

def run_js_comp_test_label_permutation_test(
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
        bin_count_mode='bin_count_occupancy',
):
    r"""Study-label Permutation-Test of the k-NN JS-Comp-Test

    Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value
    (null hypothesis: it does not). The same data, parameters and options as
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`,
    whose input layout (all studies in `expression`, NaN for a missing gene, all-NaN replicate
    slots are padding) and assignment convention apply here too.

    1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    2. `n_permutations` times, the study labels of the samples are shuffled without replacement (a
    uniformly random permutation, unbiased Fisher-Yates, from one random number generator seeded
    with `random_seed`); every study keeps its sample count. The whole JS-Comp-Test runs on the
    shuffled data; row `i` of `jsd_null` holds its values.
    3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`,
    the Laplace correction of the residual-level test, so it is never 0 and at least
    `1/(n_permutations + 1)`.

    **Rerunning the parameter search** (`rerun_parameter_search = True`): every data set, the
    observed one included, first runs the fixed-k parameter search and then the run with the
    `(n_points, n_neighbors)` it selects; `n_points` and `n_neighbors` then only size the work
    arrays. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
    runs the whole parameter search first. The routine is sequential: a caller who wants
    permutations in parallel calls
    :func:`tensor_omics.run_js_comp_test_on_sample_assignment`
    itself with its own assignments (e.g. from R), one call per assignment.

    **Failures.** A permutation whose JS-Comp-Test fails with invalid input (for example a study
    left without residuals at some reference point) gets a NaN row in `jsd_null`, is counted in
    `n_failed_permutations` and does not count in `n_valid`. No valid permutation gives NaN
    p-values (a direct call with `n_permutations = 0` too; the published entry points reject 0,
    since it is the extent of `jsd_null`). Any other error ends the test and is returned, with NaN
    p-values; the rows of the permutations run so far are kept. An error on the observed data is
    returned at once (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study
    without any sample (the call fails with invalid input). In the Fortran expert tier, where the
    caller supplies the work-array capacities, a capacity below the search's own bound fails the
    same way with the rerun option.

    **Known limitation.** A gene missing in a whole study travels with the samples, so after
    shuffling that study may have values for it.

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
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_occupancy'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_occupancy'`.

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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation::run_js_comp_test_label_permutation_test`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression = np.asfortranarray(expression, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression' must be an array of np.float64: {error}") from None
    if expression.ndim != 3:
        raise ValueError(f"'expression' must have 3 dimensions, but has {expression.ndim}")
    join_method = np.array([str(join_method).lower().encode().ljust(11)], dtype="S11")
    plateau_mode = np.array([str(plateau_mode).lower().encode().ljust(19)], dtype="S19")
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

    # what the inputs already say, rather than asking for it again
    n_studies = expression.shape[2]
    max_n_genes_all_studies = expression.shape[1]
    max_n_reps_all_studies = expression.shape[0]

    # outputs and work arrays, which the caller never sees
    global_js_divergence_observed = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values = np.empty((n_studies,), dtype=np.float64, order='C')
    jsd_null = np.empty((n_permutations, n_studies,), dtype=np.float64, order='F')
    n_failed_permutations = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_label_permutation_test_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        expression,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_permutations)),
        global_js_divergence_observed,
        p_values,
        jsd_null,
        ctypes.byref(n_failed_permutations),
        ctypes.byref(ctypes.c_bool(rerun_parameter_search)),
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        join_method,
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

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_LABEL_PERMUTATION_TEST_ARGUMENTS, _RUN_JS_COMP_TEST_LABEL_PERMUTATION_TEST_ARGUMENT_SOURCES)

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

def run_js_comp_test_adaptive_label_permutation_test(
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
        bin_count_mode='bin_count_occupancy',
):
    r"""Study-label Permutation-Test of the adaptive JS-Comp-Test

    The adaptive counterpart of
    :func:`tensor_omics.run_js_comp_test_label_permutation_test`:
    the same test (observed value, `n_permutations` shuffles of the study labels without
    replacement, p-value `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`), each
    data set run by
    :func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
    with the same data, knobs and options; its input layout and assignment convention apply here too.

    **Rerunning the parameter search** (`rerun_parameter_search = True`): every data set first
    runs the adaptive parameter search and then the construction and run with the knobs it
    selects. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
    runs the whole parameter search first. The routine is sequential: a caller who wants
    permutations in parallel calls
    :func:`tensor_omics.run_js_comp_test_adaptive_on_sample_assignment`
    itself with its own assignments (e.g. from R), one call per assignment.

    **Failures.** A permutation that fails with invalid input (a study left without residuals at
    some reference point, or a construction that builds no neighborhood or more neighborhoods
    than the work arrays hold) gets a NaN row in `jsd_null`, is counted in `n_failed_permutations`
    and does not count in `n_valid`. No valid permutation gives NaN p-values (a direct call with
    `n_permutations = 0` too; the published entry points reject 0, since it is the extent of
    `jsd_null`). Any other error ends the test and is returned, with NaN p-values; the rows of the
    permutations run so far are kept. An error on the observed data is returned at once
    (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study without any
    sample (the call fails with invalid input). In the Fortran expert tier, where the caller
    supplies the work-array capacity, a capacity below the search's own fails the same way with
    the rerun option.

    **Known limitation.** A gene missing in a whole study travels with the samples, so after
    shuffling that study may have values for it.

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
    bin_count_mode : str, one of 'bin_count_occupancy' | 'bin_count_stone_equal' | 'bin_count_stone_size_weighted' | 'bin_count_stone_pooled', optional, default 'bin_count_occupancy'
        How each reference point's histogram bin count is chosen; the Stone modes test every
        bin count from `m_min` to `m_max`
        (:func:`tensor_omics.determine_bin_count_stone`),
        with `min_residuals_per_bin` as the floor

        The default value is `'bin_count_occupancy'`.

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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_permutation::run_js_comp_test_adaptive_label_permutation_test`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression = np.asfortranarray(expression, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression' must be an array of np.float64: {error}") from None
    if expression.ndim != 3:
        raise ValueError(f"'expression' must have 3 dimensions, but has {expression.ndim}")
    join_method = np.array([str(join_method).lower().encode().ljust(11)], dtype="S11")
    plateau_mode = np.array([str(plateau_mode).lower().encode().ljust(19)], dtype="S19")
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

    # what the inputs already say, rather than asking for it again
    n_studies = expression.shape[2]
    max_n_genes_all_studies = expression.shape[1]
    max_n_reps_all_studies = expression.shape[0]

    # outputs and work arrays, which the caller never sees
    global_js_divergence_observed = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values = np.empty((n_studies,), dtype=np.float64, order='C')
    jsd_null = np.empty((n_permutations, n_studies,), dtype=np.float64, order='F')
    n_failed_permutations = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_adaptive_label_permutation_test_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        expression,
        ctypes.byref(ctypes.c_int(k_start)),
        ctypes.byref(ctypes.c_int(k_step)),
        ctypes.byref(ctypes.c_int(k_max)),
        ctypes.byref(ctypes.c_int(n_permutations)),
        global_js_divergence_observed,
        p_values,
        jsd_null,
        ctypes.byref(n_failed_permutations),
        ctypes.byref(ctypes.c_bool(rerun_parameter_search)),
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        join_method,
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

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_ADAPTIVE_LABEL_PERMUTATION_TEST_ARGUMENTS, _RUN_JS_COMP_TEST_ADAPTIVE_LABEL_PERMUTATION_TEST_ARGUMENT_SOURCES)

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
