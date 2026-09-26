r"""tox_data_integration_js_comp_test

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
allocates its own. For adaptive (heteroscedastic, Issue #217) neighborhood construction,
:func:`tensor_omics.generate_adaptive_js_comp_test_candidates`
generates the ascending `(k_start, k_step, k_max)` growth-knob candidates instead, and
`calc_adaptive_js_comp_test_bounds` recommends the reference-point capacity their search's
per-point arrays are sized by;
:func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`
is that search, the adaptive counterpart of the fixed-k parameter search, with the same
plateau criteria. Once a candidate has passed both gates, its bootstrap confidence interval
is resampled from the pooled consensus histogram by
:func:`tensor_omics.bootstrap_histogram` (heap
size recommended by
:func:`tensor_omics.calc_js_comp_test_n_top_k_jsds`).
Later stages of the port add the K-study permutation test and the top-level orchestrator that
wire these building blocks together.

Python binding, generated from tox_data_integration_js_comp_test. Do not edit.
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

_lib.generate_js_comp_test_candidates_c.restype = None
_lib.generate_js_comp_test_candidates_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GENERATE_JS_COMP_TEST_CANDIDATES_ARGUMENTS = ("max_n_genes_all_studies", "candidates_n_points_n_neighbors", "n_candidates", "ierr",)

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

_lib.construct_adaptive_neighborhoods_c.restype = None
_lib.construct_adaptive_neighborhoods_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CONSTRUCT_ADAPTIVE_NEIGHBORHOODS_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "k_start", "k_step", "k_max", "n_points", "x_star", "pooled_neighborhood_range", "n_neighbors_per_point", "stop_reason", "neighborhood_dispersion", "neighborhood_mad", "max_n_neighbors", "construction_status", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "ierr",)
#: For a derived argument, the one the caller passed it in
_CONSTRUCT_ADAPTIVE_NEIGHBORHOODS_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.check_neighborhood_overlaps_c.restype = None
_lib.check_neighborhood_overlaps_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENTS = ("neighborhood_range", "n_points", "min_neighbor_overlap", "all_have_min_neighbor_overlap", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENT_SOURCES = (None, "neighborhood_range", None, None, None,)

_lib.check_mean_pmf_min_counts_c.restype = None
_lib.check_mean_pmf_min_counts_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENTS = ("mean_pmf_counts", "n_bins", "n_bins_per_point", "n_points", "min_count", "all_bins_have_min_count", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENT_SOURCES = (None, "mean_pmf_counts", None, "mean_pmf_counts", None, None, None,)

_lib.check_plateau_condition_c.restype = None
_lib.check_plateau_condition_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_PLATEAU_CONDITION_ARGUMENTS = ("confidence_interval", "best_candidate_pair_confidence_interval", "n_studies", "best_candidate_index", "best_exceeded_ci_overlap_count", "candidate_index", "join_method", "succeeding_ci_overlap", "plateau_found", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_PLATEAU_CONDITION_ARGUMENT_SOURCES = (None, None, "confidence_interval", None, None, None, None, None, None, None,)

_lib.check_effect_size_plateau_condition_c.restype = None
_lib.check_effect_size_plateau_condition_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENTS = ("global_js_divergence", "prev_global_js_divergence", "n_studies", "has_previous", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "n_consecutive_ok", "delta", "delta_median", "delta_max", "plateau_found", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENT_SOURCES = (None, None, "global_js_divergence", None, None, None, None, None, None, None, None, None, None, None,)

_lib.create_mean_pmf_c.restype = None
_lib.create_mean_pmf_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CREATE_MEAN_PMF_ARGUMENTS = ("pmfs", "counts", "n_bins", "n_points", "n_studies", "included_n_reps", "mean_pmf", "mean_pmf_included_n_reps", "mean_pmf_counts", "ierr",)
#: For a derived argument, the one the caller passed it in
_CREATE_MEAN_PMF_ARGUMENT_SOURCES = (None, None, "pmfs", "pmfs", "pmfs", None, None, None, None, None,)

_lib.create_mean_pmf_only_c.restype = None
_lib.create_mean_pmf_only_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CREATE_MEAN_PMF_ONLY_ARGUMENTS = ("pmfs", "n_bins", "n_points", "n_studies", "mean_pmf", "ierr",)
#: For a derived argument, the one the caller passed it in
_CREATE_MEAN_PMF_ONLY_ARGUMENT_SOURCES = (None, "pmfs", "pmfs", "pmfs", None, None,)

_lib.bootstrap_histogram_c.restype = None
_lib.bootstrap_histogram_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BOOTSTRAP_HISTOGRAM_ARGUMENTS = ("n_bootstraps", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf_included_n_reps", "included_n_reps", "confidence_interval", "two_sided_bootstrapping_significance_level", "random_seed", "ierr",)
#: For a derived argument, the one the caller passed it in
_BOOTSTRAP_HISTOGRAM_ARGUMENT_SOURCES = (None, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", None, None, None, None, None, None, None,)

_lib.run_js_comp_test_c.restype = None
_lib.run_js_comp_test_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
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
_RUN_JS_COMP_TEST_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "gene_means_perms", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.run_js_comp_test_adaptive_c.restype = None
_lib.run_js_comp_test_adaptive_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
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
_RUN_JS_COMP_TEST_ADAPTIVE_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "gene_means", "residuals", "pooled_neighborhood_range", "n_neighbors_per_point", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", "pooled_neighborhood_range", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

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
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "max_n_neighbors_candidate", "n_points", "n_neighbors", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_pair_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_PARAMETER_SEARCH_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", None, None, None, None, "n_bins_per_point", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

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
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_ADAPTIVE_PARAMETER_SEARCH_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "k_start", "k_step", "k_max", "n_points", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_k_start", "trace_k_step", "trace_k_max", "trace_n_points", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "n_candidates_tried", "candidate_k_start", "candidate_n_points", "candidate_status", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_PARAMETER_SEARCH_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", None, None, None, None, "n_bins_per_point", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def estimate_bin_count(
        residuals,
        max_n_reps_all_studies,
        n_neighbors,
        shared_residual_range,
):
    r"""Estimate the histogram bin count for one (n_points, n_neighbors) candidate

    Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
    the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
    is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
    range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
    independently, each clamped on its own to at most
    ``MAX_N_BINS`` bins and returned as
    `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
    interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
    skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.

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
            to the (clamped) sturges_bins when the interquartile range is too close to zero to
            divide by (see the is_close guard below)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::estimate_bin_count`, whose argument names are
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

    Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
    the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
    is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
    range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
    independently, each clamped on its own to at most
    ``MAX_N_BINS`` bins and returned as
    `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
    interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
    skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.

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
            to the (clamped) sturges_bins when the interquartile range is too close to zero to
            divide by (see the is_close guard below)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::estimate_bin_count_expert`, whose argument names are
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

    Implements Issue #187's two-stage geometric-search-then-local-refinement algorithm for one
    neighborhood's pooled residuals (`pooled_residuals`, across all its neighbors and all
    studies): find the largest bin count `M` in `[m_min, m_max]` whose equal-width histogram
    over `[shared_residual_range_low, shared_residual_range_high]` -- this neighborhood's own
    asymmetric range, the `lower_residual_range_quantile`/`upper_residual_range_quantile`
    percentiles of its own pooled signed residuals, rather than a single dataset-wide symmetric
    range -- has every bin at or above `min_residuals_per_bin` (the occupancy criterion), rather
    than the generic
    Sturges/Freedman-Diaconis rule
    :func:`tensor_omics.estimate_bin_count` alone
    applies, which is why that routine is still called here too -- purely for the
    `sturges_bins`/`fd_bins` diagnostic outputs, never for the decision itself.

    Stage 1 grows a candidate `trial_m` geometrically from `m_min` (`next_m =
    ceiling(gamma_occupancy*trial_m)`, guaranteed to advance by at least 1 via
    `max(trial_m + 1, next_m)`, clamped to `m_max`) until occupancy first fails, recording the
    largest admissible `m_valid` and the first inadmissible `m_invalid`; reaching `m_max` while
    still admissible returns it immediately, skipping stage 2 entirely. Stage 2 then tests every
    integer strictly between `m_valid` and `m_invalid` -- not a binary search, since equal-width
    bin boundaries are recomputed for every candidate `M` and occupancy is therefore not
    guaranteed monotonic in `M` (Issue #187 is explicit about this) -- and keeps the largest one
    that still passes. Both stages reuse the occupancy diagnostics (`min`/`max_bin_occupancy`)
    computed for whichever `M` ends up selected, rather than recomputing them a second time
    afterward; `mean_bin_occupancy` needs no such bookkeeping, since every non-NaN pooled
    residual lands in exactly one bin at any `M`, so it is always exactly
    `n_pooled_residuals / selected_n_bins`.

    Per the issue's own FAILURE policy, `occupancy_failed = True` (even `m_min` bins could not
    satisfy the occupancy criterion, or every pooled residual is NaN) means the neighborhood
    should be rejected by the caller rather than built from `selected_n_bins` -- which is still
    set to `m_min` in that case, purely so a caller ignoring `occupancy_failed` has *a* value to
    build with, never as an indication the search actually found `m_min` admissible.

    The outer geometric search is a genuine sequential `do`/`exit` state machine, not `do
    concurrent`: `m_valid`/`m_invalid` accumulate across iterations and each iteration's
    continuation depends on the previous one's outcome, exactly the "loops with data-dependent
    control flow across iterations" case Fortran_Coding_Guides.pdf Sec 10 carves out as the
    deliberate exception to `do concurrent`.

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
        Largest candidate bin count the search will ever test (M_max); if a caller passes
        `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
        relying on an unconfirmed generator capability to bound one optional argument by
        another
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
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for this neighborhood's own upper residual-range bound
        (shared_residual_range_high)
        The minimum valid value is `0.0`.
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
            `True` iff even m_min bins could not satisfy the occupancy criterion (including
            the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
            caller should reject this neighborhood rather than build a histogram from
            selected_n_bins
        shared_residual_range_low : float
            This neighborhood's own lower residual-range bound (R_low): the
            lower_residual_range_quantile percentile of its own pooled signed residuals --
            replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            occupancy_failed because every pooled residual is NaN
        shared_residual_range_high : float
            This neighborhood's own upper residual-range bound (R_high): the
            upper_residual_range_quantile percentile of its own pooled signed residuals --
            replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            occupancy_failed because every pooled residual is NaN
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
            Sturges' rule estimate for this neighborhood's own pooled residuals, from
            estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        fd_bins : int
            Freedman-Diaconis rule estimate for this neighborhood's own pooled residuals, from
            estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::determine_bin_count_occupancy`, whose argument names are
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

    Implements Issue #187's two-stage geometric-search-then-local-refinement algorithm for one
    neighborhood's pooled residuals (`pooled_residuals`, across all its neighbors and all
    studies): find the largest bin count `M` in `[m_min, m_max]` whose equal-width histogram
    over `[shared_residual_range_low, shared_residual_range_high]` -- this neighborhood's own
    asymmetric range, the `lower_residual_range_quantile`/`upper_residual_range_quantile`
    percentiles of its own pooled signed residuals, rather than a single dataset-wide symmetric
    range -- has every bin at or above `min_residuals_per_bin` (the occupancy criterion), rather
    than the generic
    Sturges/Freedman-Diaconis rule
    :func:`tensor_omics.estimate_bin_count` alone
    applies, which is why that routine is still called here too -- purely for the
    `sturges_bins`/`fd_bins` diagnostic outputs, never for the decision itself.

    Stage 1 grows a candidate `trial_m` geometrically from `m_min` (`next_m =
    ceiling(gamma_occupancy*trial_m)`, guaranteed to advance by at least 1 via
    `max(trial_m + 1, next_m)`, clamped to `m_max`) until occupancy first fails, recording the
    largest admissible `m_valid` and the first inadmissible `m_invalid`; reaching `m_max` while
    still admissible returns it immediately, skipping stage 2 entirely. Stage 2 then tests every
    integer strictly between `m_valid` and `m_invalid` -- not a binary search, since equal-width
    bin boundaries are recomputed for every candidate `M` and occupancy is therefore not
    guaranteed monotonic in `M` (Issue #187 is explicit about this) -- and keeps the largest one
    that still passes. Both stages reuse the occupancy diagnostics (`min`/`max_bin_occupancy`)
    computed for whichever `M` ends up selected, rather than recomputing them a second time
    afterward; `mean_bin_occupancy` needs no such bookkeeping, since every non-NaN pooled
    residual lands in exactly one bin at any `M`, so it is always exactly
    `n_pooled_residuals / selected_n_bins`.

    Per the issue's own FAILURE policy, `occupancy_failed = True` (even `m_min` bins could not
    satisfy the occupancy criterion, or every pooled residual is NaN) means the neighborhood
    should be rejected by the caller rather than built from `selected_n_bins` -- which is still
    set to `m_min` in that case, purely so a caller ignoring `occupancy_failed` has *a* value to
    build with, never as an indication the search actually found `m_min` admissible.

    The outer geometric search is a genuine sequential `do`/`exit` state machine, not `do
    concurrent`: `m_valid`/`m_invalid` accumulate across iterations and each iteration's
    continuation depends on the previous one's outcome, exactly the "loops with data-dependent
    control flow across iterations" case Fortran_Coding_Guides.pdf Sec 10 carves out as the
    deliberate exception to `do concurrent`.

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
        Largest candidate bin count the search will ever test (M_max); if a caller passes
        `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
        relying on an unconfirmed generator capability to bound one optional argument by
        another
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
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for this neighborhood's own upper residual-range bound
        (shared_residual_range_high)
        The minimum valid value is `0.0`.
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
            `True` iff even m_min bins could not satisfy the occupancy criterion (including
            the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
            caller should reject this neighborhood rather than build a histogram from
            selected_n_bins
        shared_residual_range_low : float
            This neighborhood's own lower residual-range bound (R_low): the
            lower_residual_range_quantile percentile of its own pooled signed residuals --
            replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            occupancy_failed because every pooled residual is NaN
        shared_residual_range_high : float
            This neighborhood's own upper residual-range bound (R_high): the
            upper_residual_range_quantile percentile of its own pooled signed residuals --
            replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            occupancy_failed because every pooled residual is NaN
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
            Sturges' rule estimate for this neighborhood's own pooled residuals, from
            estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        fd_bins : int
            Freedman-Diaconis rule estimate for this neighborhood's own pooled residuals, from
            estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::determine_bin_count_occupancy_expert`, whose argument names are
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

    Tests every candidate bin count `M` in `[m_min, m_max]` independently and keeps the largest
    one whose pooled histogram satisfies the occupancy criterion, instead of
    :func:`tensor_omics.determine_bin_count_occupancy`'s
    own fast geometric-search-then-refinement. Exists purely to validate that routine's result:
    occupancy is not guaranteed monotonic in `M` once bin boundaries are recomputed per candidate
    (Issue #187 is explicit about this), so a search that stops at the first failure can in
    principle miss a larger, independently-admissible `M` the geometric ladder never tries. Takes
    `shared_residual_range_low`/`shared_residual_range_high`/`n_pooled_residuals` as direct
    inputs, already produced by a prior
    :func:`tensor_omics.determine_bin_count_occupancy`
    call, rather than re-deriving them -- this isolates the comparison to just the `M`-selection
    algorithm, uncontaminated by a second independent percentile computation.

    Every candidate is independent (no early exit, no state carried between iterations, unlike
    the production routine's own Stage 1/Stage 2), so this is a genuine `do concurrent` with
    `reduce(max:...)`, not a sequential search: `candidate_bin_counts` is declared local to the
    loop (an ordinary MAX_N_BINS-sized local, not a `tmp_` dummy -- precedented by
    `run_js_comp_test_impl`'s own `candidates_n_points_n_neighbors` local) so each concurrent
    iteration gets its own private scratch instead of racing on a shared buffer sliced by
    `trial_m`. `reduce(max:...)` has no "argmax" form, so the winning `M`'s own
    min/mean/max_bin_occupancy are recovered with one extra, ordinary (non-concurrent) call to
    `histogram_bin_counts` for `best_m` alone once the reduction is done -- trivial cost next to
    the search itself.

    `n_pooled_residuals == 0` and a degenerate zero-width range need no special-case branch here:
    `histogram_bin_counts` already guards the degenerate range internally, and an all-zero pool
    naturally resolves to `occupancy_failed` (or a trivial pass at `min_residuals_per_bin=0`,
    with `mean_bin_occupancy=0.0`, no divide-by-zero) -- intentional, not an oversight.

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
        Largest candidate bin count tested (M_max); if a caller passes `m_max < m_min`, the
        implementation clamps it up to `m_min` internally rather than relying on an
        unconfirmed generator capability to bound one optional argument by another
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
            `True` iff no candidate bin count in [m_min, m_max] satisfies the occupancy
            criterion (including the case where n_pooled_residuals is 0 and min_residuals_per_bin
            is not itself 0)
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::determine_bin_count_occupancy_exhaustive`, whose argument names are
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

    Tests every candidate bin count `M` in `[m_min, m_max]` independently and keeps the largest
    one whose pooled histogram satisfies the occupancy criterion, instead of
    :func:`tensor_omics.determine_bin_count_occupancy`'s
    own fast geometric-search-then-refinement. Exists purely to validate that routine's result:
    occupancy is not guaranteed monotonic in `M` once bin boundaries are recomputed per candidate
    (Issue #187 is explicit about this), so a search that stops at the first failure can in
    principle miss a larger, independently-admissible `M` the geometric ladder never tries. Takes
    `shared_residual_range_low`/`shared_residual_range_high`/`n_pooled_residuals` as direct
    inputs, already produced by a prior
    :func:`tensor_omics.determine_bin_count_occupancy`
    call, rather than re-deriving them -- this isolates the comparison to just the `M`-selection
    algorithm, uncontaminated by a second independent percentile computation.

    Every candidate is independent (no early exit, no state carried between iterations, unlike
    the production routine's own Stage 1/Stage 2), so this is a genuine `do concurrent` with
    `reduce(max:...)`, not a sequential search: `candidate_bin_counts` is declared local to the
    loop (an ordinary MAX_N_BINS-sized local, not a `tmp_` dummy -- precedented by
    `run_js_comp_test_impl`'s own `candidates_n_points_n_neighbors` local) so each concurrent
    iteration gets its own private scratch instead of racing on a shared buffer sliced by
    `trial_m`. `reduce(max:...)` has no "argmax" form, so the winning `M`'s own
    min/mean/max_bin_occupancy are recovered with one extra, ordinary (non-concurrent) call to
    `histogram_bin_counts` for `best_m` alone once the reduction is done -- trivial cost next to
    the search itself.

    `n_pooled_residuals == 0` and a degenerate zero-width range need no special-case branch here:
    `histogram_bin_counts` already guards the degenerate range internally, and an all-zero pool
    naturally resolves to `occupancy_failed` (or a trivial pass at `min_residuals_per_bin=0`,
    with `mean_bin_occupancy=0.0`, no divide-by-zero) -- intentional, not an oversight.

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
        Largest candidate bin count tested (M_max); if a caller passes `m_max < m_min`, the
        implementation clamps it up to `m_min` internally rather than relying on an
        unconfirmed generator capability to bound one optional argument by another
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
            `True` iff no candidate bin count in [m_min, m_max] satisfies the occupancy
            criterion (including the case where n_pooled_residuals is 0 and min_residuals_per_bin
            is not itself 0)
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::determine_bin_count_occupancy_exhaustive_expert`, whose argument names are
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

def generate_js_comp_test_candidates(
        max_n_genes_all_studies,
):
    r"""Generate the GAMMA-decay (n_points, n_neighbors) candidate grid

    Ported from the grid-building half of 125-stabilize-jscomp's
    `determine_js_comp_test_n_points_n_neighbors_helper`: starting from an initial
    `n_points_high` (clamped between MIN_POINTS and MAX_POINTS), repeatedly multiplies by
    GAMMA until it would drop below `n_points_low`, and for each distinct resulting
    `n_points` candidate pairs it with up to `size(KX_FACTORS)` distinct `n_neighbors`
    candidates. A duplicate `n_points` or `n_neighbors` value (from clamping or
    integer rounding) collapses rather than repeating -- this is real, derived behavior the
    grid depends on to avoid redundant candidates at small `max_n_genes_all_studies`, not a
    bug: a small enough `max_n_genes_all_studies` collapses the whole grid down to exactly one
    candidate.

    Issue #187: this routine no longer estimates a per-candidate histogram bin count as a side
    effect -- both
    :func:`tensor_omics.run_js_comp_test` and
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    now compute real per-neighborhood bin counts via
    :func:`tensor_omics.determine_bin_count_occupancy`
    once a candidate has passed admissibility, superseding the old global-pool
    `estimate_bin_count_impl` estimate this routine used to produce for every candidate
    regardless of admissibility.

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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::generate_js_comp_test_candidates`, whose argument names are
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

def generate_adaptive_js_comp_test_candidates(
        max_n_genes_all_studies,
        n_studies,
):
    r"""Generate the ascending adaptive (k_start, k_step, k_max) candidate sequence

    The adaptive-neighborhood counterpart (Issue #217) of
    :func:`tensor_omics.generate_js_comp_test_candidates`.
    An adaptive candidate does not fix its number of reference points: each neighborhood is
    grown from `k_start` pooled entries in rounds of `k_step` up to `k_max` while its residual
    dispersion stays stable, and the next reference point is placed beyond it, so the point
    count emerges from the growth knobs. This routine only generates those knobs. Because
    `k_start` shrinks from one candidate to the next, the candidates **ascend** in the number of
    reference points they are expected to produce.

    With `N = max_n_genes_all_studies * n_studies` (the padded pool size, not the count of
    non-NaN means, so that the capacity producer
    :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`
    and the search agree on the same number) and `floor_k = max(10, 2*n_studies)`:

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
    With them, one decay step lowers any `k_start >= 10` by at least 2, so the values are
    strictly descending without any deduplication, every `k_step` is at least 1 and every
    `k_max` at least `k_start`.

    Each `k_start_t` is computed in double precision directly from `k_start_1`, not by
    rounding the previous candidate again, so rounding never accumulates along the sequence.

    Every knob scales with `N`, so the sequence does not collapse for small data the way the
    fixed-k grid does below about 8,743 genes, where its absolute `n_points_low` floor leaves
    room for only a single `n_points` value. The first candidate is always accepted, so even a
    data set smaller than `k_start_1` (for example 3 genes in 1 study, where `k_start_1 = 10`)
    yields exactly one candidate; growing its neighborhoods then fails with a too-few-means
    status, which is the adaptive search's to report, not this routine's.

    `N` is formed in 64-bit integer arithmetic and never overflows. The generated values are
    representable as long as `k_start_1` does not exceed `huge(1)/4 = 536870911` (so
    that `k_max_1` fits a 32-bit integer), i.e. for `N` up to about 2.7e10 and `n_studies` up
    to about 2.7e8 -- far beyond any real data set. This routine does not check that bound;
    :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`,
    which the adaptive search needs for its sizing anyway, rejects inputs beyond it.

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
        Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending (so the
        emerging number of reference points ascends)
        The first `n_candidates` elements will hold the results.
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::generate_adaptive_js_comp_test_candidates`, whose argument names are
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

def construct_adaptive_neighborhoods(
        gene_means,
        residuals,
        k_start,
        k_step,
        k_max,
        tau=0.1,
        mad_distance_factor=1.0,
        max_pooled_residuals=0,
        min_study_neighbors=1,
):
    r"""Grow adaptive, dispersion-controlled neighborhoods over all studies' pooled gene means (Issue #217)

    The heteroscedastic neighborhood construction of Issue #217, after Aaron Schroeder's
    noise-model growth: instead of one fixed neighbor count, each neighborhood grows while the
    mean absolute residual of its pool stays stable, and the next reference point is placed
    beyond it. The number of reference points emerges from the growth knobs `k_start`, `k_step`
    and `k_max`, all counted in pooled entries (one entry is one gene of one study).

    **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
    n_studies` flat entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted
    ascending; NaN means sort last and take no part. A neighborhood is a contiguous run
    `a..b` of that sorted order, returned as `pooled_neighborhood_range`, so its per-study
    membership follows by decoding every position back to its study, and the studies' counts
    `n_neighbors_per_point` add up to `b - a + 1` exactly. These ranges are what
    :func:`tensor_omics.check_neighborhood_overlaps`
    takes; unlike the fixed-k ranges they are not extended over tied means, because a range is
    the exact membership here.

    **Growth of one neighborhood** around its reference point `x_star`, the mean at its seed
    position. Entries join nearest first by `|mean - x_star|`, a tie going to the lower mean.
    The dispersion `S` is the mean absolute value of every non-NaN residual of every entry in
    the pool.

    1. The first `k_start` entries join unconditionally.
    2. If they hold fewer than
    ``ADAPTIVE_MIN_VALID_RESIDUALS``
    (`10`) non-NaN residuals, the neighborhood stops there with stop reason
    ``ADAPTIVE_STOP_TOO_FEW_RESIDUALS``
    (`6`). If their residuals are all zero, it stops with
    ``ADAPTIVE_STOP_ZERO_DISPERSION``
    (`5`), since no relative change is defined.
    3. Otherwise it grows in rounds: the next `min(k_step, k_max - count)` nearest entries are
    staged and `S_new` of the pool with them is compared with the previous round's `S_old`.
    If `(S_new - S_old) / S_old > tau` the round is discarded and growth stops
    (``ADAPTIVE_STOP_TAU``, `1`);
    otherwise the round is committed and `S_old = S_new`. A round without any non-NaN
    residual leaves `S` unchanged and is committed.
    4. Growth also stops once the neighborhood holds `k_max` entries
    (``ADAPTIVE_STOP_K_MAX``, `2`),
    or once no entry is left on either side
    (``ADAPTIVE_STOP_EXHAUSTED``, `4`);
    `k_max` is checked first. A round cut short by the pool's end is still evaluated.
    5. With `max_pooled_residuals > 0`, an adaptive round stages whole entries only while the
    pool's non-NaN residuals stay at or below that cap; NaN residuals do not count toward
    it, and the `k_start` phase is not capped. A round cut short by the cap is evaluated as
    usual and, if committed, ends growth; a round that could stage nothing ends it at once.
    Both stop with
    ``ADAPTIVE_STOP_RESIDUAL_CAP`` (`3`).

    **Reference points.** The first seed is the smallest pooled mean. After neighborhood `i`
    ends at sorted position `b_i`, the next target is `mean(b_i) + mad_distance_factor *
    MAD_i`, `MAD_i` being the raw median absolute deviation of neighborhood `i`'s pooled means
    (see ``calc_sorted_slice_mad``),
    and the next seed is the position after `b_i` whose mean is nearest that target (a tie
    going to the lower mean; the largest pooled mean if the target lies beyond it); where
    several positions share that mean, the seed is the first of them. Seeds therefore
    strictly advance, which bounds the number of reference points by `N`; the
    construction ends with the neighborhood that reaches the largest pooled mean.

    **Construction status.**
    ``ADAPTIVE_STATUS_TOO_FEW_MEANS``
    (`1`) when fewer than `k_start` pooled means are non-NaN (all-NaN means included): no
    neighborhood is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction
    always completes, so that it can be inspected, and the status reports
    ``ADAPTIVE_STATUS_TOO_FEW_RESIDUALS``
    (`3`) if any neighborhood stopped with too few residuals, else
    ``ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD``
    (`2`) if any neighborhood has fewer than `min_study_neighbors` entries of some study, else
    ``ADAPTIVE_STATUS_OK`` (`0`). This
    priority is fixed, whichever neighborhood failed first.

    All outputs sized by `N` are filled only for their first `n_points` reference points.

    Parameters
    ----------
    gene_means : np.ndarray[np.float64] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Mean expression of every gene in every study, NaN for a missing gene
        NaN is permitted for this value.
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Signed residuals of every replicate of every gene in every study, NaN for a missing value
        NaN is permitted for this value.
    k_start : int
        Pooled entries every neighborhood takes unconditionally
        The minimum valid value is `1`.
    k_step : int
        Pooled entries staged per adaptive growth round
        The minimum valid value is `1`.
    k_max : int
        Largest number of pooled entries in one neighborhood
        The minimum valid value is `k_start`.
    tau : float, optional, default 0.1
        Largest relative increase of the dispersion an adaptive round may cause and still
        be committed
        The minimum valid value is `0.0`.
        The default value is `0.1`.
    mad_distance_factor : float, optional, default 1.0
        Multiple of a neighborhood's median absolute deviation that the next reference
        point's target lies beyond the neighborhood's largest mean
        The minimum valid value is `0.0`.
        The default value is `1.0`.
    max_pooled_residuals : int, optional, default 0
        Cap on the non-NaN residuals adaptive rounds may grow a neighborhood's pool to; 0
        for no cap
        The minimum valid value is `0`.
        The default value is `0`.
    min_study_neighbors : int, optional, default 1
        Fewest entries of each study every neighborhood must have for the construction
        status to stay ok
        The minimum valid value is `1`.
        The default value is `1`.

    Returns
    -------
    dict
        with keys:

        x_star : np.ndarray[np.float64] of shape (max_n_genes_all_studies*n_studies,), read-only
            Reference point of each neighborhood: the pooled mean at its seed position
            The first `n_points` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        pooled_neighborhood_range : np.ndarray[np.int32] of shape (2, max_n_genes_all_studies*n_studies,), column-major (order='F'), read-only
            For each neighborhood, its first and last position in the ascending order of the
            pooled means
            The first `n_points` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_neighbors_per_point : np.ndarray[np.int32] of shape (n_studies, max_n_genes_all_studies*n_studies,), column-major (order='F'), read-only
            For each neighborhood, how many of its pooled entries belong to each study
            The first `n_points` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        stop_reason : np.ndarray[np.int32] of shape (max_n_genes_all_studies*n_studies,), read-only
            Why each neighborhood stopped growing: `1` tau exceeded,
            `2` `k_max` reached, `3` residual cap
            reached, `4` pool exhausted, `5`
            zero dispersion, `6` too few residuals
            The first `n_points` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        neighborhood_dispersion : np.ndarray[np.float64] of shape (max_n_genes_all_studies*n_studies,), read-only
            Mean absolute non-NaN residual of each final neighborhood (the last committed
            dispersion); NaN when the neighborhood holds no non-NaN residual at all
            The first `n_points` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        neighborhood_mad : np.ndarray[np.float64] of shape (max_n_genes_all_studies*n_studies,), read-only
            Raw median absolute deviation of each final neighborhood's pooled means
            The first `n_points` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        max_n_neighbors : int
            Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        construction_status : int
            Overall outcome: `0` ok, `1` too few
            means, `2` a study neighborhood below
            `min_study_neighbors`, `3` a neighborhood with too
            few residuals

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::construct_adaptive_neighborhoods`, whose argument names are
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
    x_star = np.empty((max_n_genes_all_studies*n_studies,), dtype=np.float64, order='C')
    pooled_neighborhood_range = np.empty((2, max_n_genes_all_studies*n_studies,), dtype=np.int32, order='F')
    n_neighbors_per_point = np.empty((n_studies, max_n_genes_all_studies*n_studies,), dtype=np.int32, order='F')
    stop_reason = np.empty((max_n_genes_all_studies*n_studies,), dtype=np.int32, order='C')
    neighborhood_dispersion = np.empty((max_n_genes_all_studies*n_studies,), dtype=np.float64, order='C')
    neighborhood_mad = np.empty((max_n_genes_all_studies*n_studies,), dtype=np.float64, order='C')
    max_n_neighbors = ctypes.c_int(0)
    construction_status = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.construct_adaptive_neighborhoods_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        gene_means,
        residuals,
        ctypes.byref(ctypes.c_int(k_start)),
        ctypes.byref(ctypes.c_int(k_step)),
        ctypes.byref(ctypes.c_int(k_max)),
        ctypes.byref(n_points),
        x_star,
        pooled_neighborhood_range,
        n_neighbors_per_point,
        stop_reason,
        neighborhood_dispersion,
        neighborhood_mad,
        ctypes.byref(max_n_neighbors),
        ctypes.byref(construction_status),
        ctypes.byref(ctypes.c_double(tau)),
        ctypes.byref(ctypes.c_double(mad_distance_factor)),
        ctypes.byref(ctypes.c_int(max_pooled_residuals)),
        ctypes.byref(ctypes.c_int(min_study_neighbors)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CONSTRUCT_ADAPTIVE_NEIGHBORHOODS_ARGUMENTS, _CONSTRUCT_ADAPTIVE_NEIGHBORHOODS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    x_star.flags.writeable = False
    pooled_neighborhood_range.flags.writeable = False
    n_neighbors_per_point.flags.writeable = False
    stop_reason.flags.writeable = False
    neighborhood_dispersion.flags.writeable = False
    neighborhood_mad.flags.writeable = False

    return {
        "x_star": x_star[..., :n_points.value],
        "pooled_neighborhood_range": pooled_neighborhood_range[..., :n_points.value],
        "n_neighbors_per_point": n_neighbors_per_point[..., :n_points.value],
        "stop_reason": stop_reason[..., :n_points.value],
        "neighborhood_dispersion": neighborhood_dispersion[..., :n_points.value],
        "neighborhood_mad": neighborhood_mad[..., :n_points.value],
        "max_n_neighbors": max_n_neighbors.value,
        "construction_status": construction_status.value,
    }

def check_neighborhood_overlaps(
        neighborhood_range,
        min_neighbor_overlap,
):
    r"""Test whether every pair of consecutive neighborhoods overlaps by at least a minimum fraction

    Ported from 125-stabilize-jscomp's `test_neighborhood_overlaps_helper`: the first
    admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, using the
    `[min_idx, max_idx]` neighborhood spans
    :func:`tensor_omics.construct_neighborhoods_ranged`
    produces. Named `check_*` rather than 125's `test_*`, a deliberate deviation from the
    verbatim port: the generated R binding is published under the Fortran name, and the R test
    harness (`r/test_helpers.R`'s `run_all_tests`) discovers every `test_`-prefixed name in the
    environment as a test case to run with no arguments -- a `test_`-prefixed export would be
    swept up and fail every R suite that sources the package, not just this module's own.

    Parameters
    ----------
    neighborhood_range : np.ndarray[np.int32] of shape (2, n_points,), column-major (order='F')
        For each reference point, the `[min_idx, max_idx]` neighborhood span, as produced
        by construct_neighborhoods_ranged_impl
        The minimum valid value is `1`.
    min_neighbor_overlap : float
        Minimum fractional overlap two consecutive neighborhoods must have
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    all_have_min_neighbor_overlap : bool
        `True` if every pair of consecutive neighborhoods overlaps by at least
        `min_neighbor_overlap`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::check_neighborhood_overlaps`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_range = np.asfortranarray(neighborhood_range, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_range' must be an array of np.int32: {error}") from None
    if neighborhood_range.ndim != 2:
        raise ValueError(f"'neighborhood_range' must have 2 dimensions, but has {neighborhood_range.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = neighborhood_range.shape[1]

    # outputs and work arrays, which the caller never sees
    all_have_min_neighbor_overlap = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.check_neighborhood_overlaps_c(
        neighborhood_range,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(min_neighbor_overlap)),
        ctypes.byref(all_have_min_neighbor_overlap),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENTS, _CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENT_SOURCES)

    return all_have_min_neighbor_overlap.value

def check_mean_pmf_min_counts(
        mean_pmf_counts,
        n_bins_per_point,
        min_count,
):
    r"""Test whether every bin of a mean pmf reaches a minimum absolute count

    Ported from 125-stabilize-jscomp's `test_mean_pmf_min_counts_helper`: the second
    admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, checked once the
    first gate
    (:func:`tensor_omics.check_neighborhood_overlaps`)
    has already passed. Named `check_*` rather than 125's `test_*` for the same reason as its
    sibling above: a `test_`-prefixed R export collides with the R test harness's own
    test-discovery convention.
    Issue #187: `n_bins_per_point` scopes the reduction to each point's own valid bin range, so
    a point whose own bin count is narrower than `n_bins` (this array's second extent, i.e. the
    widest bin count any point uses) has its legitimate zero-padded columns skipped rather than
    mistaken for an occupancy failure.

    Parameters
    ----------
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_bins, n_points,), column-major (order='F')
        Absolute counts of a residual per bin for the mean pmf
        The minimum valid value is `0`.
    n_bins_per_point : np.ndarray[np.int32] of shape (n_points,)
        This point's own bin count -- only `mean_pmf_counts(1:n_bins_per_point(i_point), i_point)`
        is inspected; columns beyond it are legitimate zero-padding, not failures
        The minimum valid value is `1`.
        The maximum valid value is `n_bins`.
    min_count : int
        Minimum count each bin of the mean pmf must reach
        The minimum valid value is `0`.

    Returns
    -------
    all_bins_have_min_count : bool
        `True` if every bin within each reference point's own `n_bins_per_point`, at every
        reference point, reaches at least `min_count`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::check_mean_pmf_min_counts`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        mean_pmf_counts = np.asfortranarray(mean_pmf_counts, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_pmf_counts' must be an array of np.int32: {error}") from None
    if mean_pmf_counts.ndim != 2:
        raise ValueError(f"'mean_pmf_counts' must have 2 dimensions, but has {mean_pmf_counts.ndim}")
    try:
        n_bins_per_point = np.ascontiguousarray(n_bins_per_point, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'n_bins_per_point' must be an array of np.int32: {error}") from None
    if n_bins_per_point.ndim != 1:
        raise ValueError(f"'n_bins_per_point' must have 1 dimension, but has {n_bins_per_point.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bins = mean_pmf_counts.shape[0]
    n_points = mean_pmf_counts.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if n_bins_per_point.shape[0] != n_points:
        raise ValueError(f"'n_bins_per_point' has {n_bins_per_point.shape[0]} along axis 0, but "
            f"'mean_pmf_counts' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    all_bins_have_min_count = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.check_mean_pmf_min_counts_c(
        mean_pmf_counts,
        ctypes.byref(ctypes.c_int(n_bins)),
        n_bins_per_point,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(min_count)),
        ctypes.byref(all_bins_have_min_count),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENTS, _CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENT_SOURCES)

    return all_bins_have_min_count.value

def check_plateau_condition(
        confidence_interval,
        best_candidate_pair_confidence_interval,
        best_candidate_index,
        best_exceeded_ci_overlap_count,
        candidate_index,
        join_method,
        succeeding_ci_overlap,
):
    r"""Test one candidate pair's bootstrapped confidence intervals against the current best, and detect a JSD plateau

    Ported from 125-stabilize-jscomp's `check_plateau_condition_helper`. A search over
    candidates (finest resolution to coarsest) stops -- "plateaus" -- either when a new
    candidate is no better than the previous best (the short-circuit below: keep the previous
    best and stop searching), or once the new candidate's confidence-interval overlap with the
    previous best meets the condition `join_method` names. `join_method` replaces
    125-stabilize-jscomp's hand-rolled join-method-range validation macro entirely: the mode
    table below is itself the validation, checked against exactly the values it names.

    Parameters
    ----------
    confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F')
        JSD confidence interval `[lower, upper]` from bootstrapping, for the candidate
        pair under test
    best_candidate_pair_confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), modified in place
        JSD confidence intervals for the current best candidate pair; overwritten with
        `confidence_interval` unless the new candidate is worse
    best_candidate_index : int
        Candidate-grid index of the current best candidate pair; overwritten with
        `candidate_index` unless the new candidate is worse
    best_exceeded_ci_overlap_count : int
        Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
        best candidate pair; overwritten unless the new candidate is worse
        The minimum valid value is `0`.
    candidate_index : int
        Candidate-grid index of the candidate pair that produced `confidence_interval`
        The minimum valid value is `1`.
    join_method : str, one of 'join_min' | 'join_max' | 'join_median'
        The way to evaluate all studies' confidence-interval overlaps for the plateau
        condition: METHOD_JOIN_MIN requires every study's overlap to exceed
        `succeeding_ci_overlap`, METHOD_JOIN_MAX requires only one study's overlap to
        exceed it, and METHOD_JOIN_MEDIAN requires a majority
        (`count > (n_studies - 1) / 2`) to exceed it

    succeeding_ci_overlap : float
        Minimum fractional overlap an interval in `confidence_interval` must have with its
        respective interval in `best_candidate_pair_confidence_interval` to count as
        "exceeded"
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    dict
        with keys:

        best_candidate_index : int
            Candidate-grid index of the current best candidate pair; overwritten with
            `candidate_index` unless the new candidate is worse
        best_exceeded_ci_overlap_count : int
            Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            best candidate pair; overwritten unless the new candidate is worse
            The minimum valid value is `0`.
        plateau_found : bool
            `True` once the new candidate is no better than the previous best, or once
            `join_method`'s overlap condition is met by the new candidate

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::check_plateau_condition`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        confidence_interval = np.asfortranarray(confidence_interval, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'confidence_interval' must be an array of np.float64: {error}") from None
    if confidence_interval.ndim != 2:
        raise ValueError(f"'confidence_interval' must have 2 dimensions, but has {confidence_interval.ndim}")
    if not isinstance(best_candidate_pair_confidence_interval, np.ndarray) or best_candidate_pair_confidence_interval.dtype != np.float64:
        raise TypeError("'best_candidate_pair_confidence_interval' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if best_candidate_pair_confidence_interval.ndim != 2:
        raise ValueError(f"'best_candidate_pair_confidence_interval' must have 2 dimensions, but has {best_candidate_pair_confidence_interval.ndim}")
    if not best_candidate_pair_confidence_interval.flags.f_contiguous:
        raise ValueError("'best_candidate_pair_confidence_interval' is modified in place, so it must already be column-major (order='F')")
    best_candidate_index = ctypes.c_int(best_candidate_index)
    best_exceeded_ci_overlap_count = ctypes.c_int(best_exceeded_ci_overlap_count)
    join_method = np.array([str(join_method).lower().encode().ljust(11)], dtype="S11")

    # what the inputs already say, rather than asking for it again
    n_studies = confidence_interval.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if best_candidate_pair_confidence_interval.shape[1] != n_studies:
        raise ValueError(f"'best_candidate_pair_confidence_interval' has {best_candidate_pair_confidence_interval.shape[1]} along axis 1, but "
            f"'confidence_interval' implies n_studies == {n_studies}"
        )

    # outputs and work arrays, which the caller never sees
    plateau_found = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.check_plateau_condition_c(
        confidence_interval,
        best_candidate_pair_confidence_interval,
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(best_candidate_index),
        ctypes.byref(best_exceeded_ci_overlap_count),
        ctypes.byref(ctypes.c_int(candidate_index)),
        join_method,
        ctypes.byref(ctypes.c_double(succeeding_ci_overlap)),
        ctypes.byref(plateau_found),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CHECK_PLATEAU_CONDITION_ARGUMENTS, _CHECK_PLATEAU_CONDITION_ARGUMENT_SOURCES)

    return {
        "best_candidate_index": best_candidate_index.value,
        "best_exceeded_ci_overlap_count": best_exceeded_ci_overlap_count.value,
        "plateau_found": plateau_found.value,
    }

def check_effect_size_plateau_condition(
        global_js_divergence,
        prev_global_js_divergence,
        has_previous,
        delta_median_threshold,
        delta_max_threshold,
        delta_epsilon,
        delta_min_consecutive_transitions,
        n_consecutive_ok,
):
    r"""Test one candidate's per-study JSD against the previous admissible candidate for a relative-effect-size plateau

    Implements Issue #178's relative-effect-size plateau criterion, complementary to
    :func:`tensor_omics.check_plateau_condition`'s
    CI-overlap one: for each study `i`, the relative change in observed JSD between successive
    ADMISSIBLE parameter settings (both admissibility gates already passed),
    `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    delta_epsilon)`, summarized across studies by its median (`delta_median`, via the
    already-shipped
    :func:`tensor_omics.calc_percentile`) and maximum (`delta_max`). A
    plateau is declared once both stay under their respective thresholds for
    `delta_min_consecutive_transitions` consecutive transitions in a row -- tracked across calls
    via `n_consecutive_ok`, reset the moment either threshold is missed.

    No transition exists for the very first admissible candidate a caller ever passes in
    (`has_previous = False`): `delta`/`delta_median`/`delta_max` are all set to
    `-1.0` -- the same not-yet-computed sentinel
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    already uses for its own confidence-interval fallback, usable here for the same reason:
    every quantity this routine tracks is structurally non-negative.

    The 0.05/0.10 defaults `run_js_comp_test_parameter_search_impl` passes for
    `delta_median_threshold`/`delta_max_threshold` are Issue #178's own suggested starting
    point, explicitly not yet empirically validated -- see that routine's doc comment.

    Parameters
    ----------
    global_js_divergence : np.ndarray[np.float64] of shape (n_studies,)
        Current admissible candidate's observed global JSD per study
        The minimum valid value is `0.0`.
    prev_global_js_divergence : np.ndarray[np.float64] of shape (n_studies,)
        Previous admissible candidate's observed global JSD per study; ignored when
        `has_previous` is `False`
        The minimum valid value is `0.0`.
    has_previous : bool
        `False` for the very first admissible candidate a caller has ever passed in, where
        no transition exists to compute a relative change from
    delta_median_threshold : float
        Upper bound the median relative change across studies must stay under for a
        transition to count toward a plateau
        The minimum valid value is `above(0.0)`.
    delta_max_threshold : float
        Upper bound the largest relative change across studies must stay under for a
        transition to count toward a plateau
        The minimum valid value is `above(0.0)`.
    delta_epsilon : float
        Small constant preventing division by zero when a study's previous JSD was zero
        The minimum valid value is `above(0.0)`.
    delta_min_consecutive_transitions : int
        Number of consecutive qualifying transitions required to declare a plateau
        The minimum valid value is `1`.
    n_consecutive_ok : int
        Running count of consecutive qualifying transitions; incremented when this
        transition qualifies, reset to zero otherwise (and whenever `has_previous` is
        `False`)
        The minimum valid value is `0`.

    Returns
    -------
    dict
        with keys:

        n_consecutive_ok : int
            Running count of consecutive qualifying transitions; incremented when this
            transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            `False`)
            The minimum valid value is `0`.
        delta : np.ndarray[np.float64] of shape (n_studies,), read-only
            Per-study relative JSD change from the previous admissible candidate; `-1.0`
            throughout iff `.not. has_previous`
            A result is a value; call `.copy()` to obtain a modifiable array.
        delta_median : float
            Median of `delta` across studies; `-1.0` iff `.not. has_previous`
        delta_max : float
            Maximum of `delta` across studies; `-1.0` iff `.not. has_previous`
        plateau_found : bool
            `True` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::check_effect_size_plateau_condition`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        global_js_divergence = np.ascontiguousarray(global_js_divergence, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'global_js_divergence' must be an array of np.float64: {error}") from None
    if global_js_divergence.ndim != 1:
        raise ValueError(f"'global_js_divergence' must have 1 dimension, but has {global_js_divergence.ndim}")
    try:
        prev_global_js_divergence = np.ascontiguousarray(prev_global_js_divergence, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'prev_global_js_divergence' must be an array of np.float64: {error}") from None
    if prev_global_js_divergence.ndim != 1:
        raise ValueError(f"'prev_global_js_divergence' must have 1 dimension, but has {prev_global_js_divergence.ndim}")
    n_consecutive_ok = ctypes.c_int(n_consecutive_ok)

    # what the inputs already say, rather than asking for it again
    n_studies = global_js_divergence.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if prev_global_js_divergence.shape[0] != n_studies:
        raise ValueError(f"'prev_global_js_divergence' has {prev_global_js_divergence.shape[0]} along axis 0, but "
            f"'global_js_divergence' implies n_studies == {n_studies}"
        )

    # outputs and work arrays, which the caller never sees
    delta = np.empty((n_studies,), dtype=np.float64, order='C')
    delta_median = ctypes.c_double(0)
    delta_max = ctypes.c_double(0)
    plateau_found = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.check_effect_size_plateau_condition_c(
        global_js_divergence,
        prev_global_js_divergence,
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_bool(has_previous)),
        ctypes.byref(ctypes.c_double(delta_median_threshold)),
        ctypes.byref(ctypes.c_double(delta_max_threshold)),
        ctypes.byref(ctypes.c_double(delta_epsilon)),
        ctypes.byref(ctypes.c_int(delta_min_consecutive_transitions)),
        ctypes.byref(n_consecutive_ok),
        delta,
        ctypes.byref(delta_median),
        ctypes.byref(delta_max),
        ctypes.byref(plateau_found),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENTS, _CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    delta.flags.writeable = False

    return {
        "n_consecutive_ok": n_consecutive_ok.value,
        "delta": delta,
        "delta_median": delta_median.value,
        "delta_max": delta_max.value,
        "plateau_found": plateau_found.value,
    }

def create_mean_pmf(
        pmfs,
        counts,
        included_n_reps,
):
    r"""Build the consensus pmf and its histogram counts from all studies' pmfs

    Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_helper`.

    Known limitation: averages over all n_studies including the study being compared against it,
    rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    the fix.

    Parameters
    ----------
    pmfs : np.ndarray[np.float64] of shape (n_bins, n_points, n_studies,), column-major (order='F')
        Per-study probabilities of each bin per reference point, from
        :func:`tensor_omics.build_residual_histograms`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    counts : np.ndarray[np.int32] of shape (n_bins, n_points, n_studies,), column-major (order='F')
        Absolute counts of a residual per bin for `pmfs`
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates (included ones) per reference point, per study
        The minimum valid value is `0`.

    Returns
    -------
    dict
        with keys:

        mean_pmf : np.ndarray[np.float64] of shape (n_bins, n_points,), column-major (order='F'), read-only
            The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            known-limitation note above
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN replicates (included ones) per reference point for `mean_pmf`,
            summed across all n_studies
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_counts : np.ndarray[np.int32] of shape (n_bins, n_points,), column-major (order='F'), read-only
            Absolute counts of a residual per bin for the mean pmf -> `sum(counts, dim=3)`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::create_mean_pmf`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pmfs = np.asfortranarray(pmfs, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pmfs' must be an array of np.float64: {error}") from None
    if pmfs.ndim != 3:
        raise ValueError(f"'pmfs' must have 3 dimensions, but has {pmfs.ndim}")
    try:
        counts = np.asfortranarray(counts, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'counts' must be an array of np.int32: {error}") from None
    if counts.ndim != 3:
        raise ValueError(f"'counts' must have 3 dimensions, but has {counts.ndim}")
    try:
        included_n_reps = np.asfortranarray(included_n_reps, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'included_n_reps' must be an array of np.int32: {error}") from None
    if included_n_reps.ndim != 2:
        raise ValueError(f"'included_n_reps' must have 2 dimensions, but has {included_n_reps.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bins = pmfs.shape[0]
    n_points = pmfs.shape[1]
    n_studies = pmfs.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if counts.shape[0] != n_bins:
        raise ValueError(f"'counts' has {counts.shape[0]} along axis 0, but "
            f"'pmfs' implies n_bins == {n_bins}"
        )
    if counts.shape[1] != n_points:
        raise ValueError(f"'counts' has {counts.shape[1]} along axis 1, but "
            f"'pmfs' implies n_points == {n_points}"
        )
    if included_n_reps.shape[0] != n_points:
        raise ValueError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but "
            f"'pmfs' implies n_points == {n_points}"
        )
    if counts.shape[2] != n_studies:
        raise ValueError(f"'counts' has {counts.shape[2]} along axis 2, but "
            f"'pmfs' implies n_studies == {n_studies}"
        )
    if included_n_reps.shape[1] != n_studies:
        raise ValueError(f"'included_n_reps' has {included_n_reps.shape[1]} along axis 1, but "
            f"'pmfs' implies n_studies == {n_studies}"
        )

    # outputs and work arrays, which the caller never sees
    mean_pmf = np.empty((n_bins, n_points,), dtype=np.float64, order='F')
    mean_pmf_included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    mean_pmf_counts = np.empty((n_bins, n_points,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.create_mean_pmf_c(
        pmfs,
        counts,
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_studies)),
        included_n_reps,
        mean_pmf,
        mean_pmf_included_n_reps,
        mean_pmf_counts,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CREATE_MEAN_PMF_ARGUMENTS, _CREATE_MEAN_PMF_ARGUMENT_SOURCES)

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

    Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_only_helper`: useful where the
    mean pmf's own counts don't matter, e.g. for the bootstrap confidence interval a later
    stage of this port adds.

    Known limitation: averages over all n_studies including the study being compared against it,
    rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    the fix.

    Parameters
    ----------
    pmfs : np.ndarray[np.float64] of shape (n_bins, n_points, n_studies,), column-major (order='F')
        Per-study probabilities of each bin per reference point, from
        :func:`tensor_omics.build_residual_histograms`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    mean_pmf : np.ndarray[np.float64] of shape (n_bins, n_points,), column-major (order='F'), read-only
        The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
        known-limitation note above
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::create_mean_pmf_only`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pmfs = np.asfortranarray(pmfs, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pmfs' must be an array of np.float64: {error}") from None
    if pmfs.ndim != 3:
        raise ValueError(f"'pmfs' must have 3 dimensions, but has {pmfs.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bins = pmfs.shape[0]
    n_points = pmfs.shape[1]
    n_studies = pmfs.shape[2]

    # outputs and work arrays, which the caller never sees
    mean_pmf = np.empty((n_bins, n_points,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.create_mean_pmf_only_c(
        pmfs,
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_studies)),
        mean_pmf,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CREATE_MEAN_PMF_ONLY_ARGUMENTS, _CREATE_MEAN_PMF_ONLY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    mean_pmf.flags.writeable = False

    return mean_pmf

def bootstrap_histogram(
        n_bootstraps,
        mean_pmf_counts,
        mean_pmf_included_n_reps,
        included_n_reps,
        confidence_interval,
        two_sided_bootstrapping_significance_level=2.5,
        random_seed=42,
):
    r"""Bootstrap a confidence interval for each study's global JSD by resampling the pooled consensus histogram

    Ported from 125-stabilize-jscomp's `bootstrap_histogram_helper`. Resamples from the
    POOLED/consensus histogram counts (`mean_pmf_counts`), not from each study's own histogram
    -- this is 125's deliberate design, preserved as-is (see the plan for this port). Draws
    random numbers via
    ``random_multinomial``, so this implementation is
    deliberately impure, matching the project's existing precedent for the other
    permutation-test module's purity
    (:func:`tensor_omics.perform_permutation_test`).

    `confidence_interval` is both an input and an output: its incoming `[lower, upper]` values
    seed every slot of the top-k/bottom-k heaps (`tmp_bootstrapping_top_k_jsds`) -- ported
    verbatim from 125, which fills the whole heap with the *same* incoming reference value
    rather than the usual plus/minus-infinity heap initialization, so the observed
    (pre-bootstrap) value can only be displaced by a strictly more extreme bootstrap draw. On
    return it holds `[largest of the n_bootstrapping_top_k_jsds smallest bootstrap draws,
    smallest of the n_bootstrapping_top_k_jsds largest bootstrap draws]`.

    `mean_pmf_counts`/`tmp_pmfs`/`tmp_mean_pmf` are laid out bin-major (`(n_bins, n_points[,
    n_studies])`), matching
    :func:`tensor_omics.create_mean_pmf` and its
    own `create_mean_pmf_only_impl` above -- both ported from 125, whose own convention this
    is. The already-shipped
    :func:`tensor_omics.compute_divergence_per_reference_point`
    predates 125's own port and instead takes its pmf arguments point-major
    (`(n_points, n_bins)`); the two calls below bridge the two conventions with an explicit
    `transpose`, rather than picking one shape and silently reinterpreting the other's memory
    under it (which would scramble every non-square `(n_bins, n_points)` histogram).

    A GSL allocation failure in `create_rng` is a genuine runtime error no input check could
    have foreseen (codegen_guide.md Sec 5.14): every work array and `confidence_interval` are
    then left untouched (arrays not yet written to keep their caller-visible defined state) and
    `ierr` reports `ERR_ALLOC_FAIL`. A `random_multinomial` draw failing (which validated,
    internally-consistent inputs should never trigger) is likewise folded into `ierr`, first
    failure only, without stopping the resampling already in flight.

    Parameters
    ----------
    n_bootstraps : int
        Number of bootstrap resamples to perform
        The minimum valid value is `1`.
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_bins, n_points,), column-major (order='F')
        Absolute counts of a residual per bin for the pooled/consensus pmf, from
        create_mean_pmf_impl -- resampled with replacement each bootstrap
        The minimum valid value is `0`.
    mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates (included ones) per reference point for the pooled pmf
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates (included ones) per reference point, per study --
        how many elements are drawn (with replacement) from the pooled pool per study
        The minimum valid value is `0`.
    confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), modified in place
        Confidence interval to be bootstrapped -- incoming values are the reference values
        that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
        `[lower, upper]` interval per study
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
        otherwise used here
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::bootstrap_histogram`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        mean_pmf_counts = np.asfortranarray(mean_pmf_counts, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_pmf_counts' must be an array of np.int32: {error}") from None
    if mean_pmf_counts.ndim != 2:
        raise ValueError(f"'mean_pmf_counts' must have 2 dimensions, but has {mean_pmf_counts.ndim}")
    try:
        mean_pmf_included_n_reps = np.ascontiguousarray(mean_pmf_included_n_reps, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_pmf_included_n_reps' must be an array of np.int32: {error}") from None
    if mean_pmf_included_n_reps.ndim != 1:
        raise ValueError(f"'mean_pmf_included_n_reps' must have 1 dimension, but has {mean_pmf_included_n_reps.ndim}")
    try:
        included_n_reps = np.asfortranarray(included_n_reps, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'included_n_reps' must be an array of np.int32: {error}") from None
    if included_n_reps.ndim != 2:
        raise ValueError(f"'included_n_reps' must have 2 dimensions, but has {included_n_reps.ndim}")
    if not isinstance(confidence_interval, np.ndarray) or confidence_interval.dtype != np.float64:
        raise TypeError("'confidence_interval' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if confidence_interval.ndim != 2:
        raise ValueError(f"'confidence_interval' must have 2 dimensions, but has {confidence_interval.ndim}")
    if not confidence_interval.flags.f_contiguous:
        raise ValueError("'confidence_interval' is modified in place, so it must already be column-major (order='F')")

    # what the inputs already say, rather than asking for it again
    n_bins = mean_pmf_counts.shape[0]
    n_points = mean_pmf_counts.shape[1]
    n_studies = included_n_reps.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if mean_pmf_included_n_reps.shape[0] != n_points:
        raise ValueError(f"'mean_pmf_included_n_reps' has {mean_pmf_included_n_reps.shape[0]} along axis 0, but "
            f"'mean_pmf_counts' implies n_points == {n_points}"
        )
    if included_n_reps.shape[0] != n_points:
        raise ValueError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but "
            f"'mean_pmf_counts' implies n_points == {n_points}"
        )
    if confidence_interval.shape[1] != n_studies:
        raise ValueError(f"'confidence_interval' has {confidence_interval.shape[1]} along axis 1, but "
            f"'included_n_reps' implies n_studies == {n_studies}"
        )

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.bootstrap_histogram_c(
        ctypes.byref(ctypes.c_int(n_bootstraps)),
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_studies)),
        mean_pmf_counts,
        mean_pmf_included_n_reps,
        included_n_reps,
        confidence_interval,
        ctypes.byref(ctypes.c_double(two_sided_bootstrapping_significance_level)),
        ctypes.byref(ctypes.c_int(random_seed)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BOOTSTRAP_HISTOGRAM_ARGUMENTS, _BOOTSTRAP_HISTOGRAM_ARGUMENT_SOURCES)

    return None

def run_js_comp_test(
        n_neighbors,
        gene_means,
        gene_means_perms,
        residuals,
        x_star,
        n_permutations=1000,
        random_seed=42,
        min_residuals_per_bin=10,
        m_min=3,
        m_max=120,
        gamma_occupancy=1.25,
        lower_residual_range_quantile=0.05,
        upper_residual_range_quantile=0.95,
):
    r"""Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting

    Ported from 125-stabilize-jscomp's `js_comp_test_helper`, restructured for Issue #187's
    occupancy-constrained per-neighborhood histogram binning into three passes, mirroring
    :func:`tensor_omics.run_js_comp_test_parameter_search`'s
    own Pass A/B/C split (Issue #187's own Steps 2.5/2.6):

    - Pass A (per study): builds every study's neighborhoods
    (:func:`tensor_omics.construct_neighborhoods_ranged`),
    writing into `neighborhood_indices`/`neighborhood_range`, which already retain every
    study's own values simultaneously (both are real `intent(out)` arguments sized
    `(..., n_points, n_studies)` -- unlike `run_js_comp_test_parameter_search_impl`, no new
    buffer was needed for this). Unlike that routine, there is no admissibility gate here, so
    Pass A always runs to completion for every study.
    - Pass B (per point, sequential -- see the implementation body's own comment for why): pools
    every study's residuals for one reference point at a time
    (:func:`tensor_omics.gather_pooled_neighborhood_residuals`,
    now published as its own entry point) and runs Issue #187's occupancy-constrained
    bin-count search on the pooled result
    (:func:`tensor_omics.determine_bin_count_occupancy`),
    deciding `n_bins_per_point(i_point)` independently for every reference point, plus the
    `occupancy_failed`/`n_pooled_residuals`/`min_bin_occupancy`/`mean_bin_occupancy`/
    `max_bin_occupancy`/`sturges_bins`/`fd_bins` diagnostics. `max_n_bins_per_point`
    (`maxval(n_bins_per_point(1:n_points))`) is derived once after Pass B and replaces the old
    caller-supplied scalar `n_bins` everywhere downstream.
    - Pass C (per study): re-gathers this study's residual values from the neighbor indices Pass
    A already computed, then builds its residual histograms at the real per-point bin counts
    (:func:`tensor_omics.build_residual_histograms`).

    A neighbor gene index outside `[1, max_n_genes_all_studies]` -- which Pass A produces when
    `n_neighbors` exceeds a study's gene count -- is reported by Pass B as `ERR_INVALID_INPUT`,
    and the routine returns right there, before Pass C would read `residuals` out of bounds.

    After Pass C, the pipeline continues exactly as before: pools the per-study pmfs into the
    consensus pmf
    (:func:`tensor_omics.create_mean_pmf`), computes
    each study's observed JSD against that consensus
    (:func:`tensor_omics.compute_divergence_per_reference_point`/:func:`tensor_omics.compute_weighted_global_divergence`,
    called with the consensus pmf as the second argument), runs the permutation test
    (:func:`tensor_omics.gjct_permutation_test`), and
    finally re-derives each study's pmf/JSD/weights/global JSD from its own UNTOUCHED `counts`
    via
    :func:`tensor_omics.calc_pmf` -- `mean_pmf`/`mean_pmf_counts`
    are NOT re-derived, since they are invariant across permutations by construction (the
    permutation test above only resamples its own scratch copies, never `mean_pmf_counts`
    itself), exactly as 125 relies on.

    **Behavioral asymmetry vs.
    :func:`tensor_omics.run_js_comp_test_parameter_search`
    -- read before using this entry point where inadequately-supported neighborhoods must be
    rejected:** unlike that routine, THIS one has NO multi-candidate fallback and NO
    admissibility gate at all (no
    :func:`tensor_omics.check_neighborhood_overlaps`,
    no :func:`tensor_omics.check_mean_pmf_min_counts`,
    no early exit). A reference point whose Pass B occupancy search fails even at `m_min`
    (`occupancy_failed(i_point) = True_c_bool`) still gets a real histogram built, at
    `n_bins_per_point(i_point) == m_min`, and that point still contributes to
    `global_js_divergence` exactly like every other point -- its contribution is down-weighted
    only by `included_n_reps` (an orthogonal quantity: how many non-NaN replicates it has), never
    by bin sparsity. A caller that needs inadequately-supported neighborhoods rejected outright
    should use `run_js_comp_test_parameter_search_impl` instead, which gates on exactly this via
    `check_mean_pmf_min_counts_impl`.

    **A real, deliberate change to this routine's public array shapes (Issue #187):** the old
    mandatory scalar input `n_bins` is gone -- there is no way for a caller to know the right bin
    count in advance, since it is now genuinely computed inside this routine by Pass B's
    occupancy search, independently per reference point. Every array whose bin-sized dimension
    used to be sized by that input (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`,
    `tmp_counts_point_major`, `tmp_pmf_point_major`, `tmp_permutation_mean_pmf_counts`,
    `tmp_permutation_counts`, `tmp_permutation_pmfs`) is now sized to the fixed compile-time
    ceiling ``MAX_N_BINS`` (`256`)
    instead, exactly mirroring how `run_js_comp_test_parameter_search_impl`'s own
    `tmp_counts_point_major`/`tmp_pmf_point_major`/`tmp_pmfs`/`tmp_counts` etc. have been sized
    since Issue #187's earlier steps. The new `max_n_bins_per_point` output tells a caller how many of the
    LEADING bins/rows of each of those arrays are actually meaningful
    (`maxval(n_bins_per_point(1:n_points))`); the rest is unused padding. The generator's own
    result-size trimming directive cannot express this trim, because it only ever trims an
    array's LAST declared extent, and bins is the FIRST declared extent of every one of those
    arrays -- so a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves, exactly as a
    caller of `run_js_comp_test_parameter_search_impl`'s own jagged `trace_*` arrays already has
    to.

    `x_star` is an ordinary input here, not computed by this routine -- 125's own
    `js_comp_test_helper` takes it the same way, since a caller running several studies/several
    parameter settings is expected to compute the reference points once
    (:func:`tensor_omics.pool_means`) and reuse
    them consistently.

    `construct_neighborhoods_ranged_impl` reports neighbor gene INDICES, not gathered residual
    values (unlike its distance-sort sibling
    :func:`tensor_omics.construct_neighborhoods`),
    so Pass C gathers each neighbor's actual residual values from `residuals` itself
    (one slice of `tmp_neighborhood_residuals_gathered` per reference point, reused for every
    study) before building that point's histogram with `build_residual_histograms_impl` for
    that single point, whose bins are copied straight into the matching column of
    `counts`/`pmfs` -- no transpose is needed there. `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts`
    here are BIN-major (`(256, n_points, n_studies)`) to match
    :func:`tensor_omics.create_mean_pmf`'s own
    convention, while `compute_divergence_per_reference_point_impl` and `calc_pmf_impl` are
    POINT-major (`(n_points, max_n_bins_per_point)`): the calls to those two after Pass C bridge
    with an explicit `transpose`, exactly as
    :func:`tensor_omics.bootstrap_histogram` and
    :func:`tensor_omics.gjct_permutation_test` already do.

    Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is returned
    in `ierr`, and the routine returns right there: `pmfs`, `js_divergences`, `weights` and
    `global_js_divergence` then hold the pre-permutation values, not the final re-derived ones.

    Parameters
    ----------
    n_neighbors : int
        Number of neighbors per neighborhood
        The minimum valid value is `1`.
    gene_means : np.ndarray[np.float64] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Per-gene mean expression values for all studies
        NaN is permitted for this value.
    gene_means_perms : np.ndarray[np.int32] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Per-study sorting permutation for `gene_means` (ascending, NaN last)
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies`.
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Matrix of signed residuals per study
        NaN is permitted for this value.
    x_star : np.ndarray[np.float64] of shape (n_points,)
        Mean-expression reference points
        NaN is permitted for this value.
    n_permutations : int, optional, default 1000
        Number of permutations, forwarded to gjct_permutation_test_impl
        The minimum valid value is `0`.
        The default value is `1000`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a candidate bin count to
        be admissible in Pass B's occupancy search, forwarded to
        determine_bin_count_occupancy_impl
        The minimum valid value is `0`.
        The default value is `10`.
    m_min : int, optional, default 3
        Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
        forwarded to determine_bin_count_occupancy_impl
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count Pass B's occupancy search will ever test (M_max),
        forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
        determine_bin_count_occupancy_impl clamps it up to `m_min` internally
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    gamma_occupancy : float, optional, default 1.25
        Geometric growth factor for Pass B's occupancy search's coarse search stage,
        forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
        advances
        The minimum valid value is `above(1.0)`.
        The default value is `1.25`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for each reference point's own lower residual-range bound,
        forwarded to determine_bin_count_occupancy_impl
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for each reference point's own upper residual-range bound,
        forwarded to determine_bin_count_occupancy_impl
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    dict
        with keys:

        neighborhood_indices : np.ndarray[np.int32] of shape (n_neighbors, n_points, n_studies,), column-major (order='F'), read-only
            Gene indices of the selected neighborhood, per reference point, per study (Pass A)
            A result is a value; call `.copy()` to obtain a modifiable array.
        neighborhood_range : np.ndarray[np.int32] of shape (2, n_points, n_studies,), column-major (order='F'), read-only
            For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
            produced by construct_neighborhoods_ranged_impl (Pass A)
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_bins_per_point : np.ndarray[np.int32] of shape (n_points,), read-only
            This reference point's own selected histogram bin count (Issue #187's `M_j`), from
            Pass B's occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood
            may use a different bin count
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (n_points,), read-only
            This reference point's own lower residual-range bound (R_low), from Pass B's
            occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            different, asymmetric range (Step 3)
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (n_points,), read-only
            This reference point's own upper residual-range bound (R_high), from Pass B's
            occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            different, asymmetric range (Step 3)
            A result is a value; call `.copy()` to obtain a modifiable array.
        max_n_bins_per_point : int
            The widest `n_bins_per_point` value across all `n_points` reference points
            (`maxval(n_bins_per_point(1:n_points))`), derived once after Pass B. The number of
            leading, meaningful bins/rows in `pmfs`, `counts`, `mean_pmf` and `mean_pmf_counts`
            below (and of the bin-sized work arrays) -- those are all declared
            with a fixed 256-bin ceiling (MAX_N_BINS) rather than a caller-supplied bin count,
            since `n_bins_per_point` can no longer be known by a caller in advance. A Python/R
            caller must slice `[:max_n_bins_per_point, ...]` themselves: the generator's own result-size
            trimming directive cannot express this trim, because it only ever trims an array's
            LAST declared extent, and bins is the FIRST declared extent of every one of those
            arrays
        occupancy_failed : np.ndarray[np.bool_] of shape (n_points,), read-only
            This reference point's `occupancy_failed` flag from Pass B
            (determine_bin_count_occupancy_impl) -- `True` iff even `m_min` bins could not
            satisfy the occupancy criterion for it. See this routine's own doc block above for
            the behavioral asymmetry this implies vs. run_js_comp_test_parameter_search_impl: a
            `True` point here still gets a real histogram and still contributes to
            `global_js_divergence`, it is never rejected
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_pooled_residuals : np.ndarray[np.int32] of shape (n_points,), read-only
            This reference point's pooled residual count (N_j) from Pass B
            (determine_bin_count_occupancy_impl)
            A result is a value; call `.copy()` to obtain a modifiable array.
        min_bin_occupancy : np.ndarray[np.int32] of shape (n_points,), read-only
            This reference point's minimum bin occupancy at `n_bins_per_point`, from Pass B
            (determine_bin_count_occupancy_impl)
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_bin_occupancy : np.ndarray[np.float64] of shape (n_points,), read-only
            This reference point's mean bin occupancy at `n_bins_per_point`, from Pass B
            (determine_bin_count_occupancy_impl)
            A result is a value; call `.copy()` to obtain a modifiable array.
        max_bin_occupancy : np.ndarray[np.int32] of shape (n_points,), read-only
            This reference point's maximum bin occupancy at `n_bins_per_point`, from Pass B
            (determine_bin_count_occupancy_impl)
            A result is a value; call `.copy()` to obtain a modifiable array.
        sturges_bins : np.ndarray[np.int32] of shape (n_points,), read-only
            This reference point's Sturges' rule bin-count diagnostic from Pass B
            (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            decision
            A result is a value; call `.copy()` to obtain a modifiable array.
        fd_bins : np.ndarray[np.int32] of shape (n_points,), read-only
            This reference point's Freedman-Diaconis rule bin-count diagnostic from Pass B
            (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            decision
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmfs : np.ndarray[np.float64] of shape (256, n_points, n_studies,), column-major (order='F'), read-only
            `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`. `256`
            = MAX_N_BINS, a fixed ceiling (see `max_n_bins_per_point` above) -- only rows `1:max_n_bins_per_point`
            are meaningful; a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves
            A result is a value; call `.copy()` to obtain a modifiable array.
        counts : np.ndarray[np.int32] of shape (256, n_points, n_studies,), column-major (order='F'), read-only
            Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
            `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Count of non-NaN replicates (included ones) per reference point, per study
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf : np.ndarray[np.float64] of shape (256, n_points,), column-major (order='F'), read-only
            The consensus pmf, from create_mean_pmf_impl. `256` = MAX_N_BINS; only rows
            `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_counts : np.ndarray[np.int32] of shape (256, n_points,), column-major (order='F'), read-only
            Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS;
            only rows `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            A result is a value; call `.copy()` to obtain a modifiable array.
        js_divergences : np.ndarray[np.float64] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Per-reference-point JSD of each study against the consensus pmf
            A result is a value; call `.copy()` to obtain a modifiable array.
        weights : np.ndarray[np.float64] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Per-reference-point weights for `global_js_divergence`
            A result is a value; call `.copy()` to obtain a modifiable array.
        global_js_divergence : np.ndarray[np.float64] of shape (n_studies,), read-only
            Weighted global JSD of each study against the consensus pmf
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study from gjct_permutation_test_impl
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::run_js_comp_test`, whose argument names are
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
        gene_means_perms = np.asfortranarray(gene_means_perms, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_means_perms' must be an array of np.int32: {error}") from None
    if gene_means_perms.ndim != 2:
        raise ValueError(f"'gene_means_perms' must have 2 dimensions, but has {gene_means_perms.ndim}")
    try:
        residuals = np.asfortranarray(residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'residuals' must be an array of np.float64: {error}") from None
    if residuals.ndim != 3:
        raise ValueError(f"'residuals' must have 3 dimensions, but has {residuals.ndim}")
    try:
        x_star = np.ascontiguousarray(x_star, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x_star' must be an array of np.float64: {error}") from None
    if x_star.ndim != 1:
        raise ValueError(f"'x_star' must have 1 dimension, but has {x_star.ndim}")

    # what the inputs already say, rather than asking for it again
    n_studies = gene_means.shape[1]
    max_n_genes_all_studies = gene_means.shape[0]
    max_n_reps_all_studies = residuals.shape[0]
    n_points = x_star.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_means_perms.shape[1] != n_studies:
        raise ValueError(f"'gene_means_perms' has {gene_means_perms.shape[1]} along axis 1, but "
            f"'gene_means' implies n_studies == {n_studies}"
        )
    if residuals.shape[2] != n_studies:
        raise ValueError(f"'residuals' has {residuals.shape[2]} along axis 2, but "
            f"'gene_means' implies n_studies == {n_studies}"
        )
    if gene_means_perms.shape[0] != max_n_genes_all_studies:
        raise ValueError(f"'gene_means_perms' has {gene_means_perms.shape[0]} along axis 0, but "
            f"'gene_means' implies max_n_genes_all_studies == {max_n_genes_all_studies}"
        )
    if residuals.shape[1] != max_n_genes_all_studies:
        raise ValueError(f"'residuals' has {residuals.shape[1]} along axis 1, but "
            f"'gene_means' implies max_n_genes_all_studies == {max_n_genes_all_studies}"
        )

    # outputs and work arrays, which the caller never sees
    neighborhood_indices = np.empty((n_neighbors, n_points, n_studies,), dtype=np.int32, order='F')
    neighborhood_range = np.empty((2, n_points, n_studies,), dtype=np.int32, order='F')
    n_bins_per_point = np.empty((n_points,), dtype=np.int32, order='C')
    shared_residual_range_low = np.empty((n_points,), dtype=np.float64, order='C')
    shared_residual_range_high = np.empty((n_points,), dtype=np.float64, order='C')
    max_n_bins_per_point = ctypes.c_int(0)
    occupancy_failed = np.empty((n_points,), dtype=np.bool_, order='C')
    n_pooled_residuals = np.empty((n_points,), dtype=np.int32, order='C')
    min_bin_occupancy = np.empty((n_points,), dtype=np.int32, order='C')
    mean_bin_occupancy = np.empty((n_points,), dtype=np.float64, order='C')
    max_bin_occupancy = np.empty((n_points,), dtype=np.int32, order='C')
    sturges_bins = np.empty((n_points,), dtype=np.int32, order='C')
    fd_bins = np.empty((n_points,), dtype=np.int32, order='C')
    pmfs = np.empty((256, n_points, n_studies,), dtype=np.float64, order='F')
    counts = np.empty((256, n_points, n_studies,), dtype=np.int32, order='F')
    included_n_reps = np.empty((n_points, n_studies,), dtype=np.int32, order='F')
    mean_pmf = np.empty((256, n_points,), dtype=np.float64, order='F')
    mean_pmf_counts = np.empty((256, n_points,), dtype=np.int32, order='F')
    mean_pmf_included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    js_divergences = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    weights = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    global_js_divergence = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values = np.empty((n_studies,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        gene_means,
        gene_means_perms,
        residuals,
        x_star,
        neighborhood_indices,
        neighborhood_range,
        n_bins_per_point,
        shared_residual_range_low,
        shared_residual_range_high,
        ctypes.byref(max_n_bins_per_point),
        occupancy_failed,
        n_pooled_residuals,
        min_bin_occupancy,
        mean_bin_occupancy,
        max_bin_occupancy,
        sturges_bins,
        fd_bins,
        pmfs,
        counts,
        included_n_reps,
        mean_pmf,
        mean_pmf_counts,
        mean_pmf_included_n_reps,
        js_divergences,
        weights,
        global_js_divergence,
        p_values,
        ctypes.byref(ctypes.c_int(n_permutations)),
        ctypes.byref(ctypes.c_int(random_seed)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_double(gamma_occupancy)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_ARGUMENTS, _RUN_JS_COMP_TEST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neighborhood_indices.flags.writeable = False
    neighborhood_range.flags.writeable = False
    n_bins_per_point.flags.writeable = False
    shared_residual_range_low.flags.writeable = False
    shared_residual_range_high.flags.writeable = False
    occupancy_failed.flags.writeable = False
    n_pooled_residuals.flags.writeable = False
    min_bin_occupancy.flags.writeable = False
    mean_bin_occupancy.flags.writeable = False
    max_bin_occupancy.flags.writeable = False
    sturges_bins.flags.writeable = False
    fd_bins.flags.writeable = False
    pmfs.flags.writeable = False
    counts.flags.writeable = False
    included_n_reps.flags.writeable = False
    mean_pmf.flags.writeable = False
    mean_pmf_counts.flags.writeable = False
    mean_pmf_included_n_reps.flags.writeable = False
    js_divergences.flags.writeable = False
    weights.flags.writeable = False
    global_js_divergence.flags.writeable = False
    p_values.flags.writeable = False

    return {
        "neighborhood_indices": neighborhood_indices,
        "neighborhood_range": neighborhood_range,
        "n_bins_per_point": n_bins_per_point,
        "shared_residual_range_low": shared_residual_range_low,
        "shared_residual_range_high": shared_residual_range_high,
        "max_n_bins_per_point": max_n_bins_per_point.value,
        "occupancy_failed": occupancy_failed,
        "n_pooled_residuals": n_pooled_residuals,
        "min_bin_occupancy": min_bin_occupancy,
        "mean_bin_occupancy": mean_bin_occupancy,
        "max_bin_occupancy": max_bin_occupancy,
        "sturges_bins": sturges_bins,
        "fd_bins": fd_bins,
        "pmfs": pmfs,
        "counts": counts,
        "included_n_reps": included_n_reps,
        "mean_pmf": mean_pmf,
        "mean_pmf_counts": mean_pmf_counts,
        "mean_pmf_included_n_reps": mean_pmf_included_n_reps,
        "js_divergences": js_divergences,
        "weights": weights,
        "global_js_divergence": global_js_divergence,
        "p_values": p_values,
    }

def run_js_comp_test_adaptive(
        gene_means,
        residuals,
        pooled_neighborhood_range,
        n_permutations=1000,
        random_seed=42,
        min_residuals_per_bin=10,
        m_min=3,
        m_max=120,
        gamma_occupancy=1.25,
        lower_residual_range_quantile=0.05,
        upper_residual_range_quantile=0.95,
):
    r"""Run the JSD-Comp-Test pipeline on adaptive (Issue #217) neighborhoods

    The adaptive counterpart of
    :func:`tensor_omics.run_js_comp_test`: the
    same pipeline and the same outputs, on neighborhoods whose size varies per reference point
    and per study. Call
    :func:`tensor_omics.construct_adaptive_neighborhoods`
    first, then pass its `pooled_neighborhood_range` (trimmed to its `n_points` reference
    points) here together with the same `gene_means` and `residuals`.

    **Membership.** The pooled gene means are sorted again exactly as the construction sorts them
    (the same routine on the same input), so a pooled position means the same entry here as it
    did there. Every range `a..b` is decoded back to per-study gene lists by
    ``materialize_pooled_neighborhood``;
    the per-study counts are returned as `n_neighbors_per_point`, in the construction's own
    `(n_studies, n_points)` orientation, and add up to `b - a + 1`.

    **Pipeline, per reference point.** Pass B pools the point's residuals across all its genes
    of all studies and picks its bin count and residual range with Issue #187's occupancy
    search
    (:func:`tensor_omics.determine_bin_count_occupancy`),
    its Sturges/Freedman-Diaconis diagnostics computed with the rounded mean per-study neighbor
    count. Pass C builds each study's histogram from that study's own genes of the point. Then,
    exactly as in `run_js_comp_test_impl`: the consensus pmf
    (:func:`tensor_omics.create_mean_pmf`), each
    study's JSD against it and the weighted global JSD, the permutation test
    (:func:`tensor_omics.gjct_permutation_test`) and the
    final re-derivation from the untouched `counts`. Neighborhoods that decode to the same gene
    sets as a fixed-k run's therefore give bit-identical results. Every point is weighted by its
    non-NaN residual count, which under adaptive growth genuinely differs between points.
    Like `run_js_comp_test_impl`, there is no admissibility gate: overlap and occupancy are the
    parameter search's concern.

    **Runtime errors.** Before any histogram work, every range is checked in point order; the
    first failing one sets `ERR_INVALID_INPUT` and the routine returns at once, every output
    undefined. For each point, in this order: its first position lies after its last; its last
    position lies beyond the non-NaN pooled means (the entries a construction can ever use); a
    study has no gene in it, which would leave that study's pmf empty and the consensus pmf, an
    average over all studies, meaningless. A range a construction produced passes the first two
    checks by construction, and the third whenever its status was not an empty-study one.

    **Memory.** The work arrays are sized by the per-study gene bound `max_n_genes_all_studies`,
    the largest number of genes one study can have in one neighborhood, not by a caller-supplied
    neighbor count. The pooling buffers `tmp_pooled_residuals`/`tmp_pooled_residuals_perm`
    therefore hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies` values each,
    the size of `residuals` itself; Pass B and Pass C are sequential over points so one such
    buffer serves every point.

    The bin-sized arrays have the same fixed 256-bin
    (``MAX_N_BINS``) leading extent as in
    `run_js_comp_test_impl`, of which only the first `max_n_bins_per_point` rows are
    meaningful; a Python/R caller slices `[:max_n_bins_per_point, ...]` themselves.

    Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is returned
    in `ierr`, and the routine returns right there: `pmfs`, `js_divergences`, `weights` and
    `global_js_divergence` then hold the pre-permutation values, not the final re-derived ones.

    Parameters
    ----------
    gene_means : np.ndarray[np.float64] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Mean expression of every gene in every study, NaN for a missing gene; the same
        array the neighborhoods were constructed from
        NaN is permitted for this value.
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Signed residuals of every replicate of every gene in every study, NaN for a missing value
        NaN is permitted for this value.
    pooled_neighborhood_range : np.ndarray[np.int32] of shape (2, n_points,), column-major (order='F')
        For each neighborhood, its first and last position in the ascending order of the
        pooled means, as returned by construct_adaptive_neighborhoods
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies*n_studies`.
    n_permutations : int, optional, default 1000
        Number of permutations, forwarded to gjct_permutation_test_impl
        The minimum valid value is `0`.
        The default value is `1000`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a candidate bin count to
        be admissible in the occupancy search
        The minimum valid value is `0`.
        The default value is `10`.
    m_min : int, optional, default 3
        Smallest candidate bin count the occupancy search tests (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count the occupancy search tests (M_max); raised to `m_min`
        when smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    gamma_occupancy : float, optional, default 1.25
        Geometric growth factor of the occupancy search's coarse stage
        The minimum valid value is `above(1.0)`.
        The default value is `1.25`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for each reference point's lower residual-range bound
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for each reference point's upper residual-range bound
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.

    Returns
    -------
    dict
        with keys:

        n_neighbors_per_point : np.ndarray[np.int32] of shape (n_studies, n_points,), column-major (order='F'), read-only
            For each neighborhood, how many of its genes belong to each study (at least 1 each)
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_bins_per_point : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's selected histogram bin count, from the occupancy search
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (n_points,), read-only
            Each reference point's lower residual-range bound (R_low), from the occupancy search
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (n_points,), read-only
            Each reference point's upper residual-range bound (R_high), from the occupancy search
            A result is a value; call `.copy()` to obtain a modifiable array.
        max_n_bins_per_point : int
            The widest `n_bins_per_point` value, `maxval(n_bins_per_point)`: the number of leading,
            meaningful rows of `pmfs`, `counts`, `mean_pmf` and `mean_pmf_counts`
        occupancy_failed : np.ndarray[np.bool_] of shape (n_points,), read-only
            `True` iff even `m_min` bins could not satisfy the occupancy criterion for this
            reference point; such a point still gets a histogram at `m_min` bins and still
            contributes to `global_js_divergence`
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_pooled_residuals : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's pooled non-NaN residual count (N_j), across all its genes of
            all studies
            A result is a value; call `.copy()` to obtain a modifiable array.
        min_bin_occupancy : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's minimum bin occupancy at `n_bins_per_point`
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_bin_occupancy : np.ndarray[np.float64] of shape (n_points,), read-only
            Each reference point's mean bin occupancy at `n_bins_per_point`
            A result is a value; call `.copy()` to obtain a modifiable array.
        max_bin_occupancy : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's maximum bin occupancy at `n_bins_per_point`
            A result is a value; call `.copy()` to obtain a modifiable array.
        sturges_bins : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's Sturges' rule diagnostic, with the rounded mean per-study
            neighbor count; never part of the occupancy search's decision
            A result is a value; call `.copy()` to obtain a modifiable array.
        fd_bins : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's Freedman-Diaconis rule diagnostic, with the rounded mean
            per-study neighbor count; never part of the occupancy search's decision
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmfs : np.ndarray[np.float64] of shape (256, n_points, n_studies,), column-major (order='F'), read-only
            `counts` normalized per reference point and study. `256` = MAX_N_BINS; only rows
            `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        counts : np.ndarray[np.int32] of shape (256, n_points, n_studies,), column-major (order='F'), read-only
            Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
            `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Count of non-NaN residuals binned per reference point, per study
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf : np.ndarray[np.float64] of shape (256, n_points,), column-major (order='F'), read-only
            The consensus pmf. `256` = MAX_N_BINS; only rows `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_counts : np.ndarray[np.int32] of shape (256, n_points,), column-major (order='F'), read-only
            Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS; only
            rows `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN residuals per reference point for the consensus pmf
            A result is a value; call `.copy()` to obtain a modifiable array.
        js_divergences : np.ndarray[np.float64] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Per-reference-point JSD of each study against the consensus pmf
            A result is a value; call `.copy()` to obtain a modifiable array.
        weights : np.ndarray[np.float64] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Per-reference-point weights for `global_js_divergence`
            A result is a value; call `.copy()` to obtain a modifiable array.
        global_js_divergence : np.ndarray[np.float64] of shape (n_studies,), read-only
            Weighted global JSD of each study against the consensus pmf
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study from the permutation test
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::run_js_comp_test_adaptive`, whose argument names are
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
    try:
        pooled_neighborhood_range = np.asfortranarray(pooled_neighborhood_range, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_neighborhood_range' must be an array of np.int32: {error}") from None
    if pooled_neighborhood_range.ndim != 2:
        raise ValueError(f"'pooled_neighborhood_range' must have 2 dimensions, but has {pooled_neighborhood_range.ndim}")

    # what the inputs already say, rather than asking for it again
    n_studies = gene_means.shape[1]
    max_n_genes_all_studies = gene_means.shape[0]
    max_n_reps_all_studies = residuals.shape[0]
    n_points = pooled_neighborhood_range.shape[1]

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
    n_neighbors_per_point = np.empty((n_studies, n_points,), dtype=np.int32, order='F')
    n_bins_per_point = np.empty((n_points,), dtype=np.int32, order='C')
    shared_residual_range_low = np.empty((n_points,), dtype=np.float64, order='C')
    shared_residual_range_high = np.empty((n_points,), dtype=np.float64, order='C')
    max_n_bins_per_point = ctypes.c_int(0)
    occupancy_failed = np.empty((n_points,), dtype=np.bool_, order='C')
    n_pooled_residuals = np.empty((n_points,), dtype=np.int32, order='C')
    min_bin_occupancy = np.empty((n_points,), dtype=np.int32, order='C')
    mean_bin_occupancy = np.empty((n_points,), dtype=np.float64, order='C')
    max_bin_occupancy = np.empty((n_points,), dtype=np.int32, order='C')
    sturges_bins = np.empty((n_points,), dtype=np.int32, order='C')
    fd_bins = np.empty((n_points,), dtype=np.int32, order='C')
    pmfs = np.empty((256, n_points, n_studies,), dtype=np.float64, order='F')
    counts = np.empty((256, n_points, n_studies,), dtype=np.int32, order='F')
    included_n_reps = np.empty((n_points, n_studies,), dtype=np.int32, order='F')
    mean_pmf = np.empty((256, n_points,), dtype=np.float64, order='F')
    mean_pmf_counts = np.empty((256, n_points,), dtype=np.int32, order='F')
    mean_pmf_included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    js_divergences = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    weights = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    global_js_divergence = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values = np.empty((n_studies,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_adaptive_c(
        ctypes.byref(ctypes.c_int(n_studies)),
        ctypes.byref(ctypes.c_int(max_n_genes_all_studies)),
        ctypes.byref(ctypes.c_int(max_n_reps_all_studies)),
        ctypes.byref(ctypes.c_int(n_points)),
        gene_means,
        residuals,
        pooled_neighborhood_range,
        n_neighbors_per_point,
        n_bins_per_point,
        shared_residual_range_low,
        shared_residual_range_high,
        ctypes.byref(max_n_bins_per_point),
        occupancy_failed,
        n_pooled_residuals,
        min_bin_occupancy,
        mean_bin_occupancy,
        max_bin_occupancy,
        sturges_bins,
        fd_bins,
        pmfs,
        counts,
        included_n_reps,
        mean_pmf,
        mean_pmf_counts,
        mean_pmf_included_n_reps,
        js_divergences,
        weights,
        global_js_divergence,
        p_values,
        ctypes.byref(ctypes.c_int(n_permutations)),
        ctypes.byref(ctypes.c_int(random_seed)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_double(gamma_occupancy)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_ADAPTIVE_ARGUMENTS, _RUN_JS_COMP_TEST_ADAPTIVE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    n_neighbors_per_point.flags.writeable = False
    n_bins_per_point.flags.writeable = False
    shared_residual_range_low.flags.writeable = False
    shared_residual_range_high.flags.writeable = False
    occupancy_failed.flags.writeable = False
    n_pooled_residuals.flags.writeable = False
    min_bin_occupancy.flags.writeable = False
    mean_bin_occupancy.flags.writeable = False
    max_bin_occupancy.flags.writeable = False
    sturges_bins.flags.writeable = False
    fd_bins.flags.writeable = False
    pmfs.flags.writeable = False
    counts.flags.writeable = False
    included_n_reps.flags.writeable = False
    mean_pmf.flags.writeable = False
    mean_pmf_counts.flags.writeable = False
    mean_pmf_included_n_reps.flags.writeable = False
    js_divergences.flags.writeable = False
    weights.flags.writeable = False
    global_js_divergence.flags.writeable = False
    p_values.flags.writeable = False

    return {
        "n_neighbors_per_point": n_neighbors_per_point,
        "n_bins_per_point": n_bins_per_point,
        "shared_residual_range_low": shared_residual_range_low,
        "shared_residual_range_high": shared_residual_range_high,
        "max_n_bins_per_point": max_n_bins_per_point.value,
        "occupancy_failed": occupancy_failed,
        "n_pooled_residuals": n_pooled_residuals,
        "min_bin_occupancy": min_bin_occupancy,
        "mean_bin_occupancy": mean_bin_occupancy,
        "max_bin_occupancy": max_bin_occupancy,
        "sturges_bins": sturges_bins,
        "fd_bins": fd_bins,
        "pmfs": pmfs,
        "counts": counts,
        "included_n_reps": included_n_reps,
        "mean_pmf": mean_pmf,
        "mean_pmf_counts": mean_pmf_counts,
        "mean_pmf_included_n_reps": mean_pmf_included_n_reps,
        "js_divergences": js_divergences,
        "weights": weights,
        "global_js_divergence": global_js_divergence,
        "p_values": p_values,
    }

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
):
    r"""Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting

    Ported from 125-stabilize-jscomp's `determine_js_comp_test_n_points_n_neighbors_helper` and
    `_alloc`, merged into one implementation now that the new `_impl` rules leave no separate
    hand-written allocation layer. Pools all studies' gene means, sorts them once
    (``sort_real_heapsort_expl_size``), generates the candidate
    grid from the gene count alone
    (:func:`tensor_omics.generate_js_comp_test_candidates`
    -- Issue #187, Step 2.7: this no longer needs the pooled residuals, since real per-neighborhood
    bin counts are decided later, in Pass B below),
    then walks it from finest to coarsest resolution: for each candidate, builds every study's
    neighborhoods and checks the first admissibility gate
    (:func:`tensor_omics.check_neighborhood_overlaps`);
    once every study passes, pools the consensus pmf and checks the second gate
    (:func:`tensor_omics.check_mean_pmf_min_counts`);
    once that passes too, seeds a confidence interval with the observed JSD, bootstraps it
    (:func:`tensor_omics.bootstrap_histogram`), and
    tests it against the running best candidate for a plateau
    (:func:`tensor_omics.check_plateau_condition`).
    `plateau_mode` picks which of that CI-overlap criterion and Issue #178's complementary
    relative-effect-size one
    (:func:`tensor_omics.check_effect_size_plateau_condition`)
    governs the stop condition; both are always computed and traced (`trace_*` below) once a
    candidate is admissible, regardless of `plateau_mode`, so a caller can compare what either
    criterion would have decided. The search stops (`exit`) the moment the SELECTED criterion's
    plateau is found -- see `plateau_mode`'s own mode table below for the accepted values.

    Issue #178 also names 2 blocking dependencies for validating the effect-size thresholds
    empirically -- the KX_FACTORS default and the Freedman-Diaconis bin-count overestimate --
    both deliberately left as-is here; see the project's JSD-Comp-Test follow-up issue. (A third
    candidate blocker, a one-sided-vs-symmetric JSD formula question, was raised in the same
    follow-up issue but confirmed by the issue's own author to be a mistake in the issue text,
    not a real discrepancy -- the code's symmetric formula is correct as written.)
    `delta_median_threshold`/`delta_max_threshold` default to the issue's own suggested (not yet
    validated) 0.05/0.10.

    When `plateau_mode` selects the effect-size criterion (`MODE_PLATEAU_EFFECT_SIZE` or
    `MODE_PLATEAU_BOTH`) and it plateaus independently of the CI-overlap criterion's own running
    "best candidate" bookkeeping, `best_candidate_index`/`best_candidate_pair_confidence_interval`
    are overridden to the candidate that actually triggered the effect-size plateau, so the
    candidate this routine returns is always the one that stopped the search.

    Ported verbatim, including the
    fallback 125 relies on: if no candidate ever plateaus, the search falls back to the FIRST
    (finest-resolution) candidate and resets `best_candidate_pair_confidence_interval` to
    `-1.0`; if the grid collapsed to a single candidate (see
    :func:`tensor_omics.generate_js_comp_test_candidates`'s
    own small-N collapse note), that one candidate is used regardless of whether it plateaued or
    even passed either gate -- the plateau machinery is bypassed entirely, exactly as 125 does.

    Per the plan's work-array translation for this routine specifically: `max_n_bins_all_candidates`
    (data-dependent, not cheaply closed-form in 125) is replaced by the fixed
    ``MAX_N_BINS`` ceiling, so every
    bin-dimensioned work array below is sized to MAX_N_BINS and sliced `(1:max_n_bins, ...)` per
    candidate, rather than carrying a separate recommend-sized dimension argument for it.
    `max_n_bins` is the widest per-point bin count Issue #187's occupancy search (Pass B below)
    chose for the current candidate, `maxval(tmp_n_bins_per_point(1:n_points))` -- it replaces
    the old single scalar `n_bins` that used to come from the global-pool Sturges/FD estimate.
    `gene_means` is passed to
    ``sort_real_heapsort_expl_size``
    as its own multi-dimensional self -- that callee declares its matching dummy with an
    explicit shape, so standard Fortran sequence association reinterprets the contiguous actual
    argument as the flat 1-D array it expects, exactly as 125's own `_alloc` layer did for the
    same call.

    Impure: calls the impure
    :func:`tensor_omics.bootstrap_histogram`. A GSL
    failure it reports is folded into `ierr` (first failure only) without aborting the search,
    matching that routine's own tolerant precedent.

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
        The way to evaluate all studies' confidence-interval overlaps for the plateau
        condition, forwarded to check_plateau_condition_impl

    max_n_points_candidate : int
        Exact upper bound on the grid's first (largest) `n_points` candidate. Issue #187's
        per-point outputs below (`n_bins_per_point`, `trace_selected_n_bins`, and the other
        jagged `trace_*` arrays) are sized by this argument, so unlike before Issue #187 it
        is no longer purely an internal sizing detail the plain wrapper can compute and
        hide -- the caller must know it up front to receive those arrays, hence JUST_INFO
        rather than AUTO here now
        It is recommended to compute this argument from the `max_n_points_candidate` output produced by :func:`tensor_omics.calc_js_comp_test_candidate_bounds`.
        The minimum valid value is `1`.
    max_n_neighbors_candidate : int
        Safe upper bound on the grid's largest `n_neighbors` candidate. Also JUST_INFO, not
        because anything returned is sized by it (nothing is), but because it comes from the
        same `calc_js_comp_test_candidate_bounds` call as `max_n_points_candidate` above --
        now that that call can no longer run automatically inside this wrapper, splitting
        this one back into an AUTO call would just be a second, redundant call to the same
        routine for no benefit
        It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by :func:`tensor_omics.calc_js_comp_test_candidate_bounds`.
        The minimum valid value is `1`.
    min_residuals_per_bin : int, optional, default 10
        Minimum count each bin of the consensus pmf must reach to pass the second
        admissibility gate. Reuses Issue #187's occupancy-search default rather than an
        independently-tunable threshold of its own: once
        :func:`tensor_omics.determine_bin_count_occupancy`
        wires real per-neighborhood bin counts in, a separate laxer threshold here would
        silently let a candidate the occupancy search already marked `occupancy_failed`
        pass this gate anyway, defeating the FAILURE-detection mechanism
        The minimum valid value is `0`.
        The default value is `10`.
    min_neighbor_overlap : float, optional, default 0.1
        Minimum fractional overlap two consecutive neighborhoods must have to pass the first
        admissibility gate
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.1`.
    succeeding_ci_overlap : float, optional, default 0.9
        Minimum fractional overlap a candidate's confidence interval must have with the
        running best, per `join_method`, to plateau
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.
    plateau_mode : str, one of 'plateau_ci_overlap' | 'plateau_effect_size' | 'plateau_both', optional, default 'plateau_ci_overlap'
        Which plateau criterion decides when the search stops

        The default value is `'plateau_ci_overlap'`.
    delta_median_threshold : float, optional, default 0.05
        Upper bound the median relative JSD change across studies must stay under for a
        transition to count toward an effect-size plateau, forwarded to
        check_effect_size_plateau_condition_impl
        The minimum valid value is `above(0.0)`.
        The default value is `0.05`.
    delta_max_threshold : float, optional, default 0.1
        Upper bound the largest relative JSD change across studies must stay under for a
        transition to count toward an effect-size plateau, forwarded to
        check_effect_size_plateau_condition_impl
        The minimum valid value is `above(0.0)`.
        The default value is `0.10`.
    delta_epsilon : float, optional, default 1e-10
        Small constant preventing division by zero when a study's previous admissible
        candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
        The minimum valid value is `above(0.0)`.
        The default value is `1.0e-10`.
    delta_min_consecutive_transitions : int, optional, default 2
        Number of consecutive qualifying transitions required to declare an effect-size
        plateau, forwarded to check_effect_size_plateau_condition_impl
        The minimum valid value is `1`.
        The default value is `2`.
    m_min : int, optional, default 3
        Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
        forwarded to determine_bin_count_occupancy_impl
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count Pass B's occupancy search will ever test (M_max),
        forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
        determine_bin_count_occupancy_impl clamps it up to `m_min` internally
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    gamma_occupancy : float, optional, default 1.25
        Geometric growth factor for Pass B's occupancy search's coarse search stage,
        forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
        advances
        The minimum valid value is `above(1.0)`.
        The default value is `1.25`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for each reference point's own lower residual-range bound,
        forwarded to determine_bin_count_occupancy_impl
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for each reference point's own upper residual-range bound,
        forwarded to determine_bin_count_occupancy_impl
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
        to bootstrap_histogram_impl itself
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.

    Returns
    -------
    dict
        with keys:

        n_points : int
            The finally chosen candidate's `n_points`
        n_neighbors : int
            The finally chosen candidate's `n_neighbors`
        n_bins_per_point : np.ndarray[np.int32] of shape (max_n_points_candidate,), read-only
            The finally chosen candidate's per-point histogram bin count, one per reference
            point (Issue #187: every neighborhood may use a different bin count). Only the
            leading `n_points` entries are meaningful, mirroring how `n_points`/`n_neighbors`
            above are the finally chosen candidate's own values
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The finally chosen candidate's per-point lower residual-range bound (R_low), one per
            reference point (Step 3: every neighborhood may use a different, asymmetric range).
            Only the leading `n_points` entries are meaningful, mirroring `n_bins_per_point`
            above; `0.0` throughout in the two genuinely-degenerate cases where Pass B
            never ran for the returned candidate (see the final three-way branch's own comments)
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The finally chosen candidate's per-point upper residual-range bound (R_high),
            mirroring `shared_residual_range_low` above in every respect
            A result is a value; call `.copy()` to obtain a modifiable array.
        best_candidate_pair_confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), read-only
            The bootstrapped JSD confidence interval for the finally chosen candidate pair;
            `-1.0` throughout only when `plateau_established` is `False` and no
            smallest-bootstrap-uncertainty candidate could be substituted either (see
            `plateau_established`)
            A result is a value; call `.copy()` to obtain a modifiable array.
        plateau_established : bool
            `True` when a real plateau was found (by whichever criterion
            `plateau_mode` selected) or the candidate grid never had more than one candidate to
            begin with. `False` when the search exhausted every admissible candidate
            without ever finding one -- Issue #178's own "report that parameter stability could
            not be established". When `False` and at least one candidate was admissible, the
            routine still returns a real (non-`-1.0`) candidate and confidence interval: the
            admissible candidate with the smallest bootstrapped uncertainty, per the issue's own
            fallback recommendation -- `plateau_established` is what distinguishes that case from
            an actual plateau, not the confidence interval's sentinel value
        trace_n_points : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `n_points`, one entry per column of the other `trace_*`
            arrays. `16` = MAX_CANDIDATE_PAIRS, written as a literal for the same reason
            candidates_n_points_n_neighbors is in generate_js_comp_test_candidates_impl
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_n_neighbors : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_global_js_divergence : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
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
            (`W_{i,t} / J_{i,t}`), denominator floored at `delta_epsilon` -- the issue's own
            formula omits this floor, but the same near-zero-JSD instability that motivates
            `delta_epsilon` in the `Delta_{i,t}` formula applies here too (a near-zero `J`
            destabilizes any ratio that divides by it, whichever candidate's `J` it is)
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study relative JSD change from the previous admissible
            candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            `-1.0` throughout at the first admissible candidate specifically (no
            predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            holds a real value
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta_median : np.ndarray[np.float64] of shape (16,), read-only
            Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            `-1.0` at the first admissible candidate, see trace_delta above
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta_max : np.ndarray[np.float64] of shape (16,), read-only
            Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            `-1.0` at the first admissible candidate, see trace_delta above
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_selected_n_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point selected histogram bin count
            (Issue #187's `M_j`), from determine_bin_count_occupancy_impl. Unlike every OTHER
            trace_* array above, whose first extent is a fixed thing like `n_studies`, this
            array's first extent is `max_n_points_candidate`, NOT `n_points`, because `n_points`
            itself varies per candidate (that is why `trace_n_points(16)` exists as its own
            array): this array is genuinely JAGGED per candidate column `t` -- only rows
            `1:trace_n_points(t)` are meaningful for that column, rows beyond that are undefined
            padding. The result-size directive below only trims the LAST extent (candidates, via
            `n_admissible_evaluated`), not this row dimension, so a Python/R caller must
            additionally slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_occupancy_failed : np.ndarray[np.bool_] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point `occupancy_failed` flag from
            determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_n_pooled_residuals : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point pooled residual count (`N_j`) from
            determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_min_bin_occupancy : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point minimum bin occupancy at
            trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_mean_bin_occupancy : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point mean bin occupancy at
            trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_max_bin_occupancy : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point maximum bin occupancy at
            trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_sturges_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point Sturges' rule bin-count diagnostic
            from determine_bin_count_occupancy_impl (never part of the occupancy search's own
            decision). Jagged per candidate column exactly as trace_selected_n_bins above --
            only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R caller
            must slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_fd_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point Freedman-Diaconis rule bin-count
            diagnostic from determine_bin_count_occupancy_impl (never part of the occupancy
            search's own decision). Jagged per candidate column exactly as trace_selected_n_bins
            above -- only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R
            caller must slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point lower residual-range bound (R_low)
            from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_shared_residual_range_high : np.ndarray[np.float64] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point upper residual-range bound (R_high)
            from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::run_js_comp_test_parameter_search`, whose argument names are
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
):
    r"""Search the ascending adaptive (k_start, k_step, k_max) candidate sequence for a stable JSD parameter setting (Issue #217)

    The adaptive-neighborhood counterpart of
    :func:`tensor_omics.run_js_comp_test_parameter_search`.
    Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
    :func:`tensor_omics.generate_adaptive_js_comp_test_candidates`
    in order. Their `k_start` shrinks, so the number of reference points that emerges ascends.
    The plateau criteria, `plateau_mode` and the fallbacks are exactly the fixed-k search's.

    **Per candidate**, in order:

    1. Grow its neighborhoods over the pooled, sorted gene means. The pooled means are sorted
    once for the whole search, with the same routine on the same input as
    :func:`tensor_omics.construct_adaptive_neighborhoods`,
    and the growth is that routine's own, with the same `tau`, `mad_distance_factor`,
    `max_pooled_residuals` and `min_study_neighbors`. So a candidate's neighborhoods are
    exactly what `construct_adaptive_neighborhoods` returns for its knobs.
    2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    reference points emerged.
    3. First admissibility gate:
    :func:`tensor_omics.check_neighborhood_overlaps`
    on the pooled ranges `[a_i, b_i]`.
    4. Per point, decode the range into per-study gene lists and run Issue #187's occupancy
    search on the pooled residuals (Pass B). Then build each study's histogram per point
    (Pass C), the consensus pmf, and apply the second gate,
    :func:`tensor_omics.check_mean_pmf_min_counts`.
    This is the pipeline of
    :func:`tensor_omics.run_js_comp_test_adaptive`.
    5. From here on everything is the fixed-k search's post-gate bookkeeping, keyed on the
    candidate index: the observed JSD and the traces, the effect-size check (before the
    bootstrap, which reuses the observed-JSD buffer as scratch), the bootstrap, the
    smallest-uncertainty tracking, the CI-overlap check, the `plateau_mode` selection, the
    effect-size override, and the stop at the first plateau.

    **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
    `candidate_n_points` and `candidate_status`, in order, until the search stops. The status is

    - ``ADAPTIVE_CANDIDATE_EVALUATED``
    (`0`): passed both gates and has a trace column;
    - ``ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED``
    (`1`): construction status too few means or too
    few residuals, so the growth could not run as specified;
    - ``ADAPTIVE_CANDIDATE_EMPTY_STUDY``
    (`2`): construction status empty study neighborhood;
    - ``ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED``
    (`3`): more than `max_n_points_candidate` points;
    - ``ADAPTIVE_CANDIDATE_OVERLAP_FAILED``
    (`4`): the first gate failed;
    - ``ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED``
    (`5`): the second gate failed.

    The checks run in that order and the first failing one is logged. `candidate_n_points` is the
    number of reference points the construction produced (0 for too few means), whatever the
    status.

    **Result**, one of three cases, as in the fixed-k search:

    1. A plateau was found, or the sequence had a single candidate: that candidate's knobs,
    `plateau_established = True`.
    2. No plateau but at least one admissible candidate: the admissible candidate with the
    smallest median confidence-interval width, `plateau_established = False`.
    3. Nothing admissible: candidate 1's knobs, confidence interval `-1.0`,
    `plateau_established = False`.

    `n_points` is the selected candidate's own logged point count. Where no admissible candidate
    backs the result (case 3, or case 1 with a single, rejected candidate), every
    `n_bins_per_point` entry up to `min(n_points, max_n_points_candidate)` is `m_min` and both
    ranges are `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
    `construct_adaptive_neighborhoods` with the returned knobs and the same construction
    options, then `run_js_comp_test_adaptive` with the same occupancy options: its
    `global_js_divergence` equals the selected candidate's `trace_global_js_divergence` column.

    **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
    place of `trace_n_neighbors`. The per-point traces are jagged per column exactly as there.

    **Memory.** Beyond the fixed-k search's per-point arrays, which are sized by
    `max_n_points_candidate` (practical capacity from
    :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`),
    the construction needs buffers over all `N = max_n_genes_all_studies * n_studies` pooled
    entries, the largest being `tmp_n_neighbors_per_point` with `n_studies * N` integers. The
    pooling buffers `tmp_pooled_residuals`/`tmp_pooled_residuals_perm` hold
    `max_n_reps_all_studies * N` values each, the size of `residuals`, because one study may
    have all its genes in one neighborhood. The histogram buffers hold `256 * max_n_points_candidate
    * n_studies` values each.

    Impure: calls the impure
    :func:`tensor_omics.bootstrap_histogram`. A GSL
    failure it reports is folded into `ierr` (first failure only) without aborting the search,
    as in the fixed-k search.

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
        The way to evaluate all studies' confidence-interval overlaps for the plateau
        condition, forwarded to check_plateau_condition_impl

    max_n_points_candidate : int
        Most reference points a candidate may emerge with; a candidate with more is rejected
        with the capacity status. The per-point outputs and traces below are sized by it, so
        the caller has to know it up front
        It is recommended to compute this argument from the `max_n_points_candidate` output produced by :func:`tensor_omics.calc_adaptive_js_comp_test_bounds`.
        The minimum valid value is `1`.
    min_residuals_per_bin : int, optional, default 10
        Minimum count every bin must reach, both in the occupancy search and in the second
        admissibility gate (the consensus pmf)
        The minimum valid value is `0`.
        The default value is `10`.
    min_neighbor_overlap : float, optional, default 0.1
        Minimum fractional overlap two consecutive neighborhoods' pooled ranges must have to
        pass the first admissibility gate
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.1`.
    succeeding_ci_overlap : float, optional, default 0.9
        Minimum fractional overlap a candidate's confidence interval must have with the
        running best, per `join_method`, to plateau
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.9`.
    plateau_mode : str, one of 'plateau_ci_overlap' | 'plateau_effect_size' | 'plateau_both', optional, default 'plateau_ci_overlap'
        Which plateau criterion decides when the search stops

        The default value is `'plateau_ci_overlap'`.
    delta_median_threshold : float, optional, default 0.05
        Upper bound the median relative JSD change across studies must stay under for a
        transition to count toward an effect-size plateau
        The minimum valid value is `above(0.0)`.
        The default value is `0.05`.
    delta_max_threshold : float, optional, default 0.1
        Upper bound the largest relative JSD change across studies must stay under for a
        transition to count toward an effect-size plateau
        The minimum valid value is `above(0.0)`.
        The default value is `0.10`.
    delta_epsilon : float, optional, default 1e-10
        Small constant preventing division by zero when a study's previous admissible JSD was
        zero
        The minimum valid value is `above(0.0)`.
        The default value is `1.0e-10`.
    delta_min_consecutive_transitions : int, optional, default 2
        Number of consecutive qualifying transitions required to declare an effect-size plateau
        The minimum valid value is `1`.
        The default value is `2`.
    m_min : int, optional, default 3
        Smallest candidate bin count the occupancy search tests (M_min)
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `3`.
    m_max : int, optional, default 120
        Largest candidate bin count the occupancy search tests (M_max); raised to `m_min`
        when smaller
        The minimum valid value is `1`.
        The maximum valid value is `MAX_N_BINS`.
        The default value is `120`.
    gamma_occupancy : float, optional, default 1.25
        Geometric growth factor of the occupancy search's coarse stage
        The minimum valid value is `above(1.0)`.
        The default value is `1.25`.
    lower_residual_range_quantile : float, optional, default 0.05
        Quantile in [0,1] for each reference point's lower residual-range bound
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.05`.
    upper_residual_range_quantile : float, optional, default 0.95
        Quantile in [0,1] for each reference point's upper residual-range bound
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
        The default value is `0.95`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
        to bootstrap_histogram_impl itself
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.
    tau : float, optional, default 0.1
        Largest relative increase of the dispersion an adaptive growth round may cause and
        still be committed, as in construct_adaptive_neighborhoods
        The minimum valid value is `0.0`.
        The default value is `0.1`.
    mad_distance_factor : float, optional, default 1.0
        Multiple of a neighborhood's median absolute deviation that the next reference
        point's target lies beyond it, as in construct_adaptive_neighborhoods
        The minimum valid value is `0.0`.
        The default value is `1.0`.
    max_pooled_residuals : int, optional, default 0
        Cap on the non-NaN residuals adaptive rounds may grow a neighborhood's pool to; 0
        for no cap, as in construct_adaptive_neighborhoods
        The minimum valid value is `0`.
        The default value is `0`.
    min_study_neighbors : int, optional, default 1
        Fewest entries of each study every neighborhood must have for a candidate to be
        admissible, as in construct_adaptive_neighborhoods
        The minimum valid value is `1`.
        The default value is `1`.

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
            The finally chosen candidate's per-point histogram bin count; only the leading
            `n_points` entries are meaningful. Where no admissible candidate backs the result,
            the leading `min(n_points, max_n_points_candidate)` entries are `m_min`, since that
            candidate's point count may exceed the capacity
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The finally chosen candidate's per-point lower residual-range bound (R_low); only the
            leading `n_points` entries are meaningful, `0.0` where no admissible candidate
            backs the result
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The finally chosen candidate's per-point upper residual-range bound (R_high),
            mirroring `shared_residual_range_low`
            A result is a value; call `.copy()` to obtain a modifiable array.
        best_candidate_confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), read-only
            The bootstrapped JSD confidence interval of the finally chosen candidate;
            `-1.0` throughout where no admissible candidate backs the result
            A result is a value; call `.copy()` to obtain a modifiable array.
        plateau_established : bool
            `True` when a plateau was found (by the criterion `plateau_mode` selects) or the
            sequence had a single candidate; `False` when the search exhausted the sequence
            without one, in which case the smallest-uncertainty admissible candidate is
            returned if there is one
        trace_k_start : np.ndarray[np.int32] of shape (16,), read-only
            Per-admissible-candidate `k_start`. `16` = MAX_CANDIDATE_PAIRS, the adaptive
            sequence's cap
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
            Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
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
            Per-admissible-candidate, per-study confidence-interval width divided by the observed
            global JSD, the denominator floored at `delta_epsilon`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta : np.ndarray[np.float64] of shape (n_studies, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-study relative JSD change from the previous admissible
            candidate; `-1.0` throughout at the first one
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta_median : np.ndarray[np.float64] of shape (16,), read-only
            Per-admissible-candidate median of trace_delta across studies; `-1.0` at the
            first one
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_delta_max : np.ndarray[np.float64] of shape (16,), read-only
            Per-admissible-candidate maximum of trace_delta across studies; `-1.0` at the
            first one
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_selected_n_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point selected histogram bin count. Jagged:
            only rows `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller
            must additionally slice `[:trace_n_points[t], t]` themselves
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_occupancy_failed : np.ndarray[np.bool_] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point `occupancy_failed` flag; jagged as
            trace_selected_n_bins
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
            Per-admissible-candidate, per-reference-point Sturges' rule diagnostic, with the
            rounded mean per-study neighbor count; jagged as trace_selected_n_bins
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_fd_bins : np.ndarray[np.int32] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic, with
            the rounded mean per-study neighbor count; jagged as trace_selected_n_bins
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
            Number of reference points every tried candidate's construction produced (0 for too
            few means), whether or not the candidate was admissible
            The first `n_candidates_tried` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        candidate_status : np.ndarray[np.int32] of shape (16,), read-only
            What happened to every tried candidate: `0` evaluated,
            `1` construction failed,
            `2` a study neighborhood below `min_study_neighbors`,
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::run_js_comp_test_adaptive_parameter_search`, whose argument names are
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
