r"""tox_data_integration_js_comp_test

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
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_PLATEAU_CONDITION_ARGUMENTS = ("confidence_interval", "best_candidate_pair_confidence_interval", "n_studies", "best_candidate_index", "best_exceeded_ci_overlap_count", "candidate_index", "join_method", "succeeding_ci_overlap", "ci_overlap", "plateau_found", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_PLATEAU_CONDITION_ARGUMENT_SOURCES = (None, None, "confidence_interval", None, None, None, None, None, None, None, None,)

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
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "gene_means_perms", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

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
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_RUN_JS_COMP_TEST_ADAPTIVE_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "gene_means", "residuals", "pooled_neighborhood_range", "n_neighbors_per_point", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ADAPTIVE_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", "pooled_neighborhood_range", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

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
        Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
        with invalid input reported on upper_residual_range_quantile.
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
        Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
        with invalid input reported on upper_residual_range_quantile.
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
        Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
        with invalid input reported on upper_residual_range_quantile.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::determine_bin_count_stone`, whose argument names are
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
        Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
        with invalid input reported on upper_residual_range_quantile.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test::determine_bin_count_stone_expert`, whose argument names are
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

    Not checked here: the values fit 32-bit integers only while `k_start_1 <= huge(1)/4`;
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

    Instead of one fixed neighbor count, each neighborhood grows while the mean absolute residual
    of its pool stays stable, and the next reference point is placed beyond it; the number of
    reference points emerges. `k_start`, `k_step` and `k_max` count pooled entries (one entry is
    one gene of one study).

    **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
    n_studies` entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted ascending; NaN
    means take no part. A neighborhood is a contiguous run `a..b` of that order
    (`pooled_neighborhood_range`, the input of
    :func:`tensor_omics.check_neighborhood_overlaps`),
    and its per-study counts `n_neighbors_per_point` add up to `b - a + 1`.

    **Growth** around the reference point `x_star` (the mean at the seed position): entries join
    nearest first, ties to the lower mean; the dispersion `S` is the mean absolute non-NaN
    residual of the pool.

    1. The first `k_start` entries join unconditionally. With fewer than
    ``ADAPTIVE_MIN_VALID_RESIDUALS``
    (`10`) non-NaN residuals the neighborhood stops with
    ``ADAPTIVE_STOP_TOO_FEW_RESIDUALS``
    (`6`); with all residuals zero, with
    ``ADAPTIVE_STOP_ZERO_DISPERSION``
    (`5`).
    2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
    `(S_new - S_old) / S_old > tau` is discarded and growth stops
    (``ADAPTIVE_STOP_TAU``, `1`);
    otherwise it is committed. Growth also stops at `k_max` entries
    (``ADAPTIVE_STOP_K_MAX``, `2`,
    checked first) or when no entry is left
    (``ADAPTIVE_STOP_EXHAUSTED``, `4`).
    3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
    non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
    is still evaluated, and growth then stops with
    ``ADAPTIVE_STOP_RESIDUAL_CAP`` (`3`).

    **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
    sorted position `b`, the next seed is the first position after `b` whose mean is nearest
    `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
    absolute deviation of the neighborhood's pooled means. The construction ends with the
    neighborhood that reaches the largest pooled mean.

    **Construction status.**
    ``ADAPTIVE_STATUS_TOO_FEW_MEANS``
    (`1`) when fewer than `k_start` pooled means are non-NaN: nothing
    is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
    and reports, in this priority,
    ``ADAPTIVE_STATUS_TOO_FEW_RESIDUALS``
    (`3`) if any neighborhood stopped with too few residuals,
    ``ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD``
    (`2`) if any has fewer than `min_study_neighbors`
    entries of some study, else
    ``ADAPTIVE_STATUS_OK`` (`0`).

    Outputs sized by `N` are filled only for their first `n_points` reference points.

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
        Largest relative dispersion increase a growth round may cause and still be committed
        The minimum valid value is `0.0`.
        The default value is `0.1`.
    mad_distance_factor : float, optional, default 1.0
        Multiple of a neighborhood's MAD that the next reference point's target lies beyond
        the neighborhood's largest mean
        The minimum valid value is `0.0`.
        The default value is `1.0`.
    max_pooled_residuals : int, optional, default 0
        Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
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
            Mean absolute non-NaN residual of each final neighborhood; NaN without any
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

    The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    neighborhood spans of
    :func:`tensor_omics.construct_neighborhoods_ranged`
    (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    number of positions they share divided by the length of the earlier span; a reversed span
    (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    `test_`-prefixed export as a test case.

    Parameters
    ----------
    neighborhood_range : np.ndarray[np.int32] of shape (2, n_points,), column-major (order='F')
        For each reference point, the `[min_idx, max_idx]` neighborhood span
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

    The second admissibility gate of a parameter-search candidate, checked after
    :func:`tensor_omics.check_neighborhood_overlaps`.
    Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    them does not count as a failure.

    Parameters
    ----------
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
        Absolute counts of a residual per bin for the mean pmf
        The minimum valid value is `0`.
    n_bins_per_point : np.ndarray[np.int32] of shape (n_points,)
        Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
        is inspected
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
    n_bins = mean_pmf_counts.shape[1]
    n_points = mean_pmf_counts.shape[0]

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

    For each study, the overlap of the new candidate's interval with the best one so far is
    measured as a fraction of the new interval's width, and counted when it reaches
    `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    best (the best is then kept), or, after the new candidate became the best, when the count
    meets the condition `join_method` names.

    Parameters
    ----------
    confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F')
        Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
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
        How the studies' overlaps are joined into the plateau condition; the median method
        requires `count > (n_studies - 1) / 2`

    succeeding_ci_overlap : float
        Minimum fractional overlap with the best candidate's interval for a study to pass
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
        ci_overlap : np.ndarray[np.float64] of shape (n_studies,), read-only
            Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            it was on entry, relative to the new interval's width; always filled. Meaningless
            for the very first candidate, whose best interval is only a seed
            A result is a value; call `.copy()` to obtain a modifiable array.
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
    ci_overlap = np.empty((n_studies,), dtype=np.float64, order='C')
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
        ci_overlap,
        ctypes.byref(plateau_found),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CHECK_PLATEAU_CONDITION_ARGUMENTS, _CHECK_PLATEAU_CONDITION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    ci_overlap.flags.writeable = False

    return {
        "best_candidate_index": best_candidate_index.value,
        "best_exceeded_ci_overlap_count": best_exceeded_ci_overlap_count.value,
        "ci_overlap": ci_overlap,
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

    The relative-effect-size alternative to
    :func:`tensor_omics.check_plateau_condition`:
    for each study `i`, the relative change of the observed global JSD between successive
    admissible candidates,
    `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    delta_epsilon)`, is summarized by its median (`delta_median`) and maximum (`delta_max`). A
    plateau is found once both stay below their thresholds for
    `delta_min_consecutive_transitions` consecutive transitions, counted across calls in
    `n_consecutive_ok`.

    For the first admissible candidate (`has_previous = False`) there is no transition:
    `delta`, `delta_median` and `delta_max` are `-1.0` ("not computed").

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
        `False` for the first admissible candidate, which has no transition
    delta_median_threshold : float
        Upper bound on the median relative change for a transition to qualify
        The minimum valid value is `above(0.0)`.
    delta_max_threshold : float
        Upper bound on the largest relative change for a transition to qualify
        The minimum valid value is `above(0.0)`.
    delta_epsilon : float
        Small constant preventing division by zero when a study's previous JSD was zero
        The minimum valid value is `above(0.0)`.
    delta_min_consecutive_transitions : int
        Number of consecutive qualifying transitions required to declare a plateau
        The minimum valid value is `1`.
    n_consecutive_ok : int
        Running count of consecutive qualifying transitions; reset to zero when this one does
        not qualify or `has_previous` is `False`
        The minimum valid value is `0`.

    Returns
    -------
    dict
        with keys:

        n_consecutive_ok : int
            Running count of consecutive qualifying transitions; reset to zero when this one does
            not qualify or `has_previous` is `False`
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

    `mean_pmf` is the equal-weight mean of all studies' pmfs, `mean_pmf_counts` and
    `mean_pmf_included_n_reps` the sums of their counts. Known limitation: every study is
    compared against a consensus that includes itself, not a leave-one-out background.

    Parameters
    ----------
    pmfs : np.ndarray[np.float64] of shape (n_points, n_bins, n_studies,), column-major (order='F')
        Per-study probabilities of each bin per reference point, from
        :func:`tensor_omics.build_residual_histograms`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    counts : np.ndarray[np.int32] of shape (n_points, n_bins, n_studies,), column-major (order='F')
        Absolute counts of a residual per bin for `pmfs`
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates per reference point, per study; every study needs at
        least one at every reference point
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        mean_pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
            The consensus pmf, `sum(pmfs, dim=3) / n_studies`
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN replicates per reference point, summed over the studies
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F'), read-only
            Residual counts per bin of the consensus, `sum(counts, dim=3)`
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
    n_bins = pmfs.shape[1]
    n_points = pmfs.shape[0]
    n_studies = pmfs.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if counts.shape[1] != n_bins:
        raise ValueError(f"'counts' has {counts.shape[1]} along axis 1, but "
            f"'pmfs' implies n_bins == {n_bins}"
        )
    if counts.shape[0] != n_points:
        raise ValueError(f"'counts' has {counts.shape[0]} along axis 0, but "
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
    mean_pmf = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    mean_pmf_included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    mean_pmf_counts = np.empty((n_points, n_bins,), dtype=np.int32, order='F')
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

    As :func:`tensor_omics.create_mean_pmf`, for
    callers that need no counts (the bootstrap). Same known limitation: the consensus includes
    the study compared against it.

    Parameters
    ----------
    pmfs : np.ndarray[np.float64] of shape (n_points, n_bins, n_studies,), column-major (order='F')
        Per-study probabilities of each bin per reference point, from
        :func:`tensor_omics.build_residual_histograms`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    mean_pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
        The consensus pmf, `sum(pmfs, dim=3) / n_studies`
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
    n_bins = pmfs.shape[1]
    n_points = pmfs.shape[0]
    n_studies = pmfs.shape[2]

    # outputs and work arrays, which the caller never sees
    mean_pmf = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
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

    In every bootstrap, each study's histogram at each reference point is drawn with replacement
    (multinomial) from the pooled consensus counts `mean_pmf_counts`, with as many residuals as
    the study observed; the consensus pmf, each study's JSD against it and its weighted global
    JSD are recomputed from the draw.

    `confidence_interval` is input and output: its incoming `[lower, upper]` values seed every
    slot of the bottom-k/top-k heaps, so they are only displaced by more extreme bootstrap
    values. On return it holds `[largest of the n_bootstrapping_top_k_jsds smallest values,
    smallest of the n_bootstrapping_top_k_jsds largest values]`.

    Errors (`ierr`), with `confidence_interval` left unchanged: `ERR_INVALID_INPUT`, before any
    draw, when at some reference point the studies' counts do not add up to
    `mean_pmf_included_n_reps` or that total is not the sum of `mean_pmf_counts` (as
    :func:`tensor_omics.create_mean_pmf` produces
    them); `ERR_ALLOC_FAIL` when the random number generator cannot be created. A failed draw
    later on is recorded (first failure only) without stopping the bootstrap.

    Parameters
    ----------
    n_bootstraps : int
        Number of bootstrap resamples to perform
        The minimum valid value is `1`.
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
        Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
        The minimum valid value is `0`.
    mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates per reference point of the pool
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates per reference point, per study: the size of each study's
        draw. Every study needs at least one at every reference point
        The minimum valid value is `1`.
    confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), modified in place
        On entry the seed values of the heaps (see above); on return the bootstrapped
        `[lower, upper]` interval per study
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
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
    n_bins = mean_pmf_counts.shape[1]
    n_points = mean_pmf_counts.shape[0]
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
        bin_count_mode='bin_count_occupancy',
):
    r"""Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting

    The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:

    1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
    the caller, e.g. with
    :func:`tensor_omics.pool_means`) form its
    neighborhoods
    (:func:`tensor_omics.construct_neighborhoods_ranged`).
    2. Per reference point, all studies' residuals are pooled and the bin-count search
    `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
    :func:`tensor_omics.determine_bin_count_occupancy`
    by default, or Stone's rule
    :func:`tensor_omics.determine_bin_count_stone`).
    3. Each study's histogram and pmf is built per point
    (:func:`tensor_omics.build_residual_histograms`).
    4. The consensus pmf
    (:func:`tensor_omics.create_mean_pmf`), each
    study's per-point JSD against it and its weighted global JSD, and the permutation test
    (:func:`tensor_omics.gjct_permutation_test`) give
    `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
    `n_permutations = 0` both p-values are NaN.

    There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
    still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
    caller that needs such neighborhoods rejected should use
    :func:`tensor_omics.run_js_comp_test_parameter_search`.

    The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
    second extent ``MAX_N_BINS``
    (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
    zero. A Python/R caller slices `[:, :max_n_bins_per_point, ...]` (R:
    `[, seq_len(max_n_bins_per_point), ]`) itself.

    Errors (`ierr`): `ERR_INVALID_INPUT` when some study has no non-NaN residual at some
    reference point (for example a study whose gene means are all missing), since the
    equal-weight consensus needs every study; the consensus, JSD, weight and p-value outputs are
    then not filled. `ERR_ALLOC_FAIL` when the permutation test cannot create its random number
    generator; `pmfs`, `js_divergences`, `weights` and `global_js_divergence` then hold their
    values from before the permutation test. An out-of-range neighbor index is also reported
    as `ERR_INVALID_INPUT`, a safeguard valid input cannot reach.

    Parameters
    ----------
    n_neighbors : int
        Number of neighbors per neighborhood
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies`.
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
    n_permutations : int, optional, default 1000
        Number of permutations of the permutation test; `0` skips it (p-values NaN)
        The minimum valid value is `0`.
        The default value is `1000`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a bin count to be
        admissible, in every `bin_count_mode`
        The minimum valid value is `0`.
        The default value is `10`.
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

        neighborhood_indices : np.ndarray[np.int32] of shape (n_neighbors, n_points, n_studies,), column-major (order='F'), read-only
            Gene indices of the selected neighborhood, per reference point, per study
            A result is a value; call `.copy()` to obtain a modifiable array.
        neighborhood_range : np.ndarray[np.int32] of shape (2, n_points, n_studies,), column-major (order='F'), read-only
            For each reference point and study, the `[min_idx, max_idx]` neighborhood span
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_bins_per_point : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's histogram bin count `M_j`, from the bin-count search
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (n_points,), read-only
            Each reference point's lower residual-range bound (R_low), from the bin-count search
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (n_points,), read-only
            Each reference point's upper residual-range bound (R_high), from the bin-count search
            A result is a value; call `.copy()` to obtain a modifiable array.
        max_n_bins_per_point : int
            `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            `mean_pmf` and `mean_pmf_counts` (see above)
        occupancy_failed : np.ndarray[np.bool_] of shape (n_points,), read-only
            `True` iff the bin-count search failed for this point; it then uses `m_min` bins and
            still contributes (see above)
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_pooled_residuals : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's non-NaN pooled residual count (N_j)
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
            Each reference point's Sturges' rule diagnostic (not used for the decision)
            A result is a value; call `.copy()` to obtain a modifiable array.
        fd_bins : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmfs : np.ndarray[np.float64] of shape (n_points, 256, n_studies,), column-major (order='F'), read-only
            `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        counts : np.ndarray[np.int32] of shape (n_points, 256, n_studies,), column-major (order='F'), read-only
            Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Count of non-NaN replicates per reference point, per study
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf : np.ndarray[np.float64] of shape (n_points, 256,), column-major (order='F'), read-only
            The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, 256,), column-major (order='F'), read-only
            Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN replicates per reference point of the consensus pmf
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
            Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values_observed_consensus : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            comparison only; NaN when `n_permutations` is 0
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
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

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
    pmfs = np.empty((n_points, 256, n_studies,), dtype=np.float64, order='F')
    counts = np.empty((n_points, 256, n_studies,), dtype=np.int32, order='F')
    included_n_reps = np.empty((n_points, n_studies,), dtype=np.int32, order='F')
    mean_pmf = np.empty((n_points, 256,), dtype=np.float64, order='F')
    mean_pmf_counts = np.empty((n_points, 256,), dtype=np.int32, order='F')
    mean_pmf_included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    js_divergences = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    weights = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    global_js_divergence = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values_observed_consensus = np.empty((n_studies,), dtype=np.float64, order='C')
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
        p_values_observed_consensus,
        ctypes.byref(ctypes.c_int(n_permutations)),
        ctypes.byref(ctypes.c_int(random_seed)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_double(gamma_occupancy)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        bin_count_mode,
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
    p_values_observed_consensus.flags.writeable = False

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
        "p_values_observed_consensus": p_values_observed_consensus,
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
        bin_count_mode='bin_count_occupancy',
):
    r"""Run the JSD-Comp-Test pipeline on adaptive (Issue #217) neighborhoods

    The adaptive counterpart of
    :func:`tensor_omics.run_js_comp_test`: the
    same pipeline and outputs, on neighborhoods whose size varies per reference point and per
    study. Call
    :func:`tensor_omics.construct_adaptive_neighborhoods`
    first and pass its `pooled_neighborhood_range` (its first `n_points` columns) with the same
    `gene_means` and `residuals`. Each range is decoded back to per-study gene lists; their sizes
    are returned as `n_neighbors_per_point`. Neighborhoods that decode to the same gene sets as a
    fixed-k run give bit-identical results. There is no admissibility gate.

    Errors (`ierr`): `ERR_INVALID_INPUT` before any histogram work for the first range that is
    reversed, reaches past the non-NaN pooled means, or leaves a study without a gene (every
    output then undefined); otherwise as in `run_js_comp_test_impl` (a study without a non-NaN
    residual at some point, or `ERR_ALLOC_FAIL` from the permutation test).

    The pooling work arrays hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies`
    values each, the size of `residuals`, since one study may contribute all its genes to one
    neighborhood. The bin-sized arrays have the fixed second extent
    ``MAX_N_BINS``; only their first
    `max_n_bins_per_point` bins are meaningful (Python `[:, :max_n_bins_per_point, ...]`, R
    `[, seq_len(max_n_bins_per_point), ]`).

    Parameters
    ----------
    gene_means : np.ndarray[np.float64] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Mean expression of every gene in every study, NaN for a missing gene; the array the
        neighborhoods were constructed from
        NaN is permitted for this value.
    residuals : np.ndarray[np.float64] of shape (max_n_reps_all_studies, max_n_genes_all_studies, n_studies,), column-major (order='F')
        Signed residuals of every replicate of every gene in every study, NaN for a missing value
        NaN is permitted for this value.
    pooled_neighborhood_range : np.ndarray[np.int32] of shape (2, n_points,), column-major (order='F')
        For each neighborhood, its first and last sorted pooled position, from
        construct_adaptive_neighborhoods
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies*n_studies`.
    n_permutations : int, optional, default 1000
        Number of permutations of the permutation test; `0` skips it (p-values NaN)
        The minimum valid value is `0`.
        The default value is `1000`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.
    min_residuals_per_bin : int, optional, default 10
        Minimum number of pooled residuals every bin must reach for a bin count to be
        admissible, in every `bin_count_mode`
        The minimum valid value is `0`.
        The default value is `10`.
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

        n_neighbors_per_point : np.ndarray[np.int32] of shape (n_studies, n_points,), column-major (order='F'), read-only
            For each neighborhood, how many of its genes belong to each study (at least 1 each)
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_bins_per_point : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's selected histogram bin count, from the bin-count search
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (n_points,), read-only
            Each reference point's lower residual-range bound (R_low), from the bin-count search
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_high : np.ndarray[np.float64] of shape (n_points,), read-only
            Each reference point's upper residual-range bound (R_high), from the bin-count search
            A result is a value; call `.copy()` to obtain a modifiable array.
        max_n_bins_per_point : int
            `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            `mean_pmf` and `mean_pmf_counts`
        occupancy_failed : np.ndarray[np.bool_] of shape (n_points,), read-only
            `True` iff the bin-count search failed for this point; it then uses `m_min` bins and
            still contributes
            A result is a value; call `.copy()` to obtain a modifiable array.
        n_pooled_residuals : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's non-NaN pooled residual count (N_j)
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
            Each reference point's Sturges' rule diagnostic (not used for the decision)
            A result is a value; call `.copy()` to obtain a modifiable array.
        fd_bins : np.ndarray[np.int32] of shape (n_points,), read-only
            Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmfs : np.ndarray[np.float64] of shape (n_points, 256, n_studies,), column-major (order='F'), read-only
            `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        counts : np.ndarray[np.int32] of shape (n_points, 256, n_studies,), column-major (order='F'), read-only
            Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F'), read-only
            Count of non-NaN residuals binned per reference point, per study
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf : np.ndarray[np.float64] of shape (n_points, 256,), column-major (order='F'), read-only
            The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, 256,), column-major (order='F'), read-only
            Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            meaningful
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
            Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values_observed_consensus : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            comparison only; NaN when `n_permutations` is 0
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
    bin_count_mode = np.array([str(bin_count_mode).lower().encode().ljust(29)], dtype="S29")

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
    pmfs = np.empty((n_points, 256, n_studies,), dtype=np.float64, order='F')
    counts = np.empty((n_points, 256, n_studies,), dtype=np.int32, order='F')
    included_n_reps = np.empty((n_points, n_studies,), dtype=np.int32, order='F')
    mean_pmf = np.empty((n_points, 256,), dtype=np.float64, order='F')
    mean_pmf_counts = np.empty((n_points, 256,), dtype=np.int32, order='F')
    mean_pmf_included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    js_divergences = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    weights = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    global_js_divergence = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values = np.empty((n_studies,), dtype=np.float64, order='C')
    p_values_observed_consensus = np.empty((n_studies,), dtype=np.float64, order='C')
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
        p_values_observed_consensus,
        ctypes.byref(ctypes.c_int(n_permutations)),
        ctypes.byref(ctypes.c_int(random_seed)),
        ctypes.byref(ctypes.c_int(min_residuals_per_bin)),
        ctypes.byref(ctypes.c_int(m_min)),
        ctypes.byref(ctypes.c_int(m_max)),
        ctypes.byref(ctypes.c_double(gamma_occupancy)),
        ctypes.byref(ctypes.c_double(lower_residual_range_quantile)),
        ctypes.byref(ctypes.c_double(upper_residual_range_quantile)),
        bin_count_mode,
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
    p_values_observed_consensus.flags.writeable = False

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
        "p_values_observed_consensus": p_values_observed_consensus,
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
            are meaningful
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The chosen candidate's lower residual-range bound (R_low) per reference point; the
            first `n_points` entries are meaningful, `0.0` when no candidate was admissible
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
            `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            `[:trace_n_points[t], t]`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_occupancy_failed : np.ndarray[np.bool_] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point flag: `True` iff the bin-count search
            failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
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
            `min(n_points, max_n_points_candidate)` entries are `m_min`
            A result is a value; call `.copy()` to obtain a modifiable array.
        shared_residual_range_low : np.ndarray[np.float64] of shape (max_n_points_candidate,), read-only
            The chosen candidate's lower residual-range bound (R_low) per reference point; the
            first `n_points` entries are meaningful, `0.0` without an admissible candidate
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
            `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            `[:trace_n_points[t], t]`
            The first `n_admissible_evaluated` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        trace_occupancy_failed : np.ndarray[np.bool_] of shape (max_n_points_candidate, 16,), column-major (order='F'), read-only
            Per-admissible-candidate, per-reference-point flag: `True` iff the bin-count search
            failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
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
