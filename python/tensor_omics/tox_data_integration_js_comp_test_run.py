r"""tox_data_integration_js_comp_test_run

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Final Run

The JS-Comp-Test of K studies for neighborhood parameters already chosen: per-study
histograms at each reference point's bin count, each study's JSD against the consensus pmf,
their weighting into one global JSD per study, and the permutation test.
:func:`tensor_omics.run_js_comp_test` runs it
on fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each),
:func:`tensor_omics.run_js_comp_test_adaptive`
on adaptive neighborhoods built by
:func:`tensor_omics.construct_adaptive_neighborhoods`.

Python binding, generated from tox_data_integration_js_comp_test_run. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.run_js_comp_test_c.restype = None
_lib.run_js_comp_test_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
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
_RUN_JS_COMP_TEST_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

_lib.run_js_comp_test_expert_c.restype = None
_lib.run_js_comp_test_expert_c.argtypes = (
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
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
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
_RUN_JS_COMP_TEST_EXPERT_ARGUMENTS = ("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "n_points", "n_neighbors", "gene_means", "gene_means_perms", "residuals", "x_star", "neighborhood_indices", "neighborhood_range", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "max_n_bins_per_point", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy", "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins", "pmfs", "counts", "included_n_reps", "mean_pmf", "mean_pmf_counts", "mean_pmf_included_n_reps", "js_divergences", "weights", "global_js_divergence", "p_values", "p_values_observed_consensus", "tmp_neighborhood_residuals_gathered", "tmp_permutation_mean_pmf", "tmp_permutation_point_pool", "tmp_permutation_point_draw", "tmp_pooled_residuals", "tmp_pooled_residuals_perm", "tmp_bin_counts_search", "tmp_residual_to_study_perm", "tmp_study_block_ends", "tmp_study_block_sizes", "tmp_point_n_neighbors", "tmp_permutation_mean_pmf_counts", "tmp_permutation_counts", "tmp_permutation_pmfs", "tmp_permutation_js_divergences", "tmp_permutation_weights", "tmp_permutation_global_js_divergence", "n_permutations", "random_seed", "min_residuals_per_bin", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "bin_count_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_RUN_JS_COMP_TEST_EXPERT_ARGUMENT_SOURCES = ("gene_means", "gene_means", "residuals", "x_star", "neighborhood_indices", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

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

def run_js_comp_test(
        n_neighbors,
        gene_means,
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

    The plain `run_js_comp_test` sorts each study's gene means itself (heapsort, ascending, NaN
    last) and takes no permutation; `run_js_comp_test_expert` takes `gene_means_perms` from the
    caller instead. That order must still be ascending with NaN last (only the order among tied
    genes is the caller's choice); it is checked neither to be a permutation of
    `1..max_n_genes_all_studies` per column nor to be sorted.
    The allocating wrapper runs ``prepare_gene_means_perms`` first, which may handle the call and skip this one.

    Parameters
    ----------
    n_neighbors : int
        Number of neighbors per neighborhood
        The minimum valid value is `1`.
        The maximum valid value is `max_n_genes_all_studies`.
    gene_means : np.ndarray[np.float64] of shape (max_n_genes_all_studies, n_studies,), column-major (order='F')
        Per-gene mean expression values for all studies
        NaN is permitted for this value.
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
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_run::run_js_comp_test`, whose argument names are
    the ones an error message reports.

    This entry point computes `gene_means_perms`, `tmp_neighborhood_residuals_gathered`, `tmp_permutation_mean_pmf`, `tmp_permutation_point_pool`, `tmp_permutation_point_draw`, `tmp_pooled_residuals`, `tmp_pooled_residuals_perm`, `tmp_bin_counts_search`, `tmp_residual_to_study_perm`, `tmp_study_block_ends`, `tmp_study_block_sizes`, `tmp_point_n_neighbors`, `tmp_permutation_mean_pmf_counts`, `tmp_permutation_counts`, `tmp_permutation_pmfs`, `tmp_permutation_js_divergences`, `tmp_permutation_weights` and `tmp_permutation_global_js_divergence` for you and runs `prepare_gene_means_perms` first, which may answer the call outright.
    Call `run_js_comp_test_expert` to do that yourself.
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
    if residuals.shape[2] != n_studies:
        raise ValueError(f"'residuals' has {residuals.shape[2]} along axis 2, but "
            f"'gene_means' implies n_studies == {n_studies}"
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

def run_js_comp_test_expert(
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

    The plain `run_js_comp_test` sorts each study's gene means itself (heapsort, ascending, NaN
    last) and takes no permutation; `run_js_comp_test_expert` takes `gene_means_perms` from the
    caller instead. That order must still be ascending with NaN last (only the order among tied
    genes is the caller's choice); it is checked neither to be a permutation of
    `1..max_n_genes_all_studies` per column nor to be sorted.
    The allocating wrapper runs ``prepare_gene_means_perms`` first, which may handle the call and skip this one.

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
        Per study (column), a permutation of that study's `gene_means` in ascending order, NaN last.
        The plain entry point builds it for you (heapsort); the expert one takes yours, which
        must be sorted that way too (only the order among tied genes is free). It is checked
        neither to be a permutation of `1..max_n_genes_all_studies` per column nor to be sorted.
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
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_run::run_js_comp_test_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `gene_means_perms`, `tmp_neighborhood_residuals_gathered`, `tmp_permutation_mean_pmf`, `tmp_permutation_point_pool`, `tmp_permutation_point_draw`, `tmp_pooled_residuals`, `tmp_pooled_residuals_perm`, `tmp_bin_counts_search`, `tmp_residual_to_study_perm`, `tmp_study_block_ends`, `tmp_study_block_sizes`, `tmp_point_n_neighbors`, `tmp_permutation_mean_pmf_counts`, `tmp_permutation_counts`, `tmp_permutation_pmfs`, `tmp_permutation_js_divergences`, `tmp_permutation_weights` and `tmp_permutation_global_js_divergence` yourself.
    `run_js_comp_test` computes `gene_means_perms`, `tmp_neighborhood_residuals_gathered`, `tmp_permutation_mean_pmf`, `tmp_permutation_point_pool`, `tmp_permutation_point_draw`, `tmp_pooled_residuals`, `tmp_pooled_residuals_perm`, `tmp_bin_counts_search`, `tmp_residual_to_study_perm`, `tmp_study_block_ends`, `tmp_study_block_sizes`, `tmp_point_n_neighbors`, `tmp_permutation_mean_pmf_counts`, `tmp_permutation_counts`, `tmp_permutation_pmfs`, `tmp_permutation_js_divergences`, `tmp_permutation_weights` and `tmp_permutation_global_js_divergence` for you and runs `prepare_gene_means_perms` first, which may answer the call outright.
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
    tmp_neighborhood_residuals_gathered = np.empty((max_n_reps_all_studies, n_neighbors, n_points,), dtype=np.float64, order='F')
    tmp_permutation_mean_pmf = np.empty((n_points, 256,), dtype=np.float64, order='F')
    tmp_permutation_point_pool = np.empty((256,), dtype=np.int32, order='C')
    tmp_permutation_point_draw = np.empty((256,), dtype=np.int32, order='C')
    tmp_pooled_residuals = np.empty((max_n_reps_all_studies*n_neighbors*n_studies,), dtype=np.float64, order='C')
    tmp_pooled_residuals_perm = np.empty((max_n_reps_all_studies*n_neighbors*n_studies,), dtype=np.int32, order='C')
    tmp_bin_counts_search = np.empty((256,), dtype=np.int32, order='C')
    tmp_residual_to_study_perm = np.empty((max_n_reps_all_studies*n_neighbors*n_studies,), dtype=np.int32, order='C')
    tmp_study_block_ends = np.empty((n_studies,), dtype=np.int32, order='C')
    tmp_study_block_sizes = np.empty((n_studies,), dtype=np.int32, order='C')
    tmp_point_n_neighbors = np.empty((n_studies,), dtype=np.int32, order='C')
    tmp_permutation_mean_pmf_counts = np.empty((n_points, 256,), dtype=np.int32, order='F')
    tmp_permutation_counts = np.empty((n_points, 256,), dtype=np.int32, order='F')
    tmp_permutation_pmfs = np.empty((n_points, 256, n_studies,), dtype=np.float64, order='F')
    tmp_permutation_js_divergences = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    tmp_permutation_weights = np.empty((n_points, n_studies,), dtype=np.float64, order='F')
    tmp_permutation_global_js_divergence = np.empty((n_studies,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.run_js_comp_test_expert_c(
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
        tmp_neighborhood_residuals_gathered,
        tmp_permutation_mean_pmf,
        tmp_permutation_point_pool,
        tmp_permutation_point_draw,
        tmp_pooled_residuals,
        tmp_pooled_residuals_perm,
        tmp_bin_counts_search,
        tmp_residual_to_study_perm,
        tmp_study_block_ends,
        tmp_study_block_sizes,
        tmp_point_n_neighbors,
        tmp_permutation_mean_pmf_counts,
        tmp_permutation_counts,
        tmp_permutation_pmfs,
        tmp_permutation_js_divergences,
        tmp_permutation_weights,
        tmp_permutation_global_js_divergence,
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

    check_err_code(ierr.value, _RUN_JS_COMP_TEST_EXPERT_ARGUMENTS, _RUN_JS_COMP_TEST_EXPERT_ARGUMENT_SOURCES)

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
        Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
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
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_run::run_js_comp_test_adaptive`, whose argument names are
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
