"""Python binding to In multi-study omics analyses, it is often unclear whether biological replicates originating from different studies can be safely treated as sampling the same biological condition.

Generated from tox_data_integration. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.compute_gene_means_c.restype = None
_lib.compute_gene_means_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_GENE_MEANS_ARGUMENTS = ("n_genes", "n_reps", "expr", "means", "ierr",)

_lib.compute_residuals_c.restype = None
_lib.compute_residuals_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_RESIDUALS_ARGUMENTS = ("n_genes", "n_reps", "expr", "means", "resid", "ierr",)

_lib.pool_means_c.restype = None
_lib.pool_means_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_POOL_MEANS_ARGUMENTS = ("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "n_pool", "x_star", "ierr",)

_lib.pool_means_expert_c.restype = None
_lib.pool_means_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_POOL_MEANS_EXPERT_ARGUMENTS = ("pooled_means", "pooled_means_perm", "pool_size", "n_points", "n_pool", "x_star", "ierr",)

_lib.calc_neighborhood_size_c.restype = None
_lib.calc_neighborhood_size_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_NEIGHBORHOOD_SIZE_ARGUMENTS = ("n_pool", "n_points", "n_genes_S", "mean_S", "desired_size",)

_lib.construct_neighborhoods_c.restype = None
_lib.construct_neighborhoods_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CONSTRUCT_NEIGHBORHOODS_ARGUMENTS = ("n_points", "x_star", "n_genes_S", "mean_S", "n_reps_S", "resid_S", "neighborhood_residuals", "neighborhood_indices", "n_neighbors", "ierr",)

_lib.gjct_permutation_test_c.restype = None
_lib.gjct_permutation_test_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
)

#: The wrapped procedure's arguments, so an error can name one
_GJCT_PERMUTATION_TEST_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "ierr", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2",)

_lib.gjct_permutation_test_expert_c.restype = None
_lib.gjct_permutation_test_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    nullable(ctypes.POINTER(ctypes.c_int)),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
)

#: The wrapped procedure's arguments, so an error can name one
_GJCT_PERMUTATION_TEST_EXPERT_ARGUMENTS = ("neighborhood_residuals_S1_copy", "neighborhood_residuals_S2_copy", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "tmp_pool", "tmp_pmf_S1", "tmp_pmf_S2", "tmp_counts", "tmp_included_n_reps_S1", "tmp_included_n_reps_S2", "tmp_js_divergences", "tmp_weights", "ierr", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2",)

_lib.determine_shared_residual_range_expert_c.restype = None
_lib.determine_shared_residual_range_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS = ("abs_residual_pool", "abs_residual_pool_perm", "pool_size", "shared_residual_range", "ierr", "residual_range_quantile",)

_lib.determine_shared_residual_range_c.restype = None
_lib.determine_shared_residual_range_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
)

#: The wrapped procedure's arguments, so an error can name one
_DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "shared_residual_range", "ierr", "residual_range_quantile",)

_lib.build_residual_histograms_c.restype = None
_lib.build_residual_histograms_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    nullable(np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS')),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_RESIDUAL_HISTOGRAMS_ARGUMENTS = ("neighborhood_residuals", "n_reps", "n_neighbors", "n_points", "shared_residual_range", "n_bins", "counts", "pmf", "included_n_reps", "ierr", "neighbor_mask",)

_lib.compute_divergence_per_reference_point_c.restype = None
_lib.compute_divergence_per_reference_point_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENTS = ("pmf_S1", "pmf_S2", "n_points", "n_bins", "js_divergences", "ierr",)

_lib.compute_weighted_global_divergence_c.restype = None
_lib.compute_weighted_global_divergence_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENTS = ("js_divergences", "n_points", "included_n_reps_S1", "included_n_reps_S2", "global_js_divergence", "weights", "ierr",)

_lib.fjct_compute_jsd_c.restype = None
_lib.fjct_compute_jsd_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FJCT_COMPUTE_JSD_ARGUMENTS = ("family_idx", "gene_to_family_S1", "gene_to_family_S2", "n_genes_S1", "n_genes_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_genes_S1", "neighborhood_genes_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "ierr",)

_lib.fjct_compute_jsd_expert_c.restype = None
_lib.fjct_compute_jsd_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FJCT_COMPUTE_JSD_EXPERT_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "neighbor_mask_S1", "neighbor_mask_S2", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "pmf_S1", "pmf_S2", "tmp_counts", "ierr",)

_lib.fjct_compute_contribution_scores_c.restype = None
_lib.fjct_compute_contribution_scores_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FJCT_COMPUTE_CONTRIBUTION_SCORES_ARGUMENTS = ("global_js_divergences", "total_included_n_reps_per_f", "k_families", "support_weights", "contribution_scores", "ierr",)

def compute_gene_means(
        expr,
):
    r"""Compute per-gene mean expression, ignoring NaN values

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_reps, n_genes,), column-major (order='F')
        Expression matrix

    Returns
    -------
    means : np.ndarray[np.float64] of shape (n_genes,)
        Per-gene mean expression values

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::compute_gene_means`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expr = np.asfortranarray(expr, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expr' must be an array of np.float64: {error}") from None
    if expr.ndim != 2:
        raise ValueError(f"'expr' must have 2 dimensions, but has {expr.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_reps = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    means = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_gene_means_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_reps)),
        expr,
        means,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_GENE_MEANS_ARGUMENTS)

    return means

def compute_residuals(
        expr,
        means,
):
    r"""Compute signed residuals (centering by mean)

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_reps, n_genes,), column-major (order='F')
        Expression matrix containing
    means : np.ndarray[np.float64] of shape (n_genes,)
        Per-gene mean expression values

    Returns
    -------
    resid : np.ndarray[np.float64] of shape (n_reps, n_genes,), column-major (order='F')
        Matrix of signed residuals

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::compute_residuals`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expr = np.asfortranarray(expr, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expr' must be an array of np.float64: {error}") from None
    if expr.ndim != 2:
        raise ValueError(f"'expr' must have 2 dimensions, but has {expr.ndim}")
    try:
        means = np.ascontiguousarray(means, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'means' must be an array of np.float64: {error}") from None
    if means.ndim != 1:
        raise ValueError(f"'means' must have 1 dimension, but has {means.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_reps = expr.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if means.shape[0] != n_genes:
        raise ValueError(f"'means' has {means.shape[0]} along axis 0, but "
            f"'expr' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    resid = np.empty((n_reps, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.compute_residuals_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_reps)),
        expr,
        means,
        resid,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_RESIDUALS_ARGUMENTS)

    return resid

def pool_means(
        mean_S1,
        mean_S2,
        n_points,
):
    r"""Pool per-gene mean expression values across studies

    Parameters
    ----------
    mean_S1 : np.ndarray[np.float64] of shape (n_genes_S1,)
        Per-gene mean expression values
    mean_S2 : np.ndarray[np.float64] of shape (n_genes_S2,)
        Per-gene mean expression values
    n_points : int
        Number of reference points to define

    Returns
    -------
    dict
        with keys:

        n_pool : int
            Total number of included (non-NaN) pooled mean-expression values
        x_star : np.ndarray[np.float64] of shape (n_points,)
            Mean-expression reference points

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::pool_means_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        mean_S1 = np.ascontiguousarray(mean_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_S1' must be an array of np.float64: {error}") from None
    if mean_S1.ndim != 1:
        raise ValueError(f"'mean_S1' must have 1 dimension, but has {mean_S1.ndim}")
    try:
        mean_S2 = np.ascontiguousarray(mean_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_S2' must be an array of np.float64: {error}") from None
    if mean_S2.ndim != 1:
        raise ValueError(f"'mean_S2' must have 1 dimension, but has {mean_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes_S1 = mean_S1.shape[0]
    n_genes_S2 = mean_S2.shape[0]

    # outputs and work arrays, which the caller never sees
    n_pool = ctypes.c_int(0)
    x_star = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.pool_means_c(
        ctypes.byref(ctypes.c_int(n_genes_S1)),
        mean_S1,
        ctypes.byref(ctypes.c_int(n_genes_S2)),
        mean_S2,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(n_pool),
        x_star,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _POOL_MEANS_ARGUMENTS)

    return {
        "n_pool": n_pool.value,
        "x_star": x_star,
    }

def pool_means_expert(
        pooled_means,
        pooled_means_perm,
        n_points,
):
    r"""Pool per-gene mean expression values across studies (expert entry point)

    Parameters
    ----------
    pooled_means : np.ndarray[np.float64] of shape (pool_size,)
        Pooled means
    pooled_means_perm : np.ndarray[np.int32] of shape (pool_size,)
        Sorting permutation for `pooled_means`
    n_points : int
        Number of reference points to define

    Returns
    -------
    dict
        with keys:

        n_pool : int
            Total number of included (non-NaN) pooled mean-expression values
        x_star : np.ndarray[np.float64] of shape (n_points,)
            Mean-expression reference points

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::pool_means`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_means = np.ascontiguousarray(pooled_means, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_means' must be an array of np.float64: {error}") from None
    if pooled_means.ndim != 1:
        raise ValueError(f"'pooled_means' must have 1 dimension, but has {pooled_means.ndim}")
    try:
        pooled_means_perm = np.ascontiguousarray(pooled_means_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_means_perm' must be an array of np.int32: {error}") from None
    if pooled_means_perm.ndim != 1:
        raise ValueError(f"'pooled_means_perm' must have 1 dimension, but has {pooled_means_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    pool_size = pooled_means.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if pooled_means_perm.shape[0] != pool_size:
        raise ValueError(f"'pooled_means_perm' has {pooled_means_perm.shape[0]} along axis 0, but "
            f"'pooled_means' implies pool_size == {pool_size}"
        )

    # outputs and work arrays, which the caller never sees
    n_pool = ctypes.c_int(0)
    x_star = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.pool_means_expert_c(
        pooled_means,
        pooled_means_perm,
        ctypes.byref(ctypes.c_int(pool_size)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(n_pool),
        x_star,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _POOL_MEANS_EXPERT_ARGUMENTS)

    return {
        "n_pool": n_pool.value,
        "x_star": x_star,
    }

def calc_neighborhood_size(
        n_pool,
        n_points,
        mean_S,
        desired_size=1000,
):
    r"""Calculate the number of neighbors to be used for constructing neighborhoods

    Parameters
    ----------
    n_pool : int
        Total number of pooled mean-expression values across both studies
    n_points : int
        Number of reference points
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
    desired_size : int, optional, default 1000
        Optional desired neighborhood size.
        The default value is `1000_int32`.

    Returns
    -------
    n_neighbors : int
        Calculated neighborhood size

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::calc_neighborhood_size`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        mean_S = np.ascontiguousarray(mean_S, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_S' must be an array of np.float64: {error}") from None
    if mean_S.ndim != 1:
        raise ValueError(f"'mean_S' must have 1 dimension, but has {mean_S.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes_S = mean_S.shape[0]

    # outputs and work arrays, which the caller never sees
    n_neighbors = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_neighborhood_size_c(
        ctypes.byref(ctypes.c_int(n_pool)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_genes_S)),
        mean_S,
        ctypes.byref(ctypes.c_int(desired_size)),
        ctypes.byref(n_neighbors),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_NEIGHBORHOOD_SIZE_ARGUMENTS)

    return n_neighbors.value

def construct_neighborhoods(
        x_star,
        mean_S,
        resid_S,
        n_neighbors,
):
    r"""Construct neighborhood-based residual sets (kNN)

    Parameters
    ----------
    x_star : np.ndarray[np.float64] of shape (n_points,)
        Mean-expression reference points
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
    resid_S : np.ndarray[np.float64] of shape (n_reps_S, n_genes_S,), column-major (order='F')
        Matrix of signed residuals
    n_neighbors : int
        Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**

    Returns
    -------
    dict
        with keys:

        neighborhood_residuals : np.ndarray[np.float64] of shape (n_reps_S, n_neighbors, n_points,), column-major (order='F')
            Collection of residual vectors for each neighborhood
        neighborhood_indices : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F')
            Indices of selected neighborhood genes

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::construct_neighborhoods_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x_star = np.ascontiguousarray(x_star, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x_star' must be an array of np.float64: {error}") from None
    if x_star.ndim != 1:
        raise ValueError(f"'x_star' must have 1 dimension, but has {x_star.ndim}")
    try:
        mean_S = np.ascontiguousarray(mean_S, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_S' must be an array of np.float64: {error}") from None
    if mean_S.ndim != 1:
        raise ValueError(f"'mean_S' must have 1 dimension, but has {mean_S.ndim}")
    try:
        resid_S = np.asfortranarray(resid_S, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'resid_S' must be an array of np.float64: {error}") from None
    if resid_S.ndim != 2:
        raise ValueError(f"'resid_S' must have 2 dimensions, but has {resid_S.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = x_star.shape[0]
    n_genes_S = mean_S.shape[0]
    n_reps_S = resid_S.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if resid_S.shape[1] != n_genes_S:
        raise ValueError(f"'resid_S' has {resid_S.shape[1]} along axis 1, but "
            f"'mean_S' implies n_genes_S == {n_genes_S}"
        )

    # outputs and work arrays, which the caller never sees
    neighborhood_residuals = np.empty((n_reps_S, n_neighbors, n_points,), dtype=np.float64, order='F')
    neighborhood_indices = np.empty((n_neighbors, n_points,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.construct_neighborhoods_c(
        ctypes.byref(ctypes.c_int(n_points)),
        x_star,
        ctypes.byref(ctypes.c_int(n_genes_S)),
        mean_S,
        ctypes.byref(ctypes.c_int(n_reps_S)),
        resid_S,
        neighborhood_residuals,
        neighborhood_indices,
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CONSTRUCT_NEIGHBORHOODS_ARGUMENTS)

    return {
        "neighborhood_residuals": neighborhood_residuals,
        "neighborhood_indices": neighborhood_indices,
    }

def gjct_permutation_test(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        global_jsd_observed,
        n_bins,
        shared_residual_range,
        n_permutations,
        random_seed=None,
        neighbor_mask_S1=None,
        neighbor_mask_S2=None,
):
    r"""Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    global_jsd_observed : float
        Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
    n_bins : int
        Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
    shared_residual_range : float
        Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
    n_permutations : int
        Number of permutations to perform
    random_seed : int, optional
        Seed to use for shuffling
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        jsd_null : np.ndarray[np.float64] of shape (n_permutations,)
            Vector of global divergence values obtained under the null hypothesis
        p_value : float
            Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::gjct_permutation_test_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals_S1 = np.asfortranarray(neighborhood_residuals_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S1' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}")
    try:
        neighborhood_residuals_S2 = np.asfortranarray(neighborhood_residuals_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S2' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}")
    if neighbor_mask_S1 is not None:
        try:
            neighbor_mask_S1 = np.asfortranarray(neighbor_mask_S1, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S1' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S1.ndim != 2:
            raise ValueError(f"'neighbor_mask_S1' must have 2 dimensions, but has {neighbor_mask_S1.ndim}")
    if neighbor_mask_S2 is not None:
        try:
            neighbor_mask_S2 = np.asfortranarray(neighbor_mask_S2, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S2' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S2.ndim != 2:
            raise ValueError(f"'neighbor_mask_S2' must have 2 dimensions, but has {neighbor_mask_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    jsd_null = np.empty((n_permutations,), dtype=np.float64, order='C')
    p_value = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.gjct_permutation_test_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(global_jsd_observed)),
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(ctypes.c_int(n_permutations)),
        jsd_null,
        ctypes.byref(p_value),
        ctypes.byref(ierr),
        None if random_seed is None else ctypes.byref(ctypes.c_int(random_seed)),
        neighbor_mask_S1,
        neighbor_mask_S2,
    )

    check_err_code(ierr.value, _GJCT_PERMUTATION_TEST_ARGUMENTS)

    return {
        "jsd_null": jsd_null,
        "p_value": p_value.value,
    }

def gjct_permutation_test_expert(
        neighborhood_residuals_S1_copy,
        neighborhood_residuals_S2_copy,
        global_jsd_observed,
        n_bins,
        shared_residual_range,
        n_permutations,
        random_seed=None,
        neighbor_mask_S1=None,
        neighbor_mask_S2=None,
):
    r"""Estimates the permutation-test p-value (expert entry point with caller-provided work arrays)

    Parameters
    ----------
    neighborhood_residuals_S1_copy : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F'), modified in place
        Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
    neighborhood_residuals_S2_copy : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F'), modified in place
        Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
    global_jsd_observed : float
        Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
    n_bins : int
        Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
    shared_residual_range : float
        Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
    n_permutations : int
        Number of permutations to perform
    random_seed : int, optional
        Seed to use for shuffling
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        jsd_null : np.ndarray[np.float64] of shape (n_permutations,)
            Vector of global divergence values obtained under the null hypothesis
        p_value : float
            Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::gjct_permutation_test`.
    """
    # accept anything array-like, converting only when C needs it
    if not isinstance(neighborhood_residuals_S1_copy, np.ndarray) or neighborhood_residuals_S1_copy.dtype != np.float64:
        raise TypeError("'neighborhood_residuals_S1_copy' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if neighborhood_residuals_S1_copy.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S1_copy' must have 3 dimensions, but has {neighborhood_residuals_S1_copy.ndim}")
    if not neighborhood_residuals_S1_copy.flags.f_contiguous:
        raise ValueError("'neighborhood_residuals_S1_copy' is modified in place, so it must already be column-major (order='F')")
    if not isinstance(neighborhood_residuals_S2_copy, np.ndarray) or neighborhood_residuals_S2_copy.dtype != np.float64:
        raise TypeError("'neighborhood_residuals_S2_copy' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if neighborhood_residuals_S2_copy.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S2_copy' must have 3 dimensions, but has {neighborhood_residuals_S2_copy.ndim}")
    if not neighborhood_residuals_S2_copy.flags.f_contiguous:
        raise ValueError("'neighborhood_residuals_S2_copy' is modified in place, so it must already be column-major (order='F')")
    if neighbor_mask_S1 is not None:
        try:
            neighbor_mask_S1 = np.asfortranarray(neighbor_mask_S1, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S1' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S1.ndim != 2:
            raise ValueError(f"'neighbor_mask_S1' must have 2 dimensions, but has {neighbor_mask_S1.ndim}")
    if neighbor_mask_S2 is not None:
        try:
            neighbor_mask_S2 = np.asfortranarray(neighbor_mask_S2, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask_S2' must be an array of np.bool_: {error}") from None
        if neighbor_mask_S2.ndim != 2:
            raise ValueError(f"'neighbor_mask_S2' must have 2 dimensions, but has {neighbor_mask_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1_copy.shape[0]
    n_reps_S2 = neighborhood_residuals_S2_copy.shape[0]
    n_neighbors = neighborhood_residuals_S1_copy.shape[1]
    n_points = neighborhood_residuals_S1_copy.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2_copy.shape[1] != n_neighbors:
        raise ValueError(f"'neighborhood_residuals_S2_copy' has {neighborhood_residuals_S2_copy.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1_copy' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_residuals_S2_copy.shape[2] != n_points:
        raise ValueError(f"'neighborhood_residuals_S2_copy' has {neighborhood_residuals_S2_copy.shape[2]} along axis 2, but "
            f"'neighborhood_residuals_S1_copy' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    jsd_null = np.empty((n_permutations,), dtype=np.float64, order='C')
    p_value = ctypes.c_double(0)
    tmp_pool = np.empty((n_reps_S1 + n_reps_S2, n_neighbors,), dtype=np.float64, order='F')
    tmp_pmf_S1 = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    tmp_pmf_S2 = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    tmp_counts = np.empty((n_points, n_bins,), dtype=np.int32, order='F')
    tmp_included_n_reps_S1 = np.empty((n_points,), dtype=np.int32, order='C')
    tmp_included_n_reps_S2 = np.empty((n_points,), dtype=np.int32, order='C')
    tmp_js_divergences = np.empty((n_points,), dtype=np.float64, order='C')
    tmp_weights = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.gjct_permutation_test_expert_c(
        neighborhood_residuals_S1_copy,
        neighborhood_residuals_S2_copy,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(global_jsd_observed)),
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(ctypes.c_int(n_permutations)),
        jsd_null,
        ctypes.byref(p_value),
        tmp_pool,
        tmp_pmf_S1,
        tmp_pmf_S2,
        tmp_counts,
        tmp_included_n_reps_S1,
        tmp_included_n_reps_S2,
        tmp_js_divergences,
        tmp_weights,
        ctypes.byref(ierr),
        None if random_seed is None else ctypes.byref(ctypes.c_int(random_seed)),
        neighbor_mask_S1,
        neighbor_mask_S2,
    )

    check_err_code(ierr.value, _GJCT_PERMUTATION_TEST_EXPERT_ARGUMENTS)

    return {
        "jsd_null": jsd_null,
        "p_value": p_value.value,
    }

def determine_shared_residual_range_expert(
        abs_residual_pool,
        abs_residual_pool_perm,
        residual_range_quantile=95.0,
):
    r"""Computes the shared residual range [-R, R] from S1/S2 residuals (expert entry point)

    Parameters
    ----------
    abs_residual_pool : np.ndarray[np.float64] of shape (pool_size,)
        The absolute residual values of the concatenated S1,S2 residuals
    abs_residual_pool_perm : np.ndarray[np.int32] of shape (pool_size,)
        The permutation vector that sorts `abs_residual_pool`
    residual_range_quantile : float, optional, default 95.0
        Quantile for determining the residual range.
        The default value is `95.0_real64`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::determine_shared_residual_range`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        abs_residual_pool = np.ascontiguousarray(abs_residual_pool, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'abs_residual_pool' must be an array of np.float64: {error}") from None
    if abs_residual_pool.ndim != 1:
        raise ValueError(f"'abs_residual_pool' must have 1 dimension, but has {abs_residual_pool.ndim}")
    try:
        abs_residual_pool_perm = np.ascontiguousarray(abs_residual_pool_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'abs_residual_pool_perm' must be an array of np.int32: {error}") from None
    if abs_residual_pool_perm.ndim != 1:
        raise ValueError(f"'abs_residual_pool_perm' must have 1 dimension, but has {abs_residual_pool_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    pool_size = abs_residual_pool.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if abs_residual_pool_perm.shape[0] != pool_size:
        raise ValueError(f"'abs_residual_pool_perm' has {abs_residual_pool_perm.shape[0]} along axis 0, but "
            f"'abs_residual_pool' implies pool_size == {pool_size}"
        )

    # outputs and work arrays, which the caller never sees
    shared_residual_range = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.determine_shared_residual_range_expert_c(
        abs_residual_pool,
        abs_residual_pool_perm,
        ctypes.byref(ctypes.c_int(pool_size)),
        ctypes.byref(shared_residual_range),
        ctypes.byref(ierr),
        ctypes.byref(ctypes.c_double(residual_range_quantile)),
    )

    check_err_code(ierr.value, _DETERMINE_SHARED_RESIDUAL_RANGE_EXPERT_ARGUMENTS)

    return shared_residual_range.value

def determine_shared_residual_range(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        residual_range_quantile=95.0,
):
    r"""Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    residual_range_quantile : float, optional, default 95.0
        Quantile for determining the residual range.
        The default value is `95.0_real64`.

    Returns
    -------
    shared_residual_range : float
        Computed residual range (R)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::determine_shared_residual_range_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals_S1 = np.asfortranarray(neighborhood_residuals_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S1' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}")
    try:
        neighborhood_residuals_S2 = np.asfortranarray(neighborhood_residuals_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S2' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    shared_residual_range = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.determine_shared_residual_range_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(shared_residual_range),
        ctypes.byref(ierr),
        ctypes.byref(ctypes.c_double(residual_range_quantile)),
    )

    check_err_code(ierr.value, _DETERMINE_SHARED_RESIDUAL_RANGE_ARGUMENTS)

    return shared_residual_range.value

def build_residual_histograms(
        neighborhood_residuals,
        shared_residual_range,
        n_bins,
        neighbor_mask=None,
):
    r"""Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)

    Parameters
    ----------
    neighborhood_residuals : np.ndarray[np.float64] of shape (n_reps, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    shared_residual_range : float
        Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
    n_bins : int
        Number of equally sized histogram bins in range [-R,R]
    neighbor_mask : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
            Absolute counts of a residual per bin
        pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
            `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        included_n_reps : np.ndarray[np.int32] of shape (n_points,)
            Stores the count of non-NaN replicates (included ones)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::build_residual_histograms`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals = np.asfortranarray(neighborhood_residuals, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals' must be an array of np.float64: {error}") from None
    if neighborhood_residuals.ndim != 3:
        raise ValueError(f"'neighborhood_residuals' must have 3 dimensions, but has {neighborhood_residuals.ndim}")
    if neighbor_mask is not None:
        try:
            neighbor_mask = np.asfortranarray(neighbor_mask, dtype=np.bool_)
        except (TypeError, ValueError) as error:
            raise TypeError(f"'neighbor_mask' must be an array of np.bool_: {error}") from None
        if neighbor_mask.ndim != 2:
            raise ValueError(f"'neighbor_mask' must have 2 dimensions, but has {neighbor_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps = neighborhood_residuals.shape[0]
    n_neighbors = neighborhood_residuals.shape[1]
    n_points = neighborhood_residuals.shape[2]

    # outputs and work arrays, which the caller never sees
    counts = np.empty((n_points, n_bins,), dtype=np.int32, order='F')
    pmf = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    included_n_reps = np.empty((n_points,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.build_residual_histograms_c(
        neighborhood_residuals,
        ctypes.byref(ctypes.c_int(n_reps)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        ctypes.byref(ctypes.c_int(n_bins)),
        counts,
        pmf,
        included_n_reps,
        ctypes.byref(ierr),
        neighbor_mask,
    )

    check_err_code(ierr.value, _BUILD_RESIDUAL_HISTOGRAMS_ARGUMENTS)

    return {
        "counts": counts,
        "pmf": pmf,
        "included_n_reps": included_n_reps,
    }

def compute_divergence_per_reference_point(
        pmf_S1,
        pmf_S2,
):
    r"""Computes the Jensen-Shannon divergence per reference point from the histogram pmfs

    Parameters
    ----------
    pmf_S1 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
        Computed normalized histogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
    pmf_S2 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
        Computed normalized histogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2

    Returns
    -------
    js_divergences : np.ndarray[np.float64] of shape (n_points,)
        Jensen-Shannon divergence per reference point

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::compute_divergence_per_reference_point`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pmf_S1 = np.asfortranarray(pmf_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pmf_S1' must be an array of np.float64: {error}") from None
    if pmf_S1.ndim != 2:
        raise ValueError(f"'pmf_S1' must have 2 dimensions, but has {pmf_S1.ndim}")
    try:
        pmf_S2 = np.asfortranarray(pmf_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pmf_S2' must be an array of np.float64: {error}") from None
    if pmf_S2.ndim != 2:
        raise ValueError(f"'pmf_S2' must have 2 dimensions, but has {pmf_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = pmf_S1.shape[0]
    n_bins = pmf_S1.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if pmf_S2.shape[0] != n_points:
        raise ValueError(f"'pmf_S2' has {pmf_S2.shape[0]} along axis 0, but "
            f"'pmf_S1' implies n_points == {n_points}"
        )
    if pmf_S2.shape[1] != n_bins:
        raise ValueError(f"'pmf_S2' has {pmf_S2.shape[1]} along axis 1, but "
            f"'pmf_S1' implies n_bins == {n_bins}"
        )

    # outputs and work arrays, which the caller never sees
    js_divergences = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_divergence_per_reference_point_c(
        pmf_S1,
        pmf_S2,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_bins)),
        js_divergences,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_DIVERGENCE_PER_REFERENCE_POINT_ARGUMENTS)

    return js_divergences

def compute_weighted_global_divergence(
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
):
    r"""Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences

    Parameters
    ----------
    js_divergences : np.ndarray[np.float64] of shape (n_points,)
        Jensen-Shannon divergence per reference point, computed for studies S1 and S2
    included_n_reps_S1 : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
    included_n_reps_S2 : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])

    Returns
    -------
    dict
        with keys:

        global_js_divergence : float
            Weighted global Jensen-Shannon divergence
        weights : np.ndarray[np.float64] of shape (n_points,)
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::compute_weighted_global_divergence`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        js_divergences = np.ascontiguousarray(js_divergences, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'js_divergences' must be an array of np.float64: {error}") from None
    if js_divergences.ndim != 1:
        raise ValueError(f"'js_divergences' must have 1 dimension, but has {js_divergences.ndim}")
    try:
        included_n_reps_S1 = np.ascontiguousarray(included_n_reps_S1, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'included_n_reps_S1' must be an array of np.int32: {error}") from None
    if included_n_reps_S1.ndim != 1:
        raise ValueError(f"'included_n_reps_S1' must have 1 dimension, but has {included_n_reps_S1.ndim}")
    try:
        included_n_reps_S2 = np.ascontiguousarray(included_n_reps_S2, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'included_n_reps_S2' must be an array of np.int32: {error}") from None
    if included_n_reps_S2.ndim != 1:
        raise ValueError(f"'included_n_reps_S2' must have 1 dimension, but has {included_n_reps_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = js_divergences.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if included_n_reps_S1.shape[0] != n_points:
        raise ValueError(f"'included_n_reps_S1' has {included_n_reps_S1.shape[0]} along axis 0, but "
            f"'js_divergences' implies n_points == {n_points}"
        )
    if included_n_reps_S2.shape[0] != n_points:
        raise ValueError(f"'included_n_reps_S2' has {included_n_reps_S2.shape[0]} along axis 0, but "
            f"'js_divergences' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    global_js_divergence = ctypes.c_double(0)
    weights = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_weighted_global_divergence_c(
        js_divergences,
        ctypes.byref(ctypes.c_int(n_points)),
        included_n_reps_S1,
        included_n_reps_S2,
        ctypes.byref(global_js_divergence),
        weights,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_WEIGHTED_GLOBAL_DIVERGENCE_ARGUMENTS)

    return {
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
    }

def fjct_compute_jsd(
        family_idx,
        gene_to_family_S1,
        gene_to_family_S2,
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighborhood_genes_S1,
        neighborhood_genes_S2,
        n_bins,
        shared_residual_range,
):
    r"""Computes the family-level compatibility score for a single gene family (`family_idx`)

    Parameters
    ----------
    family_idx : int
        Index of the family that should be analyzed
    gene_to_family_S1 : np.ndarray[np.int32] of shape (n_genes_S1,)
        Mapping for study 1: Each index (gene) holds the index of its family
    gene_to_family_S2 : np.ndarray[np.int32] of shape (n_genes_S2,)
        Mapping for study 2: Each index (gene) holds the index of its family
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    neighborhood_genes_S1 : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F')
        Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
    neighborhood_genes_S2 : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F')
        Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
    n_bins : int
        Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
    shared_residual_range : float
        Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]

    Returns
    -------
    dict
        with keys:

        js_divergences : np.ndarray[np.float64] of shape (n_points,)
            Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        included_n_reps_S1 : np.ndarray[np.int32] of shape (n_points,)
            Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        included_n_reps_S2 : np.ndarray[np.int32] of shape (n_points,)
            Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        total_included_n_reps : int
            Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        global_js_divergence : float
            Weighted global Jensen-Shannon divergence
        weights : np.ndarray[np.float64] of shape (n_points,)
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::fjct_compute_jsd_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_to_family_S1 = np.ascontiguousarray(gene_to_family_S1, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_family_S1' must be an array of np.int32: {error}") from None
    if gene_to_family_S1.ndim != 1:
        raise ValueError(f"'gene_to_family_S1' must have 1 dimension, but has {gene_to_family_S1.ndim}")
    try:
        gene_to_family_S2 = np.ascontiguousarray(gene_to_family_S2, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_family_S2' must be an array of np.int32: {error}") from None
    if gene_to_family_S2.ndim != 1:
        raise ValueError(f"'gene_to_family_S2' must have 1 dimension, but has {gene_to_family_S2.ndim}")
    try:
        neighborhood_residuals_S1 = np.asfortranarray(neighborhood_residuals_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S1' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}")
    try:
        neighborhood_residuals_S2 = np.asfortranarray(neighborhood_residuals_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S2' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}")
    try:
        neighborhood_genes_S1 = np.asfortranarray(neighborhood_genes_S1, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_genes_S1' must be an array of np.int32: {error}") from None
    if neighborhood_genes_S1.ndim != 2:
        raise ValueError(f"'neighborhood_genes_S1' must have 2 dimensions, but has {neighborhood_genes_S1.ndim}")
    try:
        neighborhood_genes_S2 = np.asfortranarray(neighborhood_genes_S2, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_genes_S2' must be an array of np.int32: {error}") from None
    if neighborhood_genes_S2.ndim != 2:
        raise ValueError(f"'neighborhood_genes_S2' must have 2 dimensions, but has {neighborhood_genes_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes_S1 = gene_to_family_S1.shape[0]
    n_genes_S2 = gene_to_family_S2.shape[0]
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_genes_S1.shape[0] != n_neighbors:
        raise ValueError(f"'neighborhood_genes_S1' has {neighborhood_genes_S1.shape[0]} along axis 0, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_genes_S2.shape[0] != n_neighbors:
        raise ValueError(f"'neighborhood_genes_S2' has {neighborhood_genes_S2.shape[0]} along axis 0, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )
    if neighborhood_genes_S1.shape[1] != n_points:
        raise ValueError(f"'neighborhood_genes_S1' has {neighborhood_genes_S1.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )
    if neighborhood_genes_S2.shape[1] != n_points:
        raise ValueError(f"'neighborhood_genes_S2' has {neighborhood_genes_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    js_divergences = np.empty((n_points,), dtype=np.float64, order='C')
    included_n_reps_S1 = np.empty((n_points,), dtype=np.int32, order='C')
    included_n_reps_S2 = np.empty((n_points,), dtype=np.int32, order='C')
    total_included_n_reps = ctypes.c_int(0)
    global_js_divergence = ctypes.c_double(0)
    weights = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.fjct_compute_jsd_c(
        ctypes.byref(ctypes.c_int(family_idx)),
        gene_to_family_S1,
        gene_to_family_S2,
        ctypes.byref(ctypes.c_int(n_genes_S1)),
        ctypes.byref(ctypes.c_int(n_genes_S2)),
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighborhood_genes_S1,
        neighborhood_genes_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
        ctypes.byref(total_included_n_reps),
        ctypes.byref(global_js_divergence),
        weights,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FJCT_COMPUTE_JSD_ARGUMENTS)

    return {
        "js_divergences": js_divergences,
        "included_n_reps_S1": included_n_reps_S1,
        "included_n_reps_S2": included_n_reps_S2,
        "total_included_n_reps": total_included_n_reps.value,
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
    }

def fjct_compute_jsd_expert(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighbor_mask_S1,
        neighbor_mask_S2,
        n_bins,
        shared_residual_range,
):
    r"""Computes the compatibility score for a single sub-neighborhood/family (expert entry point with caller-provided masks and work arrays)

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F')
        Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F')
        Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
    n_bins : int
        Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
    shared_residual_range : float
        Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]

    Returns
    -------
    dict
        with keys:

        js_divergences : np.ndarray[np.float64] of shape (n_points,)
            Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        included_n_reps_S1 : np.ndarray[np.int32] of shape (n_points,)
            Count of non-NaN residuals (included ones) in study 1 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        included_n_reps_S2 : np.ndarray[np.int32] of shape (n_points,)
            Count of non-NaN residuals (included ones) in study 2 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        total_included_n_reps : int
            Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        global_js_divergence : float
            Weighted global Jensen-Shannon divergence
        weights : np.ndarray[np.float64] of shape (n_points,)
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        pmf_S1 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
            Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        pmf_S2 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
            Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::fjct_compute_jsd`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals_S1 = np.asfortranarray(neighborhood_residuals_S1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S1' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}")
    try:
        neighborhood_residuals_S2 = np.asfortranarray(neighborhood_residuals_S2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighborhood_residuals_S2' must be an array of np.float64: {error}") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise ValueError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}")
    try:
        neighbor_mask_S1 = np.asfortranarray(neighbor_mask_S1, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighbor_mask_S1' must be an array of np.bool_: {error}") from None
    if neighbor_mask_S1.ndim != 2:
        raise ValueError(f"'neighbor_mask_S1' must have 2 dimensions, but has {neighbor_mask_S1.ndim}")
    try:
        neighbor_mask_S2 = np.asfortranarray(neighbor_mask_S2, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'neighbor_mask_S2' must be an array of np.bool_: {error}") from None
    if neighbor_mask_S2.ndim != 2:
        raise ValueError(f"'neighbor_mask_S2' must have 2 dimensions, but has {neighbor_mask_S2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighbor_mask_S1.shape[0] != n_neighbors:
        raise ValueError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[0]} along axis 0, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighbor_mask_S2.shape[0] != n_neighbors:
        raise ValueError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[0]} along axis 0, but "
            f"'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}"
        )
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise ValueError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )
    if neighbor_mask_S1.shape[1] != n_points:
        raise ValueError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )
    if neighbor_mask_S2.shape[1] != n_points:
        raise ValueError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[1]} along axis 1, but "
            f"'neighborhood_residuals_S1' implies n_points == {n_points}"
        )

    # outputs and work arrays, which the caller never sees
    js_divergences = np.empty((n_points,), dtype=np.float64, order='C')
    included_n_reps_S1 = np.empty((n_points,), dtype=np.int32, order='C')
    included_n_reps_S2 = np.empty((n_points,), dtype=np.int32, order='C')
    total_included_n_reps = ctypes.c_int(0)
    global_js_divergence = ctypes.c_double(0)
    weights = np.empty((n_points,), dtype=np.float64, order='C')
    pmf_S1 = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    pmf_S2 = np.empty((n_points, n_bins,), dtype=np.float64, order='F')
    tmp_counts = np.empty((n_points, n_bins,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.fjct_compute_jsd_expert_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        ctypes.byref(ctypes.c_int(n_reps_S1)),
        ctypes.byref(ctypes.c_int(n_reps_S2)),
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ctypes.c_int(n_points)),
        neighbor_mask_S1,
        neighbor_mask_S2,
        ctypes.byref(ctypes.c_int(n_bins)),
        ctypes.byref(ctypes.c_double(shared_residual_range)),
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
        ctypes.byref(total_included_n_reps),
        ctypes.byref(global_js_divergence),
        weights,
        pmf_S1,
        pmf_S2,
        tmp_counts,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FJCT_COMPUTE_JSD_EXPERT_ARGUMENTS)

    return {
        "js_divergences": js_divergences,
        "included_n_reps_S1": included_n_reps_S1,
        "included_n_reps_S2": included_n_reps_S2,
        "total_included_n_reps": total_included_n_reps.value,
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
        "pmf_S1": pmf_S1,
        "pmf_S2": pmf_S2,
    }

def fjct_compute_contribution_scores(
        global_js_divergences,
        total_included_n_reps_per_f,
):
    r"""Computes the per-family/per-sub-neighborhood contribution score

    Parameters
    ----------
    global_js_divergences : np.ndarray[np.float64] of shape (k_families,)
        Per-sub-neighborhood weighted global JSD
    total_included_n_reps_per_f : np.ndarray[np.int32] of shape (k_families,)
        Per-sub-neighborhood `total_included_n_reps`

    Returns
    -------
    dict
        with keys:

        support_weights : np.ndarray[np.float64] of shape (k_families,)
            Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
        contribution_scores : np.ndarray[np.float64] of shape (k_families,)
            Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration::fjct_compute_contribution_scores`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        global_js_divergences = np.ascontiguousarray(global_js_divergences, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'global_js_divergences' must be an array of np.float64: {error}") from None
    if global_js_divergences.ndim != 1:
        raise ValueError(f"'global_js_divergences' must have 1 dimension, but has {global_js_divergences.ndim}")
    try:
        total_included_n_reps_per_f = np.ascontiguousarray(total_included_n_reps_per_f, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'total_included_n_reps_per_f' must be an array of np.int32: {error}") from None
    if total_included_n_reps_per_f.ndim != 1:
        raise ValueError(f"'total_included_n_reps_per_f' must have 1 dimension, but has {total_included_n_reps_per_f.ndim}")

    # what the inputs already say, rather than asking for it again
    k_families = global_js_divergences.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if total_included_n_reps_per_f.shape[0] != k_families:
        raise ValueError(f"'total_included_n_reps_per_f' has {total_included_n_reps_per_f.shape[0]} along axis 0, but "
            f"'global_js_divergences' implies k_families == {k_families}"
        )

    # outputs and work arrays, which the caller never sees
    support_weights = np.empty((k_families,), dtype=np.float64, order='C')
    contribution_scores = np.empty((k_families,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.fjct_compute_contribution_scores_c(
        global_js_divergences,
        total_included_n_reps_per_f,
        ctypes.byref(ctypes.c_int(k_families)),
        support_weights,
        contribution_scores,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _FJCT_COMPUTE_CONTRIBUTION_SCORES_ARGUMENTS)

    return {
        "support_weights": support_weights,
        "contribution_scores": contribution_scores,
    }
