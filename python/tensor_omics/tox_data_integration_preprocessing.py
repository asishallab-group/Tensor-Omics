"""Python binding to Generated from the kernel; do not edit -- regenerate instead.

Generated from tox_data_integration_preprocessing. Do not edit.
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

_lib.pool_means_c.restype = None
_lib.pool_means_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_POOL_MEANS_ARGUMENTS = ("pooled_means", "pool_size", "n_points", "n_pool", "x_star", "ierr",)

_lib.pool_study_means_expert_c.restype = None
_lib.pool_study_means_expert_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_POOL_STUDY_MEANS_EXPERT_ARGUMENTS = ("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "tmp_pooled_means", "tmp_pooled_means_perm", "n_pool", "x_star", "ierr",)

_lib.pool_study_means_c.restype = None
_lib.pool_study_means_c.argtypes = (
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
_POOL_STUDY_MEANS_ARGUMENTS = ("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "n_pool", "x_star", "ierr",)

_lib.construct_neighborhoods_expert_c.restype = None
_lib.construct_neighborhoods_expert_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CONSTRUCT_NEIGHBORHOODS_EXPERT_ARGUMENTS = ("n_points", "x_star", "n_genes_S", "mean_S", "n_reps_S", "resid_S", "tmp_distances", "tmp_distances_perm", "neighborhood_residuals", "neighborhood_indices", "n_neighbors", "ierr",)

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

def compute_gene_means(
        expr,
):
    r"""Compute per-gene mean expression, ignoring NaN values

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_reps, n_genes,), column-major (order='F')
        Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.

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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::compute_gene_means`.
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
        Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    means : np.ndarray[np.float64] of shape (n_genes,)
        Per-gene mean expression values; NaN where every replicate of a gene was NaN
        NaN is permitted for this value.

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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::compute_residuals`.
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

def pool_means_expert(
        pooled_means,
        pooled_means_perm,
        n_points,
):
    r"""Turn a sorted pool of per-gene mean expression values into reference points

    Parameters
    ----------
    pooled_means : np.ndarray[np.float64] of shape (pool_size,)
        Pooled means
        NaN is permitted for this value.
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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_means`.
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

def pool_means(
        pooled_means,
        n_points,
):
    r"""Turn a sorted pool of per-gene mean expression values into reference points

    Parameters
    ----------
    pooled_means : np.ndarray[np.float64] of shape (pool_size,)
        Pooled means
        NaN is permitted for this value.
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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_means_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        pooled_means = np.ascontiguousarray(pooled_means, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'pooled_means' must be an array of np.float64: {error}") from None
    if pooled_means.ndim != 1:
        raise ValueError(f"'pooled_means' must have 1 dimension, but has {pooled_means.ndim}")

    # what the inputs already say, rather than asking for it again
    pool_size = pooled_means.shape[0]

    # outputs and work arrays, which the caller never sees
    n_pool = ctypes.c_int(0)
    x_star = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.pool_means_c(
        pooled_means,
        ctypes.byref(ctypes.c_int(pool_size)),
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

def pool_study_means_expert(
        mean_S1,
        mean_S2,
        n_points,
):
    r"""Pool the per-gene mean expression values of two studies into reference points

    Parameters
    ----------
    mean_S1 : np.ndarray[np.float64] of shape (n_genes_S1,)
        Per-gene mean expression values of study S1
        NaN is permitted for this value.
    mean_S2 : np.ndarray[np.float64] of shape (n_genes_S2,)
        Per-gene mean expression values of study S2
        NaN is permitted for this value.
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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_study_means`.
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
    tmp_pooled_means = np.empty((n_genes_S1+n_genes_S2,), dtype=np.float64, order='C')
    tmp_pooled_means_perm = np.empty((n_genes_S1+n_genes_S2,), dtype=np.int32, order='C')
    n_pool = ctypes.c_int(0)
    x_star = np.empty((n_points,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.pool_study_means_expert_c(
        ctypes.byref(ctypes.c_int(n_genes_S1)),
        mean_S1,
        ctypes.byref(ctypes.c_int(n_genes_S2)),
        mean_S2,
        ctypes.byref(ctypes.c_int(n_points)),
        tmp_pooled_means,
        tmp_pooled_means_perm,
        ctypes.byref(n_pool),
        x_star,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _POOL_STUDY_MEANS_EXPERT_ARGUMENTS)

    return {
        "n_pool": n_pool.value,
        "x_star": x_star,
    }

def pool_study_means(
        mean_S1,
        mean_S2,
        n_points,
):
    r"""Pool the per-gene mean expression values of two studies into reference points

    Parameters
    ----------
    mean_S1 : np.ndarray[np.float64] of shape (n_genes_S1,)
        Per-gene mean expression values of study S1
        NaN is permitted for this value.
    mean_S2 : np.ndarray[np.float64] of shape (n_genes_S2,)
        Per-gene mean expression values of study S2
        NaN is permitted for this value.
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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_study_means_alloc`.
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

    _lib.pool_study_means_c(
        ctypes.byref(ctypes.c_int(n_genes_S1)),
        mean_S1,
        ctypes.byref(ctypes.c_int(n_genes_S2)),
        mean_S2,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(n_pool),
        x_star,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _POOL_STUDY_MEANS_ARGUMENTS)

    return {
        "n_pool": n_pool.value,
        "x_star": x_star,
    }

def construct_neighborhoods_expert(
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
        NaN is permitted for this value.
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
        NaN is permitted for this value.
    resid_S : np.ndarray[np.float64] of shape (n_reps_S, n_genes_S,), column-major (order='F')
        Matrix of signed residuals
        NaN is permitted for this value.
    n_neighbors : int
        Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
        cannot exceed the number of genes with a defined mean
        The minimum valid value is `1_int32`.
        The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
        It is recommended to compute this with
        [[tox_data_integration_preprocessing_kernel(module):calc_neighborhood_size(function)]].

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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::construct_neighborhoods`.
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
    tmp_distances = np.empty((n_genes_S,), dtype=np.float64, order='C')
    tmp_distances_perm = np.empty((n_genes_S,), dtype=np.int32, order='C')
    neighborhood_residuals = np.empty((n_reps_S, n_neighbors, n_points,), dtype=np.float64, order='F')
    neighborhood_indices = np.empty((n_neighbors, n_points,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.construct_neighborhoods_expert_c(
        ctypes.byref(ctypes.c_int(n_points)),
        x_star,
        ctypes.byref(ctypes.c_int(n_genes_S)),
        mean_S,
        ctypes.byref(ctypes.c_int(n_reps_S)),
        resid_S,
        tmp_distances,
        tmp_distances_perm,
        neighborhood_residuals,
        neighborhood_indices,
        ctypes.byref(ctypes.c_int(n_neighbors)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CONSTRUCT_NEIGHBORHOODS_EXPERT_ARGUMENTS)

    return {
        "neighborhood_residuals": neighborhood_residuals,
        "neighborhood_indices": neighborhood_indices,
    }

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
        NaN is permitted for this value.
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
        NaN is permitted for this value.
    resid_S : np.ndarray[np.float64] of shape (n_reps_S, n_genes_S,), column-major (order='F')
        Matrix of signed residuals
        NaN is permitted for this value.
    n_neighbors : int
        Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
        cannot exceed the number of genes with a defined mean
        The minimum valid value is `1_int32`.
        The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
        It is recommended to compute this with
        [[tox_data_integration_preprocessing_kernel(module):calc_neighborhood_size(function)]].

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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::construct_neighborhoods_alloc`.
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
