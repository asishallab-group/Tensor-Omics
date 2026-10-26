r"""tox_data_integration_preprocessing

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing

The step that turns expression vectors into the neighborhood residuals the rest of the test
consumes: gene-wise means, the signed deviation of each replicate from them, and the
neighborhoods of reference points those residuals are grouped into so the comparison is
conditioned on expression level rather than pooled across it.

Besides neighborhoods of a fixed size, the adaptive construction
:func:`tensor_omics.construct_adaptive_neighborhoods`
grows each reference point's neighborhood while its residual dispersion stays stable.

`calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own, and
`gather_pooled_neighborhood_residuals` gathers a neighborhood's residuals of all studies.

Python binding, generated from tox_data_integration_preprocessing. Do not edit.
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
#: For a derived argument, the one the caller passed it in
_COMPUTE_GENE_MEANS_ARGUMENT_SOURCES = ("expr", "expr", None, None, None,)

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
#: For a derived argument, the one the caller passed it in
_COMPUTE_RESIDUALS_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, None,)

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
#: For a derived argument, the one the caller passed it in
_POOL_MEANS_ARGUMENT_SOURCES = (None, "pooled_means", "x_star", None, None, None,)

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
#: For a derived argument, the one the caller passed it in
_POOL_MEANS_EXPERT_ARGUMENT_SOURCES = (None, None, "pooled_means", "x_star", None, None, None,)

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
#: For a derived argument, the one the caller passed it in
_POOL_STUDY_MEANS_ARGUMENT_SOURCES = ("mean_S1", None, "mean_S2", None, "x_star", None, None, None,)

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
#: For a derived argument, the one the caller passed it in
_CONSTRUCT_NEIGHBORHOODS_ARGUMENT_SOURCES = ("x_star", None, "mean_S", None, "resid_S", None, None, None, "neighborhood_residuals", None,)

_lib.construct_neighborhoods_ranged_c.restype = None
_lib.construct_neighborhoods_ranged_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CONSTRUCT_NEIGHBORHOODS_RANGED_ARGUMENTS = ("n_points", "x_star", "n_genes_S", "mean_S", "n_neighbors", "neighborhood_indices", "neighborhood_range", "ierr",)
#: For a derived argument, the one the caller passed it in
_CONSTRUCT_NEIGHBORHOODS_RANGED_ARGUMENT_SOURCES = ("x_star", None, "mean_S", None, "neighborhood_indices", None, None, None,)

_lib.construct_neighborhoods_ranged_expert_c.restype = None
_lib.construct_neighborhoods_ranged_expert_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CONSTRUCT_NEIGHBORHOODS_RANGED_EXPERT_ARGUMENTS = ("n_points", "x_star", "n_genes_S", "mean_S", "mean_S_perm", "n_neighbors", "neighborhood_indices", "neighborhood_range", "ierr",)
#: For a derived argument, the one the caller passed it in
_CONSTRUCT_NEIGHBORHOODS_RANGED_EXPERT_ARGUMENT_SOURCES = ("x_star", None, "mean_S", None, None, "neighborhood_indices", None, None, None,)

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
    means : np.ndarray[np.float64] of shape (n_genes,), read-only
        Per-gene mean expression values
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::compute_gene_means`, whose argument names are
    the ones an error message reports.
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

    check_err_code(ierr.value, _COMPUTE_GENE_MEANS_ARGUMENTS, _COMPUTE_GENE_MEANS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    means.flags.writeable = False

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
    resid : np.ndarray[np.float64] of shape (n_reps, n_genes,), column-major (order='F'), read-only
        Matrix of signed residuals
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::compute_residuals`, whose argument names are
    the ones an error message reports.
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

    check_err_code(ierr.value, _COMPUTE_RESIDUALS_ARGUMENTS, _COMPUTE_RESIDUALS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    resid.flags.writeable = False

    return resid

def pool_means(
        pooled_means,
        n_points,
):
    r"""Turn a sorted pool of per-gene mean expression values into reference points

    This takes the pool already built; `pool_study_means` pools the means of two studies
    first, if that is what is at hand.

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
        x_star : np.ndarray[np.float64] of shape (n_points,), read-only
            Mean-expression reference points. If every pooled mean is NaN, every reference
            point is NaN, and `construct_neighborhoods_ranged` rejects such points.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_means`, whose argument names are
    the ones an error message reports.

    This entry point seeds `pooled_means_perm` and sorts it by `pooled_means`.
    Call `pool_means_expert` to do that yourself.
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

    check_err_code(ierr.value, _POOL_MEANS_ARGUMENTS, _POOL_MEANS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    x_star.flags.writeable = False

    return {
        "n_pool": n_pool.value,
        "x_star": x_star,
    }

def pool_means_expert(
        pooled_means,
        pooled_means_perm,
        n_points,
):
    r"""Turn a sorted pool of per-gene mean expression values into reference points

    This takes the pool already built; `pool_study_means` pools the means of two studies
    first, if that is what is at hand.

    Parameters
    ----------
    pooled_means : np.ndarray[np.float64] of shape (pool_size,)
        Pooled means
        NaN is permitted for this value.
    pooled_means_perm : np.ndarray[np.int32] of shape (pool_size,)
        Sorting permutation for `pooled_means`
        The minimum valid value is `1`.
        The maximum valid value is `pool_size`.
    n_points : int
        Number of reference points to define

    Returns
    -------
    dict
        with keys:

        n_pool : int
            Total number of included (non-NaN) pooled mean-expression values
        x_star : np.ndarray[np.float64] of shape (n_points,), read-only
            Mean-expression reference points. If every pooled mean is NaN, every reference
            point is NaN, and `construct_neighborhoods_ranged` rejects such points.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_means_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `pooled_means_perm` yourself.
    `pool_means` seeds `pooled_means_perm` and sorts it by `pooled_means`.
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

    check_err_code(ierr.value, _POOL_MEANS_EXPERT_ARGUMENTS, _POOL_MEANS_EXPERT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    x_star.flags.writeable = False

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

    Concatenates the two studies' means, sorts the pool, and turns it into reference
    points exactly as `pool_means` does.

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
        x_star : np.ndarray[np.float64] of shape (n_points,), read-only
            Mean-expression reference points. If every mean of both studies is NaN, every
            reference point is NaN, and `construct_neighborhoods_ranged` rejects such points.
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_study_means`, whose argument names are
    the ones an error message reports.
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

    check_err_code(ierr.value, _POOL_STUDY_MEANS_ARGUMENTS, _POOL_STUDY_MEANS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    x_star.flags.writeable = False

    return {
        "n_pool": n_pool.value,
        "x_star": x_star,
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
        The minimum valid value is `1`.
        The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
        It is recommended to compute this with
        :func:`tensor_omics.calc_neighborhood_size`.

    Returns
    -------
    dict
        with keys:

        neighborhood_residuals : np.ndarray[np.float64] of shape (n_reps_S, n_neighbors, n_points,), column-major (order='F'), read-only
            Collection of residual vectors for each neighborhood
            A result is a value; call `.copy()` to obtain a modifiable array.
        neighborhood_indices : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F'), read-only
            Indices of selected neighborhood genes
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::construct_neighborhoods`, whose argument names are
    the ones an error message reports.
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

    check_err_code(ierr.value, _CONSTRUCT_NEIGHBORHOODS_ARGUMENTS, _CONSTRUCT_NEIGHBORHOODS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neighborhood_residuals.flags.writeable = False
    neighborhood_indices.flags.writeable = False

    return {
        "neighborhood_residuals": neighborhood_residuals,
        "neighborhood_indices": neighborhood_indices,
    }

def construct_neighborhoods_ranged(
        x_star,
        mean_S,
        n_neighbors,
):
    r"""Construct neighborhood-based residual sets (kNN), with a `[min_idx, max_idx]` range per reference point

    Ported from 125-stabilize-jscomp's `construct_neighborhoods_helper`: a binary-search +
    two-pointer kNN construction over a pre-sorted `mean_S`, distinct from
    :func:`tensor_omics.construct_neighborhoods`'s
    distance-sort algorithm above. In addition to the neighborhood gene indices, this also
    reports each reference point's `[min_idx, max_idx]` neighborhood span (tie-extended, so
    genes tied with the span's edge value are never split from it), which a candidate
    admissibility gate elsewhere in the JSD-Comp-Test parameter search reasons about directly,
    without needing the gathered residuals.

    Parameters
    ----------
    x_star : np.ndarray[np.float64] of shape (n_points,)
        Mean-expression reference points
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
        NaN is permitted for this value.
    n_neighbors : int
        Number of neighbors to select per reference point
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        neighborhood_indices : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F'), read-only
            Indices of selected neighborhood genes per reference point.

            All indices are in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case
            `n_genes_S` is lower than `n_neighbors`, remaining indices are filled with the
            ones from `n_genes_S+1...n_neighbors` (a documented, deliberately-preserved
            limitation, ported as-is from 125-stabilize-jscomp).
            A result is a value; call `.copy()` to obtain a modifiable array.
        neighborhood_range : np.ndarray[np.int32] of shape (2, n_points,), column-major (order='F'), read-only
            For each reference point, the `[min_idx, max_idx]` of the included genes. The
            index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))`
            would be the min value. In case of duplicate means, `min_idx` points to the first
            appearance of the value and `max_idx` to the last, so even though their related
            mean value is the min/max in the neighborhood, the actual gene might not be
            included. If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::construct_neighborhoods_ranged`, whose argument names are
    the ones an error message reports.

    This entry point seeds `mean_S_perm` and sorts it by `mean_S`.
    Call `construct_neighborhoods_ranged_expert` to do that yourself.
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

    # what the inputs already say, rather than asking for it again
    n_points = x_star.shape[0]
    n_genes_S = mean_S.shape[0]

    # outputs and work arrays, which the caller never sees
    neighborhood_indices = np.empty((n_neighbors, n_points,), dtype=np.int32, order='F')
    neighborhood_range = np.empty((2, n_points,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.construct_neighborhoods_ranged_c(
        ctypes.byref(ctypes.c_int(n_points)),
        x_star,
        ctypes.byref(ctypes.c_int(n_genes_S)),
        mean_S,
        ctypes.byref(ctypes.c_int(n_neighbors)),
        neighborhood_indices,
        neighborhood_range,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CONSTRUCT_NEIGHBORHOODS_RANGED_ARGUMENTS, _CONSTRUCT_NEIGHBORHOODS_RANGED_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neighborhood_indices.flags.writeable = False
    neighborhood_range.flags.writeable = False

    return {
        "neighborhood_indices": neighborhood_indices,
        "neighborhood_range": neighborhood_range,
    }

def construct_neighborhoods_ranged_expert(
        x_star,
        mean_S,
        mean_S_perm,
        n_neighbors,
):
    r"""Construct neighborhood-based residual sets (kNN), with a `[min_idx, max_idx]` range per reference point

    Ported from 125-stabilize-jscomp's `construct_neighborhoods_helper`: a binary-search +
    two-pointer kNN construction over a pre-sorted `mean_S`, distinct from
    :func:`tensor_omics.construct_neighborhoods`'s
    distance-sort algorithm above. In addition to the neighborhood gene indices, this also
    reports each reference point's `[min_idx, max_idx]` neighborhood span (tie-extended, so
    genes tied with the span's edge value are never split from it), which a candidate
    admissibility gate elsewhere in the JSD-Comp-Test parameter search reasons about directly,
    without needing the gathered residuals.

    Parameters
    ----------
    x_star : np.ndarray[np.float64] of shape (n_points,)
        Mean-expression reference points
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
        NaN is permitted for this value.
    mean_S_perm : np.ndarray[np.int32] of shape (n_genes_S,)
        Sorting permutation for `mean_S`
        The minimum valid value is `1`.
        The maximum valid value is `n_genes_S`.
    n_neighbors : int
        Number of neighbors to select per reference point
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        neighborhood_indices : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F'), read-only
            Indices of selected neighborhood genes per reference point.

            All indices are in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case
            `n_genes_S` is lower than `n_neighbors`, remaining indices are filled with the
            ones from `n_genes_S+1...n_neighbors` (a documented, deliberately-preserved
            limitation, ported as-is from 125-stabilize-jscomp).
            A result is a value; call `.copy()` to obtain a modifiable array.
        neighborhood_range : np.ndarray[np.int32] of shape (2, n_points,), column-major (order='F'), read-only
            For each reference point, the `[min_idx, max_idx]` of the included genes. The
            index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))`
            would be the min value. In case of duplicate means, `min_idx` points to the first
            appearance of the value and `max_idx` to the last, so even though their related
            mean value is the min/max in the neighborhood, the actual gene might not be
            included. If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::construct_neighborhoods_ranged_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `mean_S_perm` yourself.
    `construct_neighborhoods_ranged` seeds `mean_S_perm` and sorts it by `mean_S`.
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
        mean_S_perm = np.ascontiguousarray(mean_S_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_S_perm' must be an array of np.int32: {error}") from None
    if mean_S_perm.ndim != 1:
        raise ValueError(f"'mean_S_perm' must have 1 dimension, but has {mean_S_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = x_star.shape[0]
    n_genes_S = mean_S.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if mean_S_perm.shape[0] != n_genes_S:
        raise ValueError(f"'mean_S_perm' has {mean_S_perm.shape[0]} along axis 0, but "
            f"'mean_S' implies n_genes_S == {n_genes_S}"
        )

    # outputs and work arrays, which the caller never sees
    neighborhood_indices = np.empty((n_neighbors, n_points,), dtype=np.int32, order='F')
    neighborhood_range = np.empty((2, n_points,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.construct_neighborhoods_ranged_expert_c(
        ctypes.byref(ctypes.c_int(n_points)),
        x_star,
        ctypes.byref(ctypes.c_int(n_genes_S)),
        mean_S,
        mean_S_perm,
        ctypes.byref(ctypes.c_int(n_neighbors)),
        neighborhood_indices,
        neighborhood_range,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CONSTRUCT_NEIGHBORHOODS_RANGED_EXPERT_ARGUMENTS, _CONSTRUCT_NEIGHBORHOODS_RANGED_EXPERT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neighborhood_indices.flags.writeable = False
    neighborhood_range.flags.writeable = False

    return {
        "neighborhood_indices": neighborhood_indices,
        "neighborhood_range": neighborhood_range,
    }

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
    Generated from the Fortran procedure `tox_data_integration_preprocessing::construct_adaptive_neighborhoods`, whose argument names are
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
