r"""tox_data_integration_preprocessing

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing

The step that turns expression vectors into the neighborhood residuals the rest of the test
consumes: gene-wise means, the signed deviation of each replicate from them, and the
neighborhoods of reference points those residuals are grouped into so the comparison is
conditioned on expression level rather than pooled across it.

`calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own.

Python binding, generated from tox_data_integration_preprocessing. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.compute_gene_means_c.restype = None
_lib.compute_gene_means_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_GENE_MEANS_ARGUMENTS = ("n_genes", "n_reps", "expr", "means", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_GENE_MEANS_ARGUMENT_SOURCES = ("expr", "expr", None, None, None,)

_lib.compute_residuals_c.restype = None
_lib.compute_residuals_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_RESIDUALS_ARGUMENTS = ("n_genes", "n_reps", "expr", "means", "resid", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_RESIDUALS_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, None,)

_lib.pool_means_c.restype = None
_lib.pool_means_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_POOL_MEANS_ARGUMENTS = ("pooled_means", "pool_size", "n_points", "n_pool", "x_star", "ierr",)
#: For a derived argument, the one the caller passed it in
_POOL_MEANS_ARGUMENT_SOURCES = (None, "pooled_means", "x_star", None, None, None,)

_lib.pool_means_expert_c.restype = None
_lib.pool_means_expert_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_POOL_MEANS_EXPERT_ARGUMENTS = ("pooled_means", "pooled_means_perm", "pool_size", "n_points", "n_pool", "x_star", "ierr",)
#: For a derived argument, the one the caller passed it in
_POOL_MEANS_EXPERT_ARGUMENT_SOURCES = (None, None, "pooled_means", "x_star", None, None, None,)

_lib.pool_study_means_c.restype = None
_lib.pool_study_means_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_POOL_STUDY_MEANS_ARGUMENTS = ("n_genes_S1", "mean_S1", "n_genes_S2", "mean_S2", "n_points", "n_pool", "x_star", "ierr",)
#: For a derived argument, the one the caller passed it in
_POOL_STUDY_MEANS_ARGUMENT_SOURCES = ("mean_S1", None, "mean_S2", None, "x_star", None, None, None,)

_lib.construct_neighborhoods_c.restype = None
_lib.construct_neighborhoods_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CONSTRUCT_NEIGHBORHOODS_ARGUMENTS = ("n_points", "x_star", "n_genes_S", "mean_S", "n_reps_S", "resid_S", "neighborhood_residuals", "neighborhood_indices", "n_neighbors", "ierr",)
#: For a derived argument, the one the caller passed it in
_CONSTRUCT_NEIGHBORHOODS_ARGUMENT_SOURCES = ("x_star", None, "mean_S", None, "resid_S", None, None, None, "neighborhood_residuals", None,)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::compute_gene_means`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expr = _contract.reals(expr, "expr", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_reps = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    means = _np.empty((n_genes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_gene_means_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_reps)),
        expr,
        means,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_GENE_MEANS_ARGUMENTS, _COMPUTE_GENE_MEANS_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::compute_residuals`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expr = _contract.reals(expr, "expr", ndim=2, order="F")
    means = _contract.reals(means, "means", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_reps = expr.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if means.shape[0] != n_genes:
        raise _ToxArgumentError(f"'means' has {means.shape[0]} along axis 0, but 'expr' implies n_genes == {n_genes}; 'means' is declared with shape (n_genes)", argument="means")

    # outputs and work arrays, which the caller never sees
    resid = _np.empty((n_reps, n_genes,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.compute_residuals_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_reps)),
        expr,
        means,
        resid,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_RESIDUALS_ARGUMENTS, _COMPUTE_RESIDUALS_ARGUMENT_SOURCES)

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
            Mean-expression reference points
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_means`, whose argument names are
    the ones an error message reports.

    This entry point seeds `pooled_means_perm` and sorts it by `pooled_means`.
    Call `pool_means_expert` to do that yourself.
    """
    # check every input by the input contract, converting only when C needs it
    pooled_means = _contract.reals(pooled_means, "pooled_means", ndim=1)
    n_points = _contract.integer(n_points, "n_points")

    # what the inputs already say, rather than asking for it again
    pool_size = pooled_means.shape[0]

    # outputs and work arrays, which the caller never sees
    n_pool = _ctypes.c_int(0)
    x_star = _np.empty((n_points,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.pool_means_c(
        pooled_means,
        _ctypes.byref(_ctypes.c_int(pool_size)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(n_pool),
        x_star,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _POOL_MEANS_ARGUMENTS, _POOL_MEANS_ARGUMENT_SOURCES)

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
    n_points : int
        Number of reference points to define

    Returns
    -------
    dict
        with keys:

        n_pool : int
            Total number of included (non-NaN) pooled mean-expression values
        x_star : np.ndarray[np.float64] of shape (n_points,), read-only
            Mean-expression reference points
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_means_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `pooled_means_perm` yourself.
    `pool_means` seeds `pooled_means_perm` and sorts it by `pooled_means`.
    """
    # check every input by the input contract, converting only when C needs it
    pooled_means = _contract.reals(pooled_means, "pooled_means", ndim=1)
    pooled_means_perm = _contract.integers(pooled_means_perm, "pooled_means_perm", ndim=1)
    n_points = _contract.integer(n_points, "n_points")

    # what the inputs already say, rather than asking for it again
    pool_size = pooled_means.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if pooled_means_perm.shape[0] != pool_size:
        raise _ToxArgumentError(f"'pooled_means_perm' has {pooled_means_perm.shape[0]} along axis 0, but 'pooled_means' implies pool_size == {pool_size}; 'pooled_means_perm' is declared with shape (pool_size)", argument="pooled_means_perm")

    # outputs and work arrays, which the caller never sees
    n_pool = _ctypes.c_int(0)
    x_star = _np.empty((n_points,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.pool_means_expert_c(
        pooled_means,
        pooled_means_perm,
        _ctypes.byref(_ctypes.c_int(pool_size)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(n_pool),
        x_star,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _POOL_MEANS_EXPERT_ARGUMENTS, _POOL_MEANS_EXPERT_ARGUMENT_SOURCES)

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
            Mean-expression reference points
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::pool_study_means`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    mean_S1 = _contract.reals(mean_S1, "mean_S1", ndim=1)
    mean_S2 = _contract.reals(mean_S2, "mean_S2", ndim=1)
    n_points = _contract.integer(n_points, "n_points")

    # what the inputs already say, rather than asking for it again
    n_genes_S1 = mean_S1.shape[0]
    n_genes_S2 = mean_S2.shape[0]

    # outputs and work arrays, which the caller never sees
    n_pool = _ctypes.c_int(0)
    x_star = _np.empty((n_points,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.pool_study_means_c(
        _ctypes.byref(_ctypes.c_int(n_genes_S1)),
        mean_S1,
        _ctypes.byref(_ctypes.c_int(n_genes_S2)),
        mean_S2,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(n_pool),
        x_star,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _POOL_STUDY_MEANS_ARGUMENTS, _POOL_STUDY_MEANS_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing::construct_neighborhoods`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    x_star = _contract.reals(x_star, "x_star", ndim=1)
    mean_S = _contract.reals(mean_S, "mean_S", ndim=1)
    resid_S = _contract.reals(resid_S, "resid_S", ndim=2, order="F")
    n_neighbors = _contract.integer(n_neighbors, "n_neighbors")

    # what the inputs already say, rather than asking for it again
    n_points = x_star.shape[0]
    n_genes_S = mean_S.shape[0]
    n_reps_S = resid_S.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if resid_S.shape[1] != n_genes_S:
        raise _ToxArgumentError(f"'resid_S' has {resid_S.shape[1]} along axis 1, but 'mean_S' implies n_genes_S == {n_genes_S}; 'resid_S' is declared with shape (n_reps_S, n_genes_S)", argument="resid_S")

    # outputs and work arrays, which the caller never sees
    neighborhood_residuals = _np.empty((n_reps_S, n_neighbors, n_points,), dtype=_np.float64, order='F')
    neighborhood_indices = _np.empty((n_neighbors, n_points,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.construct_neighborhoods_c(
        _ctypes.byref(_ctypes.c_int(n_points)),
        x_star,
        _ctypes.byref(_ctypes.c_int(n_genes_S)),
        mean_S,
        _ctypes.byref(_ctypes.c_int(n_reps_S)),
        resid_S,
        neighborhood_residuals,
        neighborhood_indices,
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CONSTRUCT_NEIGHBORHOODS_ARGUMENTS, _CONSTRUCT_NEIGHBORHOODS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    neighborhood_residuals.flags.writeable = False
    neighborhood_indices.flags.writeable = False

    return {
        "neighborhood_residuals": neighborhood_residuals,
        "neighborhood_indices": neighborhood_indices,
    }
