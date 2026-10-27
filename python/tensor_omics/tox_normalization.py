r"""tox_normalization

Normalization of expression data: putting samples on a comparable scale before anything
distance-based is computed from them.

Several schemes, each its own entry point rather than a mode: a log2 transformation, scaling
by standard deviation (LOESS-smoothed against expression level, so the correction follows the
mean-variance trend rather than assuming one), root-mean-square scaling, quantile
normalization, and normalization to unit length.

`normalization_pipeline` chains them in the order the analysis expects, and is what most
callers want. The individual routines are published for a caller assembling their own.

Python binding, generated from tox_normalization. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.normalize_unit_length_c.restype = None
_lib.normalize_unit_length_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_UNIT_LENGTH_ARGUMENTS = ("vector", "n_dims", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZE_UNIT_LENGTH_ARGUMENT_SOURCES = (None, "vector", None,)

_lib.normalization_pipeline_c.restype = None
_lib.normalization_pipeline_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZATION_PIPELINE_ARGUMENTS = ("n_genes", "n_replicates", "expr", "log_transformed_expr", "reps_per_tissue", "n_tissues", "span", "degree", "use_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZATION_PIPELINE_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, "log_transformed_expr", None, None, None, None,)

_lib.normalize_by_std_dev_c.restype = None
_lib.normalize_by_std_dev_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_BY_STD_DEV_ARGUMENTS = ("n_genes", "n_replicates", "expr", "normalized_expr", "span", "degree", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZE_BY_STD_DEV_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, None, None,)

_lib.root_mean_sq_normalization_c.restype = None
_lib.root_mean_sq_normalization_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ROOT_MEAN_SQ_NORMALIZATION_ARGUMENTS = ("n_genes", "n_replicates", "expr", "normalized_expr", "ierr",)
#: For a derived argument, the one the caller passed it in
_ROOT_MEAN_SQ_NORMALIZATION_ARGUMENT_SOURCES = ("expr", "expr", None, None, None,)

_lib.quantile_normalization_c.restype = None
_lib.quantile_normalization_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_QUANTILE_NORMALIZATION_ARGUMENTS = ("n_genes", "n_replicates", "expr", "normalized_expr", "rank_means", "ierr",)
#: For a derived argument, the one the caller passed it in
_QUANTILE_NORMALIZATION_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, None,)

_lib.log2_transformation_c.restype = None
_lib.log2_transformation_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOG2_TRANSFORMATION_ARGUMENTS = ("n_genes", "n_tissues", "expr", "transformed_expr", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOG2_TRANSFORMATION_ARGUMENT_SOURCES = ("expr", "expr", None, None, None,)

_lib.calc_tiss_avg_c.restype = None
_lib.calc_tiss_avg_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_TISS_AVG_ARGUMENTS = ("n_genes", "n_tissues", "reps_per_tissue", "expr", "tissue_averages", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_TISS_AVG_ARGUMENT_SOURCES = ("expr", "reps_per_tissue", None, None, None, None,)

_lib.calc_fchange_c.restype = None
_lib.calc_fchange_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_FCHANGE_ARGUMENTS = ("n_genes", "n_tissues", "n_pairs", "control_tissues", "condition_tissues", "expr", "fold_changes", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_FCHANGE_ARGUMENT_SOURCES = ("expr", "expr", "control_tissues", None, None, None, None, None,)

def normalize_unit_length(
        vector,
):
    r"""Normalizes an input vector to unit length in-place

    Parameters
    ----------
    vector : np.ndarray[np.float64] of shape (n_dims,), modified in place
        Vector that will be normalized to unit length
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    None

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::normalize_unit_length`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    if not isinstance(vector, _np.ndarray) or vector.dtype != _np.float64:
        raise _ToxArgumentError("'vector' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="vector")
    if vector.ndim != 1:
        raise _ToxArgumentError(f"'vector' must have 1 dimension, but has {vector.ndim}", argument="vector")

    # what the inputs already say, rather than asking for it again
    n_dims = vector.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.normalize_unit_length_c(
        vector,
        _ctypes.byref(_ctypes.c_int(n_dims)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _NORMALIZE_UNIT_LENGTH_ARGUMENTS, _NORMALIZE_UNIT_LENGTH_ARGUMENT_SOURCES)

    return None

def normalization_pipeline(
        expr,
        reps_per_tissue,
        span=0.7,
        degree=2,
        use_quantile=False,
):
    r"""Complete normalization pipeline for gene expression data.

    Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
        Gene Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    reps_per_tissue : np.ndarray[np.int32] of shape (n_tissues,)
        Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
        e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
    span : float, optional, default 0.7
        LOESS span parameter.
        The default value is `0.7`.
    degree : int, optional, default 2
        LOESS degree parameter.
        The default value is `2`.
    use_quantile : bool, optional, default False
        Use quantile normalization.
        The default value is `False`.

    Returns
    -------
    log_transformed_expr : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F'), read-only
        Log-transformed grouped `expr`
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
    Generated from the Fortran procedure `tox_normalization::normalization_pipeline`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expr = _np.asfortranarray(expr, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expr' must be an array of np.float64: {_error}", argument="expr") from None
    if expr.ndim != 2:
        raise _ToxArgumentError(f"'expr' must have 2 dimensions, but has {expr.ndim}", argument="expr")
    try:
        reps_per_tissue = _np.ascontiguousarray(reps_per_tissue, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'reps_per_tissue' must be an array of np.int32: {_error}", argument="reps_per_tissue") from None
    if reps_per_tissue.ndim != 1:
        raise _ToxArgumentError(f"'reps_per_tissue' must have 1 dimension, but has {reps_per_tissue.ndim}", argument="reps_per_tissue")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_replicates = expr.shape[0]
    n_tissues = reps_per_tissue.shape[0]

    # outputs and work arrays, which the caller never sees
    log_transformed_expr = _np.empty((n_tissues, n_genes,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.normalization_pipeline_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_replicates)),
        expr,
        log_transformed_expr,
        reps_per_tissue,
        _ctypes.byref(_ctypes.c_int(n_tissues)),
        _ctypes.byref(_ctypes.c_double(span)),
        _ctypes.byref(_ctypes.c_int(degree)),
        _ctypes.byref(_ctypes.c_bool(use_quantile)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _NORMALIZATION_PIPELINE_ARGUMENTS, _NORMALIZATION_PIPELINE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    log_transformed_expr.flags.writeable = False

    return log_transformed_expr

def normalize_by_std_dev(
        expr,
        span=0.7,
        degree=2,
):
    r"""Normalizes each gene's expression vector using LOESS-stabilized standard deviation.

    This procedure applies a global stabilization based on the relationship between
    gene-wise mean expression and empirical standard deviation.

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
        Gene Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    span : float, optional, default 0.7
        LOESS span parameter.
        The default value is `0.7`.
    degree : int, optional, default 2
        LOESS degree parameter.
        The default value is `2`.

    Returns
    -------
    normalized_expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F'), read-only
        Normalized `expr`
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
    Generated from the Fortran procedure `tox_normalization::normalize_by_std_dev`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expr = _np.asfortranarray(expr, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expr' must be an array of np.float64: {_error}", argument="expr") from None
    if expr.ndim != 2:
        raise _ToxArgumentError(f"'expr' must have 2 dimensions, but has {expr.ndim}", argument="expr")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_replicates = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    normalized_expr = _np.empty((n_replicates, n_genes,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.normalize_by_std_dev_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_replicates)),
        expr,
        normalized_expr,
        _ctypes.byref(_ctypes.c_double(span)),
        _ctypes.byref(_ctypes.c_int(degree)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _NORMALIZE_BY_STD_DEV_ARGUMENTS, _NORMALIZE_BY_STD_DEV_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    normalized_expr.flags.writeable = False

    return normalized_expr

def root_mean_sq_normalization(
        expr,
):
    r"""Normalizes each gene's expression vector using `sqrt(mean(x^2))`

    across tissues (not classical standard deviation).

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
        Gene Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    normalized_expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F'), read-only
        Normalized `expr`
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
    Generated from the Fortran procedure `tox_normalization::root_mean_sq_normalization`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expr = _np.asfortranarray(expr, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expr' must be an array of np.float64: {_error}", argument="expr") from None
    if expr.ndim != 2:
        raise _ToxArgumentError(f"'expr' must have 2 dimensions, but has {expr.ndim}", argument="expr")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_replicates = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    normalized_expr = _np.empty((n_replicates, n_genes,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.root_mean_sq_normalization_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_replicates)),
        expr,
        normalized_expr,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _ROOT_MEAN_SQ_NORMALIZATION_ARGUMENTS, _ROOT_MEAN_SQ_NORMALIZATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    normalized_expr.flags.writeable = False

    return normalized_expr

def quantile_normalization(
        expr,
):
    r"""Quantile normalization of a gene expression matrix (F42-compliant).

    Computes average expression per rank across tissues.

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
        Gene Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    dict
        with keys:

        normalized_expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F'), read-only
            Normalized `expr`
            A result is a value; call `.copy()` to obtain a modifiable array.
        rank_means : np.ndarray[np.float64] of shape (n_genes,), read-only
            The mean of each rank across tissues, one per gene
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
    Generated from the Fortran procedure `tox_normalization::quantile_normalization`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expr = _np.asfortranarray(expr, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expr' must be an array of np.float64: {_error}", argument="expr") from None
    if expr.ndim != 2:
        raise _ToxArgumentError(f"'expr' must have 2 dimensions, but has {expr.ndim}", argument="expr")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_replicates = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    normalized_expr = _np.empty((n_replicates, n_genes,), dtype=_np.float64, order='F')
    rank_means = _np.empty((n_genes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.quantile_normalization_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_replicates)),
        expr,
        normalized_expr,
        rank_means,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _QUANTILE_NORMALIZATION_ARGUMENTS, _QUANTILE_NORMALIZATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    normalized_expr.flags.writeable = False
    rank_means.flags.writeable = False

    return {
        "normalized_expr": normalized_expr,
        "rank_means": rank_means,
    }

def log2_transformation(
        expr,
):
    r"""Apply `log2(x + 1)` transformation to each element of the input matrix.

    This subroutine performs element-wise `log2(x + 1)` transformation on a
    matrix flattened in column-major order. The `log2` is computed via:
    `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
    non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Gene Expression matrix, from :func:`tensor_omics.calc_tiss_avg`
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    transformed_expr : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F'), read-only
        Log-transformed `expr`
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
    Generated from the Fortran procedure `tox_normalization::log2_transformation`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expr = _np.asfortranarray(expr, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expr' must be an array of np.float64: {_error}", argument="expr") from None
    if expr.ndim != 2:
        raise _ToxArgumentError(f"'expr' must have 2 dimensions, but has {expr.ndim}", argument="expr")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_tissues = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    transformed_expr = _np.empty((n_tissues, n_genes,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.log2_transformation_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_tissues)),
        expr,
        transformed_expr,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _LOG2_TRANSFORMATION_ARGUMENTS, _LOG2_TRANSFORMATION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    transformed_expr.flags.writeable = False

    return transformed_expr

def calc_tiss_avg(
        reps_per_tissue,
        expr,
):
    r"""Calculate tissue averages by averaging replicates within each tissue.

    For each tissue of tissue replicates, this subroutine computes the average
    expression per gene.

    Parameters
    ----------
    reps_per_tissue : np.ndarray[np.int32] of shape (n_tissues,)
        Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
        e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
        The minimum valid value is `1`.
    expr : np.ndarray[np.float64] of shape (sum(reps_per_tissue), n_genes,), column-major (order='F')
        Gene Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    tissue_averages : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F'), read-only
        Tissue averages per gene
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
    Generated from the Fortran procedure `tox_normalization::calc_tiss_avg`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        reps_per_tissue = _np.ascontiguousarray(reps_per_tissue, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'reps_per_tissue' must be an array of np.int32: {_error}", argument="reps_per_tissue") from None
    if reps_per_tissue.ndim != 1:
        raise _ToxArgumentError(f"'reps_per_tissue' must have 1 dimension, but has {reps_per_tissue.ndim}", argument="reps_per_tissue")
    try:
        expr = _np.asfortranarray(expr, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expr' must be an array of np.float64: {_error}", argument="expr") from None
    if expr.ndim != 2:
        raise _ToxArgumentError(f"'expr' must have 2 dimensions, but has {expr.ndim}", argument="expr")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_tissues = reps_per_tissue.shape[0]
    n_expr_elements_dim_1 = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    tissue_averages = _np.empty((n_tissues, n_genes,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.calc_tiss_avg_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_tissues)),
        reps_per_tissue,
        expr,
        _ctypes.byref(_ctypes.c_int(n_expr_elements_dim_1)),
        tissue_averages,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_TISS_AVG_ARGUMENTS, _CALC_TISS_AVG_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    tissue_averages.flags.writeable = False

    return tissue_averages

def calc_fchange(
        control_tissues,
        condition_tissues,
        expr,
):
    r"""Calculate `log2 fold changes` between condition and control groups.

    For each control-condition pair, this subroutine computes the `log2 fold change`
    by subtracting the expression value in the control group from the corresponding
    value in the condition group, for all genes.

    Parameters
    ----------
    control_tissues : np.ndarray[np.int32] of shape (n_pairs,)
        Control tissue indices
        The minimum valid value is `1`.
        The maximum valid value is `n_tissues`.
    condition_tissues : np.ndarray[np.int32] of shape (n_pairs,)
        Condition tissue indices
        The minimum valid value is `1`.
        The maximum valid value is `n_tissues`.
    expr : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Gene Expression matrix, from :func:`tensor_omics.calc_tiss_avg`
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    fold_changes : np.ndarray[np.float64] of shape (n_pairs, n_genes,), column-major (order='F'), read-only
        Output matrix for fold changes
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
    Generated from the Fortran procedure `tox_normalization::calc_fchange`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        control_tissues = _np.ascontiguousarray(control_tissues, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'control_tissues' must be an array of np.int32: {_error}", argument="control_tissues") from None
    if control_tissues.ndim != 1:
        raise _ToxArgumentError(f"'control_tissues' must have 1 dimension, but has {control_tissues.ndim}", argument="control_tissues")
    try:
        condition_tissues = _np.ascontiguousarray(condition_tissues, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'condition_tissues' must be an array of np.int32: {_error}", argument="condition_tissues") from None
    if condition_tissues.ndim != 1:
        raise _ToxArgumentError(f"'condition_tissues' must have 1 dimension, but has {condition_tissues.ndim}", argument="condition_tissues")
    try:
        expr = _np.asfortranarray(expr, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expr' must be an array of np.float64: {_error}", argument="expr") from None
    if expr.ndim != 2:
        raise _ToxArgumentError(f"'expr' must have 2 dimensions, but has {expr.ndim}", argument="expr")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_tissues = expr.shape[0]
    n_pairs = control_tissues.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if condition_tissues.shape[0] != n_pairs:
        raise _ToxArgumentError(f"'condition_tissues' has {condition_tissues.shape[0]} along axis 0, but 'control_tissues' implies n_pairs == {n_pairs}; 'condition_tissues' has shape (n_pairs)", argument="condition_tissues")

    # outputs and work arrays, which the caller never sees
    fold_changes = _np.empty((n_pairs, n_genes,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.calc_fchange_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(_ctypes.c_int(n_tissues)),
        _ctypes.byref(_ctypes.c_int(n_pairs)),
        control_tissues,
        condition_tissues,
        expr,
        fold_changes,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_FCHANGE_ARGUMENTS, _CALC_FCHANGE_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    fold_changes.flags.writeable = False

    return fold_changes
