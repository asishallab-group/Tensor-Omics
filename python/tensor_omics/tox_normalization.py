"""tox_normalization

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_normalization. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.normalize_unit_length_c.restype = None
_lib.normalize_unit_length_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_UNIT_LENGTH_ARGUMENTS = ("vector", "n_dims", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZE_UNIT_LENGTH_ARGUMENT_SOURCES = (None, "vector", None,)

_lib.normalization_pipeline_c.restype = None
_lib.normalization_pipeline_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZATION_PIPELINE_ARGUMENTS = ("n_genes", "n_replicates", "expr", "log_transformed_expr", "reps_per_tissue", "n_tissues", "span", "degree", "use_quantile", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZATION_PIPELINE_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, "log_transformed_expr", None, None, None, None,)

_lib.normalize_by_std_dev_c.restype = None
_lib.normalize_by_std_dev_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_BY_STD_DEV_ARGUMENTS = ("n_genes", "n_replicates", "expr", "normalized_expr", "span", "degree", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZE_BY_STD_DEV_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, None, None,)

_lib.root_mean_sq_normalization_c.restype = None
_lib.root_mean_sq_normalization_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ROOT_MEAN_SQ_NORMALIZATION_ARGUMENTS = ("n_genes", "n_replicates", "expr", "normalized_expr", "ierr",)
#: For a derived argument, the one the caller passed it in
_ROOT_MEAN_SQ_NORMALIZATION_ARGUMENT_SOURCES = ("expr", "expr", None, None, None,)

_lib.quantile_normalization_expert_c.restype = None
_lib.quantile_normalization_expert_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_QUANTILE_NORMALIZATION_EXPERT_ARGUMENTS = ("n_genes", "n_replicates", "expr", "normalized_expr", "rank_means", "tmp_genes_row", "tmp_perm", "ierr",)
#: For a derived argument, the one the caller passed it in
_QUANTILE_NORMALIZATION_EXPERT_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, None, None, None,)

_lib.quantile_normalization_c.restype = None
_lib.quantile_normalization_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_QUANTILE_NORMALIZATION_ARGUMENTS = ("n_genes", "n_replicates", "expr", "normalized_expr", "rank_means", "ierr",)
#: For a derived argument, the one the caller passed it in
_QUANTILE_NORMALIZATION_ARGUMENT_SOURCES = ("expr", "expr", None, None, None, None,)

_lib.log2_transformation_c.restype = None
_lib.log2_transformation_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOG2_TRANSFORMATION_ARGUMENTS = ("n_genes", "n_tissues", "expr", "transformed_expr", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOG2_TRANSFORMATION_ARGUMENT_SOURCES = ("expr", "expr", None, None, None,)

_lib.calc_tiss_avg_c.restype = None
_lib.calc_tiss_avg_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_TISS_AVG_ARGUMENTS = ("n_genes", "n_tissues", "reps_per_tissue", "expr", "tissue_averages", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_TISS_AVG_ARGUMENT_SOURCES = ("expr", "reps_per_tissue", None, None, None, None,)

_lib.calc_fchange_c.restype = None
_lib.calc_fchange_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::normalize_unit_length`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    if not isinstance(vector, np.ndarray) or vector.dtype != np.float64:
        raise TypeError("'vector' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if vector.ndim != 1:
        raise ValueError(f"'vector' must have 1 dimension, but has {vector.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dims = vector.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.normalize_unit_length_c(
        vector,
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMALIZE_UNIT_LENGTH_ARGUMENTS, _NORMALIZE_UNIT_LENGTH_ARGUMENT_SOURCES)

    return None

def normalization_pipeline(
        expr,
        reps_per_tissue,
        span=0.7,
        degree=2,
        use_quantile=False,
):
    r"""Complete normalization pipeline for gene expression data.

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
    log_transformed_expr : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Log-transformed grouped `expr`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::normalization_pipeline_alloc`, whose argument names are
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
        reps_per_tissue = np.ascontiguousarray(reps_per_tissue, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'reps_per_tissue' must be an array of np.int32: {error}") from None
    if reps_per_tissue.ndim != 1:
        raise ValueError(f"'reps_per_tissue' must have 1 dimension, but has {reps_per_tissue.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_replicates = expr.shape[0]
    n_tissues = reps_per_tissue.shape[0]

    # outputs and work arrays, which the caller never sees
    log_transformed_expr = np.empty((n_tissues, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.normalization_pipeline_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_replicates)),
        expr,
        log_transformed_expr,
        reps_per_tissue,
        ctypes.byref(ctypes.c_int(n_tissues)),
        ctypes.byref(ctypes.c_double(span)),
        ctypes.byref(ctypes.c_int(degree)),
        ctypes.byref(ctypes.c_bool(use_quantile)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMALIZATION_PIPELINE_ARGUMENTS, _NORMALIZATION_PIPELINE_ARGUMENT_SOURCES)

    return log_transformed_expr

def normalize_by_std_dev(
        expr,
        span=0.7,
        degree=2,
):
    r"""Normalizes each gene's expression vector using LOESS-stabilized standard deviation.

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
    normalized_expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
        Normalized `expr`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::normalize_by_std_dev_alloc`, whose argument names are
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
    n_replicates = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    normalized_expr = np.empty((n_replicates, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.normalize_by_std_dev_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_replicates)),
        expr,
        normalized_expr,
        ctypes.byref(ctypes.c_double(span)),
        ctypes.byref(ctypes.c_int(degree)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMALIZE_BY_STD_DEV_ARGUMENTS, _NORMALIZE_BY_STD_DEV_ARGUMENT_SOURCES)

    return normalized_expr

def root_mean_sq_normalization(
        expr,
):
    r"""Normalizes each gene's expression vector using `sqrt(mean(x^2))`

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
        Gene Expression matrix
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    normalized_expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
        Normalized `expr`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::root_mean_sq_normalization`, whose argument names are
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
    n_replicates = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    normalized_expr = np.empty((n_replicates, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.root_mean_sq_normalization_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_replicates)),
        expr,
        normalized_expr,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ROOT_MEAN_SQ_NORMALIZATION_ARGUMENTS, _ROOT_MEAN_SQ_NORMALIZATION_ARGUMENT_SOURCES)

    return normalized_expr

def quantile_normalization_expert(
        expr,
):
    r"""Quantile normalization of a gene expression matrix (F42-compliant).

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

        normalized_expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
            Normalized `expr`
        rank_means : np.ndarray[np.float64] of shape (n_genes,)
            The mean of each rank across tissues, one per gene

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::quantile_normalization`, whose argument names are
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
    n_replicates = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    normalized_expr = np.empty((n_replicates, n_genes,), dtype=np.float64, order='F')
    rank_means = np.empty((n_genes,), dtype=np.float64, order='C')
    tmp_genes_row = np.empty((n_genes,), dtype=np.float64, order='C')
    tmp_perm = np.empty((n_genes,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.quantile_normalization_expert_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_replicates)),
        expr,
        normalized_expr,
        rank_means,
        tmp_genes_row,
        tmp_perm,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _QUANTILE_NORMALIZATION_EXPERT_ARGUMENTS, _QUANTILE_NORMALIZATION_EXPERT_ARGUMENT_SOURCES)

    return {
        "normalized_expr": normalized_expr,
        "rank_means": rank_means,
    }

def quantile_normalization(
        expr,
):
    r"""Quantile normalization of a gene expression matrix (F42-compliant).

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

        normalized_expr : np.ndarray[np.float64] of shape (n_replicates, n_genes,), column-major (order='F')
            Normalized `expr`
        rank_means : np.ndarray[np.float64] of shape (n_genes,)
            The mean of each rank across tissues, one per gene

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::quantile_normalization_alloc`, whose argument names are
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
    n_replicates = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    normalized_expr = np.empty((n_replicates, n_genes,), dtype=np.float64, order='F')
    rank_means = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.quantile_normalization_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_replicates)),
        expr,
        normalized_expr,
        rank_means,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _QUANTILE_NORMALIZATION_ARGUMENTS, _QUANTILE_NORMALIZATION_ARGUMENT_SOURCES)

    return {
        "normalized_expr": normalized_expr,
        "rank_means": rank_means,
    }

def log2_transformation(
        expr,
):
    r"""Apply `log2(x + 1)` transformation to each element of the input matrix.

    Parameters
    ----------
    expr : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Gene Expression matrix, from :func:`tensor_omics.calc_tiss_avg`
        NaN is permitted for this value.
        Infinite values are permitted for this value.

    Returns
    -------
    transformed_expr : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Log-transformed `expr`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::log2_transformation`, whose argument names are
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
    n_tissues = expr.shape[0]

    # outputs and work arrays, which the caller never sees
    transformed_expr = np.empty((n_tissues, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.log2_transformation_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_tissues)),
        expr,
        transformed_expr,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOG2_TRANSFORMATION_ARGUMENTS, _LOG2_TRANSFORMATION_ARGUMENT_SOURCES)

    return transformed_expr

def calc_tiss_avg(
        reps_per_tissue,
        expr,
):
    r"""Calculate tissue averages by averaging replicates within each tissue.

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
    tissue_averages : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Tissue averages per gene

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::calc_tiss_avg`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        reps_per_tissue = np.ascontiguousarray(reps_per_tissue, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'reps_per_tissue' must be an array of np.int32: {error}") from None
    if reps_per_tissue.ndim != 1:
        raise ValueError(f"'reps_per_tissue' must have 1 dimension, but has {reps_per_tissue.ndim}")
    try:
        expr = np.asfortranarray(expr, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expr' must be an array of np.float64: {error}") from None
    if expr.ndim != 2:
        raise ValueError(f"'expr' must have 2 dimensions, but has {expr.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_tissues = reps_per_tissue.shape[0]

    # outputs and work arrays, which the caller never sees
    tissue_averages = np.empty((n_tissues, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.calc_tiss_avg_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_tissues)),
        reps_per_tissue,
        expr,
        tissue_averages,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_TISS_AVG_ARGUMENTS, _CALC_TISS_AVG_ARGUMENT_SOURCES)

    return tissue_averages

def calc_fchange(
        control_tissues,
        condition_tissues,
        expr,
):
    r"""Calculate `log2 fold changes` between condition and control groups.

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
    fold_changes : np.ndarray[np.float64] of shape (n_pairs, n_genes,), column-major (order='F')
        Output matrix for fold changes

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_normalization::calc_fchange`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        control_tissues = np.ascontiguousarray(control_tissues, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'control_tissues' must be an array of np.int32: {error}") from None
    if control_tissues.ndim != 1:
        raise ValueError(f"'control_tissues' must have 1 dimension, but has {control_tissues.ndim}")
    try:
        condition_tissues = np.ascontiguousarray(condition_tissues, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'condition_tissues' must be an array of np.int32: {error}") from None
    if condition_tissues.ndim != 1:
        raise ValueError(f"'condition_tissues' must have 1 dimension, but has {condition_tissues.ndim}")
    try:
        expr = np.asfortranarray(expr, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expr' must be an array of np.float64: {error}") from None
    if expr.ndim != 2:
        raise ValueError(f"'expr' must have 2 dimensions, but has {expr.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = expr.shape[1]
    n_tissues = expr.shape[0]
    n_pairs = control_tissues.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if condition_tissues.shape[0] != n_pairs:
        raise ValueError(f"'condition_tissues' has {condition_tissues.shape[0]} along axis 0, but "
            f"'control_tissues' implies n_pairs == {n_pairs}"
        )

    # outputs and work arrays, which the caller never sees
    fold_changes = np.empty((n_pairs, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.calc_fchange_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_tissues)),
        ctypes.byref(ctypes.c_int(n_pairs)),
        control_tissues,
        condition_tissues,
        expr,
        fold_changes,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_FCHANGE_ARGUMENTS, _CALC_FCHANGE_ARGUMENT_SOURCES)

    return fold_changes
