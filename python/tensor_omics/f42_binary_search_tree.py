"""Python binding to Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.

Generated from f42_binary_search_tree. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.build_bst_index_c.restype = None
_lib.build_bst_index_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_BST_INDEX_ARGUMENTS = ("values", "n_values", "sorted_indices", "ierr",)

_lib.bst_range_query_c.restype = None
_lib.bst_range_query_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BST_RANGE_QUERY_ARGUMENTS = ("values", "sorted_indices", "n_values", "lower_bound", "upper_bound", "output_indices", "n_matches", "ierr",)

def build_bst_index(
        values,
):
    r"""Build the BST index by sorting indices using values in x

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Input real array to be indexed

    Returns
    -------
    sorted_indices : np.ndarray[np.int32] of shape (n_values,)
        Output permutation index

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_binary_search_tree::build_bst_index`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = np.ascontiguousarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values' must be an array of np.float64: {error}") from None
    if values.ndim != 1:
        raise ValueError(f"'values' must have 1 dimension, but has {values.ndim}")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # outputs and work arrays, which the caller never sees
    sorted_indices = np.empty((n_values,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.build_bst_index_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        sorted_indices,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BUILD_BST_INDEX_ARGUMENTS)

    return sorted_indices

def bst_range_query(
        values,
        sorted_indices,
        lower_bound,
        upper_bound,
):
    r"""Perform a 1D range query over the sorted index

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Input real array
    sorted_indices : np.ndarray[np.int32] of shape (n_values,)
        Permutation index array (sorted)
    lower_bound : float
        Lower bound of range (inclusive)
    upper_bound : float
        Upper bound of range (inclusive)

    Returns
    -------
    output_indices : np.ndarray[np.int32] of shape (n_values,)
        Output array of matching indices.
        The first `n_matches` elements will hold the results.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_binary_search_tree::bst_range_query`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = np.ascontiguousarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values' must be an array of np.float64: {error}") from None
    if values.ndim != 1:
        raise ValueError(f"'values' must have 1 dimension, but has {values.ndim}")
    try:
        sorted_indices = np.ascontiguousarray(sorted_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'sorted_indices' must be an array of np.int32: {error}") from None
    if sorted_indices.ndim != 1:
        raise ValueError(f"'sorted_indices' must have 1 dimension, but has {sorted_indices.ndim}")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if sorted_indices.shape[0] != n_values:
        raise ValueError(f"'sorted_indices' has {sorted_indices.shape[0]} along axis 0, but "
            f"'values' implies n_values == {n_values}"
        )

    # outputs and work arrays, which the caller never sees
    output_indices = np.empty((n_values,), dtype=np.int32, order='C')
    n_matches = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.bst_range_query_c(
        values,
        sorted_indices,
        ctypes.byref(ctypes.c_int(n_values)),
        ctypes.byref(ctypes.c_double(lower_bound)),
        ctypes.byref(ctypes.c_double(upper_bound)),
        output_indices,
        ctypes.byref(n_matches),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BST_RANGE_QUERY_ARGUMENTS)

    return output_indices[..., :n_matches.value]
