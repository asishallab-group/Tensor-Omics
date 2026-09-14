r"""f42_binary_search_tree

Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
This module provides routines to build a BST index (via sorting), access sorted values,
and perform range queries over a real-valued array using the sorted index.

Python binding, generated from f42_binary_search_tree. Do not edit.
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
#: For a derived argument, the one the caller passed it in
_BUILD_BST_INDEX_ARGUMENT_SOURCES = (None, "values", None, None,)

_lib.bst_range_query_c.restype = None
_lib.bst_range_query_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BST_RANGE_QUERY_ARGUMENTS = ("values", "n_values", "lower_bound", "upper_bound", "output_indices", "n_matches", "ierr",)
#: For a derived argument, the one the caller passed it in
_BST_RANGE_QUERY_ARGUMENT_SOURCES = (None, "values", None, None, None, None, None,)

_lib.bst_range_query_expert_c.restype = None
_lib.bst_range_query_expert_c.argtypes = (
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
_BST_RANGE_QUERY_EXPERT_ARGUMENTS = ("values", "values_perm", "n_values", "lower_bound", "upper_bound", "output_indices", "n_matches", "ierr",)
#: For a derived argument, the one the caller passed it in
_BST_RANGE_QUERY_EXPERT_ARGUMENT_SOURCES = (None, None, "values", None, None, None, None, None,)

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
    sorted_indices : np.ndarray[np.int32] of shape (n_values,), read-only
        Output permutation index
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_binary_search_tree::build_bst_index`, whose argument names are
    the ones an error message reports.
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

    check_err_code(ierr.value, _BUILD_BST_INDEX_ARGUMENTS, _BUILD_BST_INDEX_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    sorted_indices.flags.writeable = False

    return sorted_indices

def bst_range_query(
        values,
        lower_bound,
        upper_bound,
):
    r"""Perform a 1D range query over the sorted index

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Input real array
    lower_bound : float
        Lower bound of range (inclusive)
        The maximum valid value is `upper_bound`.
    upper_bound : float
        Upper bound of range (inclusive)

    Returns
    -------
    output_indices : np.ndarray[np.int32] of shape (n_values,), read-only
        Output array of matching indices.
        The first `n_matches` elements will hold the results.
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_binary_search_tree::bst_range_query`, whose argument names are
    the ones an error message reports.

    This entry point seeds `values_perm` and sorts it by `values`.
    Call `bst_range_query_expert` to do that yourself.
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
    output_indices = np.empty((n_values,), dtype=np.int32, order='C')
    n_matches = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.bst_range_query_c(
        values,
        ctypes.byref(ctypes.c_int(n_values)),
        ctypes.byref(ctypes.c_double(lower_bound)),
        ctypes.byref(ctypes.c_double(upper_bound)),
        output_indices,
        ctypes.byref(n_matches),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BST_RANGE_QUERY_ARGUMENTS, _BST_RANGE_QUERY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    output_indices.flags.writeable = False

    return output_indices[..., :n_matches.value]

def bst_range_query_expert(
        values,
        values_perm,
        lower_bound,
        upper_bound,
):
    r"""Perform a 1D range query over the sorted index

    Parameters
    ----------
    values : np.ndarray[np.float64] of shape (n_values,)
        Input real array
    values_perm : np.ndarray[np.int32] of shape (n_values,)
        Permutation of `values` in ascending order -- the BST index. The allocating entry
        point builds and heapsorts it for you; the expert one takes whatever order you
        supply, so a caller that already holds one from
        :func:`tensor_omics.build_bst_index` can reuse it across queries.
        The minimum valid value is `1`.
        The maximum valid value is `n_values`.
    lower_bound : float
        Lower bound of range (inclusive)
        The maximum valid value is `upper_bound`.
    upper_bound : float
        Upper bound of range (inclusive)

    Returns
    -------
    output_indices : np.ndarray[np.int32] of shape (n_values,), read-only
        Output array of matching indices.
        The first `n_matches` elements will hold the results.
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_binary_search_tree::bst_range_query_expert`, whose argument names are
    the ones an error message reports.

    The expert entry point: you supply `values_perm` yourself.
    `bst_range_query` seeds `values_perm` and sorts it by `values`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = np.ascontiguousarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values' must be an array of np.float64: {error}") from None
    if values.ndim != 1:
        raise ValueError(f"'values' must have 1 dimension, but has {values.ndim}")
    try:
        values_perm = np.ascontiguousarray(values_perm, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'values_perm' must be an array of np.int32: {error}") from None
    if values_perm.ndim != 1:
        raise ValueError(f"'values_perm' must have 1 dimension, but has {values_perm.ndim}")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if values_perm.shape[0] != n_values:
        raise ValueError(f"'values_perm' has {values_perm.shape[0]} along axis 0, but "
            f"'values' implies n_values == {n_values}"
        )

    # outputs and work arrays, which the caller never sees
    output_indices = np.empty((n_values,), dtype=np.int32, order='C')
    n_matches = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.bst_range_query_expert_c(
        values,
        values_perm,
        ctypes.byref(ctypes.c_int(n_values)),
        ctypes.byref(ctypes.c_double(lower_bound)),
        ctypes.byref(ctypes.c_double(upper_bound)),
        output_indices,
        ctypes.byref(n_matches),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BST_RANGE_QUERY_EXPERT_ARGUMENTS, _BST_RANGE_QUERY_EXPERT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    output_indices.flags.writeable = False

    return output_indices[..., :n_matches.value]
