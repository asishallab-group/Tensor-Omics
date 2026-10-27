r"""f42_binary_search_tree

Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
This module provides routines to build a BST index (via sorting), access sorted values,
and perform range queries over a real-valued array using the sorted index.

Python binding, generated from f42_binary_search_tree. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.build_bst_index_c.restype = None
_lib.build_bst_index_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_BST_INDEX_ARGUMENTS = ("values", "n_values", "sorted_indices", "ierr",)
#: For a derived argument, the one the caller passed it in
_BUILD_BST_INDEX_ARGUMENT_SOURCES = (None, "values", None, None,)

_lib.bst_range_query_c.restype = None
_lib.bst_range_query_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BST_RANGE_QUERY_ARGUMENTS = ("values", "n_values", "lower_bound", "upper_bound", "output_indices", "n_matches", "ierr",)
#: For a derived argument, the one the caller passed it in
_BST_RANGE_QUERY_ARGUMENT_SOURCES = (None, "values", None, None, None, None, None,)

_lib.bst_range_query_expert_c.restype = None
_lib.bst_range_query_expert_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_binary_search_tree::build_bst_index`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        values = _np.ascontiguousarray(values, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'values' must be an array of np.float64: {_error}", argument="values") from None
    if values.ndim != 1:
        raise _ToxArgumentError(f"'values' must have 1 dimension, but has {values.ndim}", argument="values")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # outputs and work arrays, which the caller never sees
    sorted_indices = _np.empty((n_values,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.build_bst_index_c(
        values,
        _ctypes.byref(_ctypes.c_int(n_values)),
        sorted_indices,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BUILD_BST_INDEX_ARGUMENTS, _BUILD_BST_INDEX_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
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
        values = _np.ascontiguousarray(values, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'values' must be an array of np.float64: {_error}", argument="values") from None
    if values.ndim != 1:
        raise _ToxArgumentError(f"'values' must have 1 dimension, but has {values.ndim}", argument="values")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # outputs and work arrays, which the caller never sees
    output_indices = _np.empty((n_values,), dtype=_np.int32, order='C')
    n_matches = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.bst_range_query_c(
        values,
        _ctypes.byref(_ctypes.c_int(n_values)),
        _ctypes.byref(_ctypes.c_double(lower_bound)),
        _ctypes.byref(_ctypes.c_double(upper_bound)),
        output_indices,
        _ctypes.byref(n_matches),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BST_RANGE_QUERY_ARGUMENTS, _BST_RANGE_QUERY_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
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
        values = _np.ascontiguousarray(values, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'values' must be an array of np.float64: {_error}", argument="values") from None
    if values.ndim != 1:
        raise _ToxArgumentError(f"'values' must have 1 dimension, but has {values.ndim}", argument="values")
    try:
        values_perm = _np.ascontiguousarray(values_perm, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'values_perm' must be an array of np.int32: {_error}", argument="values_perm") from None
    if values_perm.ndim != 1:
        raise _ToxArgumentError(f"'values_perm' must have 1 dimension, but has {values_perm.ndim}", argument="values_perm")

    # what the inputs already say, rather than asking for it again
    n_values = values.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if values_perm.shape[0] != n_values:
        raise _ToxArgumentError(f"'values_perm' has {values_perm.shape[0]} along axis 0, but 'values' implies n_values == {n_values}; 'values_perm' has shape (n_values)", argument="values_perm")

    # outputs and work arrays, which the caller never sees
    output_indices = _np.empty((n_values,), dtype=_np.int32, order='C')
    n_matches = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.bst_range_query_expert_c(
        values,
        values_perm,
        _ctypes.byref(_ctypes.c_int(n_values)),
        _ctypes.byref(_ctypes.c_double(lower_bound)),
        _ctypes.byref(_ctypes.c_double(upper_bound)),
        output_indices,
        _ctypes.byref(n_matches),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BST_RANGE_QUERY_EXPERT_ARGUMENTS, _BST_RANGE_QUERY_EXPERT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    output_indices.flags.writeable = False

    return output_indices[..., :n_matches.value]
