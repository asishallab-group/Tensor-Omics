"""tox_paralog_analysis_impl

Implementations for detecting paralog-subset expression patterns (dosage effect and subfunctionalization) relative to an ancestral ortholog.

Python binding, generated from tox_paralog_analysis_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.mask_check_state_c.restype = None
_lib.mask_check_state_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MASK_CHECK_STATE_ARGUMENTS = ("bit_mask", "i_gene",)

_lib.mask_chunk_count_c.restype = None
_lib.mask_chunk_count_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MASK_CHUNK_COUNT_ARGUMENTS = ("n_genes", "count",)

_lib.calc_work_arr_paralog_subsets_size_c.restype = None
_lib.calc_work_arr_paralog_subsets_size_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_WORK_ARR_PARALOG_SUBSETS_SIZE_ARGUMENTS = ("max_subset_size", "n_genes", "work_array_size", "filtered_paralogs_mask", "n_mask_chunks", "ierr",)
#: For a derived argument, the one the caller passed it in
_CALC_WORK_ARR_PARALOG_SUBSETS_SIZE_ARGUMENT_SOURCES = (None, None, None, None, "filtered_paralogs_mask", None,)

def mask_check_state(
        bit_mask,
        i_gene,
):
    r"""Checks the state of a bit/paralog in `bit_mask` -> True if 1 else False

    Parameters
    ----------
    bit_mask : np.ndarray[np.int32] of shape (n_bit_mask_elements,)
        chunked mask to mark active paralogs
    i_gene : int
        index of paralog to be marked active

    Returns
    -------
    state : bool
        check result

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis_impl::mask_check_state`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'bit_mask' must be an array of np.int32: {error}") from None
    if bit_mask.ndim != 1:
        raise ValueError(f"'bit_mask' must have 1 dimension, but has {bit_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bit_mask_elements = bit_mask.shape[0]

    # outputs and work arrays, which the caller never sees
    state = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.mask_check_state_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(ctypes.c_int(i_gene)),
        ctypes.byref(state),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _MASK_CHECK_STATE_ARGUMENTS)

    return state.value

def mask_chunk_count(
        n_genes,
):
    r"""Determines the needed chunk count for subset bit masks (an integer has only 32 bits)

    Parameters
    ----------
    n_genes : int
        number of genes

    Returns
    -------
    count : int
        number of 32 bit chunks a mask needs to encode `n_genes` genes

        Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis_impl::mask_chunk_count`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    count = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.mask_chunk_count_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(count),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _MASK_CHUNK_COUNT_ARGUMENTS)

    return count.value

def calc_work_arr_paralog_subsets_size(
        max_subset_size,
        n_genes,
        filtered_paralogs_mask,
):
    r"""Calculates the needed size for the paralog-subsets work array

    Parameters
    ----------
    max_subset_size : int
        maximum size that a subset must not exceed. Zero is in range and means there is
        nothing to size a work array for, which is reported back as a size of zero.
        The minimum valid value is `0`.
        If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.

        Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
    n_genes : int
        number of genes
    filtered_paralogs_mask : np.ndarray[np.int32] of shape (n_mask_chunks,)
        Output mask with all genes disabled that did not pass the filter

    Returns
    -------
    dict
        with keys:

        max_subset_size : int
            maximum size that a subset must not exceed. Zero is in range and means there is
            nothing to size a work array for, which is reported back as a size of zero.
            The minimum valid value is `0`.
            If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.

            Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
        work_array_size : int
            The calculated needed work array size in absolute worst case scenario. Look into source for details.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis_impl::calc_work_arr_paralog_subsets_size`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    max_subset_size = ctypes.c_int(max_subset_size)
    try:
        filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'filtered_paralogs_mask' must be an array of np.int32: {error}") from None
    if filtered_paralogs_mask.ndim != 1:
        raise ValueError(f"'filtered_paralogs_mask' must have 1 dimension, but has {filtered_paralogs_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # outputs and work arrays, which the caller never sees
    work_array_size = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_work_arr_paralog_subsets_size_c(
        ctypes.byref(max_subset_size),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(work_array_size),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_WORK_ARR_PARALOG_SUBSETS_SIZE_ARGUMENTS, _CALC_WORK_ARR_PARALOG_SUBSETS_SIZE_ARGUMENT_SOURCES)

    return {
        "max_subset_size": max_subset_size.value,
        "work_array_size": work_array_size.value,
    }
