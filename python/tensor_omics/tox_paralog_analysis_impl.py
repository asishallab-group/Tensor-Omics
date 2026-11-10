r"""tox_paralog_analysis_impl

Paralog-subset expression patterns -- dosage effect and subfunctionalization -- relative to an ancestral ortholog.

Candidate paralog subsets are enumerated as bitmask-encoded gene sets, built up one gene at a time
starting from single genes. At every extension step the candidate is scored against the pattern-specific
criterion (small angle plus magnitude gain for dosage effect, or bounded residual distance for
subfunctionalization); subsets that can no longer satisfy the criterion are pruned instead of being
extended further, which keeps the combinatorial subset search tractable.

There is a routine per pattern rather than one taking a pattern flag:
`detect_dosage_effect` and `detect_subfunctionalization` find the subsets, and
`filter_paralogs_by_pattern_dosage_effect` / `_subfunctionalization` reduce a set of genes to
those matching. Which criterion a result was scored against is therefore visible at the call
site, and each routine takes only the thresholds its own pattern needs.

Python binding, generated from tox_paralog_analysis_impl. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.mask_check_state_c.restype = None
_lib.mask_check_state_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MASK_CHECK_STATE_ARGUMENTS = ("bit_mask", "i_gene",)

_lib.mask_chunk_count_c.restype = None
_lib.mask_chunk_count_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MASK_CHUNK_COUNT_ARGUMENTS = ("n_genes", "count",)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis_impl::mask_check_state`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    bit_mask = _contract.integers(bit_mask, "bit_mask", ndim=1)
    i_gene = _contract.integer(i_gene, "i_gene")

    # what the inputs already say, rather than asking for it again
    n_bit_mask_elements = bit_mask.shape[0]

    # outputs and work arrays, which the caller never sees
    state = _ctypes.c_bool(0)
    ierr = _ctypes.c_int(0)

    _lib.mask_check_state_c(
        bit_mask,
        _ctypes.byref(_ctypes.c_int(n_bit_mask_elements)),
        _ctypes.byref(_ctypes.c_int(i_gene)),
        _ctypes.byref(state),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _MASK_CHECK_STATE_ARGUMENTS)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_paralog_analysis_impl::mask_chunk_count`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_genes = _contract.integer(n_genes, "n_genes")

    # outputs and work arrays, which the caller never sees
    count = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.mask_chunk_count_c(
        _ctypes.byref(_ctypes.c_int(n_genes)),
        _ctypes.byref(count),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _MASK_CHUNK_COUNT_ARGUMENTS)

    return count.value
