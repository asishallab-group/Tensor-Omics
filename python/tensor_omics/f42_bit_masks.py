r"""f42_bit_masks

A bit mask holds the same information as a `logical(c_bool)` array in an eighth of the memory,
and whole-mask operations work on 32 flags at once. `misc/bench/bit_masks/` measures both, and
records why this module is a set of procedures on plain words rather than a derived type.

**Layout.** Flag `i` (counting from 1) is bit `mod(i - 1, 32)` of word `(i - 1)/32 + 1`, the
least significant bit first, so `n_bits` flags take `ceiling(n_bits/32)` words. This is
the order numpy's `np.unpackbits(words.view(np.uint8), bitorder="little")` and R's
`intToBits(words)` decode.

**All 32 bits are used.** A word whose bit 31 is set is negative, so a word is tested against
zero with `/= 0`, never with `> 0`.

**The bits past `n_bits` in the last word are always zero.** Every procedure here keeps it that
way, and the counts, ``bitmask_all`` and
``bitmask_equal`` rely on it.

**Several masks** are the columns of a `(n_words, n_masks)` array. Every column starts on a
new word, so `bitmasks(:, i_mask)` is a mask of its own and passes to any procedure here
without a copy.

**Parallel loops write whole words.** Two iterations that set flags in the same word race, so
a `do concurrent` over flags may only read a mask; one that writes runs over words or columns.

**Hot loops** use the whole-mask procedures, or step from one set flag to the next with
``bitmask_next``. A procedure call per flag costs more
than the flag, because nothing is inlined across modules; where a loop must touch single
flags, the bit-mask macros in `src/macros.h` test a flag or locate its word inline, and one
of them sizes a mask.

Python binding, generated from f42_bit_masks. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.bitmask_n_words_c.restype = None
_lib.bitmask_n_words_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BITMASK_N_WORDS_ARGUMENTS = ("n_bits", "n_words", "ierr",)

_lib.bitmask_test_c.restype = None
_lib.bitmask_test_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BITMASK_TEST_ARGUMENTS = ("n_bits", "n_words", "bitmask", "i_bit", "is_set", "ierr",)
#: For a derived argument, the one the caller passed it in
_BITMASK_TEST_ARGUMENT_SOURCES = (None, "bitmask", None, None, None, None,)

_lib.bitmask_from_logical_c.restype = None
_lib.bitmask_from_logical_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BITMASK_FROM_LOGICAL_ARGUMENTS = ("n_bits", "n_words", "flags", "bitmask", "ierr",)
#: For a derived argument, the one the caller passed it in
_BITMASK_FROM_LOGICAL_ARGUMENT_SOURCES = ("flags", None, None, None, None,)

_lib.bitmask_to_logical_c.restype = None
_lib.bitmask_to_logical_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BITMASK_TO_LOGICAL_ARGUMENTS = ("n_bits", "n_words", "bitmask", "flags", "ierr",)
#: For a derived argument, the one the caller passed it in
_BITMASK_TO_LOGICAL_ARGUMENT_SOURCES = (None, "bitmask", None, None, None,)

_lib.bitmasks_from_logical_2D_c.restype = None
_lib.bitmasks_from_logical_2D_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BITMASKS_FROM_LOGICAL_2D_ARGUMENTS = ("n_bits", "n_masks", "n_words", "flags", "bitmasks", "ierr",)
#: For a derived argument, the one the caller passed it in
_BITMASKS_FROM_LOGICAL_2D_ARGUMENT_SOURCES = ("flags", "flags", None, None, None, None,)

_lib.bitmasks_to_logical_2D_c.restype = None
_lib.bitmasks_to_logical_2D_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BITMASKS_TO_LOGICAL_2D_ARGUMENTS = ("n_bits", "n_masks", "n_words", "bitmasks", "flags", "ierr",)
#: For a derived argument, the one the caller passed it in
_BITMASKS_TO_LOGICAL_2D_ARGUMENT_SOURCES = (None, "bitmasks", "bitmasks", None, None, None,)

def bitmask_n_words(
        n_bits,
):
    r"""The number of words a mask of `n_bits` flags takes.

    Inside Fortran, a bit-mask macro in `src/macros.h` gives the same number.

    Parameters
    ----------
    n_bits : int
        number of flags in the mask; no flags take no words
        The minimum valid value is `0`.

    Returns
    -------
    n_words : int
        number of words the mask takes

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bitmask_n_words`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_bits = _contract.integer(n_bits, "n_bits")

    # outputs and work arrays, which the caller never sees
    n_words = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.bitmask_n_words_c(
        _ctypes.byref(_ctypes.c_int(n_bits)),
        _ctypes.byref(n_words),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BITMASK_N_WORDS_ARGUMENTS)

    return n_words.value

def bitmask_test(
        n_bits,
        bitmask,
        i_bit,
):
    r"""Whether flag `i_bit` of a mask is set.

    Inside Fortran, a loop tests a flag inline with a bit-mask macro from `src/macros.h` instead.

    Parameters
    ----------
    n_bits : int
        number of flags in the mask
        The minimum valid value is `1`.
    bitmask : np.ndarray[np.int32] of shape (n_words,)
        the mask
    i_bit : int
        index of the flag, counting from 1
        The minimum valid value is `1`.
        The maximum valid value is `n_bits`.

    Returns
    -------
    is_set : bool
        `True` if the flag is set

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bitmask_test`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_bits = _contract.integer(n_bits, "n_bits")
    bitmask = _contract.integers(bitmask, "bitmask", ndim=1)
    i_bit = _contract.integer(i_bit, "i_bit")

    # what the inputs already say, rather than asking for it again
    n_words = bitmask.shape[0]

    # outputs and work arrays, which the caller never sees
    is_set = _ctypes.c_bool(0)
    ierr = _ctypes.c_int(0)

    _lib.bitmask_test_c(
        _ctypes.byref(_ctypes.c_int(n_bits)),
        _ctypes.byref(_ctypes.c_int(n_words)),
        bitmask,
        _ctypes.byref(_ctypes.c_int(i_bit)),
        _ctypes.byref(is_set),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BITMASK_TEST_ARGUMENTS, _BITMASK_TEST_ARGUMENT_SOURCES)

    return is_set.value

def bitmask_from_logical(
        n_words,
        flags,
):
    r"""Packs a `logical(c_bool)` array into a bit mask.

    Parameters
    ----------
    n_words : int
        number of words of the mask
        The minimum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
        The maximum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
    flags : np.ndarray[np.bool_] of shape (n_bits,)
        the flags to pack

    Returns
    -------
    bitmask : np.ndarray[np.int32] of shape (n_words,), read-only
        the mask, flag `i` set where `flags(i)` is `True`
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
    Generated from the Fortran procedure `f42_bit_masks::bitmask_from_logical`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_words = _contract.integer(n_words, "n_words")
    flags = _contract.logicals(flags, "flags", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_bits = flags.shape[0]

    # outputs and work arrays, which the caller never sees
    bitmask = _np.empty((n_words,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.bitmask_from_logical_c(
        _ctypes.byref(_ctypes.c_int(n_bits)),
        _ctypes.byref(_ctypes.c_int(n_words)),
        flags,
        bitmask,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BITMASK_FROM_LOGICAL_ARGUMENTS, _BITMASK_FROM_LOGICAL_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    bitmask.flags.writeable = False

    return bitmask

def bitmask_to_logical(
        n_bits,
        bitmask,
):
    r"""Unpacks a bit mask into a `logical(c_bool)` array.

    Parameters
    ----------
    n_bits : int
        number of flags
        The minimum valid value is `1`.
    bitmask : np.ndarray[np.int32] of shape (n_words,)
        the mask

    Returns
    -------
    flags : np.ndarray[np.bool_] of shape (n_bits,), read-only
        `flags(i)` is `True` where flag `i` of the mask is set
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
    Generated from the Fortran procedure `f42_bit_masks::bitmask_to_logical`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_bits = _contract.integer(n_bits, "n_bits")
    bitmask = _contract.integers(bitmask, "bitmask", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_words = bitmask.shape[0]

    # outputs and work arrays, which the caller never sees
    flags = _np.empty((n_bits,), dtype=_np.bool_, order='C')
    ierr = _ctypes.c_int(0)

    _lib.bitmask_to_logical_c(
        _ctypes.byref(_ctypes.c_int(n_bits)),
        _ctypes.byref(_ctypes.c_int(n_words)),
        bitmask,
        flags,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BITMASK_TO_LOGICAL_ARGUMENTS, _BITMASK_TO_LOGICAL_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    flags.flags.writeable = False

    return flags

def bitmasks_from_logical_2D(
        n_words,
        flags,
):
    r"""Packs every column of a `logical(c_bool)` matrix into its own bit mask.

    Parameters
    ----------
    n_words : int
        number of words per mask
        The minimum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
        The maximum valid value is `((n_bits)/32 + min(1, mod((n_bits), 32)))`.
    flags : np.ndarray[np.bool_] of shape (n_bits, n_masks,), column-major (order='F')
        the flags, one mask per column

    Returns
    -------
    bitmasks : np.ndarray[np.int32] of shape (n_words, n_masks,), column-major (order='F'), read-only
        the masks, column `j` packed from `flags(:, j)`
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
    Generated from the Fortran procedure `f42_bit_masks::bitmasks_from_logical_2D`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_words = _contract.integer(n_words, "n_words")
    flags = _contract.logicals(flags, "flags", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_bits = flags.shape[0]
    n_masks = flags.shape[1]

    # outputs and work arrays, which the caller never sees
    bitmasks = _np.empty((n_words, n_masks,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.bitmasks_from_logical_2D_c(
        _ctypes.byref(_ctypes.c_int(n_bits)),
        _ctypes.byref(_ctypes.c_int(n_masks)),
        _ctypes.byref(_ctypes.c_int(n_words)),
        flags,
        bitmasks,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BITMASKS_FROM_LOGICAL_2D_ARGUMENTS, _BITMASKS_FROM_LOGICAL_2D_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    bitmasks.flags.writeable = False

    return bitmasks

def bitmasks_to_logical_2D(
        n_bits,
        bitmasks,
):
    r"""Unpacks every column of a matrix of bit masks into a `logical(c_bool)` column.

    Parameters
    ----------
    n_bits : int
        number of flags per mask
        The minimum valid value is `1`.
    bitmasks : np.ndarray[np.int32] of shape (n_words, n_masks,), column-major (order='F')
        the masks, one per column

    Returns
    -------
    flags : np.ndarray[np.bool_] of shape (n_bits, n_masks,), column-major (order='F'), read-only
        the flags, column `j` unpacked from `bitmasks(:, j)`
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
    Generated from the Fortran procedure `f42_bit_masks::bitmasks_to_logical_2D`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_bits = _contract.integer(n_bits, "n_bits")
    bitmasks = _contract.integers(bitmasks, "bitmasks", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_masks = bitmasks.shape[1]
    n_words = bitmasks.shape[0]

    # outputs and work arrays, which the caller never sees
    flags = _np.empty((n_bits, n_masks,), dtype=_np.bool_, order='F')
    ierr = _ctypes.c_int(0)

    _lib.bitmasks_to_logical_2D_c(
        _ctypes.byref(_ctypes.c_int(n_bits)),
        _ctypes.byref(_ctypes.c_int(n_masks)),
        _ctypes.byref(_ctypes.c_int(n_words)),
        bitmasks,
        flags,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BITMASKS_TO_LOGICAL_2D_ARGUMENTS, _BITMASKS_TO_LOGICAL_2D_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    flags.flags.writeable = False

    return flags
