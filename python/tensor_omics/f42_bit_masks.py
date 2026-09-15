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
way, and the counts, ``bit_mask_all`` and
``bit_mask_equal`` rely on it.

**Several masks** are the columns of a `(n_words, n_masks)` array. Every column starts on a
new word, so `bit_masks(:, i_mask)` is a mask of its own and passes to any procedure here
without a copy.

**Parallel loops write whole words.** Two iterations that set flags in the same word race, so
a `do concurrent` over flags may only read a mask; one that writes runs over words or columns.

**Hot loops** use the whole-mask procedures, or step from one set flag to the next with
``bit_mask_next``. A procedure call per flag costs more
than the flag, because nothing is inlined across modules; where a loop must touch single
flags, the bit-mask macros in `src/macros.h` test a flag or locate its word inline, and one
of them sizes a mask.

Python binding, generated from f42_bit_masks. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.bit_mask_n_words_c.restype = None
_lib.bit_mask_n_words_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BIT_MASK_N_WORDS_ARGUMENTS = ("n_bits", "n_words", "ierr",)

_lib.bit_mask_test_c.restype = None
_lib.bit_mask_test_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BIT_MASK_TEST_ARGUMENTS = ("n_bits", "n_words", "bit_mask", "i_bit", "is_set", "ierr",)
#: For a derived argument, the one the caller passed it in
_BIT_MASK_TEST_ARGUMENT_SOURCES = (None, "bit_mask", None, None, None, None,)

_lib.bit_mask_from_logical_c.restype = None
_lib.bit_mask_from_logical_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BIT_MASK_FROM_LOGICAL_ARGUMENTS = ("n_bits", "n_words", "flags", "bit_mask", "ierr",)
#: For a derived argument, the one the caller passed it in
_BIT_MASK_FROM_LOGICAL_ARGUMENT_SOURCES = ("flags", "bit_mask", None, None, None,)

_lib.bit_mask_to_logical_c.restype = None
_lib.bit_mask_to_logical_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BIT_MASK_TO_LOGICAL_ARGUMENTS = ("n_bits", "n_words", "bit_mask", "flags", "ierr",)
#: For a derived argument, the one the caller passed it in
_BIT_MASK_TO_LOGICAL_ARGUMENT_SOURCES = ("flags", "bit_mask", None, None, None,)

_lib.bit_masks_from_logical_2D_c.restype = None
_lib.bit_masks_from_logical_2D_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BIT_MASKS_FROM_LOGICAL_2D_ARGUMENTS = ("n_bits", "n_masks", "n_words", "flags", "bit_masks", "ierr",)
#: For a derived argument, the one the caller passed it in
_BIT_MASKS_FROM_LOGICAL_2D_ARGUMENT_SOURCES = ("flags", "flags", "bit_masks", None, None, None,)

_lib.bit_masks_to_logical_2D_c.restype = None
_lib.bit_masks_to_logical_2D_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BIT_MASKS_TO_LOGICAL_2D_ARGUMENTS = ("n_bits", "n_masks", "n_words", "bit_masks", "flags", "ierr",)
#: For a derived argument, the one the caller passed it in
_BIT_MASKS_TO_LOGICAL_2D_ARGUMENT_SOURCES = ("flags", "bit_masks", "bit_masks", None, None, None,)

def bit_mask_n_words(
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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bit_mask_n_words`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    n_words = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.bit_mask_n_words_c(
        ctypes.byref(ctypes.c_int(n_bits)),
        ctypes.byref(n_words),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BIT_MASK_N_WORDS_ARGUMENTS)

    return n_words.value

def bit_mask_test(
        n_bits,
        bit_mask,
        i_bit,
):
    r"""Whether flag `i_bit` of a mask is set.

    Inside Fortran, a loop tests a flag inline with a bit-mask macro from `src/macros.h` instead.

    Parameters
    ----------
    n_bits : int
        number of flags in the mask
        The minimum valid value is `1`.
    bit_mask : np.ndarray[np.int32] of shape (n_words,)
        the mask
    i_bit : int
        index of the flag, counting from 1
        The minimum valid value is `1`.
        The maximum valid value is `n_bits`.

    Returns
    -------
    is_set : bool
        `True` if the flag is set; `False` for `i_bit > n_bits`, which the bits past
        `n_bits` would say too, as long as the mask keeps them zero

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bit_mask_test`, whose argument names are
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
    n_words = bit_mask.shape[0]

    # outputs and work arrays, which the caller never sees
    is_set = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.bit_mask_test_c(
        ctypes.byref(ctypes.c_int(n_bits)),
        ctypes.byref(ctypes.c_int(n_words)),
        bit_mask,
        ctypes.byref(ctypes.c_int(i_bit)),
        ctypes.byref(is_set),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BIT_MASK_TEST_ARGUMENTS, _BIT_MASK_TEST_ARGUMENT_SOURCES)

    return is_set.value

def bit_mask_from_logical(
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
    bit_mask : np.ndarray[np.int32] of shape (n_words,), read-only
        the mask, flag `i` set where `flags(i)` is `True`
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bit_mask_from_logical`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        flags = np.ascontiguousarray(flags, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'flags' must be an array of np.bool_: {error}") from None
    if flags.ndim != 1:
        raise ValueError(f"'flags' must have 1 dimension, but has {flags.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bits = flags.shape[0]

    # outputs and work arrays, which the caller never sees
    bit_mask = np.empty((n_words,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.bit_mask_from_logical_c(
        ctypes.byref(ctypes.c_int(n_bits)),
        ctypes.byref(ctypes.c_int(n_words)),
        flags,
        bit_mask,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BIT_MASK_FROM_LOGICAL_ARGUMENTS, _BIT_MASK_FROM_LOGICAL_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    bit_mask.flags.writeable = False

    return bit_mask

def bit_mask_to_logical(
        n_bits,
        bit_mask,
):
    r"""Unpacks a bit mask into a `logical(c_bool)` array.

    Parameters
    ----------
    n_bits : int
        number of flags
        The minimum valid value is `1`.
    bit_mask : np.ndarray[np.int32] of shape (n_words,)
        the mask

    Returns
    -------
    flags : np.ndarray[np.bool_] of shape (n_bits,), read-only
        `flags(i)` is `True` where flag `i` of the mask is set
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bit_mask_to_logical`, whose argument names are
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
    n_words = bit_mask.shape[0]

    # outputs and work arrays, which the caller never sees
    flags = np.empty((n_bits,), dtype=np.bool_, order='C')
    ierr = ctypes.c_int(0)

    _lib.bit_mask_to_logical_c(
        ctypes.byref(ctypes.c_int(n_bits)),
        ctypes.byref(ctypes.c_int(n_words)),
        bit_mask,
        flags,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BIT_MASK_TO_LOGICAL_ARGUMENTS, _BIT_MASK_TO_LOGICAL_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    flags.flags.writeable = False

    return flags

def bit_masks_from_logical_2D(
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
    bit_masks : np.ndarray[np.int32] of shape (n_words, n_masks,), column-major (order='F'), read-only
        the masks, column `j` packed from `flags(:, j)`
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bit_masks_from_logical_2D`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        flags = np.asfortranarray(flags, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'flags' must be an array of np.bool_: {error}") from None
    if flags.ndim != 2:
        raise ValueError(f"'flags' must have 2 dimensions, but has {flags.ndim}")

    # what the inputs already say, rather than asking for it again
    n_bits = flags.shape[0]
    n_masks = flags.shape[1]

    # outputs and work arrays, which the caller never sees
    bit_masks = np.empty((n_words, n_masks,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.bit_masks_from_logical_2D_c(
        ctypes.byref(ctypes.c_int(n_bits)),
        ctypes.byref(ctypes.c_int(n_masks)),
        ctypes.byref(ctypes.c_int(n_words)),
        flags,
        bit_masks,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BIT_MASKS_FROM_LOGICAL_2D_ARGUMENTS, _BIT_MASKS_FROM_LOGICAL_2D_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    bit_masks.flags.writeable = False

    return bit_masks

def bit_masks_to_logical_2D(
        n_bits,
        bit_masks,
):
    r"""Unpacks every column of a matrix of bit masks into a `logical(c_bool)` column.

    Parameters
    ----------
    n_bits : int
        number of flags per mask
        The minimum valid value is `1`.
    bit_masks : np.ndarray[np.int32] of shape (n_words, n_masks,), column-major (order='F')
        the masks, one per column

    Returns
    -------
    flags : np.ndarray[np.bool_] of shape (n_bits, n_masks,), column-major (order='F'), read-only
        the flags, column `j` unpacked from `bit_masks(:, j)`
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_bit_masks::bit_masks_to_logical_2D`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        bit_masks = np.asfortranarray(bit_masks, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'bit_masks' must be an array of np.int32: {error}") from None
    if bit_masks.ndim != 2:
        raise ValueError(f"'bit_masks' must have 2 dimensions, but has {bit_masks.ndim}")

    # what the inputs already say, rather than asking for it again
    n_masks = bit_masks.shape[1]
    n_words = bit_masks.shape[0]

    # outputs and work arrays, which the caller never sees
    flags = np.empty((n_bits, n_masks,), dtype=np.bool_, order='F')
    ierr = ctypes.c_int(0)

    _lib.bit_masks_to_logical_2D_c(
        ctypes.byref(ctypes.c_int(n_bits)),
        ctypes.byref(ctypes.c_int(n_masks)),
        ctypes.byref(ctypes.c_int(n_words)),
        bit_masks,
        flags,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _BIT_MASKS_TO_LOGICAL_2D_ARGUMENTS, _BIT_MASKS_TO_LOGICAL_2D_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    flags.flags.writeable = False

    return flags
