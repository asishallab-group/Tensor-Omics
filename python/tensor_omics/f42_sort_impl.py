r"""f42_sort_impl

Indirect sorting -- every routine reorders a permutation vector rather than the data -- plus the searches over a sorted array.

One of the modules ``f42_utils_impl`` gathers; `use f42_utils_impl` reaches all of them.

Python binding, generated from f42_sort_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.sort_real_get_perm_c.restype = None
_lib.sort_real_get_perm_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SORT_REAL_GET_PERM_ARGUMENTS = ("array", "n", "perm", "ierr",)
#: For a derived argument, the one the caller passed it in
_SORT_REAL_GET_PERM_ARGUMENT_SOURCES = (None, "array", None, None,)

def sort_real_get_perm(
        array,
):
    r"""Sort a real vector, returning the permutation (ascending, NaN last)

    Parameters
    ----------
    array : np.ndarray[np.float64] of shape (n,)
        values to sort

    Returns
    -------
    perm : np.ndarray[np.int32] of shape (n,), read-only
        permutation that sorts `array` ascending, NaN last
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_sort_impl::sort_real_get_perm`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        array = np.ascontiguousarray(array, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'array' must be an array of np.float64: {error}") from None
    if array.ndim != 1:
        raise ValueError(f"'array' must have 1 dimension, but has {array.ndim}")

    # what the inputs already say, rather than asking for it again
    n = array.shape[0]

    # outputs and work arrays, which the caller never sees
    perm = np.empty((n,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.sort_real_get_perm_c(
        array,
        ctypes.byref(ctypes.c_int(n)),
        perm,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SORT_REAL_GET_PERM_ARGUMENTS, _SORT_REAL_GET_PERM_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    perm.flags.writeable = False

    return perm
