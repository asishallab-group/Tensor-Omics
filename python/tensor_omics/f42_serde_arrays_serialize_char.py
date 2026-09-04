r"""f42_serde_arrays_serialize_char

Module for serializing character arrays into files

Python binding, generated from f42_serde_arrays_serialize_char. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.serialize_char_helper_c.restype = None
_lib.serialize_char_helper_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_CHAR_HELPER_ARGUMENTS = ("arr", "n_strings", "arr_shape", "filename", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_CHAR_HELPER_ARGUMENT_SOURCES = (None, "arr", "arr", None, None,)

def serialize_char_helper(
        arr,
        filename,
):
    r"""Subroutine to serialize a flat character array into a file

    Parameters
    ----------
    arr : sequence of str, of length *
        Array to be serialized
    filename : str
        Name of the file to write to

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_serde_arrays_serialize_char::serialize_char_helper`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    arr_shape = np.ascontiguousarray(np.shape(arr), dtype=np.int32)
    try:
        _arr_bytes = [str(_s).encode() for _s in np.asarray(arr).ravel(order='F')]
        _arr_width = max(map(len, _arr_bytes), default=0) or 1
        arr = np.asarray([_b.ljust(_arr_width) for _b in _arr_bytes], dtype="S")
    except TypeError as error:
        raise TypeError(f"'arr' must be a sequence of strings: {error}") from None
    filename = np.array([str(filename).encode().ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    arr_strlen = arr.itemsize
    n_strings = arr.size
    n_arr_shape_elements = arr_shape.shape[0]
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.serialize_char_helper_c(
        arr,
        ctypes.byref(ctypes.c_int(arr_strlen)),
        ctypes.byref(ctypes.c_int(n_strings)),
        arr_shape,
        ctypes.byref(ctypes.c_int(n_arr_shape_elements)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SERIALIZE_CHAR_HELPER_ARGUMENTS, _SERIALIZE_CHAR_HELPER_ARGUMENT_SOURCES)

    return None
