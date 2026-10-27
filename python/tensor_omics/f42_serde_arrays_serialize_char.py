r"""f42_serde_arrays_serialize_char

Module for serializing character arrays into files

Python binding, generated from f42_serde_arrays_serialize_char. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.serialize_char_helper_c.restype = None
_lib.serialize_char_helper_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_serde_arrays_serialize_char::serialize_char_helper`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    arr_shape = _np.ascontiguousarray(_np.shape(arr), dtype=_np.int32)
    try:
        _arr_bytes = [str(_s).encode() for _s in _np.asarray(arr).ravel(order='F')]
        _arr_width = max(map(len, _arr_bytes), default=0) or 1
        arr = _np.asarray([_b.ljust(_arr_width) for _b in _arr_bytes], dtype="S")
    except TypeError as _error:
        raise _ToxArgumentError(f"'arr' must be a sequence of strings: {_error}", argument="arr") from None
    filename = _np.array([str(filename).encode().ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    arr_strlen = arr.itemsize
    n_strings = arr.size
    n_arr_shape_elements = arr_shape.shape[0]
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.serialize_char_helper_c(
        arr,
        _ctypes.byref(_ctypes.c_int(arr_strlen)),
        _ctypes.byref(_ctypes.c_int(n_strings)),
        arr_shape,
        _ctypes.byref(_ctypes.c_int(n_arr_shape_elements)),
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _SERIALIZE_CHAR_HELPER_ARGUMENTS, _SERIALIZE_CHAR_HELPER_ARGUMENT_SOURCES)

    return None
