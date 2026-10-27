r"""f42_serde_arrays_serialize_logical

Module for serializing logical arrays into files

Python binding, generated from f42_serde_arrays_serialize_logical. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.serialize_logical_helper_c.restype = None
_lib.serialize_logical_helper_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_LOGICAL_HELPER_ARGUMENTS = ("arr", "n_elements", "arr_shape", "filename", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_LOGICAL_HELPER_ARGUMENT_SOURCES = (None, "arr", "arr", None, None,)

def serialize_logical_helper(
        arr,
        filename,
):
    r"""Subroutine to serialize a flat logical array into a file

    Parameters
    ----------
    arr : np.ndarray[np.bool_] of shape (*,)
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
    Generated from the Fortran procedure `f42_serde_arrays_serialize_logical::serialize_logical_helper`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    arr_shape = _np.ascontiguousarray(_np.shape(arr), dtype=_np.int32)
    try:
        arr = _np.asfortranarray(arr, dtype=_np.bool_)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'arr' must be an array of np.bool_: {_error}", argument="arr") from None
    filename = _np.array([str(filename).encode().ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    n_elements = arr.size
    n_arr_shape_elements = arr_shape.shape[0]
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.serialize_logical_helper_c(
        arr.ravel(order='F'),
        _ctypes.byref(_ctypes.c_int(n_elements)),
        arr_shape,
        _ctypes.byref(_ctypes.c_int(n_arr_shape_elements)),
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _SERIALIZE_LOGICAL_HELPER_ARGUMENTS, _SERIALIZE_LOGICAL_HELPER_ARGUMENT_SOURCES)

    return None
