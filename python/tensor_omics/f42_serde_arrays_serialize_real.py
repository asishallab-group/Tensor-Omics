r"""f42_serde_arrays_serialize_real

Module for serializing real arrays into files

Python binding, generated from f42_serde_arrays_serialize_real. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.serialize_real_helper_c.restype = None
_lib.serialize_real_helper_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_REAL_HELPER_ARGUMENTS = ("arr", "n_elements", "arr_shape", "filename", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_REAL_HELPER_ARGUMENT_SOURCES = (None, "arr", "arr", None, None,)
#: For an array the C wrapper compares with its declaration, that shape
_SERIALIZE_REAL_HELPER_ARGUMENT_SHAPES = ("n_elements", None, None, None, None,)

def serialize_real_helper(
        arr,
        filename,
):
    r"""Subroutine to serialize a flat real array into a file

    Parameters
    ----------
    arr : np.ndarray[np.float64] of shape (*,)
        Array to be serialized
    filename : str
        Name of the file to write to
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

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
    Generated from the Fortran procedure `f42_serde_arrays_serialize_real::serialize_real_helper`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    arr_shape = _np.ascontiguousarray(_np.shape(arr), dtype=_np.int32)
    arr = _contract.reals(arr, "arr", order="F")
    filename = _np.array([_contract.string(filename, "filename").ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    n_elements = arr.size
    n_arr_shape_elements = arr_shape.shape[0]
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.serialize_real_helper_c(
        arr.ravel(order='F'),
        _ctypes.byref(_ctypes.c_int(n_elements)),
        arr_shape,
        _ctypes.byref(_ctypes.c_int(n_arr_shape_elements)),
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _SERIALIZE_REAL_HELPER_ARGUMENTS, _SERIALIZE_REAL_HELPER_ARGUMENT_SOURCES, shapes=_SERIALIZE_REAL_HELPER_ARGUMENT_SHAPES)

    return None
