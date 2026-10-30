r"""f42_serde_arrays_deserialize_char

Module for deserializing character arrays from files

Python binding, generated from f42_serde_arrays_deserialize_char. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.deserialize_char_helper_c.restype = None
_lib.deserialize_char_helper_c.argtypes = (
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
_DESERIALIZE_CHAR_HELPER_ARGUMENTS = ("arr", "n_strings", "strlen", "arr_shape", "filename", "ierr",)
#: For a derived argument, the one the caller passed it in
_DESERIALIZE_CHAR_HELPER_ARGUMENT_SOURCES = (None, "arr", None, "arr", None, None,)
#: For an array the C wrapper compares with its declaration, that shape
_DESERIALIZE_CHAR_HELPER_ARGUMENT_SHAPES = ("n_strings", None, None, None, None, None,)

def deserialize_char_helper(
        filename,
):
    r"""Subroutine to deserialize a flat character array from a file

    Parameters
    ----------
    filename : str
        Name of the file
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

    Returns
    -------
    arr : sequence of str, of length *
        Pre-allocated array to read the data into

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_serde_arrays_deserialize_char::deserialize_char_helper`, whose argument names are
    the ones an error message reports.
    """
    # kept before conversion, for the producers called below
    _filename_raw = filename

    # check every input by the input contract, converting only when C needs it
    filename = _np.array([_contract.string(filename, "filename").ljust(1)], dtype="S")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .f42_serde_arrays_utils import get_array_metadata
    _get_array_metadata_result = get_array_metadata(filename=_filename_raw, dims_out_capacity=5)
    strlen = _get_array_metadata_result["type_code"]
    from .f42_serde_arrays_utils import get_array_metadata
    arr_shape = _get_array_metadata_result["dims_out"]

    # what the inputs already say, rather than asking for it again
    n_strings = int(_np.prod(arr_shape))
    n_arr_shape_elements = arr_shape.shape[0]

    # outputs and work arrays, which the caller never sees
    arr = _np.zeros((int(_np.prod(arr_shape)),), dtype=f"S{strlen}")
    ierr = _ctypes.c_int(0)

    _lib.deserialize_char_helper_c(
        arr,
        _ctypes.byref(_ctypes.c_int(n_strings)),
        _ctypes.byref(_ctypes.c_int(strlen)),
        arr_shape,
        _ctypes.byref(_ctypes.c_int(n_arr_shape_elements)),
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _DESERIALIZE_CHAR_HELPER_ARGUMENTS, _DESERIALIZE_CHAR_HELPER_ARGUMENT_SOURCES, shapes=_DESERIALIZE_CHAR_HELPER_ARGUMENT_SHAPES)

    return _np.asarray([_s.decode().rstrip(' ') for _s in arr], dtype=str).reshape(tuple(arr_shape), order='F')
