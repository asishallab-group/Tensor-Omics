"""f42_serde_arrays_serialize_complex

Module for serializing complex arrays into files

Python binding, generated from f42_serde_arrays_serialize_complex. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.serialize_complex_helper_c.restype = None
_lib.serialize_complex_helper_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_COMPLEX_HELPER_ARGUMENTS = ("arr", "n_elements", "arr_shape", "filename", "ierr",)

def serialize_complex_helper(
        arr,
        filename,
):
    r"""Subroutine to serialize a flat complex array into a file

    Parameters
    ----------
    arr : np.ndarray[np.complex128] of shape (*,)
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
    Generated from the Fortran procedure `f42_serde_arrays_serialize_complex::serialize_complex_helper`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    arr_shape = np.ascontiguousarray(np.shape(arr), dtype=np.int32)
    try:
        arr = np.asfortranarray(arr, dtype=np.complex128)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'arr' must be an array of np.complex128: {error}") from None
    filename = np.array([str(filename).encode()], dtype="S")

    # what the inputs already say, rather than asking for it again
    n_elements = arr.size
    n_arr_shape_elements = arr_shape.shape[0]
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.serialize_complex_helper_c(
        arr.ravel(order='F'),
        ctypes.byref(ctypes.c_int(n_elements)),
        arr_shape,
        ctypes.byref(ctypes.c_int(n_arr_shape_elements)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SERIALIZE_COMPLEX_HELPER_ARGUMENTS)

    return None
