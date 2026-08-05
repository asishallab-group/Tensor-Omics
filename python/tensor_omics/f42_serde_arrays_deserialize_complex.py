"""f42_serde_arrays_deserialize_complex

Module for deserializing complex arrays from files

Python binding, generated from f42_serde_arrays_deserialize_complex. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.deserialize_complex_helper_c.restype = None
_lib.deserialize_complex_helper_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DESERIALIZE_COMPLEX_HELPER_ARGUMENTS = ("arr", "n_elements", "arr_shape", "filename", "ierr",)
#: For a derived argument, the one the caller passed it in
_DESERIALIZE_COMPLEX_HELPER_ARGUMENT_SOURCES = (None, "arr", "arr", None, None,)

def deserialize_complex_helper(
        filename,
):
    r"""Deserialize a flat complex array from a file

    Parameters
    ----------
    filename : str
        Name of the file

    Returns
    -------
    arr : np.ndarray[np.complex128] of shape (*,)
        Pre-allocated array to read the data into

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_serde_arrays_deserialize_complex::deserialize_complex_helper`, whose argument names are
    the ones an error message reports.
    """
    # kept before conversion, for the producers called below
    _filename_raw = filename

    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .f42_serde_arrays_utils import get_array_metadata
    _get_array_metadata_result = get_array_metadata(filename=_filename_raw, dims_out_capacity=5)
    arr_shape = _get_array_metadata_result["dims_out"]

    # what the inputs already say, rather than asking for it again
    n_elements = int(np.prod(arr_shape))
    n_arr_shape_elements = arr_shape.shape[0]

    # outputs and work arrays, which the caller never sees
    arr = np.empty((int(np.prod(arr_shape)),), dtype=np.complex128, order='C')
    ierr = ctypes.c_int(0)

    _lib.deserialize_complex_helper_c(
        arr.ravel(order='F'),
        ctypes.byref(ctypes.c_int(n_elements)),
        arr_shape,
        ctypes.byref(ctypes.c_int(n_arr_shape_elements)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DESERIALIZE_COMPLEX_HELPER_ARGUMENTS, _DESERIALIZE_COMPLEX_HELPER_ARGUMENT_SOURCES)

    return arr.reshape(tuple(arr_shape), order='F')
