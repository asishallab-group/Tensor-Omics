"""Python binding to Module for array utilities.

Generated from f42_serde_arrays_utils. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.get_array_metadata_c.restype = None
_lib.get_array_metadata_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GET_ARRAY_METADATA_ARGUMENTS = ("filename", "dims_out", "dims_out_capacity", "ndims", "type_code", "ierr",)

def get_array_metadata(
        filename,
        dims_out_capacity,
):
    r"""Get the metadata of an array file

    Parameters
    ----------
    filename : str
        Name of the file
    dims_out_capacity : int
        Capacity of the dims_out array

    Returns
    -------
    dict
        with keys:

        dims_out : np.ndarray[np.int32] of shape (dims_out_capacity,)
            Array to store output dimensions
            The first `ndims` elements will hold the results.
        type_code : int
            Type code of the serialized array

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_serde_arrays_utils::get_array_metadata`.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    dims_out = np.empty((dims_out_capacity,), dtype=np.int32, order='C')
    ndims = ctypes.c_int(0)
    type_code = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.get_array_metadata_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        dims_out,
        ctypes.byref(ctypes.c_int(dims_out_capacity)),
        ctypes.byref(ndims),
        ctypes.byref(type_code),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GET_ARRAY_METADATA_ARGUMENTS)

    return {
        "dims_out": dims_out[..., :ndims.value],
        "type_code": type_code.value,
    }
