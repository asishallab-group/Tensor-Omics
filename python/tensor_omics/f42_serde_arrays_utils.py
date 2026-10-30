r"""f42_serde_arrays_utils

Module for array utilities.

Defines the shared on-disk binary layout used by all typed array
serialize/deserialize modules (int/real/complex/logical/char) and the
header read/write/validate helpers that implement it. The file header is
a fixed sequence of unformatted stream records, written and read in this
order: magic number (``ARRAY_FILE_MAGIC``),
type code, number of dimensions `ndim`, then `ndim` dimension sizes. The
raw array payload follows immediately after the header, written as one
contiguous block by the type-specific serializers.

The header does NOT record the width of an element, so the payload is only
readable by a build that agrees with the writer on the storage size of the
type code. Logical arrays changed width once: they are now written as
`logical(c_bool)`, one byte per element, where earlier builds wrote the
default logical kind at four. A logical `.bin` file written before that
change therefore decodes to garbage here rather than failing, and one
written here does the same there. No other type code has moved.

Python binding, generated from f42_serde_arrays_utils. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.get_array_metadata_c.restype = None
_lib.get_array_metadata_c.argtypes = (
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GET_ARRAY_METADATA_ARGUMENTS = ("filename", "dims_out", "dims_out_capacity", "ndims", "type_code", "ierr",)
#: For a derived argument, the one the caller passed it in
_GET_ARRAY_METADATA_ARGUMENT_SOURCES = (None, None, "dims_out", None, None, None,)

def get_array_metadata(
        filename,
        dims_out_capacity,
):
    r"""Get the metadata of an array file

    Parameters
    ----------
    filename : str
        Name of the file
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    dims_out_capacity : int
        Capacity of the dims_out array

    Returns
    -------
    dict
        with keys:

        dims_out : np.ndarray[np.int32] of shape (dims_out_capacity,), read-only
            Array to store output dimensions
            The first `ndims` elements will hold the results.
            A result is a value; call `.copy()` to obtain a modifiable array.
        type_code : int
            Type code of the serialized array

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `f42_serde_arrays_utils::get_array_metadata`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    filename = _np.array([_contract.string(filename, "filename").ljust(1)], dtype="S")
    dims_out_capacity = _contract.integer(dims_out_capacity, "dims_out_capacity")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize

    # outputs and work arrays, which the caller never sees
    dims_out = _np.empty((dims_out_capacity,), dtype=_np.int32, order='C')
    ndims = _ctypes.c_int(0)
    type_code = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.get_array_metadata_c(
        filename,
        _ctypes.byref(_ctypes.c_int(filename_strlen)),
        dims_out,
        _ctypes.byref(_ctypes.c_int(dims_out_capacity)),
        _ctypes.byref(ndims),
        _ctypes.byref(type_code),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GET_ARRAY_METADATA_ARGUMENTS, _GET_ARRAY_METADATA_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    dims_out.flags.writeable = False

    return {
        "dims_out": dims_out[..., :ndims.value],
        "type_code": type_code.value,
    }
