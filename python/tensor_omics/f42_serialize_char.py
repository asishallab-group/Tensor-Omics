from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def serialize_char_1d(
        arr,
        filename
        ):
    """
    Parameters
    ----------
    arr : np.ndarray[f"S{arr_strlen}"] of shape (n_arr_elements,) in column-major layout (order='F')
        array to save
    filename : str
        output filename

    Returns
    -------
    None

    Notes
    -----
    Serialize a 1D character array to a binary file.The file will contain a magic number, type code, dimension, shape, character length, and the array data.
    """

    # ensure all array inputs are numpy arrays
    arr = np.asarray(arr)
    filename = np.asarray(filename)

    # extract dimension arguments
    arr_strlen = arr.dtype.itemsize // arr.dtype.alignment
    n_arr_elements = arr.shape[0]
    filename_strlen = filename.dtype.itemsize // filename.dtype.alignment


    # Create temporaries and/or outputs
    arr = arr.astype(f"S{arr_strlen}", order="F")
    filename = filename.astype(f"S{filename_strlen}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.serialize_char_1d_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='F_CONTIGUOUS', dtype=f"S{arr_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{filename_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.serialize_char_1d_c.restype = None

    tox.serialize_char_1d_c(
        arr,
        ctypes.byref(ctypes.c_int(arr_strlen)),
        ctypes.byref(ctypes.c_int(n_arr_elements)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return None
