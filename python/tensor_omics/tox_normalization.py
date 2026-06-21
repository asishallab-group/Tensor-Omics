from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def normalize_unit_length(
        vector
        ):
    """
    Parameters
    ----------
    vector : np.ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F'), modified in-place
        Vector that will be normalized to unit length

    Returns
    -------
    None

    Notes
    -----
    Normalizes an input vector to unit length in-place
    """

    # ensure all array inputs are numpy arrays
    assert type(vector) is np.ndarray and vector.flags.f_contiguous and vector.dtype == np.float64, "'vector' must be column-major numpy array (order='F')"

    # extract dimension arguments
    n_dims = vector.shape[0]


    # Create temporaries and/or outputs
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.normalize_unit_length_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.normalize_unit_length_c.restype = None

    tox.normalize_unit_length_c(
        vector,
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    vector.setflags(write=False)

    return None
