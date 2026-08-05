"""tox_loess_kernel

Kernels for LOESS (netlib `dloess`/`lowesd` family) local polynomial regression smoothing.

Python binding, generated from tox_loess_kernel. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.tox_loess_required_workspace_c.restype = None
_lib.tox_loess_required_workspace_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_TOX_LOESS_REQUIRED_WORKSPACE_ARGUMENTS = ("n_dim", "max_neighborhood_size", "int_workspace_size", "real_workspace_size", "save_factorization",)

def tox_loess_required_workspace(
        n_dim,
        max_neighborhood_size,
        save_factorization,
):
    r"""Recommend workspace sizes based on Netlib exact formulas

    Parameters
    ----------
    n_dim : int
        Dimensionality of the data
    max_neighborhood_size : int
        Maximum neighborhood size
    save_factorization : bool
        Save matrix factorization flag

    Returns
    -------
    dict
        with keys:

        int_workspace_size : int
            Required size of the integer workspace array
        real_workspace_size : int
            Required size of the real workspace array

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess_kernel::tox_loess_required_workspace`.
    """
    # outputs and work arrays, which the caller never sees
    int_workspace_size = ctypes.c_int(0)
    real_workspace_size = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.tox_loess_required_workspace_c(
        ctypes.byref(ctypes.c_int(n_dim)),
        ctypes.byref(ctypes.c_int(max_neighborhood_size)),
        ctypes.byref(int_workspace_size),
        ctypes.byref(real_workspace_size),
        ctypes.byref(ctypes.c_bool(save_factorization)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _TOX_LOESS_REQUIRED_WORKSPACE_ARGUMENTS)

    return {
        "int_workspace_size": int_workspace_size.value,
        "real_workspace_size": real_workspace_size.value,
    }
