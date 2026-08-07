"""tox_shape_truthful_clustering_accept_kernel

# Shape Truthful Clustering (STC): Accept

Python binding, generated from tox_shape_truthful_clustering_accept_kernel. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.tox_stc_accept_ensemble_svd_workspace_c.restype = None
_lib.tox_stc_accept_ensemble_svd_workspace_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_TOX_STC_ACCEPT_ENSEMBLE_SVD_WORKSPACE_ARGUMENTS = ("d_t", "d_tp1", "lwork",)

def tox_stc_accept_ensemble_svd_workspace(
        d_t,
        d_tp1,
):
    r"""Recommend LAPACK dgesvd workspace size for accept_ensemble's principal-angle SVD

    Parameters
    ----------
    d_t : int
        Ensemble's intrinsic dimension at t
    d_tp1 : int
        Ensemble's intrinsic dimension at t+1

    Returns
    -------
    lwork : int
        Recommended size of the real LAPACK workspace

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_accept_kernel::tox_stc_accept_ensemble_svd_workspace`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    lwork = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.tox_stc_accept_ensemble_svd_workspace_c(
        ctypes.byref(ctypes.c_int(d_t)),
        ctypes.byref(ctypes.c_int(d_tp1)),
        ctypes.byref(lwork),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _TOX_STC_ACCEPT_ENSEMBLE_SVD_WORKSPACE_ARGUMENTS)

    return lwork.value
