"""tox_shape_truthful_clustering_parameter_estimation_kernel

# Shape Truthful Clustering (STC): Parameter Estimation

Python binding, generated from tox_shape_truthful_clustering_parameter_estimation_kernel. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.tox_stc_estimate_parameters_svd_workspace_c.restype = None
_lib.tox_stc_estimate_parameters_svd_workspace_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_TOX_STC_ESTIMATE_PARAMETERS_SVD_WORKSPACE_ARGUMENTS = ("n_dimensions", "n_vectors", "lwork_observable", "iwork_size", "lwork_angle",)

def tox_stc_estimate_parameters_svd_workspace(
        n_dimensions,
        n_vectors,
):
    r"""Recommend LAPACK workspace sizes for estimate_stc_parameters' SVD calls

    Parameters
    ----------
    n_dimensions : int
        Ambient dimension D
    n_vectors : int
        Number of input vectors N

    Returns
    -------
    dict
        with keys:

        lwork_observable : int
            Recommended size of observable's real LAPACK workspace (worst case)
        iwork_size : int
            Recommended size of observable's integer LAPACK workspace (worst case)
        lwork_angle : int
            Recommended size of the pairwise principal-angle LAPACK workspace (worst case)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_parameter_estimation_kernel::tox_stc_estimate_parameters_svd_workspace`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    lwork_observable = ctypes.c_int(0)
    iwork_size = ctypes.c_int(0)
    lwork_angle = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.tox_stc_estimate_parameters_svd_workspace_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(lwork_observable),
        ctypes.byref(iwork_size),
        ctypes.byref(lwork_angle),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _TOX_STC_ESTIMATE_PARAMETERS_SVD_WORKSPACE_ARGUMENTS)

    return {
        "lwork_observable": lwork_observable.value,
        "iwork_size": iwork_size.value,
        "lwork_angle": lwork_angle.value,
    }
