r"""tox_shape_truthful_clustering_observable_impl

# Shape Truthful Clustering (STC): Observable

`observable`: the tuple (U, d, G, mu, normal_error, tangent_scales) for an ensemble,
obtained from the economy-mode singular value decomposition (LAPACK `dgesdd`) of its
centered member vectors -- never an eigendecomposition of an explicitly formed
covariance matrix (see `misc/mod_STC.md`, "Numerical Linear Algebra"). `normal_error` and
`tangent_scales` are simple, dependency-free reductions over the eigenvalues `observable`
computes. See `misc/mod_STC.md`, SKG `observable`/`normal_error`/`tangent_scales`, for the
full algorithm definitions.

Python binding, generated from tox_shape_truthful_clustering_observable_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.tox_stc_observable_svd_workspace_c.restype = None
_lib.tox_stc_observable_svd_workspace_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_TOX_STC_OBSERVABLE_SVD_WORKSPACE_ARGUMENTS = ("n_dimensions", "n_selected_member", "lwork", "iwork_size",)

def tox_stc_observable_svd_workspace(
        n_dimensions,
        n_selected_member,
):
    r"""Recommend LAPACK dgesdd workspace sizes for observable's economy-mode SVD

    Parameters
    ----------
    n_dimensions : int
        Ambient dimension D
    n_selected_member : int
        Number of selected ensemble members

    Returns
    -------
    dict
        with keys:

        lwork : int
            Recommended size of the real LAPACK workspace
        iwork_size : int
            Recommended size of the integer LAPACK workspace

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_observable_impl::tox_stc_observable_svd_workspace`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    lwork = ctypes.c_int(0)
    iwork_size = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.tox_stc_observable_svd_workspace_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(n_selected_member)),
        ctypes.byref(lwork),
        ctypes.byref(iwork_size),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _TOX_STC_OBSERVABLE_SVD_WORKSPACE_ARGUMENTS)

    return {
        "lwork": lwork.value,
        "iwork_size": iwork_size.value,
    }
