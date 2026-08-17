r"""tox_shape_truthful_clustering_accept_impl

# Shape Truthful Clustering (STC): Accept

`accept_ensemble`: whether a grown ensemble at t+1 is still compatible with its own growth
trajectory, judged by four criteria -- tangent-space drift (chordal distance, compared
against a reference set: the bootstrap iteration plus the trailing o-window, not just the
immediately preceding iteration), change in intrinsic dimension (against both the bootstrap
iteration and the immediately preceding one), relative change in spectral gap, and relative
change in residual (RMSE), both against the immediately preceding iteration only. See
`misc/mod_STC.md`, SKG `accept_ensemble`, for the full algorithm definition and the
"no cumulative-rotation budget" rationale for comparing against a reference set rather than
a single previous state. This compares the SAME ensemble across one growth step -- not two
different ensembles/anchors at a possible junction -- so, unlike
`misc/STC_for_LoManLe.md` section 4's explicit "angle never gates a junction" rule, a
tangent-space-drift mismatch here legitimately contributes to rejection.

Python binding, generated from tox_shape_truthful_clustering_accept_impl. Do not edit.
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
)

#: The wrapped procedure's arguments, so an error can name one
_TOX_STC_ACCEPT_ENSEMBLE_SVD_WORKSPACE_ARGUMENTS = ("n_dimensions", "lwork",)

def tox_stc_accept_ensemble_svd_workspace(
        n_dimensions,
):
    r"""Recommend LAPACK dgesvd workspace size for accept_ensemble's principal-angle SVDs

    Parameters
    ----------
    n_dimensions : int
        Ambient dimension D

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
    Generated from the Fortran procedure `tox_shape_truthful_clustering_accept_impl::tox_stc_accept_ensemble_svd_workspace`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    lwork = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.tox_stc_accept_ensemble_svd_workspace_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(lwork),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _TOX_STC_ACCEPT_ENSEMBLE_SVD_WORKSPACE_ARGUMENTS)

    return lwork.value
