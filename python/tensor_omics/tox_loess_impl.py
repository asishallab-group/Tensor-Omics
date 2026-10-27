r"""tox_loess_impl

LOESS local polynomial regression smoothing, over the netlib `dloess`/`lowesd` family.

Two fits: `loess_fit_plain` and `loess_fit_robust`, the latter reweighting against outliers over
a number of iterations. There is deliberately no combined entry point dispatching on a mode --
a caller chooses the routine, and supplies the weights and the evaluation points it wants.

A sample too degenerate to fit (fewer distinct points than the degree needs) is answered
directly from the observations rather than handed to netlib, which cannot fit one and dies
inside its decomposition instead of reporting. Nothing about that is visible to a caller
beyond the result.

`tox_loess_required_workspace` sizes the workspace for a caller that wants to own it; the
mode and iteration constants and `EPS_LOESS` are here for the same reason. The netlib
routines themselves are not re-documented beyond their calling convention.

Python binding, generated from tox_loess_impl. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.tox_loess_required_workspace_c.restype = None
_lib.tox_loess_required_workspace_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_TOX_LOESS_REQUIRED_WORKSPACE_ARGUMENTS = ("n_dim", "max_neighborhood_size", "int_workspace_size", "real_workspace_size", "save_factorization",)

def tox_loess_required_workspace(
        n_dim,
        max_neighborhood_size,
        save_factorization,
):
    r"""Recommend workspace sizes based on Netlib exact formulas

    Computes the required sizes for integer and real workspace arrays.
    These sizes depend on the dimensionality of the data and the maximum neighborhood size.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess_impl::tox_loess_required_workspace`, whose argument names are
    the ones an error message reports.
    """
    # outputs and work arrays, which the caller never sees
    int_workspace_size = _ctypes.c_int(0)
    real_workspace_size = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.tox_loess_required_workspace_c(
        _ctypes.byref(_ctypes.c_int(n_dim)),
        _ctypes.byref(_ctypes.c_int(max_neighborhood_size)),
        _ctypes.byref(int_workspace_size),
        _ctypes.byref(real_workspace_size),
        _ctypes.byref(_ctypes.c_bool(save_factorization)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _TOX_LOESS_REQUIRED_WORKSPACE_ARGUMENTS)

    return {
        "int_workspace_size": int_workspace_size.value,
        "real_workspace_size": real_workspace_size.value,
    }
