"""tox_shape_truthful_clustering_accept

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_shape_truthful_clustering_accept. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.accept_ensemble_c.restype = None
_lib.accept_ensemble_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ACCEPT_ENSEMBLE_ARGUMENTS = ("n_dimensions", "U_t", "d_t", "G_t", "U_tp1", "d_tp1", "G_tp1", "alpha_max", "d_max", "G_max", "is_accepted", "ierr",)
#: For a derived argument, the one the caller passed it in
_ACCEPT_ENSEMBLE_ARGUMENT_SOURCES = ("U_t", None, None, None, None, None, None, None, None, None, None, None,)

def accept_ensemble(
        U_t,
        d_t,
        G_t,
        U_tp1,
        d_tp1,
        G_tp1,
        alpha_max,
        d_max,
        G_max,
):
    r"""Whether a grown ensemble at t+1 is still compatible with its own state at t

    Parameters
    ----------
    U_t : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions,), column-major (order='F')
        Ensemble's tangent+normal basis at t, see `observable`
    d_t : int
        Ensemble's intrinsic dimension at t
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    G_t : float
        Ensemble's spectral gap at t
        The minimum valid value is `above(0.0)`.
    U_tp1 : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions,), column-major (order='F')
        Ensemble's tangent+normal basis at t+1
    d_tp1 : int
        Ensemble's intrinsic dimension at t+1
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    G_tp1 : float
        Ensemble's spectral gap at t+1
        The minimum valid value is `above(0.0)`.
    alpha_max : float
        Maximum tolerated principal angle (radians)
        The minimum valid value is `0.0`.
        The maximum valid value is `2.0 * atan(1.0)`.
    d_max : int
        Maximum tolerated change in intrinsic dimension
        The minimum valid value is `0`.
    G_max : float
        Maximum tolerated |log(G_tp1/G_t)|
        The minimum valid value is `0.0`.

    Returns
    -------
    is_accepted : bool
        True if all three acceptance criteria are satisfied

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shape_truthful_clustering_accept::accept_ensemble_alloc`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        U_t = np.asfortranarray(U_t, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'U_t' must be an array of np.float64: {error}") from None
    if U_t.ndim != 2:
        raise ValueError(f"'U_t' must have 2 dimensions, but has {U_t.ndim}")
    try:
        U_tp1 = np.asfortranarray(U_tp1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'U_tp1' must be an array of np.float64: {error}") from None
    if U_tp1.ndim != 2:
        raise ValueError(f"'U_tp1' must have 2 dimensions, but has {U_tp1.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = U_t.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if U_tp1.shape[0] != n_dimensions:
        raise ValueError(f"'U_tp1' has {U_tp1.shape[0]} along axis 0, but "
            f"'U_t' implies n_dimensions == {n_dimensions}"
        )

    # outputs and work arrays, which the caller never sees
    is_accepted = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.accept_ensemble_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        U_t,
        ctypes.byref(ctypes.c_int(d_t)),
        ctypes.byref(ctypes.c_double(G_t)),
        U_tp1,
        ctypes.byref(ctypes.c_int(d_tp1)),
        ctypes.byref(ctypes.c_double(G_tp1)),
        ctypes.byref(ctypes.c_double(alpha_max)),
        ctypes.byref(ctypes.c_int(d_max)),
        ctypes.byref(ctypes.c_double(G_max)),
        ctypes.byref(is_accepted),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ACCEPT_ENSEMBLE_ARGUMENTS, _ACCEPT_ENSEMBLE_ARGUMENT_SOURCES)

    return is_accepted.value
