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
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_ACCEPT_ENSEMBLE_ARGUMENTS = ("n_dimensions", "o", "U_first", "d_first", "U_history", "d_history", "history_len", "G_t", "normal_error_t", "U_tp1", "d_tp1", "G_tp1", "normal_error_tp1", "chordal_dist_max_as_prcnt_of_range", "d_max", "G_max", "RMSE_change_max", "is_accepted", "ierr",)
#: For a derived argument, the one the caller passed it in
_ACCEPT_ENSEMBLE_ARGUMENT_SOURCES = ("U_first", "U_history", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,)

def accept_ensemble(
        U_first,
        d_first,
        U_history,
        d_history,
        history_len,
        G_t,
        normal_error_t,
        U_tp1,
        d_tp1,
        G_tp1,
        normal_error_tp1,
        chordal_dist_max_as_prcnt_of_range,
        d_max,
        G_max,
        RMSE_change_max,
):
    r"""Whether a grown ensemble at t+1 is still compatible with its own growth trajectory

    Parameters
    ----------
    U_first : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions,), column-major (order='F')
        Ensemble's tangent+normal basis at its bootstrap iteration (iteration 1), see
        `misc/mod_STC.md`, "Output"
    d_first : int
        Ensemble's intrinsic dimension at its bootstrap iteration
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    U_history : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions, o,), column-major (order='F')
        Trailing tangent+normal bases, oldest to newest; only columns 1:history_len are
        valid (see `history_len`)
    d_history : np.ndarray[np.int32] of shape (o,)
        Trailing intrinsic dimensions, one per column of U_history; only entries
        1:history_len are valid
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    history_len : int
        Number of valid columns in U_history/d_history; column history_len is the most
        recently accepted iteration
        The minimum valid value is `1`.
        The maximum valid value is `o`.
    G_t : float
        Ensemble's spectral gap at the most recently accepted iteration
        The minimum valid value is `above(0.0)`.
    normal_error_t : float
        Ensemble's normal_error at the most recently accepted iteration, see
        `normal_error`. Zero is valid -- a perfectly flat/collinear ensemble has no
        noise at all in its normal directions -- unlike G_t, which is a ratio already
        protected by its own +epsilon denominator (see `observable`) and so is
        required to be strictly positive.
        The minimum valid value is `0.0`.
    U_tp1 : np.ndarray[np.float64] of shape (n_dimensions, n_dimensions,), column-major (order='F')
        Candidate ensemble's tangent+normal basis
    d_tp1 : int
        Candidate ensemble's intrinsic dimension
        The minimum valid value is `0`.
        The maximum valid value is `n_dimensions`.
    G_tp1 : float
        Candidate ensemble's spectral gap
        The minimum valid value is `above(0.0)`.
    normal_error_tp1 : float
        Candidate ensemble's normal_error. Zero is valid, see `normal_error_t`.
        The minimum valid value is `0.0`.
    chordal_dist_max_as_prcnt_of_range : float
        Maximum tolerated chordal distance between tangent bases, as a fraction of its
        own [0, sqrt(d)] range, see `accept_ensemble`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    d_max : int
        Maximum tolerated change in intrinsic dimension, see `accept_ensemble`
        The minimum valid value is `0`.
    G_max : float
        Maximum tolerated |log(G_tp1/G_t)|
        The minimum valid value is `0.0`.
    RMSE_change_max : float
        Maximum tolerated |log(RMSE_tp1/RMSE_t)|
        The minimum valid value is `0.0`.

    Returns
    -------
    is_accepted : bool
        True if all four acceptance criteria are satisfied

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
        U_first = np.asfortranarray(U_first, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'U_first' must be an array of np.float64: {error}") from None
    if U_first.ndim != 2:
        raise ValueError(f"'U_first' must have 2 dimensions, but has {U_first.ndim}")
    try:
        U_history = np.asfortranarray(U_history, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'U_history' must be an array of np.float64: {error}") from None
    if U_history.ndim != 3:
        raise ValueError(f"'U_history' must have 3 dimensions, but has {U_history.ndim}")
    try:
        d_history = np.ascontiguousarray(d_history, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'d_history' must be an array of np.int32: {error}") from None
    if d_history.ndim != 1:
        raise ValueError(f"'d_history' must have 1 dimension, but has {d_history.ndim}")
    try:
        U_tp1 = np.asfortranarray(U_tp1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'U_tp1' must be an array of np.float64: {error}") from None
    if U_tp1.ndim != 2:
        raise ValueError(f"'U_tp1' must have 2 dimensions, but has {U_tp1.ndim}")

    # what the inputs already say, rather than asking for it again
    n_dimensions = U_first.shape[0]
    o = U_history.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if U_history.shape[0] != n_dimensions:
        raise ValueError(f"'U_history' has {U_history.shape[0]} along axis 0, but "
            f"'U_first' implies n_dimensions == {n_dimensions}"
        )
    if U_tp1.shape[0] != n_dimensions:
        raise ValueError(f"'U_tp1' has {U_tp1.shape[0]} along axis 0, but "
            f"'U_first' implies n_dimensions == {n_dimensions}"
        )
    if d_history.shape[0] != o:
        raise ValueError(f"'d_history' has {d_history.shape[0]} along axis 0, but "
            f"'U_history' implies o == {o}"
        )

    # outputs and work arrays, which the caller never sees
    is_accepted = ctypes.c_bool(0)
    ierr = ctypes.c_int(0)

    _lib.accept_ensemble_c(
        ctypes.byref(ctypes.c_int(n_dimensions)),
        ctypes.byref(ctypes.c_int(o)),
        U_first,
        ctypes.byref(ctypes.c_int(d_first)),
        U_history,
        d_history,
        ctypes.byref(ctypes.c_int(history_len)),
        ctypes.byref(ctypes.c_double(G_t)),
        ctypes.byref(ctypes.c_double(normal_error_t)),
        U_tp1,
        ctypes.byref(ctypes.c_int(d_tp1)),
        ctypes.byref(ctypes.c_double(G_tp1)),
        ctypes.byref(ctypes.c_double(normal_error_tp1)),
        ctypes.byref(ctypes.c_double(chordal_dist_max_as_prcnt_of_range)),
        ctypes.byref(ctypes.c_int(d_max)),
        ctypes.byref(ctypes.c_double(G_max)),
        ctypes.byref(ctypes.c_double(RMSE_change_max)),
        ctypes.byref(is_accepted),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _ACCEPT_ENSEMBLE_ARGUMENTS, _ACCEPT_ENSEMBLE_ARGUMENT_SOURCES)

    return is_accepted.value
