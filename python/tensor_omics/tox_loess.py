r"""tox_loess

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

Python binding, generated from tox_loess. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.loess_fit_plain_c.restype = None
_lib.loess_fit_plain_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_bool),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_FIT_PLAIN_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "fitted_values", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOESS_FIT_PLAIN_ARGUMENT_SOURCES = ("x", None, None, None, None, None, None, None, None, None, None, None,)

_lib.loess_fit_robust_c.restype = None
_lib.loess_fit_robust_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_FIT_ROBUST_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "n_iters", "fitted_values", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOESS_FIT_ROBUST_ARGUMENT_SOURCES = ("x", None, None, None, None, None, None, None, None, None, None, None, None,)

def loess_fit_plain(
        x,
        y,
        weights,
        eval_points,
        span,
        degree,
        max_neighborhood_size,
        compute_influence=False,
        save_factorization=False,
):
    r"""Perform plain LOESS fitting

    Data too degenerate to fit is answered directly, by the observations themselves; see
    ``loess_degenerate_fit``.
    Fits a LOESS model to the data using the specified smoothing parameter and outputs the smoothed
    response array.

    Parameters
    ----------
    x : np.ndarray[np.float64] of shape (n,)
        Predictor variable array
    y : np.ndarray[np.float64] of shape (n,)
        Response variable array
    weights : np.ndarray[np.float64] of shape (n,)
        Weight array for data points
    eval_points : np.ndarray[np.float64] of shape (n, 1,), column-major (order='F')
        Evaluation points (x values at which the fitted curve is computed)
    span : float
        Smoothing parameter for LOESS
        The minimum valid value is `EPS_LOESS`.
        The maximum valid value is `1.0`.
    degree : int
        Degree of the LOESS polynomial
        The minimum valid value is `0`.
        The maximum valid value is `2`.
    max_neighborhood_size : int
        Maximum neighborhood size
    compute_influence : bool, optional, default False
        Influence calculation flag
        The default value is `False`.
    save_factorization : bool, optional, default False
        Save matrix factorization flag
        The default value is `False`.

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (n,), read-only
        Fitted (smoothed) values of y at the evaluation points
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_plain`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x = _np.ascontiguousarray(x, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'x' must be an array of np.float64: {_error}", argument="x") from None
    if x.ndim != 1:
        raise _ToxArgumentError(f"'x' must have 1 dimension, but has {x.ndim}", argument="x")
    try:
        y = _np.ascontiguousarray(y, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'y' must be an array of np.float64: {_error}", argument="y") from None
    if y.ndim != 1:
        raise _ToxArgumentError(f"'y' must have 1 dimension, but has {y.ndim}", argument="y")
    try:
        weights = _np.ascontiguousarray(weights, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'weights' must be an array of np.float64: {_error}", argument="weights") from None
    if weights.ndim != 1:
        raise _ToxArgumentError(f"'weights' must have 1 dimension, but has {weights.ndim}", argument="weights")
    try:
        eval_points = _np.asfortranarray(eval_points, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'eval_points' must be an array of np.float64: {_error}", argument="eval_points") from None
    if eval_points.ndim != 2:
        raise _ToxArgumentError(f"'eval_points' must have 2 dimensions, but has {eval_points.ndim}", argument="eval_points")

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if y.shape[0] != n:
        raise _ToxArgumentError(f"'y' has {y.shape[0]} along axis 0, but 'x' implies n == {n}; 'y' is declared with shape (n)", argument="y")
    if weights.shape[0] != n:
        raise _ToxArgumentError(f"'weights' has {weights.shape[0]} along axis 0, but 'x' implies n == {n}; 'weights' is declared with shape (n)", argument="weights")
    if eval_points.shape[0] != n:
        raise _ToxArgumentError(f"'eval_points' has {eval_points.shape[0]} along axis 0, but 'x' implies n == {n}; 'eval_points' is declared with shape (n, 1)", argument="eval_points")
    if eval_points.shape[1] != 1:
        raise _ToxArgumentError(f"'eval_points' has {eval_points.shape[1]} along axis 1, but its declared shape fixes it at {1}; 'eval_points' is declared with shape (n, 1)", argument="eval_points")

    # outputs and work arrays, which the caller never sees
    fitted_values = _np.empty((n,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.loess_fit_plain_c(
        _ctypes.byref(_ctypes.c_int(n)),
        x,
        y,
        weights,
        eval_points,
        _ctypes.byref(_ctypes.c_double(span)),
        _ctypes.byref(_ctypes.c_int(degree)),
        _ctypes.byref(_ctypes.c_int(max_neighborhood_size)),
        _ctypes.byref(_ctypes.c_bool(compute_influence)),
        _ctypes.byref(_ctypes.c_bool(save_factorization)),
        fitted_values,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _LOESS_FIT_PLAIN_ARGUMENTS, _LOESS_FIT_PLAIN_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    fitted_values.flags.writeable = False

    return fitted_values

def loess_fit_robust(
        x,
        y,
        weights,
        eval_points,
        span,
        degree,
        max_neighborhood_size,
        compute_influence=False,
        save_factorization=False,
        n_iters=3,
):
    r"""Perform robust LOESS fitting with bisquare reweighting

    Data too degenerate to fit is answered directly, by the observations themselves; see
    ``loess_degenerate_fit``.
    Fits a LOESS model to the data using robust iterations to handle outliers.
    The robust fitting process iterates n_iters times, each iteration:
    - Combines original weights with robust weights (down-weights from previous iteration)
    - Runs LOESS fitting with combined weights
    - Computes residuals (y - fitted values)
    - Updates robust weights using bisquare function (suppresses large residuals)

    Parameters
    ----------
    x : np.ndarray[np.float64] of shape (n,)
        Predictor variable array
    y : np.ndarray[np.float64] of shape (n,)
        Response variable array
    weights : np.ndarray[np.float64] of shape (n,)
        Weight array for data points
    eval_points : np.ndarray[np.float64] of shape (n, 1,), column-major (order='F')
        Evaluation points (x values at which the fitted curve is computed)
    span : float
        Smoothing parameter for LOESS
        The minimum valid value is `EPS_LOESS`.
        The maximum valid value is `1.0`.
    degree : int
        Degree of the LOESS polynomial
        The minimum valid value is `0`.
        The maximum valid value is `2`.
    max_neighborhood_size : int
        Maximum neighborhood size
    compute_influence : bool, optional, default False
        Influence calculation flag
        The default value is `False`.
    save_factorization : bool, optional, default False
        Save matrix factorization flag
        The default value is `False`.
    n_iters : int, optional, default 3
        Number of robust iterations
        The minimum valid value is `1`.
        The default value is `3`.

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (n,), read-only
        Fitted (smoothed) values of y at the evaluation points
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_robust`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x = _np.ascontiguousarray(x, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'x' must be an array of np.float64: {_error}", argument="x") from None
    if x.ndim != 1:
        raise _ToxArgumentError(f"'x' must have 1 dimension, but has {x.ndim}", argument="x")
    try:
        y = _np.ascontiguousarray(y, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'y' must be an array of np.float64: {_error}", argument="y") from None
    if y.ndim != 1:
        raise _ToxArgumentError(f"'y' must have 1 dimension, but has {y.ndim}", argument="y")
    try:
        weights = _np.ascontiguousarray(weights, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'weights' must be an array of np.float64: {_error}", argument="weights") from None
    if weights.ndim != 1:
        raise _ToxArgumentError(f"'weights' must have 1 dimension, but has {weights.ndim}", argument="weights")
    try:
        eval_points = _np.asfortranarray(eval_points, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'eval_points' must be an array of np.float64: {_error}", argument="eval_points") from None
    if eval_points.ndim != 2:
        raise _ToxArgumentError(f"'eval_points' must have 2 dimensions, but has {eval_points.ndim}", argument="eval_points")

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if y.shape[0] != n:
        raise _ToxArgumentError(f"'y' has {y.shape[0]} along axis 0, but 'x' implies n == {n}; 'y' is declared with shape (n)", argument="y")
    if weights.shape[0] != n:
        raise _ToxArgumentError(f"'weights' has {weights.shape[0]} along axis 0, but 'x' implies n == {n}; 'weights' is declared with shape (n)", argument="weights")
    if eval_points.shape[0] != n:
        raise _ToxArgumentError(f"'eval_points' has {eval_points.shape[0]} along axis 0, but 'x' implies n == {n}; 'eval_points' is declared with shape (n, 1)", argument="eval_points")
    if eval_points.shape[1] != 1:
        raise _ToxArgumentError(f"'eval_points' has {eval_points.shape[1]} along axis 1, but its declared shape fixes it at {1}; 'eval_points' is declared with shape (n, 1)", argument="eval_points")

    # outputs and work arrays, which the caller never sees
    fitted_values = _np.empty((n,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.loess_fit_robust_c(
        _ctypes.byref(_ctypes.c_int(n)),
        x,
        y,
        weights,
        eval_points,
        _ctypes.byref(_ctypes.c_double(span)),
        _ctypes.byref(_ctypes.c_int(degree)),
        _ctypes.byref(_ctypes.c_int(max_neighborhood_size)),
        _ctypes.byref(_ctypes.c_bool(compute_influence)),
        _ctypes.byref(_ctypes.c_bool(save_factorization)),
        _ctypes.byref(_ctypes.c_int(n_iters)),
        fitted_values,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _LOESS_FIT_ROBUST_ARGUMENTS, _LOESS_FIT_ROBUST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    fitted_values.flags.writeable = False

    return fitted_values
