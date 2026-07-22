"""Python interface to Wraps the netlib LOESS (`dloess`/`lowesd` family) Fortran routines for local polynomial regression smoothing.

Generated from tox_loess. Do not edit.
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

_lib.loess_fit_plain_c.restype = None
_lib.loess_fit_plain_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_bool),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_FIT_PLAIN_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "int_workspace", "int_workspace_size", "real_workspace", "real_workspace_size", "hat_diag", "fitted_values", "ierr",)

_lib.loess_fit_robust_c.restype = None
_lib.loess_fit_robust_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_FIT_ROBUST_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "n_iters", "int_workspace", "int_workspace_size", "real_workspace", "real_workspace_size", "hat_diag", "robust_weights", "combined_weights", "residuals", "permutation_indices", "fitted_values", "ierr",)

_lib.loess_c.restype = None
_lib.loess_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_ARGUMENTS = ("x", "y", "span", "degree", "fitted_values", "mode", "n_iters", "ierr",)

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
    Generated from the Fortran procedure `tox_loess::tox_loess_required_workspace`.
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

def loess_fit_plain(
        x,
        y,
        weights,
        eval_points,
        span,
        degree,
        max_neighborhood_size,
        compute_influence,
        save_factorization,
        int_workspace,
        real_workspace,
        hat_diag,
):
    r"""Perform plain LOESS fitting

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
    degree : int
        Degree of the LOESS polynomial
    max_neighborhood_size : int
        Maximum neighborhood size
    compute_influence : bool
        Influence calculation flag
    save_factorization : bool
        Save matrix factorization flag
    int_workspace : np.ndarray[np.int32] of shape (int_workspace_size,), modified in place
        Integer workspace array
    real_workspace : np.ndarray[np.float64] of shape (real_workspace_size,), modified in place
        Real workspace array
    hat_diag : np.ndarray[np.float64] of shape (n,), modified in place
        Diagonal elements of the hat matrix

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (n,)
        Fitted (smoothed) values of y at the evaluation points

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_plain`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x = np.ascontiguousarray(x, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x' must be an array of np.float64: {error}") from None
    if x.ndim != 1:
        raise ValueError(f"'x' must have 1 dimension, but has {x.ndim}")
    try:
        y = np.ascontiguousarray(y, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'y' must be an array of np.float64: {error}") from None
    if y.ndim != 1:
        raise ValueError(f"'y' must have 1 dimension, but has {y.ndim}")
    try:
        weights = np.ascontiguousarray(weights, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'weights' must be an array of np.float64: {error}") from None
    if weights.ndim != 1:
        raise ValueError(f"'weights' must have 1 dimension, but has {weights.ndim}")
    try:
        eval_points = np.asfortranarray(eval_points, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'eval_points' must be an array of np.float64: {error}") from None
    if eval_points.ndim != 2:
        raise ValueError(f"'eval_points' must have 2 dimensions, but has {eval_points.ndim}")
    if not isinstance(int_workspace, np.ndarray) or int_workspace.dtype != np.int32:
        raise TypeError("'int_workspace' is modified in place, so it must already be a numpy array of {}".format(np.int32))
    if int_workspace.ndim != 1:
        raise ValueError(f"'int_workspace' must have 1 dimension, but has {int_workspace.ndim}")
    if not isinstance(real_workspace, np.ndarray) or real_workspace.dtype != np.float64:
        raise TypeError("'real_workspace' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if real_workspace.ndim != 1:
        raise ValueError(f"'real_workspace' must have 1 dimension, but has {real_workspace.ndim}")
    if not isinstance(hat_diag, np.ndarray) or hat_diag.dtype != np.float64:
        raise TypeError("'hat_diag' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if hat_diag.ndim != 1:
        raise ValueError(f"'hat_diag' must have 1 dimension, but has {hat_diag.ndim}")

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]
    int_workspace_size = int_workspace.shape[0]
    real_workspace_size = real_workspace.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if y.shape[0] != n:
        raise ValueError(f"'y' has {y.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if weights.shape[0] != n:
        raise ValueError(f"'weights' has {weights.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if eval_points.shape[0] != n:
        raise ValueError(f"'eval_points' has {eval_points.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if hat_diag.shape[0] != n:
        raise ValueError(f"'hat_diag' has {hat_diag.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )

    # outputs and work arrays, which the caller never sees
    fitted_values = np.empty((n,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.loess_fit_plain_c(
        ctypes.byref(ctypes.c_int(n)),
        x,
        y,
        weights,
        eval_points,
        ctypes.byref(ctypes.c_double(span)),
        ctypes.byref(ctypes.c_int(degree)),
        ctypes.byref(ctypes.c_int(max_neighborhood_size)),
        ctypes.byref(ctypes.c_bool(compute_influence)),
        ctypes.byref(ctypes.c_bool(save_factorization)),
        int_workspace,
        ctypes.byref(ctypes.c_int(int_workspace_size)),
        real_workspace,
        ctypes.byref(ctypes.c_int(real_workspace_size)),
        hat_diag,
        fitted_values,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_FIT_PLAIN_ARGUMENTS)

    return fitted_values

def loess_fit_robust(
        x,
        y,
        weights,
        eval_points,
        span,
        degree,
        max_neighborhood_size,
        compute_influence,
        save_factorization,
        n_iters,
        int_workspace,
        real_workspace,
        hat_diag,
        robust_weights,
        combined_weights,
        residuals,
        permutation_indices,
):
    r"""Perform robust LOESS fitting with bisquare reweighting

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
    degree : int
        Degree of the LOESS polynomial
    max_neighborhood_size : int
        Maximum neighborhood size
    compute_influence : bool
        Influence calculation flag
    save_factorization : bool
        Save matrix factorization flag
    n_iters : int
        Number of robust iterations
    int_workspace : np.ndarray[np.int32] of shape (int_workspace_size,), modified in place
        Integer workspace array
    real_workspace : np.ndarray[np.float64] of shape (real_workspace_size,), modified in place
        Real workspace array
    hat_diag : np.ndarray[np.float64] of shape (n,), modified in place
        Diagonal elements of the hat matrix
    robust_weights : np.ndarray[np.float64] of shape (n,), modified in place
        Robust bisquare weights (updated each iteration, initialized to 1.0)
    combined_weights : np.ndarray[np.float64] of shape (n,), modified in place
        Combined weights: product of user weights and robust weights (weights(i) * robust_weights(i))
    residuals : np.ndarray[np.float64] of shape (n,), modified in place
        Residuals (y - fitted_values), used to compute bisquare robust weights
    permutation_indices : np.ndarray[np.int32] of shape (n,), modified in place
        Permutation indices array (from NetLib bisquare weight computation)

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (n,)
        Fitted (smoothed) values of y at the evaluation points

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_robust`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x = np.ascontiguousarray(x, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x' must be an array of np.float64: {error}") from None
    if x.ndim != 1:
        raise ValueError(f"'x' must have 1 dimension, but has {x.ndim}")
    try:
        y = np.ascontiguousarray(y, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'y' must be an array of np.float64: {error}") from None
    if y.ndim != 1:
        raise ValueError(f"'y' must have 1 dimension, but has {y.ndim}")
    try:
        weights = np.ascontiguousarray(weights, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'weights' must be an array of np.float64: {error}") from None
    if weights.ndim != 1:
        raise ValueError(f"'weights' must have 1 dimension, but has {weights.ndim}")
    try:
        eval_points = np.asfortranarray(eval_points, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'eval_points' must be an array of np.float64: {error}") from None
    if eval_points.ndim != 2:
        raise ValueError(f"'eval_points' must have 2 dimensions, but has {eval_points.ndim}")
    if not isinstance(int_workspace, np.ndarray) or int_workspace.dtype != np.int32:
        raise TypeError("'int_workspace' is modified in place, so it must already be a numpy array of {}".format(np.int32))
    if int_workspace.ndim != 1:
        raise ValueError(f"'int_workspace' must have 1 dimension, but has {int_workspace.ndim}")
    if not isinstance(real_workspace, np.ndarray) or real_workspace.dtype != np.float64:
        raise TypeError("'real_workspace' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if real_workspace.ndim != 1:
        raise ValueError(f"'real_workspace' must have 1 dimension, but has {real_workspace.ndim}")
    if not isinstance(hat_diag, np.ndarray) or hat_diag.dtype != np.float64:
        raise TypeError("'hat_diag' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if hat_diag.ndim != 1:
        raise ValueError(f"'hat_diag' must have 1 dimension, but has {hat_diag.ndim}")
    if not isinstance(robust_weights, np.ndarray) or robust_weights.dtype != np.float64:
        raise TypeError("'robust_weights' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if robust_weights.ndim != 1:
        raise ValueError(f"'robust_weights' must have 1 dimension, but has {robust_weights.ndim}")
    if not isinstance(combined_weights, np.ndarray) or combined_weights.dtype != np.float64:
        raise TypeError("'combined_weights' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if combined_weights.ndim != 1:
        raise ValueError(f"'combined_weights' must have 1 dimension, but has {combined_weights.ndim}")
    if not isinstance(residuals, np.ndarray) or residuals.dtype != np.float64:
        raise TypeError("'residuals' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if residuals.ndim != 1:
        raise ValueError(f"'residuals' must have 1 dimension, but has {residuals.ndim}")
    if not isinstance(permutation_indices, np.ndarray) or permutation_indices.dtype != np.int32:
        raise TypeError("'permutation_indices' is modified in place, so it must already be a numpy array of {}".format(np.int32))
    if permutation_indices.ndim != 1:
        raise ValueError(f"'permutation_indices' must have 1 dimension, but has {permutation_indices.ndim}")

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]
    int_workspace_size = int_workspace.shape[0]
    real_workspace_size = real_workspace.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if y.shape[0] != n:
        raise ValueError(f"'y' has {y.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if weights.shape[0] != n:
        raise ValueError(f"'weights' has {weights.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if eval_points.shape[0] != n:
        raise ValueError(f"'eval_points' has {eval_points.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if hat_diag.shape[0] != n:
        raise ValueError(f"'hat_diag' has {hat_diag.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if robust_weights.shape[0] != n:
        raise ValueError(f"'robust_weights' has {robust_weights.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if combined_weights.shape[0] != n:
        raise ValueError(f"'combined_weights' has {combined_weights.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if residuals.shape[0] != n:
        raise ValueError(f"'residuals' has {residuals.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )
    if permutation_indices.shape[0] != n:
        raise ValueError(f"'permutation_indices' has {permutation_indices.shape[0]} along axis 0, but "
            f"'x' implies n == {n}"
        )

    # outputs and work arrays, which the caller never sees
    fitted_values = np.empty((n,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.loess_fit_robust_c(
        ctypes.byref(ctypes.c_int(n)),
        x,
        y,
        weights,
        eval_points,
        ctypes.byref(ctypes.c_double(span)),
        ctypes.byref(ctypes.c_int(degree)),
        ctypes.byref(ctypes.c_int(max_neighborhood_size)),
        ctypes.byref(ctypes.c_bool(compute_influence)),
        ctypes.byref(ctypes.c_bool(save_factorization)),
        ctypes.byref(ctypes.c_int(n_iters)),
        int_workspace,
        ctypes.byref(ctypes.c_int(int_workspace_size)),
        real_workspace,
        ctypes.byref(ctypes.c_int(real_workspace_size)),
        hat_diag,
        robust_weights,
        combined_weights,
        residuals,
        permutation_indices,
        fitted_values,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_FIT_ROBUST_ARGUMENTS)

    return fitted_values

def loess(
        x,
        y,
        span,
        degree,
        mode,
        n_iters=3,
):
    r"""Wrapper subroutine for LOESS fitting (plain or robust)

    Parameters
    ----------
    x : np.ndarray[np.float64] of shape (n_x_elements,)
        Predictor variable array
    y : np.ndarray[np.float64] of shape (n_y_elements,)
        Response variable array
    span : float
        Smoothing parameter for LOESS
    degree : int
        Degree of the LOESS polynomial
    mode : str, one of 'plain' | 'robust'
        Mode of operation

    n_iters : int, optional, default 3
        Number of robust iterations, ignored in [[tox_loess(module):MODE_PLAIN(variable)]].
        The default value is `3_int32`.

    Returns
    -------
    fitted_values : np.ndarray[np.float64] of shape (size(y),)
        Fitted (smoothed) values of y

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        x = np.ascontiguousarray(x, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'x' must be an array of np.float64: {error}") from None
    if x.ndim != 1:
        raise ValueError(f"'x' must have 1 dimension, but has {x.ndim}")
    try:
        y = np.ascontiguousarray(y, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'y' must be an array of np.float64: {error}") from None
    if y.ndim != 1:
        raise ValueError(f"'y' must have 1 dimension, but has {y.ndim}")
    mode = np.array([str(mode).lower().encode()], dtype="S6")

    # what the inputs already say, rather than asking for it again
    n_x_elements = x.shape[0]
    n_y_elements = y.shape[0]

    # outputs and work arrays, which the caller never sees
    fitted_values = np.empty((y.size,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.loess_c(
        x,
        ctypes.byref(ctypes.c_int(n_x_elements)),
        y,
        ctypes.byref(ctypes.c_int(n_y_elements)),
        ctypes.byref(ctypes.c_double(span)),
        ctypes.byref(ctypes.c_int(degree)),
        fitted_values,
        mode,
        ctypes.byref(ctypes.c_int(n_iters)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_ARGUMENTS)

    return fitted_values
