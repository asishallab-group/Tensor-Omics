"""tox_loess

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_loess. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.loess_fit_plain_expert_c.restype = None
_lib.loess_fit_plain_expert_c.argtypes = (
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
_LOESS_FIT_PLAIN_EXPERT_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "fitted_values", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOESS_FIT_PLAIN_EXPERT_ARGUMENT_SOURCES = ("x", None, None, None, None, None, None, None, None, None, None, "tmp_int_workspace", None, "tmp_real_workspace", None, None, None,)

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
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_FIT_PLAIN_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "fitted_values", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOESS_FIT_PLAIN_ARGUMENT_SOURCES = ("x", None, None, None, None, None, None, None, None, None, None, None,)

_lib.loess_fit_robust_expert_c.restype = None
_lib.loess_fit_robust_expert_c.argtypes = (
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
_LOESS_FIT_ROBUST_EXPERT_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "n_iters", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "tmp_robust_weights", "tmp_combined_weights", "tmp_residuals", "tmp_permutation_indices", "fitted_values", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOESS_FIT_ROBUST_EXPERT_ARGUMENT_SOURCES = ("x", None, None, None, None, None, None, None, None, None, None, None, "tmp_int_workspace", None, "tmp_real_workspace", None, None, None, None, None, None, None,)

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
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LOESS_FIT_ROBUST_ARGUMENTS = ("n", "x", "y", "weights", "eval_points", "span", "degree", "max_neighborhood_size", "compute_influence", "save_factorization", "n_iters", "fitted_values", "ierr",)
#: For a derived argument, the one the caller passed it in
_LOESS_FIT_ROBUST_ARGUMENT_SOURCES = ("x", None, None, None, None, None, None, None, None, None, None, None, None,)

def loess_fit_plain_expert(
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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_plain`, whose argument names are
    the ones an error message reports.
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

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_loess_kernel import tox_loess_required_workspace
    _tox_loess_required_workspace_result = tox_loess_required_workspace(n_dim=1, max_neighborhood_size=max_neighborhood_size, save_factorization=save_factorization)
    int_workspace_size = _tox_loess_required_workspace_result["int_workspace_size"]
    from .tox_loess_kernel import tox_loess_required_workspace
    real_workspace_size = _tox_loess_required_workspace_result["real_workspace_size"]

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

    # outputs and work arrays, which the caller never sees
    tmp_int_workspace = np.empty((int_workspace_size,), dtype=np.int32, order='C')
    tmp_real_workspace = np.empty((real_workspace_size,), dtype=np.float64, order='C')
    tmp_hat_diag = np.empty((n,), dtype=np.float64, order='C')
    fitted_values = np.empty((n,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.loess_fit_plain_expert_c(
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
        tmp_int_workspace,
        ctypes.byref(ctypes.c_int(int_workspace_size)),
        tmp_real_workspace,
        ctypes.byref(ctypes.c_int(real_workspace_size)),
        tmp_hat_diag,
        fitted_values,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_FIT_PLAIN_EXPERT_ARGUMENTS, _LOESS_FIT_PLAIN_EXPERT_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    fitted_values.flags.writeable = False

    return fitted_values

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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_plain_alloc`, whose argument names are
    the ones an error message reports.
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

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]

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
        fitted_values,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_FIT_PLAIN_ARGUMENTS, _LOESS_FIT_PLAIN_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    fitted_values.flags.writeable = False

    return fitted_values

def loess_fit_robust_expert(
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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_robust`, whose argument names are
    the ones an error message reports.
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

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]

    # work out what other procedures must supply, per DM_OUTPUT_FROM
    from .tox_loess_kernel import tox_loess_required_workspace
    _tox_loess_required_workspace_result = tox_loess_required_workspace(n_dim=1, max_neighborhood_size=max_neighborhood_size, save_factorization=save_factorization)
    int_workspace_size = _tox_loess_required_workspace_result["int_workspace_size"]
    from .tox_loess_kernel import tox_loess_required_workspace
    real_workspace_size = _tox_loess_required_workspace_result["real_workspace_size"]

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

    # outputs and work arrays, which the caller never sees
    tmp_int_workspace = np.empty((int_workspace_size,), dtype=np.int32, order='C')
    tmp_real_workspace = np.empty((real_workspace_size,), dtype=np.float64, order='C')
    tmp_hat_diag = np.empty((n,), dtype=np.float64, order='C')
    tmp_robust_weights = np.empty((n,), dtype=np.float64, order='C')
    tmp_combined_weights = np.empty((n,), dtype=np.float64, order='C')
    tmp_residuals = np.empty((n,), dtype=np.float64, order='C')
    tmp_permutation_indices = np.empty((n,), dtype=np.int32, order='C')
    fitted_values = np.empty((n,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.loess_fit_robust_expert_c(
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
        tmp_int_workspace,
        ctypes.byref(ctypes.c_int(int_workspace_size)),
        tmp_real_workspace,
        ctypes.byref(ctypes.c_int(real_workspace_size)),
        tmp_hat_diag,
        tmp_robust_weights,
        tmp_combined_weights,
        tmp_residuals,
        tmp_permutation_indices,
        fitted_values,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_FIT_ROBUST_EXPERT_ARGUMENTS, _LOESS_FIT_ROBUST_EXPERT_ARGUMENT_SOURCES)

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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_loess::loess_fit_robust_alloc`, whose argument names are
    the ones an error message reports.
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

    # what the inputs already say, rather than asking for it again
    n = x.shape[0]

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
        fitted_values,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LOESS_FIT_ROBUST_ARGUMENTS, _LOESS_FIT_ROBUST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    fitted_values.flags.writeable = False

    return fitted_values
