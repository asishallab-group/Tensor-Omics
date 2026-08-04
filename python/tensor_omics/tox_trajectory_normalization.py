"""Python binding to Generated from the kernel; do not edit -- regenerate instead.

Generated from tox_trajectory_normalization. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.normalize_variable_timeseries_c.restype = None
_lib.normalize_variable_timeseries_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_VARIABLE_TIMESERIES_ARGUMENTS = ("v", "v_norm", "n_points", "status", "ierr",)

_lib.normalize_single_trajectory_c.restype = None
_lib.normalize_single_trajectory_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_SINGLE_TRAJECTORY_ARGUMENTS = ("trajectory", "trajectory_norm", "n_factors", "n_timepoints", "status", "ierr",)

_lib.normalize_all_trajectories_expert_c.restype = None
_lib.normalize_all_trajectories_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_ALL_TRAJECTORIES_EXPERT_ARGUMENTS = ("trajectories", "trajectories_norm", "n_factors", "n_samples", "n_timepoints", "tmp_series", "tmp_series_norm", "status", "ierr",)

_lib.normalize_all_trajectories_c.restype = None
_lib.normalize_all_trajectories_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_ALL_TRAJECTORIES_ARGUMENTS = ("trajectories", "trajectories_norm", "n_factors", "n_samples", "n_timepoints", "status", "ierr",)

def normalize_variable_timeseries(
        v,
):
    r"""Normalize a single variable across time using min-max scaling

    Parameters
    ----------
    v : np.ndarray[np.float64] of shape (n_points,)
        Original time series

    Returns
    -------
    dict
        with keys:

        v_norm : np.ndarray[np.float64] of shape (n_points,)
            Normalized time series
        status : int
            Status code for specific warnings

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_normalization::normalize_variable_timeseries`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        v = np.ascontiguousarray(v, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'v' must be an array of np.float64: {error}") from None
    if v.ndim != 1:
        raise ValueError(f"'v' must have 1 dimension, but has {v.ndim}")

    # what the inputs already say, rather than asking for it again
    n_points = v.shape[0]

    # outputs and work arrays, which the caller never sees
    v_norm = np.empty((n_points,), dtype=np.float64, order='C')
    status = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.normalize_variable_timeseries_c(
        v,
        v_norm,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(status),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMALIZE_VARIABLE_TIMESERIES_ARGUMENTS)

    return {
        "v_norm": v_norm,
        "status": status.value,
    }

def normalize_single_trajectory(
        trajectory,
):
    r"""Normalize all factors in a single trajectory independently across time

    Parameters
    ----------
    trajectory : np.ndarray[np.float64] of shape (n_timepoints, n_factors,), column-major (order='F')
        Original trajectory for one sample

    Returns
    -------
    dict
        with keys:

        trajectory_norm : np.ndarray[np.float64] of shape (n_timepoints, n_factors,), column-major (order='F')
            Normalized trajectory for one sample
        status : np.ndarray[np.int32] of shape (n_factors,)
            Status code for specific warnings, one per factor

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_normalization::normalize_single_trajectory`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectory = np.asfortranarray(trajectory, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectory' must be an array of np.float64: {error}") from None
    if trajectory.ndim != 2:
        raise ValueError(f"'trajectory' must have 2 dimensions, but has {trajectory.ndim}")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectory.shape[1]
    n_timepoints = trajectory.shape[0]

    # outputs and work arrays, which the caller never sees
    trajectory_norm = np.empty((n_timepoints, n_factors,), dtype=np.float64, order='F')
    status = np.empty((n_factors,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.normalize_single_trajectory_c(
        trajectory,
        trajectory_norm,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        status,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMALIZE_SINGLE_TRAJECTORY_ARGUMENTS)

    return {
        "trajectory_norm": trajectory_norm,
        "status": status,
    }

def normalize_all_trajectories_expert(
        trajectories,
):
    r"""Normalize all trajectories across multiple entities

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        Original trajectories

    Returns
    -------
    dict
        with keys:

        trajectories_norm : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
            Normalized trajectories
        status : np.ndarray[np.int32] of shape (n_factors, n_samples,), column-major (order='F')
            Status code for specific warnings, one per factor per sample

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_normalization::normalize_all_trajectories`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectories' must be an array of np.float64: {error}") from None
    if trajectories.ndim != 3:
        raise ValueError(f"'trajectories' must have 3 dimensions, but has {trajectories.ndim}")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    trajectories_norm = np.empty((n_factors, n_samples, n_timepoints,), dtype=np.float64, order='F')
    tmp_series = np.empty((n_timepoints,), dtype=np.float64, order='C')
    tmp_series_norm = np.empty((n_timepoints,), dtype=np.float64, order='C')
    status = np.empty((n_factors, n_samples,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.normalize_all_trajectories_expert_c(
        trajectories,
        trajectories_norm,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        tmp_series,
        tmp_series_norm,
        status,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMALIZE_ALL_TRAJECTORIES_EXPERT_ARGUMENTS)

    return {
        "trajectories_norm": trajectories_norm,
        "status": status,
    }

def normalize_all_trajectories(
        trajectories,
):
    r"""Normalize all trajectories across multiple entities

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        Original trajectories

    Returns
    -------
    dict
        with keys:

        trajectories_norm : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
            Normalized trajectories
        status : np.ndarray[np.int32] of shape (n_factors, n_samples,), column-major (order='F')
            Status code for specific warnings, one per factor per sample

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_normalization::normalize_all_trajectories_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectories' must be an array of np.float64: {error}") from None
    if trajectories.ndim != 3:
        raise ValueError(f"'trajectories' must have 3 dimensions, but has {trajectories.ndim}")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    trajectories_norm = np.empty((n_factors, n_samples, n_timepoints,), dtype=np.float64, order='F')
    status = np.empty((n_factors, n_samples,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    _lib.normalize_all_trajectories_c(
        trajectories,
        trajectories_norm,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        status,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _NORMALIZE_ALL_TRAJECTORIES_ARGUMENTS)

    return {
        "trajectories_norm": trajectories_norm,
        "status": status,
    }
