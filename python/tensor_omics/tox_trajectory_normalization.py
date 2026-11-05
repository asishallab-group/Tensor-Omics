r"""tox_trajectory_normalization

Min-max normalization of factor trajectories over time.

Each factor's time series is rescaled to `[0,1]` independently, per sample -- so trajectories
of very different magnitudes can be compared by shape. `normalize_single_trajectory` does one;
`normalize_all_trajectories` does a whole tensor of them in one call.

Python binding, generated from tox_trajectory_normalization. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.normalize_variable_timeseries_c.restype = None
_lib.normalize_variable_timeseries_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_VARIABLE_TIMESERIES_ARGUMENTS = ("v", "v_norm", "n_points", "status", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZE_VARIABLE_TIMESERIES_ARGUMENT_SOURCES = (None, None, "v", None, None,)

_lib.normalize_single_trajectory_c.restype = None
_lib.normalize_single_trajectory_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_SINGLE_TRAJECTORY_ARGUMENTS = ("trajectory", "trajectory_norm", "n_factors", "n_timepoints", "status", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZE_SINGLE_TRAJECTORY_ARGUMENT_SOURCES = (None, None, "trajectory", "trajectory", None, None,)

_lib.normalize_all_trajectories_c.restype = None
_lib.normalize_all_trajectories_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_NORMALIZE_ALL_TRAJECTORIES_ARGUMENTS = ("trajectories", "trajectories_norm", "n_factors", "n_samples", "n_timepoints", "status", "ierr",)
#: For a derived argument, the one the caller passed it in
_NORMALIZE_ALL_TRAJECTORIES_ARGUMENT_SOURCES = (None, None, "trajectories", "trajectories", "trajectories", None, None,)

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

        v_norm : np.ndarray[np.float64] of shape (n_points,), read-only
            Normalized time series
            A result is a value; call `.copy()` to obtain a modifiable array.
        status : int
            Status code for specific warnings

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_normalization::normalize_variable_timeseries`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    v = _contract.reals(v, "v", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_points = v.shape[0]

    # outputs and work arrays, which the caller never sees
    v_norm = _np.empty((n_points,), dtype=_np.float64, order='C')
    status = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.normalize_variable_timeseries_c(
        v,
        v_norm,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(status),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _NORMALIZE_VARIABLE_TIMESERIES_ARGUMENTS, _NORMALIZE_VARIABLE_TIMESERIES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    v_norm.flags.writeable = False

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

        trajectory_norm : np.ndarray[np.float64] of shape (n_timepoints, n_factors,), column-major (order='F'), read-only
            Normalized trajectory for one sample
            A result is a value; call `.copy()` to obtain a modifiable array.
        status : np.ndarray[np.int32] of shape (n_factors,), read-only
            Status code for specific warnings, one per factor
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
    Generated from the Fortran procedure `tox_trajectory_normalization::normalize_single_trajectory`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectory = _contract.reals(trajectory, "trajectory", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectory.shape[1]
    n_timepoints = trajectory.shape[0]

    # outputs and work arrays, which the caller never sees
    trajectory_norm = _np.empty((n_timepoints, n_factors,), dtype=_np.float64, order='F')
    status = _np.empty((n_factors,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.normalize_single_trajectory_c(
        trajectory,
        trajectory_norm,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        status,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _NORMALIZE_SINGLE_TRAJECTORY_ARGUMENTS, _NORMALIZE_SINGLE_TRAJECTORY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    trajectory_norm.flags.writeable = False
    status.flags.writeable = False

    return {
        "trajectory_norm": trajectory_norm,
        "status": status,
    }

def normalize_all_trajectories(
        trajectories,
):
    r"""Normalize all trajectories across multiple entities

    independently across time for each sample.

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        Original trajectories

    Returns
    -------
    dict
        with keys:

        trajectories_norm : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F'), read-only
            Normalized trajectories
            A result is a value; call `.copy()` to obtain a modifiable array.
        status : np.ndarray[np.int32] of shape (n_factors, n_samples,), column-major (order='F'), read-only
            Status code for specific warnings, one per factor per sample
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
    Generated from the Fortran procedure `tox_trajectory_normalization::normalize_all_trajectories`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectories = _contract.reals(trajectories, "trajectories", ndim=3, order="F")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    trajectories_norm = _np.empty((n_factors, n_samples, n_timepoints,), dtype=_np.float64, order='F')
    status = _np.empty((n_factors, n_samples,), dtype=_np.int32, order='F')
    ierr = _ctypes.c_int(0)

    _lib.normalize_all_trajectories_c(
        trajectories,
        trajectories_norm,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        status,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _NORMALIZE_ALL_TRAJECTORIES_ARGUMENTS, _NORMALIZE_ALL_TRAJECTORIES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    trajectories_norm.flags.writeable = False
    status.flags.writeable = False

    return {
        "trajectories_norm": trajectories_norm,
        "status": status,
    }
