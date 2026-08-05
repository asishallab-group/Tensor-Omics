"""tox_trajectory_contribution_analysis_kernel

Module for quantifying how much one trajectory (a "factor") contributes to another (a "dependent") over time.

Python binding, generated from tox_trajectory_contribution_analysis_kernel. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.compute_all_contributions_kernel_c.restype = None
_lib.compute_all_contributions_kernel_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ALL_CONTRIBUTIONS_KERNEL_ARGUMENTS = ("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_indices", "n_selected_factors", "dependent_indices", "n_selected_dependents", "baseline_mode", "local_contributions", "total_contributions", "tmp_factors", "tmp_dependent", "ierr",)

_lib.compute_baselines_factor_dependent_kernel_c.restype = None
_lib.compute_baselines_factor_dependent_kernel_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_BASELINES_FACTOR_DEPENDENT_KERNEL_ARGUMENTS = ("n_timepoints", "factor", "dependent", "baseline_mode", "factor_baseline", "dependent_baseline", "ierr",)

_lib.compute_velocity_trajectory_kernel_c.restype = None
_lib.compute_velocity_trajectory_kernel_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_VELOCITY_TRAJECTORY_KERNEL_ARGUMENTS = ("trajectory", "velocity", "n_timepoints",)

_lib.compute_acceleration_from_velocity_trajectory_kernel_c.restype = None
_lib.compute_acceleration_from_velocity_trajectory_kernel_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_KERNEL_ARGUMENTS = ("velocity", "acceleration", "n_timepoints",)

_lib.compute_velocity_acceleration_contributions_kernel_c.restype = None
_lib.compute_velocity_acceleration_contributions_kernel_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_KERNEL_ARGUMENTS = ("trajectories", "n_factors", "n_samples", "n_timepoints", "baseline_mode", "tmp_factors", "tmp_dependent", "tmp_contributions", "contrib_velocity", "velocity_contribution_series", "contrib_acceleration", "acceleration_contribution_series", "ierr",)

def compute_all_contributions_kernel(
        trajectories,
        factor_indices,
        dependent_indices,
        baseline_mode,
):
    r"""Contribution analysis for every selected factor-dependent pair

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        expression vectors across different samples over time
    factor_indices : np.ndarray[np.int32] of shape (n_selected_factors,)
        indices of factors to compute the contributions for
        The minimum valid value is `1`.
        The maximum valid value is `n_factors`.
    dependent_indices : np.ndarray[np.int32] of shape (n_selected_dependents,)
        indices of dependents to compute the contributions for
        The minimum valid value is `1`.
        The maximum valid value is `n_factors`.
    baseline_mode : str, one of 'raw' | 'mean' | 'min'

    Returns
    -------
    dict
        with keys:

        local_contributions : np.ndarray[np.float64] of shape (n_timepoints, n_selected_factors, n_selected_dependents, n_samples,), column-major (order='F')
            Per-timepoint contributions per sample-dependent-factor combination
        total_contributions : np.ndarray[np.float64] of shape (n_selected_factors, n_selected_dependents, n_samples,), column-major (order='F')
            Total contribution (`sum(local_contributions)`) per sample-dependent-factor combination

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis_kernel::compute_all_contributions_kernel`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectories' must be an array of np.float64: {error}") from None
    if trajectories.ndim != 3:
        raise ValueError(f"'trajectories' must have 3 dimensions, but has {trajectories.ndim}")
    try:
        factor_indices = np.ascontiguousarray(factor_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'factor_indices' must be an array of np.int32: {error}") from None
    if factor_indices.ndim != 1:
        raise ValueError(f"'factor_indices' must have 1 dimension, but has {factor_indices.ndim}")
    try:
        dependent_indices = np.ascontiguousarray(dependent_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dependent_indices' must be an array of np.int32: {error}") from None
    if dependent_indices.ndim != 1:
        raise ValueError(f"'dependent_indices' must have 1 dimension, but has {dependent_indices.ndim}")
    baseline_mode = np.array([str(baseline_mode).lower().encode()], dtype="S4")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]
    n_selected_factors = factor_indices.shape[0]
    n_selected_dependents = dependent_indices.shape[0]

    # outputs and work arrays, which the caller never sees
    local_contributions = np.empty((n_timepoints, n_selected_factors, n_selected_dependents, n_samples,), dtype=np.float64, order='F')
    total_contributions = np.empty((n_selected_factors, n_selected_dependents, n_samples,), dtype=np.float64, order='F')
    tmp_factors = np.empty((n_timepoints, n_selected_factors,), dtype=np.float64, order='F')
    tmp_dependent = np.empty((n_timepoints,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_all_contributions_kernel_c(
        trajectories,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        factor_indices,
        ctypes.byref(ctypes.c_int(n_selected_factors)),
        dependent_indices,
        ctypes.byref(ctypes.c_int(n_selected_dependents)),
        baseline_mode,
        local_contributions,
        total_contributions,
        tmp_factors,
        tmp_dependent,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ALL_CONTRIBUTIONS_KERNEL_ARGUMENTS)

    return {
        "local_contributions": local_contributions,
        "total_contributions": total_contributions,
    }

def compute_baselines_factor_dependent_kernel(
        factor,
        dependent,
        baseline_mode,
):
    r"""Compute scalar baselines for a factor and dependent variable time series

    Parameters
    ----------
    factor : np.ndarray[np.float64] of shape (n_timepoints,)
        Factor time series, length n_timepoints
    dependent : np.ndarray[np.float64] of shape (n_timepoints,)
        Dependent variable time series, length n_timepoints
    baseline_mode : str, one of 'raw' | 'mean' | 'min'

    Returns
    -------
    dict
        with keys:

        factor_baseline : float
            Computed baseline for factor
        dependent_baseline : float
            Computed baseline for dependent variable

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis_kernel::compute_baselines_factor_dependent_kernel`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        factor = np.ascontiguousarray(factor, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'factor' must be an array of np.float64: {error}") from None
    if factor.ndim != 1:
        raise ValueError(f"'factor' must have 1 dimension, but has {factor.ndim}")
    try:
        dependent = np.ascontiguousarray(dependent, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'dependent' must be an array of np.float64: {error}") from None
    if dependent.ndim != 1:
        raise ValueError(f"'dependent' must have 1 dimension, but has {dependent.ndim}")
    baseline_mode = np.array([str(baseline_mode).lower().encode()], dtype="S4")

    # what the inputs already say, rather than asking for it again
    n_timepoints = factor.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if dependent.shape[0] != n_timepoints:
        raise ValueError(f"'dependent' has {dependent.shape[0]} along axis 0, but "
            f"'factor' implies n_timepoints == {n_timepoints}"
        )

    # outputs and work arrays, which the caller never sees
    factor_baseline = ctypes.c_double(0)
    dependent_baseline = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.compute_baselines_factor_dependent_kernel_c(
        ctypes.byref(ctypes.c_int(n_timepoints)),
        factor,
        dependent,
        baseline_mode,
        ctypes.byref(factor_baseline),
        ctypes.byref(dependent_baseline),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_BASELINES_FACTOR_DEPENDENT_KERNEL_ARGUMENTS)

    return {
        "factor_baseline": factor_baseline.value,
        "dependent_baseline": dependent_baseline.value,
    }

def compute_velocity_trajectory_kernel(
        trajectory,
):
    r"""Compute velocity trajectory from a single position trajectory

    Parameters
    ----------
    trajectory : np.ndarray[np.float64] of shape (n_timepoints,)
        input position trajectory

    Returns
    -------
    velocity : np.ndarray[np.float64] of shape (max(0, n_timepoints - 1),)
        output velocity trajectory

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis_kernel::compute_velocity_trajectory_kernel`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectory = np.ascontiguousarray(trajectory, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectory' must be an array of np.float64: {error}") from None
    if trajectory.ndim != 1:
        raise ValueError(f"'trajectory' must have 1 dimension, but has {trajectory.ndim}")

    # what the inputs already say, rather than asking for it again
    n_timepoints = trajectory.shape[0]

    # outputs and work arrays, which the caller never sees
    velocity = np.empty((max(0, n_timepoints - 1),), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_velocity_trajectory_kernel_c(
        trajectory,
        velocity,
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_VELOCITY_TRAJECTORY_KERNEL_ARGUMENTS)

    return velocity

def compute_acceleration_from_velocity_trajectory_kernel(
        velocity,
        n_timepoints,
):
    r"""Compute acceleration trajectory from a single velocity trajectory

    Parameters
    ----------
    velocity : np.ndarray[np.float64] of shape (max(0, n_timepoints - 1),)
        velocity trajectory
    n_timepoints : int
        number of timepoints

    Returns
    -------
    acceleration : np.ndarray[np.float64] of shape (max(0, n_timepoints - 2),)
        acceleration trajectory

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis_kernel::compute_acceleration_from_velocity_trajectory_kernel`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        velocity = np.ascontiguousarray(velocity, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'velocity' must be an array of np.float64: {error}") from None
    if velocity.ndim != 1:
        raise ValueError(f"'velocity' must have 1 dimension, but has {velocity.ndim}")

    # outputs and work arrays, which the caller never sees
    acceleration = np.empty((max(0, n_timepoints - 2),), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_acceleration_from_velocity_trajectory_kernel_c(
        velocity,
        acceleration,
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_KERNEL_ARGUMENTS)

    return acceleration

def compute_velocity_acceleration_contributions_kernel(
        trajectories,
        baseline_mode,
):
    r"""Compute velocity and acceleration contributions for all variable pairs

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        input position trajectories
    baseline_mode : str, one of 'raw' | 'mean' | 'min'

    Returns
    -------
    dict
        with keys:

        contrib_velocity : np.ndarray[np.float64] of shape (n_factors, n_factors, n_samples,), column-major (order='F')
            output velocity contributions
        velocity_contribution_series : np.ndarray[np.float64] of shape (n_timepoints, n_factors, n_factors, n_samples,), column-major (order='F')
            output velocity contribution series
        contrib_acceleration : np.ndarray[np.float64] of shape (n_factors, n_factors, n_samples,), column-major (order='F')
            output acceleration contributions
        acceleration_contribution_series : np.ndarray[np.float64] of shape (n_timepoints, n_factors, n_factors, n_samples,), column-major (order='F')
            output acceleration contribution series

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis_kernel::compute_velocity_acceleration_contributions_kernel`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectories' must be an array of np.float64: {error}") from None
    if trajectories.ndim != 3:
        raise ValueError(f"'trajectories' must have 3 dimensions, but has {trajectories.ndim}")
    baseline_mode = np.array([str(baseline_mode).lower().encode()], dtype="S4")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    tmp_factors = np.empty((n_timepoints - 1, n_factors,), dtype=np.float64, order='F')
    tmp_dependent = np.empty((n_timepoints - 1,), dtype=np.float64, order='C')
    tmp_contributions = np.empty((n_timepoints - 1,), dtype=np.float64, order='C')
    contrib_velocity = np.empty((n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    velocity_contribution_series = np.empty((n_timepoints, n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    contrib_acceleration = np.empty((n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    acceleration_contribution_series = np.empty((n_timepoints, n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.compute_velocity_acceleration_contributions_kernel_c(
        trajectories,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        baseline_mode,
        tmp_factors,
        tmp_dependent,
        tmp_contributions,
        contrib_velocity,
        velocity_contribution_series,
        contrib_acceleration,
        acceleration_contribution_series,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_KERNEL_ARGUMENTS)

    return {
        "contrib_velocity": contrib_velocity,
        "velocity_contribution_series": velocity_contribution_series,
        "contrib_acceleration": contrib_acceleration,
        "acceleration_contribution_series": acceleration_contribution_series,
    }
