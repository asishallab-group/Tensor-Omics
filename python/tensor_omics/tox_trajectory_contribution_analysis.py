r"""tox_trajectory_contribution_analysis

Module for quantifying how much one trajectory (a "factor") contributes to another (a "dependent") over time.

Contributions are computed per timepoint as the product of both series' deviations from a chosen
baseline, for raw expression trajectories as well as for their velocity (first difference) and
acceleration (second difference) derivatives. Statistical significance of an observed contribution can
be assessed via a permutation test that recomputes the same contribution against a randomly chosen
other sample.

Python binding, generated from tox_trajectory_contribution_analysis. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.perform_permutation_test_c.restype = None
_lib.perform_permutation_test_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    nullable(ctypes.POINTER(ctypes.c_int)),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_PERFORM_PERMUTATION_TEST_ARGUMENTS = ("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_idx", "dependent_idx", "sample_idx", "baseline_mode", "n_permutations", "local_contributions", "total_contributions", "random_seed", "ierr",)
#: For a derived argument, the one the caller passed it in
_PERFORM_PERMUTATION_TEST_ARGUMENT_SOURCES = (None, "trajectories", "trajectories", "trajectories", None, None, None, None, "local_contributions", None, None, None, None,)

_lib.compute_p_values_c.restype = None
_lib.compute_p_values_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_P_VALUES_ARGUMENTS = ("local_contributions_observed", "total_contribution_observed", "local_contributions_perm_test", "total_contributions_perm_test", "n_timepoints", "n_permutations", "local_p_values", "total_p_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_P_VALUES_ARGUMENT_SOURCES = (None, None, None, None, "local_contributions_observed", "local_contributions_perm_test", None, None, None,)

_lib.compute_contributions_c.restype = None
_lib.compute_contributions_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_CONTRIBUTIONS_ARGUMENTS = ("factor", "dependent", "n_dims", "baseline_mode", "local_contributions", "total_contribution", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_CONTRIBUTIONS_ARGUMENT_SOURCES = (None, None, "factor", None, None, None, None,)

_lib.compute_all_contributions_c.restype = None
_lib.compute_all_contributions_c.argtypes = (
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
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ALL_CONTRIBUTIONS_ARGUMENTS = ("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_indices", "n_selected_factors", "dependent_indices", "n_selected_dependents", "baseline_mode", "local_contributions", "total_contributions", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ALL_CONTRIBUTIONS_ARGUMENT_SOURCES = (None, "trajectories", "trajectories", "trajectories", None, "factor_indices", None, "dependent_indices", None, None, None, None,)

_lib.compute_baselines_factor_dependent_c.restype = None
_lib.compute_baselines_factor_dependent_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENTS = ("n_timepoints", "factor", "dependent", "baseline_mode", "factor_baseline", "dependent_baseline", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENT_SOURCES = ("factor", None, None, None, None, None, None,)

_lib.compute_velocity_trajectory_c.restype = None
_lib.compute_velocity_trajectory_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_VELOCITY_TRAJECTORY_ARGUMENTS = ("trajectory", "velocity", "n_timepoints", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_VELOCITY_TRAJECTORY_ARGUMENT_SOURCES = (None, None, "trajectory", None,)

_lib.compute_acceleration_from_velocity_trajectory_c.restype = None
_lib.compute_acceleration_from_velocity_trajectory_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_ARGUMENTS = ("velocity", "acceleration", "n_timepoints", "ierr",)

_lib.compute_velocity_trajectories_c.restype = None
_lib.compute_velocity_trajectories_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENTS = ("trajectories", "velocity", "n_factors", "n_samples", "n_timepoints", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENT_SOURCES = (None, None, "trajectories", "trajectories", "trajectories", None,)

_lib.compute_acceleration_from_velocity_c.restype = None
_lib.compute_acceleration_from_velocity_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENTS = ("velocity", "acceleration", "n_factors", "n_samples", "n_timepoints", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENT_SOURCES = (None, None, "velocity", "velocity", None, None,)

_lib.compute_velocity_acceleration_contributions_c.restype = None
_lib.compute_velocity_acceleration_contributions_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=4, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_ARGUMENTS = ("trajectories", "n_factors", "n_samples", "n_timepoints", "baseline_mode", "contrib_velocity", "velocity_contribution_series", "contrib_acceleration", "acceleration_contribution_series", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_ARGUMENT_SOURCES = (None, "trajectories", "trajectories", "trajectories", None, None, None, None, None, None,)

def perform_permutation_test(
        trajectories,
        factor_idx,
        dependent_idx,
        sample_idx,
        baseline_mode,
        n_permutations,
        random_seed=None,
):
    r"""For a factor-dependent pair, calculates the contributions against the same dependent taken from a random different sample

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        expression vectors across different samples over time
    factor_idx : int
        index of factor to compute the permutation contributions for
        The minimum valid value is `1`.
        The maximum valid value is `n_factors`.
    dependent_idx : int
        index of dependent to compute the permutation contributions for
        The minimum valid value is `1`.
        The maximum valid value is `n_factors`.
    sample_idx : int
        index of sample to compute the permutation contributions for
        The minimum valid value is `1`.
        The maximum valid value is `n_samples`.
    baseline_mode : str, one of 'raw' | 'mean' | 'min'
    n_permutations : int
        number of permutations to perform
    random_seed : int, optional
        Seed to use for random number generation.

    Returns
    -------
    dict
        with keys:

        local_contributions : np.ndarray[np.float64] of shape (n_timepoints, n_permutations,), column-major (order='F'), read-only
            Per-timepoint contributions per permutation
            A result is a value; call `.copy()` to obtain a modifiable array.
        total_contributions : np.ndarray[np.float64] of shape (n_permutations,), read-only
            Total contribution (`sum(local_contributions)`) per permutation
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::perform_permutation_test`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectories' must be an array of np.float64: {error}") from None
    if trajectories.ndim != 3:
        raise ValueError(f"'trajectories' must have 3 dimensions, but has {trajectories.ndim}")
    baseline_mode = np.array([str(baseline_mode).lower().encode().ljust(4)], dtype="S4")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    local_contributions = np.empty((n_timepoints, n_permutations,), dtype=np.float64, order='F')
    total_contributions = np.empty((n_permutations,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.perform_permutation_test_c(
        trajectories,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ctypes.c_int(factor_idx)),
        ctypes.byref(ctypes.c_int(dependent_idx)),
        ctypes.byref(ctypes.c_int(sample_idx)),
        baseline_mode,
        ctypes.byref(ctypes.c_int(n_permutations)),
        local_contributions,
        total_contributions,
        None if random_seed is None else ctypes.byref(ctypes.c_int(random_seed)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _PERFORM_PERMUTATION_TEST_ARGUMENTS, _PERFORM_PERMUTATION_TEST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    local_contributions.flags.writeable = False
    total_contributions.flags.writeable = False

    return {
        "local_contributions": local_contributions,
        "total_contributions": total_contributions,
    }

def compute_p_values(
        local_contributions_observed,
        total_contribution_observed,
        local_contributions_perm_test,
        total_contributions_perm_test,
):
    r"""Calculates the p values for the contributions once the permutation tests are done

    Parameters
    ----------
    local_contributions_observed : np.ndarray[np.float64] of shape (n_timepoints,)
        Per-timepoint contributions for the observed factor-dependent-sample combination
    total_contribution_observed : float
        Total contribution (`sum(local_contributions)`) for the observed factor-dependent-sample combination
    local_contributions_perm_test : np.ndarray[np.float64] of shape (n_timepoints, n_permutations,), column-major (order='F')
        Per-timepoint contributions for the factor-dependent-random_sample combinations from :func:`tensor_omics.perform_permutation_test`
    total_contributions_perm_test : np.ndarray[np.float64] of shape (n_permutations,)
        Total contribution (`sum(local_contributions)`) for the factor-dependent-random_sample combinations from :func:`tensor_omics.perform_permutation_test`

    Returns
    -------
    dict
        with keys:

        local_p_values : np.ndarray[np.float64] of shape (n_timepoints,), read-only
            calculated p values for local contributions, like: `(local_contributions_perm_test >= local_contributions_observed)/n_permutations`
            A result is a value; call `.copy()` to obtain a modifiable array.
        total_p_value : float
            calculated p values for total contributions, like: `(total_contributions_perm_test >= total_contribution_observed)/n_permutations`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_p_values`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        local_contributions_observed = np.ascontiguousarray(local_contributions_observed, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'local_contributions_observed' must be an array of np.float64: {error}") from None
    if local_contributions_observed.ndim != 1:
        raise ValueError(f"'local_contributions_observed' must have 1 dimension, but has {local_contributions_observed.ndim}")
    try:
        local_contributions_perm_test = np.asfortranarray(local_contributions_perm_test, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'local_contributions_perm_test' must be an array of np.float64: {error}") from None
    if local_contributions_perm_test.ndim != 2:
        raise ValueError(f"'local_contributions_perm_test' must have 2 dimensions, but has {local_contributions_perm_test.ndim}")
    try:
        total_contributions_perm_test = np.ascontiguousarray(total_contributions_perm_test, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'total_contributions_perm_test' must be an array of np.float64: {error}") from None
    if total_contributions_perm_test.ndim != 1:
        raise ValueError(f"'total_contributions_perm_test' must have 1 dimension, but has {total_contributions_perm_test.ndim}")

    # what the inputs already say, rather than asking for it again
    n_timepoints = local_contributions_observed.shape[0]
    n_permutations = local_contributions_perm_test.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if local_contributions_perm_test.shape[0] != n_timepoints:
        raise ValueError(f"'local_contributions_perm_test' has {local_contributions_perm_test.shape[0]} along axis 0, but "
            f"'local_contributions_observed' implies n_timepoints == {n_timepoints}"
        )
    if total_contributions_perm_test.shape[0] != n_permutations:
        raise ValueError(f"'total_contributions_perm_test' has {total_contributions_perm_test.shape[0]} along axis 0, but "
            f"'local_contributions_perm_test' implies n_permutations == {n_permutations}"
        )

    # outputs and work arrays, which the caller never sees
    local_p_values = np.empty((n_timepoints,), dtype=np.float64, order='C')
    total_p_value = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.compute_p_values_c(
        local_contributions_observed,
        ctypes.byref(ctypes.c_double(total_contribution_observed)),
        local_contributions_perm_test,
        total_contributions_perm_test,
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ctypes.c_int(n_permutations)),
        local_p_values,
        ctypes.byref(total_p_value),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_P_VALUES_ARGUMENTS, _COMPUTE_P_VALUES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    local_p_values.flags.writeable = False

    return {
        "local_p_values": local_p_values,
        "total_p_value": total_p_value.value,
    }

def compute_contributions(
        factor,
        dependent,
        baseline_mode,
):
    r"""Performs contribution analysis for a specific factor-dependent pair

    Parameters
    ----------
    factor : np.ndarray[np.float64] of shape (n_dims,)
        Factor time series, length n_timepoints
    dependent : np.ndarray[np.float64] of shape (n_dims,)
        Dependent variable time series, length n_timepoints
    baseline_mode : str, one of 'raw' | 'mean' | 'min'

    Returns
    -------
    dict
        with keys:

        local_contributions : np.ndarray[np.float64] of shape (n_dims,), read-only
            Per-element contributions
            A result is a value; call `.copy()` to obtain a modifiable array.
        total_contribution : float
            Total contribution (`sum(local_contributions)`)

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_contributions`, whose argument names are
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
    baseline_mode = np.array([str(baseline_mode).lower().encode().ljust(4)], dtype="S4")

    # what the inputs already say, rather than asking for it again
    n_dims = factor.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if dependent.shape[0] != n_dims:
        raise ValueError(f"'dependent' has {dependent.shape[0]} along axis 0, but "
            f"'factor' implies n_dims == {n_dims}"
        )

    # outputs and work arrays, which the caller never sees
    local_contributions = np.empty((n_dims,), dtype=np.float64, order='C')
    total_contribution = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.compute_contributions_c(
        factor,
        dependent,
        ctypes.byref(ctypes.c_int(n_dims)),
        baseline_mode,
        local_contributions,
        ctypes.byref(total_contribution),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_CONTRIBUTIONS_ARGUMENTS, _COMPUTE_CONTRIBUTIONS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    local_contributions.flags.writeable = False

    return {
        "local_contributions": local_contributions,
        "total_contribution": total_contribution.value,
    }

def compute_all_contributions(
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

        local_contributions : np.ndarray[np.float64] of shape (n_timepoints, n_selected_factors, n_selected_dependents, n_samples,), column-major (order='F'), read-only
            Per-timepoint contributions per sample-dependent-factor combination
            A result is a value; call `.copy()` to obtain a modifiable array.
        total_contributions : np.ndarray[np.float64] of shape (n_selected_factors, n_selected_dependents, n_samples,), column-major (order='F'), read-only
            Total contribution (`sum(local_contributions)`) per sample-dependent-factor combination
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_all_contributions`, whose argument names are
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
    baseline_mode = np.array([str(baseline_mode).lower().encode().ljust(4)], dtype="S4")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]
    n_selected_factors = factor_indices.shape[0]
    n_selected_dependents = dependent_indices.shape[0]

    # outputs and work arrays, which the caller never sees
    local_contributions = np.empty((n_timepoints, n_selected_factors, n_selected_dependents, n_samples,), dtype=np.float64, order='F')
    total_contributions = np.empty((n_selected_factors, n_selected_dependents, n_samples,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.compute_all_contributions_c(
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
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ALL_CONTRIBUTIONS_ARGUMENTS, _COMPUTE_ALL_CONTRIBUTIONS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    local_contributions.flags.writeable = False
    total_contributions.flags.writeable = False

    return {
        "local_contributions": local_contributions,
        "total_contributions": total_contributions,
    }

def compute_baselines_factor_dependent(
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
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_baselines_factor_dependent`, whose argument names are
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
    baseline_mode = np.array([str(baseline_mode).lower().encode().ljust(4)], dtype="S4")

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

    _lib.compute_baselines_factor_dependent_c(
        ctypes.byref(ctypes.c_int(n_timepoints)),
        factor,
        dependent,
        baseline_mode,
        ctypes.byref(factor_baseline),
        ctypes.byref(dependent_baseline),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENTS, _COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENT_SOURCES)

    return {
        "factor_baseline": factor_baseline.value,
        "dependent_baseline": dependent_baseline.value,
    }

def compute_velocity_trajectory(
        trajectory,
):
    r"""Compute velocity trajectory from a single position trajectory

    Parameters
    ----------
    trajectory : np.ndarray[np.float64] of shape (n_timepoints,)
        input position trajectory

    Returns
    -------
    velocity : np.ndarray[np.float64] of shape (max(0, n_timepoints - 1),), read-only
        output velocity trajectory
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_velocity_trajectory`, whose argument names are
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

    _lib.compute_velocity_trajectory_c(
        trajectory,
        velocity,
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_VELOCITY_TRAJECTORY_ARGUMENTS, _COMPUTE_VELOCITY_TRAJECTORY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    velocity.flags.writeable = False

    return velocity

def compute_acceleration_from_velocity_trajectory(
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
    acceleration : np.ndarray[np.float64] of shape (max(0, n_timepoints - 2),), read-only
        acceleration trajectory
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_acceleration_from_velocity_trajectory`, whose argument names are
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

    _lib.compute_acceleration_from_velocity_trajectory_c(
        velocity,
        acceleration,
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_ARGUMENTS)

    # a result is a value: modify a copy, not this
    acceleration.flags.writeable = False

    return acceleration

def compute_velocity_trajectories(
        trajectories,
):
    r"""Computes velocity trajectories from position trajectories

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        input position trajectories

    Returns
    -------
    velocity : np.ndarray[np.float64] of shape (max(0, n_timepoints - 1), n_factors, n_samples,), column-major (order='F'), read-only
        output velocity trajectories
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_velocity_trajectories`, whose argument names are
    the ones an error message reports.
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
    velocity = np.empty((max(0, n_timepoints - 1), n_factors, n_samples,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.compute_velocity_trajectories_c(
        trajectories,
        velocity,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENTS, _COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    velocity.flags.writeable = False

    return velocity

def compute_acceleration_from_velocity(
        velocity,
        n_timepoints,
):
    r"""Computes acceleration trajectories from velocity trajectories

    Parameters
    ----------
    velocity : np.ndarray[np.float64] of shape (max(0, n_timepoints - 1), n_factors, n_samples,), column-major (order='F')
        input velocity trajectories
    n_timepoints : int
        number of timepoints

    Returns
    -------
    acceleration : np.ndarray[np.float64] of shape (max(0, n_timepoints - 2), n_factors, n_samples,), column-major (order='F'), read-only
        output acceleration trajectories
        A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_acceleration_from_velocity`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        velocity = np.asfortranarray(velocity, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'velocity' must be an array of np.float64: {error}") from None
    if velocity.ndim != 3:
        raise ValueError(f"'velocity' must have 3 dimensions, but has {velocity.ndim}")

    # what the inputs already say, rather than asking for it again
    n_factors = velocity.shape[1]
    n_samples = velocity.shape[2]

    # outputs and work arrays, which the caller never sees
    acceleration = np.empty((max(0, n_timepoints - 2), n_factors, n_samples,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.compute_acceleration_from_velocity_c(
        velocity,
        acceleration,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENTS, _COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    acceleration.flags.writeable = False

    return acceleration

def compute_velocity_acceleration_contributions(
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

        contrib_velocity : np.ndarray[np.float64] of shape (n_factors, n_factors, n_samples,), column-major (order='F'), read-only
            output velocity contributions
            A result is a value; call `.copy()` to obtain a modifiable array.
        velocity_contribution_series : np.ndarray[np.float64] of shape (n_timepoints, n_factors, n_factors, n_samples,), column-major (order='F'), read-only
            output velocity contribution series
            A result is a value; call `.copy()` to obtain a modifiable array.
        contrib_acceleration : np.ndarray[np.float64] of shape (n_factors, n_factors, n_samples,), column-major (order='F'), read-only
            output acceleration contributions
            A result is a value; call `.copy()` to obtain a modifiable array.
        acceleration_contribution_series : np.ndarray[np.float64] of shape (n_timepoints, n_factors, n_factors, n_samples,), column-major (order='F'), read-only
            output acceleration contribution series
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_velocity_acceleration_contributions`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectories' must be an array of np.float64: {error}") from None
    if trajectories.ndim != 3:
        raise ValueError(f"'trajectories' must have 3 dimensions, but has {trajectories.ndim}")
    baseline_mode = np.array([str(baseline_mode).lower().encode().ljust(4)], dtype="S4")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    contrib_velocity = np.empty((n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    velocity_contribution_series = np.empty((n_timepoints, n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    contrib_acceleration = np.empty((n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    acceleration_contribution_series = np.empty((n_timepoints, n_factors, n_factors, n_samples,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.compute_velocity_acceleration_contributions_c(
        trajectories,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        baseline_mode,
        contrib_velocity,
        velocity_contribution_series,
        contrib_acceleration,
        acceleration_contribution_series,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_ARGUMENTS, _COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    contrib_velocity.flags.writeable = False
    velocity_contribution_series.flags.writeable = False
    contrib_acceleration.flags.writeable = False
    acceleration_contribution_series.flags.writeable = False

    return {
        "contrib_velocity": contrib_velocity,
        "velocity_contribution_series": velocity_contribution_series,
        "contrib_acceleration": contrib_acceleration,
        "acceleration_contribution_series": acceleration_contribution_series,
    }
