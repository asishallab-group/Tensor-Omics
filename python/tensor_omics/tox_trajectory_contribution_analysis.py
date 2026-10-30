r"""tox_trajectory_contribution_analysis

Module for quantifying how much one trajectory (a "factor") contributes to another (a "dependent") over time.

Contributions are computed per timepoint as the product of both series' deviations from a chosen
baseline, for raw expression trajectories as well as for their velocity (first difference) and
acceleration (second difference) derivatives. Statistical significance of an observed contribution can
be assessed via a permutation test that recomputes the same contribution against a randomly chosen
other sample.

Python binding, generated from tox_trajectory_contribution_analysis. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.perform_permutation_test_c.restype = None
_lib.perform_permutation_test_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _nullable(_ctypes.POINTER(_ctypes.c_int)),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_PERFORM_PERMUTATION_TEST_ARGUMENTS = ("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_idx", "dependent_idx", "sample_idx", "baseline_mode", "n_permutations", "local_contributions", "total_contributions", "random_seed", "ierr",)
#: For a derived argument, the one the caller passed it in
_PERFORM_PERMUTATION_TEST_ARGUMENT_SOURCES = (None, "trajectories", "trajectories", "trajectories", None, None, None, None, "local_contributions", None, None, None, None,)

_lib.compute_p_values_c.restype = None
_lib.compute_p_values_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_P_VALUES_ARGUMENTS = ("local_contributions_observed", "total_contribution_observed", "local_contributions_perm_test", "total_contributions_perm_test", "n_timepoints", "n_permutations", "local_p_values", "total_p_value", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_P_VALUES_ARGUMENT_SOURCES = (None, None, None, None, "local_contributions_observed", "local_contributions_perm_test", None, None, None,)

_lib.compute_contributions_c.restype = None
_lib.compute_contributions_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_CONTRIBUTIONS_ARGUMENTS = ("factor", "dependent", "n_dims", "baseline_mode", "local_contributions", "total_contribution", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_CONTRIBUTIONS_ARGUMENT_SOURCES = (None, None, "factor", None, None, None, None,)

_lib.compute_all_contributions_c.restype = None
_lib.compute_all_contributions_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=4, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ALL_CONTRIBUTIONS_ARGUMENTS = ("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_indices", "n_selected_factors", "dependent_indices", "n_selected_dependents", "baseline_mode", "local_contributions", "total_contributions", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ALL_CONTRIBUTIONS_ARGUMENT_SOURCES = (None, "trajectories", "trajectories", "trajectories", None, "factor_indices", None, "dependent_indices", None, None, None, None,)

_lib.compute_baselines_factor_dependent_c.restype = None
_lib.compute_baselines_factor_dependent_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENTS = ("n_timepoints", "factor", "dependent", "baseline_mode", "factor_baseline", "dependent_baseline", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENT_SOURCES = ("factor", None, None, None, None, None, None,)

_lib.compute_velocity_trajectory_c.restype = None
_lib.compute_velocity_trajectory_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_VELOCITY_TRAJECTORY_ARGUMENTS = ("trajectory", "velocity", "n_timepoints", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_VELOCITY_TRAJECTORY_ARGUMENT_SOURCES = (None, None, "trajectory", None,)

_lib.compute_acceleration_from_velocity_trajectory_c.restype = None
_lib.compute_acceleration_from_velocity_trajectory_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_ARGUMENTS = ("velocity", "acceleration", "n_timepoints", "ierr",)
#: For an array the C wrapper compares with its declaration, that shape
_COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_ARGUMENT_SHAPES = ("max(0, n_timepoints - 1)", None, None, None,)

_lib.compute_velocity_trajectories_c.restype = None
_lib.compute_velocity_trajectories_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENTS = ("trajectories", "velocity", "n_factors", "n_samples", "n_timepoints", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENT_SOURCES = (None, None, "trajectories", "trajectories", "trajectories", None,)

_lib.compute_acceleration_from_velocity_c.restype = None
_lib.compute_acceleration_from_velocity_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENTS = ("velocity", "acceleration", "n_factors", "n_samples", "n_timepoints", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENT_SOURCES = (None, None, "velocity", "velocity", None, None,)
#: For an array the C wrapper compares with its declaration, that shape
_COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENT_SHAPES = ("max(0, n_timepoints - 1), n_factors, n_samples", None, None, None, None, None,)

_lib.compute_velocity_acceleration_contributions_c.restype = None
_lib.compute_velocity_acceleration_contributions_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=4, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=4, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
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
    baseline_mode : str, one of 'baseline_raw' | 'baseline_mean' | 'baseline_min'
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::perform_permutation_test`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectories = _contract.reals(trajectories, "trajectories", ndim=3, order="F")
    factor_idx = _contract.integer(factor_idx, "factor_idx")
    dependent_idx = _contract.integer(dependent_idx, "dependent_idx")
    sample_idx = _contract.integer(sample_idx, "sample_idx")
    baseline_mode = _np.array([_contract.mode(baseline_mode, "baseline_mode", ("baseline_raw", "baseline_mean", "baseline_min",)).ljust(13)], dtype="S13")
    n_permutations = _contract.integer(n_permutations, "n_permutations")
    if random_seed is not None:
        random_seed = _contract.integer(random_seed, "random_seed")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    local_contributions = _np.empty((n_timepoints, n_permutations,), dtype=_np.float64, order='F')
    total_contributions = _np.empty((n_permutations,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.perform_permutation_test_c(
        trajectories,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        _ctypes.byref(_ctypes.c_int(factor_idx)),
        _ctypes.byref(_ctypes.c_int(dependent_idx)),
        _ctypes.byref(_ctypes.c_int(sample_idx)),
        baseline_mode,
        _ctypes.byref(_ctypes.c_int(n_permutations)),
        local_contributions,
        total_contributions,
        None if random_seed is None else _ctypes.byref(_ctypes.c_int(random_seed)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _PERFORM_PERMUTATION_TEST_ARGUMENTS, _PERFORM_PERMUTATION_TEST_ARGUMENT_SOURCES)

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

    Given the permutation tests (:func:`tensor_omics.perform_permutation_test`),
    this counts how many of the permutation contributions were at least as high as the real ones.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_p_values`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    local_contributions_observed = _contract.reals(local_contributions_observed, "local_contributions_observed", ndim=1)
    total_contribution_observed = _contract.real(total_contribution_observed, "total_contribution_observed")
    local_contributions_perm_test = _contract.reals(local_contributions_perm_test, "local_contributions_perm_test", ndim=2, order="F")
    total_contributions_perm_test = _contract.reals(total_contributions_perm_test, "total_contributions_perm_test", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_timepoints = local_contributions_observed.shape[0]
    n_permutations = local_contributions_perm_test.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if local_contributions_perm_test.shape[0] != n_timepoints:
        raise _ToxArgumentError(f"'local_contributions_perm_test' has {local_contributions_perm_test.shape[0]} along axis 0, but 'local_contributions_observed' implies n_timepoints == {n_timepoints}; 'local_contributions_perm_test' is declared with shape (n_timepoints, n_permutations)", argument="local_contributions_perm_test")
    if total_contributions_perm_test.shape[0] != n_permutations:
        raise _ToxArgumentError(f"'total_contributions_perm_test' has {total_contributions_perm_test.shape[0]} along axis 0, but 'local_contributions_perm_test' implies n_permutations == {n_permutations}; 'total_contributions_perm_test' is declared with shape (n_permutations)", argument="total_contributions_perm_test")

    # outputs and work arrays, which the caller never sees
    local_p_values = _np.empty((n_timepoints,), dtype=_np.float64, order='C')
    total_p_value = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.compute_p_values_c(
        local_contributions_observed,
        _ctypes.byref(_ctypes.c_double(total_contribution_observed)),
        local_contributions_perm_test,
        total_contributions_perm_test,
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        _ctypes.byref(_ctypes.c_int(n_permutations)),
        local_p_values,
        _ctypes.byref(total_p_value),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_P_VALUES_ARGUMENTS, _COMPUTE_P_VALUES_ARGUMENT_SOURCES)

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
    baseline_mode : str, one of 'baseline_raw' | 'baseline_mean' | 'baseline_min'
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_contributions`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    factor = _contract.reals(factor, "factor", ndim=1)
    dependent = _contract.reals(dependent, "dependent", ndim=1)
    baseline_mode = _np.array([_contract.mode(baseline_mode, "baseline_mode", ("baseline_raw", "baseline_mean", "baseline_min",)).ljust(13)], dtype="S13")

    # what the inputs already say, rather than asking for it again
    n_dims = factor.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dependent.shape[0] != n_dims:
        raise _ToxArgumentError(f"'dependent' has {dependent.shape[0]} along axis 0, but 'factor' implies n_dims == {n_dims}; 'dependent' is declared with shape (n_dims)", argument="dependent")

    # outputs and work arrays, which the caller never sees
    local_contributions = _np.empty((n_dims,), dtype=_np.float64, order='C')
    total_contribution = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.compute_contributions_c(
        factor,
        dependent,
        _ctypes.byref(_ctypes.c_int(n_dims)),
        baseline_mode,
        local_contributions,
        _ctypes.byref(total_contribution),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_CONTRIBUTIONS_ARGUMENTS, _COMPUTE_CONTRIBUTIONS_ARGUMENT_SOURCES)

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
    baseline_mode : str, one of 'baseline_raw' | 'baseline_mean' | 'baseline_min'
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_all_contributions`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectories = _contract.reals(trajectories, "trajectories", ndim=3, order="F")
    factor_indices = _contract.integers(factor_indices, "factor_indices", ndim=1)
    dependent_indices = _contract.integers(dependent_indices, "dependent_indices", ndim=1)
    baseline_mode = _np.array([_contract.mode(baseline_mode, "baseline_mode", ("baseline_raw", "baseline_mean", "baseline_min",)).ljust(13)], dtype="S13")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]
    n_selected_factors = factor_indices.shape[0]
    n_selected_dependents = dependent_indices.shape[0]

    # outputs and work arrays, which the caller never sees
    local_contributions = _np.empty((n_timepoints, n_selected_factors, n_selected_dependents, n_samples,), dtype=_np.float64, order='F')
    total_contributions = _np.empty((n_selected_factors, n_selected_dependents, n_samples,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.compute_all_contributions_c(
        trajectories,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        factor_indices,
        _ctypes.byref(_ctypes.c_int(n_selected_factors)),
        dependent_indices,
        _ctypes.byref(_ctypes.c_int(n_selected_dependents)),
        baseline_mode,
        local_contributions,
        total_contributions,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_ALL_CONTRIBUTIONS_ARGUMENTS, _COMPUTE_ALL_CONTRIBUTIONS_ARGUMENT_SOURCES)

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
    baseline_mode : str, one of 'baseline_raw' | 'baseline_mean' | 'baseline_min'
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_baselines_factor_dependent`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    factor = _contract.reals(factor, "factor", ndim=1)
    dependent = _contract.reals(dependent, "dependent", ndim=1)
    baseline_mode = _np.array([_contract.mode(baseline_mode, "baseline_mode", ("baseline_raw", "baseline_mean", "baseline_min",)).ljust(13)], dtype="S13")

    # what the inputs already say, rather than asking for it again
    n_timepoints = factor.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dependent.shape[0] != n_timepoints:
        raise _ToxArgumentError(f"'dependent' has {dependent.shape[0]} along axis 0, but 'factor' implies n_timepoints == {n_timepoints}; 'dependent' is declared with shape (n_timepoints)", argument="dependent")

    # outputs and work arrays, which the caller never sees
    factor_baseline = _ctypes.c_double(0)
    dependent_baseline = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.compute_baselines_factor_dependent_c(
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        factor,
        dependent,
        baseline_mode,
        _ctypes.byref(factor_baseline),
        _ctypes.byref(dependent_baseline),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENTS, _COMPUTE_BASELINES_FACTOR_DEPENDENT_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_velocity_trajectory`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectory = _contract.reals(trajectory, "trajectory", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_timepoints = trajectory.shape[0]

    # outputs and work arrays, which the caller never sees
    velocity = _np.empty((max(0, n_timepoints - 1),), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_velocity_trajectory_c(
        trajectory,
        velocity,
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_VELOCITY_TRAJECTORY_ARGUMENTS, _COMPUTE_VELOCITY_TRAJECTORY_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_acceleration_from_velocity_trajectory`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    velocity = _contract.reals(velocity, "velocity", ndim=1)
    n_timepoints = _contract.integer(n_timepoints, "n_timepoints")

    # what the inputs already say, rather than asking for it again
    n_velocity_elements = velocity.shape[0]

    # outputs and work arrays, which the caller never sees
    acceleration = _np.empty((max(0, n_timepoints - 2),), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_acceleration_from_velocity_trajectory_c(
        velocity,
        _ctypes.byref(_ctypes.c_int(n_velocity_elements)),
        acceleration,
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_ARGUMENTS, shapes=_COMPUTE_ACCELERATION_FROM_VELOCITY_TRAJECTORY_ARGUMENT_SHAPES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_velocity_trajectories`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectories = _contract.reals(trajectories, "trajectories", ndim=3, order="F")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    velocity = _np.empty((max(0, n_timepoints - 1), n_factors, n_samples,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.compute_velocity_trajectories_c(
        trajectories,
        velocity,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENTS, _COMPUTE_VELOCITY_TRAJECTORIES_ARGUMENT_SOURCES)

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_acceleration_from_velocity`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    velocity = _contract.reals(velocity, "velocity", ndim=3, order="F")
    n_timepoints = _contract.integer(n_timepoints, "n_timepoints")

    # what the inputs already say, rather than asking for it again
    n_velocity_elements_dim_1 = velocity.shape[0]
    n_factors = velocity.shape[1]
    n_samples = velocity.shape[2]

    # outputs and work arrays, which the caller never sees
    acceleration = _np.empty((max(0, n_timepoints - 2), n_factors, n_samples,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.compute_acceleration_from_velocity_c(
        velocity,
        _ctypes.byref(_ctypes.c_int(n_velocity_elements_dim_1)),
        acceleration,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENTS, _COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENT_SOURCES, shapes=_COMPUTE_ACCELERATION_FROM_VELOCITY_ARGUMENT_SHAPES)

    # a result is a value: modify a copy, not this
    acceleration.flags.writeable = False

    return acceleration

def compute_velocity_acceleration_contributions(
        trajectories,
        baseline_mode,
):
    r"""Compute velocity and acceleration contributions for all variable pairs

    Performance layout:

    `trajectories` uses `(n_factors, n_samples, n_timepoints)`.
    Velocity and acceleration use time-first layouts:

    - `velocity`     -> `(max(0, n_timepoints-1), n_factors, n_samples)`
    - `acceleration` -> `(max(0, n_timepoints-2), n_factors, n_samples)`

    This keeps slices like `velocity(:, factor, sample)` contiguous,
    avoids expensive tmporaries, and improves cache efficiency.

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        input position trajectories
    baseline_mode : str, one of 'baseline_raw' | 'baseline_mean' | 'baseline_min'
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_trajectory_contribution_analysis::compute_velocity_acceleration_contributions`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectories = _contract.reals(trajectories, "trajectories", ndim=3, order="F")
    baseline_mode = _np.array([_contract.mode(baseline_mode, "baseline_mode", ("baseline_raw", "baseline_mean", "baseline_min",)).ljust(13)], dtype="S13")

    # what the inputs already say, rather than asking for it again
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # outputs and work arrays, which the caller never sees
    contrib_velocity = _np.empty((n_factors, n_factors, n_samples,), dtype=_np.float64, order='F')
    velocity_contribution_series = _np.empty((n_timepoints, n_factors, n_factors, n_samples,), dtype=_np.float64, order='F')
    contrib_acceleration = _np.empty((n_factors, n_factors, n_samples,), dtype=_np.float64, order='F')
    acceleration_contribution_series = _np.empty((n_timepoints, n_factors, n_factors, n_samples,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.compute_velocity_acceleration_contributions_c(
        trajectories,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        baseline_mode,
        contrib_velocity,
        velocity_contribution_series,
        contrib_acceleration,
        acceleration_contribution_series,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_ARGUMENTS, _COMPUTE_VELOCITY_ACCELERATION_CONTRIBUTIONS_ARGUMENT_SOURCES)

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
