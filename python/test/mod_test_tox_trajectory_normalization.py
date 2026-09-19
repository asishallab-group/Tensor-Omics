"""
The `tox_trajectory_normalization` suite: each published procedure can be called from Python,
returns the documented type and shape, and raises the documented error. Whether the numbers are
right is the Fortran suite's job (test/mod_test_tox_trajectory_normalization.F90), as the coding
guide's "Where to Test What" asks.

A series that cannot be scaled (constant over time) is not an error: the procedures return its
per-series `status`, so those tests check that the status reaches the caller, at the right index.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    normalize_all_trajectories,
    normalize_single_trajectory,
    normalize_variable_timeseries,
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import (
    ERR_DIVISION_BY_ZERO,
    ERR_EMPTY_INPUT,
    ERR_NAN_INF,
    ERR_OK,
)

N_FACTORS, N_SAMPLES, N_TIMEPOINTS = 2, 3, 4


def _assert_result_array(result, shape, dtype, name, fortran_order):
    assert isinstance(result, np.ndarray), f"{name}: expected a numpy array, got {type(result).__name__}"
    assert result.dtype == dtype, f"{name}: expected {np.dtype(dtype)}, got {result.dtype}"
    assert result.shape == shape, f"{name}: expected shape {shape}, got {result.shape}"
    assert not result.flags.writeable, f"{name}: a result is a value, it should be read-only"
    if fortran_order:
        assert result.flags.f_contiguous, f"{name}: expected column-major (order='F')"
    else:
        assert result.flags.c_contiguous, f"{name}: expected a contiguous array"


def _assert_keys(result, keys):
    assert isinstance(result, dict), f"expected a dict, got {type(result).__name__}"
    assert set(result) == keys, f"unexpected keys {sorted(result)}"


def _trajectories():
    """(factors, samples, timepoints), every series varying over time."""
    return np.fromfunction(
        lambda i_factor, i_sample, i_timepoint: (i_factor + 1) * (i_sample + 1) * (i_timepoint + 1.0),
        (N_FACTORS, N_SAMPLES, N_TIMEPOINTS))


# -----------------------------------------------------------------------------------------------
# normalize_variable_timeseries
# -----------------------------------------------------------------------------------------------

def test_normalize_variable_timeseries():
    result = normalize_variable_timeseries([3.0, 1.0, 2.0])
    _assert_keys(result, {"v_norm", "status"})
    _assert_result_array(result["v_norm"], (3,), np.float64, "v_norm", fortran_order=False)
    assert isinstance(result["status"], int), f"status: expected an int, got {type(result['status']).__name__}"
    assert result["status"] == ERR_OK, f"status: expected ERR_OK, got {result['status']}"
    # the one value: the binding hands the points over in order, the maximum first
    assert np.array_equal(result["v_norm"], [1.0, 0.0, 0.5]), f"v_norm: got {result['v_norm']}"


def test_normalize_variable_timeseries_reports_a_constant_series_in_status():
    # not raised: the status is returned next to the result
    result = normalize_variable_timeseries(np.full(3, 7.0))
    assert result["status"] == ERR_DIVISION_BY_ZERO, f"status: expected ERR_DIVISION_BY_ZERO, got {result['status']}"


def test_normalize_variable_timeseries_rejects_an_empty_series():
    assert_error(lambda: normalize_variable_timeseries(np.empty(0)), "no points", ERR_EMPTY_INPUT)


def test_normalize_variable_timeseries_rejects_nan():
    assert_error(lambda: normalize_variable_timeseries([1.0, np.nan, 3.0]), "a NaN point", ERR_NAN_INF)


def test_normalize_variable_timeseries_rejects_a_matrix():
    # the binding's own rank check, before the library is called
    assert_error(lambda: normalize_variable_timeseries(np.ones((2, 2))), "v must be 1-D")


# -----------------------------------------------------------------------------------------------
# normalize_single_trajectory
# -----------------------------------------------------------------------------------------------

def test_normalize_single_trajectory():
    # (timepoints, factors), given in C order so the binding has to convert it
    trajectory = np.array([[0.0, 2.0],
                           [1.0, 1.0],
                           [2.0, 0.0]])
    result = normalize_single_trajectory(trajectory)
    _assert_keys(result, {"trajectory_norm", "status"})
    _assert_result_array(result["trajectory_norm"], (3, 2), np.float64, "trajectory_norm", fortran_order=True)
    _assert_result_array(result["status"], (2,), np.int32, "status", fortran_order=False)
    assert np.all(result["status"] == ERR_OK), f"status: got {result['status']}"
    # the one value: time runs down the rows, so the first factor rises and the second falls
    assert np.array_equal(result["trajectory_norm"], [[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]]), \
        f"trajectory_norm: got {result['trajectory_norm']}"


def test_normalize_single_trajectory_reports_a_constant_factor_in_status():
    # the second factor is constant over time; only its status says so
    result = normalize_single_trajectory(np.array([[0.0, 5.0], [1.0, 5.0], [2.0, 5.0]]))
    assert list(result["status"]) == [ERR_OK, ERR_DIVISION_BY_ZERO], f"status: got {result['status']}"


def test_normalize_single_trajectory_rejects_no_factors():
    assert_error(lambda: normalize_single_trajectory(np.empty((3, 0))), "no factors", ERR_EMPTY_INPUT)


def test_normalize_single_trajectory_rejects_infinity():
    trajectory = np.ones((3, 2))
    trajectory[1, 1] = np.inf
    assert_error(lambda: normalize_single_trajectory(trajectory), "an infinite point", ERR_NAN_INF)


def test_normalize_single_trajectory_rejects_a_vector():
    # the binding's own rank check, before the library is called
    assert_error(lambda: normalize_single_trajectory(np.ones(3)), "trajectory must be 2-D")


# -----------------------------------------------------------------------------------------------
# normalize_all_trajectories
# -----------------------------------------------------------------------------------------------

def test_normalize_all_trajectories():
    trajectories = _trajectories()
    trajectories[0, 1, :] = [3.0, 2.0, 1.0, 0.0]
    result = normalize_all_trajectories(trajectories)
    _assert_keys(result, {"trajectories_norm", "status"})
    _assert_result_array(result["trajectories_norm"], (N_FACTORS, N_SAMPLES, N_TIMEPOINTS), np.float64,
                         "trajectories_norm", fortran_order=True)
    _assert_result_array(result["status"], (N_FACTORS, N_SAMPLES), np.int32, "status", fortran_order=True)
    assert np.all(result["status"] == ERR_OK), f"status: got {result['status']}"
    # the one value: time is the last axis, and this series falls over it
    assert np.allclose(result["trajectories_norm"][0, 1, :], [1.0, 2.0 / 3.0, 1.0 / 3.0, 0.0]), \
        f"trajectories_norm[0, 1, :]: got {result['trajectories_norm'][0, 1, :]}"


def test_normalize_all_trajectories_reports_a_constant_series_in_status():
    # one (factor, sample) constant over time: its status, and only its, says so
    trajectories = _trajectories()
    trajectories[1, 2, :] = 5.0
    status = normalize_all_trajectories(trajectories)["status"]
    expected = np.full((N_FACTORS, N_SAMPLES), ERR_OK)
    expected[1, 2] = ERR_DIVISION_BY_ZERO
    assert np.array_equal(status, expected), f"status: got {status}"


def test_normalize_all_trajectories_rejects_no_samples():
    assert_error(lambda: normalize_all_trajectories(np.empty((N_FACTORS, 0, N_TIMEPOINTS))),
                 "no samples", ERR_EMPTY_INPUT)


def test_normalize_all_trajectories_rejects_nan():
    trajectories = _trajectories()
    trajectories[1, 0, 2] = np.nan
    assert_error(lambda: normalize_all_trajectories(trajectories), "a NaN point", ERR_NAN_INF)


def test_normalize_all_trajectories_rejects_a_matrix():
    # the binding's own rank check, before the library is called
    assert_error(lambda: normalize_all_trajectories(np.ones((2, 3))), "trajectories must be 3-D")


if __name__ == '__main__':
    run_all_tests(globals().values())
