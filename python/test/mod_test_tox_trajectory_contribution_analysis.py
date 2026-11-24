"""
The `tox_trajectory_contribution_analysis` suite: each published procedure can be called from
Python, returns the documented type and shape, and raises the documented error, naming the
argument. Whether the numbers are right is the Fortran suite's job
(test/mod_test_tox_trajectory_contribution_analysis/), as the coding guide's "Where to Test What"
asks; the one value a test checks here is there to show that the binding hands the arrays over
in the right order and the indices 1-based.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    compute_acceleration_from_velocity,
    compute_acceleration_from_velocity_trajectory,
    compute_all_contributions,
    compute_baselines_factor_dependent,
    compute_contributions,
    compute_p_values,
    compute_velocity_acceleration_contributions,
    compute_velocity_trajectories,
    compute_velocity_trajectory,
    perform_permutation_test,
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import (
    ERR_EMPTY_INPUT,
    ERR_INVALID_INPUT,
    ERR_NAN_INF,
    ERR_SIZE_MISMATCH,
)

N_FACTORS, N_SAMPLES, N_TIMEPOINTS = 2, 3, 4


def _assert_result_array(result, shape, name, fortran_order):
    assert isinstance(result, np.ndarray), f"{name}: expected a numpy array, got {type(result).__name__}"
    assert result.dtype == np.float64, f"{name}: expected float64, got {result.dtype}"
    assert result.shape == shape, f"{name}: expected shape {shape}, got {result.shape}"
    assert not result.flags.writeable, f"{name}: a result is a value, it should be read-only"
    if fortran_order:
        assert result.flags.f_contiguous, f"{name}: expected column-major (order='F')"
    else:
        assert result.flags.c_contiguous, f"{name}: expected a contiguous array"


def _assert_float(value, name):
    assert isinstance(value, float), f"{name}: expected a float, got {type(value).__name__}"


def _assert_keys(result, keys):
    assert isinstance(result, dict), f"expected a dict, got {type(result).__name__}"
    assert set(result) == keys, f"unexpected keys {sorted(result)}"


def _trajectories():
    """(factors, samples, timepoints), every series varying over time."""
    return np.fromfunction(
        lambda i_factor, i_sample, i_timepoint: (i_factor + 1) * (i_sample + 1) * (i_timepoint + 1.0) ** 2,
        (N_FACTORS, N_SAMPLES, N_TIMEPOINTS))


def _sample_marking(n_samples):
    """Factor 1 is [1, 0] in every sample, factor 2 of sample s is [s, 0]: a raw permutation's
    total is the number of the sample it drew."""
    trajectories = np.zeros((2, n_samples, 2))
    trajectories[0, :, 0] = 1.0
    trajectories[1, :, 0] = np.arange(1, n_samples + 1)
    return trajectories


# -----------------------------------------------------------------------------------------------
# perform_permutation_test
# -----------------------------------------------------------------------------------------------

def test_perform_permutation_test():
    result = perform_permutation_test(_trajectories(), 1, 2, 2, "baseline_mean", 5, random_seed=42)
    _assert_keys(result, {"local_contributions", "total_contributions"})
    _assert_result_array(result["local_contributions"], (N_TIMEPOINTS, 5), "local_contributions", fortran_order=True)
    _assert_result_array(result["total_contributions"], (5,), "total_contributions", fortran_order=False)


def test_perform_permutation_test_draws_the_other_sample():
    # the one value: of two samples, sample 1 (1-based) can only draw sample 2
    result = perform_permutation_test(_sample_marking(2), 1, 2, 1, "baseline_raw", 4, random_seed=42)
    assert np.array_equal(result["total_contributions"], [2.0] * 4), f"got {result['total_contributions']}"


def test_perform_permutation_test_random_seed_reaches_the_library():
    # the same seed twice draws the same: the binding hands it over, it does not drop it
    trajectories = _sample_marking(5)
    first = perform_permutation_test(trajectories, 1, 2, 3, "baseline_raw", 64, random_seed=7)
    second = perform_permutation_test(trajectories, 1, 2, 3, "baseline_raw", 64, random_seed=7)
    assert np.array_equal(first["total_contributions"], second["total_contributions"]), "the same seed drew otherwise"


def test_perform_permutation_test_without_a_seed():
    result = perform_permutation_test(_sample_marking(3), 1, 2, 1, "baseline_raw", 3)
    assert set(result["total_contributions"]) <= {2.0, 3.0}, f"got {result['total_contributions']}"


def test_perform_permutation_test_rejects_a_sample_index_out_of_range():
    assert_error(lambda: perform_permutation_test(_trajectories(), 1, 2, N_SAMPLES + 1, "baseline_raw", 2,
                                                  random_seed=42),
                 "sample_idx past n_samples", ERR_INVALID_INPUT, argument="sample_idx")


def test_perform_permutation_test_rejects_zero_permutations():
    assert_error(lambda: perform_permutation_test(_trajectories(), 1, 2, 1, "baseline_raw", 0, random_seed=42),
                 "no permutations", ERR_EMPTY_INPUT, argument="n_permutations")


def test_perform_permutation_test_rejects_a_negative_permutation_count():
    # the binding's own extent check, before it allocates by it
    assert_error(lambda: perform_permutation_test(_trajectories(), 1, 2, 1, "baseline_raw", -1, random_seed=42),
                 "a negative n_permutations", argument="n_permutations")


def test_perform_permutation_test_rejects_nan():
    trajectories = _trajectories()
    trajectories[1, 2, 3] = np.nan
    assert_error(lambda: perform_permutation_test(trajectories, 1, 2, 1, "baseline_raw", 2, random_seed=42),
                 "a NaN point", ERR_NAN_INF, argument="trajectories")


def test_perform_permutation_test_rejects_an_unknown_mode():
    # the binding's own check, before the library is called
    assert_error(lambda: perform_permutation_test(_trajectories(), 1, 2, 1, "baseline_median", 2, random_seed=42),
                 "an unknown mode", argument="baseline_mode")


# -----------------------------------------------------------------------------------------------
# compute_p_values
# -----------------------------------------------------------------------------------------------

def test_compute_p_values():
    # (timepoints, permutations): the columns are the permutations
    local_perm = np.array([[0.0, 1.0, 2.0, 3.0],
                           [5.0, 5.0, 5.0, 5.0]])
    result = compute_p_values([2.0, 1.0], 3.0, local_perm, [1.0, 2.0, 3.0, 4.0])
    _assert_keys(result, {"local_p_values", "total_p_value"})
    _assert_result_array(result["local_p_values"], (2,), "local_p_values", fortran_order=False)
    _assert_float(result["total_p_value"], "total_p_value")
    # the one value: point 1 is beaten by permutations 3 and 4, point 2 by all, the total by 3 and 4
    assert np.array_equal(result["local_p_values"], [0.5, 1.0]), f"local_p_values: got {result['local_p_values']}"
    assert result["total_p_value"] == 0.5, f"total_p_value: got {result['total_p_value']}"


def test_compute_p_values_rejects_nan():
    assert_error(lambda: compute_p_values([2.0, np.nan], 3.0, np.ones((2, 4)), np.ones(4)),
                 "a NaN observation", ERR_NAN_INF, argument="local_contributions_observed")


def test_compute_p_values_rejects_infinity_in_the_permutations():
    local_perm = np.ones((2, 4))
    local_perm[1, 3] = np.inf
    assert_error(lambda: compute_p_values([2.0, 1.0], 3.0, local_perm, np.ones(4)),
                 "an infinite permutation", ERR_NAN_INF, argument="local_contributions_perm_test")


def test_compute_p_values_rejects_mismatched_permutations():
    # the binding's own extent checks, before the library is called
    assert_error(lambda: compute_p_values([2.0, 1.0], 3.0, np.ones((3, 4)), np.ones(4)),
                 "a permutation of three points against two", argument="local_contributions_perm_test")
    assert_error(lambda: compute_p_values([2.0, 1.0], 3.0, np.ones((2, 4)), np.ones(3)),
                 "three totals for four permutations", argument="total_contributions_perm_test")


# -----------------------------------------------------------------------------------------------
# compute_contributions
# -----------------------------------------------------------------------------------------------

def test_compute_contributions():
    result = compute_contributions([1.0, 2.0, 3.0], [3.0, 1.0, 2.0], "baseline_raw")
    _assert_keys(result, {"local_contributions", "total_contribution"})
    _assert_result_array(result["local_contributions"], (3,), "local_contributions", fortran_order=False)
    _assert_float(result["total_contribution"], "total_contribution")
    # the one value: the points stay in order
    assert np.array_equal(result["local_contributions"], [3.0, 2.0, 6.0]), f"got {result['local_contributions']}"


def test_compute_contributions_rejects_no_points():
    assert_error(lambda: compute_contributions(np.empty(0), np.empty(0), "baseline_raw"),
                 "no points", ERR_EMPTY_INPUT, argument="factor")


def test_compute_contributions_rejects_nan():
    assert_error(lambda: compute_contributions([1.0, 2.0], [np.nan, 1.0], "baseline_raw"),
                 "a NaN point", ERR_NAN_INF, argument="dependent")


def test_compute_contributions_rejects_mismatched_series():
    # the binding's own extent check, before the library is called
    assert_error(lambda: compute_contributions([1.0, 2.0, 3.0], [1.0, 2.0], "baseline_raw"),
                 "a dependent shorter than the factor", argument="dependent")


def test_compute_contributions_rejects_an_unknown_mode():
    assert_error(lambda: compute_contributions([1.0, 2.0], [3.0, 4.0], "mean"),
                 "an unknown mode", argument="baseline_mode")


# -----------------------------------------------------------------------------------------------
# compute_all_contributions
# -----------------------------------------------------------------------------------------------

def test_compute_all_contributions():
    result = compute_all_contributions(_trajectories(), [2, 1, 2], [1, 2], "baseline_min")
    _assert_keys(result, {"local_contributions", "total_contributions"})
    _assert_result_array(result["local_contributions"], (N_TIMEPOINTS, 3, 2, N_SAMPLES), "local_contributions",
                         fortran_order=True)
    _assert_result_array(result["total_contributions"], (3, 2, N_SAMPLES), "total_contributions",
                         fortran_order=True)


def test_compute_all_contributions_takes_indices_1_based():
    # the one value: factor 1 against factor 2, raw, in the one sample: [1*3, 2*4] = [3, 8]
    trajectories = np.array([[[1.0, 2.0]], [[3.0, 4.0]]])
    result = compute_all_contributions(trajectories, [1], [2], "baseline_raw")
    assert np.array_equal(result["local_contributions"][:, 0, 0, 0], [3.0, 8.0]), \
        f"got {result['local_contributions'][:, 0, 0, 0]}"


def test_compute_all_contributions_rejects_an_index_out_of_range():
    assert_error(lambda: compute_all_contributions(_trajectories(), [1, N_FACTORS + 1], [1], "baseline_raw"),
                 "a factor index past n_factors", ERR_INVALID_INPUT, argument="factor_indices")
    assert_error(lambda: compute_all_contributions(_trajectories(), [1], [0], "baseline_raw"),
                 "a dependent index 0", ERR_INVALID_INPUT, argument="dependent_indices")


def test_compute_all_contributions_rejects_an_empty_selection():
    assert_error(lambda: compute_all_contributions(_trajectories(), np.empty(0, dtype=np.int32), [1], "baseline_raw"),
                 "no factors selected", ERR_EMPTY_INPUT, argument="factor_indices")


def test_compute_all_contributions_rejects_infinity():
    trajectories = _trajectories()
    trajectories[0, 1, 2] = -np.inf
    assert_error(lambda: compute_all_contributions(trajectories, [1], [2], "baseline_raw"),
                 "an infinite point", ERR_NAN_INF, argument="trajectories")


# -----------------------------------------------------------------------------------------------
# compute_baselines_factor_dependent
# -----------------------------------------------------------------------------------------------

def test_compute_baselines_factor_dependent():
    result = compute_baselines_factor_dependent([3.0, -1.0, 4.0, 2.0], [5.0, 9.0, 2.0, 6.0], "baseline_mean")
    _assert_keys(result, {"factor_baseline", "dependent_baseline"})
    _assert_float(result["factor_baseline"], "factor_baseline")
    _assert_float(result["dependent_baseline"], "dependent_baseline")
    # the one value: each series has its own baseline
    assert (result["factor_baseline"], result["dependent_baseline"]) == (2.0, 5.5), f"got {result}"


def test_compute_baselines_factor_dependent_rejects_no_points():
    assert_error(lambda: compute_baselines_factor_dependent(np.empty(0), np.empty(0), "baseline_min"),
                 "no points", ERR_EMPTY_INPUT, argument="factor")


def test_compute_baselines_factor_dependent_rejects_nan():
    assert_error(lambda: compute_baselines_factor_dependent([1.0, np.nan], [1.0, 2.0], "baseline_raw"),
                 "a NaN point", ERR_NAN_INF, argument="factor")


def test_compute_baselines_factor_dependent_rejects_mismatched_series():
    # the binding's own extent check, before the library is called
    assert_error(lambda: compute_baselines_factor_dependent([1.0, 2.0, 3.0], [1.0, 2.0], "baseline_raw"),
                 "a dependent shorter than the factor", argument="dependent")


def test_compute_baselines_factor_dependent_rejects_an_unknown_mode():
    # it reached the C wrapper cut to the buffer's 13 bytes before (#221 G-I1, G-I2)
    assert_error(lambda: compute_baselines_factor_dependent([1.0, 2.0], [3.0, 4.0], "baseline_unknown_mode"),
                 "an unknown mode", argument="baseline_mode")


# -----------------------------------------------------------------------------------------------
# compute_velocity_trajectory
# -----------------------------------------------------------------------------------------------

def test_compute_velocity_trajectory():
    velocity = compute_velocity_trajectory([1.0, 2.0, 4.0, 7.0])
    _assert_result_array(velocity, (3,), "velocity", fortran_order=False)
    # the one value: the steps stay in order
    assert np.array_equal(velocity, [1.0, 2.0, 3.0]), f"velocity: got {velocity}"


def test_compute_velocity_trajectory_of_a_single_point_is_empty():
    _assert_result_array(compute_velocity_trajectory([5.0]), (0,), "velocity", fortran_order=False)


def test_compute_velocity_trajectory_rejects_no_points():
    assert_error(lambda: compute_velocity_trajectory(np.empty(0)), "no points", ERR_EMPTY_INPUT,
                 argument="trajectory")


def test_compute_velocity_trajectory_rejects_nan():
    assert_error(lambda: compute_velocity_trajectory([1.0, np.nan, 3.0]), "a NaN point", ERR_NAN_INF,
                 argument="trajectory")


def test_compute_velocity_trajectory_rejects_a_matrix():
    # the binding's own rank check, before the library is called
    assert_error(lambda: compute_velocity_trajectory(np.ones((1, 4))), "a matrix", argument="trajectory")


# -----------------------------------------------------------------------------------------------
# compute_acceleration_from_velocity_trajectory
# -----------------------------------------------------------------------------------------------

def test_compute_acceleration_from_velocity_trajectory():
    acceleration = compute_acceleration_from_velocity_trajectory([1.0, 3.0, 0.0], 4)
    _assert_result_array(acceleration, (2,), "acceleration", fortran_order=False)
    # the one value: the steps stay in order
    assert np.array_equal(acceleration, [2.0, -3.0]), f"acceleration: got {acceleration}"


def test_compute_acceleration_from_velocity_trajectory_rejects_too_few_velocities():
    # `velocity` has n_timepoints - 1 steps; 3 against 6 timepoints used to be read past (#221)
    assert_error(lambda: compute_acceleration_from_velocity_trajectory([1.0, 2.0, 3.0], 6),
                 "too few velocity steps", ERR_SIZE_MISMATCH, argument="velocity")


def test_compute_acceleration_from_velocity_trajectory_rejects_nan():
    assert_error(lambda: compute_acceleration_from_velocity_trajectory([1.0, np.nan, 3.0], 4),
                 "a NaN velocity", ERR_NAN_INF, argument="velocity")


def test_compute_acceleration_from_velocity_trajectory_rejects_a_matrix():
    # the binding's own rank check, before the library is called
    assert_error(lambda: compute_acceleration_from_velocity_trajectory(np.ones((1, 3)), 4), "a matrix",
                 argument="velocity")


# -----------------------------------------------------------------------------------------------
# compute_velocity_trajectories
# -----------------------------------------------------------------------------------------------

def test_compute_velocity_trajectories():
    trajectories = _trajectories()
    velocity = compute_velocity_trajectories(trajectories)
    _assert_result_array(velocity, (N_TIMEPOINTS - 1, N_FACTORS, N_SAMPLES), "velocity", fortran_order=True)
    # the one value: time comes first in the result, last in the input; factor 2 of sample 3 is
    # 6 * [1, 4, 9, 16]
    assert np.array_equal(velocity[:, 1, 2], [18.0, 30.0, 42.0]), f"velocity[:, 1, 2]: got {velocity[:, 1, 2]}"


def test_compute_velocity_trajectories_rejects_no_samples():
    assert_error(lambda: compute_velocity_trajectories(np.empty((N_FACTORS, 0, N_TIMEPOINTS))),
                 "no samples", ERR_EMPTY_INPUT, argument="trajectories")


def test_compute_velocity_trajectories_rejects_nan():
    trajectories = _trajectories()
    trajectories[1, 2, 3] = np.nan
    assert_error(lambda: compute_velocity_trajectories(trajectories), "a NaN point", ERR_NAN_INF,
                 argument="trajectories")


def test_compute_velocity_trajectories_rejects_a_matrix():
    # the binding's own rank check, before the library is called
    assert_error(lambda: compute_velocity_trajectories(np.ones((2, 3))), "a matrix", argument="trajectories")


# -----------------------------------------------------------------------------------------------
# compute_acceleration_from_velocity
# -----------------------------------------------------------------------------------------------

def test_compute_acceleration_from_velocity():
    # (n_timepoints - 1, factors, samples)
    velocity = np.zeros((3, 1, 2))
    velocity[:, 0, 1] = [1.0, 3.0, 0.0]
    acceleration = compute_acceleration_from_velocity(velocity, 4)
    _assert_result_array(acceleration, (2, 1, 2), "acceleration", fortran_order=True)
    # the one value: time is the first axis
    assert np.array_equal(acceleration[:, 0, 1], [2.0, -3.0]), f"acceleration[:, 0, 1]: got {acceleration[:, 0, 1]}"


def test_compute_acceleration_from_velocity_rejects_too_few_velocities():
    # `velocity` has n_timepoints - 1 steps; 3 against 6 timepoints used to be read past (#221)
    assert_error(lambda: compute_acceleration_from_velocity(np.ones((3, 1, 2)), 6),
                 "too few velocity steps", ERR_SIZE_MISMATCH, argument="velocity")


def test_compute_acceleration_from_velocity_rejects_infinity():
    velocity = np.ones((3, 1, 2))
    velocity[2, 0, 1] = np.inf
    assert_error(lambda: compute_acceleration_from_velocity(velocity, 4), "an infinite velocity", ERR_NAN_INF,
                 argument="velocity")


def test_compute_acceleration_from_velocity_rejects_no_factors():
    assert_error(lambda: compute_acceleration_from_velocity(np.empty((3, 0, 2)), 4), "no factors", ERR_EMPTY_INPUT,
                 argument="velocity")


# -----------------------------------------------------------------------------------------------
# compute_velocity_acceleration_contributions
# -----------------------------------------------------------------------------------------------

def test_compute_velocity_acceleration_contributions():
    result = compute_velocity_acceleration_contributions(_trajectories(), "baseline_raw")
    _assert_keys(result, {"contrib_velocity", "velocity_contribution_series",
                          "contrib_acceleration", "acceleration_contribution_series"})
    _assert_result_array(result["contrib_velocity"], (N_FACTORS, N_FACTORS, N_SAMPLES), "contrib_velocity",
                         fortran_order=True)
    _assert_result_array(result["velocity_contribution_series"], (N_TIMEPOINTS, N_FACTORS, N_FACTORS, N_SAMPLES),
                         "velocity_contribution_series", fortran_order=True)
    _assert_result_array(result["contrib_acceleration"], (N_FACTORS, N_FACTORS, N_SAMPLES), "contrib_acceleration",
                         fortran_order=True)
    _assert_result_array(result["acceleration_contribution_series"],
                         (N_TIMEPOINTS, N_FACTORS, N_FACTORS, N_SAMPLES), "acceleration_contribution_series",
                         fortran_order=True)


def test_compute_velocity_acceleration_contributions_series_start_late():
    # the one value: the velocity series starts at the second time point, the acceleration
    # series at the third. x1 [1, 3, 6, 10] has the velocity [2, 3, 4] and the acceleration
    # [1, 1], x2 [1, 2, 2, 1] has [1, 0, -1] and [-1, -1]
    trajectories = np.array([[[1.0, 3.0, 6.0, 10.0]], [[1.0, 2.0, 2.0, 1.0]]])
    result = compute_velocity_acceleration_contributions(trajectories, "baseline_raw")
    assert np.array_equal(result["velocity_contribution_series"][:, 0, 1, 0], [0.0, 2.0, 0.0, -4.0]), \
        f"velocity_contribution_series: got {result['velocity_contribution_series'][:, 0, 1, 0]}"
    assert np.array_equal(result["acceleration_contribution_series"][:, 0, 1, 0], [0.0, 0.0, -1.0, -1.0]), \
        f"acceleration_contribution_series: got {result['acceleration_contribution_series'][:, 0, 1, 0]}"


def test_compute_velocity_acceleration_contributions_rejects_no_timepoints():
    assert_error(lambda: compute_velocity_acceleration_contributions(np.empty((N_FACTORS, N_SAMPLES, 0)),
                                                                     "baseline_raw"),
                 "no time points", ERR_EMPTY_INPUT, argument="trajectories")


def test_compute_velocity_acceleration_contributions_rejects_nan():
    trajectories = _trajectories()
    trajectories[0, 0, 0] = np.nan
    assert_error(lambda: compute_velocity_acceleration_contributions(trajectories, "baseline_mean"),
                 "a NaN point", ERR_NAN_INF, argument="trajectories")


def test_compute_velocity_acceleration_contributions_rejects_an_unknown_mode():
    assert_error(lambda: compute_velocity_acceleration_contributions(_trajectories(), "baseline"),
                 "an unknown mode", argument="baseline_mode")


if __name__ == '__main__':
    run_all_tests(globals().values())
