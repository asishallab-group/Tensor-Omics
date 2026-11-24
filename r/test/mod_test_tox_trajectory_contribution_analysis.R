# The `tox_trajectory_contribution_analysis` suite: each published procedure can be called from R,
# returns the documented type and shape, and raises the documented error, naming the argument.
# Whether the numbers are right is the Fortran suite's job
# (test/mod_test_tox_trajectory_contribution_analysis/), as the coding guide's "Where to Test What"
# asks; the one value a test checks here is there to show that the binding hands the arrays over
# in the right order and the indices 1-based.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_FACTORS <- 2L
N_SAMPLES <- 3L
N_TIMEPOINTS <- 4L

# (factors, samples, timepoints), every series varying over time.
.trajectories <- function() {
  array(outer(outer(seq_len(N_FACTORS), seq_len(N_SAMPLES)), seq_len(N_TIMEPOINTS)^2),
        dim = c(N_FACTORS, N_SAMPLES, N_TIMEPOINTS))
}

# Factor 1 is [1, 0] in every sample, factor 2 of sample s is [s, 0]: a raw permutation's total is
# the number of the sample it drew.
.sample_marking <- function(n_samples) {
  trajectories <- array(0, dim = c(2L, n_samples, 2L))
  trajectories[1, , 1] <- 1
  trajectories[2, , 1] <- seq_len(n_samples)
  trajectories
}

.assert_names <- function(result, names) {
  assert_true(is.list(result), "expected a list")
  assert_true(identical(sort(names(result)), sort(names)), paste("unexpected names", toString(names(result))))
}

.assert_double_array <- function(result, dims, name) {
  assert_true(is.double(result), sprintf("%s: expected doubles", name))
  actual <- if (is.null(dim(result))) length(result) else dim(result)
  assert_true(identical(as.integer(actual), as.integer(dims)),
              sprintf("%s: expected dim %s, got %s", name, toString(dims), toString(actual)))
}

.assert_double_scalar <- function(result, name) {
  assert_true(is.double(result) && length(result) == 1L && is.null(dim(result)),
              sprintf("%s: expected a double scalar", name))
}

# The expected error -- `code` NULL for one the binding raises itself, before the library is
# called -- naming `argument`.
.assert_refused <- function(expr, msg, code = NULL, argument) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  if (is.null(err)) stop(msg, ": nothing was raised", call. = FALSE)
  assert_error(stop(err), msg, code)
  assert_true(identical(err$argument, argument),
              sprintf("%s: expected '%s' named, got '%s'", msg, argument, toString(err$argument)))
}

# -----------------------------------------------------------------------------------------------
# perform_permutation_test
# -----------------------------------------------------------------------------------------------

test_perform_permutation_test <- function() {
  result <- perform_permutation_test(.trajectories(), 1L, 2L, 2L, "baseline_mean", 5L, random_seed = 42L)
  .assert_names(result, c("local_contributions", "total_contributions"))
  assert_true(is.matrix(result$local_contributions), "local_contributions: expected a matrix")
  .assert_double_array(result$local_contributions, c(N_TIMEPOINTS, 5L), "local_contributions")
  .assert_double_array(result$total_contributions, 5L, "total_contributions")
}

test_perform_permutation_test_draws_the_other_sample <- function() {
  # the one value: of two samples, sample 1 can only draw sample 2
  result <- perform_permutation_test(.sample_marking(2L), 1L, 2L, 1L, "baseline_raw", 4L, random_seed = 42L)
  assert_true(identical(result$total_contributions, rep(2, 4)),
              paste("total_contributions: got", toString(result$total_contributions)))
}

test_perform_permutation_test_random_seed_reaches_the_library <- function() {
  # the same seed twice draws the same: the binding hands it over, it does not drop it
  trajectories <- .sample_marking(5L)
  first <- perform_permutation_test(trajectories, 1L, 2L, 3L, "baseline_raw", 64L, random_seed = 7L)
  second <- perform_permutation_test(trajectories, 1L, 2L, 3L, "baseline_raw", 64L, random_seed = 7L)
  assert_true(identical(first$total_contributions, second$total_contributions), "the same seed drew otherwise")
}

test_perform_permutation_test_without_a_seed <- function() {
  result <- perform_permutation_test(.sample_marking(3L), 1L, 2L, 1L, "baseline_raw", 3L)
  assert_true(all(result$total_contributions %in% c(2, 3)),
              paste("total_contributions: got", toString(result$total_contributions)))
}

test_perform_permutation_test_rejects_a_sample_index_out_of_range <- function() {
  .assert_refused(perform_permutation_test(.trajectories(), 1L, 2L, N_SAMPLES + 1L, "baseline_raw", 2L,
                                           random_seed = 42L),
                  "sample_idx past n_samples", ERR_INVALID_INPUT, "sample_idx")
}

test_perform_permutation_test_rejects_a_single_sample <- function() {
  # n_samples has the minimum 2: a permutation draws another sample; the extent is read off
  # `trajectories`, so that is the argument named
  .assert_refused(perform_permutation_test(.sample_marking(1L), 1L, 2L, 1L, "baseline_raw", 2L, random_seed = 42L),
                  "a single sample", ERR_INVALID_INPUT, "trajectories")
}

test_perform_permutation_test_rejects_zero_permutations <- function() {
  .assert_refused(perform_permutation_test(.trajectories(), 1L, 2L, 1L, "baseline_raw", 0L, random_seed = 42L),
                  "no permutations", ERR_EMPTY_INPUT, "n_permutations")
}

test_perform_permutation_test_rejects_a_negative_permutation_count <- function() {
  # the binding's own extent check, before it allocates by it
  .assert_refused(perform_permutation_test(.trajectories(), 1L, 2L, 1L, "baseline_raw", -1L, random_seed = 42L),
                  "a negative n_permutations", argument = "n_permutations")
}

test_perform_permutation_test_rejects_na <- function() {
  # R's missing value is a NaN to the library
  trajectories <- .trajectories()
  trajectories[2, 3, 4] <- NA
  .assert_refused(perform_permutation_test(trajectories, 1L, 2L, 1L, "baseline_raw", 2L, random_seed = 42L),
                  "a missing point", ERR_NAN_INF, "trajectories")
}

test_perform_permutation_test_rejects_an_unknown_mode <- function() {
  # the binding's own check, before the library is called
  .assert_refused(perform_permutation_test(.trajectories(), 1L, 2L, 1L, "baseline_median", 2L, random_seed = 42L),
                  "an unknown mode", argument = "baseline_mode")
}

# -----------------------------------------------------------------------------------------------
# compute_p_values
# -----------------------------------------------------------------------------------------------

test_compute_p_values <- function() {
  # (timepoints, permutations): the columns are the permutations
  local_perm <- rbind(c(0, 1, 2, 3), c(5, 5, 5, 5))
  result <- compute_p_values(c(2, 1), 3, local_perm, c(1, 2, 3, 4))
  .assert_names(result, c("local_p_values", "total_p_value"))
  .assert_double_array(result$local_p_values, 2L, "local_p_values")
  .assert_double_scalar(result$total_p_value, "total_p_value")
  # the one value: point 1 is beaten by permutations 3 and 4, point 2 by all, the total by 3 and 4
  assert_true(identical(result$local_p_values, c(0.5, 1)),
              paste("local_p_values: got", toString(result$local_p_values)))
  assert_true(identical(result$total_p_value, 0.5), paste("total_p_value: got", result$total_p_value))
}

test_compute_p_values_rejects_nan <- function() {
  .assert_refused(compute_p_values(c(2, NaN), 3, matrix(1, 2, 4), rep(1, 4)),
                  "a NaN observation", ERR_NAN_INF, "local_contributions_observed")
}

test_compute_p_values_rejects_infinity_in_the_permutations <- function() {
  local_perm <- matrix(1, 2, 4)
  local_perm[2, 4] <- Inf
  .assert_refused(compute_p_values(c(2, 1), 3, local_perm, rep(1, 4)),
                  "an infinite permutation", ERR_NAN_INF, "local_contributions_perm_test")
}

test_compute_p_values_rejects_mismatched_permutations <- function() {
  # the binding's own extent checks, before the library is called
  .assert_refused(compute_p_values(c(2, 1), 3, matrix(1, 3, 4), rep(1, 4)),
                  "a permutation of three points against two", argument = "local_contributions_perm_test")
  .assert_refused(compute_p_values(c(2, 1), 3, matrix(1, 2, 4), rep(1, 3)),
                  "three totals for four permutations", argument = "total_contributions_perm_test")
}

# -----------------------------------------------------------------------------------------------
# compute_contributions
# -----------------------------------------------------------------------------------------------

test_compute_contributions <- function() {
  result <- compute_contributions(c(1, 2, 3), c(3, 1, 2), "baseline_raw")
  .assert_names(result, c("local_contributions", "total_contribution"))
  .assert_double_array(result$local_contributions, 3L, "local_contributions")
  .assert_double_scalar(result$total_contribution, "total_contribution")
  # the one value: the points stay in order
  assert_true(identical(result$local_contributions, c(3, 2, 6)),
              paste("local_contributions: got", toString(result$local_contributions)))
}

test_compute_contributions_rejects_no_points <- function() {
  .assert_refused(compute_contributions(numeric(0), numeric(0), "baseline_raw"),
                  "no points", ERR_EMPTY_INPUT, "factor")
}

test_compute_contributions_rejects_na <- function() {
  .assert_refused(compute_contributions(c(1, 2), c(NA, 1), "baseline_raw"),
                  "a missing point", ERR_NAN_INF, "dependent")
}

test_compute_contributions_rejects_mismatched_series <- function() {
  # the binding's own extent check, before the library is called
  .assert_refused(compute_contributions(c(1, 2, 3), c(1, 2), "baseline_raw"),
                  "a dependent shorter than the factor", argument = "dependent")
}

test_compute_contributions_rejects_an_unknown_mode <- function() {
  .assert_refused(compute_contributions(c(1, 2), c(3, 4), "mean"), "an unknown mode", argument = "baseline_mode")
}

# -----------------------------------------------------------------------------------------------
# compute_all_contributions
# -----------------------------------------------------------------------------------------------

test_compute_all_contributions <- function() {
  result <- compute_all_contributions(.trajectories(), c(2L, 1L, 2L), c(1L, 2L), "baseline_min")
  .assert_names(result, c("local_contributions", "total_contributions"))
  .assert_double_array(result$local_contributions, c(N_TIMEPOINTS, 3L, 2L, N_SAMPLES), "local_contributions")
  .assert_double_array(result$total_contributions, c(3L, 2L, N_SAMPLES), "total_contributions")
}

test_compute_all_contributions_takes_indices_1_based <- function() {
  # the one value: factor 1 against factor 2, raw, in the one sample: [1*3, 2*4] = [3, 8]
  trajectories <- array(c(1, 3, 2, 4), dim = c(2L, 1L, 2L))
  result <- compute_all_contributions(trajectories, 1L, 2L, "baseline_raw")
  assert_true(identical(as.vector(result$local_contributions[, 1, 1, 1]), c(3, 8)),
              paste("local_contributions: got", toString(result$local_contributions)))
}

test_compute_all_contributions_rejects_an_index_out_of_range <- function() {
  .assert_refused(compute_all_contributions(.trajectories(), c(1L, N_FACTORS + 1L), 1L, "baseline_raw"),
                  "a factor index past n_factors", ERR_INVALID_INPUT, "factor_indices")
  .assert_refused(compute_all_contributions(.trajectories(), 1L, 0L, "baseline_raw"),
                  "a dependent index 0", ERR_INVALID_INPUT, "dependent_indices")
}

test_compute_all_contributions_rejects_an_empty_selection <- function() {
  .assert_refused(compute_all_contributions(.trajectories(), integer(0), 1L, "baseline_raw"),
                  "no factors selected", ERR_EMPTY_INPUT, "factor_indices")
}

test_compute_all_contributions_rejects_infinity <- function() {
  trajectories <- .trajectories()
  trajectories[1, 2, 3] <- -Inf
  .assert_refused(compute_all_contributions(trajectories, 1L, 2L, "baseline_raw"),
                  "an infinite point", ERR_NAN_INF, "trajectories")
}

# -----------------------------------------------------------------------------------------------
# compute_baselines_factor_dependent
# -----------------------------------------------------------------------------------------------

test_compute_baselines_factor_dependent <- function() {
  result <- compute_baselines_factor_dependent(c(3, -1, 4, 2), c(5, 9, 2, 6), "baseline_mean")
  .assert_names(result, c("factor_baseline", "dependent_baseline"))
  .assert_double_scalar(result$factor_baseline, "factor_baseline")
  .assert_double_scalar(result$dependent_baseline, "dependent_baseline")
  # the one value: each series has its own baseline
  assert_true(identical(c(result$factor_baseline, result$dependent_baseline), c(2, 5.5)),
              paste("got", result$factor_baseline, result$dependent_baseline))
}

test_compute_baselines_factor_dependent_rejects_no_points <- function() {
  .assert_refused(compute_baselines_factor_dependent(numeric(0), numeric(0), "baseline_min"),
                  "no points", ERR_EMPTY_INPUT, "factor")
}

test_compute_baselines_factor_dependent_rejects_na <- function() {
  .assert_refused(compute_baselines_factor_dependent(c(1, NA), c(1, 2), "baseline_raw"),
                  "a missing point", ERR_NAN_INF, "factor")
}

test_compute_baselines_factor_dependent_rejects_mismatched_series <- function() {
  # the binding's own extent check, before the library is called
  .assert_refused(compute_baselines_factor_dependent(c(1, 2, 3), c(1, 2), "baseline_raw"),
                  "a dependent shorter than the factor", argument = "dependent")
}

test_compute_baselines_factor_dependent_rejects_an_unknown_mode <- function() {
  .assert_refused(compute_baselines_factor_dependent(c(1, 2), c(3, 4), "baseline_unknown_mode"),
                  "an unknown mode", argument = "baseline_mode")
}

# -----------------------------------------------------------------------------------------------
# compute_velocity_trajectory
# -----------------------------------------------------------------------------------------------

test_compute_velocity_trajectory <- function() {
  velocity <- compute_velocity_trajectory(c(1, 2, 4, 7))
  .assert_double_array(velocity, 3L, "velocity")
  # the one value: the steps stay in order
  assert_true(identical(velocity, c(1, 2, 3)), paste("velocity: got", toString(velocity)))
}

test_compute_velocity_trajectory_of_a_single_point_is_empty <- function() {
  .assert_double_array(compute_velocity_trajectory(5), 0L, "velocity")
}

test_compute_velocity_trajectory_rejects_no_points <- function() {
  .assert_refused(compute_velocity_trajectory(numeric(0)), "no points", ERR_EMPTY_INPUT, "trajectory")
}

test_compute_velocity_trajectory_rejects_na <- function() {
  .assert_refused(compute_velocity_trajectory(c(1, NA, 3)), "a missing point", ERR_NAN_INF, "trajectory")
}

test_compute_velocity_trajectory_rejects_text <- function() {
  # the binding's own type check, before the library is called
  .assert_refused(compute_velocity_trajectory(c("1", "2")), "text", argument = "trajectory")
}

# -----------------------------------------------------------------------------------------------
# compute_acceleration_from_velocity_trajectory
# -----------------------------------------------------------------------------------------------

test_compute_acceleration_from_velocity_trajectory <- function() {
  acceleration <- compute_acceleration_from_velocity_trajectory(c(1, 3, 0), 4L)
  .assert_double_array(acceleration, 2L, "acceleration")
  # the one value: the steps stay in order
  assert_true(identical(acceleration, c(2, -3)), paste("acceleration: got", toString(acceleration)))
}

test_compute_acceleration_from_velocity_trajectory_rejects_too_few_velocities <- function() {
  # `velocity` has n_timepoints - 1 steps; 3 against 6 timepoints used to be read past (#221)
  .assert_refused(compute_acceleration_from_velocity_trajectory(c(1, 2, 3), 6L),
                  "too few velocity steps", ERR_SIZE_MISMATCH, "velocity")
}

test_compute_acceleration_from_velocity_trajectory_rejects_na <- function() {
  .assert_refused(compute_acceleration_from_velocity_trajectory(c(1, NA, 3), 4L),
                  "a missing velocity", ERR_NAN_INF, "velocity")
}

test_a_computed_axis_cannot_overflow_into_an_allocation <- function() {
  # acceleration has max(0, n_timepoints - 2) elements. In C int arithmetic -2147483647 - 2
  # overflowed (undefined, in practice 2147483647: a 16 GB vector); the shim now computes every
  # axis in double, so it is 0, as Fortran sizes it, and nothing is allocated (#215)
  r <- compute_acceleration_from_velocity_trajectory(numeric(0), -2147483647L)
  assert_true(is.double(r) && length(r) == 0L, "an empty acceleration, not an allocation")
}

# -----------------------------------------------------------------------------------------------
# compute_velocity_trajectories
# -----------------------------------------------------------------------------------------------

test_compute_velocity_trajectories <- function() {
  velocity <- compute_velocity_trajectories(.trajectories())
  .assert_double_array(velocity, c(N_TIMEPOINTS - 1L, N_FACTORS, N_SAMPLES), "velocity")
  # the one value: time comes first in the result, last in the input; factor 2 of sample 3 is
  # 6 * [1, 4, 9, 16]
  assert_true(identical(as.vector(velocity[, 2, 3]), c(18, 30, 42)),
              paste("velocity[, 2, 3]: got", toString(velocity[, 2, 3])))
}

test_compute_velocity_trajectories_rejects_no_samples <- function() {
  .assert_refused(compute_velocity_trajectories(array(numeric(0), dim = c(N_FACTORS, 0L, N_TIMEPOINTS))),
                  "no samples", ERR_EMPTY_INPUT, "trajectories")
}

test_compute_velocity_trajectories_rejects_na <- function() {
  trajectories <- .trajectories()
  trajectories[2, 3, 4] <- NA
  .assert_refused(compute_velocity_trajectories(trajectories), "a missing point", ERR_NAN_INF, "trajectories")
}

test_compute_velocity_trajectories_rejects_a_matrix <- function() {
  # the binding's own rank check, before the library is called
  .assert_refused(compute_velocity_trajectories(matrix(1, 2, 3)), "a matrix", argument = "trajectories")
}

# -----------------------------------------------------------------------------------------------
# compute_acceleration_from_velocity
# -----------------------------------------------------------------------------------------------

test_compute_acceleration_from_velocity <- function() {
  # (n_timepoints - 1, factors, samples)
  velocity <- array(0, dim = c(3L, 1L, 2L))
  velocity[, 1, 2] <- c(1, 3, 0)
  acceleration <- compute_acceleration_from_velocity(velocity, 4L)
  .assert_double_array(acceleration, c(2L, 1L, 2L), "acceleration")
  # the one value: time is the first axis
  assert_true(identical(as.vector(acceleration[, 1, 2]), c(2, -3)),
              paste("acceleration[, 1, 2]: got", toString(acceleration[, 1, 2])))
}

test_compute_acceleration_from_velocity_rejects_too_few_velocities <- function() {
  # `velocity` has n_timepoints - 1 steps; 3 against 6 timepoints used to be read past (#221)
  .assert_refused(compute_acceleration_from_velocity(array(1, dim = c(3L, 1L, 2L)), 6L),
                  "too few velocity steps", ERR_SIZE_MISMATCH, "velocity")
}

test_compute_acceleration_from_velocity_rejects_infinity <- function() {
  velocity <- array(1, dim = c(3L, 1L, 2L))
  velocity[3, 1, 2] <- Inf
  .assert_refused(compute_acceleration_from_velocity(velocity, 4L), "an infinite velocity", ERR_NAN_INF, "velocity")
}

test_compute_acceleration_from_velocity_rejects_no_factors <- function() {
  .assert_refused(compute_acceleration_from_velocity(array(numeric(0), dim = c(3L, 0L, 2L)), 4L),
                  "no factors", ERR_EMPTY_INPUT, "velocity")
}

# -----------------------------------------------------------------------------------------------
# compute_velocity_acceleration_contributions
# -----------------------------------------------------------------------------------------------

test_compute_velocity_acceleration_contributions <- function() {
  result <- compute_velocity_acceleration_contributions(.trajectories(), "baseline_raw")
  .assert_names(result, c("contrib_velocity", "velocity_contribution_series",
                          "contrib_acceleration", "acceleration_contribution_series"))
  .assert_double_array(result$contrib_velocity, c(N_FACTORS, N_FACTORS, N_SAMPLES), "contrib_velocity")
  .assert_double_array(result$velocity_contribution_series, c(N_TIMEPOINTS, N_FACTORS, N_FACTORS, N_SAMPLES),
                       "velocity_contribution_series")
  .assert_double_array(result$contrib_acceleration, c(N_FACTORS, N_FACTORS, N_SAMPLES), "contrib_acceleration")
  .assert_double_array(result$acceleration_contribution_series, c(N_TIMEPOINTS, N_FACTORS, N_FACTORS, N_SAMPLES),
                       "acceleration_contribution_series")
}

test_compute_velocity_acceleration_contributions_series_start_late <- function() {
  # the one value: the velocity series starts at the second time point, the acceleration series
  # at the third. x1 [1, 3, 6, 10] has the velocity [2, 3, 4] and the acceleration [1, 1], x2
  # [1, 2, 2, 1] has [1, 0, -1] and [-1, -1]
  trajectories <- array(c(1, 1, 3, 2, 6, 2, 10, 1), dim = c(2L, 1L, 4L))
  result <- compute_velocity_acceleration_contributions(trajectories, "baseline_raw")
  assert_true(identical(as.vector(result$velocity_contribution_series[, 1, 2, 1]), c(0, 2, 0, -4)),
              paste("velocity_contribution_series: got", toString(result$velocity_contribution_series[, 1, 2, 1])))
  assert_true(identical(as.vector(result$acceleration_contribution_series[, 1, 2, 1]), c(0, 0, -1, -1)),
              paste("acceleration_contribution_series: got",
                    toString(result$acceleration_contribution_series[, 1, 2, 1])))
}

test_compute_velocity_acceleration_contributions_rejects_no_timepoints <- function() {
  .assert_refused(compute_velocity_acceleration_contributions(array(numeric(0), dim = c(N_FACTORS, N_SAMPLES, 0L)),
                                                              "baseline_raw"),
                  "no time points", ERR_EMPTY_INPUT, "trajectories")
}

test_compute_velocity_acceleration_contributions_rejects_na <- function() {
  trajectories <- .trajectories()
  trajectories[1, 1, 1] <- NA
  .assert_refused(compute_velocity_acceleration_contributions(trajectories, "baseline_mean"),
                  "a missing point", ERR_NAN_INF, "trajectories")
}

test_compute_velocity_acceleration_contributions_rejects_an_unknown_mode <- function() {
  .assert_refused(compute_velocity_acceleration_contributions(.trajectories(), "baseline"),
                  "an unknown mode", argument = "baseline_mode")
}

run_all_tests()
