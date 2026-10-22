# The `tox_trajectory_normalization` suite: each published procedure can be called from R, returns
# the documented type and shape, and raises the documented error. Whether the numbers are right is
# the Fortran suite's job (test/mod_test_tox_trajectory_normalization.F90), as the coding guide's
# "Where to Test What" asks.
#
# A series that cannot be scaled (constant over time) is not an error: the procedures return its
# per-series `status`, so those tests check that the status reaches the caller, at the right index.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_FACTORS <- 2L
N_SAMPLES <- 3L
N_TIMEPOINTS <- 4L

# (factors, samples, timepoints), every series varying over time.
.trajectories <- function() {
  array(outer(outer(seq_len(N_FACTORS), seq_len(N_SAMPLES)), seq_len(N_TIMEPOINTS)),
        dim = c(N_FACTORS, N_SAMPLES, N_TIMEPOINTS))
}

.assert_names <- function(result, names) {
  assert_true(is.list(result), "expected a list")
  assert_true(identical(sort(names(result)), sort(names)), paste("unexpected names", toString(names(result))))
}

.assert_dims <- function(result, dims, name) {
  assert_true(identical(dim(result), as.integer(dims)),
              sprintf("%s: expected dim %s, got %s", name, toString(dims), toString(dim(result))))
}

# -----------------------------------------------------------------------------------------------
# normalize_variable_timeseries
# -----------------------------------------------------------------------------------------------

test_normalize_variable_timeseries <- function() {
  result <- normalize_variable_timeseries(c(3, 1, 2))
  .assert_names(result, c("v_norm", "status"))
  assert_true(is.double(result$v_norm) && is.null(dim(result$v_norm)) && length(result$v_norm) == 3L,
              "v_norm: expected a double vector of length 3")
  assert_true(is.integer(result$status) && length(result$status) == 1L, "status: expected an integer scalar")
  assert_true(result$status == ERR_OK, paste("status: expected ERR_OK, got", result$status))
  # the one value: the binding hands the points over in order, the maximum first
  assert_true(identical(result$v_norm, c(1, 0, 0.5)), paste("v_norm: got", toString(result$v_norm)))
}

test_normalize_variable_timeseries_reports_a_constant_series_in_status <- function() {
  # not raised: the status is returned next to the result
  result <- normalize_variable_timeseries(rep(7, 3))
  assert_true(result$status == ERR_DIVISION_BY_ZERO,
              paste("status: expected ERR_DIVISION_BY_ZERO, got", result$status))
}

test_normalize_variable_timeseries_rejects_an_empty_series <- function() {
  assert_error(normalize_variable_timeseries(numeric(0)), "no points", ERR_EMPTY_INPUT)
}

test_normalize_variable_timeseries_rejects_na <- function() {
  # R's missing value is a NaN to the library
  assert_error(normalize_variable_timeseries(c(1, NA, 3)), "a missing point", ERR_NAN_INF)
}

test_normalize_variable_timeseries_rejects_text <- function() {
  # the binding's own type check, before the library is called
  assert_error(normalize_variable_timeseries(c("1", "2")), "v must be numeric")
}

# -----------------------------------------------------------------------------------------------
# normalize_single_trajectory
# -----------------------------------------------------------------------------------------------

test_normalize_single_trajectory <- function() {
  # (timepoints, factors)
  trajectory <- cbind(c(0, 1, 2), c(2, 1, 0))
  result <- normalize_single_trajectory(trajectory)
  .assert_names(result, c("trajectory_norm", "status"))
  assert_true(is.matrix(result$trajectory_norm) && is.double(result$trajectory_norm),
              "trajectory_norm: expected a double matrix")
  .assert_dims(result$trajectory_norm, c(3L, 2L), "trajectory_norm")
  assert_true(is.integer(result$status) && is.null(dim(result$status)) && length(result$status) == 2L,
              "status: expected an integer vector, one per factor")
  assert_true(all(result$status == ERR_OK), paste("status: got", toString(result$status)))
  # the one value: time runs down the rows, so the first factor rises and the second falls
  assert_true(identical(result$trajectory_norm, cbind(c(0, 0.5, 1), c(1, 0.5, 0))),
              paste("trajectory_norm: got", toString(result$trajectory_norm)))
}

test_normalize_single_trajectory_reports_a_constant_factor_in_status <- function() {
  # the second factor is constant over time; only its status says so
  result <- normalize_single_trajectory(cbind(c(0, 1, 2), c(5, 5, 5)))
  assert_true(identical(result$status, as.integer(c(ERR_OK, ERR_DIVISION_BY_ZERO))),
              paste("status: got", toString(result$status)))
}

test_normalize_single_trajectory_rejects_no_factors <- function() {
  assert_error(normalize_single_trajectory(matrix(numeric(0), nrow = 3, ncol = 0)), "no factors", ERR_EMPTY_INPUT)
}

test_normalize_single_trajectory_rejects_infinity <- function() {
  trajectory <- matrix(1, nrow = 3, ncol = 2)
  trajectory[2, 2] <- Inf
  assert_error(normalize_single_trajectory(trajectory), "an infinite point", ERR_NAN_INF)
}

test_normalize_single_trajectory_rejects_a_vector <- function() {
  # the binding's own type check, before the library is called
  assert_error(normalize_single_trajectory(c(1, 2, 3)), "trajectory must be a matrix")
}

# -----------------------------------------------------------------------------------------------
# normalize_all_trajectories
# -----------------------------------------------------------------------------------------------

test_normalize_all_trajectories <- function() {
  trajectories <- .trajectories()
  trajectories[1, 2, ] <- c(3, 2, 1, 0)
  result <- normalize_all_trajectories(trajectories)
  .assert_names(result, c("trajectories_norm", "status"))
  assert_true(is.array(result$trajectories_norm) && is.double(result$trajectories_norm),
              "trajectories_norm: expected a double array")
  .assert_dims(result$trajectories_norm, c(N_FACTORS, N_SAMPLES, N_TIMEPOINTS), "trajectories_norm")
  assert_true(is.matrix(result$status) && is.integer(result$status), "status: expected an integer matrix")
  .assert_dims(result$status, c(N_FACTORS, N_SAMPLES), "status")
  assert_true(all(result$status == ERR_OK), paste("status: got", toString(result$status)))
  # the one value: time is the last axis, and this series falls over it
  assert_true(isTRUE(all.equal(result$trajectories_norm[1, 2, ], c(1, 2 / 3, 1 / 3, 0))),
              paste("trajectories_norm[1, 2, ]: got", toString(result$trajectories_norm[1, 2, ])))
}

test_normalize_all_trajectories_reports_a_constant_series_in_status <- function() {
  # one (factor, sample) constant over time: its status, and only its, says so
  trajectories <- .trajectories()
  trajectories[2, 3, ] <- 5
  status <- normalize_all_trajectories(trajectories)$status
  expected <- matrix(as.integer(ERR_OK), nrow = N_FACTORS, ncol = N_SAMPLES)
  expected[2, 3] <- as.integer(ERR_DIVISION_BY_ZERO)
  assert_true(identical(status, expected), paste("status: got", toString(status)))
}

test_normalize_all_trajectories_rejects_no_samples <- function() {
  assert_error(normalize_all_trajectories(array(numeric(0), dim = c(N_FACTORS, 0L, N_TIMEPOINTS))),
               "no samples", ERR_EMPTY_INPUT)
}

test_normalize_all_trajectories_rejects_nan <- function() {
  trajectories <- .trajectories()
  trajectories[2, 1, 3] <- NaN
  assert_error(normalize_all_trajectories(trajectories), "a NaN point", ERR_NAN_INF)
}

test_normalize_all_trajectories_rejects_a_matrix <- function() {
  # the binding's own rank check, before the library is called
  assert_error(normalize_all_trajectories(matrix(1, nrow = 2, ncol = 3)), "trajectories must have rank 3")
}

run_all_tests()
