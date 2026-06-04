source("rcpp/tensoromics_functions.R")

cat("=== Testing trajectory normalization R wrappers ===\n")

TOL <- 1e-10

# -----------------------------------------------------
# tox_normalize_variable_timeseries
# -----------------------------------------------------
test_normalize_variable_timeseries <- function() {
  x <- c(1, 2, 3, 4, 5)
  res <- tox_normalize_variable_timeseries(x)
  expected <- c(0, 0.25, 0.5, 0.75, 1.0)

  stopifnot(length(res$v_norm) == length(x))
  stopifnot(all(abs(res$v_norm - expected) < TOL))
  cat("test_normalize_variable_timeseries passed ✓\n")
}

# -----------------------------------------------------
# tox_normalize_single_trajectory
# -----------------------------------------------------
test_normalize_single_trajectory <- function() {
  trajectory <- matrix(c(
    10, 20,
    11, 21,
    12, 22,
    13, 23
  ), nrow = 4, byrow = TRUE)

  res <- tox_normalize_single_trajectory(trajectory)

  expected_col <- c(0, 1/3, 2/3, 1)
  stopifnot(all(dim(res$traj_norm) == dim(trajectory)))
  stopifnot(all(abs(res$traj_norm[, 1] - expected_col) < TOL))
  stopifnot(all(abs(res$traj_norm[, 2] - expected_col) < TOL))
  cat("test_normalize_single_trajectory passed ✓\n")
}

# -----------------------------------------------------
# tox_normalize_all_trajectories
# -----------------------------------------------------
test_normalize_all_trajectories <- function() {
  trajectories <- array(0, dim = c(2, 2, 4))
  trajectories[1, 1, ] <- c(1, 2, 3, 4)
  trajectories[2, 1, ] <- c(10, 20, 30, 40)
  trajectories[1, 2, ] <- c(5, 6, 7, 8)
  trajectories[2, 2, ] <- c(2, 4, 6, 8)

  res <- tox_normalize_all_trajectories(trajectories)

  stopifnot(all(dim(res$traj_norm) == dim(trajectories)))

  for (i_factor in 1:dim(res$traj_norm)[1]) {
    for (i_sample in 1:dim(res$traj_norm)[2]) {
      ts <- res$traj_norm[i_factor, i_sample, ]
      stopifnot(abs(min(ts) - 0) < TOL)
      stopifnot(abs(max(ts) - 1) < TOL)
    }
  }

  cat("test_normalize_all_trajectories passed ✓\n")
}

cat("\nRunning trajectory normalization tests...\n")
test_normalize_variable_timeseries()
test_normalize_single_trajectory()
test_normalize_all_trajectories()
cat("All trajectory normalization R tests passed.\n")
