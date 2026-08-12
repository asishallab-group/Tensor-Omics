

# Set library path and compile
source("r/load_tensor_omics.R")
source("r/test_helpers.R")


# Constants
TOL <- 1e-12

# =====================================================
# Test: compute_baselines_factor_dependent
# =====================================================

test_compute_baselines_factor_dependent <- function() {
  
  factor <- c(1.0, 3.0, 2.0, 4.0)
  dependent <- c(5.0, 7.0, 6.0, 8.0)
  
  # RAW baseline_mode = > zero baselines
  res_raw <- compute_baselines_factor_dependent(factor, dependent, baseline_mode = "raw")
  assert_true(abs(res_raw$factor_baseline - 0.0) < TOL)
  assert_true(abs(res_raw$dependent_baseline - 0.0) < TOL)
  
  # MIN baseline_mode = > min values
  res_min <- compute_baselines_factor_dependent(factor, dependent, baseline_mode = "min")
  assert_true(abs(res_min$factor_baseline - min(factor)) < TOL)
  assert_true(abs(res_min$dependent_baseline - min(dependent)) < TOL)
  
  # MEAN baseline_mode = > arithmetic mean
  res_mean <- compute_baselines_factor_dependent(factor, dependent, baseline_mode = "mean")
  assert_true(abs(res_mean$factor_baseline - mean(factor)) < TOL)
  assert_true(abs(res_mean$dependent_baseline - mean(dependent)) < TOL)
  
  # Mismatched lengths should raise error
  assert_error(compute_baselines_factor_dependent(factor, dependent[-length(dependent)], baseline_mode = "raw"), "Should have raised error for mismatched lengths")
  
  # Invalid mode should raise error
  assert_error(compute_baselines_factor_dependent(factor, dependent, baseline_mode = "unknown_mode"), "Should have raised error for invalid mode")
  
}

# =====================================================
# Test: compute_contributions
# =====================================================

test_compute_contributions <- function() {
  
  # Case 1: RAW baseline
  factor <- c(1.0, 2.0, 3.0, 4.0)
  dependent <- c(2.0, 1.0, 0.0, -1.0)
  result <- compute_contributions(factor, dependent, baseline_mode = "raw")
  
  expected_local <- factor * dependent
  expected_total <- sum(expected_local)
  assert_true(all(abs(result$local_contributions - expected_local) < TOL))
  assert_true(abs(result$total_contribution - expected_total) < TOL)
  
  # Case 2: MIN baseline
  factor <- c(3.0, 5.0, 2.0, 4.0)
  dependent <- c(1.0, 2.0, 0.0, -1.0)
  result <- compute_contributions(factor, dependent, baseline_mode = "min")
  
  expected_local <- (factor - min(factor)) * (dependent - min(dependent))
  expected_total <- sum(expected_local)
  assert_true(all(abs(result$local_contributions - expected_local) < TOL))
  assert_true(abs(result$total_contribution - expected_total) < TOL)
  
  # Case 3: MEAN baseline
  factor <- c(1.0, 2.0, 3.0, 4.0)
  dependent <- c(4.0, 3.0, 2.0, 1.0)
  result <- compute_contributions(factor, dependent, baseline_mode = "mean")
  
  expected_local <- (factor - mean(factor)) * (dependent - mean(dependent))
  expected_total <- sum(expected_local)
  assert_true(all(abs(result$local_contributions - expected_local) < TOL))
  assert_true(abs(result$total_contribution - expected_total) < TOL)
  
}

# =====================================================
# Test: compute_all_contributions
# =====================================================

test_compute_all_contributions <- function() {
  
  # Case 1: MEAN baseline with 2 factors and 1 sample
  trajectories <- array(c(1.0, 4.0, 2.0, 5.0, 3.0, 6.0), dim = c(2L, 1L, 3L))  # 3D array
  factor_indices <- c(1L)
  dependent_indices <- c(2L)
  
  result <- compute_all_contributions(trajectories, factor_indices, dependent_indices, baseline_mode = "mean")
  expected_local <- c(1.0, 0.0, 1.0)
  expected_total <- 2.0
  # Update indexing to match actual shape
  dims <- dim(result$local_contributions)
  if (length(dims) == 4) {
    assert_true(all(abs(result$local_contributions[,1,1,1] - expected_local) < TOL))
    assert_true(abs(result$total_contributions[1,1,1] - expected_total) < TOL)
  } else if (length(dims) == 3) {
    assert_true(all(abs(result$local_contributions[,1,1] - expected_local) < TOL))
    assert_true(abs(result$total_contributions[1,1] - expected_total) < TOL)
  } else if (length(dims) == 2) {
    assert_true(all(abs(result$local_contributions[,1] - expected_local) < TOL))
    assert_true(abs(result$total_contributions[1] - expected_total) < TOL)
  } else {
    stop("Unexpected number of dimensions in local_contributions: ", length(dims))
  }
  
  # Case 2: MIN baseline
  trajectories <- array(c(2.0, 5.0, 4.0, 3.0, 6.0, 5.0), dim = c(2L, 1L, 3L))  # 3D array
  result <- compute_all_contributions(trajectories, factor_indices, dependent_indices, baseline_mode = "min")
  expected_local <- c(0.0, 0.0, 8.0)
  expected_total <- 8.0
  assert_true(all(abs(result$local_contributions[,1,1,1] - expected_local) < TOL))
  assert_true(abs(result$total_contributions[1,1,1] - expected_total) < TOL)
  
}

# =====================================================
# Test: perform_permutation_test
# =====================================================

test_perform_permutation_test <- function() {
  
  n_factors <- 2L
  n_samples <- 3L
  n_timepoints <- 3L
  n_permutations <- 3L
  trajectories <- numeric(n_factors * n_samples * n_timepoints)
  trajectories[1:3] <- c(1.0, 2.0, 3.0)
  trajectories[4:6] <- c(2.0, 4.0, 6.0)
  trajectories[7:9] <- c(3.0, 6.0, 9.0)
  trajectories[10:12] <- c(4.0, 5.0, 6.0)
  trajectories[13:15] <- c(1.0, 3.0, 5.0)
  trajectories[16:18] <- c(2.0, 4.0, 6.0)
  dims <- c(n_factors, n_samples, n_timepoints)
  factor_idx <- 1L
  dependent_idx <- 2L
  sample_idx <- 1L
  mode <- "mean"  # MEAN baseline
  random_seed <- 12345L
  trajectories <- array(trajectories, dim = dims)
  result <- perform_permutation_test(
    trajectories, factor_idx, dependent_idx, sample_idx,
    baseline_mode = mode, 
    n_permutations = n_permutations, 
    random_seed = random_seed
  )
  local <- result$local_contributions
  total <- result$total_contributions
  assert_true(all(is.finite(local)))
  assert_true(all(is.finite(total)))
  assert_true(!all(abs(local) < TOL))
}

# =====================================================
# Test: compute_p_values
# =====================================================

test_compute_p_values <- function() {
  # Case 1: Valid inputs
  n_timepoints <- 3L
  n_permutations <- 4L
  local_obs <- c(2.0, 0.0, 2.0)
  total_obs <- 4.0
  local_perm <- matrix(c(
    1.0, 0.0, 1.0,  # permutation 1
    2.0, 0.0, 2.0,  # permutation 2
    3.0, 1.0, 3.0,  # permutation 3
    0.0, 0.0, 0.0   # permutation 4
  ), nrow = n_timepoints, ncol = n_permutations, byrow = FALSE)
  total_perm <- c(2.0, 4.0, 7.0, 0.0)
  dims <- c(n_timepoints, n_permutations)
  result <- compute_p_values(local_obs, total_obs, local_perm, total_perm)
  local_p <- result$local_p_values
  total_p <- result$total_p_value
  expected_local <- c(0.5, 1.0, 0.5)
  expected_total <- 0.5
  assert_true(all(abs(local_p - expected_local) < TOL))
  assert_true(abs(total_p - expected_total) < TOL)
  # Case 2: NaN in observed contributions
  local_obs_nan <- c(2.0, 0.0, NaN)
  assert_error(compute_p_values(local_obs_nan, total_obs, local_perm, total_perm), "Should have raised error for NaN input", ERR_NAN_INF)
  # Case 3: Inf in permutation contributions
  local_perm_inf <- local_perm
  local_perm_inf[3, 4] <- Inf
  assert_error(compute_p_values(local_obs, total_obs, local_perm_inf, total_perm), "Should have raised error for Inf input", ERR_NAN_INF)
}


# =====================================================
# Test: compute_velocity_trajectories
# =====================================================
test_compute_velocity_trajectories <- function() {

  # Input trajectories shape:
  # (n_factors=1, n_samples=2, n_timepoints=4)
  trajectories <- array(0.0, dim = c(1L, 2L, 4L))

  # factor 1, sample 1
  trajectories[1, 1, ] <- c(1.0, 2.0, 4.0, 7.0)

  # factor 1, sample 2
  trajectories[1, 2, ] <- c(0.0, -1.0, -1.0, 0.0)


  velocity <- compute_velocity_trajectories(trajectories)

  # Expected output shape:
  # (n_timepoints-1, n_factors, n_samples) = (3, 1, 2)
  expected_velocity <- array(0.0, dim = c(3L, 1L, 2L))

  # sample 1: diff(c(1,2,4,7)) = c(1,2,3)
  expected_velocity[, 1, 1] <- c(1.0, 2.0, 3.0)

  # sample 2: diff(c(0,-1,-1,0)) = c(-1,0,1)
  expected_velocity[, 1, 2] <- c(-1.0, 0.0, 1.0)

  assert_true(identical(dim(velocity), c(3L, 1L, 2L)))
  assert_true(all(is.finite(velocity)))
  assert_true(all(abs(velocity - expected_velocity) < TOL))

}

# =====================================================
# Test: compute_acceleration_from_velocity
# =====================================================
test_compute_acceleration_from_velocity <- function() {

  # Input velocity shape:
  # (n_timepoints-1, n_factors, n_samples) = (3, 1, 2)
  velocity <- array(0.0, dim = c(3L, 1L, 2L))

  # sample 1 velocity: c(1,2,3)
  velocity[, 1, 1] <- c(1.0, 2.0, 3.0)

  # sample 2 velocity: c(-1,0,1)
  velocity[, 1, 2] <- c(-1.0, 0.0, 1.0)


  acceleration <- compute_acceleration_from_velocity(velocity, dim(velocity)[1] + 1L)

  # Since n_timepoints = n_vel + 1 = 4,
  # expected output shape = (n_timepoints-2, n_factors, n_samples) = (2, 1, 2)
  expected_acceleration <- array(0.0, dim = c(2L, 1L, 2L))

  # sample 1: diff(c(1,2,3)) = c(1,1)
  expected_acceleration[, 1, 1] <- c(1.0, 1.0)

  # sample 2: diff(c(-1,0,1)) = c(1,1)
  expected_acceleration[, 1, 2] <- c(1.0, 1.0)

  assert_true(identical(dim(acceleration), c(2L, 1L, 2L)))
  assert_true(all(is.finite(acceleration)))
  assert_true(all(abs(acceleration - expected_acceleration) < TOL))

}

# =====================================================
# Test: compute_velocity_trajectory
# =====================================================
test_compute_velocity_trajectory <- function() {

  trajectory <- c(1.0, 2.0, 4.0, 7.0)
  velocity <- compute_velocity_trajectory(trajectory)

  # Expected raw output length = n_timepoints - 1 = 3
  expected_velocity <- c(1.0, 2.0, 3.0)

  assert_true(length(velocity) == 3L)
  assert_true(all(is.finite(velocity)))
  assert_true(all(abs(velocity - expected_velocity) < TOL))

}

# =====================================================
# Test: compute_acceleration_from_velocity_trajectory
# =====================================================
test_compute_acceleration_from_velocity_trajectory <- function() {

  velocity <- c(1.0, 2.0, 3.0)
  acceleration <- compute_acceleration_from_velocity_trajectory(velocity, length(velocity) + 1L)

  # Expected raw output length = n_timepoints - 2 = 2
  # Since original trajectory length would be 4
  expected_acceleration <- c(1.0, 1.0)

  assert_true(length(acceleration) == 2L)
  assert_true(all(is.finite(acceleration)))
  assert_true(all(abs(acceleration - expected_acceleration) < TOL))

}
# =====================================================
# Test: compute_velocity_acceleration_contributions
# =====================================================
test_compute_velocity_acceleration_contributions <- function() {

  # Input shape: (n_factors=2, n_samples=1, n_timepoints=4)
  trajectories <- array(0.0, dim = c(2L, 1L, 4L))

  # Factor 1, Sample 1
  trajectories[1, 1, ] <- c(1.0, 3.0, 6.0, 10.0)

  # Factor 2, Sample 1
  trajectories[2, 1, ] <- c(1.0, 2.0, 2.0, 1.0)

  mode <- "raw"

  # Defensive shape and value checks
  assert_true(length(dim(trajectories)) == 3)
  assert_true(all(is.finite(trajectories)))
  # Ensure Fortran order (column-major)
  trajectories <- aperm(trajectories, c(1,2,3))
  result <- compute_velocity_acceleration_contributions(trajectories, baseline_mode = mode)

  contrib_velocity <- result$contrib_velocity
  velocity_contribution_series <- result$velocity_contribution_series
  contrib_acceleration <- result$contrib_acceleration
  acceleration_contribution_series <- result$acceleration_contribution_series

  # Shape checks based on Fortran declarations
  assert_true(identical(dim(contrib_velocity), c(2L, 2L, 1L)))
  assert_true(identical(dim(contrib_acceleration), c(2L, 2L, 1L)))

  # Series arrays are expected from your wrapper as:
  # (n_timepoints, n_factors, n_factors, n_samples) = (4, 2, 2, 1)
  assert_true(identical(dim(velocity_contribution_series), c(4L, 2L, 2L, 1L)))
  assert_true(identical(dim(acceleration_contribution_series), c(4L, 2L, 2L, 1L)))

  # Finite checks
  assert_true(all(is.finite(contrib_velocity)))
  assert_true(all(is.finite(contrib_acceleration)))
  assert_true(all(is.finite(velocity_contribution_series)))
  assert_true(all(is.finite(acceleration_contribution_series)))

  # Optional numerical sanity check for one pair:
  # raw mode uses plain products of velocity / acceleration series
  #
  # trajectories:
  # factor 1: [1, 3, 6, 10] -> velocity [2, 3, 4] -> acceleration [1, 1]
  # factor 2: [1, 2, 2, 1]  -> velocity [1, 0,-1] -> acceleration [-1,-1]
  #
  # velocity pair contribution (1 -> 2): [2*1, 3*0, 4*(-1)] = [2, 0, -4], total = -2
  # acceleration pair contribution (1 -> 2): [1*(-1), 1*(-1)] = [-1, -1], total = -2
  #
  # If your backend stores time-aligned padded series of length 4:
  # velocity expected ~= [0, 2, 0, -4]
  # acceleration expected ~= [0, 0, -1, -1]

  expected_vel_total_12 <- -2.0
  expected_acc_total_12 <- -2.0

  assert_true(abs(contrib_velocity[1, 2, 1] - expected_vel_total_12) < TOL)
  assert_true(abs(contrib_acceleration[1, 2, 1] - expected_acc_total_12) < TOL)

  expected_vel_series_12 <- c(0.0, 2.0, 0.0, -4.0)
  expected_acc_series_12 <- c(0.0, 0.0, -1.0, -1.0)

  assert_true(all(abs(velocity_contribution_series[, 1, 2, 1] - expected_vel_series_12) < TOL))
  assert_true(all(abs(acceleration_contribution_series[, 1, 2, 1] - expected_acc_series_12) < TOL))

}

# =====================================================
# Test: compute_velocity_acceleration_contributions
# =====================================================
test_compute_velocity_acceleration_contributions_plain <- function() {

  # Same input as the expert test so outputs should match
  trajectories <- array(0.0, dim = c(2L, 1L, 4L))

  # Factor 1, Sample 1
  trajectories[1, 1, ] <- c(1.0, 3.0, 6.0, 10.0)

  # Factor 2, Sample 1
  trajectories[2, 1, ] <- c(1.0, 2.0, 2.0, 1.0)

  mode <- "raw"


  result <- compute_velocity_acceleration_contributions(trajectories, baseline_mode = mode)

  contrib_velocity <- result$contrib_velocity
  velocity_contribution_series <- result$velocity_contribution_series
  contrib_acceleration <- result$contrib_acceleration
  acceleration_contribution_series <- result$acceleration_contribution_series

  # Shape checks
  assert_true(identical(dim(contrib_velocity), c(2L, 2L, 1L)))
  assert_true(identical(dim(contrib_acceleration), c(2L, 2L, 1L)))
  assert_true(identical(dim(velocity_contribution_series), c(4L, 2L, 2L, 1L)))
  assert_true(identical(dim(acceleration_contribution_series), c(4L, 2L, 2L, 1L)))

  # Finite checks
  assert_true(all(is.finite(contrib_velocity)))
  assert_true(all(is.finite(contrib_acceleration)))
  assert_true(all(is.finite(velocity_contribution_series)))
  assert_true(all(is.finite(acceleration_contribution_series)))

  # Same numerical check as above
  expected_vel_total_12 <- -2.0
  expected_acc_total_12 <- -2.0

  assert_true(abs(contrib_velocity[1, 2, 1] - expected_vel_total_12) < TOL)
  assert_true(abs(contrib_acceleration[1, 2, 1] - expected_acc_total_12) < TOL)

  expected_vel_series_12 <- c(0.0, 2.0, 0.0, -4.0)
  expected_acc_series_12 <- c(0.0, 0.0, -1.0, -1.0)

  assert_true(all(abs(velocity_contribution_series[, 1, 2, 1] - expected_vel_series_12) < TOL))
  assert_true(all(abs(acceleration_contribution_series[, 1, 2, 1] - expected_acc_series_12) < TOL))

}

run_all_tests()
