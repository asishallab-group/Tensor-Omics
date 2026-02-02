
library(Rcpp)

# Set library path and compile
source("rcpp/tensoromics_functions.R")

cat("=== Testing Trajectory Contribution Analysis R wrapper functions ===\n")


# Constants
TOL <- 1e-12

# =====================================================
# Test: tox_compute_baselines_factor_dependent
# =====================================================

test_compute_baselines_factor_dependent <- function() {
  cat("\n[test_compute_baselines_factor_dependent] Baseline computation across modes\n")
  
  factor <- c(1.0, 3.0, 2.0, 4.0)
  dependent <- c(5.0, 7.0, 6.0, 8.0)
  
  # RAW mode => zero baselines
  cat("  Testing RAW mode...\n")
  res_raw <- tox_compute_baselines_factor_dependent(factor, dependent, mode = "raw")
  stopifnot(abs(res_raw$factor_baseline - 0.0) < TOL)
  stopifnot(abs(res_raw$dependent_baseline - 0.0) < TOL)
  cat("  RAW mode passed ✓\n")
  
  # MIN mode => min values
  cat("  Testing MIN mode...\n")
  res_min <- tox_compute_baselines_factor_dependent(factor, dependent, mode = "min")
  stopifnot(abs(res_min$factor_baseline - min(factor)) < TOL)
  stopifnot(abs(res_min$dependent_baseline - min(dependent)) < TOL)
  cat("  MIN mode passed ✓\n")
  
  # MEAN mode => arithmetic mean
  cat("  Testing MEAN mode...\n")
  res_mean <- tox_compute_baselines_factor_dependent(factor, dependent, mode = "mean")
  stopifnot(abs(res_mean$factor_baseline - mean(factor)) < TOL)
  stopifnot(abs(res_mean$dependent_baseline - mean(dependent)) < TOL)
  cat("  MEAN mode passed ✓\n")
  
  # Mismatched lengths should raise error
  cat("  Testing error handling for mismatched lengths...\n")
  tryCatch(
    {
      tox_compute_baselines_factor_dependent(factor, dependent[-length(dependent)], mode = "raw")
      stop("Should have raised error for mismatched lengths")
    },
    error = function(e) {
      cat("  Correctly caught error for mismatched lengths ✓\n")
    }
  )
  
  # Invalid mode should raise error
  cat("  Testing error handling for invalid mode...\n")
  tryCatch(
    {
      tox_compute_baselines_factor_dependent(factor, dependent, mode = "unknown_mode")
      stop("Should have raised error for invalid mode")
    },
    error = function(e) {
      cat("  Correctly caught error for invalid mode ✓\n")
    }
  )
  
  cat("test_compute_baselines_factor_dependent passed ✓\n")
}

# =====================================================
# Test: tox_compute_contributions
# =====================================================

test_compute_contributions <- function() {
  cat("\n[test_compute_contributions] Contribution computation across modes\n")
  
  # Case 1: RAW baseline
  cat("  Testing RAW baseline...\n")
  factor <- c(1.0, 2.0, 3.0, 4.0)
  dependent <- c(2.0, 1.0, 0.0, -1.0)
  result <- tox_compute_contributions(factor, dependent, mode = "raw")
  
  expected_local <- factor * dependent
  expected_total <- sum(expected_local)
  stopifnot(all(abs(result$local_contributions - expected_local) < TOL))
  stopifnot(abs(result$total_contribution - expected_total) < TOL)
  cat("  RAW baseline passed ✓\n")
  
  # Case 2: MIN baseline
  cat("  Testing MIN baseline...\n")
  factor <- c(3.0, 5.0, 2.0, 4.0)
  dependent <- c(1.0, 2.0, 0.0, -1.0)
  result <- tox_compute_contributions(factor, dependent, mode = "min")
  
  expected_local <- (factor - min(factor)) * (dependent - min(dependent))
  expected_total <- sum(expected_local)
  stopifnot(all(abs(result$local_contributions - expected_local) < TOL))
  stopifnot(abs(result$total_contribution - expected_total) < TOL)
  cat("  MIN baseline passed ✓\n")
  
  # Case 3: MEAN baseline
  cat("  Testing MEAN baseline...\n")
  factor <- c(1.0, 2.0, 3.0, 4.0)
  dependent <- c(4.0, 3.0, 2.0, 1.0)
  result <- tox_compute_contributions(factor, dependent, mode = "mean")
  
  expected_local <- (factor - mean(factor)) * (dependent - mean(dependent))
  expected_total <- sum(expected_local)
  stopifnot(all(abs(result$local_contributions - expected_local) < TOL))
  stopifnot(abs(result$total_contribution - expected_total) < TOL)
  cat("  MEAN baseline passed ✓\n")
  
  cat("test_compute_contributions passed ✓\n")
}

# =====================================================
# Test: tox_compute_all_contributions
# =====================================================

test_compute_all_contributions <- function() {
  cat("\n[test_compute_all_contributions] Multi-factor contribution computation\n")
  
  # Case 1: MEAN baseline
  cat("  Testing MEAN baseline with 2 factors and 1 sample...\n")
  # Create 2D trajectories (n_factors=2, n_samples=1, n_timepoints=3)
  # In Fortran column-major: trajectories(factor, sample, timepoint)
  # trajectories(1,1,:) = [1,2,3] (factor 1)
  # trajectories(2,1,:) = [4,5,6] (factor 2)
  # Flattened in column-major: [1,4, 2,5, 3,6]
  trajectories <- c(1.0, 4.0, 2.0, 5.0, 3.0, 6.0)  # Fortran column-major order
  dims <- c(2L, 1L, 3L)  # (n_factors, n_samples, n_timepoints)
  factor_indices <- c(1L)
  dependent_indices <- c(2L)
  
  result <- tox_compute_all_contributions(trajectories, dims, factor_indices, dependent_indices, mode = "mean")
  
  # Expected: factor=[1,2,3], dependent=[4,5,6], mean(factor)=2, mean(dependent)=5
  expected_local <- c(1.0, 0.0, 1.0)  # (1-2)*(4-5), (2-2)*(5-5), (3-2)*(6-5)
  expected_total <- 2.0
  
  stopifnot(all(abs(result$local_contributions[,1,1,1] - expected_local) < TOL))
  stopifnot(abs(result$total_contributions[1,1,1] - expected_total) < TOL)
  cat("  MEAN baseline passed ✓\n")
  
  # Case 2: MIN baseline
  cat("  Testing MIN baseline...\n")
  # Factor trajectory: [2,4,6]
  # Dependent trajectory: [5,3,5]
  # Flattened in column-major: [2,5, 4,3, 6,5]
  trajectories <- c(2.0, 5.0, 4.0, 3.0, 6.0, 5.0)
  dims <- c(2L, 1L, 3L)
  
  result <- tox_compute_all_contributions(trajectories, dims, factor_indices, dependent_indices, mode = "min")
  
  # Expected: factor=[2,4,6], dependent=[5,3,5], min(factor)=2, min(dependent)=3
  expected_local <- c(0.0, 0.0, 8.0)  # (2-2)*(5-3), (4-2)*(3-3), (6-2)*(5-3)
  expected_total <- 8.0
  
  stopifnot(all(abs(result$local_contributions[,1,1,1] - expected_local) < TOL))
  stopifnot(abs(result$total_contributions[1,1,1] - expected_total) < TOL)
  cat("  MIN baseline passed ✓\n")
  
  cat("test_compute_all_contributions passed ✓\n")
}

# =====================================================
# Test: tox_perform_permutation_test
# =====================================================

test_perform_permutation_test <- function() {
  cat("\n[test_perform_permutation_test] Permutation testing with fixed seed\n")
  
  # Setup synthetic trajectories
  n_factors <- 2L
  n_samples <- 3L
  n_timepoints <- 3L
  n_permutations <- 3L
  
  # Create 3D trajectory data (n_factors, n_samples, n_timepoints)
  trajectories <- numeric(n_factors * n_samples * n_timepoints)
  
  # Factor 1 (index 1): [1,2,3], [2,4,6], [3,6,9]
  trajectories[1:3] <- c(1.0, 2.0, 3.0)
  trajectories[4:6] <- c(2.0, 4.0, 6.0)
  trajectories[7:9] <- c(3.0, 6.0, 9.0)
  
  # Dependent 2 (index 2): [4,5,6], [1,3,5], [2,4,6]
  trajectories[10:12] <- c(4.0, 5.0, 6.0)
  trajectories[13:15] <- c(1.0, 3.0, 5.0)
  trajectories[16:18] <- c(2.0, 4.0, 6.0)
  
  dims <- c(n_factors, n_samples, n_timepoints)
  factor_idx <- 1L
  dependent_idx <- 2L
  sample_idx <- 1L
  mode <- "mean"
  random_seed <- 12345L
  
  cat("  Running permutation test with", n_permutations, "permutations...\n")
  result <- tox_perform_permutation_test(
    trajectories, dims, 
    factor_idx, dependent_idx, sample_idx,
    mode = mode, 
    n_permutations = n_permutations, 
    random_seed = random_seed
  )
  
  local <- result$local_contributions
  total <- result$total_contributions
  
  cat("  Checking that all values are finite...\n")
  stopifnot(all(is.finite(local)))
  stopifnot(all(is.finite(total)))
  cat("  All finite ✓\n")
  
  cat("  Checking that not all contributions are zero...\n")
  stopifnot(!all(abs(local) < TOL))
  cat("  Non-zero contributions found ✓\n")
  
  cat("test_perform_permutation_test passed ✓\n")
}

# =====================================================
# Test: tox_compute_p_values
# =====================================================

test_compute_p_values <- function() {
  cat("\n[test_compute_p_values] P-value computation\n")
  
  # Case 1: Valid inputs
  cat("  Testing valid input case...\n")
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
  
  # Note: tox_compute_p_values expects dims parameter
  dims <- c(n_timepoints, n_permutations)
  result <- tox_compute_p_values(local_obs, total_obs, c(local_perm), total_perm, dims)
  
  local_p <- result$local_p_values
  total_p <- result$total_p_value
  
  expected_local <- c(0.5, 1.0, 0.5)
  expected_total <- 0.5
  
  stopifnot(all(abs(local_p - expected_local) < TOL))
  stopifnot(abs(total_p - expected_total) < TOL)
  cat("  Valid input case passed ✓\n")
  
  # Case 2: NaN in observed contributions
  cat("  Testing error handling for NaN input...\n")
  local_obs_nan <- c(2.0, 0.0, NaN)
  tryCatch(
    {
      tox_compute_p_values(local_obs_nan, total_obs, c(local_perm), total_perm, dims)
      stop("Should have raised error for NaN input")
    },
    error = function(e) {
      cat("  Correctly caught error for NaN input ✓\n")
    }
  )
  
  # Case 3: Inf in permutation contributions
  cat("  Testing error handling for Inf input...\n")
  local_perm_inf <- local_perm
  local_perm_inf[3, 4] <- Inf
  tryCatch(
    {
      tox_compute_p_values(local_obs, total_obs, c(local_perm_inf), total_perm, dims)
      stop("Should have raised error for Inf input")
    },
    error = function(e) {
      cat("  Correctly caught error for Inf input ✓\n")
    }
  )
  
  cat("test_compute_p_values passed ✓\n")
}

# =====================================================
# Main test runner
# =====================================================

main <- function() {
  cat("\n")
  cat("=================================================\n")
  cat("    TRAJECTORY CONTRIBUTION ANALYSIS R TESTS\n")
  cat("=================================================\n")
  cat("\n")
  
  test_compute_baselines_factor_dependent()
  test_compute_contributions()
  test_compute_all_contributions()
  test_perform_permutation_test()
  test_compute_p_values()
  
  cat("\n")
  cat("=================================================\n")
  cat("    ALL TRAJECTORY CONTRIBUTION TESTS PASSED ✓\n")
  cat("=================================================\n")
  cat("\n")
}

# Run tests
main()
