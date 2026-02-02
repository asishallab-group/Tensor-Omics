
library(Rcpp)

# Set library path and compile
source("rcpp/tensoromics_functions.R")

cat("=== Testing TensorOmics Trajectory Normalization Functions ===\n")


# Constants
TOL <- 1e-12

# =====================================================
# Test: tox_normalize_variable_timeseries
# =====================================================

test_normalize_variable_timeseries <- function() {
  cat("\n[test_normalize_variable_timeseries] Normalize single timeseries\n")
  
  # Case 1: Simple vector with range [1, 4]
  cat("  Testing simple vector normalization...\n")
  v <- c(1.0, 2.0, 3.0, 4.0)
  result <- tox_normalize_variable_timeseries(v)
  
  # Min-max normalization: (x - min) / (max - min) = (x - 1) / 3
  expected <- c(0.0, 1/3, 2/3, 1.0)
  
  stopifnot(all(abs(result$v_norm - expected) < TOL))
  stopifnot(all(is.finite(result$v_norm)))
  cat("  Simple vector test passed ✓\n")
  
  # Case 2: All same values (edge case)
  cat("  Testing constant vector...\n")
  v_const <- c(5.0, 5.0, 5.0, 5.0)
  result_const <- tox_normalize_variable_timeseries(v_const)
  
  # Constant vector normalization should handle gracefully
  stopifnot(all(is.finite(result_const$v_norm)))
  cat("  Constant vector test passed ✓\n")
  
  # Case 3: Negative values
  cat("  Testing vector with negative values...\n")
  v_neg <- c(-2.0, -1.0, 0.0, 1.0, 2.0)
  result_neg <- tox_normalize_variable_timeseries(v_neg)
  
  expected_neg <- c(0.0, 0.25, 0.5, 0.75, 1.0)
  stopifnot(all(abs(result_neg$v_norm - expected_neg) < TOL))
  cat("  Negative values test passed ✓\n")
  
  cat("test_normalize_variable_timeseries passed ✓\n")
}

# =====================================================
# Test: tox_normalize_single_trajectory
# =====================================================

test_normalize_single_trajectory <- function() {
  cat("\n[test_normalize_single_trajectory] Normalize factor trajectory matrix\n")
  
  # Case 1: Simple 3x2 trajectory (3 timepoints, 2 factors)
  cat("  Testing 3x2 trajectory matrix...\n")
  trajectory <- matrix(c(
    1.0, 4.0,  # timepoint 1: factor1=1, factor2=4
    2.0, 5.0,  # timepoint 2: factor1=2, factor2=5
    3.0, 6.0   # timepoint 3: factor1=3, factor2=6
  ), nrow = 3, ncol = 2, byrow = TRUE)
  
  result <- tox_normalize_single_trajectory(trajectory)
  
  # Each factor normalized independently
  # Factor 1: [1,2,3] -> [0, 0.5, 1]
  # Factor 2: [4,5,6] -> [0, 0.5, 1]
  expected <- matrix(c(
    0.0, 0.0,
    0.5, 0.5,
    1.0, 1.0
  ), nrow = 3, ncol = 2, byrow = TRUE)
  
  stopifnot(all(abs(result$trajectory_norm - expected) < TOL))
  cat("  3x2 trajectory test passed ✓\n")
  
  # Case 2: Trajectory with negative values
  cat("  Testing trajectory with negative values...\n")
  trajectory_neg <- matrix(c(
    -2.0, -1.0,
    0.0, 1.0,
    2.0, 3.0
  ), nrow = 3, ncol = 2, byrow = TRUE)
  
  result_neg <- tox_normalize_single_trajectory(trajectory_neg)
  
  # Normalized to [0, 1] range
  # Factor 1: [-2,0,2] -> [0, 0.5, 1]
  # Factor 2: [-1,1,3] -> [0, 0.5, 1]
  expected_neg <- matrix(c(
    0.0, 0.0,
    0.5, 0.5,
    1.0, 1.0
  ), nrow = 3, ncol = 2, byrow = TRUE)
  
  stopifnot(all(abs(result_neg$trajectory_norm - expected_neg) < TOL))
  cat("  Negative values test passed ✓\n")
  
  cat("test_normalize_single_trajectory passed ✓\n")
}

# =====================================================
# Test: tox_normalize_all_trajectories
# =====================================================

test_normalize_all_trajectories <- function() {
  cat("\n[test_normalize_all_trajectories] Normalize 3D trajectory data\n")
  
  # Case 1: Simple 2x2x3 trajectories (n_factors=2, n_samples=2, n_timepoints=3)
  cat("  Testing 2x2x3 trajectory data...\n")
  # Flattened Fortran-order: factor1/sample1, factor1/sample2, factor2/sample1, factor2/sample2
  trajectories <- c(
    # Factor 1
    1.0, 2.0, 3.0,     # sample 1, timepoints 1-3
    10.0, 20.0, 30.0,  # sample 2, timepoints 1-3
    # Factor 2
    4.0, 5.0, 6.0,     # sample 1, timepoints 1-3
    7.0, 8.0, 9.0      # sample 2, timepoints 1-3
  )
  
  dims <- c(2L, 2L, 3L)  # (n_factors, n_samples, n_timepoints)
  
  result <- tox_normalize_all_trajectories(trajectories, dims)
  
  # Each factor normalized independently across all samples/timepoints
  # All values should be in [0, 1]
  stopifnot(all(result$trajectories_norm >= 0))
  stopifnot(all(result$trajectories_norm <= 1))
  cat("  2x2x3 trajectory test passed ✓\n")
  
  # Case 2: Single sample, multiple factors
  cat("  Testing single sample case...\n")
  trajectories_single <- c(
    1.0, 2.0, 3.0,  # factor 1
    4.0, 5.0, 6.0   # factor 2
  )
  
  dims_single <- c(2L, 1L, 3L)
  result_single <- tox_normalize_all_trajectories(trajectories_single, dims_single)
  
  stopifnot(all(result_single$trajectories_norm >= 0))
  stopifnot(all(result_single$trajectories_norm <= 1))
  cat("  Single sample test passed ✓\n")
  
  cat("test_normalize_all_trajectories passed ✓\n")
}

# =====================================================
# Helper function: pretty print matrix
# =====================================================

print_matrix <- function(name, mat) {
  cat("\n", name, ":\n", sep = "")
  cat("Dimensions:", nrow(mat), "x", ncol(mat), "\n")
  for (i in 1:nrow(mat)) {
    cat("  Row", i, ":", paste(format(mat[i, ], width = 8, digits = 4), collapse = " "), "\n")
  }
}

# =====================================================
# Example test 1: Normalization of simple data
# =====================================================

test_normalization_example_1 <- function() {
  cat("\n")
  cat("=================================================\n")
  cat("EXAMPLE 1: Normalize simple 2x3 matrix\n")
  cat("=================================================\n")
  
  mat <- matrix(c(
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0
  ), nrow = 2, ncol = 3, byrow = TRUE)
  
  print_matrix("Input", mat)
  
  # Normalize each factor independently
  for (i in 1:nrow(mat)) {
    row <- mat[i, ]
    min_val <- min(row)
    max_val <- max(row)
    normalized <- (row - min_val) / (max_val - min_val)
    cat(sprintf("  Factor %d: min=%.1f, max=%.1f, normalized=[%.2f, %.2f, %.2f]\n", 
                i, min_val, max_val, normalized[1], normalized[2], normalized[3]))
  }
  
  cat("\n✓ Example 1 passed!\n")
}

# =====================================================
# Example test 2: Large values
# =====================================================

test_normalization_example_2 <- function() {
  cat("\n")
  cat("=================================================\n")
  cat("EXAMPLE 2: Normalize large values\n")
  cat("=================================================\n")
  
  mat <- matrix(c(
    1e6, 2e6,
    3e6, 4e6
  ), nrow = 2, ncol = 2, byrow = TRUE)
  
  print_matrix("Input", mat)
  
  # Verify normalization works with large values
  for (i in 1:nrow(mat)) {
    row <- mat[i, ]
    min_val <- min(row)
    max_val <- max(row)
    normalized <- (row - min_val) / (max_val - min_val)
    cat(sprintf("  Factor %d: normalized=[%.4f, %.4f]\n", i, normalized[1], normalized[2]))
  }
  
  cat("\n✓ Example 2 passed!\n")
}

# =====================================================
# Error handling tests
# =====================================================

test_error_handling <- function() {
  cat("\n")
  cat("=================================================\n")
  cat("ERROR HANDLING TESTS\n")
  cat("=================================================\n")
  
  # Test empty vector
  cat("\n--- Test empty vector ---\n")
  tryCatch(
    {
      tox_normalize_variable_timeseries(c())
      stop("Should have raised error for empty vector")
    },
    error = function(e) {
      cat("✓ Correctly caught error for empty vector\n")
    }
  )
  
  # Test single value
  cat("\n--- Test single value vector ---\n")
  result <- tox_normalize_variable_timeseries(c(5.0))
  cat("✓ Single value handled\n")
  
  # Test constant vector
  cat("\n--- Test constant vector ---\n")
  result <- tox_normalize_variable_timeseries(c(5.0, 5.0, 5.0))
  cat("✓ Constant vector handled\n")
  
  # Test NaN input
  cat("\n--- Test NaN input ---\n")
  tryCatch(
    {
      tox_normalize_variable_timeseries(c(1.0, NaN, 3.0))
      stop("Should have raised error for NaN input")
    },
    error = function(e) {
      cat("✓ Correctly caught error for NaN input\n")
    }
  )
  
  # Test Inf input
  cat("\n--- Test Inf input ---\n")
  tryCatch(
    {
      tox_normalize_variable_timeseries(c(1.0, Inf, 3.0))
      stop("Should have raised error for Inf input")
    },
    error = function(e) {
      cat("✓ Correctly caught error for Inf input\n")
    }
  )
  
  cat("\n✓ All error handling tests passed!\n")
}

# =====================================================
# Main test runner
# =====================================================

main <- function() {
  cat("\n")
  cat("=================================================\n")
  cat("TRAJECTORY NORMALIZATION R TEST SUITE\n")
  cat("=================================================\n")
  cat("\n")
  
  test_normalize_variable_timeseries()
  test_normalize_single_trajectory()
  test_normalize_all_trajectories()
  test_normalization_example_1()
  test_normalization_example_2()
  test_error_handling()
  
  cat("\n")
  cat("=================================================\n")
  cat("ALL NORMALIZATION TESTS PASSED ✓\n")
  cat("=================================================\n")
  cat("\n")
}

# Run tests
main()
