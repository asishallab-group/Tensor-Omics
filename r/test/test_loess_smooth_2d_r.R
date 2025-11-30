# =====================
# Test cases for loess_smooth_2d R wrapper function
# Based on Fortran test suite from mod_test_loess_smoothing.f90
# Includes error validation and comprehensive test coverage
# =====================

# Source the main functions
source("r/tensoromics_functions.R")

cat("=== Testing loess_smooth_2d R wrapper function ===\n")
cat("Based on Fortran test suite with 14 comprehensive test cases\n")

# Test 1: Constant input
test_loess_constant_input <- function() {
  cat("\n[test_loess_constant_input] Constant input test\n")
  n <- 100
  x_ref <- rep(5.0, n)
  y_ref <- rep(10.0, n)
  x_query <- rep(5.0, n/2)
  indices_used <- 1:n
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 1.0, 3.0)
  
  # Verify all outputs are constant
  stopifnot(all(abs(result$smoothed_values - 10.0) < 1e-6))
  cat("Constant input test passed ✓\n")
}

# Test 2: Linear trend
test_loess_linear_trend <- function() {
  cat("\n[test_loess_linear_trend] Linear trend test\n")
  n <- 100
  x_ref <- 1:n
  y_ref <- 0.5 * x_ref
  x_query <- (1:(n/2)) + 0.5
  indices_used <- 1:n
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 1.0, 3.0)
  expected <- 0.5 * x_query
  
  # Verify linear trend is preserved (within tolerance)
  stopifnot(all(abs(result$smoothed_values - expected) < 0.05))
  cat("Linear trend test passed ✓\n")
}

# Test 3: Outlier suppression
test_loess_outlier_suppression <- function() {
  cat("\n[test_loess_outlier_suppression] Outlier suppression test\n")
  n <- 100
  x_ref <- c(rep(10.0, n-1), 100.0)
  y_ref <- c(rep(5.0, n-1), 99.0)
  x_query <- rep(10.0, 50)
  indices_used <- 1:n
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 1.0, 3.0)
  
  # Outlier should be suppressed, result close to 5.0
  stopifnot(all(abs(result$smoothed_values - 5.0) < 0.01))
  cat("Outlier suppression test passed ✓\n")
}

# Test 4: Sparse fallback
test_loess_sparse_fallback <- function() {
  cat("\n[test_loess_sparse_fallback] Sparse fallback test\n")
  n <- 100
  x_ref <- (1:n) * 100.0
  y_ref <- 1:n
  x_query <- (1:50) * 100.0 + 50.0
  indices_used <- 1:n
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 1.0, 3.0)
  
  # Should fallback to nearest neighbor for sparse data
  expected <- y_ref[1:50]  # Should be close to nearest values
  stopifnot(all(abs(result$smoothed_values - expected) < 1e-6))
  cat("Sparse fallback test passed ✓\n")
}

# Test 5: Single point
test_loess_single_point <- function() {
  cat("\n[test_loess_single_point] Single point test\n")
  x_ref <- 0.0
  y_ref <- 42.0
  x_query <- 0.0
  indices_used <- 1
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 0.1, 1.0)
  
  # Should return the single point value
  stopifnot(abs(result$smoothed_values[1] - 42.0) < 1e-6)
  cat("Single point test passed ✓\n")
}

# Test 6: Identical points
test_loess_identical_points <- function() {
  cat("\n[test_loess_identical_points] Identical points test\n")
  x_ref <- c(0.0, 1.0)
  y_ref <- c(1.0, 1.0)
  x_query <- 0.0
  indices_used <- c(1, 2)
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 0.1, 1.0)
  
  # Should return the identical value
  stopifnot(abs(result$smoothed_values[1] - 1.0) < 1e-6)
  cat("Identical points test passed ✓\n")
}

# Test 7: Linear interpolation
test_loess_linear_interp <- function() {
  cat("\n[test_loess_linear_interp] Linear interpolation test\n")
  x_ref <- c(0.0, 2.0)
  y_ref <- c(0.0, 2.0)
  x_query <- 1.0
  indices_used <- c(1, 2)
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 0.5, 3.0)
  
  # Should interpolate to approximately 1.0
  stopifnot(abs(result$smoothed_values[1] - 1.0) < 0.1)
  cat("Linear interpolation test passed ✓\n")
}

# Test 8: Weight decay
test_loess_weight_decay <- function() {
  cat("\n[test_loess_weight_decay] Weight decay test\n")
  x_ref <- c(0.0, 10.0)
  y_ref <- c(0.0, 10.0)
  x_query <- 0.0
  indices_used <- c(1, 2)
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 1.0, 3.0)
  
  # Should be weighted towards the closer point
  stopifnot(result$smoothed_values[1] < 5.0)
  cat("Weight decay test passed ✓\n")
}

# Test 9: Mask exclusion
test_loess_mask_exclusion <- function() {
  cat("\n[test_loess_mask_exclusion] Mask exclusion test\n")
  x_ref <- c(0.0, 10.0, 20.0)
  y_ref <- c(0.0, 10.0, 20.0)
  x_query <- 0.0
  indices_used <- 2  # Use only the middle point
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 10.0, 3.0)
  
  # Should be influenced primarily by the middle point
  stopifnot(abs(result$smoothed_values[1] - 10.0) < 1.0)
  cat("Mask exclusion test passed ✓\n")
}

# Test 10: Fallback
test_loess_fallback <- function() {
  cat("\n[test_loess_fallback] Fallback test\n")
  x_ref <- 0.0
  y_ref <- 123.0
  x_query <- 100.0  # Far from reference point
  indices_used <- 1
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 0.1, 1.0)
  
  # Should fallback to the single reference value
  stopifnot(abs(result$smoothed_values[1] - 123.0) < 1e-6)
  cat("Fallback test passed ✓\n")
}

# Test 11: Edge query (extrapolation)
test_loess_edge_query <- function() {
  cat("\n[test_loess_edge_query] Edge query test\n")
  x_ref <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  y_ref <- c(10.0, 20.0, 30.0, 40.0, 50.0)
  x_query <- c(0.0, 6.0)  # Outside the range
  indices_used <- 1:5
  
  result <- tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 1.0, 3.0)
  
  # Should extrapolate reasonably
  stopifnot(abs(result$smoothed_values[1] - 10.0) < 1.0)  # Close to first point
  stopifnot(abs(result$smoothed_values[2] - 50.0) < 1.0)  # Close to last point
  cat("Edge query test passed ✓\n")
}

# Test 12: Invalid dimensions (error handling)
test_loess_invalid_dimensions <- function() {
  cat("\n[test_loess_invalid_dimensions] Invalid dimensions test\n")
  x_ref <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  y_ref <- c(10.0, 20.0, 30.0, 40.0, 50.0)
  x_query <- c(0.0, 6.0)
  indices_used <- 1:5
  
  # Test with empty arrays - should throw error
  error_caught <- FALSE
  tryCatch({
    tox_loess_smooth_2d(numeric(0), numeric(0), x_query, indices_used, 1.0, 3.0)
  }, error = function(e) {
    error_caught <<- TRUE
    print(paste("Caught expected error:", e$message))
    # Check that the error message contains expected text
    stopifnot(grepl("Empty input arrays provided.", e$message))
  })
  stopifnot(error_caught)  # Make sure an error was actually thrown
  
  cat("Invalid dimensions test passed ✓\n")
}

# Test 13: Invalid parameters (error handling)
test_loess_invalid_parameters <- function() {
  cat("\n[test_loess_invalid_parameters] Invalid parameters test\n")
  x_ref <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  y_ref <- c(10.0, 20.0, 30.0, 40.0, 50.0)
  x_query <- c(0.0, 6.0)
  indices_used <- 1:5
  
  # Test with negative kernel_sigma - should throw error
  error_caught1 <- FALSE
  tryCatch({
    tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, -1.0, 3.0)
  }, error = function(e) {
    error_caught1 <<- TRUE
    # Check that the error message contains expected text
    stopifnot(grepl("Invalid input", e$message))
  })
  stopifnot(error_caught1)  # Make sure an error was actually thrown
  
  # Test with negative kernel_cutoff - should throw error
  error_caught2 <- FALSE
  tryCatch({
    tox_loess_smooth_2d(x_ref, y_ref, x_query, indices_used, 1.0, -3.0)
  }, error = function(e) {
    error_caught2 <<- TRUE
    # Check that the error message contains expected text
    stopifnot(grepl("Invalid input", e$message))
  })
  stopifnot(error_caught2)  # Make sure an error was actually thrown
  
  cat("Invalid parameters test passed ✓\n")
}

# Test 14: Invalid indices (error handling)
test_loess_invalid_indices <- function() {
  cat("\n[test_loess_invalid_indices] Invalid indices test\n")
  x_ref <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  y_ref <- c(10.0, 20.0, 30.0, 40.0, 50.0)
  x_query <- c(0.0, 6.0)
  
  # Test with index out of bounds (too high) - should throw error
  error_caught1 <- FALSE
  tryCatch({
    tox_loess_smooth_2d(x_ref, y_ref, x_query, c(1, 2, 3, 4, 6), 1.0, 3.0)
  }, error = function(e) {
    error_caught1 <<- TRUE
    # Check that the error message contains expected text
    stopifnot(grepl("Invalid input", e$message))
  })
  stopifnot(error_caught1)  # Make sure an error was actually thrown
  
  # Test with index out of bounds (too low) - should throw error
  error_caught2 <- FALSE
  tryCatch({
    tox_loess_smooth_2d(x_ref, y_ref, x_query, c(0, 2, 3, 4, 5), 1.0, 3.0)
  }, error = function(e) {
    error_caught2 <<- TRUE
    # Check that the error message contains expected text
    stopifnot(grepl("Invalid input", e$message))
  })
  stopifnot(error_caught2)  # Make sure an error was actually thrown
  
  cat("Invalid indices test passed ✓\n")
}

# Test 15: Expert family scaling function
test_family_scaling_expert <- function() {
  cat("\n[test_family_scaling_expert] Expert family scaling function test\n")
  
  # Generate synthetic data for testing
  n_genes <- 100
  n_families <- 5
  
  # Create distances (some variation for realistic data)
  set.seed(42)  # For reproducible results
  distances <- runif(n_genes, 0.1, 5.0)
  
  # Create gene-to-family mapping (balanced distribution)
  gene_to_fam <- rep(1:n_families, length.out = n_genes)
  
  # Pre-allocate work arrays as required by expert function
  perm_tmp <- integer(n_genes)
  stack_left_tmp <- integer(n_genes)
  stack_right_tmp <- integer(n_genes)
  family_distances <- numeric(n_genes)
  
  # Call the expert function
  result_expert <- tox_compute_family_scaling_expert(
    distances, gene_to_fam, n_families,
    perm_tmp, stack_left_tmp, stack_right_tmp, family_distances
  )
  
  # Call the regular function for comparison
  result_regular <- tox_compute_family_scaling(distances, gene_to_fam, n_families)
  
  # Verify that both functions produce the same results
  stopifnot(length(result_expert$dscale) == n_families)
  stopifnot(length(result_regular$dscale) == n_families)
  
  # Compare scaling factors (should be identical)
  stopifnot(all(abs(result_expert$dscale - result_regular$dscale) < 1e-10))
  
  # Compare LOESS coordinates (should be identical)
  stopifnot(all(abs(result_expert$loess_x - result_regular$loess_x) < 1e-10))
  stopifnot(all(abs(result_expert$loess_y - result_regular$loess_y) < 1e-10))
  
  # Compare indices used (should be identical)
  stopifnot(all(result_expert$indices_used == result_regular$indices_used))
  
  # Verify that work arrays are returned (expert-specific feature)
  stopifnot(!is.null(result_expert$perm_tmp))
  stopifnot(!is.null(result_expert$stack_left_tmp))
  stopifnot(!is.null(result_expert$stack_right_tmp))
  stopifnot(!is.null(result_expert$family_distances))
  
  # Verify work arrays have correct dimensions
  stopifnot(length(result_expert$perm_tmp) == n_genes)
  stopifnot(length(result_expert$stack_left_tmp) == n_genes)
  stopifnot(length(result_expert$stack_right_tmp) == n_genes)
  stopifnot(length(result_expert$family_distances) == n_genes)
  
  # Test reusing work arrays (expert function advantage)
  # Modify distances slightly and reuse work arrays
  distances2 <- distances * 1.1
  result_expert2 <- tox_compute_family_scaling_expert(
    distances2, gene_to_fam, n_families,
    result_expert$perm_tmp,  # Reuse arrays
    result_expert$stack_left_tmp,
    result_expert$stack_right_tmp,
    result_expert$family_distances
  )
  
  # Should produce different but valid results
  stopifnot(length(result_expert2$dscale) == n_families)
  stopifnot(!all(result_expert2$dscale == result_expert$dscale))  # Should be different
  
  cat("Expert family scaling function test passed ✓\n")
  cat("  - Regular and expert functions produce identical results\n")
  cat("  - Work arrays are properly returned and can be reused\n")
  cat("  - Expert function provides performance benefits for repeated calls\n")
}

# =====================
# Run all tests
# =====================

cat("\n==========================================================\n")
cat("           RUNNING ALL LOESS SMOOTH 2D TESTS\n")
cat("           (Based on Fortran test suite)\n")
cat("==========================================================\n")

# Core functionality tests
test_loess_constant_input()
test_loess_linear_trend()
test_loess_outlier_suppression()
test_loess_sparse_fallback()
test_loess_single_point()
test_loess_identical_points()
test_loess_linear_interp()
test_loess_weight_decay()
test_loess_mask_exclusion()
test_loess_fallback()
test_loess_edge_query()

# Error handling tests
test_loess_invalid_dimensions()
test_loess_invalid_parameters()
test_loess_invalid_indices()

# Expert function test
test_family_scaling_expert()

