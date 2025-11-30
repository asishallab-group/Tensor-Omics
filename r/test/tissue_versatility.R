# Comprehensive R test suite for tissue versatility (mirrors Fortran unit tests)
# Source the main functions
source("r/tensoromics_functions.R")

# 1. Uniform expression (should yield TV=0)
test_uniform_expression <- function() {
  expr <- matrix(2, nrow=3, ncol=1)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE), c(TRUE,TRUE,TRUE))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)
  cat("test_uniform_expression passed\n")
}

# 2. Single axis expression (should yield TV=1)
test_single_axis_expression <- function() {
  expr <- matrix(c(0,0,5), nrow=3, ncol=1)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE), c(TRUE,TRUE,TRUE))
  stopifnot(abs(res$tissue_versatilities[1] - 1) < 1e-12)
  stopifnot(res$tissue_angles_deg[1] > 0)
  cat("test_single_axis_expression passed\n")
}

# 3. Null vector (should yield TV=1, angle=90)
test_null_vector <- function() {
  expr <- matrix(0, nrow=3, ncol=1)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE), c(TRUE,TRUE,TRUE))
  stopifnot(abs(res$tissue_versatilities[1] - 1) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1] - 90) < 1e-12)
  cat("test_null_vector passed\n")
}

# 4. Partial axis selection (subspace)
test_partial_axis_selection <- function() {
  expr <- matrix(c(1,2,3), nrow=3, ncol=1)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE), c(TRUE,FALSE,TRUE))
  stopifnot(res$tissue_versatilities[1] >= 0 && res$tissue_versatilities[1] <= 1)
  stopifnot(res$tissue_angles_deg[1] >= 0 && res$tissue_angles_deg[1] <= 90)
  cat("test_partial_axis_selection passed\n")
}

# 5. Mixed vectors (uniform, single axis, null)
test_mixed_vectors <- function() {
  expr <- matrix(c(1,1,1, 0,0,2, 0,0,0), nrow=3, ncol=3)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE,TRUE,TRUE), c(TRUE,TRUE,TRUE))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res$tissue_versatilities[2] - 1) < 1e-12)
  stopifnot(abs(res$tissue_versatilities[3] - 1) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)
  stopifnot(res$tissue_angles_deg[2] > 0)
  stopifnot(abs(res$tissue_angles_deg[3] - 90) < 1e-12)
  cat("test_mixed_vectors passed\n")
}

# 6. Angle output in degrees for a known case (should be 45)
test_angle_degrees <- function() {
  expr <- matrix(c(1,0), nrow=2, ncol=1)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE), c(TRUE,TRUE))
  stopifnot(abs(res$tissue_angles_deg[1] - 45) < 1e-12)
  cat("test_angle_degrees passed\n")
}

# 7. Multiple vectors selection
test_multiple_vectors_selection <- function() {
  expr <- matrix(c(1,1, 0,2, 0,0), nrow=2, ncol=3)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE,FALSE,TRUE), c(TRUE,TRUE))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res$tissue_versatilities[2] - 1) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[2] - 90) < 1e-12)
  cat("test_multiple_vectors_selection passed\n")
}

# 8. High-dimensional vectors (4D, 5D)
test_high_dimensional_vectors <- function() {
  expr4 <- matrix(1, nrow=4, ncol=1)
  expr5 <- matrix(2, nrow=5, ncol=1)
  res4 <- tox_calculate_tissue_versatility(expr4, c(TRUE), rep(TRUE,4))
  res5 <- tox_calculate_tissue_versatility(expr5, c(TRUE), rep(TRUE,5))
  stopifnot(abs(res4$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res4$tissue_angles_deg[1]) < 1e-12)
  stopifnot(abs(res5$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res5$tissue_angles_deg[1]) < 1e-12)
  cat("test_high_dimensional_vectors passed\n")
}

# 9. Randomized vectors and axes
test_randomized_vectors_axes <- function() {
  set.seed(42)
  n_axes <- 5
  n_vecs <- 4
  expr <- matrix(runif(n_axes * n_vecs), nrow=n_axes, ncol=n_vecs)
  res <- tox_calculate_tissue_versatility(expr, rep(TRUE,n_vecs), c(TRUE,FALSE,TRUE,FALSE,TRUE))
  stopifnot(all(res$tissue_versatilities >= 0 & res$tissue_versatilities <= 1))
  stopifnot(all(res$tissue_angles_deg >= 0 & res$tissue_angles_deg <= 90))
  cat("test_randomized_vectors_axes passed\n")
}

# 10. Numerical stability (comprehensive edge cases - mirrors Fortran test)
test_numerical_stability <- function() {
  # Mirror the test_comprehensive_edge_cases from Fortran
  # Case 1: Large numbers (should work normally - uniform → TV=0)  
  # Case 2: Small numbers above threshold (should work normally - uniform → TV=0)
  expr <- matrix(c(1e15,1e15,1e15, 1e-4,1e-4,1e-4), nrow=3, ncol=2)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE,TRUE), c(TRUE,TRUE,TRUE))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)  # Large uniform → TV=0
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)    # Large uniform → angle=0
  stopifnot(abs(res$tissue_versatilities[2]) < 1e-12)  # Small uniform → TV=0  
  stopifnot(abs(res$tissue_angles_deg[2]) < 1e-12)    # Small uniform → angle=0
  cat("test_numerical_stability passed\n")
}

# 10b. Epsilon threshold protection (mirrors Fortran epsilon stability test)
test_epsilon_threshold_protection <- function() {
  # Mirror test_epsilon_threshold_stability from Fortran exactly
  # eps_sqrt ≈ 1.49e-8, large_component = 1e-5
  eps_sqrt <- sqrt(.Machine$double.eps)  # R equivalent of sqrt(epsilon(1.0_real64))
  large_component <- 1e-5
  
  # Test case 2: Vector with norm slightly below sqrt(epsilon) threshold → should get TV=1, angle=90°
  # Test case 4: Underflow case → should get TV=1, angle=90°
  expr <- matrix(c(eps_sqrt*0.5/sqrt(3), eps_sqrt*0.5/sqrt(3), eps_sqrt*0.5/sqrt(3),   # Case 2
                   1e-200, 1e-200, 1e-200), nrow=3, ncol=2)                              # Case 4
  res <- tox_calculate_tissue_versatility(expr, c(TRUE,TRUE), c(TRUE,TRUE,TRUE))
  stopifnot(abs(res$tissue_versatilities[1] - 1) < 1e-12)  # Below threshold → TV=1
  stopifnot(abs(res$tissue_angles_deg[1] - 90) < 1e-12)   # Below threshold → angle=90°
  stopifnot(abs(res$tissue_versatilities[2] - 1) < 1e-12)  # Underflow → TV=1
  stopifnot(abs(res$tissue_angles_deg[2] - 90) < 1e-12)   # Underflow → angle=90°
  cat("test_epsilon_threshold_protection passed\n")
}

# 11. Invalid input: Empty input arrays provided. (should throw error with code 202)
test_invalid_input_no_axes <- function() {
  expr <- matrix(c(1,2,3), nrow=3, ncol=1)
  error_caught <- FALSE
  tryCatch({
    tox_calculate_tissue_versatility(expr, c(TRUE), c(FALSE,FALSE,FALSE))
  }, error = function(e) {
    error_caught <<- TRUE
    # Check that the error message contains the expected text
    stopifnot(grepl("Empty input arrays provided.", e$message))
  })
  stopifnot(error_caught)  # Make sure an error was actually thrown
  cat("test_invalid_input_no_axes passed\n")
}

# 12. Multiple selection, partial axes
test_multiple_selection_partial_axes <- function() {
  expr <- matrix(c(1,2, 3,4, 5,6), nrow=2, ncol=3)
  res <- tox_calculate_tissue_versatility(expr, c(TRUE,FALSE,TRUE), c(TRUE,FALSE))
  stopifnot(length(res$tissue_versatilities) == 2)
  stopifnot(all(res$tissue_versatilities >= 0 & res$tissue_versatilities <= 1))
  cat("test_multiple_selection_partial_axes passed\n")
}

# Run all tests
cat("=================================================\n")
cat("    TISSUE VERSATILITY FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

test_uniform_expression()
test_single_axis_expression()
test_null_vector()
test_partial_axis_selection()
test_mixed_vectors()
test_angle_degrees()
test_multiple_vectors_selection()
test_high_dimensional_vectors()
test_randomized_vectors_axes()
test_numerical_stability()
test_epsilon_threshold_protection()
test_invalid_input_no_axes()
test_multiple_selection_partial_axes()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all tissue versatility R interface tests passed! ✓\n")
