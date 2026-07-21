# Comprehensive R test suite for tissue versatility (mirrors Fortran unit tests)
# Source the main functions
source("rcpp/load_tensor_omics.R")
source("rcpp/test_helpers.R")

# The two selection counts are derived from the masks, but the Fortran takes them.
tissue_versatility <- function(expression_vectors, exp_vecs_selection_index, axes_selection) {
  compute_tissue_versatility(expression_vectors, exp_vecs_selection_index,
                             sum(exp_vecs_selection_index), axes_selection,
                             sum(axes_selection))
}

# 1. Uniform expression (should yield TV=0)
test_uniform_expression <- function() {
  expr <- matrix(2, nrow=3, ncol=1)
  res <- tissue_versatility(expr, c(TRUE), c(TRUE,TRUE,TRUE))
  assert_true(abs(res$tissue_versatilities[1]) < 1e-12)
  assert_true(abs(res$tissue_angles_deg[1]) < 1e-12)
}

# 2. Single axis expression (should yield TV=1)
test_single_axis_expression <- function() {
  expr <- matrix(c(0,0,5), nrow=3, ncol=1)
  res <- tissue_versatility(expr, c(TRUE), c(TRUE,TRUE,TRUE))
  assert_true(abs(res$tissue_versatilities[1] - 1) < 1e-12)
  assert_true(res$tissue_angles_deg[1] > 0)
}

# 3. Null vector (should yield TV=1, angle=90)
test_null_vector <- function() {
  expr <- matrix(0, nrow=3, ncol=1)
  res <- tissue_versatility(expr, c(TRUE), c(TRUE,TRUE,TRUE))
  assert_true(abs(res$tissue_versatilities[1] - 1) < 1e-12)
  assert_true(abs(res$tissue_angles_deg[1] - 90) < 1e-12)
}

# 4. Partial axis selection (subspace)
test_partial_axis_selection <- function() {
  expr <- matrix(c(1,2,3), nrow=3, ncol=1)
  res <- tissue_versatility(expr, c(TRUE), c(TRUE,FALSE,TRUE))
  assert_true(res$tissue_versatilities[1] >= 0 && res$tissue_versatilities[1] <= 1)
  assert_true(res$tissue_angles_deg[1] >= 0 && res$tissue_angles_deg[1] <= 90)
}

# 5. Mixed vectors (uniform, single axis, null)
test_mixed_vectors <- function() {
  expr <- matrix(c(1,1,1, 0,0,2, 0,0,0), nrow=3, ncol=3)
  res <- tissue_versatility(expr, c(TRUE,TRUE,TRUE), c(TRUE,TRUE,TRUE))
  assert_true(abs(res$tissue_versatilities[1]) < 1e-12)
  assert_true(abs(res$tissue_versatilities[2] - 1) < 1e-12)
  assert_true(abs(res$tissue_versatilities[3] - 1) < 1e-12)
  assert_true(abs(res$tissue_angles_deg[1]) < 1e-12)
  assert_true(res$tissue_angles_deg[2] > 0)
  assert_true(abs(res$tissue_angles_deg[3] - 90) < 1e-12)
}

# 6. Angle output in degrees for a known case (should be 45)
test_angle_degrees <- function() {
  expr <- matrix(c(1,0), nrow=2, ncol=1)
  res <- tissue_versatility(expr, c(TRUE), c(TRUE,TRUE))
  assert_true(abs(res$tissue_angles_deg[1] - 45) < 1e-12)
}

# 7. Multiple vectors selection
test_multiple_vectors_selection <- function() {
  expr <- matrix(c(1,1, 0,2, 0,0), nrow=2, ncol=3)
  res <- tissue_versatility(expr, c(TRUE,FALSE,TRUE), c(TRUE,TRUE))
  assert_true(abs(res$tissue_versatilities[1]) < 1e-12)
  assert_true(abs(res$tissue_versatilities[2] - 1) < 1e-12)
  assert_true(abs(res$tissue_angles_deg[1]) < 1e-5)
  assert_true(abs(res$tissue_angles_deg[2] - 90) < 1e-12)
}

# 8. High-dimensional vectors (4D, 5D)
test_high_dimensional_vectors <- function() {
  expr4 <- matrix(1, nrow=4, ncol=1)
  expr5 <- matrix(2, nrow=5, ncol=1)
  res4 <- tissue_versatility(expr4, c(TRUE), rep(TRUE,4))
  res5 <- tissue_versatility(expr5, c(TRUE), rep(TRUE,5))
  assert_true(abs(res4$tissue_versatilities[1]) < 1e-12)
  assert_true(abs(res4$tissue_angles_deg[1]) < 1e-12)
  assert_true(abs(res5$tissue_versatilities[1]) < 1e-12)
  assert_true(abs(res5$tissue_angles_deg[1]) < 1e-5)
}

# 9. Randomized vectors and axes
test_randomized_vectors_axes <- function() {
  set.seed(42)
  n_axes <- 5
  n_vecs <- 4
  expr <- matrix(runif(n_axes * n_vecs), nrow=n_axes, ncol=n_vecs)
  res <- tissue_versatility(expr, rep(TRUE,n_vecs), c(TRUE,FALSE,TRUE,FALSE,TRUE))
  assert_true(all(res$tissue_versatilities >= 0 & res$tissue_versatilities <= 1))
  assert_true(all(res$tissue_angles_deg >= 0 & res$tissue_angles_deg <= 90))
}

# 10. Numerical stability (comprehensive edge cases - mirrors Fortran test)
test_numerical_stability <- function() {
  # Mirror the test_comprehensive_edge_cases from Fortran
  # Case 1: Large numbers (should work normally - uniform → TV=0)  
  # Case 2: Small numbers above threshold (should work normally - uniform → TV=0)
  expr <- matrix(c(1e15,1e15,1e15, 1e-4,1e-4,1e-4), nrow=3, ncol=2)
  res <- tissue_versatility(expr, c(TRUE,TRUE), c(TRUE,TRUE,TRUE))
  assert_true(abs(res$tissue_versatilities[1]) < 1e-12)  # Large uniform → TV=0
  assert_true(abs(res$tissue_angles_deg[1]) < 1e-12)    # Large uniform → angle=0
  assert_true(abs(res$tissue_versatilities[2]) < 1e-12)  # Small uniform → TV=0  
  assert_true(abs(res$tissue_angles_deg[2]) < 1e-12)    # Small uniform → angle=0
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
  res <- tissue_versatility(expr, c(TRUE,TRUE), c(TRUE,TRUE,TRUE))
  assert_true(abs(res$tissue_versatilities[1] - 1) < 1e-12)  # Below threshold → TV=1
  assert_true(abs(res$tissue_angles_deg[1] - 90) < 1e-12)   # Below threshold → angle=90°
  assert_true(abs(res$tissue_versatilities[2] - 1) < 1e-12)  # Underflow → TV=1
  assert_true(abs(res$tissue_angles_deg[2] - 90) < 1e-12)   # Underflow → angle=90°
}

# 11. Invalid input: Empty input arrays provided. (should throw error with code 202)
test_invalid_input_no_axes <- function() {
  expr <- matrix(c(1,2,3), nrow=3, ncol=1)
  error_caught <- FALSE
  tryCatch({
    tissue_versatility(expr, c(TRUE), c(FALSE,FALSE,FALSE))
  }, error = function(e) {
    error_caught <<- TRUE
    # Check that the error message contains the expected text
    assert_true(grepl("empty input", e$message, ignore.case = TRUE))
  })
  assert_true(error_caught)  # Make sure an error was actually thrown
}

# 12. Multiple selection, partial axes
test_multiple_selection_partial_axes <- function() {
  expr <- matrix(c(1,2, 3,4, 5,6), nrow=2, ncol=3)
  res <- tissue_versatility(expr, c(TRUE,FALSE,TRUE), c(TRUE,FALSE))
  assert_true(length(res$tissue_versatilities) == 2)
  assert_true(all(res$tissue_versatilities >= 0 & res$tissue_versatilities <= 1))
}

run_all_tests()
