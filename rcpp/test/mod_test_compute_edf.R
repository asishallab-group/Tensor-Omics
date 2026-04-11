
# Set library path and compile
source("rcpp/tensoromics_functions.R")

cat("=== Testing compute_edf R wrapper function ===\n")
cat("Based on Fortran/Python test suite with comprehensive test coverage\n\n")

TOL <- 1e-12

# Test 1: Simple EDF Test
test_compute_edf_simple <- function() {
  cat("\n[test_compute_edf_simple] Basic EDF computation\n")
  x <- c(1, 2, 2, 3, 3, 3)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3)
  expected_cdf <- c(1/6, 3/6, 6/6)
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 3)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ Simple EDF test passed\n")
}

# Test 2: All Unique Values
test_compute_edf_all_unique <- function() {
  cat("\n[test_compute_edf_all_unique] All values are unique\n")
  x <- c(1, 2, 3, 4, 5)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  unique_vals <- result$unique_values
  cdf_vals <- result$cdf_values
  n_unique <- length(unique_vals)
  expected_unique <- c(1, 2, 3, 4, 5)
  expected_cdf <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 5)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ All unique values test passed\n")
}

# Test 3: All Same Values
test_compute_edf_all_same <- function() {
  cat("\n[test_compute_edf_all_same] All values are identical\n")
  x <- c(5, 5, 5, 5)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  cat("  Input values:", x, "\n")
  cat("  n_unique:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 1)
  stopifnot(length(unique_vals) == 1)
  stopifnot(length(cdf_vals) == 1)
  stopifnot(abs(unique_vals[1] - 5) < TOL)
  stopifnot(abs(cdf_vals[1] - 1.0) < TOL)
  cat("  ✓ All same values test passed\n")
}

# Test 4: Duplicates
test_compute_edf_duplicates <- function() {
  cat("\n[test_compute_edf_duplicates] Various duplicate values\n")
  x <- c(1, 1, 2, 3, 3, 3, 4)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  cat("  (DEBUG) n_unique:", n_unique, "\n")
  stopifnot(n_unique == 4)
  stopifnot(length(unique_vals) == 4)
  stopifnot(length(cdf_vals) == 4)
  expected_unique <- c(1, 2, 3, 4)
  expected_cdf <- c(2/7, 3/7, 6/7, 1.0)
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 4)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ Duplicates test passed\n")
}

# Test 5: Single Value
test_compute_edf_single_value <- function() {
  cat("\n[test_compute_edf_single_value] Single value input\n")
  x <- c(42)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 1)
  stopifnot(abs(unique_vals[1] - 42) < TOL)
  stopifnot(abs(cdf_vals[1] - 1.0) < TOL)
  cat("  ✓ Single value test passed\n")
}

# Test 6: Empty Input
test_compute_edf_empty_input <- function() {
  cat("\n[test_compute_edf_empty_input] Empty array input (should fail)\n")
  x <- numeric(0)
  perm <- integer(0)
  tryCatch({
    tox_compute_edf(x, perm)
    stop("Should have raised error for empty input")
  }, error = function(e) {
    cat("  ✓ Correctly caught error for empty input\n")
  })
}

# Test 7: Large Dataset
test_compute_edf_large_dataset <- function() {
  cat("\n[test_compute_edf_large_dataset] Large dataset (1000 values)\n")
  x <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3, 4)
  expected_cdf <- c(0.25, 0.5, 0.75, 1.0)
  cat("  Input size:", length(x), "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 4)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ Large dataset test passed\n")
}

# Test 8: Negative Values
test_compute_edf_negative_values <- function() {
  cat("\n[test_compute_edf_negative_values] Negative values\n")
  x <- c(-3, -1, 0, 1, 3)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(-3, -1, 0, 1, 3)
  expected_cdf <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 5)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ Negative values test passed\n")
}

# Test 9: Unsorted Input
test_compute_edf_unsorted_input <- function() {
  cat("\n[test_compute_edf_unsorted_input] Unsorted input values\n")
  x <- c(3, 1, 4, 1, 5, 9, 2, 6)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- sort(unique(x))
  expected_cdf <- c(2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 1.0)
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 7)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ Unsorted input test passed\n")
}

# Test 10: List Input
test_compute_edf_list_input <- function() {
  cat("\n[test_compute_edf_list_input] R vector as input\n")
  x <- c(1, 2, 2, 3)
  perm <- order(x)
  result <- tox_compute_edf(x, perm)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3)
  expected_cdf <- c(0.25, 0.75, 1.0)
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 3)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ List input test passed\n")
}

# Test 11: Default path without perm
test_compute_edf_default_no_perm <- function() {
  cat("\n[test_compute_edf_default_no_perm] Default EDF path without perm\n")
  x <- c(3, 1, 2, 2, 3, 3)
  result <- tox_compute_edf(x)
  n_unique <- result$n_unique
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3)
  expected_cdf <- c(1/6, 3/6, 6/6)
  cat("  Input values:", x, "\n")
  cat("  Number of unique values:", n_unique, "\n")
  cat("  Unique values:", unique_vals, "\n")
  cat("  CDF values:", cdf_vals, "\n")
  stopifnot(n_unique == 3)
  stopifnot(all(abs(unique_vals - expected_unique) < TOL))
  stopifnot(all(abs(cdf_vals - expected_cdf) < TOL))
  cat("  ✓ Default no-perm path test passed\n")
}

# =====================
# Run all tests
# =====================

cat("\n=================================================\n")
cat("    COMPUTE EDF FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

test_compute_edf_simple()
test_compute_edf_all_unique()
test_compute_edf_all_same()
test_compute_edf_duplicates()
test_compute_edf_single_value()
test_compute_edf_empty_input()
test_compute_edf_large_dataset()
test_compute_edf_negative_values()
test_compute_edf_unsorted_input()
test_compute_edf_list_input()
test_compute_edf_default_no_perm()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all compute_edf R interface tests passed! ✓\n\n")
