
# Set library path and compile
source("r/load_tensor_omics.R")
source("r/test_helpers.R")


TOL <- 1e-12

# Test 1: Simple EDF Test
test_compute_edf_simple <- function() {
  x <- c(1, 2, 2, 3, 3, 3)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3)
  expected_cdf <- c(1/6, 3/6, 6/6)
  assert_true(n_unique == 3)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

# Test 2: All Unique Values
test_compute_edf_all_unique <- function() {
  x <- c(1, 2, 3, 4, 5)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  unique_vals <- result$unique_values
  cdf_vals <- result$cdf_values
  n_unique <- length(unique_vals)
  expected_unique <- c(1, 2, 3, 4, 5)
  expected_cdf <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  assert_true(n_unique == 5)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

# Test 3: All Same Values
test_compute_edf_all_same <- function() {
  x <- c(5, 5, 5, 5)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  assert_true(n_unique == 1)
  assert_true(length(unique_vals) == 1)
  assert_true(length(cdf_vals) == 1)
  assert_true(abs(unique_vals[1] - 5) < TOL)
  assert_true(abs(cdf_vals[1] - 1.0) < TOL)
}

test_compute_edf_duplicates <- function() {
  x <- c(1, 1, 2, 3, 3, 3, 4)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  assert_true(n_unique == 4)
  assert_true(length(unique_vals) == 4)
  assert_true(length(cdf_vals) == 4)
  expected_unique <- c(1, 2, 3, 4)
  expected_cdf <- c(2/7, 3/7, 6/7, 1.0)
  assert_true(n_unique == 4)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

# Test 5: Single Value
test_compute_edf_single_value <- function() {
  x <- c(42)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  assert_true(n_unique == 1)
  assert_true(abs(unique_vals[1] - 42) < TOL)
  assert_true(abs(cdf_vals[1] - 1.0) < TOL)
}

# Test 6: Empty Input
test_compute_edf_empty_input <- function() {
  x <- numeric(0)
  perm <- integer(0)
  assert_error(compute_edf_expert(x, perm), "Should have raised error for empty input", ERR_EMPTY_INPUT)
}

# Test 7: Large Dataset
test_compute_edf_large_dataset <- function() {
  x <- c(rep(1, 250), rep(2, 250), rep(3, 250), rep(4, 250))
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3, 4)
  expected_cdf <- c(0.25, 0.5, 0.75, 1.0)
  assert_true(n_unique == 4)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

# Test 8: Negative Values
test_compute_edf_negative_values <- function() {
  x <- c(-3, -1, 0, 1, 3)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(-3, -1, 0, 1, 3)
  expected_cdf <- c(0.2, 0.4, 0.6, 0.8, 1.0)
  assert_true(n_unique == 5)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

# Test 9: Unsorted Input
test_compute_edf_unsorted_input <- function() {
  x <- c(3, 1, 4, 1, 5, 9, 2, 6)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- sort(unique(x))
  expected_cdf <- c(2/8, 3/8, 4/8, 5/8, 6/8, 7/8, 1.0)
  assert_true(n_unique == 7)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

# Test 10: List Input
test_compute_edf_list_input <- function() {
  x <- c(1, 2, 2, 3)
  perm <- order(x)
  result <- compute_edf_expert(x, perm)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3)
  expected_cdf <- c(0.25, 0.75, 1.0)
  assert_true(n_unique == 3)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

# Test 11: Default path without perm
test_compute_edf_default_no_perm <- function() {
  x <- c(3, 1, 2, 2, 3, 3)
  result <- compute_edf(x)
  n_unique <- length(result$unique_values)
  unique_vals <- result$unique_values[seq_len(n_unique)]
  cdf_vals <- result$cdf_values[seq_len(n_unique)]
  expected_unique <- c(1, 2, 3)
  expected_cdf <- c(1/6, 3/6, 6/6)
  assert_true(n_unique == 3)
  assert_true(all(abs(unique_vals - expected_unique) < TOL))
  assert_true(all(abs(cdf_vals - expected_cdf) < TOL))
}

run_all_tests()
