# =====================
# Comprehensive R test suite for outlier detection 
# Uses tensoromics_functions.R wrapper functions
# =====================

# Source the main functions
source("rcpp/load_tensor_omics.R")
source("rcpp/test_helpers.R")

# The generated interface exposes the Fortran routine as it is; the sort prep the old
# wrapper did is here, so the assertions below stay about the numerics.
empirical_p_values <- function(distribution, c_const) {
  if (length(distribution) == 0) return(numeric(0))
  sorted_rdi <- distribution
  sorted_rdi[sorted_rdi < 0] <- 0
  perm <- order(sorted_rdi)
  compute_empirical_p_values(distribution, sorted_rdi, perm, c_const)
}

# =====================
# Tests for compute_empirical_p_values
# =====================

# Test 1: Basic empirical p-values calculation
test_empirical_p_values_basic <- function() {
  distribution <- c(0.5, 1.2, 0.8, 0.3)
  c_const <- 1.0

  p_values <- empirical_p_values(distribution, c_const)

  # Verify p-values are within [0, 1]
  assert_true(all(p_values >= 0 & p_values <= 1))

}

# Test 2: All zeros distribution
test_empirical_p_values_all_zeros <- function() {
  distribution <- c(0, 0, 0, 0, 0)
  c_const <- 1.0

  p_values <- empirical_p_values(distribution, c_const)

  # Verify all p-values are 1
  assert_true(all(p_values == 1))

}

# Test 3: Negative values in distribution
test_empirical_p_values_negative_values <- function() {
  distribution <- c(-0.5, 1.2, -0.8, 0.3)
  c_const <- 1.0

  p_values <- empirical_p_values(distribution, c_const)

  # Verify p-values for negative values are 1
  assert_true(all(p_values[distribution < 0] == 1))

}

# Test 4: Large distribution
test_empirical_p_values_large_distribution <- function() {
  set.seed(42)  # For reproducibility
  distribution <- runif(1000, 0, 10)  # Large distribution
  c_const <- 1.0

  p_values <- empirical_p_values(distribution, c_const)

  # Verify p-values are within [0, 1]
  assert_true(all(p_values >= 0 & p_values <= 1))

}

# =====================
# Run all tests
# =====================

run_all_tests()
