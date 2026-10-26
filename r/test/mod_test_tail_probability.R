# =====================
# Comprehensive R test suite for outlier detection
# Uses tensor_omics.R wrapper functions
# =====================

# Source the main functions
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# The entry point sorts the distribution itself; only the clamping of invalid (negative)
# RDIs is left here, so the assertions below stay about the numerics.
scaled_distance_tail_probability <- function(distribution, c_const) {
  if (length(distribution) == 0) return(numeric(0))
  sorted_rdi <- distribution
  sorted_rdi[sorted_rdi < 0] <- 0
  compute_scaled_distance_tail_probability(distribution, sorted_rdi, c_const)
}

# =====================
# Tests for compute_scaled_distance_tail_probability
# =====================

# Test 1: Basic empirical tail probability calculation
test_tail_probability_basic <- function() {
  distribution <- c(0.5, 1.2, 0.8, 0.3)
  c_const <- 1.0

  tail_probability <- scaled_distance_tail_probability(distribution, c_const)

  # Verify tail_probability is within [0, 1]
  assert_true(all(tail_probability >= 0 & tail_probability <= 1))

}

# Test 2: All zeros distribution
test_tail_probability_all_zeros <- function() {
  distribution <- c(0, 0, 0, 0, 0)
  c_const <- 1.0

  tail_probability <- scaled_distance_tail_probability(distribution, c_const)

  # Verify all tail_probability values are 1
  assert_true(all(tail_probability == 1))

}

# Test 3: Negative values in distribution
test_tail_probability_negative_values <- function() {
  distribution <- c(-0.5, 1.2, -0.8, 0.3)
  c_const <- 1.0

  tail_probability <- scaled_distance_tail_probability(distribution, c_const)

  # Verify tail_probability for negative values are 1
  assert_true(all(tail_probability[distribution < 0] == 1))

}

# Test 4: Large distribution
test_tail_probability_large_distribution <- function() {
  set.seed(42)  # For reproducibility
  distribution <- runif(1000, 0, 10)  # Large distribution
  c_const <- 1.0

  tail_probability <- scaled_distance_tail_probability(distribution, c_const)

  # Verify tail_probability is within [0, 1]
  assert_true(all(tail_probability >= 0 & tail_probability <= 1))

}

# =====================
# Run all tests
# =====================

run_all_tests()
