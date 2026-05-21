# =====================
# Comprehensive R test suite for outlier detection 
# Uses tensoromics_functions.R wrapper functions
# =====================

# Source the main functions
source("rcpp/tensoromics_functions.R")

cat("=== Testing outlier detection R wrapper functions ===\n")
cat("Based on Fortran test suite with comprehensive test coverage\n")

# =====================
# Tests for compute_empirical_p_values
# =====================

# Test 1: Basic empirical p-values calculation
test_empirical_p_values_basic <- function() {
  cat("\n[test_empirical_p_values_basic] Basic empirical p-values calculation test\n")
  distribution <- c(0.5, 1.2, 0.8, 0.3)
  c_const <- 1.0

  p_values <- compute_empirical_p_values(distribution, c_const)

  # Verify p-values are within [0, 1]
  stopifnot(all(p_values >= 0 & p_values <= 1))

  cat("Basic empirical p-values calculation test passed ✓\n")
}

# Test 2: All zeros distribution
test_empirical_p_values_all_zeros <- function() {
  cat("\n[test_empirical_p_values_all_zeros] All zeros distribution test\n")
  distribution <- c(0, 0, 0, 0, 0)
  c_const <- 1.0

  p_values <- compute_empirical_p_values(distribution, c_const)

  # Verify all p-values are 1
  stopifnot(all(p_values == 1))

  cat("All zeros distribution test passed ✓\n")
}

# Test 3: Negative values in distribution
test_empirical_p_values_negative_values <- function() {
  cat("\n[test_empirical_p_values_negative_values] Negative values in distribution test\n")
  distribution <- c(-0.5, 1.2, -0.8, 0.3)
  c_const <- 1.0

  p_values <- compute_empirical_p_values(distribution, c_const)

  # Verify p-values for negative values are 1
  stopifnot(all(p_values[distribution < 0] == 1))

  cat("Negative values in distribution test passed ✓\n")
}

# Test 4: Large distribution
test_empirical_p_values_large_distribution <- function() {
  cat("\n[test_empirical_p_values_large_distribution] Large distribution test\n")
  set.seed(42)  # For reproducibility
  distribution <- runif(1000, 0, 10)  # Large distribution
  c_const <- 1.0

  p_values <- compute_empirical_p_values(distribution, c_const)

  # Verify p-values are within [0, 1]
  stopifnot(all(p_values >= 0 & p_values <= 1))

  cat("Large distribution test passed ✓\n")
}

# =====================
# Run all tests
# =====================

cat("\n=================================================\n")
cat("    EMPIRICAL P VALUE TESTS\n")
cat("=================================================\n\n")

# compute_empirical_p_values tests
test_empirical_p_values_basic()
test_empirical_p_values_all_zeros()
test_empirical_p_values_negative_values()
test_empirical_p_values_large_distribution()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all outlier detection R interface tests passed! ✓\n")
