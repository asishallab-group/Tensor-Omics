# test_tox_jensen_shannon_test.R
source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

# Test 1: Basic compute_gene_means
test_compute_gene_means_basic <- function() {
  
  # Create simple expression matrix: 3 replicates, 4 genes
  expr <- matrix(c(
    10, 20, 30, 40,
    12, 22, 32, 42,
    14, 24, 34, 44
  ), nrow = 3, ncol = 4, byrow = TRUE)
  
  result <- tox_compute_gene_means(expr)
  
  # Verify calculations
  expected <- c(12, 22, 32, 42)  # (10+12+14)/3 = 12, etc.
  assert_true(all(abs(result - expected) < 1e-12))
  
}

# Test 2: compute_gene_means with NA values
test_compute_gene_means_with_na <- function() {
  
  expr <- matrix(c(
    10,  20,  NA,  40,
    12,  NA,  32,  42,
    14,  24,  34,  NA
  ), nrow = 3, ncol = 4, byrow = TRUE)
  
  result <- tox_compute_gene_means(expr)
  
  # Verify calculations (skip NA values)
  # Gene 1: (10+12+14)/3 = 12
  # Gene 2: (20+24)/2 = 22 (NA excluded)
  # Gene 3: (32+34)/2 = 33 (NA excluded)
  # Gene 4: (40+42)/2 = 41 (NA excluded)
  expected <- c(12, 22, 33, 41)
  assert_true(all(abs(result - expected) < 1e-12, na.rm = TRUE))
  
}

# Test 3: Basic compute_residuals
test_compute_residuals_basic <- function() {
  
  expr <- matrix(c(
    10, 20, 30,
    12, 22, 32,
    14, 24, 34
  ), nrow = 3, ncol = 3, byrow = TRUE)
  
  means <- c(12, 22, 32)  # Known means
  result <- tox_compute_residuals(expr, means)
  
  # Verify calculations: residuals = expr - means
  expected <- matrix(c(
    10-12, 20-22, 30-32,
    12-12, 22-22, 32-32,
    14-12, 24-22, 34-32
  ), nrow = 3, ncol = 3, byrow = TRUE)
  
  assert_true(all(abs(result - expected) < 1e-12))
  
}

# Test 4: Basic pool_means
test_pool_means_basic <- function() {
  
  mean_S1 <- c(10, 12, 14, NA, 18)
  mean_S2 <- c(20, 22, NA, 26, 28)
  n_points <- 3
  
  result <- tox_pool_means(mean_S1, mean_S2, n_points)
  
  # Verify
  # Non-NA values: 10,12,14,18 from S1 and 20,22,26,28 from S2 = 8 total
  assert_true(result$n_pool == 8)
  assert_true(length(result$x_star) == n_points)
  
}

# Test 5: Basic construct_neighborhoods
test_construct_neighborhoods_basic <- function() {
  
  # Small test case
  n_points <- 3
  x_star <- c(10, 20, 30)  # Reference points
  
  # 5 genes, 2 replicates
  mean_S <- c(8, 12, 18, 22, 28)
  resid_S <- matrix(c(
    1, -1, 2, -2, 3,
    -1, 1, -2, 2, -3
  ), nrow = 2, ncol = 5, byrow = TRUE)
  
  n_pool <- 10  # Arbitrary value for test
  
  result <- tox_construct_neighborhoods(x_star, n_pool, mean_S, resid_S, n_pool)
  force(result)
}

# Test 7: Error handling test
test_error_handling <- function() {
  
  # Test with invalid dimensions
  expr <- matrix(1:6, nrow = 2, ncol = 3)
  means_wrong <- c(1, 2)  # Wrong length
  
  assert_error(tox_compute_residuals(expr, means_wrong), "Test failed - should have caught dimension mismatch")
  
  # Test with invalid n_points
  mean_S1 <- c(1, 2, 3)
  mean_S2 <- c(4, 5, 6)
  
  assert_error(tox_pool_means(mean_S1, mean_S2, n_points = 0), "Test failed - should have caught invalid n_points")
  
}

# Test with explicit neighborhood size
test_construct_neighborhoods_explicit <- function() {
  
  n_points <- 3
  x_star <- c(10, 20, 30)
  mean_S <- c(8, 12, 18, 22, 28, 8, 12, 18, 22, 28, 8, 12, 18, 22, 28, 8, 12, 18, 22, 28)
  resid_S <- matrix(c(
    1, -1, 2, -2, 3, 1, -1, 2, -2, 3, 1, -1, 2, -2, 3, 1, -1, 2, -2, 3,
    -1, 1, -2, 2, -3, -1, 1, -2, 2, -3, -1, 1, -2, 2, -3, -1, 1, -2, 2, -3
  ), nrow = 2, ncol = 20, byrow = TRUE)
  n_pool <- 20
  explicit_size <- 5
  
  # Test with explicit size
  result <- tox_construct_neighborhoods(x_star, n_pool, mean_S, resid_S, explicit_size)
  
  assert_true(dim(result$neighborhood_residuals)[2] == explicit_size)  # Should use the explicit size
  
}

run_all_tests()
