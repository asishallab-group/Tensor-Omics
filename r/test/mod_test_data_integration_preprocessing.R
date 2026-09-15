# test_tox_jensen_shannon_test.R
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# Test 1: Basic compute_gene_means
test_compute_gene_means_basic <- function() {
  
  # Create simple expression matrix: 3 replicates, 4 genes
  expr <- matrix(c(
    10, 20, 30, 40,
    12, 22, 32, 42,
    14, 24, 34, 44
  ), nrow = 3, ncol = 4, byrow = TRUE)
  
  result <- compute_gene_means(expr)
  
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
  
  result <- compute_gene_means(expr)
  
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
  result <- compute_residuals(expr, means)
  
  # Verify calculations: residuals = expr - means
  expected <- matrix(c(
    10-12, 20-22, 30-32,
    12-12, 22-22, 32-32,
    14-12, 24-22, 34-32
  ), nrow = 3, ncol = 3, byrow = TRUE)
  
  assert_true(all(abs(result - expected) < 1e-12))
  
}

# Test 4: Basic pool_study_means
test_pool_study_means_basic <- function() {

  mean_S1 <- c(10, 12, 14, NA, 18)
  mean_S2 <- c(20, 22, NA, 26, 28)
  n_points <- 3

  result <- pool_study_means(mean_S1, mean_S2, n_points)

  # Verify
  # Non-NA values: 10,12,14,18 from S1 and 20,22,26,28 from S2 = 8 total
  assert_true(result$n_pool == 8)
  assert_true(length(result$x_star) == n_points)

}

# Test 4b: pool_means (the plain, auto-sorting tier) directly.
# n_points=647, pool_size=1000, i_point=216 (1-based) is a confirmed case where 125's
# percentage-round-trip quantile formula and a naive direct-fraction formula round to different
# doubles, which floor() to different ranks -- a genuine index flip, not epsilon noise (see the
# Fortran regression test test_pool_means_matches_125_quantile_rounding in
# test/mod_test_data_integration.F90). With a sequential pool [1.0, ..., 1000.0], the interpolated
# x_star equals the rank itself, so the fixed formula's expected value is exactly
# 333.99999999999994, not the direct-fraction formula's (wrong) 334.0.
test_pool_means_matches_125_quantile_rounding <- function() {

  pooled_means <- as.numeric(1:1000)
  n_points <- 647

  result <- pool_means(pooled_means, n_points)

  assert_true(result$n_pool == 1000)
  assert_true(length(result$x_star) == n_points)
  assert_true(abs(result$x_star[216] - 333.99999999999994) < 1e-12)

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
  
  # the routine takes the actual neighbour count now, not a desired size
  n_neighbors <- calc_neighborhood_size(n_pool, length(x_star), mean_S, n_pool)
  result <- construct_neighborhoods(x_star, mean_S, resid_S, n_neighbors)
  force(result)
}

# Test: construct_neighborhoods_ranged / _expert basic two-point case, hand-computed from a
# sorted mean_S with no ties.
test_construct_neighborhoods_ranged_basic <- function() {

  x_star <- c(2.0, 10.0)
  mean_S <- c(1.0, 2.5, 9.0, 10.5, 20.0)
  mean_S_perm <- as.integer(c(1, 2, 3, 4, 5))  # already ascending, no ties
  n_neighbors <- 2L

  # The plain entry point sorts mean_S itself; mean_S_perm is unused here.
  plain <- function(x_star, mean_S, mean_S_perm, n_neighbors) {
    construct_neighborhoods_ranged(x_star, mean_S, n_neighbors)
  }

  for (func in list(construct_neighborhoods_ranged_expert, plain)) {
    result <- func(x_star, mean_S, mean_S_perm, n_neighbors)

    # x_star[1]=2.0: distances to [1,2.5,9,10.5,20] are [1,0.5,7,8.5,18] -> nearest are gene 2, then gene 1
    assert_true(all(result$neighborhood_indices[, 1] == c(2, 1)), "indices for point 1")
    assert_true(all(result$neighborhood_range[, 1] == c(1, 2)), "range for point 1")

    # x_star[2]=10.0: distances to [1,2.5,9,10.5,20] are [9,7.5,1,0.5,10] -> nearest are gene 4, then gene 3
    assert_true(all(result$neighborhood_indices[, 2] == c(4, 3)), "indices for point 2")
    assert_true(all(result$neighborhood_range[, 2] == c(3, 4)), "range for point 2")
  }
}

# Test: three genes tie on the same mean at the boundary of a single-neighbor neighborhood: only
# one of them is picked as the neighbor, but neighborhood_range must still be extended to span
# all three, because a later admissibility gate reasons about the range, not just the picked
# neighbor.
test_construct_neighborhoods_ranged_tie_extends_range <- function() {

  x_star <- c(3.0)
  mean_S <- c(1.0, 3.0, 3.0, 3.0, 5.0)
  mean_S_perm <- as.integer(c(1, 2, 3, 4, 5))  # already ascending
  n_neighbors <- 1L

  result <- construct_neighborhoods_ranged_expert(x_star, mean_S, mean_S_perm, n_neighbors)

  # Binary search lands on the first "3.0" (gene 2); it is closer than gene 1 (distance 0 vs 2),
  # so gene 2 is the sole neighbor -- but genes 3 and 4 share its mean and must extend the range.
  assert_true(all(result$neighborhood_indices[, 1] == c(2)), "picked neighbor")
  assert_true(all(result$neighborhood_range[, 1] == c(2, 4)), "range must include all ties")
}

# Test: when every gene's mean is NA, the reference point's own value is irrelevant: the routine
# falls back to the first min(n_genes_S, n_neighbors) positions.
test_construct_neighborhoods_ranged_all_nan_means <- function() {

  x_star <- c(5.0)
  mean_S <- c(NA_real_, NA_real_, NA_real_)
  mean_S_perm <- as.integer(c(1, 2, 3))
  n_neighbors <- 2L

  result <- construct_neighborhoods_ranged_expert(x_star, mean_S, mean_S_perm, n_neighbors)

  assert_true(all(result$neighborhood_indices[, 1] == c(1, 2)), "fallback indices")
  assert_true(all(result$neighborhood_range[, 1] == c(1, 2)), "fallback range")
}

# Test: n_neighbors exceeding the number of genes is a documented, deliberately-preserved
# limitation: the leftover slots are filled with out-of-range sentinel indices
# (n_genes_S+1, n_genes_S+2, ...) rather than left undefined.
test_construct_neighborhoods_ranged_fewer_genes_than_neighbors <- function() {

  x_star <- c(2.0)
  mean_S <- c(1.0, 2.0)
  mean_S_perm <- as.integer(c(1, 2))
  n_neighbors <- 3L

  result <- construct_neighborhoods_ranged_expert(x_star, mean_S, mean_S_perm, n_neighbors)

  # gene 2 (distance 0), then gene 1 (distance 1), then both pointers are exhausted -> sentinel 3
  assert_true(all(result$neighborhood_indices[, 1] == c(2, 1, 3)), "indices")
  assert_true(all(result$neighborhood_range[, 1] == c(1, 2)), "range")
}

# Test 7: Error handling test
test_error_handling <- function() {
  
  # Test with invalid dimensions
  expr <- matrix(1:6, nrow = 2, ncol = 3)
  means_wrong <- c(1, 2)  # Wrong length
  
  assert_error(compute_residuals(expr, means_wrong), "Test failed - should have caught dimension mismatch")
  
  # Test with invalid n_points
  mean_S1 <- c(1, 2, 3)
  mean_S2 <- c(4, 5, 6)
  
  assert_error(pool_study_means(mean_S1, mean_S2, n_points = 0), "Test failed - should have caught invalid n_points", ERR_EMPTY_INPUT)
  
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
  n_neighbors <- calc_neighborhood_size(n_pool, n_points, mean_S, explicit_size)
  result <- construct_neighborhoods(x_star, mean_S, resid_S, n_neighbors)
  
  assert_true(dim(result$neighborhood_residuals)[2] == explicit_size)  # Should use the explicit size
  
}

run_all_tests()
