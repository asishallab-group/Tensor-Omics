
# Set library path and compile
source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

# Test 1: Simple 3D vectors
test_euclidean_distance_3d <- function() {
  vec1 <- c(1.0, 2.0, 3.0)
  vec2 <- c(4.0, 5.0, 6.0)
  
  result <- tox_euclidean_distance(vec1, vec2)$euclidean_distance
  expected <- sqrt(sum((vec1 - vec2)^2))
  
  
  # Verify calculation
  assert_true(abs(result - expected) < 1e-12)
  
}

# Test 2: 2D vector to origin (3-4-5 triangle)
test_euclidean_distance_to_origin <- function() {
  vec1 <- c(3.0, 4.0)
  vec2 <- c(0.0, 0.0)
  
  result <- tox_euclidean_distance(vec1, vec2)$euclidean_distance
  expected <- 5.0
  
  
  # Verify 3-4-5 triangle
  assert_true(abs(result - expected) < 1e-12)
  
}

# Test 3: Identical vectors
test_euclidean_distance_identical <- function() {
  vec1 <- c(1.5, 2.7, 3.9)
  vec2 <- vec1  # identical
  
  result <- tox_euclidean_distance(vec1, vec2)$euclidean_distance
  expected <- 0.0
  
  
  # Verify zero distance
  assert_true(abs(result - expected) < 1e-15)
  
}

# Test 4: High-dimensional vectors
test_euclidean_distance_high_dimensional <- function() {
  # 100-dimensional vectors
  d <- 100
  vec1 <- 1:d  # [1, 2, 3, ..., 100]
  vec2 <- 2:(d+1)  # [2, 3, 4, ..., 101] (shift by 1)
  
  result <- tox_euclidean_distance(vec1, vec2)$euclidean_distance
  expected <- sqrt(d)  # sqrt(100 * 1^2) = 10
  
  
  # Verify high-dimensional calculation
  assert_true(abs(result - expected) < 1e-12)
  
}

# Test 5: Invalid inputs (should throw errors)
test_euclidean_distance_invalid_inputs <- function() {
  # Test different lengths
  assert_error(tox_euclidean_distance(c(1, 2), c(1, 2, 3)), "case: same length")
  
  # Test empty vectors
  assert_error(tox_euclidean_distance(numeric(0), numeric(0)), "case: cannot be empty")
  
  # Test non-numeric input
  assert_error(tox_euclidean_distance(c("a", "b"), c(1, 2)), "case: must be numeric")
}

# =====================
# Tests for distance_to_centroid
# =====================

# Test 6: Distance to centroid functionality
test_distance_to_centroid_basic <- function() {
  
  # Gene expression data (genes as columns, dimensions as rows)
  # Gene 1: [1, 0, 0] - Family 1
  # Gene 2: [0, 1, 0] - Family 1  
  # Gene 3: [3, 0, 0] - Family 2
  # Gene 4: [0, 3, 0] - Family 2
  genes <- matrix(c(1.0, 0.0, 0.0,  # Gene 1
             0.0, 1.0, 0.0,  # Gene 2
             3.0, 0.0, 0.0,  # Gene 3
             0.0, 3.0, 0.0), ncol=4L, nrow=3L, byrow=FALSE)  # Gene 4
  
  # Family centroids
  # Family 1 centroid: [0.5, 0.5, 0.0]
  # Family 2 centroid: [1.5, 1.5, 0.0]
  centroids <- matrix(c(0.5, 0.5, 0.0,  # Family 1
                 1.5, 1.5, 0.0), ncol=2L, nrow=3L, byrow=FALSE)  # Family 2
  
  # Gene-to-family mapping (1-based)
  gene_to_fam <- c(1L, 1L, 2L, 2L)
  
  result <- tox_distance_to_centroid(genes, centroids, gene_to_fam)$distance
  
  # Expected distances
  # Gene 1: [1,0,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
  # Gene 2: [0,1,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
  # Gene 3: [3,0,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
  # Gene 4: [0,3,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
  expected <- c(sqrt(0.5^2 + 0.5^2), sqrt(0.5^2 + 0.5^2), 
                sqrt(1.5^2 + 1.5^2), sqrt(1.5^2 + 1.5^2))
  
  for(i in 1:length(result)) {
    assert_true(abs(result[i] - expected[i]) < 1e-12)
  }
  
}

# Test 7: Handling invalid family indices (should return -1 for invalid genes)
test_distance_to_centroid_invalid_families <- function() {
  
  # Gene data
  genes <- matrix(c(1.0, 2.0,  # Gene 1
             3.0, 4.0,  # Gene 2
             5.0, 6.0), nrow=2, ncol=3, byrow=FALSE)  # Gene 3
  
  # Centroids
  centroids <- matrix(c(0.0, 0.0,  # Family 1
                 1.0, 1.0), nrow=2, ncol=2, byrow=FALSE)  # Family 2
  
  # Mixed family assignments: valid (1), invalid (3), no family (0)
  gene_to_fam <- c(1L, 0L, 0L)  # family 0 = no assignment
  
  result <- tox_distance_to_centroid(genes, centroids, gene_to_fam)$distance
  
  
  # Verify handling of invalid indices
  assert_true(result[1] > 0)      # Gene 1 has valid family, should have positive distance
  assert_true(result[2] == -1)    # Gene 2 has invalid family, should be -1
  assert_true(result[3] == -1)    # Gene 3 has no family, should be -1
  
}

# Test 8: Performance test with realistic genomic data size
test_distance_to_centroid_performance <- function() {
  
  n_genes <- 1000
  n_families <- 50
  d <- 20  # 20 tissue types
  
  # Generate random-like data
  set.seed(12345)
  genes <- matrix(rnorm(d * n_genes), ncol=n_genes, nrow=d)
  centroids <- matrix(rnorm(d * n_families), ncol=n_families, nrow=d)
  gene_to_fam <- sample(1:n_families, n_genes, replace = TRUE)
  
  # Time the operation
  start_time <- Sys.time()
  
  result <- tox_distance_to_centroid(genes, centroids, gene_to_fam)$distance
  
  end_time <- Sys.time()
  elapsed <- as.numeric(end_time - start_time, units = "secs")
  
  # Check results
  valid_distances <- sum(result > 0)
  
  # Verify all distances are positive
  assert_true(all(result > 0))
  assert_true(length(result) == n_genes)
  
}

# Test 9: Input validation for distance_to_centroid
test_distance_to_centroid_input_validation <- function() {
  # Test gene_to_fam length mismatch
  assert_error(tox_distance_to_centroid(c(1, 2, 3, 4), c(1, 2), c(1, 2, 3)), "wrong gene_to_fam length")
  
  # Test negative family indices (should throw error)
  assert_error(tox_distance_to_centroid(c(1, 2, 3, 4), c(1, 2), c(1, -1)), "negative family index")
}

# Test 10: Single-dimensional vectors
test_euclidean_distance_1d <- function() {
  vec1 <- c(5.0)
  vec2 <- c(2.0)
  
  result <- tox_euclidean_distance(vec1, vec2)$euclidean_distance
  expected <- 3.0
  
  
  # Verify 1D calculation
  assert_true(abs(result - expected) < 1e-12)
  
}

run_all_tests()
