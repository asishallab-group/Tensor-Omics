# =====================
# Comprehensive R test suite for Euclidean distance functions (mirrors Fortran unit tests)
# Uses tensoromics_functions.R wrapper functions
# =====================

# Source the main functions
source("r/tensoromics_functions.R")

cat("=== Testing euclidean distance R wrapper functions ===\n")
cat("Based on Fortran test suite with comprehensive test coverage\n")

# =====================
# Tests for euclidean_distance
# =====================

# Test 1: Simple 3D vectors
test_euclidean_distance_3d <- function() {
  cat("\n[test_euclidean_distance_3d] Simple 3D vectors test\n")
  vec1 <- c(1.0, 2.0, 3.0)
  vec2 <- c(4.0, 5.0, 6.0)
  
  result <- tox_euclidean_distance(vec1, vec2)
  expected <- sqrt(sum((vec1 - vec2)^2))
  
  cat("  Vectors [1,2,3] vs [4,5,6]\n")
  cat("  Result:", result, "\n")
  cat("  Expected:", expected, "\n")
  
  # Verify calculation
  stopifnot(abs(result - expected) < 1e-12)
  
  cat("Simple 3D vectors test passed ✓\n")
}

# Test 2: 2D vector to origin (3-4-5 triangle)
test_euclidean_distance_to_origin <- function() {
  cat("\n[test_euclidean_distance_to_origin] 2D vector to origin test\n")
  vec1 <- c(3.0, 4.0)
  vec2 <- c(0.0, 0.0)
  
  result <- tox_euclidean_distance(vec1, vec2)
  expected <- 5.0
  
  cat("  Vector [3,4] to origin\n")
  cat("  Result:", result, "\n")
  cat("  Expected:", expected, "\n")
  
  # Verify 3-4-5 triangle
  stopifnot(abs(result - expected) < 1e-12)
  
  cat("2D vector to origin test passed ✓\n")
}

# Test 3: Identical vectors
test_euclidean_distance_identical <- function() {
  cat("\n[test_euclidean_distance_identical] Identical vectors test\n")
  vec1 <- c(1.5, 2.7, 3.9)
  vec2 <- vec1  # identical
  
  result <- tox_euclidean_distance(vec1, vec2)
  expected <- 0.0
  
  cat("  Identical vectors\n")
  cat("  Result:", result, "\n")
  cat("  Expected:", expected, "\n")
  
  # Verify zero distance
  stopifnot(abs(result - expected) < 1e-15)
  
  cat("Identical vectors test passed ✓\n")
}

# Test 4: High-dimensional vectors
test_euclidean_distance_high_dimensional <- function() {
  cat("\n[test_euclidean_distance_high_dimensional] High-dimensional test\n")
  # 100-dimensional vectors
  d <- 100
  vec1 <- 1:d  # [1, 2, 3, ..., 100]
  vec2 <- 2:(d+1)  # [2, 3, 4, ..., 101] (shift by 1)
  
  result <- tox_euclidean_distance(vec1, vec2)
  expected <- sqrt(d)  # sqrt(100 * 1^2) = 10
  
  cat("  100D vectors with unit shift\n")
  cat("  Result:", result, "\n")
  cat("  Expected:", expected, "\n")
  
  # Verify high-dimensional calculation
  stopifnot(abs(result - expected) < 1e-12)
  
  cat("High-dimensional test passed ✓\n")
}

# Test 5: Invalid inputs (should throw errors)
test_euclidean_distance_invalid_inputs <- function() {
  cat("\n[test_euclidean_distance_invalid_inputs] Invalid inputs test\n")
  
  # Test different lengths
  error_caught1 <- FALSE
  tryCatch({
    tox_euclidean_distance(c(1, 2), c(1, 2, 3))
  }, error = function(e) {
    error_caught1 <<- TRUE
    stopifnot(grepl("same length", e$message))
  })
  stopifnot(error_caught1)
  
  # Test empty vectors
  error_caught2 <- FALSE
  tryCatch({
    tox_euclidean_distance(numeric(0), numeric(0))
  }, error = function(e) {
    error_caught2 <<- TRUE
    stopifnot(grepl("cannot be empty", e$message))
  })
  stopifnot(error_caught2)
  
  # Test non-numeric input
  error_caught3 <- FALSE
  tryCatch({
    tox_euclidean_distance(c("a", "b"), c(1, 2))
  }, error = function(e) {
    error_caught3 <<- TRUE
    stopifnot(grepl("must be numeric", e$message))
  })
  stopifnot(error_caught3)
  
  cat("Invalid inputs test passed ✓\n")
}

# =====================
# Tests for distance_to_centroid
# =====================

# Test 6: Distance to centroid functionality
test_distance_to_centroid_basic <- function() {
  cat("\n[test_distance_to_centroid_basic] Basic distance to centroid test\n")
  
  # Setup test data
  d <- 3
  
  # Gene expression data (genes as columns, dimensions as rows)
  # Gene 1: [1, 0, 0] - Family 1
  # Gene 2: [0, 1, 0] - Family 1  
  # Gene 3: [3, 0, 0] - Family 2
  # Gene 4: [0, 3, 0] - Family 2
  genes <- c(1.0, 0.0, 0.0,  # Gene 1
             0.0, 1.0, 0.0,  # Gene 2
             3.0, 0.0, 0.0,  # Gene 3
             0.0, 3.0, 0.0)  # Gene 4
  
  # Family centroids
  # Family 1 centroid: [0.5, 0.5, 0.0]
  # Family 2 centroid: [1.5, 1.5, 0.0]
  centroids <- c(0.5, 0.5, 0.0,  # Family 1
                 1.5, 1.5, 0.0)  # Family 2
  
  # Gene-to-family mapping (1-based)
  gene_to_fam <- c(1, 1, 2, 2)
  
  result <- tox_distance_to_centroid(genes, centroids, gene_to_fam, d)
  
  # Expected distances
  # Gene 1: [1,0,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
  # Gene 2: [0,1,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
  # Gene 3: [3,0,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
  # Gene 4: [0,3,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
  expected <- c(sqrt(0.5^2 + 0.5^2), sqrt(0.5^2 + 0.5^2), 
                sqrt(1.5^2 + 1.5^2), sqrt(1.5^2 + 1.5^2))
  
  cat("  Distance to centroid results:\n")
  for(i in 1:length(result)) {
    cat(sprintf("    Gene %d: Result=%.6f, Expected=%.6f\n", 
                i, result[i], expected[i]))
    stopifnot(abs(result[i] - expected[i]) < 1e-12)
  }
  
  cat("Basic distance to centroid test passed ✓\n")
}

# Test 7: Handling invalid family indices (should return -1 for invalid genes)
test_distance_to_centroid_invalid_families <- function() {
  cat("\n[test_distance_to_centroid_invalid_families] Invalid family indices test\n")
  
  d <- 2
  
  # Gene data
  genes <- c(1.0, 2.0,  # Gene 1
             3.0, 4.0,  # Gene 2
             5.0, 6.0)  # Gene 3
  
  # Centroids
  centroids <- c(0.0, 0.0,  # Family 1
                 1.0, 1.0)  # Family 2
  
  # Mixed family assignments: valid (1), invalid (3), no family (0)
  gene_to_fam <- c(1, 3, 0)  # Family 3 doesn't exist, family 0 = no assignment
  
  result <- tox_distance_to_centroid(genes, centroids, gene_to_fam, d)
  
  cat("  Gene 1 (family 1):", result[1], "(should be valid distance)\n")
  cat("  Gene 2 (family 3):", result[2], "(should be -1, invalid family)\n")
  cat("  Gene 3 (family 0):", result[3], "(should be -1, no family)\n")
  
  # Verify handling of invalid indices
  stopifnot(result[1] > 0)      # Gene 1 has valid family, should have positive distance
  stopifnot(result[2] == -1)    # Gene 2 has invalid family, should be -1
  stopifnot(result[3] == -1)    # Gene 3 has no family, should be -1
  
  cat("Invalid family indices test passed ✓\n")
}

# Test 8: Performance test with realistic genomic data size
test_distance_to_centroid_performance <- function() {
  cat("\n[test_distance_to_centroid_performance] Performance test\n")
  
  n_genes <- 1000
  n_families <- 50
  d <- 20  # 20 tissue types
  
  cat(sprintf("  Testing with %d genes, %d families, %d dimensions\n", 
              n_genes, n_families, d))
  
  # Generate random-like data
  set.seed(12345)
  genes <- rnorm(d * n_genes)
  centroids <- rnorm(d * n_families)
  gene_to_fam <- sample(1:n_families, n_genes, replace = TRUE)
  
  # Time the operation
  start_time <- Sys.time()
  
  result <- tox_distance_to_centroid(genes, centroids, gene_to_fam, d)
  
  end_time <- Sys.time()
  elapsed <- as.numeric(end_time - start_time, units = "secs")
  
  # Check results
  valid_distances <- sum(result > 0)
  cat(sprintf("  Completed in %.6f seconds\n", elapsed))
  cat(sprintf("  Valid distances computed: %d/%d\n", valid_distances, n_genes))
  cat(sprintf("  Mean distance: %.6f\n", mean(result)))
  
  # Verify all distances are positive
  stopifnot(all(result > 0))
  stopifnot(length(result) == n_genes)
  
  cat("Performance test passed ✓\n")
}

# Test 9: Input validation for distance_to_centroid
test_distance_to_centroid_input_validation <- function() {
  cat("\n[test_distance_to_centroid_input_validation] Input validation test\n")
  
  # Test dimension mismatch
  error_caught1 <- FALSE
  tryCatch({
    tox_distance_to_centroid(c(1, 2, 3), c(1, 2), c(1), 2)  # genes not divisible by d
  }, error = function(e) {
    error_caught1 <<- TRUE
    stopifnot(grepl("divisible by d", e$message))
  })
  stopifnot(error_caught1)
  
  # Test gene_to_fam length mismatch
  error_caught2 <- FALSE
  tryCatch({
    tox_distance_to_centroid(c(1, 2, 3, 4), c(1, 2), c(1, 2, 3), 2)  # wrong gene_to_fam length
  }, error = function(e) {
    error_caught2 <<- TRUE
    stopifnot(grepl("equal number of genes", e$message))
  })
  stopifnot(error_caught2)
  
  # Test negative family indices (should throw error)
  error_caught3 <- FALSE
  tryCatch({
    tox_distance_to_centroid(c(1, 2, 3, 4), c(1, 2), c(1, -1), 2)  # negative family index
  }, error = function(e) {
    error_caught3 <<- TRUE
    stopifnot(grepl("must be between 0 and", e$message))
  })
  stopifnot(error_caught3)
  
  cat("Input validation test passed ✓\n")
}

# Test 10: Single-dimensional vectors
test_euclidean_distance_1d <- function() {
  cat("\n[test_euclidean_distance_1d] Single-dimensional vectors test\n")
  vec1 <- c(5.0)
  vec2 <- c(2.0)
  
  result <- tox_euclidean_distance(vec1, vec2)
  expected <- 3.0
  
  cat("  1D vectors: 5.0 vs 2.0\n")
  cat("  Result:", result, "\n")
  cat("  Expected:", expected, "\n")
  
  # Verify 1D calculation
  stopifnot(abs(result - expected) < 1e-12)
  
  cat("Single-dimensional vectors test passed ✓\n")
}

# =====================
# Run all tests
# =====================

cat("\n=================================================\n")
cat("    EUCLIDEAN DISTANCE FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

# euclidean_distance tests
test_euclidean_distance_3d()
test_euclidean_distance_to_origin()
test_euclidean_distance_identical()
test_euclidean_distance_high_dimensional()
test_euclidean_distance_invalid_inputs()
test_euclidean_distance_1d()

# distance_to_centroid tests
test_distance_to_centroid_basic()
test_distance_to_centroid_invalid_families()
test_distance_to_centroid_performance()
test_distance_to_centroid_input_validation()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all euclidean distance R interface tests passed! ✓\n")
