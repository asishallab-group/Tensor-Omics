# Comprehensive R test suite for shift vector field (mirrors Fortran unit tests)
# Source the main functions
source("rcpp/load_tensor_omics.R")
source("rcpp/test_helpers.R")

# 1. Test correct mapping between families and genes
test_correct_family_mapping <- function() {
  expression_vectors <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow=3, ncol=5)
  family_centroids <- matrix(c(5,4,3,2,1,0, -1,-2,-3), nrow=3, ncol=3)
  gene_to_centroid <- c(2,3,1,3,1) # 1-based
  res <- compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  shift_vectors <- matrix(res, nrow=6, ncol=5)
  # Expected: rows 1:3 = centroid, rows 4:6 = shift
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  assert_true(all(dim(shift_vectors) == c(6,5)))
  assert_true(all(abs(shift_vectors - expected) < 1e-12))
}

# 2. Test for invalid family id mapping raising error
test_invalid_family_mapping <- function() {
  expression_vectors <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  family_centroids <- matrix(c(5,4,3,2,1,0, -1,-2,-3), nrow=3, ncol=3)
  gene_to_centroid <- c(3,4) # 4 is invalid
  assert_error(compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid), "expected error for invalid family")
}

# 3. Test for zero distance between paralog and centroid
test_zero_distance <- function() {
  expression_vectors <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  family_centroids <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  gene_to_centroid <- c(1,2)
  res <- compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  shift_vectors <- matrix(res, nrow=6, ncol=2)
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  assert_true(all(abs(shift_vectors - expected) < 1e-12))
}

# 4. Test for multiple genes per family centroid
test_multiple_genes_per_family <- function() {
  expression_vectors <- matrix(1:8, nrow=2, ncol=4)
  family_centroids <- matrix(c(10,20,30,40), nrow=2, ncol=2)
  gene_to_centroid <- c(1,2,1,2)
  res <- compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  shift_vectors <- matrix(res, nrow=4, ncol=4)
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  assert_true(all(abs(shift_vectors - expected) < 1e-12))
}

# 5. Test for single gene per family centroid
test_single_gene_per_family <- function() {
  expression_vectors <- matrix(1:8, nrow=2, ncol=4)
  family_centroids <- matrix(seq(10,80,10), nrow=2, ncol=4)
  gene_to_centroid <- c(1,2,3,4)
  res <- compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  shift_vectors <- matrix(res, nrow=4, ncol=4)
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  assert_true(all(abs(shift_vectors - expected) < 1e-12))
}

# 6. Test for dimension edge cases (0 genes with dimension 1 and 1 family)
test_dimension_edge_cases <- function() {
  expression_vectors <- matrix(numeric(0), nrow=1, ncol=0)
  family_centroids <- matrix(0, nrow=1, ncol=1)
  gene_to_centroid <- integer(0)
  assert_error(compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid), "expected error for empty input")
}

run_all_tests()