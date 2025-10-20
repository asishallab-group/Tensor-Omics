# Comprehensive R test suite for shift vector field (mirrors Fortran unit tests)
# Source the main functions
source("r/tensoromics_functions.R")

# 1. Test correct mapping between families and genes
test_correct_family_mapping <- function() {
  expression_vectors <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), nrow=3, ncol=5)
  family_centroids <- matrix(c(5,4,3,2,1,0, -1,-2,-3), nrow=3, ncol=3)
  gene_to_centroid <- c(2,3,1,3,1) # 1-based
  shift_vectors <- tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  # Expected: rows 1:3 = centroid, rows 4:6 = shift
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  stopifnot(all(dim(shift_vectors) == c(6,5)))
  stopifnot(all(abs(shift_vectors - expected) < 1e-12))
  cat("test_correct_family_mapping passed\n")
}

# 2. Test for invalid family id mapping raising error
test_invalid_family_mapping <- function() {
  expression_vectors <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  family_centroids <- matrix(c(5,4,3,2,1,0, -1,-2,-3), nrow=3, ncol=3)
  gene_to_centroid <- c(3,4) # 4 is invalid
  error_caught <- FALSE
  tryCatch({
    tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  }, error = function(e) {
    error_caught <<- TRUE
    stopifnot(grepl("Invalid input", e$message) || grepl("invalid", e$message))
  })
  stopifnot(error_caught)
  cat("test_invalid_family_mapping passed\n")
}

# 3. Test for zero distance between paralog and centroid
test_zero_distance <- function() {
  expression_vectors <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  family_centroids <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
  gene_to_centroid <- c(1,2)
  shift_vectors <- tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  stopifnot(all(abs(shift_vectors - expected) < 1e-12))
  cat("test_zero_distance passed\n")
}

# 4. Test for multiple genes per family centroid
test_multiple_genes_per_family <- function() {
  expression_vectors <- matrix(1:8, nrow=2, ncol=4)
  family_centroids <- matrix(c(10,20,30,40), nrow=2, ncol=2)
  gene_to_centroid <- c(1,2,1,2)
  shift_vectors <- tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  stopifnot(all(abs(shift_vectors - expected) < 1e-12))
  cat("test_multiple_genes_per_family passed\n")
}

# 5. Test for single gene per family centroid
test_single_gene_per_family <- function() {
  expression_vectors <- matrix(1:8, nrow=2, ncol=4)
  family_centroids <- matrix(seq(10,80,10), nrow=2, ncol=4)
  gene_to_centroid <- c(1,2,3,4)
  shift_vectors <- tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  expected_centroids <- sapply(gene_to_centroid, function(idx) family_centroids[, idx])
  expected_shifts <- expression_vectors - expected_centroids
  expected <- rbind(expected_centroids, expected_shifts)
  stopifnot(all(abs(shift_vectors - expected) < 1e-12))
  cat("test_single_gene_per_family passed\n")
}

# 6. Test for dimension edge cases (0 genes with dimension 1 and 1 family)
test_dimension_edge_cases <- function() {
  expression_vectors <- matrix(numeric(0), nrow=1, ncol=0)
  family_centroids <- matrix(0, nrow=1, ncol=1)
  gene_to_centroid <- integer(0)
  error_caught <- FALSE
  tryCatch({
    tox_compute_shift_vector_field(expression_vectors, family_centroids, gene_to_centroid)
  }, error = function(e) {
    error_caught <<- TRUE
    stopifnot(grepl("empty input", tolower(e$message)) || grepl("zero", tolower(e$message)))
  })
  stopifnot(error_caught)
  cat("test_dimension_edge_cases passed\n")
}

# Run all tests
cat("=================================================\n")
cat("    SHIFT VECTOR FIELD FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

test_correct_family_mapping()
test_invalid_family_mapping()
test_zero_distance()
test_multiple_genes_per_family()
test_single_gene_per_family()
test_dimension_edge_cases()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all shift vector R interface tests passed successfully!\n")