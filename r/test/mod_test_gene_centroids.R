# Comprehensive R test suite for gene centroids binding functions

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# 1. Basic functionality in 'all' mode
test_basic_all_mode <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  ortholog_set <- rep(FALSE, n_genes)
  expected <- matrix(c(3,3,15,15), nrow=n_axes, ncol=n_families)
  centroids <- group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_all')
  assert_true(all(abs(centroids - expected) < 1e-12))
}

# 2. Basic functionality in 'orthologs' mode
test_basic_ortho_mode <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  ortholog_set <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
  expected <- matrix(c(3,3,15,15), nrow=n_axes, ncol=n_families)
  centroids <- group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_orthologs')
  assert_true(all(abs(centroids - expected) < 1e-12))
}

# 3. A family exists but has no genes assigned to it
test_empty_family <- function() {
  n_axes <- 3; n_genes <- 4; n_families <- 2
  vectors <- matrix(1.0, nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(rep(1, n_genes))
  ortholog_set <- rep(TRUE, n_genes)
  expected <- matrix(0.0, nrow=n_axes, ncol=n_families)
  expected[,1] <- 1.0
  centroids <- group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_all')
  assert_true(all(abs(centroids - expected) < 1e-12))
}

# 4. 'orthologs' mode, but a family has no orthologs
test_no_matching_orthologs <- function() {
  n_axes <- 2; n_genes <- 3; n_families <- 1
  vectors <- matrix(c(10,10,20,20,30,30), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(rep(1, n_genes))
  ortholog_set <- rep(FALSE, n_genes)
  expected <- matrix(0.0, nrow=n_axes, ncol=n_families)
  centroids <- group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_orthologs')
  assert_true(all(abs(centroids - expected) < 1e-12))
}

# 5. A family contains only a single gene
test_single_gene_family <- function() {
  n_axes <- 3; n_genes <- 1; n_families <- 1
  vectors <- matrix(c(12.3, -4.5, 6.7), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(1)
  ortholog_set <- TRUE
  centroids <- group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_all')
  assert_true(all(abs(centroids - vectors) < 1e-12))
}

# 6. Input vectors with extreme values
test_extreme_values <- function() {
  n_axes <- 2; n_genes <- 4; n_families <- 1
  vectors <- matrix(0.0, nrow=n_axes, ncol=n_genes)
  vectors[,1] <- c(1e12, -1e-12)
  vectors[,2] <- c(-1e12, 1e-12)
  vectors[,3] <- c(0, 5)
  vectors[,4] <- c(0, -5)
  gene_to_family <- as.integer(rep(1, n_genes))
  ortholog_set <- rep(TRUE, n_genes)
  expected <- matrix(0.0, nrow=n_axes, ncol=n_families)
  centroids <- group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_all')
  assert_true(all(abs(centroids - expected) < 1e-12))
}

# 7. Higher dimensional data
test_higher_dimensions <- function() {
  n_axes <- 10; n_genes <- 100; n_families <- 5
  vectors <- matrix(0.0, nrow=n_axes, ncol=n_genes)
  gene_to_family <- integer(n_genes)
  for (i in 1:n_genes) {
    vectors[,i] <- i
    gene_to_family[i] <- as.integer(((i-1) %% n_families) + 1)
  }
  ortholog_set <- rep(TRUE, n_genes)
  centroids <- group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_all')
  # Only check family 1 centroid
  idxs <- which(gene_to_family == 1)
  expected <- apply(vectors[,idxs,drop=FALSE], 1, mean)
  assert_true(all(abs(centroids[,1] - expected) < 1e-12))
}

# 8. Ensure result is invariant to order of genes
test_gene_order_invariance <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors1 <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family1 <- as.integer(c(1,1,2,2,1))
  ortholog_set1 <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
  vectors2 <- matrix(c(5,5,10,10,1,1,3,3,20,20), nrow=n_axes, ncol=n_genes)
  gene_to_family2 <- as.integer(c(1,2,1,1,2))
  ortholog_set2 <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  centroids1 <- group_centroid(vectors1, gene_to_family1, n_families, ortholog_set1, mode='group_orthologs')
  centroids2 <- group_centroid(vectors2, gene_to_family2, n_families, ortholog_set2, mode='group_orthologs')
  assert_true(all(abs(centroids1 - centroids2) < 1e-12))
}

# 9. Test for invalid input arguments
test_invalid_input_arguments <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  ortholog_set <- rep(FALSE, n_genes)
  # Invalid n_axes
  assert_error(group_centroid(matrix(numeric(0), nrow=0, ncol=n_genes), gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_all'), "Expected error for invalid n_axes")
  # Invalid n_genes
  assert_error(group_centroid(matrix(numeric(0), nrow=n_axes, ncol=0), integer(0), n_families, logical(0), mode='group_all'), "Expected error for invalid n_genes")
  # Invalid n_families
  assert_error(group_centroid(vectors, gene_to_family, 0, ortholog_set=ortholog_set, mode='group_all'), "Expected error for invalid n_families")
}

# 10. Test for invalid family mapping
test_invalid_family_mapping <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,3,1)) # 3 is invalid
  ortholog_set <- rep(FALSE, n_genes)
  assert_error(group_centroid(vectors, gene_to_family, n_families, ortholog_set=ortholog_set, mode='group_all'), "Expected error for invalid families")
}

run_all_tests()
