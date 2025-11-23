# Comprehensive R test suite for gene centroids interface functions

source("r/tensoromics_functions.R")

# 1. Basic functionality in 'all' mode
test_basic_all_mode <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  mode <- 'all'
  expected <- matrix(c(3,3,15,15), nrow=n_axes, ncol=n_families)
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode)
  stopifnot(all(abs(centroids - expected) < 1e-12))
  cat("test_basic_all_mode passed\n")
}

# 2. Basic functionality in 'ortho' mode
test_basic_ortho_mode <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  mode <- 'orthologs'
  ortholog_set <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
  expected <- matrix(c(3,3,15,15), nrow=n_axes, ncol=n_families)
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode, ortholog_set)
  stopifnot(all(abs(centroids - expected) < 1e-12))
  cat("test_basic_ortho_mode passed\n")
}

# 3. A family exists but has no genes assigned to it
test_empty_family <- function() {
  n_axes <- 3; n_genes <- 4; n_families <- 2
  vectors <- matrix(1.0, nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(rep(1, n_genes))
  mode <- 'all'
  expected <- matrix(0.0, nrow=n_axes, ncol=n_families)
  expected[,1] <- 1.0
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode)
  stopifnot(all(abs(centroids - expected) < 1e-12))
  cat("test_empty_family passed\n")
}

# 4. 'orthologs' mode, but a family has no orthologs
test_no_matching_orthologs <- function() {
  n_axes <- 2; n_genes <- 3; n_families <- 1
  vectors <- matrix(c(10,10,20,20,30,30), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(rep(1, n_genes))
  mode <- 'orthologs'
  ortholog_set <- rep(FALSE, n_genes)
  expected <- matrix(0.0, nrow=n_axes, ncol=n_families)
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode, ortholog_set)
  stopifnot(all(abs(centroids - expected) < 1e-12))
  cat("test_no_matching_orthologs passed\n")
}

# 5. A family contains only a single gene
test_single_gene_family <- function() {
  n_axes <- 3; n_genes <- 1; n_families <- 1
  vectors <- matrix(c(12.3, -4.5, 6.7), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(1)
  mode <- 'all'
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode)
  stopifnot(all(abs(centroids - vectors) < 1e-12))
  cat("test_single_gene_family passed\n")
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
  mode <- 'all'
  expected <- matrix(0.0, nrow=n_axes, ncol=n_families)
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode)
  stopifnot(all(abs(centroids - expected) < 1e-12))
  cat("test_extreme_values passed\n")
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
  mode <- 'all'
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode)
  # Only check family 1 centroid
  idxs <- which(gene_to_family == 1)
  expected <- apply(vectors[,idxs,drop=FALSE], 1, mean)
  stopifnot(all(abs(centroids[,1] - expected) < 1e-12))
  cat("test_higher_dimensions passed\n")
}

# 8. Ensure result is invariant to order of genes
test_gene_order_invariance <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors1 <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family1 <- as.integer(c(1,1,2,2,1))
  mode <- 'orthologs'
  ortholog_set1 <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
  vectors2 <- matrix(c(5,5,10,10,1,1,3,3,20,20), nrow=n_axes, ncol=n_genes)
  gene_to_family2 <- as.integer(c(1,2,1,1,2))
  ortholog_set2 <- c(TRUE, TRUE, TRUE, FALSE, TRUE)
  centroids1 <- tox_group_centroid(vectors1, gene_to_family1, n_families, mode, ortholog_set1)
  centroids2 <- tox_group_centroid(vectors2, gene_to_family2, n_families, mode, ortholog_set2)
  stopifnot(all(abs(centroids1 - centroids2) < 1e-12))
  cat("test_gene_order_invariance passed\n")
}

# 9. Test for invalid input arguments
test_invalid_input_arguments <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  mode <- 'all'
  # Invalid n_axes
  error_caught <- FALSE
  tryCatch({
    tox_group_centroid(matrix(numeric(0), nrow=0, ncol=n_genes), gene_to_family, n_families, mode)
  }, error=function(e) { error_caught <<- TRUE })
  stopifnot(error_caught)
  # Invalid n_genes
  error_caught <- FALSE
  tryCatch({
    tox_group_centroid(matrix(numeric(0), nrow=n_axes, ncol=0), integer(0), n_families, mode)
  }, error=function(e) { error_caught <<- TRUE })
  stopifnot(error_caught)
  # Invalid n_families
  error_caught <- FALSE
  tryCatch({
    tox_group_centroid(vectors, gene_to_family, 0, mode)
  }, error=function(e) { error_caught <<- TRUE })
  stopifnot(error_caught)
  cat("test_invalid_input_arguments passed\n")
}

# 10. Test for invalid family mapping
test_invalid_family_mapping <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,3,1)) # 3 is invalid
  mode <- 'all'
  error_caught <- FALSE
  tryCatch({
    tox_group_centroid(vectors, gene_to_family, n_families, mode)
  }, error=function(e) { error_caught <<- TRUE })
  stopifnot(error_caught)
  cat("test_invalid_family_mapping passed\n")
}

# 11. Test for invalid mode string
test_invalid_mode_string <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  mode <- ''
  error_caught <- FALSE
  tryCatch({
    tox_group_centroid(vectors, gene_to_family, n_families, mode)
  }, error=function(e) { error_caught <<- TRUE })
  stopifnot(error_caught)

  error_caught <- FALSE
  mode <- 'invalid_mode'
  tryCatch({
    tox_group_centroid(vectors, gene_to_family, n_families, mode)
  }, error=function(e) { error_caught <<- TRUE })
  stopifnot(error_caught)  
  cat("test_invalid_mode_string passed\n")
}

# 12. Test for missing ortholog set in 'orthologs' mode
test_missing_ortholog_set <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  mode <- 'orthologs'
  error_caught <- FALSE
  tryCatch({
    tox_group_centroid(vectors, gene_to_family, n_families, mode)
  }, error=function(e) { error_caught <<- TRUE })
  stopifnot(error_caught)
  cat("test_missing_ortholog_set passed\n")
}

#13. Test for providing ortholog set in 'all' mode
test_present_ortholog_set_in_all_mode <- function() {
  n_axes <- 2; n_genes <- 5; n_families <- 2
  vectors <- matrix(c(1,1,3,3,10,10,20,20,5,5), nrow=n_axes, ncol=n_genes)
  gene_to_family <- as.integer(c(1,1,2,2,1))
  ortholog_set <- c(TRUE, FALSE, TRUE, TRUE, TRUE)
  mode <- 'all'
  expected <- matrix(c(3,3,15,15), nrow=n_axes, ncol=n_families)
  centroids <- tox_group_centroid(vectors, gene_to_family, n_families, mode, ortholog_set)
  stopifnot(all(abs(centroids - expected) < 1e-12))
  cat("test_present_ortholog_set_in_all_mode passed\n")
}

# Run all tests
cat("=================================================\n")
cat("    GENE CENTROIDS FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

test_basic_all_mode()
test_basic_ortho_mode()
test_empty_family()
test_no_matching_orthologs()
test_single_gene_family()
test_extreme_values()
test_higher_dimensions()
test_gene_order_invariance()
test_invalid_input_arguments()
test_invalid_family_mapping()
test_invalid_mode_string()
test_missing_ortholog_set()
test_present_ortholog_set_in_all_mode()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all gene centroids R interface tests passed successfully!\n")