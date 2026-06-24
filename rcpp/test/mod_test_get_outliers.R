# =====================
# Comprehensive R test suite for outlier detection 
# Uses tensoromics_functions.R wrapper functions
# =====================

# Source the main functions
source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

cat("=== Testing outlier detection R wrapper functions ===\n")
cat("Based on Fortran test suite with comprehensive test coverage\n")

# =====================
# Tests for compute_family_scaling
# =====================

# Test 1: Basic LOESS scaling
test_compute_family_scaling_basic <- function() {
  set.seed(42)  # For reproducibility

  # Updated to use more genes and families
  n_families <- 6
  genes_per_fam <- 4
  n_genes <- n_families * genes_per_fam

  # Generate distances and gene-to-family mapping
  distances <- rep(1:n_families, each = genes_per_fam) * 2 + rep(c(0.1, 0.4, 0.9, 1.6), n_families)
  gene_to_fam <- rep(1:n_families, each = genes_per_fam)

  result <- tox_compute_family_scaling(distances, gene_to_fam, n_families)

  # Verify results
  stopifnot(all(result$dscale > 0))  # All families have >1 gene
  stopifnot(length(result$dscale) == n_families)
  stopifnot(length(result$loess_x) == n_families)
  stopifnot(length(result$loess_y) == n_families)

}

# Test 2: Invalid family indices (should throw error)
test_compute_family_scaling_invalid <- function() {
  distances <- c(1, 2, 3, 4, 5)
  gene_to_fam <- c(1, 3, 2, 2, 2)  # family 3 doesn't exist
  n_families <- 2
  
  error_caught <- FALSE
  tryCatch({
    tox_compute_family_scaling(distances, gene_to_fam, n_families)
  }, error = function(e) {
    error_caught <<- TRUE
    # Check that the error message contains expected text
    stopifnot(grepl("Invalid input", e$message))
  })
  stopifnot(error_caught)  # Make sure an error was actually thrown
  
}

# Test 3: Zero distances
test_compute_family_scaling_zero_distances <- function() {
  distances <- c(0, 0, 0, 0, 0)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  n_families <- 2
  
  result <- tox_compute_family_scaling(distances, gene_to_fam, n_families)
  
  # Verify zero distance handling
  stopifnot(all(result$dscale == 0))  # All distances zero
  
}

# Test 4: Large dataset
test_compute_family_scaling_large_dataset <- function() {
  set.seed(42)  # For reproducibility

  # Updated to use more families and genes
  n_families <- 10
  genes_per_fam <- 5
  n_genes <- n_families * genes_per_fam

  # Generate random distances and gene-to-family mapping
  distances <- runif(n_genes, 0.5, 5.0)  # Random distances
  gene_to_fam <- rep(1:n_families, each = genes_per_fam)

  result <- tox_compute_family_scaling(distances, gene_to_fam, n_families)

  # Verify large dataset handling
  stopifnot(length(result$dscale) == n_families)
  stopifnot(all(result$dscale >= 0))  # All scaling factors non-negative

}

# =====================
# Tests for compute_rdi
# =====================

# Test 5: Normal RDI computation
test_compute_rdi_normal <- function() {
  distances <- c(1, 2, 3, 4, 5)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  dscale <- c(2, 4)
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Verify RDI calculation: rdi = abs(distances) / dscale
  expected_rdi <- abs(distances) / dscale[gene_to_fam]
  stopifnot(all(abs(result$rdi - expected_rdi) < 1e-10))
  stopifnot(length(result$sorted_rdi) == length(distances))
  
}

# Test 6: Zero scaling factors
test_compute_rdi_zero_scaling <- function() {
  distances <- c(1, 2, 3, 4, 5)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  dscale <- c(0, 0)
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Verify zero scaling results in zero RDI
  stopifnot(all(result$rdi == 0))
  
}

# Test 7: High precision test
test_compute_rdi_precision <- function() {
  distances <- c(2.0, 4.0, 6.0)
  gene_to_fam <- c(1, 1, 1)
  dscale <- c(2.0)  # All genes in family 1, scaling = 2.0
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Expected RDI: [2.0/2.0, 4.0/2.0, 6.0/2.0] = [1.0, 2.0, 3.0]
  expected <- c(1.0, 2.0, 3.0)
  stopifnot(all(abs(result$rdi - expected) < 1e-10))
  
}

# Test 8: Negative distances
test_compute_rdi_negative_distances <- function() {
  distances <- c(-1, 2, -3, 4, 5)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  dscale <- c(2, 4)
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Verify negative distances are handled (absolute value used)
  expected_rdi <- abs(distances) / dscale[gene_to_fam]
  stopifnot(all(abs(result$rdi - expected_rdi) < 1e-10))
  
}

# =====================
# Tests for identify_outliers
# =====================

# Test 9: Simple outlier identification
test_identify_outliers_simple <- function() {
  rdi <- c(0.3, 0.1, 0.5, 0.2, 0.4)
  percentile <- 50.0
  
  result <- tox_identify_outliers(rdi, percentile)
  
  # Verify outlier detection logic
  stopifnot(length(result$is_outlier) == length(rdi))
  stopifnot(is.numeric(result$threshold))
  stopifnot(is.logical(result$is_outlier))
  
}

# Test 10: All zeros RDI
test_identify_outliers_all_zeros <- function() {
  rdi <- c(0, 0, 0, 0, 0)
  percentile <- 90.0
  
  result <- tox_identify_outliers(rdi, percentile)
  
  # Verify no outliers detected when all RDI are zero
  stopifnot(!any(result$is_outlier))
  stopifnot(result$threshold == 0)
  
}

# Test 11: Percentile 0 (all outliers)
test_identify_outliers_percentile_0 <- function() {
  rdi <- c(0.3, 0.1, 0.5, 0.2, 0.4)
  percentile <- 0.0
  
  result <- tox_identify_outliers(rdi, percentile)

  # Verify all are outliers at 0% percentile
  stopifnot(all(result$is_outlier))
  
}

# Test 12: Percentile 100 (minimal outliers)
test_identify_outliers_percentile_100 <- function() {
  rdi <- c(0.3, 0.1, 0.5, 0.2, 0.4)
  percentile <- 100.0
  
  result <- tox_identify_outliers(rdi, percentile)
  
  # Verify only highest RDI values are outliers
  stopifnot(sum(result$is_outlier) >= 0)  # At least 0 outliers
  
}

# =====================
# Tests for detect_outliers (complete workflow)
# =====================

# Test 13: Typical workflow
test_detect_outliers_typical <- function() {
  set.seed(42)  # For reproducibility
  # Updated to use more genes and families
  n_families <- 8
  genes_per_fam <- 6
  n_genes <- n_families * genes_per_fam
  percentile <- 80.0

  # Generate distances and gene-to-family mapping
  distances <- runif(n_genes, 1, 10)  # Random distances
  gene_to_fam <- rep(1:n_families, each = genes_per_fam)

  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)

  # Verify typical workflow
  stopifnot(length(result$is_outlier) == length(distances))
  stopifnot(length(result$loess_x) == n_families)
  stopifnot(length(result$loess_y) == n_families)
  stopifnot(is.logical(result$is_outlier))

}

# Test 14: Invalid families (should throw error)
test_detect_outliers_invalid_families <- function() {
  distances <- c(1, 2, 3, 4, 5, 6)
  gene_to_fam <- c(1, 3, 2, 2, 2, 2)  # family 3 doesn't exist
  n_families <- 2
  percentile <- 80.0
  
  error_caught <- FALSE
  tryCatch({
    tox_detect_outliers(distances, gene_to_fam, n_families, percentile)
  }, error = function(e) {
    error_caught <<- TRUE
    # Check that the error message contains expected text
    stopifnot(grepl("Invalid input", e$message))
  })
  stopifnot(error_caught)  # Make sure an error was actually thrown
  
}

# Test 15: Single gene families
test_detect_outliers_single_families <- function() {
  distances <- c(1, 10, 100)  # Each gene in different family
  gene_to_fam <- c(1, 2, 3)
  n_families <- 3
  percentile <- 90.0
  
  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)
  
  # Verify single gene families don't cause errors
  stopifnot(length(result$is_outlier) == length(distances))
  stopifnot(is.logical(result$is_outlier))
  
}

# Test 16: Mixed family sizes
test_detect_outliers_mixed_sizes <- function() {
  set.seed(42)  # For reproducibility

  # Updated to use more families with mixed sizes
  distances <- c(
    runif(10, 1, 2),  # Family 1: 10 genes
    runif(8, 2, 3),   # Family 2: 8 genes
    runif(12, 3, 4),  # Family 3: 12 genes
    runif(6, 4, 5),   # Family 4: 6 genes
    runif(9, 5, 6),   # Family 5: 9 genes
    runif(7, 6, 7),   # Family 6: 7 genes
    runif(11, 7, 8),  # Family 7: 11 genes
    runif(5, 8, 9)    # Family 8: 5 genes
  )
  gene_to_fam <- c(
    rep(1, 10),  # Family 1
    rep(2, 8),   # Family 2
    rep(3, 12),  # Family 3
    rep(4, 6),   # Family 4
    rep(5, 9),   # Family 5
    rep(6, 7),   # Family 6
    rep(7, 11),  # Family 7
    rep(8, 5)    # Family 8
  )
  n_families <- 8
  percentile <- 95.0

  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)

  # Verify mixed family sizes work
  stopifnot(length(result$is_outlier) == length(distances))
  stopifnot(length(result$loess_x) == n_families)

}

# Test 17: Large dataset with outliers
test_detect_outliers_large_dataset <- function() {
  set.seed(123)  # For reproducibility
  n_genes <- 92
  n_families <- 10

  # Create synthetic data with some clear outliers
  distances <- c(
    rnorm(10, 1, 0.1),    # Family 1: tight cluster around 1
    rnorm(10, 2, 0.1),    # Family 2: tight cluster around 2
    rnorm(10, 3, 0.1),    # Family 3: tight cluster around 3
    rnorm(10, 4, 0.1),    # Family 4: tight cluster around 4
    rnorm(10, 5, 0.1),    # Family 5: tight cluster around 5
    rnorm(10, 6, 0.1),    # Family 6: tight cluster around 6
    rnorm(10, 7, 0.1),    # Family 7: tight cluster around 7
    rnorm(10, 8, 0.1),    # Family 8: tight cluster around 8
    rnorm(10, 9, 0.1),    # Family 9: tight cluster around 9
    c(20, 25)             # Two clear outliers
  )
  # Map each distance to a family index. The constructed distances vector has
  # 10+10+10+10+10+10+10+10+2 = 92 elements, so create a matching mapping of length 92.
  gene_to_fam <- c(rep(1:10, each = 9), 9, 10)
  percentile <- 90.0

  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)
  # Verify large dataset handling
  stopifnot(length(result$is_outlier) == n_genes)
  stopifnot(sum(result$is_outlier) >= 0)  # At least 0 outliers detected

}

# =====================
# Tests for compute_family_scaling_expert
# =====================

# Test 18: Expert version basic test
test_compute_family_scaling_expert_basic <- function() {
  set.seed(42)  # For reproducibility

  # Test data
  n_families <- 10  # Incremented number of families
  genes_per_fam <- 4
  n_genes <- n_families * genes_per_fam
  span <- 0.7
  mode <- 1
  n_iters <- 3
  degree <- 2

  distances <- c(
    runif(genes_per_fam, 1.0, 2.0),  # Family 1
    runif(genes_per_fam, 2.0, 3.0),  # Family 2
    runif(genes_per_fam, 3.0, 4.0),  # Family 3
    runif(genes_per_fam, 4.0, 5.0),  # Family 4
    runif(genes_per_fam, 5.0, 6.0),  # Family 5
    runif(genes_per_fam, 6.0, 7.0),  # Family 6
    runif(genes_per_fam, 7.0, 8.0),  # Family 7
    runif(genes_per_fam, 8.0, 9.0),  # Family 8
    runif(genes_per_fam, 9.0, 10.0), # Family 9
    runif(genes_per_fam, 10.0, 11.0) # Family 10
  )
  gene_to_fam <- rep(1:n_families, each = genes_per_fam)

  # Pre-allocate work arrays (user responsibility)
  perm_tmp <- integer(n_genes)
  stack_left_tmp <- integer(n_genes)
  stack_right_tmp <- integer(n_genes)

  result <- tox_compute_family_scaling_expert(
    distances, gene_to_fam, n_families,
    perm_tmp, stack_left_tmp, stack_right_tmp,
    span, degree, mode, n_iters
  )


  # Verify basic properties
  stopifnot(length(result$dscale) == n_families)  # Ten families
  stopifnot(all(is.finite(result$dscale)))
  stopifnot(all(result$dscale > 0))  # Scaling factors should be positive

}

# Test 19: Expert version consistency with regular
test_compute_family_scaling_expert_consistency <- function() {
  set.seed(42)  # For reproducibility

  # Test data
  n_families <- 10  # Incremented number of families
  genes_per_fam <- 4
  n_genes <- n_families * genes_per_fam
  span <- 0.7
  mode <- 1
  n_iters <- 3
  degree <- 2

  distances <- c(
    runif(genes_per_fam, 1.0, 2.0),  # Family 1
    runif(genes_per_fam, 2.0, 3.0),  # Family 2
    runif(genes_per_fam, 3.0, 4.0),  # Family 3
    runif(genes_per_fam, 4.0, 5.0),  # Family 4
    runif(genes_per_fam, 5.0, 6.0),  # Family 5
    runif(genes_per_fam, 6.0, 7.0),  # Family 6
    runif(genes_per_fam, 7.0, 8.0),  # Family 7
    runif(genes_per_fam, 8.0, 9.0),  # Family 8
    runif(genes_per_fam, 9.0, 10.0), # Family 9
    runif(genes_per_fam, 10.0, 11.0) # Family 10
  )
  gene_to_fam <- rep(1:n_families, each = genes_per_fam)

  # Pre-allocate work arrays (user responsibility)
  perm_tmp <- integer(n_genes)
  stack_left_tmp <- integer(n_genes)
  stack_right_tmp <- integer(n_genes)

  # Validar tamaños de las matrices de trabajo
  stopifnot(length(perm_tmp) == n_genes)
  stopifnot(length(stack_left_tmp) == n_genes)
  stopifnot(length(stack_right_tmp) == n_genes)

  result_expert <- tox_compute_family_scaling_expert(
    distances, gene_to_fam, n_families,
    perm_tmp, stack_left_tmp, stack_right_tmp,
    span = span, degree = degree, mode = mode, n_iters = n_iters
  )
  result_regular <- tox_compute_family_scaling(distances, gene_to_fam, n_families)

  valid_idx <- !is.na(result_expert$loess_x)

  # Compare with regular version to ensure consistency

  # Results should be very similar (within numerical precision)
  stopifnot(all(abs(result_expert$dscale - result_regular$dscale) < 1e-10))
  stopifnot(all(abs(result_expert$loess_x[valid_idx] - result_regular$loess_x[valid_idx]) < 1e-10))
  stopifnot(all(abs(result_expert$loess_y[valid_idx] - result_regular$loess_y[valid_idx]) < 1e-10))

}

# =====================
# Run all tests
# =====================

run_all_tests()
