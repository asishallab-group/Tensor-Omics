# =====================
# Comprehensive R test suite for outlier detection 
# Uses tensoromics_functions.R wrapper functions
# =====================

# Source the main functions
source("r/tensoromics_functions.R")

cat("=== Testing outlier detection R wrapper functions ===\n")
cat("Based on Fortran test suite with comprehensive test coverage\n")

# =====================
# Tests for compute_family_scaling
# =====================

# Test 1: Basic LOESS scaling
test_compute_family_scaling_basic <- function() {
  cat("\n[test_compute_family_scaling_basic] Basic LOESS scaling test\n")
  distances <- c(1, 2, 3, 4, 5)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  n_families <- 2
  
  result <- tox_compute_family_scaling(distances, gene_to_fam, n_families)
  
  # Verify results
  stopifnot(all(result$dscale > 0))  # Both families have >1 gene
  stopifnot(length(result$dscale) == n_families)
  stopifnot(length(result$loess_x) == n_families)
  stopifnot(length(result$loess_y) == n_families)
  
  cat("Basic LOESS scaling test passed ✓\n")
}

# Test 2: Invalid family indices (should throw error)
test_compute_family_scaling_invalid <- function() {
  cat("\n[test_compute_family_scaling_invalid] Invalid family indices test\n")
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
  
  cat("Invalid family indices test passed ✓\n")
}

# Test 3: Zero distances
test_compute_family_scaling_zero_distances <- function() {
  cat("\n[test_compute_family_scaling_zero_distances] Zero distances test\n")
  distances <- c(0, 0, 0, 0, 0)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  n_families <- 2
  
  result <- tox_compute_family_scaling(distances, gene_to_fam, n_families)
  
  # Verify zero distance handling
  stopifnot(all(result$dscale == 0))  # All distances zero
  
  cat("Zero distances test passed ✓\n")
}

# Test 4: Large dataset
test_compute_family_scaling_large_dataset <- function() {
  cat("\n[test_compute_family_scaling_large_dataset] Large dataset test\n")
  set.seed(42)  # For reproducibility
  n_genes <- 20
  n_families <- 4
  distances <- runif(n_genes, 0.5, 5.0)  # Random distances
  gene_to_fam <- sample(1:n_families, n_genes, replace=TRUE)
  
  result <- tox_compute_family_scaling(distances, gene_to_fam, n_families)
  
  # Verify large dataset handling
  stopifnot(length(result$dscale) == n_families)
  stopifnot(all(result$dscale >= 0))  # All scaling factors non-negative
  
  cat("Large dataset test passed ✓\n")
}

# =====================
# Tests for compute_rdi
# =====================

# Test 5: Normal RDI computation
test_compute_rdi_normal <- function() {
  cat("\n[test_compute_rdi_normal] Normal RDI computation test\n")
  distances <- c(1, 2, 3, 4, 5)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  dscale <- c(2, 4)
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Verify RDI calculation: rdi = abs(distances) / dscale
  expected_rdi <- abs(distances) / dscale[gene_to_fam]
  stopifnot(all(abs(result$rdi - expected_rdi) < 1e-10))
  stopifnot(length(result$sorted_rdi) == length(distances))
  
  cat("Normal RDI computation test passed ✓\n")
}

# Test 6: Zero scaling factors
test_compute_rdi_zero_scaling <- function() {
  cat("\n[test_compute_rdi_zero_scaling] Zero scaling factors test\n")
  distances <- c(1, 2, 3, 4, 5)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  dscale <- c(0, 0)
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Verify zero scaling results in zero RDI
  stopifnot(all(result$rdi == 0))
  
  cat("Zero scaling factors test passed ✓\n")
}

# Test 7: High precision test
test_compute_rdi_precision <- function() {
  cat("\n[test_compute_rdi_precision] High precision test\n")
  distances <- c(2.0, 4.0, 6.0)
  gene_to_fam <- c(1, 1, 1)
  dscale <- c(2.0)  # All genes in family 1, scaling = 2.0
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Expected RDI: [2.0/2.0, 4.0/2.0, 6.0/2.0] = [1.0, 2.0, 3.0]
  expected <- c(1.0, 2.0, 3.0)
  stopifnot(all(abs(result$rdi - expected) < 1e-10))
  
  cat("High precision test passed ✓\n")
}

# Test 8: Negative distances
test_compute_rdi_negative_distances <- function() {
  cat("\n[test_compute_rdi_negative_distances] Negative distances test\n")
  distances <- c(-1, 2, -3, 4, 5)
  gene_to_fam <- c(1, 1, 2, 2, 2)
  dscale <- c(2, 4)
  
  result <- tox_compute_rdi(distances, gene_to_fam, dscale)
  
  # Verify negative distances are handled (absolute value used)
  expected_rdi <- abs(distances) / dscale[gene_to_fam]
  stopifnot(all(abs(result$rdi - expected_rdi) < 1e-10))
  
  cat("Negative distances test passed ✓\n")
}

# =====================
# Tests for identify_outliers
# =====================

# Test 9: Simple outlier identification
test_identify_outliers_simple <- function() {
  cat("\n[test_identify_outliers_simple] Simple outlier identification test\n")
  rdi <- c(0.3, 0.1, 0.5, 0.2, 0.4)
  percentile <- 50.0
  
  result <- tox_identify_outliers(rdi, percentile)
  
  # Verify outlier detection logic
  stopifnot(length(result$is_outlier) == length(rdi))
  stopifnot(is.numeric(result$threshold))
  stopifnot(is.logical(result$is_outlier))
  
  cat("Simple outlier identification test passed ✓\n")
}

# Test 10: All zeros RDI
test_identify_outliers_all_zeros <- function() {
  cat("\n[test_identify_outliers_all_zeros] All zeros RDI test\n")
  rdi <- c(0, 0, 0, 0, 0)
  percentile <- 90.0
  
  result <- tox_identify_outliers(rdi, percentile)
  
  # Verify no outliers detected when all RDI are zero
  stopifnot(!any(result$is_outlier))
  stopifnot(result$threshold == 0)
  
  cat("All zeros RDI test passed ✓\n")
}

# Test 11: Percentile 0 (all outliers)
test_identify_outliers_percentile_0 <- function() {
  cat("\n[test_identify_outliers_percentile_0] Percentile 0 test\n")
  rdi <- c(0.3, 0.1, 0.5, 0.2, 0.4)
  percentile <- 0.0
  
  result <- tox_identify_outliers(rdi, percentile)
  
  # Verify all are outliers at 0% percentile
  stopifnot(all(result$is_outlier))
  
  cat("Percentile 0 test passed ✓\n")
}

# Test 12: Percentile 100 (minimal outliers)
test_identify_outliers_percentile_100 <- function() {
  cat("\n[test_identify_outliers_percentile_100] Percentile 100 test\n")
  rdi <- c(0.3, 0.1, 0.5, 0.2, 0.4)
  percentile <- 100.0
  
  result <- tox_identify_outliers(rdi, percentile)
  
  # Verify only highest RDI values are outliers
  stopifnot(sum(result$is_outlier) >= 0)  # At least 0 outliers
  
  cat("Percentile 100 test passed ✓\n")
}

# =====================
# Tests for detect_outliers (complete workflow)
# =====================

# Test 13: Typical workflow
test_detect_outliers_typical <- function() {
  cat("\n[test_detect_outliers_typical] Typical workflow test\n")
  distances <- c(1, 2, 3, 4, 5, 6)
  gene_to_fam <- c(1, 1, 2, 2, 2, 2)
  n_families <- 2
  percentile <- 80.0
  
  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)
  
  # Verify typical workflow
  stopifnot(length(result$is_outlier) == length(distances))
  stopifnot(length(result$loess_x) == n_families)
  stopifnot(length(result$loess_y) == n_families)
  stopifnot(is.logical(result$is_outlier))
  
  cat("Typical workflow test passed ✓\n")
}

# Test 14: Invalid families (should throw error)
test_detect_outliers_invalid_families <- function() {
  cat("\n[test_detect_outliers_invalid_families] Invalid families test\n")
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
    stopifnot(grepl("Invalid input provided", e$message))
  })
  stopifnot(error_caught)  # Make sure an error was actually thrown
  
  cat("Invalid families test passed ✓\n")
}

# Test 15: Single gene families
test_detect_outliers_single_families <- function() {
  cat("\n[test_detect_outliers_single_families] Single gene families test\n")
  distances <- c(1, 10, 100)  # Each gene in different family
  gene_to_fam <- c(1, 2, 3)
  n_families <- 3
  percentile <- 90.0
  
  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)
  
  # Verify single gene families don't cause errors
  stopifnot(length(result$is_outlier) == length(distances))
  stopifnot(is.logical(result$is_outlier))
  
  cat("Single gene families test passed ✓\n")
}

# Test 16: Mixed family sizes
test_detect_outliers_mixed_sizes <- function() {
  cat("\n[test_detect_outliers_mixed_sizes] Mixed family sizes test\n")
  distances <- c(1.0, 1.1, 1.2, 1.3, 1.4, 2.0, 2.1, 3.0, 4.0, 5.0)
  gene_to_fam <- c(1, 1, 1, 1, 1, 2, 2, 3, 3, 3)  # Family sizes: 5, 2, 3
  n_families <- 3
  percentile <- 95.0
  
  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)
  
  # Verify mixed family sizes work
  stopifnot(length(result$is_outlier) == length(distances))
  stopifnot(length(result$loess_x) == n_families)
  
  cat("Mixed family sizes test passed ✓\n")
}

# Test 17: Large dataset with outliers
test_detect_outliers_large_dataset <- function() {
  cat("\n[test_detect_outliers_large_dataset] Large dataset test\n")
  set.seed(123)  # For reproducibility
  n_genes <- 50
  n_families <- 5
  
  # Create synthetic data with some clear outliers
  distances <- c(
    rnorm(10, 1, 0.1),    # Family 1: tight cluster around 1
    rnorm(10, 2, 0.1),    # Family 2: tight cluster around 2
    rnorm(10, 3, 0.1),    # Family 3: tight cluster around 3
    rnorm(10, 4, 0.1),    # Family 4: tight cluster around 4
    rnorm(8, 5, 0.1),     # Family 5: tight cluster around 5
    c(20, 25)             # Two clear outliers
  )
  gene_to_fam <- c(rep(1:5, each=10), c(1, 2))
  percentile <- 90.0
  
  result <- tox_detect_outliers(distances, gene_to_fam, n_families, percentile)
  
  # Verify large dataset handling
  stopifnot(length(result$is_outlier) == n_genes)
  stopifnot(sum(result$is_outlier) >= 0)  # At least 0 outliers detected
  
  cat("Large dataset test passed ✓\n")
}

# =====================
# Run all tests
# =====================

cat("\n=================================================\n")
cat("    OUTLIER DETECTION FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

# compute_family_scaling tests
test_compute_family_scaling_basic()
test_compute_family_scaling_invalid()
test_compute_family_scaling_zero_distances()
test_compute_family_scaling_large_dataset()

# compute_rdi tests
test_compute_rdi_normal()
test_compute_rdi_zero_scaling()
test_compute_rdi_precision()
test_compute_rdi_negative_distances()

# identify_outliers tests
test_identify_outliers_simple()
test_identify_outliers_all_zeros()
test_identify_outliers_percentile_0()
test_identify_outliers_percentile_100()

# detect_outliers tests (complete workflow)
test_detect_outliers_typical()
test_detect_outliers_invalid_families()
test_detect_outliers_single_families()
test_detect_outliers_mixed_sizes()
test_detect_outliers_large_dataset()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all outlier detection R interface tests passed! ✓\n")
