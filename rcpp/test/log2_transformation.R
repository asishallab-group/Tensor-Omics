# =====================================================================
#   TESTS FOR log2_transformation_rcpp
# =====================================================================

library(Rcpp)
source("rcpp/tensoromics_functions.R")

cat("=== TESTING: log2_transformation_rcpp ===\n")

# Helper function for expected output
expected_log2 <- function(x) {
  log(x + 1, base = 2)
}

# --------------------------------------------------------------
# Test 1: Basic vector
# --------------------------------------------------------------
test_log2_basic <- function() {
  cat("\n[test_log2_basic] Basic test\n")

  input <- matrix(c(0, 1, 3, 7), ncol = 1)
  result <- log2_transformation_rcpp(input)$output
  expected <- expected_log2(input)

  print(result)
  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Basic log2 test passed ✓\n")
}

# --------------------------------------------------------------
# Test 2: Matrix test
# --------------------------------------------------------------
test_log2_matrix <- function() {
  cat("\n[test_log2_matrix] Matrix test\n")

  input <- matrix(c(0,1,2,3,  4,5,6,7), nrow = 4, ncol = 2)
  result <- log2_transformation_rcpp(input)$output
  expected <- expected_log2(input)

  print(result)
  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Matrix test passed ✓\n")
}

# --------------------------------------------------------------
# Test 3: Large data performance
# --------------------------------------------------------------
test_log2_performance <- function() {
  cat("\n[test_log2_performance] Performance test\n")

  set.seed(123)
  input <- matrix(runif(500000), ncol = 50)

  start <- Sys.time()
  result <- log2_transformation_rcpp(input)$output
  elapsed <- as.numeric(Sys.time() - start, "secs")

  cat("Elapsed:", elapsed, "seconds\n")
  stopifnot(elapsed < 1.0)

  cat("Performance test passed ✓\n")
}

# --------------------------------------------------------------
# Test 4: Negative values should produce NaN (Fortran does not error)
# --------------------------------------------------------------
test_log2_negative_values <- function() {
  cat("\n[test_log2_negative_values] Negative values test\n")

  input <- matrix(c(1, 2, -5, 3), ncol = 1)
  result <- log2_transformation_rcpp(input)$output

  print(result)

  # Expected: log2(x+1) where x = -5 → log2(-4) = NaN
  stopifnot(is.nan(result[3, 1]))

  cat("Negative-value NaN behavior test passed ✓\n")
}


# --------------------------------------------------------------
# RUN ALL TESTS
# --------------------------------------------------------------

cat("\n=================================================\n")
cat(" RUNNING LOG2 TRANSFORMATION TESTS\n")
cat("=================================================\n")

test_log2_basic()
test_log2_matrix()
test_log2_performance()
test_log2_negative_values()

cat("\n=================================================\n")
cat("   ALL TESTS PASSED ✓\n")
cat("=================================================\n")
