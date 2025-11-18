# =====================================================================
#   TESTS FOR tox_quantile_normalization
# =====================================================================

library(Rcpp)
source("rcpp/tensoromics_functions.R")

cat("=== TESTING: tox_quantile_normalization ===\n")

# --------------------------------------------------------------
# Helper: pure R quantile normalization
# --------------------------------------------------------------
quantile_normalize_R <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)

  # 1. sort each column
  sorted <- apply(mat, 2, sort)

  # 2. compute rank means
  rank_means <- rowMeans(sorted)

  # 3. reconstruct by ordering each column back
  res <- mat
  for (j in 1:nc) {
    ord <- order(mat[, j])
    res[ord, j] <- rank_means
  }

  res
}

# --------------------------------------------------------------
# Test 1: Basic correctness
# --------------------------------------------------------------
test_qn_basic <- function() {
  cat("\n[test_qn_basic] Basic quantile normalization\n")

  input <- matrix(c(5,4,3,1,
                    2,10,6,7,
                    9,8,11,12), nrow=3, byrow=TRUE)

  result <- tox_quantile_normalization(input)
  expected <- quantile_normalize_R(input)

  print(result)
  print(expected)

  stopifnot(all(abs(result - expected) < 1e-12))
  cat("Basic QN test passed ✓\n")
}

# --------------------------------------------------------------
# Test 2: Identical columns → output must remain identical
# --------------------------------------------------------------
test_qn_identical <- function() {
  cat("\n[test_qn_identical] Identical columns\n")

  input <- matrix(rep(c(1,5,10), 4), nrow=3, ncol=4)

  result <- tox_quantile_normalization(input)
  expected <- quantile_normalize_R(input)

  print(result)
  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Identical columns test passed ✓\n")
}

# --------------------------------------------------------------
# Test 3: Single column → sorting only, but values remain same
# --------------------------------------------------------------
test_qn_single_column <- function() {
  cat("\n[test_qn_single_column]\n")

  input <- matrix(c(5,1,9,3), ncol=1)
  result <- tox_quantile_normalization(input)

  # Fortran behavior: single column remains unchanged
  expected <- input

  print(result)
  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Single-column test passed ✓\n")
}


# --------------------------------------------------------------
# Test 4: Performance test
# --------------------------------------------------------------
test_qn_performance <- function() {
  cat("\n[test_qn_performance] Performance test\n")

  set.seed(123)
  input <- matrix(rnorm(300000), nrow = 10000, ncol = 30)

  start <- Sys.time()
  result <- tox_quantile_normalization(input)
  elapsed <- as.numeric(Sys.time() - start, "secs")

  cat("Elapsed:", elapsed, "seconds\n")
  stopifnot(elapsed < 2.5)

  cat("Performance test passed ✓\n")
}

# --------------------------------------------------------------
# RUN ALL TESTS
# --------------------------------------------------------------

cat("\n=================================================\n")
cat(" RUNNING QUANTILE NORMALIZATION TESTS\n")
cat("=================================================\n")

test_qn_basic()
test_qn_identical()
test_qn_single_column()
test_qn_performance()

cat("\n=================================================\n")
cat("   ALL TESTS PASSED ✓\n")
cat("=================================================\n")
