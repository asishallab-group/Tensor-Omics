# =====================================================================
#   TESTS FOR calc_tiss_avg_rcpp
# =====================================================================

library(Rcpp)
source("rcpp/tensoromics_functions.R")

cat("=== TESTING: calc_tiss_avg_rcpp ===\n")

# --------------------------------------------------------------
# Helper: R implementation of group averages
# --------------------------------------------------------------
r_avg <- function(mat, group_s, group_c) {
  n_genes <- nrow(mat)
  n_groups <- length(group_s)
  out <- matrix(0, n_genes, n_groups)

  for (g in 1:n_groups) {
    start <- group_s[g]
    count <- group_c[g]
    cols <- start:(start + count - 1)
    out[, g] <- rowMeans(mat[, cols, drop=FALSE])
  }
  out
}

# --------------------------------------------------------------
# Test 1: Basic averaging
# --------------------------------------------------------------
test_avg_basic <- function() {
  cat("\n[test_avg_basic] Basic averaging\n")

  input <- matrix(
    c(1, 2, 3,
      4, 5, 6,
      7, 8, 9),
    nrow = 3, byrow = TRUE
  )

  group_s <- c(1, 3)   # group1 = cols 1–2, group2 = col3 only
  group_c <- c(2, 1)

  result <- calc_tiss_avg_rcpp(input, group_s, group_c)$output
  expected <- r_avg(input, group_s, group_c)

  print(result)
  print(expected)

  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Basic group averaging test passed ✓\n")
}

# --------------------------------------------------------------
# Test 2: Single-column groups
# --------------------------------------------------------------
test_avg_single_column <- function() {
  cat("\n[test_avg_single_column]\n")

  input <- matrix(
    c(10, 11,
      20, 21,
      30, 31),
    nrow = 3, byrow = TRUE
  )

  group_s <- c(1,2)
  group_c <- c(1,1)

  result <- calc_tiss_avg_rcpp(input, group_s, group_c)$output
  expected <- input

  print(result)

  stopifnot(all(result == expected))

  cat("Single-column groups test passed ✓\n")
}

# --------------------------------------------------------------
# Test 3: Multiple groups of different sizes
# --------------------------------------------------------------
test_avg_multiple_groups <- function() {
  cat("\n[test_avg_multiple_groups]\n")

  input <- matrix(
    c(1,2,3,4,
      5,6,7,8,
      9,10,11,12),
    nrow = 3, byrow = TRUE
  )

  group_s <- c(1, 3, 4)
  group_c <- c(2, 1, 1)

  result <- calc_tiss_avg_rcpp(input, group_s, group_c)$output
  expected <- r_avg(input, group_s, group_c)

  print(result)
  print(expected)

  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Multiple groups test passed ✓\n")
}

# --------------------------------------------------------------
# Test 4: Rows with all zeros
# --------------------------------------------------------------
test_avg_zero_rows <- function() {
  cat("\n[test_avg_zero_rows]\n")

  input <- matrix(
    c(0,0,0,
      5,5,5,
      0,10,20),
    nrow = 3, byrow = TRUE
  )

  group_s <- c(1, 3)
  group_c <- c(2, 1)

  result <- calc_tiss_avg_rcpp(input, group_s, group_c)$output
  expected <- r_avg(input, group_s, group_c)

  print(result)

  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Zero-rows test passed ✓\n")
}

# --------------------------------------------------------------
# Test 5: Performance test
# --------------------------------------------------------------
test_avg_performance <- function() {
  cat("\n[test_avg_performance]\n")

  set.seed(123)
  input <- matrix(rnorm(200000), nrow=10000, ncol=20)

  # 4 groups
  group_s <- c(1, 6, 11, 16)
  group_c <- c(5, 5, 5, 5)

  start <- Sys.time()
  result <- calc_tiss_avg_rcpp(input, group_s, group_c)$output
  elapsed <- as.numeric(Sys.time() - start, "secs")

  cat("Elapsed:", elapsed, "seconds\n")

  stopifnot(elapsed < 1.5)

  cat("Performance test passed ✓\n")
}

# --------------------------------------------------------------
# No invalid-input tests (Fortran does NOT validate inputs)
# --------------------------------------------------------------

# --------------------------------------------------------------
# RUN ALL TESTS
# --------------------------------------------------------------
cat("\n=================================================\n")
cat(" RUNNING calc_tiss_avg_rcpp TESTS\n")
cat("=================================================\n")

test_avg_basic()
test_avg_single_column()
test_avg_multiple_groups()
test_avg_zero_rows()
test_avg_performance()

cat("\n=================================================\n")
cat("   ALL TESTS PASSED ✓\n")
cat("=================================================\n")
