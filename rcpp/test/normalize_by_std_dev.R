# =====================================================================
#   TESTS FOR tox_normalize_by_std_dev (RMS per-gene normalization)
# =====================================================================

library(Rcpp)
source("rcpp/tensoromics_functions.R")

cat("=== TESTING: tox_normalize_by_std_dev ===\n")

# --------------------------------------------------------------
# Helper: compute RMS per row (gene)
# --------------------------------------------------------------
rms_normalize <- function(mat) {
  apply(mat, 1, function(row) {
    rms <- sqrt(mean(row^2))
    if (rms == 0) return(rep(0, length(row)))
    row / rms
  }) |> t()
}

# --------------------------------------------------------------
# Test 1: Basic RMS normalization
# --------------------------------------------------------------
test_rms_basic <- function() {
  cat("\n[test_rms_basic] Basic RMS normalization\n")

  input <- matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)

  result <- tox_normalize_by_std_dev(input)   # <-- FIXED
  expected <- rms_normalize(input)

  print(result)
  print(expected)

  stopifnot(all(abs(result - expected) < 1e-12))
  cat("Basic RMS test passed ✓\n")
}

# --------------------------------------------------------------
# Test 2: Single-column → result = all 1
# --------------------------------------------------------------
test_single_column <- function() {
  cat("\n[test_single_column]\n")

  input <- matrix(c(1,2,3,4,5), ncol = 1)
  result <- tox_normalize_by_std_dev(input)

  expected <- matrix(rep(1, 5), ncol = 1)

  print(result)
  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Single-column test passed ✓\n")
}

# --------------------------------------------------------------
# Test 3: Zero-row normalization
# --------------------------------------------------------------
test_zero_row <- function() {
  cat("\n[test_zero_row]\n")

  input <- matrix(c(0,0,0,  1,2,2), nrow=2, byrow=TRUE)
  result <- tox_normalize_by_std_dev(input)
  expected <- rms_normalize(input)

  print(result)
  stopifnot(all(abs(result - expected) < 1e-12))

  cat("Zero-row test passed ✓\n")
}

# --------------------------------------------------------------
# Test 4: Performance
# --------------------------------------------------------------
test_rms_performance <- function() {
  cat("\n[test_rms_performance]\n")

  set.seed(42)
  input <- matrix(rnorm(200000), nrow=10000, ncol=20)

  start <- Sys.time()
  result <- tox_normalize_by_std_dev(input)
  elapsed <- as.numeric(Sys.time() - start, "secs")

  cat("Elapsed:", elapsed, "seconds\n")
  stopifnot(elapsed < 1.5)

  cat("Performance test passed ✓\n")
}

# --------------------------------------------------------------
# RUN ALL TESTS
# --------------------------------------------------------------

cat("\n=================================================\n")
cat(" RUNNING RMS NORMALIZATION TESTS\n")
cat("=================================================\n")

test_rms_basic()
test_single_column()
test_zero_row()
test_rms_performance()

cat("\n=================================================\n")
cat("   ALL TESTS PASSED ✓\n")
cat("=================================================\n")
