# =====================================================================
#   TESTS FOR calc_fchange_rcpp — VALID BEHAVIOR ONLY
# =====================================================================

library(Rcpp)
source("rcpp/tensoromics_functions.R")

cat("=== TESTING: calc_fchange_rcpp ===\n")

# --------------------------------------------------------------
# Test 1: Basic test
# --------------------------------------------------------------
test_fc_basic <- function() {
  cat("\n[test_fc_basic] Basic test\n")

  input <- matrix(c(1,2,3,
                    4,5,6), nrow=3)

  result <- calc_fchange_rcpp(input, c(1), c(2))$output
  expected <- c(3,3,3)

  print(result)
  print(expected)

  stopifnot(all(result == expected))

  cat("Basic fold-change test passed ✓\n")
}

# --------------------------------------------------------------
# Test 2: Multiple pairs
# --------------------------------------------------------------
test_fc_multiple_pairs <- function() {
  cat("\n[test_fc_multiple_pairs]\n")

  input <- matrix(c(1,2,
                    4,5), nrow=2)

  result <- calc_fchange_rcpp(input, c(1,1), c(2,2))$output
  expected <- matrix(c(3,3,3,3), nrow=2)

  print(result)
  print(expected)

  stopifnot(all(result == expected))

  cat("Multiple pairs test passed ✓\n")
}

# --------------------------------------------------------------
# Test 3: Negative values
# --------------------------------------------------------------
test_fc_negative_values <- function() {
  cat("\n[test_fc_negative_values]\n")

  input <- matrix(c(-1, 2, -3,
                     4, -2, 1), nrow=3)

  result <- calc_fchange_rcpp(input, c(1), c(2))$output
  expected <- input[,2] - input[,1]

  print(result)
  print(expected)

  stopifnot(all(result == expected))

  cat("Negative values test passed ✓\n")
}

# --------------------------------------------------------------
# Test 4: Invalid inputs — only check that it does NOT crash
# --------------------------------------------------------------
test_fc_invalid_inputs <- function() {
  cat("\n[test_fc_invalid_inputs]\n")

  input <- matrix(1:9, nrow=3)

  # Case 1: uneven length (should NOT crash)
  calc_fchange_rcpp(input, c(1), c(2,3))

  # Case 2: out-of-bound (should NOT crash)
  calc_fchange_rcpp(input, c(10), c(1))

  # Case 3: zero-length (should NOT crash)
  calc_fchange_rcpp(input, integer(0), integer(0))

  cat("Invalid inputs test passed (function does not crash) ✓\n")
}

# --------------------------------------------------------------
# RUN ALL TESTS
# --------------------------------------------------------------

cat("\n=================================================\n")
cat(" RUNNING calc_fchange TESTS\n")
cat("=================================================\n")

test_fc_basic()
test_fc_multiple_pairs()
test_fc_negative_values()
test_fc_invalid_inputs()

cat("\n=================================================\n")
cat("   ALL TESTS PASSED ✓\n")
cat("=================================================\n")
