## ==================== TESTS FOR RELATIVE AXIS CONTRIBUTIONS ====================
# Load the compiled Fortran library
source("r/tensoromics_functions.R")

#' Assert that two numeric vectors are approximately equal
assert_vec_equal <- function(actual, expected, tolerance = 1e-12, message = "") {
  if (length(actual) != length(expected)) {
    stop(sprintf("Test failed: %s. Length mismatch: expected %d, got %d", message, length(expected), length(actual)))
  }
  diffs <- abs(actual - expected)
  if (any(diffs > tolerance)) {
    stop(sprintf("Test failed: %s. Max diff: %e", message, max(diffs)))
  }
  cat(sprintf("✓ %s\n", message))
}

#' Assert that all values in a vector are finite
assert_all_finite <- function(vec, message = "") {
  if (any(!is.finite(vec))) {
    stop(sprintf("Test failed: %s. Non-finite values detected.", message))
  }
  cat(sprintf("✓ %s\n", message))
}

#' Assert that all values in a vector are in [0, 1]
assert_in_range <- function(vec, min_val = 0, max_val = 1, message = "") {
  if (any(vec < min_val | vec > max_val)) {
    stop(sprintf("Test failed: %s. Values out of range [%g, %g]", message, min_val, max_val))
  }
  cat(sprintf("✓ %s\n", message))
}

#' Assert that the sum of a vector is approximately equal to a value
assert_sum_equal <- function(vec, expected_sum, tolerance = 1e-12, message = "") {
  actual_sum <- sum(vec)
  if (abs(actual_sum - expected_sum) > tolerance) {
    stop(sprintf("Test failed: %s. Expected sum %f, got %f", message, expected_sum, actual_sum))
  }
  cat(sprintf("✓ %s\n", message))
}


#' Test: relative_axes_changes_from_shift_vector
test_shift_positive_vector <- function() {
  vec <- c(1, 2, 3)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift positive vector: sum")
  assert_all_finite(contrib, "shift positive vector: finite")
  assert_in_range(contrib, message = "shift positive vector: range")
}

test_shift_negative_vector <- function() {
  vec <- c(-1, -2, -3)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift negative vector: sum")
  assert_all_finite(contrib, "shift negative vector: finite")
  assert_in_range(contrib, message = "shift negative vector: range")
}

test_shift_mixed_vector <- function() {
  vec <- c(2, -2, 4)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift mixed vector: sum")
  assert_all_finite(contrib, "shift mixed vector: finite")
  assert_in_range(contrib, message = "shift mixed vector: range")
}

test_shift_zero_vector <- function() {
  vec <- c(0, 0, 0)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 0.0, message = "shift zero vector: sum")
  assert_all_finite(contrib, "shift zero vector: finite")
  assert_in_range(contrib, message = "shift zero vector: range")
}

test_shift_one_nonzero_axis <- function() {
  vec <- c(0, 5, 0)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift one nonzero axis: sum")
  assert_all_finite(contrib, "shift one nonzero axis: finite")
  assert_in_range(contrib, message = "shift one nonzero axis: range")
}

test_shift_all_equal <- function() {
  vec <- c(2, 2, 2)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift all equal: sum")
  assert_all_finite(contrib, "shift all equal: finite")
  assert_in_range(contrib, message = "shift all equal: range")
}

test_shift_large_vector <- function() {
  vec <- rep(1, 100)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift large vector: sum")
  assert_all_finite(contrib, "shift large vector: finite")
  assert_in_range(contrib, message = "shift large vector: range")
}

#' Test: relative_axes_expression_from_expression_vector
test_expr_positive_vector <- function() {
  vec <- c(1, 2, 3)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr positive vector: sum")
  assert_all_finite(contrib, "expr positive vector: finite")
  assert_in_range(contrib, message = "expr positive vector: range")
}

test_expr_negative_vector <- function() {
  vec <- c(-1, -2, -3)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr negative vector: sum")
  assert_all_finite(contrib, "expr negative vector: finite")
  assert_in_range(contrib, message = "expr negative vector: range")
}

test_expr_mixed_vector <- function() {
  vec <- c(2, -2, 4)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr mixed vector: sum")
  assert_all_finite(contrib, "expr mixed vector: finite")
  assert_in_range(contrib, message = "expr mixed vector: range")
}

test_expr_zero_vector <- function() {
  vec <- c(0, 0, 0)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 0.0, message = "expr zero vector: sum")
  assert_all_finite(contrib, "expr zero vector: finite")
  assert_in_range(contrib, message = "expr zero vector: range")
}

test_expr_one_nonzero_axis <- function() {
  vec <- c(0, 5, 0)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr one nonzero axis: sum")
  assert_all_finite(contrib, "expr one nonzero axis: finite")
  assert_in_range(contrib, message = "expr one nonzero axis: range")
}

test_expr_all_equal <- function() {
  vec <- c(2, 2, 2)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr all equal: sum")
  assert_all_finite(contrib, "expr all equal: finite")
  assert_in_range(contrib, message = "expr all equal: range")
}

test_expr_large_vector <- function() {
  vec <- rep(1, 100)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr large vector: sum")
  assert_all_finite(contrib, "expr large vector: finite")
  assert_in_range(contrib, message = "expr large vector: range")
}


cat("Running R relative axis contribution tests...\n")
cat("============================================\n")
test_functions <- list(
# shift vector tests
test_shift_positive_vector,
test_shift_negative_vector,
test_shift_mixed_vector,
test_shift_zero_vector,
test_shift_one_nonzero_axis,
test_shift_all_equal,
test_shift_large_vector,
# expression vector tests
test_expr_positive_vector,
test_expr_negative_vector,
test_expr_mixed_vector,
test_expr_zero_vector,
test_expr_one_nonzero_axis,
test_expr_all_equal,
test_expr_large_vector
)
test_names <- c(
"test_shift_positive_vector",
"test_shift_negative_vector",
"test_shift_mixed_vector",
"test_shift_zero_vector",
"test_shift_one_nonzero_axis",
"test_shift_all_equal",
"test_shift_large_vector",
"test_expr_positive_vector",
"test_expr_negative_vector",
"test_expr_mixed_vector",
"test_expr_zero_vector",
"test_expr_one_nonzero_axis",
"test_expr_all_equal",
"test_expr_large_vector"
)
passed <- 0
failed <- 0
for (i in seq_along(test_functions)) {
cat(sprintf("\n--- %s ---\n", test_names[i]))
tryCatch({
    test_functions[[i]]()
    passed <- passed + 1
}, error = function(e) {
    cat(sprintf("✗ FAILED: %s\n", e$message))
    failed <- failed + 1
})
}
cat("\n============================================\n")
cat(sprintf("Tests completed: %d passed, %d failed\n", passed, failed))
if (failed == 0) {
  cat("🎉 All R relative axis contribution tests passed!\n")
} else {
  cat("❌ Some tests failed. Please check the output above.\n")
}

