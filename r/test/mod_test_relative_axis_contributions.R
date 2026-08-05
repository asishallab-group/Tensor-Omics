## ==================== TESTS FOR RELATIVE AXIS CONTRIBUTIONS ====================
# Load the compiled Fortran library
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

#' Assert that all values in a vector are in [0, 1]
assert_in_range <- function(vec, min_val = 0, max_val = 1, message = "") {
  assert_false(any(vec < min_val | vec > max_val), sprintf("Test failed: %s. Values out of range [%g, %g]", message, min_val, max_val))
}

#' Assert that the sum of a vector is approximately equal to a value
assert_sum_equal <- function(vec, expected_sum, tolerance = 1e-12, message = "") {
  actual_sum <- sum(vec)
  assert_equal_numeric(actual_sum, expected_sum, sprintf("Test failed: %s. Expected sum %f, got %f", message, expected_sum, actual_sum))
}


#' Test: relative_axes_changes_from_shift_vector
test_shift_positive_vector <- function() {
  vec <- c(1, 2, 3)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift positive vector: sum")
  assert_in_range(contrib, message = "shift positive vector: range")
}

test_shift_negative_vector <- function() {
  vec <- c(-1, -2, -3)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift negative vector: sum")
  assert_in_range(contrib, message = "shift negative vector: range")
}

test_shift_mixed_vector <- function() {
  vec <- c(2, -2, 4)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift mixed vector: sum")
  assert_in_range(contrib, message = "shift mixed vector: range")
}

test_shift_zero_vector <- function() {
  # as above: no direction to apportion, so this is an error, not zeros
  vec <- c(0, 0, 0)
  assert_error(relative_axes_changes_from_shift_vector(vec),
               "Expected an error for a zero shift vector", ERR_DIVISION_BY_ZERO)
}

test_shift_one_nonzero_axis <- function() {
  vec <- c(0, 5, 0)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift one nonzero axis: sum")
  assert_in_range(contrib, message = "shift one nonzero axis: range")
}

test_shift_all_equal <- function() {
  vec <- c(2, 2, 2)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift all equal: sum")
  assert_in_range(contrib, message = "shift all equal: range")
}

test_shift_large_vector <- function() {
  vec <- rep(1, 100)
  contrib <- relative_axes_changes_from_shift_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "shift large vector: sum")
  assert_in_range(contrib, message = "shift large vector: range")
}

#' Test: relative_axes_expression_from_expression_vector
test_expr_positive_vector <- function() {
  vec <- c(1, 2, 3)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr positive vector: sum")
  assert_in_range(contrib, message = "expr positive vector: range")
}

test_expr_negative_vector <- function() {
  vec <- c(-1, -2, -3)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr negative vector: sum")
  assert_in_range(contrib, message = "expr negative vector: range")
}

test_expr_mixed_vector <- function() {
  vec <- c(2, -2, 4)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr mixed vector: sum")
  assert_in_range(contrib, message = "expr mixed vector: range")
}

test_expr_zero_vector <- function() {
  # a zero vector has no direction to apportion, so the Fortran reports a division by
  # zero rather than returning zeros. The Python suite asserts the same.
  vec <- c(0, 0, 0)
  assert_error(relative_axes_expression_from_expression_vector(vec),
               "Expected an error for a zero expression vector", ERR_DIVISION_BY_ZERO)
}

test_expr_one_nonzero_axis <- function() {
  vec <- c(0, 5, 0)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr one nonzero axis: sum")
  assert_in_range(contrib, message = "expr one nonzero axis: range")
}

test_expr_all_equal <- function() {
  vec <- c(2, 2, 2)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr all equal: sum")
  assert_in_range(contrib, message = "expr all equal: range")
}

test_expr_large_vector <- function() {
  vec <- rep(1, 100)
  contrib <- relative_axes_expression_from_expression_vector(vec)
  assert_sum_equal(contrib, 1.0, message = "expr large vector: sum")
  assert_in_range(contrib, message = "expr large vector: range")
}

run_all_tests()
