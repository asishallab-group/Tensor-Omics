source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

#' Comprehensive test function
test_bst_functions <- function() {
  # Test BST
  x <- c(3.0, 1.0, 4.0, 2.0)
  bst_ix <- build_bst_index(x)$bst_index
  assert_equal_int(bst_ix, c(2L, 4L, 1L, 3L))
  
  # Test get_sorted_value
  sorted_val <- get_sorted_value(x, bst_ix, 2)
  assert_equal_numeric(sorted_val, 2.0, 0.0, "Expected output match exactly the second smallest value")
  
  # Test BST range query
  range_result <- bst_range_query(x, bst_ix, 1.5, 3.5)
  slice <- x[range_result$indices[seq_len(range_result$count)]]
  assert_equal_numeric(slice, c(2.0, 3.0), 0.0, "Expected output match exactly the values in range 1.5,3.5")
}

test_kd_functions <- function() {
  # Test KD-Tree
  X <- matrix(c(1.0, 4.0, 2.0, 5.0, 3.0, 6.0), nrow = 2, ncol = 3)
  kd_ix <- build_kd_index(X, c(1, 2))
  
  # Test get_kd_point
  point <- get_kd_point(X, kd_ix, 2)
  
  # Test Spherical KD-Tree
  V <- matrix(rnorm(6), nrow = 2, ncol = 3)
  V <- V / sqrt(colSums(V^2))  # Normalize to unit length
  sphere_ix <- build_spherical_kd(V, c(1, 2))
  
  # Test get_kd_point for spherical
  sphere_point <- get_kd_point(V, sphere_ix, 1)
  force(sphere_point)
  
  # Empty cases
  assert_error(build_bst_index(numeric(0)), "Expected an error for empty bst input")
  
  # Single element
  single_bst <- build_bst_index(c(42.0))$bst_index
  assert_equal_int(single_bst, c(1L), "Single element should have only 1 in index")
}

run_all_tests()
