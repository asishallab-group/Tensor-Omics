source("r/tensoromics_functions.R")

#' Comprehensive test function
test_tree_functions <- function() {
  cat("=== Testing BST Functions ===\n")
  
  # Test BST
  x <- c(3.0, 1.0, 4.0, 2.0)
  bst_ix <- build_bst_index(x)
  cat("BST indices:", bst_ix, "\n")
  cat("Sorted values:", x[bst_ix], "\n")
  
  # Test get_sorted_value
  sorted_val <- get_sorted_value(x, bst_ix, 2)
  cat("2nd smallest value:", sorted_val, "\n")
  
  # Test BST range query
  range_result <- bst_range_query(x, bst_ix, 1.5, 3.5)
  cat("Range [1.5, 3.5] matches:", range_result$indices, "\n")
  cat("Matching values:", x[range_result$indices], "\n")
  
  cat("\n=== Testing KD-Tree Functions ===\n")
  
  # Test KD-Tree
  X <- matrix(c(1.0, 4.0, 2.0, 5.0, 3.0, 6.0), nrow = 2, ncol = 3)
  kd_ix <- build_kd_index(X, c(1, 2))
  cat("KD-Tree indices:", kd_ix, "\n")
  cat("Points in order:\n")
  print(X[, kd_ix, drop = FALSE])
  
  # Test get_kd_point
  point <- get_kd_point(X, kd_ix, 2)
  cat("2nd point in KD order:", point, "\n")
  
  cat("\n=== Testing Spherical KD-Tree Functions ===\n")
  
  # Test Spherical KD-Tree
  V <- matrix(rnorm(6), nrow = 2, ncol = 3)
  V <- V / sqrt(colSums(V^2))  # Normalize to unit length
  sphere_ix <- build_spherical_kd(V, c(1, 2))
  cat("Spherical KD-Tree indices:", sphere_ix, "\n")
  cat("Vectors in spherical order:\n")
  print(V[, sphere_ix, drop = FALSE])
  
  # Test get_kd_point for spherical
  sphere_point <- get_kd_point(V, sphere_ix, 1)
  cat("1st spherical vector:", sphere_point, "\n")
  cat("Norm:", sqrt(sum(sphere_point^2)), "\n")
  
  cat("\n=== Testing Edge Cases ===\n")
  
  # Empty cases
  empty_bst <- tryCatch(build_bst_index(numeric(0)), error = function(e) e)
  cat("Empty BST test:", ifelse(inherits(empty_bst, "error"), "Error", "Success"), "\n")
  
  # Single element
  single_bst <- build_bst_index(c(42.0))
  cat("Single element BST:", single_bst, "\n")
}

#' Memory cleanup function (optional)
cleanup_tree_memory <- function() {
  # Force garbage collection to free memory
  gc()
}

# Export functions for package use
if (sys.nframe() == 0) {
  test_tree_functions()
}
test_tree_functions()