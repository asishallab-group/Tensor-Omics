source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

set.seed(1)

test_omics_vector_RAP_projection <- function() {
  expr <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)

  vecs_selection_mask <- c(1L, 1L, 1L)                 # length 3 = ncol(expr)
  axes_selection_mask <- c(rep(1L, 3), rep(0L, 7))      # length 10 = nrow(expr)

  assert_true(length(vecs_selection_mask) == ncol(expr))
  assert_true(length(axes_selection_mask) == nrow(expr))
  assert_true(sum(vecs_selection_mask) == 3L)
  assert_true(sum(axes_selection_mask) == 3L)

  res <- tox_omics_vector_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

  assert_true(is.matrix(res))
  assert_true(nrow(res) == 3)  # selected axes
  assert_true(ncol(res) == 3)  # selected vectors

}

test_omics_field_RAP_projection <- function() {
  expr <- matrix(rnorm(6 * 3), nrow = 6, ncol = 3)

  vecs_selection_mask <- c(1L, 1L, 1L)  # length 3 = ncol(expr)
  axes_selection_mask <- c(1L, 1L, 1L)  # length 3 = n_axes = nrow(expr)/2

  assert_true(nrow(expr) %% 2 == 0)
  assert_true(length(vecs_selection_mask) == ncol(expr))
  assert_true(length(axes_selection_mask) == (nrow(expr) / 2))
  assert_true(sum(vecs_selection_mask) == 3L)
  assert_true(sum(axes_selection_mask) == 3L)

  res <- tox_omics_field_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

  assert_true(is.matrix(res))
  assert_true(nrow(res) == 3)  # selected axes
  assert_true(ncol(res) == 3)  # selected vectors

}

test_partial_selection <- function() {
  expr <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)

  vecs_selection_mask <- c(1L, 1L, 0L, 0L, 0L)          # select 2 vectors
  axes_selection_mask <- c(rep(1L, 4), rep(0L, 6))      # select 4 axes

  assert_true(length(vecs_selection_mask) == ncol(expr))
  assert_true(length(axes_selection_mask) == nrow(expr))
  assert_true(sum(vecs_selection_mask) == 2L)
  assert_true(sum(axes_selection_mask) == 4L)

  res <- tox_omics_vector_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

  assert_true(is.matrix(res))
  assert_true(nrow(res) == 4)
  assert_true(ncol(res) == 2)

}

run_all_tests()
