source("rcpp/tensoromics_functions.R")

test_omics_vector_RAP_projection_call <- function() {
  set.seed(1)

  test_omics_vector_RAP_projection <- function() {
    expr <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)

    vecs_selection_mask <- c(1L, 1L, 1L)                 # length 3 = ncol(expr)
    axes_selection_mask <- c(rep(1L, 3), rep(0L, 7))      # length 10 = nrow(expr)

    stopifnot(length(vecs_selection_mask) == ncol(expr))
    stopifnot(length(axes_selection_mask) == nrow(expr))
    stopifnot(sum(vecs_selection_mask) == 3L)
    stopifnot(sum(axes_selection_mask) == 3L)

    res <- tox_omics_vector_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

    stopifnot(is.matrix(res))
    stopifnot(nrow(res) == 3)  # selected axes
    stopifnot(ncol(res) == 3)  # selected vectors

    cat("test_omics_vector_RAP_projection passed\n")
  }

  test_omics_field_RAP_projection <- function() {
    expr <- matrix(rnorm(6 * 3), nrow = 6, ncol = 3)

    vecs_selection_mask <- c(1L, 1L, 1L)  # length 3 = ncol(expr)
    axes_selection_mask <- c(1L, 1L, 1L)  # length 3 = n_axes = nrow(expr)/2

    stopifnot(nrow(expr) %% 2 == 0)
    stopifnot(length(vecs_selection_mask) == ncol(expr))
    stopifnot(length(axes_selection_mask) == (nrow(expr) / 2))
    stopifnot(sum(vecs_selection_mask) == 3L)
    stopifnot(sum(axes_selection_mask) == 3L)

    res <- tox_omics_field_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

    stopifnot(is.matrix(res))
    stopifnot(nrow(res) == 3)  # selected axes
    stopifnot(ncol(res) == 3)  # selected vectors

    cat("test_omics_field_RAP_projection passed\n")
  }

  test_partial_selection <- function() {
    expr <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)

    vecs_selection_mask <- c(1L, 1L, 0L, 0L, 0L)          # select 2 vectors
    axes_selection_mask <- c(rep(1L, 4), rep(0L, 6))      # select 4 axes

    stopifnot(length(vecs_selection_mask) == ncol(expr))
    stopifnot(length(axes_selection_mask) == nrow(expr))
    stopifnot(sum(vecs_selection_mask) == 2L)
    stopifnot(sum(axes_selection_mask) == 4L)

    res <- tox_omics_vector_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

    stopifnot(is.matrix(res))
    stopifnot(nrow(res) == 4)
    stopifnot(ncol(res) == 2)

    cat("test_partial_selection passed\n")
  }

  cat("=================================================\n")
  cat("    RELATIVE AXIS PLANE R INTERFACE TESTS\n")
  cat("=================================================\n\n")

  test_omics_vector_RAP_projection()
  test_omics_field_RAP_projection()
  test_partial_selection()

  cat("=================================================\n")
  cat("             ALL TESTS COMPLETED\n")
  cat("=================================================\n")
  cat("All RAP R interface tests passed! ✓\n")
}

test_omics_vector_RAP_projection_call()