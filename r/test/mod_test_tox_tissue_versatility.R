# The `tox_tissue_versatility` suite: each published procedure can be called from R, returns the
# documented type and shape, and raises the documented error. Whether the numbers are right is the
# Fortran suite's job (test/mod_test_tox_tissue_versatility.F90), as the coding guide's "Where to
# Test What" asks.
#
# The one value needs no arithmetic beyond the trivial: over four selected axes a uniform vector
# scores 0 at an angle of 0 and a single-axis vector scores 1, both exactly. It only shows that the
# vectors mask selects and keeps the order, and that the axes mask reaches the library the right
# way round: the fifth axis, which the mask leaves out, would break every uniform vector.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_AXES <- 5L
N_VECTORS <- 3L

# every axis but the last
AXES_SELECTION_MASK <- c(TRUE, TRUE, TRUE, TRUE, FALSE)

# (axes, vectors). Over the first four axes, vectors 1 and 3 are uniform and vector 2 sits on axis
# 4 alone; the fifth axis breaks the uniform ones, and a reversed axes mask would read it.
.expression_vectors <- function() {
  cbind(c(1, 1, 1, 1, 9),
        c(0, 0, 0, 3, 0),
        c(4, 4, 4, 4, 1))
}

# -----------------------------------------------------------------------------------------------
# compute_tissue_versatility
# -----------------------------------------------------------------------------------------------

test_compute_tissue_versatility <- function() {
  result <- compute_tissue_versatility(.expression_vectors(), c(FALSE, TRUE, TRUE), AXES_SELECTION_MASK)
  assert_true(is.list(result) && identical(names(result), c("tissue_versatilities", "tissue_angles_deg")),
              "expected a list of tissue_versatilities and tissue_angles_deg")
  for (name in names(result))
    assert_true(is.double(result[[name]]) && is.null(dim(result[[name]])) && length(result[[name]]) == 2L,
                sprintf("%s: expected a double vector, one per selected vector", name))
  # the one value: vector 2 (on one axis) then vector 3 (uniform). A reversed vectors mask would
  # give c(0, 1), a reversed axes mask or none at all a nonzero second entry
  assert_true(identical(result$tissue_versatilities, c(1, 0)),
              paste("tissue_versatilities: got", toString(result$tissue_versatilities)))
  assert_true(identical(result$tissue_angles_deg[2], 0), paste("tissue_angles_deg: got", toString(result$tissue_angles_deg)))
}

test_compute_tissue_versatility_rejects_no_selected_vectors <- function() {
  # an empty selection is an error, not an empty result; the error names the output
  assert_error(compute_tissue_versatility(.expression_vectors(), rep(FALSE, N_VECTORS), AXES_SELECTION_MASK),
               "no selected vectors", ERR_EMPTY_INPUT)
}

test_compute_tissue_versatility_rejects_no_selected_axes <- function() {
  assert_error(compute_tissue_versatility(.expression_vectors(), rep(TRUE, N_VECTORS), rep(FALSE, N_AXES)),
               "no selected axes", ERR_INVALID_INPUT)
}

test_compute_tissue_versatility_rejects_no_axes <- function() {
  assert_error(compute_tissue_versatility(matrix(numeric(0), nrow = 0, ncol = N_VECTORS), rep(TRUE, N_VECTORS),
                                          logical(0)),
               "no axes", ERR_EMPTY_INPUT)
}

test_compute_tissue_versatility_rejects_no_vectors <- function() {
  assert_error(compute_tissue_versatility(matrix(numeric(0), nrow = N_AXES, ncol = 0), logical(0),
                                          AXES_SELECTION_MASK),
               "no vectors", ERR_EMPTY_INPUT)
}

test_compute_tissue_versatility_rejects_na <- function() {
  # R's missing value is a NaN to the library; it sits in vector 1, which the mask leaves out
  expression_vectors <- .expression_vectors()
  expression_vectors[3, 1] <- NA
  assert_error(compute_tissue_versatility(expression_vectors, c(FALSE, TRUE, TRUE), AXES_SELECTION_MASK),
               "a missing expression", ERR_NAN_INF)
}

test_compute_tissue_versatility_rejects_infinity <- function() {
  # the infinity sits on axis 5, which the mask leaves out: the whole matrix is checked
  expression_vectors <- .expression_vectors()
  expression_vectors[5, 3] <- Inf
  assert_error(compute_tissue_versatility(expression_vectors, c(FALSE, TRUE, TRUE), AXES_SELECTION_MASK),
               "an infinite expression", ERR_NAN_INF)
}

test_compute_tissue_versatility_rejects_a_missing_vector_flag <- function() {
  # the binding's own NA check, before the library is called
  assert_error(compute_tissue_versatility(.expression_vectors(), c(TRUE, NA, TRUE), AXES_SELECTION_MASK),
               "an NA in vectors_selection_mask")
}

test_compute_tissue_versatility_rejects_a_missing_axis_flag <- function() {
  # the binding's own NA check, before the library is called
  assert_error(compute_tissue_versatility(.expression_vectors(), rep(TRUE, N_VECTORS),
                                          c(TRUE, TRUE, NA, TRUE, FALSE)),
               "an NA in axes_selection_mask")
}

test_compute_tissue_versatility_rejects_a_short_vectors_selection_mask <- function() {
  # the binding's own extent check, before the library is called
  assert_error(compute_tissue_versatility(.expression_vectors(), rep(TRUE, N_VECTORS - 1L), AXES_SELECTION_MASK),
               "vectors_selection_mask shorter than n_vectors")
}

test_compute_tissue_versatility_rejects_a_long_axes_selection_mask <- function() {
  # the binding's own extent check, before the library is called
  assert_error(compute_tissue_versatility(.expression_vectors(), rep(TRUE, N_VECTORS), rep(TRUE, N_AXES + 1L)),
               "axes_selection_mask longer than n_axes")
}

test_compute_tissue_versatility_rejects_a_vector <- function() {
  # the binding's own type check, before the library is called
  assert_error(compute_tissue_versatility(c(1, 2, 3), TRUE, c(TRUE, TRUE, TRUE)),
               "expression_vectors must be a matrix")
}

run_all_tests()
