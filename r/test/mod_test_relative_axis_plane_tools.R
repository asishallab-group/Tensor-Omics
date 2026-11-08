source("r/load_tensor_omics.R")
source("r/test_helpers.R")

set.seed(1)

test_omics_vector_RAP_projection <- function() {
  expr <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)

  vecs_selection_mask <- c(TRUE, TRUE, TRUE)                 # length 3 = ncol(expr)
  axes_selection_mask <- c(rep(TRUE, 3), rep(FALSE, 7))      # length 10 = nrow(expr)

  assert_true(length(vecs_selection_mask) == ncol(expr))
  assert_true(length(axes_selection_mask) == nrow(expr))
  assert_true(sum(vecs_selection_mask) == 3L)
  assert_true(sum(axes_selection_mask) == 3L)

  res <- omics_vector_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

  assert_true(is.matrix(res))
  assert_true(nrow(res) == 3)  # selected axes
  assert_true(ncol(res) == 3)  # selected vectors

}

test_omics_field_RAP_projection <- function() {
  # fields are (n_axes, 2, n_fields): [, 1, ] origins, [, 2, ] targets
  expr <- array(rnorm(3 * 2 * 3), dim = c(3L, 2L, 3L))

  vecs_selection_mask <- c(TRUE, TRUE, TRUE)  # length 3 = ncol(expr)
  axes_selection_mask <- c(TRUE, TRUE, TRUE)  # length 3 = n_axes = nrow(expr)/2

  assert_true(length(vecs_selection_mask) == dim(expr)[3])
  assert_true(length(axes_selection_mask) == dim(expr)[1])
  assert_true(sum(vecs_selection_mask) == 3L)
  assert_true(sum(axes_selection_mask) == 3L)

  res <- omics_field_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

  assert_true(is.matrix(res))
  assert_true(nrow(res) == 3)  # selected axes
  assert_true(ncol(res) == 3)  # selected vectors

}

test_partial_selection <- function() {
  expr <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)

  vecs_selection_mask <- c(TRUE, TRUE, FALSE, FALSE, FALSE)          # select 2 vectors
  axes_selection_mask <- c(rep(TRUE, 4), rep(FALSE, 6))      # select 4 axes

  assert_true(length(vecs_selection_mask) == ncol(expr))
  assert_true(length(axes_selection_mask) == nrow(expr))
  assert_true(sum(vecs_selection_mask) == 2L)
  assert_true(sum(axes_selection_mask) == 4L)

  res <- omics_vector_RAP_projection(expr, vecs_selection_mask, axes_selection_mask)

  assert_true(is.matrix(res))
  assert_true(nrow(res) == 4)
  assert_true(ncol(res) == 2)

}

# `fields` is (n_axes, 2, n_fields); a middle extent of 1 used to be read past (#221)
test_omics_field_RAP_projection_refuses_fields_without_two_points <- function() {
  assert_error(omics_field_RAP_projection(array(1, dim = c(3L, 1L, 2L)), c(TRUE, TRUE),
                                          c(TRUE, TRUE, TRUE)),
               "a field without its target must be refused")
}

# Vectors whose coordinates are all equal. Many of them used to project to a residue of rounding
# noise instead of zeros, rep(0.7, 6) to about -1.1e-16 per coordinate.
.smallest_subnormal <- 2^-1074
.equal_coordinate_cases <- list(
  list(name = "rep(0.7, 6)", value = 0.7, n_axes = 6L),
  list(name = "rep(0.1, 3)", value = 0.1, n_axes = 3L),
  list(name = "rep(0.1, 7)", value = 0.1, n_axes = 7L),
  list(name = "rep(3.3, 13)", value = 3.3, n_axes = 13L),
  list(name = "rep(1e300, 3)", value = 1e300, n_axes = 3L),
  list(name = "rep(-2.5, 4)", value = -2.5, n_axes = 4L),
  list(name = "rep(tiny, 5)", value = .Machine$double.xmin, n_axes = 5L),
  list(name = "rep(tiny / 4, 3), subnormal", value = .Machine$double.xmin / 4, n_axes = 3L),
  list(name = "rep(smallest subnormal, 7)", value = .smallest_subnormal, n_axes = 7L),
  list(name = "5 on a single axis", value = 5, n_axes = 1L)
)

.assert_exact_zeros_that_normalization_refuses <- function(projection, case_name) {
  assert_true(all(projection == 0), paste0(case_name, ": projection is not exact zeros"))
  assert_error(normalize_unit_length(projection),
               paste0(case_name, ": normalize_unit_length must refuse the zero vector"),
               code = ERR_DIVISION_BY_ZERO)
}

test_omics_vector_RAP_projection_equal_coordinates_give_exact_zeros <- function() {
  for (case in .equal_coordinate_cases) {
    vecs <- matrix(case$value, nrow = case$n_axes, ncol = 1L)
    res <- omics_vector_RAP_projection(vecs, TRUE, rep(TRUE, case$n_axes))
    .assert_exact_zeros_that_normalization_refuses(res[, 1], case$name)
  }
}

test_omics_field_RAP_projection_equal_shifts_give_exact_zeros <- function() {
  for (case in .equal_coordinate_cases) {
    # the origin holds `value`, the target 0, so every shift is `value`
    fields <- array(0, dim = c(case$n_axes, 2L, 1L))
    fields[, 1L, 1L] <- case$value
    res <- omics_field_RAP_projection(fields, TRUE, rep(TRUE, case$n_axes))
    .assert_exact_zeros_that_normalization_refuses(res[, 1], paste("shift", case$name))
  }
}

run_all_tests()
