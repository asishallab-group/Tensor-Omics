# The `tox_relative_axis_plane_tools` suite: each published procedure can be called from R,
# returns the documented type and shape, and raises the documented error. Whether the numbers are
# right is the Fortran suite's job (test/mod_test_tox_relative_axis_plane_tools.F90), as the
# coding guide's "Where to Test What" asks. The one hand-derived value per procedure below proves
# only that the binding passes its data through unchanged: argument order, layout, selection.
#
# The binding derives every `n_selected_*` from its selection mask, so the Fortran check that a
# count agrees with its mask cannot be reached from here; an empty selection can.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_AXES <- 4L
N_VECS <- 5L
AXES_SELECTION_MASK <- c(TRUE, FALSE, TRUE, TRUE)
VECS_SELECTION_MASK <- c(TRUE, TRUE, FALSE, TRUE, FALSE)
N_SELECTED_AXES <- 3L
N_SELECTED_VECS <- 3L

# the three contribution procedures share one signature and one contract
CONTRIBUTION_PROCEDURES <- list(
  compute_relative_axis_contributions = compute_relative_axis_contributions,
  relative_axes_changes_from_shift_vector = relative_axes_changes_from_shift_vector,
  relative_axes_expression_from_expression_vector = relative_axes_expression_from_expression_vector
)

# what the Fortran rejects in every real input; NA_real_ is a NaN payload, which the binding
# leaves to the Fortran on purpose
NON_FINITE_VALUES <- c(NaN, Inf, NA_real_)

# Vectors whose coordinates are all equal. Many of them used to project to a residue of rounding
# noise instead of zeros, rep(0.1, 7) to about -1.4e-17 per coordinate.
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

# (axes, vectors), finite and all different, so nothing degenerates
.vecs <- function() {
  outer(seq_len(N_AXES), seq_len(N_VECS), function(i_axis, i_vec) i_axis * (i_vec + 1) + (i_axis - 1)^2)
}

# (axes, origin/shift, fields): origins from `.vecs`, and shifts off them
.fields <- function() {
  fields <- array(0, dim = c(N_AXES, 2L, N_VECS))
  fields[, 1, ] <- .vecs()
  fields[, 2, ] <- .vecs()^1.5
  fields
}

# `x` with its first element replaced by `value`
.with <- function(x, value) {
  x[1] <- value
  x
}

.assert_matrix <- function(result, dims, name) {
  assert_true(is.matrix(result), paste0(name, ": expected a matrix"))
  assert_true(is.double(result), paste0(name, ": expected a double matrix"))
  assert_true(identical(dim(result), as.integer(dims)),
              sprintf("%s: expected dim %s, got %s", name, toString(dims), toString(dim(result))))
  assert_true(all(is.finite(result)), paste0(name, ": non-finite values in the result"))
}

.assert_vector <- function(result, len, name) {
  assert_true(is.double(result) && is.null(dim(result)), paste0(name, ": expected a double vector"))
  assert_true(length(result) == len, sprintf("%s: expected length %d, got %d", name, len, length(result)))
  assert_true(all(is.finite(result)), paste0(name, ": non-finite values in the result"))
}

.assert_exact_zeros_that_normalization_refuses <- function(projection, case_name) {
  assert_true(all(projection == 0), paste0(case_name, ": projection is not exact zeros"))
  assert_error(normalize_unit_length(projection),
               paste0(case_name, ": normalize_unit_length must refuse the zero vector"),
               code = ERR_DIVISION_BY_ZERO)
}

# ---- omics_vector_RAP_projection -------------------------------------------------------------

test_omics_vector_RAP_projection <- function() {
  .assert_matrix(omics_vector_RAP_projection(.vecs(), VECS_SELECTION_MASK, AXES_SELECTION_MASK),
                 c(N_SELECTED_AXES, N_SELECTED_VECS), "omics_vector_RAP_projection")
}

test_omics_vector_RAP_projection_keeps_axes_as_rows <- function() {
  # pass-through only: the projection subtracts the selected column's mean, 30, from
  # (10, 20, 60); a transposed matrix or a misapplied mask would give other numbers
  vecs <- matrix(c(1, 2, 3, 10, 20, 60), nrow = 3)
  result <- omics_vector_RAP_projection(vecs, c(FALSE, TRUE), c(TRUE, TRUE, TRUE))
  assert_equal_numeric(as.vector(result), c(-20, -10, 30), msg = paste("got", toString(result)))
}

test_omics_vector_RAP_projection_rejects_a_mask_of_the_wrong_length <- function() {
  # the binding's own shape check, before the library is called
  assert_error(omics_vector_RAP_projection(.vecs(), VECS_SELECTION_MASK, AXES_SELECTION_MASK[-1]),
               "vecs has four axes, the mask three")
}

test_omics_vector_RAP_projection_rejects_a_mask_that_is_not_logical <- function() {
  # the binding's own type check: only the numbers 0 and 1 are accepted as a mask
  assert_error(omics_vector_RAP_projection(.vecs(), as.integer(VECS_SELECTION_MASK) + 1L, AXES_SELECTION_MASK),
               "a mask of ones and twos")
  assert_error(omics_vector_RAP_projection(.vecs(), .with(VECS_SELECTION_MASK, NA), AXES_SELECTION_MASK),
               "a mask with NA")
}

test_omics_vector_RAP_projection_rejects_non_finite_vecs <- function() {
  for (value in NON_FINITE_VALUES) {
    assert_error(omics_vector_RAP_projection(.with(.vecs(), value), VECS_SELECTION_MASK, AXES_SELECTION_MASK),
                 paste("vecs holds", value), ERR_NAN_INF)
  }
}

test_omics_vector_RAP_projection_rejects_an_empty_selection <- function() {
  assert_error(omics_vector_RAP_projection(.vecs(), VECS_SELECTION_MASK, rep(FALSE, N_AXES)),
               "no axis selected", ERR_EMPTY_INPUT)
}

test_omics_vector_RAP_projection_equal_coordinates_give_exact_zeros <- function() {
  for (case in .equal_coordinate_cases) {
    vecs <- matrix(case$value, nrow = case$n_axes, ncol = 1L)
    res <- omics_vector_RAP_projection(vecs, TRUE, rep(TRUE, case$n_axes))
    .assert_exact_zeros_that_normalization_refuses(res[, 1], case$name)
  }
}

# ---- omics_field_RAP_projection --------------------------------------------------------------

test_omics_field_RAP_projection <- function() {
  .assert_matrix(omics_field_RAP_projection(.fields(), VECS_SELECTION_MASK, AXES_SELECTION_MASK),
                 c(N_SELECTED_AXES, N_SELECTED_VECS), "omics_field_RAP_projection")
}

test_omics_field_RAP_projection_projects_the_selected_shift <- function() {
  # pass-through only: each field is [origin, shift], and the shift is what is projected. The
  # selected second field's shift (1, 2, 6) is mean-free (-2, -1, 3), whatever its origin (100
  # here); the first field is never selected
  fields <- array(0, dim = c(3L, 2L, 2L))
  fields[, , 1] <- 100
  fields[, 1, 2] <- 100
  fields[, 2, 2] <- c(1, 2, 6)
  result <- as.vector(omics_field_RAP_projection(fields, c(FALSE, TRUE), c(TRUE, TRUE, TRUE)))
  expected <- c(-2, -1, 3)
  assert_true(isTRUE(all.equal(result, expected)), paste("expected -2, -1, 3, got", toString(result)))
}

test_omics_field_RAP_projection_rejects_a_mask_of_the_wrong_length <- function() {
  # the binding's own shape check, before the library is called
  assert_error(omics_field_RAP_projection(.fields(), VECS_SELECTION_MASK[-1], AXES_SELECTION_MASK),
               "fields has five fields, the mask four")
}

# `fields` is (n_axes, 2, n_fields); a middle extent of 1 used to be read past (#221)
test_omics_field_RAP_projection_rejects_fields_without_two_points <- function() {
  assert_error(omics_field_RAP_projection(.fields()[, 1, , drop = FALSE], VECS_SELECTION_MASK, AXES_SELECTION_MASK),
               "a field without its shift must be refused")
}

test_omics_field_RAP_projection_rejects_non_finite_fields <- function() {
  for (value in NON_FINITE_VALUES) {
    assert_error(omics_field_RAP_projection(.with(.fields(), value), VECS_SELECTION_MASK, AXES_SELECTION_MASK),
                 paste("fields holds", value), ERR_NAN_INF)
  }
}

test_omics_field_RAP_projection_rejects_an_empty_selection <- function() {
  assert_error(omics_field_RAP_projection(.fields(), rep(FALSE, N_VECS), AXES_SELECTION_MASK),
               "no field selected", ERR_EMPTY_INPUT)
}

test_omics_field_RAP_projection_equal_shifts_give_exact_zeros <- function() {
  for (case in .equal_coordinate_cases) {
    fields <- array(0, dim = c(case$n_axes, 2L, 1L))
    fields[, 1L, 1L] <- seq_len(case$n_axes)  # an origin off the diagonal, which does not enter
    fields[, 2L, 1L] <- case$value  # the shift
    res <- omics_field_RAP_projection(fields, TRUE, rep(TRUE, case$n_axes))
    .assert_exact_zeros_that_normalization_refuses(res[, 1], paste("shift", case$name))
  }
}

# ---- clock_hand_angle_between_vectors --------------------------------------------------------

test_clock_hand_angle_between_vectors <- function() {
  result <- clock_hand_angle_between_vectors(c(1, 0, 0), c(0.6, 0.8, 0), c(0, 1, 0))
  .assert_vector(result, 1L, "clock_hand_angle_between_vectors")
}

test_clock_hand_angle_between_vectors_keeps_the_argument_order <- function() {
  # pass-through only: the turn from e1 to e2 moves along e2, which the reference (-1, 1)
  # counts positive, so +pi/2; swapped vectors turn along e1, which it counts negative
  result <- clock_hand_angle_between_vectors(c(1, 0), c(0, 1), c(-1, 1))
  assert_equal_numeric(result, pi / 2, msg = paste("expected pi/2, got", result))
}

test_clock_hand_angle_between_vectors_rejects_vectors_of_different_length <- function() {
  # the binding's own shape check, before the library is called
  assert_error(clock_hand_angle_between_vectors(c(1, 0), c(0, 1, 0), c(0, 1)), "v1 has two dimensions, v2 three")
}

test_clock_hand_angle_between_vectors_rejects_a_non_finite_reference <- function() {
  for (value in NON_FINITE_VALUES) {
    assert_error(clock_hand_angle_between_vectors(c(1, 0), c(0, 1), c(value, 1)),
                 paste("the reference holds", value), ERR_NAN_INF)
  }
}

test_clock_hand_angle_between_vectors_rejects_empty_vectors <- function() {
  assert_error(clock_hand_angle_between_vectors(numeric(0), numeric(0), numeric(0)), "no dimensions", ERR_EMPTY_INPUT)
}

test_clock_hand_angle_between_vectors_rejects_a_reference_that_orients_nothing <- function() {
  # the documented ERR_INVALID_INPUT: the turn is in the (1, 2) plane, and axis 3 cannot sign it
  assert_error(clock_hand_angle_between_vectors(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)),
               "a reference orthogonal to the rotation", ERR_INVALID_INPUT)
}

# ---- clock_hand_angles_for_shift_vectors -----------------------------------------------------

test_clock_hand_angles_for_shift_vectors <- function() {
  .assert_vector(clock_hand_angles_for_shift_vectors(.fields(), VECS_SELECTION_MASK, seq_len(N_AXES)),
                 N_SELECTED_VECS, "clock_hand_angles_for_shift_vectors")
}

test_clock_hand_angles_for_shift_vectors_turns_from_origin_to_origin_plus_shift <- function() {
  # pass-through only, the rule of the single-pair test above: each field is [origin, shift] and
  # turns from o to o + s. The selected second field turns from e1 to e1 + (-1, 1) = e2, +pi/2;
  # reading the shift as the target would give 3pi/4, and the unselected first field turns to
  # e1 + (-1, -1) = -e2, -pi/2
  fields <- array(0, dim = c(2L, 2L, 2L))
  fields[, 1, ] <- c(1, 0)
  fields[, 2, 1] <- c(-1, -1)
  fields[, 2, 2] <- c(-1, 1)
  result <- clock_hand_angles_for_shift_vectors(fields, c(FALSE, TRUE), c(-1, 1))
  assert_equal_numeric(result, pi / 2, msg = paste("expected pi/2, got", toString(result)))
}

test_clock_hand_angles_for_shift_vectors_rejects_a_mask_of_the_wrong_length <- function() {
  # the binding's own shape check, before the library is called
  assert_error(clock_hand_angles_for_shift_vectors(.fields(), VECS_SELECTION_MASK[-1], rep(1, N_AXES)),
               "fields has five fields, the mask four")
}

# `fields` is (n_dims, 2, n_fields); a middle extent of 1 used to be read past (#221)
test_clock_hand_angles_for_shift_vectors_rejects_fields_without_two_points <- function() {
  assert_error(clock_hand_angles_for_shift_vectors(.fields()[, 1, , drop = FALSE], VECS_SELECTION_MASK, rep(1, N_AXES)),
               "a field without its shift must be refused")
}

test_clock_hand_angles_for_shift_vectors_rejects_non_finite_fields <- function() {
  for (value in NON_FINITE_VALUES) {
    assert_error(clock_hand_angles_for_shift_vectors(.with(.fields(), value), VECS_SELECTION_MASK, rep(1, N_AXES)),
                 paste("fields holds", value), ERR_NAN_INF)
  }
}

test_clock_hand_angles_for_shift_vectors_rejects_an_empty_selection <- function() {
  assert_error(clock_hand_angles_for_shift_vectors(.fields(), rep(FALSE, N_VECS), rep(1, N_AXES)),
               "no field selected", ERR_EMPTY_INPUT)
}

test_clock_hand_angles_for_shift_vectors_rejects_a_reference_that_orients_nothing <- function() {
  # the documented ERR_INVALID_INPUT, as for a single pair
  fields <- array(0, dim = c(3L, 2L, 1L))
  fields[1, 1, 1] <- 1
  fields[2, 2, 1] <- 1
  assert_error(clock_hand_angles_for_shift_vectors(fields, TRUE, c(0, 0, 1)),
               "a reference orthogonal to the rotation", ERR_INVALID_INPUT)
}

# ---- the three contribution procedures -------------------------------------------------------

test_contributions <- function() {
  for (name in names(CONTRIBUTION_PROCEDURES)) {
    .assert_vector(CONTRIBUTION_PROCEDURES[[name]](seq_len(N_AXES)), N_AXES, name)
  }
}

test_contributions_keep_the_axis_order <- function() {
  # pass-through only: |1| and |-3| are a quarter and three quarters of 4; reversed, the other way
  for (name in names(CONTRIBUTION_PROCEDURES)) {
    result <- CONTRIBUTION_PROCEDURES[[name]](c(1, -3))
    assert_equal_numeric(result, c(0.25, 0.75), msg = paste0(name, ": got ", toString(result)))
  }
}

test_contributions_reject_a_zero_vector <- function() {
  for (name in names(CONTRIBUTION_PROCEDURES)) {
    assert_error(CONTRIBUTION_PROCEDURES[[name]](c(0, 0, 0)), paste0(name, ": nothing to apportion"),
                 ERR_DIVISION_BY_ZERO)
  }
}

test_contributions_reject_non_finite_values <- function() {
  for (name in names(CONTRIBUTION_PROCEDURES)) {
    for (value in NON_FINITE_VALUES) {
      assert_error(CONTRIBUTION_PROCEDURES[[name]](c(1, value)), paste0(name, ": vec holds ", value), ERR_NAN_INF)
    }
  }
}

test_contributions_reject_an_empty_vector <- function() {
  for (name in names(CONTRIBUTION_PROCEDURES)) {
    assert_error(CONTRIBUTION_PROCEDURES[[name]](numeric(0)), paste0(name, ": no axes"), ERR_EMPTY_INPUT)
  }
}

run_all_tests()
