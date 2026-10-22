# Generated. Do not edit.

#' Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
#'
#' Generated from the Fortran procedure \code{tox_relative_axis_plane_tools::omics_vector_RAP_projection}, whose argument names
#' are the ones an error message reports.
#'
#' @param vecs a numeric matrix. matrix with expression vectors
#' @param vecs_selection_mask a logical vector. `TRUE` for vectors where projection is to be computed
#' @param axes_selection_mask a logical vector. `TRUE` for axes to be included in RAP
#' @return a numeric matrix. projected vectors
#' @export
omics_vector_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
    vecs <- .tox_as_double_matrix(vecs, "vecs")
    vecs_selection_mask <- .tox_as_logical_vector(vecs_selection_mask, "vecs_selection_mask")
    axes_selection_mask <- .tox_as_logical_vector(axes_selection_mask, "axes_selection_mask")
    if (length(axes_selection_mask) != dim(vecs)[1])
        .tox_shape_error("axes_selection_mask", length(axes_selection_mask), "vecs", dim(vecs)[1])
    if (length(vecs_selection_mask) != dim(vecs)[2])
        .tox_shape_error("vecs_selection_mask", length(vecs_selection_mask), "vecs", dim(vecs)[2])

    .result <- .Call("omics_vector_RAP_projection_call", vecs, vecs_selection_mask, axes_selection_mask)
    .arguments <- c("vecs", "n_axes", "n_vecs", "vecs_selection_mask", "n_selected_vecs", "axes_selection_mask", "n_selected_axes", "projections", "ierr")
    .sources <- c(NA_character_, "vecs", "vecs", NA_character_, "projections", NA_character_, "projections", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$projections
}

#' Project the shifts of selected vector fields onto the RAP constructed from a selected set of axes.
#'
#' Generated from the Fortran procedure \code{tox_relative_axis_plane_tools::omics_field_RAP_projection}, whose argument names
#' are the ones an error message reports.
#'
#' @param fields a numeric array of rank 3. matrix with vector fields; each field holds two vectors, its origin first (e.g. a
#'   family centroid) and the shift from it second (e.g. paralog minus centroid), as
#'   \code{\link{compute_shift_vector_field}} stores
#'   them. The shift is projected; the origin does not enter the projection.
#' @param fields_selection_mask a logical vector. `TRUE` for fields where projection is to be computed
#' @param axes_selection_mask a logical vector. `TRUE` for axes to be included in RAP
#' @return a numeric matrix. projected vectors
#' @export
omics_field_RAP_projection <- function(fields, fields_selection_mask, axes_selection_mask) {
    fields <- .tox_as_double_array(fields, "fields", 3L)
    fields_selection_mask <- .tox_as_logical_vector(fields_selection_mask, "fields_selection_mask")
    axes_selection_mask <- .tox_as_logical_vector(axes_selection_mask, "axes_selection_mask")
    if (length(axes_selection_mask) != dim(fields)[1])
        .tox_shape_error("axes_selection_mask", length(axes_selection_mask), "fields", dim(fields)[1])
    if (length(fields_selection_mask) != dim(fields)[3])
        .tox_shape_error("fields_selection_mask", length(fields_selection_mask), "fields", dim(fields)[3])

    .result <- .Call("omics_field_RAP_projection_call", fields, fields_selection_mask, axes_selection_mask)
    .arguments <- c("fields", "n_axes", "n_fields", "fields_selection_mask", "n_selected_fields", "axes_selection_mask", "n_selected_axes", "projections", "ierr")
    .sources <- c(NA_character_, "fields", "fields", NA_character_, "projections", NA_character_, "projections", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$projections
}

#' Compute the signed clock hand angle between two RAP-projected vectors.
#'
#' The unsigned angle is the one between the two directions; their lengths do not matter, so
#' the vectors need not be normalized. `orientation_reference` supplies the sign by saying
#' which way round the plane the two vectors span counts as positive. Reports
#' `ERR_DIVISION_BY_ZERO` when `v1` or `v2` is the zero vector, which has no direction, and
#' `ERR_INVALID_INPUT` when the reference is orthogonal to the rotation and so orients nothing.
#'
#' Generated from the Fortran procedure \code{tox_relative_axis_plane_tools::clock_hand_angle_between_vectors}, whose argument names
#' are the ones an error message reports.
#'
#' @param v1 a numeric vector. First vector in RAP space, of any length but zero
#' @param v2 a numeric vector. Second vector in RAP space, of any length but zero
#' @param orientation_reference a numeric vector. Orients the plane the rotation happens in, so the angle can carry a sign. A
#'   rotation from one vector to another has no inherent direction above two
#'   dimensions -- and in RAP space not even in two, since the axes are tissues or
#'   factors and carry no handedness -- so the caller states which way round counts
#'   as positive. The sign is that of this vector's component along the rotation.
#'   For three selected tissues, `d x v1`, with `d` the space diagonal, reproduces the
#'   determinant rule `sign(det[d, v1, v2])`; a fixed vector, such as the anchor axis
#'   projected onto the RAP, gives every angle the same sense of clockwise.
#' @return a numeric scalar. Signed angle between vectors in radians [-pi, pi]
#' @export
clock_hand_angle_between_vectors <- function(v1, v2, orientation_reference) {
    v1 <- .tox_as_double_vector(v1, "v1")
    v2 <- .tox_as_double_vector(v2, "v2")
    orientation_reference <- .tox_as_double_vector(orientation_reference, "orientation_reference")
    if (length(v2) != length(v1))
        .tox_shape_error("v2", length(v2), "v1", length(v1))
    if (length(orientation_reference) != length(v1))
        .tox_shape_error("orientation_reference", length(orientation_reference), "v1", length(v1))

    .result <- .Call("clock_hand_angle_between_vectors_call", v1, v2, orientation_reference)
    .arguments <- c("v1", "v2", "n_dims", "orientation_reference", "signed_angle", "ierr")
    .sources <- c(NA_character_, NA_character_, "v1", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$signed_angle
}

#' Compute the signed clock hand angle of every selected shift, from its origin to the point it reaches
#'
#' Each selected field, an origin `o` and a shift `s` from it, turns from `o` to `o + s` -- from
#' a family centroid to its paralog, for the fields
#' \code{\link{compute_shift_vector_field}} stores -- by the
#' rule of
#' \code{\link{clock_hand_angle_between_vectors}},
#' with one `orientation_reference` shared by the whole batch. The rule angles RAP-space
#' vectors, so project origins and shifts first; projection is linear, so `o + s` of the
#' projected pair is the projected paralog. A single selected field whose origin or `o + s`
#' is zero fails the call with `ERR_DIVISION_BY_ZERO`, and one whose rotation the reference
#' fails to orient with `ERR_INVALID_INPUT`.
#'
#' Generated from the Fortran procedure \code{tox_relative_axis_plane_tools::clock_hand_angles_for_shift_vectors}, whose argument names
#' are the ones an error message reports.
#'
#' @param fields a numeric array of rank 3. matrix with vector fields; each field holds two vectors, its origin first and the
#'   shift from it second, as
#'   \code{\link{compute_shift_vector_field}} stores them
#' @param fields_selection_mask a logical vector. TRUE for vector pairs where angle should be computed
#' @param orientation_reference a numeric vector. Orients the plane the rotation happens in, so the angle can carry a sign. A
#'   rotation from one vector to another has no inherent direction above two
#'   dimensions -- and in RAP space not even in two, since the axes are tissues or
#'   factors and carry no handedness -- so the caller states which way round counts
#'   as positive. The sign is that of this vector's component along the rotation.
#'   For three selected tissues, `d x v1`, with `d` the space diagonal, reproduces the
#'   determinant rule `sign(det[d, v1, v2])`; a fixed vector, such as the anchor axis
#'   projected onto the RAP, gives every angle the same sense of clockwise.
#' @return a numeric vector. Signed rotation angles between vector pairs in radians [-pi, pi]
#' @export
clock_hand_angles_for_shift_vectors <- function(fields, fields_selection_mask, orientation_reference) {
    fields <- .tox_as_double_array(fields, "fields", 3L)
    fields_selection_mask <- .tox_as_logical_vector(fields_selection_mask, "fields_selection_mask")
    orientation_reference <- .tox_as_double_vector(orientation_reference, "orientation_reference")
    if (length(orientation_reference) != dim(fields)[1])
        .tox_shape_error("orientation_reference", length(orientation_reference), "fields", dim(fields)[1])
    if (length(fields_selection_mask) != dim(fields)[3])
        .tox_shape_error("fields_selection_mask", length(fields_selection_mask), "fields", dim(fields)[3])

    .result <- .Call("clock_hand_angles_for_shift_vectors_call", fields, fields_selection_mask, orientation_reference)
    .arguments <- c("fields", "n_dims", "n_fields", "fields_selection_mask", "n_selected_fields", "orientation_reference", "signed_angles", "ierr")
    .sources <- c(NA_character_, "fields", "fields", NA_character_, "signed_angles", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$signed_angles
}

#' Compute the fractional contribution of each axis to a RAP-projected vector
#'
#' Shared utility: the shift-vector and expression-vector entry points below both drive it.
#' The shares depend on the vector's direction alone, so it need not be normalized; the zero
#' vector, which has no direction, is `ERR_DIVISION_BY_ZERO`.
#'
#' Generated from the Fortran procedure \code{tox_relative_axis_plane_tools::compute_relative_axis_contributions}, whose argument names
#' are the ones an error message reports.
#'
#' @param vec a numeric vector. RAP-projected vector (expression or shift), of any length but zero
#' @return a numeric vector. Fractional contribution of each axis, values in [0,1], sum to 1
#' @export
compute_relative_axis_contributions <- function(vec) {
    vec <- .tox_as_double_vector(vec, "vec")
    .result <- .Call("compute_relative_axis_contributions_call", vec)
    .arguments <- c("vec", "n_axes", "contributions", "ierr")
    .sources <- c(NA_character_, "vec", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$contributions
}

#' Compute fractional contribution of each axis to a RAP-projected shift vector.
#'
#' Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
#'
#' Generated from the Fortran procedure \code{tox_relative_axis_plane_tools::relative_axes_changes_from_shift_vector}, whose argument names
#' are the ones an error message reports.
#'
#' @param vec a numeric vector. RAP-projected shift vector, of any length but zero
#' @return a numeric vector. Fractional contribution of each axis, values in [0,1], sum to 1
#' @export
relative_axes_changes_from_shift_vector <- function(vec) {
    vec <- .tox_as_double_vector(vec, "vec")
    .result <- .Call("relative_axes_changes_from_shift_vector_call", vec)
    .arguments <- c("vec", "n_axes", "contributions", "ierr")
    .sources <- c(NA_character_, "vec", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$contributions
}

#' Compute fractional contribution of each axis to a RAP-projected expression vector.
#'
#' Wrapper for single RAP-projected expression vectors
#'
#' Generated from the Fortran procedure \code{tox_relative_axis_plane_tools::relative_axes_expression_from_expression_vector}, whose argument names
#' are the ones an error message reports.
#'
#' @param vec a numeric vector. RAP-projected expression vector, of any length but zero
#' @return a numeric vector. Fractional contribution of each axis, values in [0,1], sum to 1
#' @export
relative_axes_expression_from_expression_vector <- function(vec) {
    vec <- .tox_as_double_vector(vec, "vec")
    .result <- .Call("relative_axes_expression_from_expression_vector_call", vec)
    .arguments <- c("vec", "n_axes", "contributions", "ierr")
    .sources <- c(NA_character_, "vec", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$contributions
}
