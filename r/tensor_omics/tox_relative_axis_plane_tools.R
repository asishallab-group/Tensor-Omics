# Generated. Do not edit.

#' Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
#'
#' Generated from the Fortran module \code{tox_relative_axis_plane_tools}.
#'
#' @param vecs a numeric matrix. matrix with expression vectors
#' @param vecs_selection_mask a logical vector. `TRUE` for vectors where projection is to be computed
#' @param axes_selection_mask a logical vector. `TRUE` for axes to be included in RAP
#' @return a numeric matrix. projected vectors
#' @export
omics_vector_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
    vecs <- .tox_as_double_matrix(vecs, "vecs")
    vecs_selection_mask <- .tox_as_logical(vecs_selection_mask, "vecs_selection_mask")
    axes_selection_mask <- .tox_as_logical(axes_selection_mask, "axes_selection_mask")
    if (length(axes_selection_mask) != dim(vecs)[1])
        .tox_shape_error("axes_selection_mask", length(axes_selection_mask), "vecs", dim(vecs)[1])
    if (length(vecs_selection_mask) != dim(vecs)[2])
        .tox_shape_error("vecs_selection_mask", length(vecs_selection_mask), "vecs", dim(vecs)[2])

    .result <- .Call("omics_vector_RAP_projection_call", vecs, vecs_selection_mask, axes_selection_mask)
    .arguments <- c("vecs", "n_axes", "n_vecs", "vecs_selection_mask", "n_selected_vecs", "axes_selection_mask", "n_selected_axes", "projections", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$projections
}

#' Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
#'
#' Generated from the Fortran module \code{tox_relative_axis_plane_tools}.
#'
#' @param fields a numeric array of rank 3. matrix with vector fields; each field holds two vectors, the origin first and the target second
#' @param fields_selection_mask a logical vector. `TRUE` for vectors where projection is to be computed
#' @param axes_selection_mask a logical vector. `TRUE` for axes to be included in RAP
#' @return a numeric matrix. projected vectors
#' @export
omics_field_RAP_projection <- function(fields, fields_selection_mask, axes_selection_mask) {
    fields <- .tox_as_double_array(fields, "fields", 3L)
    fields_selection_mask <- .tox_as_logical(fields_selection_mask, "fields_selection_mask")
    axes_selection_mask <- .tox_as_logical(axes_selection_mask, "axes_selection_mask")
    if (length(axes_selection_mask) != dim(fields)[1])
        .tox_shape_error("axes_selection_mask", length(axes_selection_mask), "fields", dim(fields)[1])
    if (length(fields_selection_mask) != dim(fields)[3])
        .tox_shape_error("fields_selection_mask", length(fields_selection_mask), "fields", dim(fields)[3])

    .result <- .Call("omics_field_RAP_projection_call", fields, fields_selection_mask, axes_selection_mask)
    .arguments <- c("fields", "n_axes", "n_fields", "fields_selection_mask", "n_selected_fields", "axes_selection_mask", "n_selected_axes", "projections", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$projections
}

#' Compute the signed clock hand angle between two RAP-projected and normalized vectors.
#'
#' Calculates the signed rotation angle between two normalized vectors in RAP space.
#' For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
#'
#' Generated from the Fortran module \code{tox_relative_axis_plane_tools}.
#'
#' @param v1 a numeric vector. First normalized vector in RAP space
#' @param v2 a numeric vector. Second normalized vector in RAP space
#' @param selected_axes_for_signed a integer vector. Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
#' @return a numeric scalar. Signed angle between vectors in radians [-π, π]
#' @export
clock_hand_angle_between_vectors <- function(v1, v2, selected_axes_for_signed) {
    v1 <- .tox_as_double_vector(v1, "v1")
    v2 <- .tox_as_double_vector(v2, "v2")
    selected_axes_for_signed <- .tox_as_integer_vector(selected_axes_for_signed, "selected_axes_for_signed")
    if (length(v2) != length(v1))
        .tox_shape_error("v2", length(v2), "v1", length(v1))

    .result <- .Call("clock_hand_angle_between_vectors_call", v1, v2, selected_axes_for_signed)
    .arguments <- c("v1", "v2", "n_dims", "signed_angle", "selected_axes_for_signed", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$signed_angle
}

#' Compute signed rotation angles between for shift vectors, so between their origin and target
#'
#' Generated from the Fortran module \code{tox_relative_axis_plane_tools}.
#'
#' @param fields a numeric array of rank 3. matrix with vector fields; each field holds two vectors, the origin first and the target second
#' @param fields_selection_mask a logical vector. TRUE for vector pairs where angle should be computed
#' @param selected_axes_for_signed a integer vector. Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
#' @return a numeric vector. Signed rotation angles between vector pairs in radians [-π, π]
#' @export
clock_hand_angles_for_shift_vectors <- function(fields, fields_selection_mask, selected_axes_for_signed) {
    fields <- .tox_as_double_array(fields, "fields", 3L)
    fields_selection_mask <- .tox_as_logical(fields_selection_mask, "fields_selection_mask")
    selected_axes_for_signed <- .tox_as_integer_vector(selected_axes_for_signed, "selected_axes_for_signed")
    if (length(fields_selection_mask) != dim(fields)[3])
        .tox_shape_error("fields_selection_mask", length(fields_selection_mask), "fields", dim(fields)[3])

    .result <- .Call("clock_hand_angles_for_shift_vectors_call", fields, fields_selection_mask, selected_axes_for_signed)
    .arguments <- c("fields", "n_dims", "n_fields", "fields_selection_mask", "n_selected_fields", "selected_axes_for_signed", "signed_angles", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$signed_angles
}

#' Compute the fractional contribution of each axis to a RAP-projected and normalized vector
#'
#' Shared utility: the shift-vector and expression-vector entry points below both drive it.
#'
#' Generated from the Fortran module \code{tox_relative_axis_plane_tools}.
#'
#' @param vec a numeric vector. RAP-projected and normalized vector (expression or shift)
#' @return a numeric vector. Fractional contribution of each axis (output), values in [0,1], sum to 1
#' @export
compute_relative_axis_contributions <- function(vec) {
    vec <- .tox_as_double_vector(vec, "vec")
    .result <- .Call("compute_relative_axis_contributions_call", vec)
    .arguments <- c("vec", "n_axes", "contributions", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$contributions
}

#' Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
#'
#' Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
#'
#' Generated from the Fortran module \code{tox_relative_axis_plane_tools}.
#'
#' @param vec a numeric vector. RAP-projected and normalized shift vector
#' @return a numeric vector. Fractional contribution of each axis (output), values in [0,1], sum to 1
#' @export
relative_axes_changes_from_shift_vector <- function(vec) {
    vec <- .tox_as_double_vector(vec, "vec")
    .result <- .Call("relative_axes_changes_from_shift_vector_call", vec)
    .arguments <- c("vec", "n_axes", "contributions", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$contributions
}

#' Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
#'
#' Wrapper for single RAP-projected expression vectors
#'
#' Generated from the Fortran module \code{tox_relative_axis_plane_tools}.
#'
#' @param vec a numeric vector. RAP-projected and normalized expression vector
#' @return a numeric vector. Fractional contribution of each axis (output), values in [0,1], sum to 1
#' @export
relative_axes_expression_from_expression_vector <- function(vec) {
    vec <- .tox_as_double_vector(vec, "vec")
    .result <- .Call("relative_axes_expression_from_expression_vector_call", vec)
    .arguments <- c("vec", "n_axes", "contributions", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$contributions
}
