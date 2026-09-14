# Generated. Do not edit.

#' Build a k-d tree index using a stack-based, non-recursive approach
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::build_kd_index}, whose argument names
#' are the ones an error message reports.
#'
#' @param points a numeric matrix. Data points
#' @param dimension_order a integer vector. Dimension order (by variance)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @return a integer vector. Output index array (k-d tree order)
#' @export
build_kd_index <- function(points, dimension_order) {
    points <- .tox_as_double_matrix(points, "points")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    if (length(dimension_order) != dim(points)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "points", dim(points)[1])

    .result <- .Call("build_kd_index_call", points, dimension_order)
    .arguments <- c("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr")
    .sources <- c(NA_character_, "points", "points", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$kd_indices
}

#' Build a k-d tree index over points on the unit sphere (unit vectors)
#'
#' This is a thin, semantically-named wrapper: partitioning is identical to
#' \code{\link{build_kd_index}} (plain per-axis median splits);
#' callers are responsible for ensuring `points` are actually unit-normalized beforehand.
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::build_spherical_kd}, whose argument names
#' are the ones an error message reports.
#'
#' @param points a numeric matrix. Data points
#' @param dimension_order a integer vector. Dimension order (by variance)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @return a integer vector. Output index array (k-d tree order)
#' @export
build_spherical_kd <- function(points, dimension_order) {
    points <- .tox_as_double_matrix(points, "points")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    if (length(dimension_order) != dim(points)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "points", dim(points)[1])

    .result <- .Call("build_spherical_kd_call", points, dimension_order)
    .arguments <- c("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr")
    .sources <- c(NA_character_, "points", "points", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$kd_indices
}
