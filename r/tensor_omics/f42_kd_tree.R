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

#' Find reference points within a radius around a query point.
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::vicinity_vectors}, whose argument names
#' are the ones an error message reports.
#'
#' @param query_point a numeric vector. Coordinate vector used as the center of the search
#' @param points a numeric matrix. Ambient point matrix
#' @param r a numeric scalar. Search radius
#'   The minimum valid value is `0.0`.
#' @param dimension_order a integer vector. Sequence of k-d tree split dimensions
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param kd_indices a integer vector. K-d tree index sequence
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_points`.
#' @return a logical vector. Mask indicating points within the search radius
#' @export
vicinity_vectors <- function(query_point, points, r, dimension_order, kd_indices) {
    query_point <- .tox_as_double_vector(query_point, "query_point")
    points <- .tox_as_double_matrix(points, "points")
    r <- .tox_as_double_scalar(r, "r")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    if (dim(points)[1] != length(query_point))
        .tox_shape_error("points", dim(points)[1], "query_point", length(query_point))
    if (length(dimension_order) != length(query_point))
        .tox_shape_error("dimension_order", length(dimension_order), "query_point", length(query_point))
    if (length(kd_indices) != dim(points)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "points", dim(points)[2])

    .result <- .Call("vicinity_vectors_call", query_point, points, r, dimension_order, kd_indices)
    .arguments <- c("query_point", "points", "n_dimensions", "n_points", "r", "dimension_order", "kd_indices", "vicinity_mask", "ierr")
    .sources <- c(NA_character_, NA_character_, "query_point", "points", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$vicinity_mask
}

#' Count reference points within a radius around a query point.
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::vicinity_vectors_count}, whose argument names
#' are the ones an error message reports.
#'
#' @param query_point a numeric vector. Coordinate vector used as the center of the search
#' @param points a numeric matrix. Ambient point matrix
#' @param r a numeric scalar. Search radius
#'   The minimum valid value is `0.0`.
#' @param dimension_order a integer vector. Sequence of k-d tree split dimensions
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param kd_indices a integer vector. K-d tree index sequence
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_points`.
#' @return a integer scalar. Number of points within the search radius
#' @export
vicinity_vectors_count <- function(query_point, points, r, dimension_order, kd_indices) {
    query_point <- .tox_as_double_vector(query_point, "query_point")
    points <- .tox_as_double_matrix(points, "points")
    r <- .tox_as_double_scalar(r, "r")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    if (dim(points)[1] != length(query_point))
        .tox_shape_error("points", dim(points)[1], "query_point", length(query_point))
    if (length(dimension_order) != length(query_point))
        .tox_shape_error("dimension_order", length(dimension_order), "query_point", length(query_point))
    if (length(kd_indices) != dim(points)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "points", dim(points)[2])

    .result <- .Call("vicinity_vectors_count_call", query_point, points, r, dimension_order, kd_indices)
    .arguments <- c("query_point", "points", "n_dimensions", "n_points", "r", "dimension_order", "kd_indices", "n_neighbors", "ierr")
    .sources <- c(NA_character_, NA_character_, "query_point", "points", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$n_neighbors
}
