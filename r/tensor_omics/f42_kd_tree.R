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

#' Find the k nearest neighbors of a query point in a pre-built k-d tree
#'
#' Via a bounded max-heap kept directly in `neighbors`/`distances` and splitting-plane
#' pruning. Requires `k_neighbors <= n_points`: every point is then guaranteed visited
#' before the heap can still have room, so no fallback for "fewer than k found" is needed.
#' Does not guarantee `neighbors`/`distances` are sorted nearest-to-farthest (max-heap
#' order internally).
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::kd_knn_query}, whose argument names
#' are the ones an error message reports.
#'
#' @param points a numeric matrix. Original points dataset
#' @param kd_indices a integer vector. Pre-built k-d tree index, see \code{\link{build_kd_index}}
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_points`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param query_point a numeric vector. Query point coordinates
#' @param k_neighbors a integer scalar. Number of neighbors to find
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_points`.
#' @return a named list with elements:
#'   \item{neighbors}{a integer vector. Output: indices of the k nearest neighbors (nearest-to-farthest order not guaranteed, max-heap order internally)}
#'   \item{distances}{a numeric vector. Output: Euclidean distances to the k nearest neighbors}
#' @export
kd_knn_query <- function(points, kd_indices, dimension_order, query_point, k_neighbors) {
    points <- .tox_as_double_matrix(points, "points")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    query_point <- .tox_as_double_vector(query_point, "query_point")
    k_neighbors <- .tox_as_integer_scalar(k_neighbors, "k_neighbors")
    if (length(dimension_order) != dim(points)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "points", dim(points)[1])
    if (length(query_point) != dim(points)[1])
        .tox_shape_error("query_point", length(query_point), "points", dim(points)[1])
    if (length(kd_indices) != dim(points)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "points", dim(points)[2])

    .result <- .Call("kd_knn_query_call", points, kd_indices, dimension_order, query_point, k_neighbors)
    .arguments <- c("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "k_neighbors", "neighbors", "distances", "ierr")
    .sources <- c(NA_character_, "points", "points", NA_character_, NA_character_, NA_character_, "neighbors", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighbors = .result$neighbors,
        distances = .result$distances
    )
}

#' Mark every point within `radius` of a query point in a pre-built k-d tree
#'
#' Same iterative, stack-based traversal and splitting-plane pruning as
#' \code{\link{kd_knn_query}} (the near side of each split is
#' always explored, the far side only when it is still within `radius` of the splitting
#' plane), with a fixed radius bound instead of a k-nearest-neighbor heap. Compares
#' squared distances against a precomputed `radius**2` (no `sqrt` per node visited).
#'
#' Fits a caller that already does an O(n_points) pass over the result (e.g. merging it
#' into an existing coverage mask via `.or.`). A caller issuing many independent range
#' queries per outer step (e.g. one per candidate point in a greedy loop) should use
#' \code{\link{kd_range_query_list}} or
#' \code{\link{kd_range_query_count}} instead, since the
#' `in_radius_mask = FALSE` reset here costs O(n_points) on every call regardless of how
#' few points are actually found.
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::kd_range_query_mask}, whose argument names
#' are the ones an error message reports.
#'
#' @param points a numeric matrix. Original points dataset
#' @param kd_indices a integer vector. Pre-built k-d tree index
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_points`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param query_point a numeric vector. Query point coordinates
#' @param radius a numeric scalar. Search radius
#'   The minimum valid value is `0.0`.
#' @return a logical vector. Output: TRUE for points within `radius` of `query_point`
#' @export
kd_range_query_mask <- function(points, kd_indices, dimension_order, query_point, radius) {
    points <- .tox_as_double_matrix(points, "points")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    query_point <- .tox_as_double_vector(query_point, "query_point")
    radius <- .tox_as_double_scalar(radius, "radius")
    if (length(dimension_order) != dim(points)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "points", dim(points)[1])
    if (length(query_point) != dim(points)[1])
        .tox_shape_error("query_point", length(query_point), "points", dim(points)[1])
    if (length(kd_indices) != dim(points)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "points", dim(points)[2])

    .result <- .Call("kd_range_query_mask_call", points, kd_indices, dimension_order, query_point, radius)
    .arguments <- c("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "in_radius_mask", "ierr")
    .sources <- c(NA_character_, "points", "points", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$in_radius_mask
}

#' List every point within `radius` of a query point in a pre-built k-d tree
#'
#' Same traversal and pruning as \code{\link{kd_range_query_mask}},
#' but writes matches into a caller-provided compact index buffer (`neighbors(1:n_found)`)
#' instead of a full-size logical mask, so repeated calls -- e.g. once per candidate point
#' in an outer greedy loop -- don't each pay an O(n_points) reset. A caller that only
#' needs the count, not the identities, should use
#' \code{\link{kd_range_query_count}} instead, to skip the
#' index-buffer writes entirely.
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::kd_range_query_list}, whose argument names
#' are the ones an error message reports.
#'
#' @param points a numeric matrix. Original points dataset
#' @param kd_indices a integer vector. Pre-built k-d tree index
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_points`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param query_point a numeric vector. Query point coordinates
#' @param radius a numeric scalar. Search radius
#'   The minimum valid value is `0.0`.
#' @return a named list with elements:
#'   \item{neighbors}{a integer vector. Output: indices within `radius`, valid in `neighbors(1:n_found)`}
#'   \item{n_found}{a integer scalar. Output: number of points within `radius`}
#' @export
kd_range_query_list <- function(points, kd_indices, dimension_order, query_point, radius) {
    points <- .tox_as_double_matrix(points, "points")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    query_point <- .tox_as_double_vector(query_point, "query_point")
    radius <- .tox_as_double_scalar(radius, "radius")
    if (length(dimension_order) != dim(points)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "points", dim(points)[1])
    if (length(query_point) != dim(points)[1])
        .tox_shape_error("query_point", length(query_point), "points", dim(points)[1])
    if (length(kd_indices) != dim(points)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "points", dim(points)[2])

    .result <- .Call("kd_range_query_list_call", points, kd_indices, dimension_order, query_point, radius)
    .arguments <- c("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "neighbors", "n_found", "ierr")
    .sources <- c(NA_character_, "points", "points", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        neighbors = .result$neighbors,
        n_found = .result$n_found
    )
}

#' Count the points within `radius` of a query point in a pre-built k-d tree
#'
#' Same traversal and pruning as
#' \code{\link{kd_range_query_mask}}, but writes no index
#' buffer at all -- only a scalar count. The right choice when the identities of the
#' points found are never needed, e.g. a per-point local-density label.
#'
#' Generated from the Fortran procedure \code{f42_kd_tree::kd_range_query_count}, whose argument names
#' are the ones an error message reports.
#'
#' @param points a numeric matrix. Original points dataset
#' @param kd_indices a integer vector. Pre-built k-d tree index
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_points`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param query_point a numeric vector. Query point coordinates
#' @param radius a numeric scalar. Search radius
#'   The minimum valid value is `0.0`.
#' @return a integer scalar. Output: number of points within `radius`
#' @export
kd_range_query_count <- function(points, kd_indices, dimension_order, query_point, radius) {
    points <- .tox_as_double_matrix(points, "points")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    query_point <- .tox_as_double_vector(query_point, "query_point")
    radius <- .tox_as_double_scalar(radius, "radius")
    if (length(dimension_order) != dim(points)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "points", dim(points)[1])
    if (length(query_point) != dim(points)[1])
        .tox_shape_error("query_point", length(query_point), "points", dim(points)[1])
    if (length(kd_indices) != dim(points)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "points", dim(points)[2])

    .result <- .Call("kd_range_query_count_call", points, kd_indices, dimension_order, query_point, radius)
    .arguments <- c("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "query_point", "radius", "neighbor_count", "ierr")
    .sources <- c(NA_character_, "points", "points", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$neighbor_count
}
