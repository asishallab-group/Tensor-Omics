# Generated. Do not edit.

#' Local-density search radius, a percentile of mean-to-vector distances
#'
#' Computes the mean vector of `vectors`, the Euclidean distance from every vector to it,
#' and returns the `mean_to_other_vecs_dist_quant` percentile of those distances -- the
#' radius later used by `density_labels` to measure local density around each vector.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_seeding::calculate_density_radius_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param mean_to_other_vecs_dist_quant a numeric scalar. Percentile (0.0 to 1.0) of mean-to-vector distances used as the density radius
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.15`.
#' @return a numeric scalar. Resulting density search radius
#' @export
calculate_density_radius <- function(vectors, mean_to_other_vecs_dist_quant = 0.15) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    mean_to_other_vecs_dist_quant <- .tox_as_double_scalar(mean_to_other_vecs_dist_quant, "mean_to_other_vecs_dist_quant")
    .result <- .Call("calculate_density_radius_call", vectors, mean_to_other_vecs_dist_quant)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "mean_to_other_vecs_dist_quant", "radius", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$radius
}

#' Per-vector density label, the count of vectors (including itself) within `radius`
#'
#' $\rho_i = \sum_j \mathbf{1}(d(v_i, v_j) \le radius)$, via a k-d tree range-count query
#' per vector -- see \code{\link{kd_range_query_count}}.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_seeding::density_labels_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param radius a numeric scalar. Density search radius, see `calculate_density_radius`
#'   The minimum valid value is `0.0`.
#' @return a numeric vector. Per-vector density label
#' @export
density_labels <- function(vectors, kd_indices, dimension_order, radius) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    radius <- .tox_as_double_scalar(radius, "radius")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])

    .result <- .Call("density_labels_call", vectors, kd_indices, dimension_order, radius)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "radius", "labels", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$labels
}

#' Select seed points via greedy, density-ranked, coverage-based selection
#'
#' Ranks vectors by density label, descending (see `density_labels`). Starting with the
#' highest-density unvisited vector, marks it a seed, marks every vector within the
#' density radius of it as visited, and continues with the next-highest-density
#' unvisited vector until none remain -- so only genuinely uncovered regions can seed
#' another ensemble.
#'
#' Generated from the Fortran procedure \code{tox_shape_truthful_clustering_seeding::seeds_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param vectors a numeric matrix. Input data matrix
#' @param kd_indices a integer vector. Pre-built k-d tree index over `vectors`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_vectors`.
#' @param dimension_order a integer vector. Dimension order used to build `kd_indices`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_dimensions`.
#' @param mean_to_other_vecs_dist_quant a numeric scalar. Density-radius percentile (0.0 to 1.0), see `calculate_density_radius`
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.15`.
#' @return a logical vector. TRUE for points selected as seeds
#' @export
seeds <- function(vectors, kd_indices, dimension_order, mean_to_other_vecs_dist_quant = 0.15) {
    vectors <- .tox_as_double_matrix(vectors, "vectors")
    kd_indices <- .tox_as_integer_vector(kd_indices, "kd_indices")
    dimension_order <- .tox_as_integer_vector(dimension_order, "dimension_order")
    mean_to_other_vecs_dist_quant <- .tox_as_double_scalar(mean_to_other_vecs_dist_quant, "mean_to_other_vecs_dist_quant")
    if (length(dimension_order) != dim(vectors)[1])
        .tox_shape_error("dimension_order", length(dimension_order), "vectors", dim(vectors)[1])
    if (length(kd_indices) != dim(vectors)[2])
        .tox_shape_error("kd_indices", length(kd_indices), "vectors", dim(vectors)[2])

    .result <- .Call("seeds_call", vectors, kd_indices, dimension_order, mean_to_other_vecs_dist_quant)
    .arguments <- c("vectors", "n_dimensions", "n_vectors", "kd_indices", "dimension_order", "mean_to_other_vecs_dist_quant", "is_seed_mask", "ierr")
    .sources <- c(NA_character_, "vectors", "vectors", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$is_seed_mask
}
