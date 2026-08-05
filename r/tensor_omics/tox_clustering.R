# Generated. Do not edit.

#' Performs k-means clustering on factor trajectories, so factor evolution over time
#'
#' @param trajectories a numeric array of rank 3. matrix with data points to cluster
#' @param centroids a numeric matrix. matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
#'   The centroids should be unique. This is not checked in this routine.
#'
#'   The final values will be the final centroids of the clusters
#' @param max_iterations a integer scalar. number of maximum iterations of the clustering
#' @return a named list with elements `centroids`, `labels`, `label_counts`.
#'
#' Generated from the Fortran procedure \code{tox_clustering::cluster_factor_trajectories_k_means}.
#' @export
cluster_factor_trajectories_k_means <- function(trajectories, centroids, max_iterations) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    centroids <- .tox_as_double_matrix(centroids, "centroids")
    max_iterations <- .tox_as_integer_scalar(max_iterations, "max_iterations")
    if (dim(centroids)[1] != dim(trajectories)[1])
        .tox_shape_error("centroids", dim(centroids)[1], "trajectories", dim(trajectories)[1])

    .result <- .Call("cluster_factor_trajectories_k_means_call", trajectories, centroids, max_iterations)
    .arguments <- c("n_clusters", "trajectories", "n_factors", "n_samples", "n_timepoints", "centroids", "labels", "label_counts", "max_iterations", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        centroids = .result$centroids,
        labels = .result$labels,
        label_counts = .result$label_counts
    )
}

#' k-means clustering algorithm
#'
#' 1. Assigns each data point to one of `k` clusters whose centroid is clostest
#' 2. Recalculates the centroids using the mean of its assigned points
#' 3. repeat 1-2 until assignment remains unchanged
#'
#' @param data_points a numeric matrix. matrix with data points to cluster
#' @param centroids a numeric matrix. matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
#'   The centroids should be unique. This is not checked in this routine.
#'
#'   The final values will be the final centroids of the clusters
#' @param max_iterations a integer scalar. number of maximum iterations of the clustering.
#'   The default value is `300`.
#' @return a named list with elements `centroids`, `labels`, `label_counts`.
#'
#' Generated from the Fortran procedure \code{tox_clustering::k_means_clustering}.
#' @export
k_means_clustering <- function(data_points, centroids, max_iterations = 300L) {
    data_points <- .tox_as_double_matrix(data_points, "data_points")
    centroids <- .tox_as_double_matrix(centroids, "centroids")
    max_iterations <- .tox_as_integer_scalar(max_iterations, "max_iterations")
    if (dim(centroids)[1] != dim(data_points)[1])
        .tox_shape_error("centroids", dim(centroids)[1], "data_points", dim(data_points)[1])

    .result <- .Call("k_means_clustering_call", data_points, centroids, max_iterations)
    .arguments <- c("n_clusters", "data_points", "n_points", "n_dims", "centroids", "labels", "label_counts", "max_iterations", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        centroids = .result$centroids,
        labels = .result$labels,
        label_counts = .result$label_counts
    )
}

#' Perform linkage clustering on a distance matrix.
#'
#' @note
#' This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
#' So there is no need to copy an existing distance matrix, just pass the original.
#' @endnote
#'
#' @param distances a numeric matrix. symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
#'
#'   @note
#'   This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
#'   So there is no need to copy an existing distance matrix, just pass the original.
#'   @endnote
#'
#'   Its structure (symmetry, non-negativity, zero diagonal) is validated by the
#'   distance-matrix naming convention in the generated wrapper.
#' @param method a string, one of "average", "weighted", "ward". used algorithm
#'   The minimum valid value is `0`.
#'   The maximum valid value is `2`.
#' @return a named list with elements `distances`, `merge_i`, `merge_j`, `heights`, `cluster_sizes`.
#'
#' Generated from the Fortran procedure \code{tox_clustering::linkage_clustering}.
#' @export
linkage_clustering <- function(distances, method) {
    distances <- .tox_as_double_matrix(distances, "distances")
    method <- .tox_as_mode(method, "method", c("average", "weighted", "ward"))
    .result <- .Call("linkage_clustering_call", distances, method)
    .arguments <- c("distances", "n_points", "merge_i", "merge_j", "heights", "cluster_sizes", "method", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        distances = .result$distances,
        merge_i = .result$merge_i,
        merge_j = .result$merge_j,
        heights = .result$heights,
        cluster_sizes = .result$cluster_sizes
    )
}
