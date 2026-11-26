# Generated. Do not edit.

#' Performs k-means clustering on whole factor trajectories, so which samples evolve alike
#'
#' One sample is one point: its trajectory is flattened over factors and time into a single
#' vector, `factor1_t1, factor2_t1, ..., factor1_t2, ...`, so the clustering answers which
#' samples follow a similar course. Clustering the individual `(sample, timepoint)` states
#' instead is \code{\link{k_means_clustering}} on the same
#' array read as `(n_factors, n_samples*n_timepoints)`, which is a different question. The
#' distances and means are scaled as that procedure documents, with the same limits.
#'
#' Generated from the Fortran procedure \code{tox_clustering::cluster_factor_trajectories_k_means}, whose argument names
#' are the ones an error message reports.
#'
#' @param trajectories a numeric array of rank 3. matrix with the trajectories to cluster
#' @param centroids a numeric matrix. matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
#'   The centroids should be unique. This is not checked in this routine.
#'
#'   One centroid is a whole flattened trajectory, so it has `n_factors*n_timepoints` entries.
#'
#'   The final values will be the final centroids of the clusters. A cluster that ends up
#'   without samples (its `label_counts` is 0) keeps the centroid it had when it lost its
#'   last sample, or its starting centroid if it never had one.
#' @param max_iterations a integer scalar. number of maximum iterations of the clustering. Unlike
#'   \code{\link{k_means_clustering}}, which defaults to
#'   300, this one has to be given.
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{centroids}{a numeric matrix. matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
#'     The centroids should be unique. This is not checked in this routine.
#'
#'     One centroid is a whole flattened trajectory, so it has `n_factors*n_timepoints` entries.
#'
#'     The final values will be the final centroids of the clusters. A cluster that ends up
#'     without samples (its `label_counts` is 0) keeps the centroid it had when it lost its
#'     last sample, or its starting centroid if it never had one.}
#'   \item{labels}{a integer vector. array of labels, one per sample, so the first label is the first sample's.
#'
#'     each label is the index of its related cluster -> `1<=label<=n_clusters=k`}
#'   \item{label_counts}{a integer vector. holds the number of samples having the respective label assigned}
#' @export
cluster_factor_trajectories_k_means <- function(trajectories, centroids, max_iterations) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    centroids <- .tox_as_double_matrix(centroids, "centroids")
    max_iterations <- .tox_as_integer_scalar(max_iterations, "max_iterations")
    .result <- .Call("cluster_factor_trajectories_k_means_call", trajectories, centroids, max_iterations)
    .arguments <- c("n_clusters", "trajectories", "n_factors", "n_samples", "n_timepoints", "centroids", "labels", "label_counts", "max_iterations", "ierr")
    .sources <- c("centroids", NA_character_, "trajectories", "trajectories", "trajectories", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .shapes <- c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "n_factors*n_timepoints, n_clusters", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources, shapes = .shapes)

    list(
        centroids = .result$centroids,
        labels = .result$labels,
        label_counts = .result$label_counts
    )
}

#' k-means clustering algorithm
#'
#' 1. Assigns each data point to one of `k` clusters whose centroid is closest
#' 2. Recalculates the centroids using the mean of its assigned points
#' 3. repeat 1-2 until assignment remains unchanged
#'
#' A point exactly as close to several centroids goes to the last of them, the one with the
#' highest cluster index.
#'
#' Each point's squared distances are compared scaled by a power of two picked for that point,
#' by f42's rule for sums of squares, from its largest coordinate difference to the nearest
#' centroid in that measure, and each coordinate of a centroid is the mean of its cluster's
#' values in that dimension, summed scaled by the power of two picked for that set. So for
#' finite input the squared distances that decide a label neither overflow nor vanish beside
#' larger values, and no centroid's sum overflows; only the distances to centroids far from a
#' point may overflow, and those centroids lose. One limit stays: where the build or the host
#' flushes subnormal numbers to zero, coordinates closer than `2**-1022` read as equal.
#'
#' Ordinary data needs no scaling and gets the arithmetic of the plain formulas: its results
#' are bit-identical to those of the code before this scaling in builds that do not contract a
#' multiplication and an addition into one fused operation (FMA); where a build does, a near
#' tie can fall the other way.
#'
#' Generated from the Fortran procedure \code{tox_clustering::k_means_clustering}, whose argument names
#' are the ones an error message reports.
#'
#' @param data_points a numeric matrix. matrix with data points to cluster
#' @param centroids a numeric matrix. matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
#'   The centroids should be unique. This is not checked in this routine.
#'
#'   The final values will be the final centroids of the clusters. A cluster that ends up
#'   without points (its `label_counts` is 0) keeps the centroid it had when it lost its
#'   last point, or its starting centroid if it never had one.
#' @param max_iterations a integer scalar. number of maximum iterations of the clustering. Defaults to 300, while
#'   \code{\link{cluster_factor_trajectories_k_means}}
#'   requires it.
#'   The minimum valid value is `1`.
#'   The default value is `300`.
#' @return a named list with elements:
#'   \item{centroids}{a numeric matrix. matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
#'     The centroids should be unique. This is not checked in this routine.
#'
#'     The final values will be the final centroids of the clusters. A cluster that ends up
#'     without points (its `label_counts` is 0) keeps the centroid it had when it lost its
#'     last point, or its starting centroid if it never had one.}
#'   \item{labels}{a integer vector. array of labels, each index corresponds to the respective point's index, so first label is first point's label.
#'
#'     each label is the index of its related cluster -> `1<=label<=n_clusters=k`}
#'   \item{label_counts}{a integer vector. holds the number of points having the respective label assigned}
#' @export
k_means_clustering <- function(data_points, centroids, max_iterations = 300L) {
    data_points <- .tox_as_double_matrix(data_points, "data_points")
    centroids <- .tox_as_double_matrix(centroids, "centroids")
    max_iterations <- .tox_as_integer_scalar(max_iterations, "max_iterations")
    if (dim(centroids)[1] != dim(data_points)[1])
        .tox_shape_error("centroids", 1L, dim(centroids)[1], "'data_points' implies n_dims == ", dim(data_points)[1], "n_dims, n_clusters")

    .result <- .Call("k_means_clustering_call", data_points, centroids, max_iterations)
    .arguments <- c("n_clusters", "data_points", "n_points", "n_dims", "centroids", "labels", "label_counts", "max_iterations", "ierr")
    .sources <- c("centroids", NA_character_, "data_points", "data_points", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        centroids = .result$centroids,
        labels = .result$labels,
        label_counts = .result$label_counts
    )
}

#' Perform linkage clustering on a distance matrix.
#'
#' The bottom triangle is used as scratch and restored from the top triangle before
#' returning, on success or on error, so the matrix comes back unchanged. There is no
#' need to copy it before calling.
#'
#' UPGMA and WPGMA average each pair of distances scaled by the power of two f42's rule picks
#' for the larger, and Ward combines its three distances scaled for the largest of them, so for
#' finite input no update overflows, and small distances beside large ones keep their value.
#' Ward's distances can still exceed every input distance, as they grow with the clusters'
#' sizes: a merge whose true height passes `huge` stops the run with `ERR_NAN_INF`. When an
#' input exceeds `2**1007`, Ward keeps its scratch distances scaled by `2**-16`, so that one
#' passing `huge` on the way is still compared right; inputs below `2**-1006` are then
#' subnormal in that scratch, and read as 0 where subnormals are flushed.
#'
#' Every update is clamped where rounding can leave the range it lies in mathematically: an
#' average into the pair it averages, a Ward update to at least the merge's height. So the
#' heights never decrease, and equal distances average to themselves. Distances that need no
#' scaling get the plain formulas: their results are bit-identical to those of the code before
#' this scaling in builds that do not contract a multiplication and an addition into one fused
#' operation (FMA), except where the clamps apply.
#'
#' Generated from the Fortran procedure \code{tox_clustering::linkage_clustering}, whose argument names
#' are the ones an error message reports.
#'
#' @param distances a numeric matrix. symmetric distance matrix, holding the non-negative distances between points. Distance of X->X is always zero.
#'
#'   The bottom triangle is used as scratch and restored from the top triangle before
#'   returning, on success or on error, so the matrix comes back unchanged. There is no
#'   need to copy it before calling.
#'
#'   Its structure (symmetry, non-negativity, zero diagonal) is validated by the
#'   distance-matrix naming convention in the generated wrapper.
#' @param method a string, one of "average", "weighted", "ward". used algorithm
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @return a named list with elements:
#'   \item{distances}{a numeric matrix. symmetric distance matrix, holding the non-negative distances between points. Distance of X->X is always zero.
#'
#'     The bottom triangle is used as scratch and restored from the top triangle before
#'     returning, on success or on error, so the matrix comes back unchanged. There is no
#'     need to copy it before calling.
#'
#'     Its structure (symmetry, non-negativity, zero diagonal) is validated by the
#'     distance-matrix naming convention in the generated wrapper.}
#'   \item{merge_i}{a integer vector. the first of the two clusters merged at merge k: a positive label is a point, its
#'     1-based index; a negative label `-j` is the cluster made by merge j. Of the two, this
#'     is the one whose lowest point index is the higher.
#'
#'     The signs are the reverse of R's `hclust`, whose `merge` holds points as negative
#'     and merged clusters as positive numbers.}
#'   \item{merge_j}{a integer vector. the second of the two clusters merged at merge k, labelled as in `merge_i`: the one
#'     whose lowest point index is the lower.
#'
#'     Two pairs at the same distance merge in the order of their lowest point indices: the
#'     pair whose lower cluster has the lowest point index first, then the one whose higher
#'     cluster has.}
#'   \item{heights}{a numeric vector. the linkage distance at which merge k joins its two clusters, by the method's
#'     criterion. The heights never decrease from one merge to the next.}
#'   \item{cluster_sizes}{a integer vector. number of points in the cluster made by merge k}
#' @export
linkage_clustering <- function(distances, method) {
    distances <- .tox_as_double_matrix(distances, "distances")
    method <- .tox_as_mode(method, "method", c("average", "weighted", "ward"))
    if (dim(distances)[2] != dim(distances)[1])
        .tox_shape_error("distances", 2L, dim(distances)[2], "'distances' implies n_points == ", dim(distances)[1], "n_points, n_points")

    .result <- .Call("linkage_clustering_call", distances, method)
    .arguments <- c("distances", "n_points", "merge_i", "merge_j", "heights", "cluster_sizes", "method", "ierr")
    .sources <- c(NA_character_, "distances", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        distances = .result$distances,
        merge_i = .result$merge_i,
        merge_j = .result$merge_j,
        heights = .result$heights,
        cluster_sizes = .result$cluster_sizes
    )
}
