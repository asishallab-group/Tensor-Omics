# The `tox_clustering` suite: `cluster_factor_trajectories_k_means`, `k_means_clustering` and
# `linkage_clustering` can be called from R, return the documented types and shapes, and raise the
# documented errors, naming the argument. Whether the numbers are right is the Fortran suite's job
# (test/mod_test_tox_clustering/), as the coding guide's "Where to Test What" asks; the one value
# each procedure is checked for here is there to show that the binding hands the arrays over in
# the right order, that labels are 1-based, and that `max_iterations` and `method` arrive.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

METHODS <- c("average", "weighted", "ward")

# (dims, points): (1, 1), (2, 2), (8, 8), (11, 11).
.points <- function() matrix(c(1, 1, 2, 2, 8, 8, 11, 11), nrow = 2L)

# (dims, clusters): the first two points.
.centroids <- function() matrix(c(1, 1, 2, 2), nrow = 2L)

# The four samples flattened (factor1_t1, factor2_t1, factor1_t2, ...), one per column; samples 1
# and 3 evolve alike, and 2 and 4.
.flat <- function() matrix(c(1, 10, 3, 30, 5, 50, 20, 2, 10, 0, 0, 1, 3, 12, 5, 32, 7, 52, 21, 3, 11, 3, 4, 3),
                           nrow = 6L)

# (factors, samples, timepoints) = (2, 4, 3): element (f, s, t) is .flat()[2*(t - 1) + f, s].
.trajectories <- function() aperm(array(.flat(), dim = c(2L, 3L, 4L)), c(1L, 3L, 2L))

# (factors*timepoints, clusters): samples 1 and 2 flattened.
.trajectory_centroids <- function() .flat()[, 1:2]

# The example of https://en.wikipedia.org/wiki/UPGMA, points a..e.
.wikipedia_distances <- function() matrix(c(0, 17, 21, 31, 23, 17, 0, 30, 34, 21, 21, 30, 0, 28, 39,
                                            31, 34, 28, 0, 43, 23, 21, 39, 43, 0), nrow = 5L)

# assert_error, and that the error names `argument`. On main, assert_error cannot check the
# argument yet; once it can (its `argument =`), this is assert_error(expr, msg, code, argument = ).
.assert_error_naming <- function(expr, msg, code, argument) {
  err <- tryCatch({ expr; NULL }, error = function(e) e)
  assert_error(if (!is.null(err)) stop(err), msg, code)
  assert_true(identical(err$argument, argument),
              paste0(msg, ": expected '", argument, "' named, but got ",
                     if (is.null(err$argument)) "none" else sprintf("'%s'", toString(err$argument))))
}

.assert_integer_vector <- function(result, name, length) {
  value <- result[[name]]
  assert_true(is.integer(value) && is.null(dim(value)) && length(value) == length,
              paste0(name, ": expected an integer vector of ", length, ", got ", class(value)[1], " of ",
                     length(value)))
}

# =====================
# cluster_factor_trajectories_k_means
# =====================

test_cluster_factor_trajectories_k_means <- function() {
  result <- cluster_factor_trajectories_k_means(.trajectories(), .trajectory_centroids(), 300L)
  assert_true(setequal(names(result), c("centroids", "labels", "label_counts")),
              paste("got", toString(names(result))))
  assert_true(is.double(result$centroids) && identical(dim(result$centroids), c(6L, 2L)),
              paste("centroids: got dim", toString(dim(result$centroids))))
  .assert_integer_vector(result, "labels", 4L)
  .assert_integer_vector(result, "label_counts", 2L)
}

test_cluster_factor_trajectories_k_means_reads_the_trajectories_as_documented <- function() {
  # the one value: the samples cluster as [1, 2, 1, 2] and the centroids end at the flattened
  # means only if the 3-D array arrives in column-major order and is flattened as documented
  result <- cluster_factor_trajectories_k_means(.trajectories(), .trajectory_centroids(), 300L)
  assert_true(identical(result$labels, c(1L, 2L, 1L, 2L)), paste("got", toString(result$labels)))
  assert_true(identical(result$centroids[, 1], c(2, 11, 4, 31, 6, 51)), paste("got", toString(result$centroids[, 1])))
}

test_cluster_factor_trajectories_k_means_rejects_centroids_of_another_length <- function() {
  # the C shim compares the centroids' rows with n_factors*n_timepoints
  .assert_error_naming(cluster_factor_trajectories_k_means(.trajectories(), matrix(0, 5L, 2L), 300L),
                       "centroids of 5 rows for 2 factors and 3 timepoints", ERR_SIZE_MISMATCH, "centroids")
}

test_cluster_factor_trajectories_k_means_rejects_more_clusters_than_samples <- function() {
  # n_clusters is read off centroids, so that is the argument named
  .assert_error_naming(cluster_factor_trajectories_k_means(.trajectories(), matrix(0, 6L, 5L), 300L),
                       "5 clusters for 4 samples", ERR_INVALID_INPUT, "centroids")
}

test_cluster_factor_trajectories_k_means_rejects_nan_trajectories <- function() {
  trajectories <- .trajectories()
  trajectories[2, 3, 2] <- NaN
  .assert_error_naming(cluster_factor_trajectories_k_means(trajectories, .trajectory_centroids(), 300L),
                       "a NaN trajectory", ERR_NAN_INF, "trajectories")
}

test_cluster_factor_trajectories_k_means_rejects_infinite_centroids <- function() {
  centroids <- .trajectory_centroids()
  centroids[6, 2] <- -Inf
  .assert_error_naming(cluster_factor_trajectories_k_means(.trajectories(), centroids, 300L),
                       "an infinite centroid", ERR_NAN_INF, "centroids")
}

test_cluster_factor_trajectories_k_means_rejects_no_factors <- function() {
  # n_factors is read off trajectories
  .assert_error_naming(cluster_factor_trajectories_k_means(array(numeric(0), dim = c(0L, 4L, 3L)),
                                                           matrix(numeric(0), 0L, 2L), 300L),
                       "no factors", ERR_EMPTY_INPUT, "trajectories")
}

test_cluster_factor_trajectories_k_means_rejects_a_matrix_of_trajectories <- function() {
  .assert_error_naming(cluster_factor_trajectories_k_means(matrix(0, 2L, 4L), .trajectory_centroids(), 300L),
                       "2-D trajectories", NULL, "trajectories")
}

test_cluster_factor_trajectories_k_means_rejects_no_iterations <- function() {
  .assert_error_naming(cluster_factor_trajectories_k_means(.trajectories(), .trajectory_centroids(), 0L),
                       "max_iterations 0", ERR_INVALID_INPUT, "max_iterations")
}

test_cluster_factor_trajectories_k_means_rejects_a_fractional_max_iterations <- function() {
  .assert_error_naming(cluster_factor_trajectories_k_means(.trajectories(), .trajectory_centroids(), 2.5),
                       "max_iterations 2.5", NULL, "max_iterations")
}

# =====================
# k_means_clustering
# =====================

test_k_means_clustering <- function() {
  result <- k_means_clustering(.points(), .centroids())
  assert_true(setequal(names(result), c("centroids", "labels", "label_counts")),
              paste("got", toString(names(result))))
  assert_true(is.double(result$centroids) && identical(dim(result$centroids), c(2L, 2L)),
              paste("centroids: got dim", toString(dim(result$centroids))))
  .assert_integer_vector(result, "labels", 4L)
  .assert_integer_vector(result, "label_counts", 2L)
}

test_k_means_clustering_labels_are_1_based <- function() {
  # the one value: two groups, labelled 1 and 2, with their means as the final centroids
  result <- k_means_clustering(.points(), .centroids())
  assert_true(identical(result$labels, c(1L, 1L, 2L, 2L)), paste("got", toString(result$labels)))
  assert_true(identical(result$centroids, matrix(c(1.5, 1.5, 9.5, 9.5), nrow = 2L)),
              paste("got", toString(result$centroids)))
}

test_k_means_clustering_takes_max_iterations <- function() {
  # one iteration stops after the first assignment, before point 2 moves over
  result <- k_means_clustering(.points(), .centroids(), max_iterations = 1L)
  assert_true(identical(result$labels, c(1L, 2L, 2L, 2L)), paste("got", toString(result$labels)))
}

test_k_means_clustering_rejects_more_clusters_than_points <- function() {
  # n_clusters is read off centroids, so that is the argument named
  .assert_error_naming(k_means_clustering(.points(), matrix(0, 2L, 5L)), "5 clusters for 4 points",
                       ERR_INVALID_INPUT, "centroids")
}

test_k_means_clustering_rejects_no_points <- function() {
  # no points is refused by n_clusters' bound [1, n_points], which is checked first, so the
  # argument named is centroids, where n_clusters comes from
  .assert_error_naming(k_means_clustering(matrix(numeric(0), 2L, 0L), .centroids()), "no points",
                       ERR_INVALID_INPUT, "centroids")
}

test_k_means_clustering_rejects_no_dimensions <- function() {
  .assert_error_naming(k_means_clustering(matrix(numeric(0), 0L, 4L), matrix(numeric(0), 0L, 2L)),
                       "no dimensions", ERR_EMPTY_INPUT, "data_points")
}

test_k_means_clustering_rejects_nan_points <- function() {
  points <- .points()
  points[2, 4] <- NaN
  .assert_error_naming(k_means_clustering(points, .centroids()), "a NaN point", ERR_NAN_INF, "data_points")
}

test_k_means_clustering_rejects_infinite_centroids <- function() {
  centroids <- .centroids()
  centroids[1, 1] <- Inf
  .assert_error_naming(k_means_clustering(.points(), centroids), "an infinite centroid", ERR_NAN_INF, "centroids")
}

test_k_means_clustering_rejects_centroids_of_other_dimensions <- function() {
  # the binding's own shape check, before the library is called
  .assert_error_naming(k_means_clustering(.points(), matrix(0, 3L, 2L)), "centroids of 3 dimensions for points of 2",
                       NULL, "centroids")
}

test_k_means_clustering_rejects_non_numeric_points <- function() {
  .assert_error_naming(k_means_clustering(matrix(c("a", "b", "c", "d"), 2L), .centroids()), "strings", NULL,
                       "data_points")
}

test_k_means_clustering_rejects_no_iterations <- function() {
  .assert_error_naming(k_means_clustering(.points(), .centroids(), max_iterations = 0L), "max_iterations 0",
                       ERR_INVALID_INPUT, "max_iterations")
}

test_k_means_clustering_rejects_a_fractional_max_iterations <- function() {
  .assert_error_naming(k_means_clustering(.points(), .centroids(), max_iterations = 1.5), "max_iterations 1.5", NULL,
                       "max_iterations")
}

# =====================
# linkage_clustering
# =====================

test_linkage_clustering <- function() {
  for (method in METHODS) {
    result <- linkage_clustering(.wikipedia_distances(), method)
    assert_true(setequal(names(result), c("distances", "merge_i", "merge_j", "heights", "cluster_sizes")),
                paste(method, ": got", toString(names(result))))
    for (name in c("merge_i", "merge_j", "cluster_sizes")) .assert_integer_vector(result, name, 4L)
    assert_true(is.double(result$heights) && is.null(dim(result$heights)) && length(result$heights) == 4L,
                paste(method, ": heights should be a double vector of 4"))
    assert_true(identical(result$distances, .wikipedia_distances()),
                paste(method, ": the matrix should come back unchanged"))
  }
}

test_linkage_clustering_takes_the_method <- function() {
  # the one value: the last merge's height tells the methods apart (UPGMA 33, WPGMA 35, Ward
  # sqrt(28876/15) = 43.88), and the first merge joins points 1 and 2, labelled 1-based
  expected_last <- c(average = 33, weighted = 35, ward = sqrt(28876 / 15))
  for (method in METHODS) {
    result <- linkage_clustering(.wikipedia_distances(), method)
    last <- result$heights[4]
    assert_true(abs(last - expected_last[[method]]) <= 1e-15 * expected_last[[method]],
                paste(method, ": got", toString(result$heights)))
    assert_true(setequal(c(result$merge_i[1], result$merge_j[1]), c(1L, 2L)),
                paste(method, ": first merge", result$merge_i[1], result$merge_j[1]))
  }
}

test_linkage_clustering_of_one_point_merges_nothing <- function() {
  merges <- linkage_clustering(matrix(0, 1, 1), "average")
  merged <- merges[c("merge_i", "merge_j", "heights", "cluster_sizes")]
  assert_true(all(lengths(merged) == 0L), "one point merges nothing")
}

test_linkage_clustering_of_no_points_is_the_librarys_empty_input <- function() {
  # the merges of n points are n - 1 long, which the shim clamps at zero as Fortran sizes
  # them: R refused -1 itself before the library could report the empty input (#215)
  .assert_error_naming(linkage_clustering(matrix(0, 0, 0), "average"), "no points", ERR_EMPTY_INPUT, "distances")
}

test_linkage_clustering_rejects_an_unknown_method <- function() {
  .assert_error_naming(linkage_clustering(.wikipedia_distances(), "single"), "method 'single'", NULL, "method")
}

test_linkage_clustering_rejects_a_non_square_matrix <- function() {
  .assert_error_naming(linkage_clustering(matrix(0, 2L, 3L), "average"), "a 2 x 3 matrix", NULL, "distances")
}

test_linkage_clustering_rejects_nan <- function() {
  distances <- .wikipedia_distances()
  distances[4, 2] <- distances[2, 4] <- NaN
  .assert_error_naming(linkage_clustering(distances, "ward"), "a NaN distance", ERR_NAN_INF, "distances")
}

test_linkage_clustering_rejects_inf <- function() {
  distances <- .wikipedia_distances()
  distances[5, 1] <- distances[1, 5] <- Inf
  .assert_error_naming(linkage_clustering(distances, "weighted"), "an infinite distance", ERR_NAN_INF, "distances")
}

test_linkage_clustering_rejects_a_negative_distance <- function() {
  distances <- .wikipedia_distances()
  distances[3, 1] <- distances[1, 3] <- -1
  .assert_error_naming(linkage_clustering(distances, "average"), "a negative distance", ERR_INVALID_INPUT, "distances")
}

test_linkage_clustering_rejects_an_asymmetric_matrix <- function() {
  distances <- .wikipedia_distances()
  distances[3, 1] <- 5
  .assert_error_naming(linkage_clustering(distances, "average"), "an asymmetric matrix", ERR_INVALID_INPUT,
                       "distances")
}

test_linkage_clustering_rejects_a_non_zero_diagonal <- function() {
  distances <- .wikipedia_distances()
  distances[2, 2] <- 1
  .assert_error_naming(linkage_clustering(distances, "average"), "a non-zero diagonal", ERR_INVALID_INPUT,
                       "distances")
}

run_all_tests()
