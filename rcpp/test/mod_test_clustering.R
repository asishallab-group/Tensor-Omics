source("rcpp/tensoromics_functions.R")
cat(" Tensoromics clustering R wrappers test\n")
cat("=== Testing tox clustering R wrappers ===\n")
# tox_k_means_clustering returns valid output
set.seed(42)
n_clusters <- 2
n_points <- 6
n_dims <- 2
data_points <- matrix(c(
  1, 2, 1, 2, 8, 9,
  1, 2, 2, 1, 8, 9
), nrow = n_dims, ncol = n_points, byrow = TRUE)
centroids <- matrix(c(1, 1, 8, 8), nrow = n_dims, ncol = n_clusters)
max_iterations <- 10
res <- tox_k_means_clustering(
  n_clusters = n_clusters,
  data_points = data_points,
  n_points = n_points,
  n_dims = n_dims,
  centroids = centroids,
  max_iterations = max_iterations
)
stopifnot(isTRUE(all.equal(dim(res$centroids), c(n_dims, n_clusters))))
stopifnot(length(res$labels) == n_points)
stopifnot(length(res$label_counts) == n_clusters)
stopifnot(res$ierr == 0)

# Example: cluster_factor_trajectories_k_means (minimal smoke test)

# tox_cluster_factor_trajectories_k_means runs without error
set.seed(1)
n_clusters <- 2
n_factors <- 2
n_samples <- 2
n_timepoints <- 2
trajectories <- array(c(
  1, 2, 3, 4, 5, 6, 7, 8
), dim = c(n_factors, n_samples, n_timepoints))
centroids <- matrix(runif(n_factors * n_clusters), nrow = n_factors, ncol = n_clusters)
max_iterations <- 5
res <- tox_cluster_factor_trajectories_k_means(
  n_clusters = n_clusters,
  trajectories = as.numeric(trajectories),
  n_factors = n_factors,
  n_samples = n_samples,
  n_timepoints = n_timepoints,
  centroids = centroids,
  max_iterations = max_iterations
)
cat("dim(res$centroids):", paste(dim(res$centroids), collapse=", "), "\n")
cat("n_factors:", n_factors, "\n")
cat("n_clusters:", n_clusters, "\n")
stopifnot(isTRUE(all.equal(dim(res$centroids), c(n_factors, n_clusters))))
stopifnot(length(res$labels) == n_samples * n_timepoints)
stopifnot(length(res$label_counts) == n_clusters)
stopifnot(res$ierr == 0)

cat("\nAll tox_clustering tests completed successfully.\n")
