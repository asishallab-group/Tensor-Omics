source("r/load_tensor_omics.R")
source("r/test_helpers.R")

test_k_means_clustering <- function() {
  # k_means_clustering returns valid output
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
  # the extents come from the arrays themselves now
  res <- k_means_clustering(
    data_points = data_points,
    centroids = centroids,
    max_iterations = max_iterations
  )
  assert_true(isTRUE(all.equal(dim(res$centroids), c(n_dims, n_clusters))))
  assert_true(length(res$labels) == n_points)
  assert_true(length(res$label_counts) == n_clusters)

  # Example: cluster_factor_trajectories_k_means (minimal smoke test)

  # cluster_factor_trajectories_k_means runs without error
  set.seed(1)
  n_clusters <- 2
  n_factors <- 2
  n_samples <- 2
  n_timepoints <- 2
  trajectories <- array(c(
    1, 2, 3, 4, 5, 6, 7, 8
  ), dim = c(n_factors, n_samples, n_timepoints))
  # one sample is one point, its trajectory flattened over factors and time
  n_traj_dims <- n_factors * n_timepoints
  centroids <- matrix(runif(n_traj_dims * n_clusters), nrow = n_traj_dims, ncol = n_clusters)
  max_iterations <- 5
  # the extents come from the arrays, so the trajectories keep their shape
  res <- cluster_factor_trajectories_k_means(
    trajectories = trajectories,
    centroids = centroids,
    max_iterations = max_iterations
  )
  assert_true(isTRUE(all.equal(dim(res$centroids), c(n_traj_dims, n_clusters))))
  assert_true(length(res$labels) == n_samples)
  assert_true(length(res$label_counts) == n_clusters)
}

run_all_tests()
