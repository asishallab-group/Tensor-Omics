#!/usr/bin/env python3
"""Group points, or whole trajectories, by k-means and by hierarchical linkage.

Two long-standing, ordinary clustering methods, with Tensor Omics' own
conventions around them:

  k_means_clustering                 partitions points in n dimensions
  cluster_factor_trajectories_k_means the same, with one whole trajectory as
                                     one point, so samples that evolve alike
                                     land together
  linkage_clustering                 agglomerative merging of a distance
                                     matrix you already computed

k-means needs you to supply the initial centroids, which is where k comes from
and where the seeding decision lives -- there is no k-means++ here. Linkage
takes a distance matrix rather than points, so the same matrix can be
re-clustered under a different criterion without recomputing it.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/clustering.py

The data is synthetic and planted: three separated blobs for k-means, and four
trajectories of which two rise and two fall.
"""

import sys

import numpy as np

from tensor_omics import (cluster_factor_trajectories_k_means, k_means_clustering,
                          linkage_clustering)


def show_tree(merge_i, merge_j, heights, sizes):
    """Print the merge sequence, decoding the sign convention.

    A positive entry is a 1-based *data point*; a negative entry is the inner
    node created at that step, so -2 is whatever step 2 produced. Note this is
    the opposite of R's hclust, where negatives are the singletons.
    """
    def name(value):
        return f"point {value}" if value > 0 else f"node {-value}"

    for step, (i, j, height, size) in enumerate(zip(merge_i, merge_j, heights, sizes), start=1):
        print(f"  step {step}: {name(i):9} + {name(j):9} -> node {step}"
              f"   height {height:.3f}, size {size}")


def main(argv=None):
    rng = np.random.default_rng(2)

    # --- k-means over points ------------------------------------------------
    points = np.hstack([
        rng.normal(np.array(centre)[:, None], 0.3, size=(2, 12))
        for centre in ([0, 0], [5, 5], [0, 5])
    ])
    data = np.asfortranarray(points)

    # k is the number of columns here, and these are the seeds. They are
    # modified in place and come back as the final centroids.
    centroids = np.asfortranarray(np.array([[0.0, 5.0, 0.0], [0.0, 5.0, 5.0]]))
    result = k_means_clustering(data, centroids, max_iterations=300)
    print(f"k-means over {data.shape[1]} points in {data.shape[0]}D, k={centroids.shape[1]}")
    print(f"  cluster sizes: {np.asarray(result['label_counts'])}")
    print(f"  final centroids:\n{np.round(centroids, 2)}")

    # --- k-means over whole trajectories ------------------------------------
    n_factors, n_timepoints = 2, 3
    rising = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
    falling = np.array([[3.0, 2.0, 1.0], [3.0, 2.0, 1.0]])
    trajectories = np.zeros((n_factors, 4, n_timepoints), order="F")
    for i_sample, course in enumerate([rising, rising * 1.05, falling, falling * 0.95]):
        trajectories[:, i_sample, :] = course

    # One centroid is a whole flattened trajectory: factor-major within each
    # timepoint, so factor1_t1, factor2_t1, factor1_t2, ...
    seeds = np.asfortranarray(np.column_stack([
        np.array([1.0, 1.0, 2.0, 2.0, 3.0, 3.0]),
        np.array([3.0, 3.0, 2.0, 2.0, 1.0, 1.0]),
    ]))
    grouped = cluster_factor_trajectories_k_means(trajectories, seeds, 100)
    labels = np.asarray(grouped["labels"])
    print(f"\nk-means over {trajectories.shape[1]} trajectories "
          f"({n_factors} factors x {n_timepoints} timepoints each)")
    print(f"  labels: {labels}  (one per sample)")
    print(f"  the two rising samples agree: {labels[0] == labels[1]}, "
          f"the two falling ones agree: {labels[2] == labels[3]}, "
          f"and the groups differ: {labels[0] != labels[2]}")

    # --- hierarchical linkage over a distance matrix -------------------------
    subset = points.T[:8]
    distances = np.sqrt(((subset[:, None, :] - subset[None, :, :]) ** 2).sum(-1))
    print(f"\nlinkage over {len(subset)} points, from a precomputed distance matrix")
    for method in ("average", "weighted", "ward"):
        # The matrix is declared as modified in place, but the lower triangle is
        # only scratch and is restored before returning -- no copy needed.
        matrix = np.asfortranarray(distances.copy())
        tree = linkage_clustering(matrix, method)
        unchanged = np.array_equal(matrix, distances)
        print(f"  {method:9} last merge at height "
              f"{np.asarray(tree['heights'])[-1]:.3f}; matrix unchanged: {unchanged}")

    tree = linkage_clustering(np.asfortranarray(distances.copy()), "average")
    print("\naverage linkage, merge by merge:")
    show_tree(np.asarray(tree["merge_i"]), np.asarray(tree["merge_j"]),
              np.asarray(tree["heights"]), np.asarray(tree["cluster_sizes"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
