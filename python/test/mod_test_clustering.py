"""
Test script for clock hand angle functions
Python equivalent of the R and Fortran clock hand angle tests
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_k_means_clustering, tox_linkage_clustering


def test_tox_k_means_clustering():
    """
    Test tox_k_means_clustering using known 2D points and initial centroids.
    """

    # Define 2D points: two clusters
    # Cluster 1: (1,1), (2,2)
    # Cluster 2: (8,8), (10,10)
    data_points = np.array([
        [1.0, 2.0, 8.0, 10.0],  # x
        [1.0, 2.0, 8.0, 10.0]   # y
    ], dtype=np.float64)

    # Initial centroids: first and last points
    centroids_init = np.array([
        [1.0, 10.0],
        [1.0, 10.0]
    ], dtype=np.float64)

    # Expected final centroids
    expected_centroids = np.array([
        [1.5, 9.0],
        [1.5, 9.0]
    ], dtype=np.float64)

    # Expected labels and counts
    expected_labels = np.array([1, 1, 2, 2], dtype=np.int32)
    expected_counts = np.array([2, 2], dtype=np.int32)

    # Run clustering
    centroids, labels, label_counts = tox_k_means_clustering(data_points, centroids_init, max_iter=10).values()

    # Validate centroids
    assert np.allclose(centroids, expected_centroids, atol=1e-12), "tox_k_means_clustering: centroid mismatch"

    # Validate labels
    assert np.array_equal(labels, expected_labels), "tox_k_means_clustering: label assignment mismatch"

    # Validate label counts
    assert np.array_equal(label_counts, expected_counts), "tox_k_means_clustering: label_counts mismatch"

    print("✅ tox_k_means_clustering passed all tests.")


def test_linkage_methods():
    import numpy as np

    TOL = 1e-12
    method_names = ["average", "weighted", "ward"]
    method_labels = ["UPGMA", "WPGMA", "Ward"]

    expected_merge_i = np.array([
        [1, -1, 3, -3],
        [1, -1, 3, -3],
        [1, -1, 3, -3]
    ])
    expected_merge_j = np.array([
        [2, 5, 4, -2],
        [2, 5, 4, -2],
        [2, 5, 4, -2]
    ])
    expected_heights = np.array([
        [17.0, 22.0, 28.0, 33.0],
        [17.0, 22.0, 28.0, 35.0],
        [17.0, 23.459184413217212, 28.0, 43.87558166755931]
    ])
    expected_cluster_sizes = np.array([
        [2, 3, 2, 5],
        [2, 3, 2, 5],
        [2, 3, 2, 5]
    ])

    for i_method, method in enumerate(method_names):
        try:
            label = method_labels[i_method]

            # Case 1: Reference Wikipedia example
            dist = np.array([
                [0.0, 17.0, 21.0, 31.0, 23.0],
                [17.0, 0.0, 30.0, 34.0, 21.0],
                [21.0, 30.0, 0.0, 28.0, 39.0],
                [31.0, 34.0, 28.0, 0.0, 43.0],
                [23.0, 21.0, 39.0, 43.0, 0.0]
            ])
            dist_copy = dist.copy()
            merge_i, merge_j, heights, cluster_sizes = tox_linkage_clustering(dist_copy, method).values()

            assert np.allclose(dist_copy, dist, atol=0.0), f"{label}: reference output matrix mismatch"
            for i in range(4):
                mi, mj = merge_i[i], merge_j[i]
                assert min(mi, mj) == expected_merge_i[i_method, i], f"{label}: merge_i[{i}] mismatch"
                assert max(mi, mj) == expected_merge_j[i_method, i], f"{label}: merge_j[{i}] mismatch"
                assert mi != mj, f"{label}: merge_i == merge_j at step {i}"
            assert np.allclose(heights, expected_heights[i_method], atol=TOL), f"{label}: reference heights mismatch"
            assert np.array_equal(cluster_sizes, expected_cluster_sizes[i_method]), f"{label}: reference cluster_sizes mismatch"

            # Case 2: Equal distances
            dist = np.ones((4, 4)) - np.eye(4)
            dist_copy = dist.copy()
            merge_i, merge_j, heights, cluster_sizes = tox_linkage_clustering(dist_copy, method).values()
            assert np.allclose(dist_copy, dist, atol=0.0), f"{label}: equal-distance output matrix mismatch"
            assert np.allclose(heights, [1.0, 1.0, 1.0], atol=TOL), f"{label}: equal-distance heights mismatch"

            # Case 3: Two points
            dist = np.array([
                [0.0, 5.0],
                [5.0, 0.0]
            ])
            dist_copy = dist.copy()
            merge_i, merge_j, heights, cluster_sizes = tox_linkage_clustering(dist_copy, method).values()
            assert np.allclose(dist_copy, dist, atol=0.0), f"{label}: two-point output matrix mismatch"
            assert np.allclose(heights, [5.0], atol=TOL), f"{label}: two-point height mismatch"

            # Case 4: Single point
            dist = np.array([[0.0]])
            merge_i, merge_j, heights, cluster_sizes = tox_linkage_clustering(dist.copy(), method).values()
            assert merge_i.size == 0, f"{label}: single-point merge_i should be empty"
            assert merge_j.size == 0, f"{label}: single-point merge_j should be empty"
            assert heights.size == 0, f"{label}: single-point heights should be empty"
            assert cluster_sizes.size == 0, f"{label}: single-point cluster_sizes should be empty"

            # Case 5: NaN in distance matrix
            dist = np.array([
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
                [np.nan, 1.0, 0.0]
            ])
            try:
                tox_linkage_clustering(dist.copy(), method)
                raise AssertionError(f"{label}: NaN case should raise an error")
            except RuntimeError:
                pass  # Expected
        except RuntimeError:
            raise AssertionError(f"{label}: Unexpected Error")

    print("✅ tox_linkage_clustering passed all tests.")


def main():
    print("=================================================")
    print("    CLUSTERING PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_tox_k_means_clustering()
    test_linkage_methods()


if __name__ == "__main__":
    main()
