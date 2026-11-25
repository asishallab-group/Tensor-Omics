"""
The `tox_clustering` suite: `cluster_factor_trajectories_k_means`, `k_means_clustering` and
`linkage_clustering` can be called from Python, return the documented types and shapes, and raise
the documented errors, naming the argument. Whether the numbers are right is the Fortran suite's
job (test/mod_test_tox_clustering/), as the coding guide's "Where to Test What" asks; the one
value each procedure is checked for here is there to show that the binding hands the arrays over
in the right order, that labels are 1-based, and that `max_iterations` and `method` arrive.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import cluster_factor_trajectories_k_means, k_means_clustering, linkage_clustering
from test_helpers import run_all_tests, assert_error, assert_refused_by_name
from tensor_omics.error_handling import (ERR_EMPTY_INPUT, ERR_INVALID_INPUT, ERR_NAN_INF,
                                         ERR_SIZE_MISMATCH)

METHODS = ("average", "weighted", "ward")


def _points():
    """(dims, points): (1, 1), (2, 2), (8, 8), (11, 11)."""
    return np.array([[1.0, 2.0, 8.0, 11.0], [1.0, 2.0, 8.0, 11.0]], order='F')


def _centroids():
    """(dims, clusters): the first two points, a fresh array each time, as it is overwritten."""
    return np.array([[1.0, 2.0], [1.0, 2.0]], order='F')


def _trajectories():
    """(factors, samples, timepoints) = (2, 4, 3); samples 1 and 3 evolve alike, and 2 and 4."""
    flat = np.array([[1, 10, 3, 30, 5, 50], [20, 2, 10, 0, 0, 1], [3, 12, 5, 32, 7, 52],
                     [21, 3, 11, 3, 4, 3]], dtype=np.float64)
    # sample s flattened is factor1_t1, factor2_t1, factor1_t2, ...: element (f, s, t) is
    # flat[s, 2*t + f]
    return np.asfortranarray(flat.reshape(4, 3, 2).transpose(2, 0, 1))


def _trajectory_centroids():
    """(factors*timepoints, clusters): samples 1 and 2 flattened."""
    return np.array([[1, 10, 3, 30, 5, 50], [20, 2, 10, 0, 0, 1]], dtype=np.float64).T.copy(order='F')


def _wikipedia_distances():
    """The example of https://en.wikipedia.org/wiki/UPGMA, points a..e."""
    return np.array([[0, 17, 21, 31, 23], [17, 0, 30, 34, 21], [21, 30, 0, 28, 39],
                     [31, 34, 28, 0, 43], [23, 21, 39, 43, 0]], dtype=np.float64, order='F')


def _assert_read_only_int32(result, name, shape):
    array = result[name]
    assert isinstance(array, np.ndarray), f"{name}: expected a numpy array, got {type(array).__name__}"
    assert array.dtype == np.int32, f"{name}: expected int32, got {array.dtype}"
    assert array.shape == shape, f"{name}: expected shape {shape}, got {array.shape}"
    assert not array.flags.writeable, f"{name}: a result is a value, it should be read-only"


# =====================
# cluster_factor_trajectories_k_means
# =====================

def test_cluster_factor_trajectories_k_means():
    centroids = _trajectory_centroids()
    result = cluster_factor_trajectories_k_means(_trajectories(), centroids, 300)
    assert set(result) == {"labels", "label_counts"}, f"got keys {sorted(result)}"
    _assert_read_only_int32(result, "labels", (4,))
    _assert_read_only_int32(result, "label_counts", (2,))
    assert centroids.shape == (6, 2) and centroids.dtype == np.float64, "centroids keep their shape and type"


def test_cluster_factor_trajectories_k_means_reads_the_trajectories_as_documented():
    # the one value: the samples cluster as [1, 2, 1, 2] and the centroids end at the flattened
    # means only if the 3-D array arrives in column-major order and is flattened as documented
    centroids = _trajectory_centroids()
    result = cluster_factor_trajectories_k_means(_trajectories(), centroids, 300)
    assert np.array_equal(result["labels"], [1, 2, 1, 2]), f"got {result['labels']}"
    assert np.array_equal(centroids[:, 0], [2, 11, 4, 31, 6, 51]), f"got {centroids[:, 0]}"


def test_cluster_factor_trajectories_k_means_rejects_centroids_of_another_length():
    # the C wrapper compares the centroids' rows with n_factors*n_timepoints
    assert_error(lambda: cluster_factor_trajectories_k_means(_trajectories(), np.zeros((5, 2), order='F'), 300),
                 "centroids of 5 rows for 2 factors and 3 timepoints", ERR_SIZE_MISMATCH, argument="centroids")


def test_cluster_factor_trajectories_k_means_rejects_more_clusters_than_samples():
    # n_clusters is read off centroids, so that is the argument named
    assert_error(lambda: cluster_factor_trajectories_k_means(_trajectories(), np.zeros((6, 5), order='F'), 300),
                 "5 clusters for 4 samples", ERR_INVALID_INPUT, argument="centroids")


def test_cluster_factor_trajectories_k_means_rejects_nan_trajectories():
    trajectories = _trajectories()
    trajectories[1, 2, 1] = np.nan
    assert_error(lambda: cluster_factor_trajectories_k_means(trajectories, _trajectory_centroids(), 300),
                 "a NaN trajectory", ERR_NAN_INF, argument="trajectories")


def test_cluster_factor_trajectories_k_means_rejects_infinite_centroids():
    centroids = _trajectory_centroids()
    centroids[-1, -1] = -np.inf
    assert_error(lambda: cluster_factor_trajectories_k_means(_trajectories(), centroids, 300),
                 "an infinite centroid", ERR_NAN_INF, argument="centroids")


def test_cluster_factor_trajectories_k_means_rejects_no_factors():
    # n_factors is read off trajectories
    assert_error(lambda: cluster_factor_trajectories_k_means(np.empty((0, 4, 3), order='F'),
                                                             np.empty((0, 2), order='F'), 300),
                 "no factors", ERR_EMPTY_INPUT, argument="trajectories")


def test_cluster_factor_trajectories_k_means_rejects_a_matrix_of_trajectories():
    assert_refused_by_name(lambda: cluster_factor_trajectories_k_means(np.zeros((2, 4), order='F'),
                                                                       _trajectory_centroids(), 300),
                           "2-D trajectories", "trajectories")


def test_cluster_factor_trajectories_k_means_rejects_row_major_centroids():
    # centroids is written in place, so the binding cannot convert it
    assert_refused_by_name(lambda: cluster_factor_trajectories_k_means(_trajectories(),
                                                                       np.ascontiguousarray(_trajectory_centroids()), 300),
                           "row-major centroids", "centroids")


def test_cluster_factor_trajectories_k_means_rejects_a_fractional_max_iterations():
    assert_refused_by_name(lambda: cluster_factor_trajectories_k_means(_trajectories(), _trajectory_centroids(), 2.5),
                           "max_iterations 2.5", "max_iterations")


# =====================
# k_means_clustering
# =====================

def test_k_means_clustering():
    centroids = _centroids()
    result = k_means_clustering(_points(), centroids)
    assert set(result) == {"labels", "label_counts"}, f"got keys {sorted(result)}"
    _assert_read_only_int32(result, "labels", (4,))
    _assert_read_only_int32(result, "label_counts", (2,))
    assert centroids.shape == (2, 2) and centroids.dtype == np.float64, "centroids keep their shape and type"


def test_k_means_clustering_labels_are_1_based_and_centroids_written_in_place():
    # the one value: two groups, labelled 1 and 2, whose means land in the caller's array
    centroids = _centroids()
    result = k_means_clustering(_points(), centroids)
    assert np.array_equal(result["labels"], [1, 1, 2, 2]), f"got {result['labels']}"
    assert np.array_equal(centroids, [[1.5, 9.5], [1.5, 9.5]]), f"got {centroids}"


def test_k_means_clustering_takes_max_iterations():
    # one iteration stops after the first assignment, before point 2 moves over
    result = k_means_clustering(_points(), _centroids(), max_iterations=1)
    assert np.array_equal(result["labels"], [1, 2, 2, 2]), f"got {result['labels']}"


def test_k_means_clustering_rejects_more_clusters_than_points():
    # n_clusters is read off centroids, so that is the argument named
    assert_error(lambda: k_means_clustering(_points(), np.zeros((2, 5), order='F')),
                 "5 clusters for 4 points", ERR_INVALID_INPUT, argument="centroids")


def test_k_means_clustering_rejects_no_points():
    # no points is refused by n_clusters' bound [1, n_points], which is checked first, so the
    # argument named is centroids, where n_clusters comes from
    assert_error(lambda: k_means_clustering(np.empty((2, 0), order='F'), _centroids()),
                 "no points", ERR_INVALID_INPUT, argument="centroids")


def test_k_means_clustering_rejects_no_dimensions():
    assert_error(lambda: k_means_clustering(np.empty((0, 4), order='F'), np.empty((0, 2), order='F')),
                 "no dimensions", ERR_EMPTY_INPUT, argument="data_points")


def test_k_means_clustering_rejects_nan_points():
    points = _points()
    points[1, 3] = np.nan
    assert_error(lambda: k_means_clustering(points, _centroids()), "a NaN point", ERR_NAN_INF, argument="data_points")


def test_k_means_clustering_rejects_infinite_centroids():
    centroids = _centroids()
    centroids[0, 0] = np.inf
    assert_error(lambda: k_means_clustering(_points(), centroids), "an infinite centroid", ERR_NAN_INF,
                 argument="centroids")


def test_k_means_clustering_rejects_centroids_of_other_dimensions():
    # the binding's own shape check, before the library is called
    assert_refused_by_name(lambda: k_means_clustering(_points(), np.zeros((3, 2), order='F')),
                           "centroids of 3 dimensions for points of 2", "centroids")


def test_k_means_clustering_rejects_centroids_it_cannot_write_in_place():
    assert_refused_by_name(lambda: k_means_clustering(_points(), [[1.0, 2.0], [1.0, 2.0]]),
                           "a list of centroids", "centroids")
    assert_refused_by_name(lambda: k_means_clustering(_points(), np.array([[1, 2], [1, 2]], order='F')),
                           "integer centroids", "centroids")
    assert_refused_by_name(lambda: k_means_clustering(_points(), np.ascontiguousarray([[1.0, 2.0], [3.0, 4.0]])),
                           "row-major centroids", "centroids")


def test_k_means_clustering_rejects_non_numeric_points():
    assert_refused_by_name(lambda: k_means_clustering(np.array([["a", "b"], ["c", "d"]]), _centroids()),
                           "strings", "data_points")


def test_k_means_clustering_rejects_a_fractional_max_iterations():
    assert_refused_by_name(lambda: k_means_clustering(_points(), _centroids(), max_iterations=1.5),
                           "max_iterations 1.5", "max_iterations")


# =====================
# linkage_clustering
# =====================

def test_linkage_clustering():
    for method in METHODS:
        distances = _wikipedia_distances()
        result = linkage_clustering(distances, method)
        assert set(result) == {"merge_i", "merge_j", "heights", "cluster_sizes"}, f"{method}: got keys {sorted(result)}"
        for name in ("merge_i", "merge_j", "cluster_sizes"):
            _assert_read_only_int32(result, name, (4,))
        heights = result["heights"]
        assert heights.dtype == np.float64 and heights.shape == (4,), f"{method}: heights {heights.dtype} {heights.shape}"
        assert not heights.flags.writeable, f"{method}: heights should be read-only"
        assert np.array_equal(distances, _wikipedia_distances()), f"{method}: the matrix should come back unchanged"


def test_linkage_clustering_takes_the_method():
    # the one value: the last merge's height tells the methods apart (UPGMA 33, WPGMA 35, Ward
    # sqrt(28876/15) = 43.88), and the first merge joins points 1 and 2, labelled 1-based
    expected_last = {"average": 33.0, "weighted": 35.0, "ward": np.sqrt(28876.0 / 15.0)}
    for method, last in expected_last.items():
        result = linkage_clustering(_wikipedia_distances(), method)
        assert np.isclose(result["heights"][-1], last, rtol=1e-15, atol=0.0), f"{method}: got {result['heights']}"
        assert {result["merge_i"][0], result["merge_j"][0]} == {1, 2}, \
            f"{method}: first merge {result['merge_i'][0]}, {result['merge_j'][0]}"


def test_linkage_clustering_of_one_point_merges_nothing():
    merges = linkage_clustering(np.zeros((1, 1), order='F'), "average")
    assert all(merges[name].shape == (0,) for name in merges), "one point merges nothing"


def test_linkage_clustering_of_no_points_is_the_librarys_empty_input():
    # the merges of n points are n - 1 long, which the binding clamps at zero as Fortran sizes
    # them: numpy refused -1 itself before the library could report the empty input (#215)
    assert_error(lambda: linkage_clustering(np.empty((0, 0), order='F'), "average"),
                 "no points", ERR_EMPTY_INPUT, argument="distances")


def test_linkage_clustering_rejects_an_unknown_method():
    assert_refused_by_name(lambda: linkage_clustering(_wikipedia_distances(), "single"), "method 'single'", "method")


def test_linkage_clustering_rejects_a_non_square_matrix():
    assert_refused_by_name(lambda: linkage_clustering(np.zeros((2, 3), order='F'), "average"),
                           "a 2 x 3 matrix", "distances")


def test_linkage_clustering_rejects_a_row_major_matrix():
    # distances is used as scratch in place, so the binding cannot convert it
    assert_refused_by_name(lambda: linkage_clustering(np.ascontiguousarray(_wikipedia_distances()), "average"),
                           "a row-major matrix", "distances")


def test_linkage_clustering_rejects_nan():
    distances = _wikipedia_distances()
    distances[3, 1] = distances[1, 3] = np.nan
    assert_error(lambda: linkage_clustering(distances, "ward"), "a NaN distance", ERR_NAN_INF, argument="distances")


def test_linkage_clustering_rejects_inf():
    distances = _wikipedia_distances()
    distances[4, 0] = distances[0, 4] = np.inf
    assert_error(lambda: linkage_clustering(distances, "weighted"), "an infinite distance", ERR_NAN_INF,
                 argument="distances")


def test_linkage_clustering_rejects_a_negative_distance():
    distances = _wikipedia_distances()
    distances[2, 0] = distances[0, 2] = -1.0
    assert_error(lambda: linkage_clustering(distances, "average"), "a negative distance", ERR_INVALID_INPUT,
                 argument="distances")


def test_linkage_clustering_rejects_an_asymmetric_matrix():
    distances = _wikipedia_distances()
    distances[2, 0] = 5.0
    assert_error(lambda: linkage_clustering(distances, "average"), "an asymmetric matrix", ERR_INVALID_INPUT,
                 argument="distances")


def test_linkage_clustering_rejects_a_non_zero_diagonal():
    distances = _wikipedia_distances()
    distances[1, 1] = 1.0
    assert_error(lambda: linkage_clustering(distances, "average"), "a non-zero diagonal", ERR_INVALID_INPUT,
                 argument="distances")


if __name__ == "__main__":
    run_all_tests(globals().values())
