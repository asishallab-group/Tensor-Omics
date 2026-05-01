import time
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import build_bst_index, build_kd_index, bst_range_query, build_spherical_kd
from test_helpers import run_all_tests


# --- Test Cases ---
def test_bst():
    x = np.array([3.0, 1.0, 4.0, 2.0], dtype=np.float64)

    # Build BST index using the wrapper function
    ix = build_bst_index(x)

    assert all(x[ix - 1] == sorted(x)), "expected x to be sorted"

    # Range query using the wrapper function
    res = bst_range_query(x, ix, 1.5, 3.5)
    matching_indices = res["matching_indices"]

    assert all(x[matching_indices - 1] == [2.0, 3.0]), "expected range values are wrong"


def test_kdtree():
    X = np.array([[4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0]], dtype=np.float64).transpose().copy(order="F")
    for i in range(X.shape[1]):
        X[:, i] = [5.0, 6.0, 3.0]
        kd_ix = build_kd_index(X, np.array([1, 2, 3], dtype=np.int32))
        assert kd_ix[-1] == i + 1, "test vector should be always last"
        X[:, i] = [4.0, 2.0, 1.0]


def test_spherical_kdtree():

    # Create some unit vectors on a sphere
    np.random.seed(42)  # For reproducible results
    vectors = np.random.randn(3, 10).astype(np.float64, order="F")  # 3D, 10 vectors
    norms = np.linalg.norm(vectors, axis=0)
    unit_vectors = vectors / norms

    X = np.array([[4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0], [4.0, 2.0, 1.0]], dtype=np.float64).transpose().copy(order="F")
    X /= np.linalg.norm(X, axis=0)

    for i in range(X.shape[1]):
        X[:, i] = [5.0, 6.0, 3.0] / np.linalg.norm([5.0, 6.0, 3.0])
        kd_ix = build_kd_index(X, np.array([1, 2, 3], dtype=np.int32))
        assert kd_ix[1] == i + 1, "test vector should be always second"
        X[:, i] = [4.0, 2.0, 1.0] / np.linalg.norm([4.0, 2.0, 1.0])


def test_spherical_kdtree_specific_cases():
    # Test case 1: Points on a sphere (more realistic spherical data)
    theta = np.linspace(0, 2*np.pi, 8)
    phi = np.linspace(0, np.pi, 4)

    sphere_points = []
    for p in phi:
        for t in theta:
            x = np.sin(p) * np.cos(t)
            y = np.sin(p) * np.sin(t)
            z = np.cos(p)
            sphere_points.append([x, y, z])

    sphere_points = np.array([[
        np.sin(p) * np.cos(t),
        np.sin(p) * np.sin(t),
        np.cos(p)
    ] for t in theta for p in phi], dtype=np.float64).T.copy(order="F")

    sphere_ix = build_spherical_kd(sphere_points, np.array([1, 2, 3], dtype=np.int32))
    hemisphere_points = sphere_points[:, sphere_points[2, :] >= 0].copy(order="F")  # Upper hemisphere
    hemisphere_ix = build_spherical_kd(hemisphere_points, np.array([1, 2, 3], dtype=np.int32))


def test_bst_edge_cases():
    # Empty array
    try:
        x = np.array([], dtype=np.float64, order="F")
        ix = build_bst_index(x)
        assert False, "Expected error for empty bst index input"
    except Exception as e:
        pass

    # Single element
    x = np.array([42.0], dtype=np.float64)
    ix = build_bst_index(x)


def test_kdtree_edge_cases():
    # Empty matrix
    try:
        X = np.empty((2, 0), dtype=np.float64, order='F')
        kd_ix = build_kd_index(X)
        assert False, "Expected error for empty kd tree input"
    except Exception as e:
        pass

    # Single point
    X = np.array([[1.0, 2.0]], dtype=np.float64).T.copy(order='F')
    kd_ix = build_kd_index(X, np.array([1, 2], dtype=np.int32))


def test_spherical_kdtree_edge_cases():
    # Empty spherical data
    try:
        empty_vectors = np.empty((3, 0), dtype=np.float64, order='F')
        sphere_ix = build_spherical_kd(empty_vectors) - 1
        assert False, "Expected error for empty spherical kd tree input"
    except Exception as e:
        pass

    # Single vector on sphere
    single_vector = np.array([[0.0, 0.0, 1.0]], dtype=np.float64).T.copy(order='F')  # North pole
    sphere_ix = build_spherical_kd(single_vector, np.array([1, 2, 3], dtype=np.int32)) - 1

    # 2D case (not truly spherical but should work)
    circle_vectors = np.array([[1.0, 0.0, -1.0, 0.0], 
                             [0.0, 1.0, 0.0, -1.0]], dtype=np.float64, order='F')
    circle_vectors = circle_vectors / np.linalg.norm(circle_vectors, axis=0)  # Normalize
    circle_ix = build_spherical_kd(circle_vectors, np.array([1, 2], dtype=np.int32))


def test_performance():
    # Test with larger datasets
    # Larger BST
    large_values = np.random.rand(1000).astype(np.float64)
    bst_indices = build_bst_index(large_values)

    # Larger spherical KD-Tree
    large_vectors = np.random.randn(3, 500).astype(np.float64, order="F")  # 3D, 500 vectors
    large_vectors = large_vectors / np.linalg.norm(large_vectors, axis=0)
    start_time = time.time()
    sphere_indices = build_spherical_kd(large_vectors)
    end_time = time.time()

    assert end_time - start_time < 1, "performance very bad, should be faster than a second"


if __name__ == "__main__":
    run_all_tests(globals().values())
