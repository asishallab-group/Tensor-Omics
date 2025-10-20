import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import build_bst_index, build_kd_index, bst_range_query, build_spherical_kd

# --- Test Cases ---
def test_bst():
    print("=== Testing BST ===")
    x = np.array([3.0, 1.0, 4.0, 2.0], dtype=np.float64)
    
    # Build BST index using the wrapper function
    ix = build_bst_index(x)

    print(f"BST indices (Python 0-based): {ix}")
    print(f"Sorted values: {x[ix]}")

    # Range query using the wrapper function
    matching_indices, count = bst_range_query(x, ix, 1.5, 3.5)
    print(f"Range [1.5, 3.5] matches: {matching_indices} (values: {x[matching_indices]})")

def test_kdtree():
    print("\n=== Testing KD-Tree ===")
    X = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float64, order='F')
    
    # Build KD-Tree using the wrapper function
    kd_ix = build_kd_index(X, np.array([1, 2], dtype=np.int32))
    print(f"KD-Tree indices (Python 0-based): {kd_ix}")
    print(f"Points in order:\n{X[:, kd_ix]}")

def test_spherical_kdtree():
    print("\n=== Testing Spherical KD-Tree ===")
    
    # Create some unit vectors on a sphere
    np.random.seed(42)  # For reproducible results
    vectors = np.random.randn(3, 10)  # 3D, 10 vectors
    # Normalize to unit length
    norms = np.linalg.norm(vectors, axis=0)
    unit_vectors = vectors / norms
    
    print(f"Unit vectors shape: {unit_vectors.shape}")
    print(f"Norms: {np.linalg.norm(unit_vectors, axis=0)}")  # Should all be ~1.0
    
    # Build spherical KD-Tree using the wrapper function
    sphere_ix = build_spherical_kd(unit_vectors, np.array([1, 2, 3], dtype=np.int32))
    print(f"Spherical KD-Tree indices (Python 0-based): {sphere_ix}")
    print(f"Vectors in spherical order:\n{unit_vectors[:, sphere_ix]}")
    
    # Verify that the vectors are still unit length after processing
    processed_norms = np.linalg.norm(unit_vectors[:, sphere_ix], axis=0)
    print(f"Norms after processing: {processed_norms}")
    assert np.allclose(processed_norms, 1.0), "Vectors should remain unit length"

def test_spherical_kdtree_specific_cases():
    print("\n=== Testing Spherical KD-Tree Specific Cases ===")
    
    # Test case 1: Points on a sphere (more realistic spherical data)
    print("Test 1: Points on sphere surface")
    theta = np.linspace(0, 2*np.pi, 8)
    phi = np.linspace(0, np.pi, 4)
    
    sphere_points = []
    for p in phi:
        for t in theta:
            x = np.sin(p) * np.cos(t)
            y = np.sin(p) * np.sin(t)
            z = np.cos(p)
            sphere_points.append([x, y, z])
    
    sphere_points = np.array(sphere_points).T.astype(np.float64)  # 3 x 32
    sphere_points = np.asfortranarray(sphere_points)
    
    print(f"Sphere points shape: {sphere_points.shape}")
    print(f"Sample point norms: {np.linalg.norm(sphere_points[:, :3], axis=0)}")
    
    sphere_ix = build_spherical_kd(sphere_points, np.array([1, 2, 3], dtype=np.int32))
    print(f"Spherical indices for sphere points: {sphere_ix[:10]}...")  # Show first 10
    
    # Test case 2: Hemisphere (common in many applications)
    print("\nTest 2: Hemisphere points")
    hemisphere_points = sphere_points[:, sphere_points[2, :] >= 0]  # Upper hemisphere
    hemisphere_points = np.asfortranarray(hemisphere_points)
    
    print(f"Hemisphere points shape: {hemisphere_points.shape}")
    hemisphere_ix = build_spherical_kd(hemisphere_points, np.array([1, 2, 3], dtype=np.int32))
    print(f"Spherical indices for hemisphere: {hemisphere_ix[:10]}...")

def test_bst_edge_cases():
    print("\n=== BST Edge Cases ===")
    # Empty array
    try:
        x = np.array([], dtype=np.float64)
        ix = build_bst_index(x)
        print("Empty array: No error (expected behavior)")
    except Exception as e:
        print(f"Empty array: Exception caught: {e}")

    # Single element
    try:
        x = np.array([42.0], dtype=np.float64)
        ix = build_bst_index(x)
        print(f"Single element BST indices: {ix}")
    except Exception as e:
        print(f"Single element test failed: {e}")

def test_kdtree_edge_cases():
    print("\n=== KD-Tree Edge Cases ===")
    # Empty matrix
    try:
        X = np.empty((2, 0), dtype=np.float64, order='F')
        kd_ix = build_kd_index(X)
        print("Empty matrix: No error (expected behavior)")
    except Exception as e:
        print(f"Empty matrix: Exception caught: {e}")

    # Single point
    try:
        X = np.array([[1.0], [2.0]], dtype=np.float64, order='F')
        kd_ix = build_kd_index(X, np.array([1, 2], dtype=np.int32))
        print(f"Single point KD-Tree indices: {kd_ix}")
    except Exception as e:
        print(f"Single point test failed: {e}")

def test_spherical_kdtree_edge_cases():
    print("\n=== Spherical KD-Tree Edge Cases ===")
    
    # Empty spherical data
    try:
        empty_vectors = np.empty((3, 0), dtype=np.float64, order='F')
        sphere_ix = build_spherical_kd(empty_vectors)
        print("Empty spherical data: No error (expected behavior)")
    except Exception as e:
        print(f"Empty spherical data: Exception caught: {e}")
    
    # Single vector on sphere
    try:
        single_vector = np.array([[0.0], [0.0], [1.0]], dtype=np.float64, order='F')  # North pole
        sphere_ix = build_spherical_kd(single_vector, np.array([1, 2, 3], dtype=np.int32))
        print(f"Single vector spherical indices: {sphere_ix}")
    except Exception as e:
        print(f"Single vector test failed: {e}")
    
    # 2D case (not truly spherical but should work)
    try:
        circle_vectors = np.array([[1.0, 0.0, -1.0, 0.0], 
                                 [0.0, 1.0, 0.0, -1.0]], dtype=np.float64, order='F')
        circle_vectors = circle_vectors / np.linalg.norm(circle_vectors, axis=0)  # Normalize
        circle_ix = build_spherical_kd(circle_vectors, np.array([1, 2], dtype=np.int32))
        print(f"2D circle indices: {circle_ix}")
    except Exception as e:
        print(f"2D circle test failed: {e}")

def test_performance():
    print("\n=== Performance Testing ===")
    
    # Test with larger datasets
    print("Testing with larger datasets...")
    
    # Larger BST
    large_values = np.random.rand(1000).astype(np.float64)
    bst_indices = build_bst_index(large_values)
    print(f"Large BST built successfully for {len(large_values)} values")
    
    # Larger spherical KD-Tree
    large_vectors = np.random.randn(3, 500).astype(np.float64)  # 3D, 500 vectors
    large_vectors = large_vectors / np.linalg.norm(large_vectors, axis=0)
    large_vectors = np.asfortranarray(large_vectors)
    
    sphere_indices = build_spherical_kd(large_vectors)
    print(f"Large spherical KD-Tree built successfully for {large_vectors.shape[1]} vectors")

def run_all_tests():
    """Run all test functions"""
    print("Running all tree structure tests...\n")
    
    test_bst()
    test_kdtree()
    test_spherical_kdtree()
    test_spherical_kdtree_specific_cases()
    test_bst_edge_cases()
    test_kdtree_edge_cases()
    test_spherical_kdtree_edge_cases()
    test_performance()
    
    print("\n=== All tests completed ===")

if __name__ == "__main__":
    run_all_tests()