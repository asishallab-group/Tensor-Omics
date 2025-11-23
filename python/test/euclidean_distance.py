#!/usr/bin/env python3
"""
Comprehensive Python test suite for Euclidean distance functions
Uses tensoromics_functions.py wrapper functions (mirrors R euclidean_distance.R tests)
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_euclidean_distance, tox_distance_to_centroid

print("=== Testing euclidean distance Python wrapper functions ===")
print("Based on Fortran test suite with comprehensive test coverage")

# =====================
# Tests for euclidean_distance
# =====================

def test_euclidean_distance_3d():
    """Test simple 3D vectors"""
    print("\n[test_euclidean_distance_3d] Simple 3D vectors test")
    vec1 = np.array([1.0, 2.0, 3.0])
    vec2 = np.array([4.0, 5.0, 6.0])
    
    result = tox_euclidean_distance(vec1, vec2)
    expected = np.sqrt(np.sum((vec1 - vec2)**2))
    
    print("  Vectors [1,2,3] vs [4,5,6]")
    print(f"  Result: {result}")
    print(f"  Expected: {expected}")
    
    # Verify calculation
    assert abs(result - expected) < 1e-12
    
    print("Simple 3D vectors test passed ✓")

def test_euclidean_distance_to_origin():
    """Test 2D vector to origin (3-4-5 triangle)"""
    print("\n[test_euclidean_distance_to_origin] 2D vector to origin test")
    vec1 = np.array([3.0, 4.0])
    vec2 = np.array([0.0, 0.0])
    
    result = tox_euclidean_distance(vec1, vec2)
    expected = 5.0
    
    print("  Vector [3,4] to origin")
    print(f"  Result: {result}")
    print(f"  Expected: {expected}")
    
    # Verify 3-4-5 triangle
    assert abs(result - expected) < 1e-12
    
    print("2D vector to origin test passed ✓")

def test_euclidean_distance_identical():
    """Test identical vectors"""
    print("\n[test_euclidean_distance_identical] Identical vectors test")
    vec1 = np.array([1.5, 2.7, 3.9])
    vec2 = vec1.copy()  # identical
    
    result = tox_euclidean_distance(vec1, vec2)
    expected = 0.0
    
    print("  Identical vectors")
    print(f"  Result: {result}")
    print(f"  Expected: {expected}")
    
    # Verify zero distance
    assert abs(result - expected) < 1e-15
    
    print("Identical vectors test passed ✓")

def test_euclidean_distance_high_dimensional():
    """Test high-dimensional vectors"""
    print("\n[test_euclidean_distance_high_dimensional] High-dimensional test")
    # 100-dimensional vectors
    d = 100
    vec1 = np.arange(1, d+1, dtype=np.float64)  # [1, 2, 3, ..., 100]
    vec2 = np.arange(2, d+2, dtype=np.float64)  # [2, 3, 4, ..., 101] (shift by 1)
    
    result = tox_euclidean_distance(vec1, vec2)
    expected = np.sqrt(d)  # sqrt(100 * 1^2) = 10
    
    print("  100D vectors with unit shift")
    print(f"  Result: {result}")
    print(f"  Expected: {expected}")
    
    # Verify high-dimensional calculation
    assert abs(result - expected) < 1e-12
    
    print("High-dimensional test passed ✓")

def test_euclidean_distance_invalid_inputs():
    """Test invalid inputs (should throw errors)"""
    print("\n[test_euclidean_distance_invalid_inputs] Invalid inputs test")
    
    # Test different lengths
    error_caught1 = False
    try:
        tox_euclidean_distance(np.array([1, 2]), np.array([1, 2, 3]))
    except ValueError as e:
        error_caught1 = True
        assert "same length" in str(e)
    assert error_caught1
    
    # Test empty vectors
    error_caught2 = False
    try:
        tox_euclidean_distance(np.array([]), np.array([]))
    except ValueError as e:
        error_caught2 = True
        assert "cannot be empty" in str(e)
    assert error_caught2
    
    # Test non-numeric input
    error_caught3 = False
    try:
        tox_euclidean_distance(np.array(["a", "b"]), np.array([1, 2]))
    except (ValueError, TypeError) as e:
        error_caught3 = True
        # Could be either ValueError from wrapper or TypeError from numpy conversion
    assert error_caught3
    
    print("Invalid inputs test passed ✓")

def test_euclidean_distance_1d():
    """Test single-dimensional vectors"""
    print("\n[test_euclidean_distance_1d] Single-dimensional vectors test")
    vec1 = np.array([5.0])
    vec2 = np.array([2.0])
    
    result = tox_euclidean_distance(vec1, vec2)
    expected = 3.0
    
    print("  1D vectors: 5.0 vs 2.0")
    print(f"  Result: {result}")
    print(f"  Expected: {expected}")
    
    # Verify 1D calculation
    assert abs(result - expected) < 1e-12
    
    print("Single-dimensional vectors test passed ✓")

# =====================
# Tests for distance_to_centroid  
# =====================

def test_distance_to_centroid_basic():
    """Test basic distance to centroid functionality"""
    print("\n[test_distance_to_centroid_basic] Basic distance to centroid test")
    
    # Setup test data
    d = 3
    
    # Gene expression data (genes as columns, dimensions as rows)
    # Gene 1: [1, 0, 0] - Family 1
    # Gene 2: [0, 1, 0] - Family 1  
    # Gene 3: [3, 0, 0] - Family 2
    # Gene 4: [0, 3, 0] - Family 2
    genes = np.array([1.0, 0.0, 0.0,  # Gene 1
                      0.0, 1.0, 0.0,  # Gene 2
                      3.0, 0.0, 0.0,  # Gene 3
                      0.0, 3.0, 0.0])  # Gene 4
    
    # Family centroids
    # Family 1 centroid: [0.5, 0.5, 0.0]
    # Family 2 centroid: [1.5, 1.5, 0.0]
    centroids = np.array([0.5, 0.5, 0.0,  # Family 1
                          1.5, 1.5, 0.0])  # Family 2
    
    # Gene-to-family mapping (1-based)
    gene_to_fam = np.array([1, 1, 2, 2], dtype=np.int32)
    
    result = tox_distance_to_centroid(genes, centroids, gene_to_fam, d)
    
    # Expected distances
    # Gene 1: [1,0,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
    # Gene 2: [0,1,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
    # Gene 3: [3,0,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
    # Gene 4: [0,3,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
    expected = np.array([np.sqrt(0.5**2 + 0.5**2), np.sqrt(0.5**2 + 0.5**2), 
                         np.sqrt(1.5**2 + 1.5**2), np.sqrt(1.5**2 + 1.5**2)])
    
    print("  Distance to centroid results:")
    for i in range(len(result)):
        print(f"    Gene {i+1}: Result={result[i]:.6f}, Expected={expected[i]:.6f}")
        assert abs(result[i] - expected[i]) < 1e-12
    
    print("Basic distance to centroid test passed ✓")

def test_distance_to_centroid_invalid_families():
    """Test handling invalid family indices (should return -1 for invalid genes)"""
    print("\n[test_distance_to_centroid_invalid_families] Invalid family indices test")
    
    d = 2
    
    # Gene data
    genes = np.array([1.0, 2.0,  # Gene 1
                      3.0, 4.0,  # Gene 2
                      5.0, 6.0])  # Gene 3
    
    # Centroids
    centroids = np.array([0.0, 0.0,  # Family 1
                          1.0, 1.0])  # Family 2
    
    # Mixed family assignments: valid (1), invalid (3), no family (0)
    gene_to_fam = np.array([1, 3, 0], dtype=np.int32)  # Family 3 doesn't exist, family 0 = no assignment
    
    result = tox_distance_to_centroid(genes, centroids, gene_to_fam, d)
    
    print(f"  Gene 1 (family 1): {result[0]} (should be valid distance)")
    print(f"  Gene 2 (family 3): {result[1]} (should be -1, invalid family)")
    print(f"  Gene 3 (family 0): {result[2]} (should be -1, no family)")
    
    # Verify handling of invalid indices
    assert result[0] > 0      # Gene 1 has valid family, should have positive distance
    assert result[1] == -1    # Gene 2 has invalid family, should be -1
    assert result[2] == -1    # Gene 3 has no family, should be -1
    
    print("Invalid family indices test passed ✓")

def test_distance_to_centroid_performance():
    """Test performance with realistic genomic data size"""
    print("\n[test_distance_to_centroid_performance] Performance test")
    
    n_genes = 1000
    n_families = 50
    d = 20  # 20 tissue types
    
    print(f"  Testing with {n_genes} genes, {n_families} families, {d} dimensions")
    
    # Generate random-like data
    np.random.seed(12345)
    genes = np.random.randn(d * n_genes)
    centroids = np.random.randn(d * n_families)
    gene_to_fam = np.random.randint(1, n_families + 1, n_genes, dtype=np.int32)
    
    # Time the operation
    import time
    start_time = time.time()
    
    result = tox_distance_to_centroid(genes, centroids, gene_to_fam, d)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    # Check results
    valid_distances = np.sum(result > 0)
    print(f"  Completed in {elapsed:.6f} seconds")
    print(f"  Valid distances computed: {valid_distances}/{n_genes}")
    print(f"  Mean distance: {np.mean(result):.6f}")
    
    # Verify all distances are positive
    assert np.all(result > 0)
    assert len(result) == n_genes
    
    print("Performance test passed ✓")

def test_distance_to_centroid_input_validation():
    """Test input validation for distance_to_centroid"""
    print("\n[test_distance_to_centroid_input_validation] Input validation test")
    
    # Test dimension mismatch
    error_caught1 = False
    try:
        tox_distance_to_centroid(np.array([1, 2, 3]), np.array([1, 2]), np.array([1]), 2)  # genes not divisible by d
    except ValueError as e:
        error_caught1 = True
        assert "divisible by d" in str(e)
    assert error_caught1
    
    # Test gene_to_fam length mismatch
    error_caught2 = False
    try:
        tox_distance_to_centroid(np.array([1, 2, 3, 4]), np.array([1, 2]), np.array([1, 2, 3]), 2)  # wrong gene_to_fam length
    except ValueError as e:
        error_caught2 = True
        assert "equal number of genes" in str(e)
    assert error_caught2
    
    # Test negative family indices (should throw error)
    error_caught3 = False
    try:
        tox_distance_to_centroid(np.array([1, 2, 3, 4]), np.array([1, 2]), np.array([1, -1]), 2)  # negative family index
    except ValueError as e:
        error_caught3 = True
        assert "must be between 0 and" in str(e)
    assert error_caught3
    
    print("Input validation test passed ✓")


# =====================
# Run all tests
# =====================

if __name__ == "__main__":
    print("\n=================================================")
    print("    EUCLIDEAN DISTANCE FULL PYTHON INTERFACE TESTS")
    print("=================================================\n")

    # euclidean_distance tests
    test_euclidean_distance_3d()
    test_euclidean_distance_to_origin()
    test_euclidean_distance_identical()
    test_euclidean_distance_high_dimensional()
    test_euclidean_distance_invalid_inputs()
    test_euclidean_distance_1d()

    # distance_to_centroid tests
    test_distance_to_centroid_basic()
    test_distance_to_centroid_invalid_families()
    test_distance_to_centroid_performance()
    test_distance_to_centroid_input_validation()

    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all euclidean distance Python interface tests passed! ✓")
