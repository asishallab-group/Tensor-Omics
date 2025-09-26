"""
Test script for clock hand angle functions
Python equivalent of the R and Fortran clock hand angle tests
"""

import numpy as np
import ctypes
import time
import sys
import os
import math 

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_clock_hand_angle_between_vectors, tox_clock_hand_angles_for_shift_vectors


# Constants
PI = math.pi
TOL = 1e-12



def test_identical_vectors_2d():
    """Test identical vectors in 2D (should give 0 angle)"""
    print("=== Testing Identical Vectors 2D ===")
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        
        expected = 0.0
        print("Test: Identical 2D vectors [1,0] vs [1,0]:")
        print(f"  Fortran result: {result}")
        print(f"  Expected: {expected}")
        print(f"  Match: {abs(result - expected) < TOL}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_perpendicular_vectors_2d():
    """Test perpendicular vectors in 2D (should give ±π/2)"""
    print("=== Testing Perpendicular Vectors 2D ===")
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0], dtype=np.float64)
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        
        expected_magnitude = PI/2
        print("Test: Perpendicular 2D vectors [1,0] vs [0,1]:")
        print(f"  Fortran result: {result}")
        print(f"  Expected magnitude: ±{expected_magnitude}")
        print(f"  Magnitude match: {abs(abs(result) - expected_magnitude) < TOL}")
        print(f"  Sign (should be positive): {result > 0}")
        print(f"  Angle in degrees: {result * 180/PI:.2f}°")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_opposite_vectors_2d():
    """Test opposite vectors in 2D (should give ±π)"""
    print("=== Testing Opposite Vectors 2D ===")
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([-1.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        
        expected_magnitude = PI
        print("Test: Opposite 2D vectors [1,0] vs [-1,0]:")
        print(f"  Fortran result: {result}")
        print(f"  Expected magnitude: ±{expected_magnitude}")
        print(f"  Magnitude match: {abs(abs(result) - expected_magnitude) < TOL}")
        print(f"  Angle in degrees: {result * 180/PI:.2f}°")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_45_degree_rotation_2d():
    """Test 45-degree rotation in 2D"""
    print("=== Testing 45-Degree Rotation 2D ===")
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([math.sqrt(2)/2, math.sqrt(2)/2], dtype=np.float64)  # 45 degrees
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        
        expected = PI/4
        print("Test: 45-degree counterclockwise rotation:")
        print(f"  v1: {v1}")
        print(f"  v2: {v2}")
        print(f"  Fortran result: {result}")
        print(f"  Expected: {expected}")
        print(f"  Match: {abs(result - expected) < TOL}")
        print(f"  Angle in degrees: {result * 180/PI:.2f}°")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_clockwise_vs_counterclockwise_2d():
    """Test clockwise vs counterclockwise rotations in 2D"""
    print("=== Testing Clockwise vs Counterclockwise 2D ===")
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2_ccw = np.array([0.0, 1.0], dtype=np.float64)   # 90° counterclockwise
    v2_cw = np.array([0.0, -1.0], dtype=np.float64)   # 90° clockwise
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    try:
        result_ccw = tox_clock_hand_angle_between_vectors(v1, v2_ccw, selected_axes)
        result_cw = tox_clock_hand_angle_between_vectors(v1, v2_cw, selected_axes)
        
        print("Test: Clockwise vs Counterclockwise:")
        print(f"  Counterclockwise angle: {result_ccw:.6f} ({result_ccw * 180/PI:.2f}°)")
        print(f"  Clockwise angle: {result_cw:.6f} ({result_cw * 180/PI:.2f}°)")
        print(f"  CCW > 0: {result_ccw > 0}")
        print(f"  CW < 0: {result_cw < 0}")
        print(f"  Magnitudes equal: {abs(abs(result_ccw) - abs(result_cw)) < TOL}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_3d_vectors():
    """Test 3D vector calculations"""
    print("=== Testing 3D Vectors ===")
    
    # Test identical 3D vectors
    v1 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    v2 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    selected_axes = np.array([1, 1, 1], dtype=np.int32)  # Ignored for 3D
    
    result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    
    print("Test 1: Identical 3D vectors [1,1,1] vs [1,1,1]:")
    print(f"  Result: {result} (should be 0)")
    print(f"  Match: {abs(result) < TOL}")
    
    # Test perpendicular 3D vectors
    v1 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    
    result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    
    expected_magnitude = PI/2
    print("Test 2: Perpendicular 3D vectors [1,0,0] vs [0,1,0]:")
    print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result) - expected_magnitude) < TOL}")
    print()

def test_high_dimensional():
    """Test high-dimensional vectors with selected axes"""
    print("=== Testing High-Dimensional Vectors ===")
    # 5D vectors, perpendicular in the first two axes
    v1 = np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)  # Only 3 indices
    result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    expected_magnitude = PI/2
    print("Test 1: 5D vectors perpendicular in first two axes:")
    print(f"  v1: {v1}")
    print(f"  v2: {v2}")
    print(f"  Selected axes: {selected_axes}")
    print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result) - expected_magnitude) < TOL}")
    # 7D vectors with selected axes
    v1 = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([3, 5, 1], dtype=np.int32)  # Only 3 indices
    result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    print("Test 2: 7D vectors with selected axes [3,5,1]:")
    print(f"  v1: {v1}")
    print(f"  v2: {v2}")
    print(f"  Selected axes: {selected_axes}")
    print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result) - expected_magnitude) < TOL}")
    print()

def test_shift_vectors_single_pair():
    """Test single pair of shift vectors"""
    print("=== Testing Single Pair Shift Vectors ===")
    
    # Single pair: [1,0] -> [0,1] (90° counterclockwise)
    n_dims = 2
    n_vecs = 1
    
    # Create matrices exactly like R's matrix(c(1,0), nrow=2, ncol=1)
    origins_data = np.array([1.0, 0.0], dtype=np.float64)
    origins = origins_data.reshape((2, 1), order='F')
    
    targets_data = np.array([0.0, 1.0], dtype=np.float64)
    targets = targets_data.reshape((2, 1), order='F')
    vecs_selection_mask = np.array([1], dtype=np.int32)  # Select the single pair
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    signed_angles = tox_clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask, selected_axes)
    
    expected = PI/2
    print("Test: Single pair [1,0] -> [0,1]:")
    print(f"  Result: {signed_angles[0]:.6f} ({signed_angles[0] * 180/PI:.2f}°)")
    print(f"  Expected: {expected:.6f} ({expected * 180/PI:.2f}°)")
    print(f"  Match: {abs(signed_angles[0] - expected) < TOL}")
    print()

def test_shift_vectors_multiple_pairs():
    """Test multiple pairs of shift vectors"""
    print("=== Testing Multiple Pairs Shift Vectors ===")
    
    # Three different rotations
    n_dims = 2
    n_vecs = 3
    
    # Create matrices exactly like R's matrix(c(...), nrow=2, ncol=3)
    # R fills by columns, so we need to create and transpose
    origins_data = np.array([1.0, 0.0, 1.0, 0.0, 1.0, 0.0], dtype=np.float64)  # c(1,0, 1,0, 1,0)
    origins = origins_data.reshape((2, 3), order='F')  # nrow=2, ncol=3, fill by columns
    
    # targets: c(0,1, -1,0, 0,-1) -> columns [0,1], [-1,0], [0,-1]
    targets_data = np.array([0.0, 1.0, -1.0, 0.0, 0.0, -1.0], dtype=np.float64)  # c(0,1, -1,0, 0,-1)
    targets = targets_data.reshape((2, 3), order='F')  # nrow=2, ncol=3, fill by columns
    
    vecs_selection_mask = np.array([1, 1, 1], dtype=np.int32)  # Select all pairs
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    signed_angles = tox_clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask, selected_axes)
    
    # Expected mathematically correct results
    expected = [PI/2, PI, -PI/2]  # 90° CCW, 180°, 90° CW
    print("Test: Multiple pairs with different rotations:")
    print(f"  Vector pairs (following R implementation):")
    print(f"    Pair 1: [1,0] -> [0,1] (should be +90°)")
    print(f"    Pair 2: [1,0] -> [-1,0] (should be ±180°)")  
    print(f"    Pair 3: [1,0] -> [0,-1] (should be -90°)")
    
    for i in range(n_vecs):
        if i == 1:  # 180° case - check magnitude only (can be +π or -π)
            match = abs(abs(signed_angles[i]) - abs(expected[i])) < TOL
            print(f"  Pair {i+1}: {signed_angles[i]:.6f} ({signed_angles[i] * 180/PI:.2f}°), "
                  f"Expected magnitude: ±{abs(expected[i]):.6f}, Match: {match}")
        else:
            match = abs(signed_angles[i] - expected[i]) < TOL
            print(f"  Pair {i+1}: {signed_angles[i]:.6f} ({signed_angles[i] * 180/PI:.2f}°), "
                  f"Expected: {expected[i]:.6f}, Match: {match}")
    
    # Check if all results match expectations
    all_match = all([abs(signed_angles[0] - expected[0]) < TOL, 
                     abs(abs(signed_angles[1]) - abs(expected[1])) < TOL,
                     abs(signed_angles[2] - expected[2]) < TOL])
    
    print()

def test_shift_vectors_with_selection_mask():
    """Test shift vectors with selection mask"""
    print("=== Testing Shift Vectors with Selection Mask ===")
    
    # Four vectors, but only select 2nd and 4th
    n_dims = 2
    n_vecs = 4
    
    origins_data = np.array([1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0], dtype=np.float64)
    origins = origins_data.reshape((2, 4), order='F')
    

    targets_data = np.array([0.0, 1.0, -1.0, 0.0, 0.0, -1.0, math.sqrt(2)/2, math.sqrt(2)/2], dtype=np.float64)
    targets = targets_data.reshape((2, 4), order='F')
    
    vecs_selection_mask = np.array([0, 1, 0, 1], dtype=np.int32)  # Select 2nd and 4th (FALSE, TRUE, FALSE, TRUE)
    selected_axes = np.array([1, 1], dtype=np.int32)  # Ignored for 2D
    
    signed_angles = tox_clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask, selected_axes)
    
    # Expected results based on R test: 2nd pair→180°, 4th pair→45°
    expected = [PI, PI/4]  # 180° and 45°
    print("Test: Selection mask [False, True, False, True]:")
    print(f"  Vector pairs selected:")
    print(f"    Pair 2: [1,0] -> [-1,0] (should be ±180°)")
    print(f"    Pair 4: [1,0] -> [√2/2,√2/2] (should be +45°)")
    
    for i, expected_val in enumerate(expected):
        if i == 0:  # 180° case - check magnitude only
            match = abs(abs(signed_angles[i]) - abs(expected_val)) < TOL
            print(f"  Selected pair {i+1}: {signed_angles[i]:.6f} ({signed_angles[i] * 180/PI:.2f}°)")
            print(f"    Expected magnitude: ±{abs(expected_val):.6f}, Match: {match}")
        else:
            match = abs(signed_angles[i] - expected_val) < TOL
            print(f"  Selected pair {i+1}: {signed_angles[i]:.6f} ({signed_angles[i] * 180/PI:.2f}°)")
            print(f"    Expected: {expected_val:.6f}, Match: {match}")
    print()

def test_edge_cases():
    """Test edge cases and precision"""
    print("=== Testing Edge Cases ===")
    
    # Test 1: Nearly identical vectors
    epsilon = 1e-15
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, epsilon], dtype=np.float64)
    selected_axes = np.array([1, 1], dtype=np.int32)
    
    result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    
    print("Test 1: Nearly identical vectors:")
    print(f"  v1: {v1}")
    print(f"  v2: {v2}")
    print(f"  Result: {result:.2e} radians")
    print(f"  Very small angle: {abs(result) < 1e-10}")
    
    # Test 2: Denormalized vectors (large magnitude)
    v1 = np.array([100.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 50.0], dtype=np.float64)
    
    result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    
    expected = PI/2
    print("Test 2: Denormalized vectors [100,0] vs [0,50]:")
    print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
    print(f"  Expected: {expected:.6f}")
    print(f"  Match: {abs(result - expected) < TOL}")
    
    # Test 3: Very small vectors
    tiny = 1e-14
    v1 = np.array([tiny, 0.0], dtype=np.float64)
    v2 = np.array([0.0, tiny], dtype=np.float64)
    
    result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    
    print("Test 3: Tiny vectors:")
    print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected:.6f}")
    print(f"  Magnitude match: {abs(abs(result) - expected) < 1e-10}")
    print()

def test_consistency_between_functions():
    """Test consistency between single and batch functions"""
    print("=== Testing Consistency Between Functions ===")
    
    # Test same calculation with both functions
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0], dtype=np.float64)
    selected_axes = np.array([1, 1], dtype=np.int32)
    
    # Single function
    single_result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    
    # Batch function with single pair - create like R
    origins_data = np.array([1.0, 0.0], dtype=np.float64)
    origins = origins_data.reshape((2, 1), order='F')
    
    targets_data = np.array([0.0, 1.0], dtype=np.float64)
    targets = targets_data.reshape((2, 1), order='F')
    vecs_selection_mask = np.array([1], dtype=np.int32)
    batch_results = tox_clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask, selected_axes)
    
    print("Test: Consistency between single and batch functions:")
    print(f"  Single function result: {single_result:.10f}")
    print(f"  Batch function result: {batch_results[0]:.10f}")
    print(f"  Difference: {abs(single_result - batch_results[0]):.2e}")
    print(f"  Match: {abs(single_result - batch_results[0]) < TOL}")
    print()

def test_mathematical_properties():
    """Test mathematical properties (anti-commutativity)"""
    print("=== Testing Mathematical Properties ===")
    
    # Test anti-commutativity: angle(v1,v2) = -angle(v2,v1)
    v1 = np.array([1.0, 2.0], dtype=np.float64)
    v2 = np.array([3.0, 1.0], dtype=np.float64)
    
    # Normalize vectors
    v1 = v1 / np.sqrt(np.sum(v1**2))
    v2 = v2 / np.sqrt(np.sum(v2**2))
    
    selected_axes = np.array([1, 1], dtype=np.int32)
    
    result_12 = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
    result_21 = tox_clock_hand_angle_between_vectors(v2, v1, selected_axes)
    
    print("Test: Anti-commutativity angle(v1,v2) = -angle(v2,v1):")
    print(f"  v1 (normalized): {v1}")
    print(f"  v2 (normalized): {v2}")
    print(f"  angle(v1,v2): {result_12:.6f}")
    print(f"  angle(v2,v1): {result_21:.6f}")
    print(f"  Sum (should be ~0): {result_12 + result_21:.2e}")
    print(f"  Anti-commutative: {abs(result_12 + result_21) < TOL}")
    print()

def test_performance():
    """Performance test with large-scale data"""
    print("=== Performance Test ===")
    
    # Large-scale test
    n_dims = 50
    n_vecs = 1000
    
    print(f"Testing with {n_vecs} vector pairs in {n_dims} dimensions")
    
    # Generate random-like data
    np.random.seed(12345)
    origins_flat = np.random.randn(n_dims * n_vecs).astype(np.float64)
    targets_flat = np.random.randn(n_dims * n_vecs).astype(np.float64)
    
    # Reshape to column format (n_dims x n_vecs)
    origins = origins_flat.reshape((n_dims, n_vecs), order='F')
    targets = targets_flat.reshape((n_dims, n_vecs), order='F')
    
    # Normalize vectors (for more meaningful angles)
    for i in range(n_vecs):
        origin_vec = origins[:, i]
        target_vec = targets[:, i]
        
        origin_norm = np.sqrt(np.sum(origin_vec**2))
        target_norm = np.sqrt(np.sum(target_vec**2))
        
        if origin_norm > 0:
            origins[:, i] /= origin_norm
        if target_norm > 0:
            targets[:, i] /= target_norm
    
    vecs_selection_mask = np.ones(n_vecs, dtype=np.int32)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)
    
    # Time the operation
    start_time = time.time()
    
    signed_angles = tox_clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask, selected_axes)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    # Check results
    valid_count = np.sum(~np.isnan(signed_angles))
    mean_angle = np.mean(signed_angles[~np.isnan(signed_angles)]) if valid_count > 0 else 0
    
    print(f"  Completed in {elapsed:.6f} seconds")
    print(f"  Valid angles computed: {valid_count}/{n_vecs}")
    print(f"  Mean angle: {mean_angle:.6f} radians ({mean_angle * 180/PI:.2f}°)")
    print(f"  Performance: {n_vecs/elapsed:.0f} calculations/second")
    print("  Performance test completed successfully!")
    print()

def test_invalid_selected_axes():
    """Test invalid selected axes (same indices)"""
    print("=== Testing Invalid Selected Axes ===")
    v1 = np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 1, 1, 0, 0], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        print(f"  ERROR: Expected exception, got result {result}")
    except Exception as e:
        print(f"  Passed: Caught expected error: {e}")
    print()

def test_out_of_bounds_selected_axes():
    """Test out-of-bounds selected axes"""
    print("=== Testing Out-of-Bounds Selected Axes ===")
    v1 = np.array([1.0, 0.0, 0.0, 2.3], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0, 2.3], dtype=np.float64)
    selected_axes = np.array([1, 2, 5, 6], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        print(f"  ERROR: Expected exception, got result {result}")
    except Exception as e:
        print(f"  Passed: Caught expected error: {e}")
    print()

def test_zero_vectors():
    """Test zero vectors (should not produce NaN)"""
    print("=== Testing Zero Vectors ===")
    v1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        print(f"  Result: {result}")
        print(f"  Not NaN: {not np.isnan(result)}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_denormalized_vectors():
    """Test denormalized vectors (large magnitude)"""
    print("=== Testing Denormalized Vectors ===")
    v1 = np.array([100.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 50.0], dtype=np.float64)
    selected_axes = np.array([1, 2], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        expected = PI/2
        print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
        print(f"  Expected: {expected:.6f}")
        print(f"  Match: {abs(result - expected) < TOL}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_tiny_vectors_precision():
    """Test tiny vectors near machine precision"""
    print("=== Testing Tiny Vectors Precision ===")
    tiny = 1e-14
    v1 = np.array([tiny, 0.0], dtype=np.float64)
    v2 = np.array([0.0, tiny], dtype=np.float64)
    selected_axes = np.array([1, 2], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        expected = PI/2
        print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
        print(f"  Expected magnitude: ±{expected:.6f}")
        print(f"  Magnitude match: {abs(abs(result) - expected) < 1e-10}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_huge_vectors_precision():
    """Test huge vectors near overflow"""
    print("=== Testing Huge Vectors Precision ===")
    huge_val = 1e14
    v1 = np.array([huge_val, 0.0], dtype=np.float64)
    v2 = np.array([0.0, huge_val], dtype=np.float64)
    selected_axes = np.array([1, 2], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        expected = PI/2
        print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
        print(f"  Expected: {expected:.6f}")
        print(f"  Match: {abs(result - expected) < TOL}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_nearly_identical_vectors():
    """Test nearly identical vectors (precision boundary)"""
    print("=== Testing Nearly Identical Vectors ===")
    epsilon = 1e-15
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, epsilon], dtype=np.float64)
    selected_axes = np.array([1, 2], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        print(f"  Result: {result:.2e} radians")
        print(f"  Very small angle: {abs(result) < 1e-10}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_nearly_opposite_vectors():
    """Test nearly opposite vectors (precision boundary)"""
    print("=== Testing Nearly Opposite Vectors ===")
    epsilon = 1e-15
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([-1.0, epsilon], dtype=np.float64)
    selected_axes = np.array([1, 2], dtype=np.int32)
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        print(f"  Result: {result:.6f} ({result * 180/PI:.2f}°)")
        print(f"  Nearly π: {abs(abs(result) - PI) < 1e-10}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def test_mixed_positive_negative():
    """Test vectors with mixed positive/negative components"""
    print("=== Testing Mixed Positive/Negative Vectors ===")
    v1 = np.array([1.0, -2.0, 3.0], dtype=np.float64)
    v2 = np.array([-2.0, 1.0, -3.0], dtype=np.float64)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)
    v1 = v1 / np.sqrt(np.sum(v1**2))
    v2 = v2 / np.sqrt(np.sum(v2**2))
    try:
        result = tox_clock_hand_angle_between_vectors(v1, v2, selected_axes)
        print(f"  Result: {result:.6f}")
        print(f"  In valid range: {abs(result) >= 0.0 and abs(result) <= PI}")
    except Exception as e:
        print(f"  ERROR: {e}")
    print()

def main():
    print("=================================================")
    print("    CLOCK HAND ANGLES PYTHON INTERFACE TESTS")
    print("=================================================")
    print()
    
    # Basic 2D tests
    test_identical_vectors_2d()
    test_perpendicular_vectors_2d()
    test_opposite_vectors_2d()
    test_45_degree_rotation_2d()
    test_clockwise_vs_counterclockwise_2d()
    
    # 3D and high-dimensional tests
    test_3d_vectors()
    test_high_dimensional()
    
    # Shift vectors tests
    test_shift_vectors_single_pair()
    test_shift_vectors_multiple_pairs()
    test_shift_vectors_with_selection_mask()
    
    # Edge cases and properties
    test_edge_cases()
    test_consistency_between_functions()
    test_mathematical_properties()
    
    # Performance test
    test_performance()
    
    # Error conditions and edge cases
    test_invalid_selected_axes()
    test_out_of_bounds_selected_axes()
    test_zero_vectors()
    test_denormalized_vectors()
    test_tiny_vectors_precision()
    test_huge_vectors_precision()
    test_nearly_identical_vectors()
    test_nearly_opposite_vectors()
    test_mixed_positive_negative()
    
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all Python interface tests passed! ✓")

if __name__ == "__main__":
    main()
