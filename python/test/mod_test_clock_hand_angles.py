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

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import clock_hand_angle_between_vectors, clock_hand_angles_for_shift_vectors
from test_helpers import *


# Constants
PI = math.pi
TOL = 1e-12
# selected_axes_for_signed is ignored for n_dims <= 3, so any valid triple will do there
ANY_AXES = np.array([1, 2, 3], dtype=np.int32)


def test_identical_vectors_2d():
    """Test identical vectors in 2D (should give 0 angle)"""

    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, 0.0], dtype=np.float64)

    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)

    expected = 0.0
    assert np.isclose(result, expected), "expected zero angle"


def test_perpendicular_vectors_2d():
    """Test perpendicular vectors in 2D (should give ±π/2)"""

    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0], dtype=np.float64)

    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)

    expected_magnitude = PI/2
    assert np.isclose(result, expected_magnitude), "expected PI/2 angle"


def test_opposite_vectors_2d():
    """Test opposite vectors in 2D (should give ±π)"""

    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([-1.0, 0.0], dtype=np.float64)

    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)

    expected_magnitude = PI
    assert np.isclose(result, expected_magnitude), "expected PI angle"


def test_45_degree_rotation_2d():
    """Test 45-degree rotation in 2D"""

    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([math.sqrt(2)/2, math.sqrt(2)/2], dtype=np.float64)  # 45 degrees

    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)

    expected = PI/4
    assert np.isclose(result, expected), "expected PI/4 angle"


def test_clockwise_vs_counterclockwise_2d():
    """Test clockwise vs counterclockwise rotations in 2D"""

    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2_ccw = np.array([0.0, 1.0], dtype=np.float64)   # 90° counterclockwise
    v2_cw = np.array([0.0, -1.0], dtype=np.float64)   # 90° clockwise

    result_ccw = clock_hand_angle_between_vectors(v1, v2_ccw, ANY_AXES)
    result_cw = clock_hand_angle_between_vectors(v1, v2_cw, ANY_AXES)

    assert np.isclose(result_ccw, PI/2), "ccw expected PI/2 angle"
    assert np.isclose(result_cw, -PI/2), "cw expected PI/2 angle"


def test_3d_vectors():
    """Test 3D vector calculations"""

    # Test identical 3D vectors
    v1 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    v2 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    selected_axes = np.array([1, 1, 1], dtype=np.int32)  # Ignored for 3D

    result = clock_hand_angle_between_vectors(v1, v2, selected_axes)

    assert np.isclose(result, 0), "expected zero angle"

    # Test perpendicular 3D vectors
    v1 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0], dtype=np.float64)

    result = clock_hand_angle_between_vectors(v1, v2, selected_axes)

    expected_magnitude = PI/2
    assert np.isclose(result, expected_magnitude), "expected PI/2 angle"


def test_high_dimensional():
    """Test high-dimensional vectors with selected axes"""

    # 5D vectors, perpendicular in the first two axes
    v1 = np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)  # Only 3 indices
    result = clock_hand_angle_between_vectors(v1, v2, selected_axes)
    expected_magnitude = PI/2
    assert np.isclose(result, expected_magnitude), "5D expected PI/2 angle"

    # 7D vectors with selected axes
    v1 = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([3, 5, 1], dtype=np.int32)  # Only 3 indices
    result = clock_hand_angle_between_vectors(v1, v2, selected_axes)
    assert np.isclose(result, expected_magnitude), "7D expected PI/2 angle"


def test_shift_vectors_single_pair():
    """Test single pair of shift vectors"""

    # Single pair: [1,0] -> [0,1] (90° counterclockwise)
    n_dims = 2
    n_vecs = 1

    fields = np.array([[
        [1.0, 0.0],  # origin
        [0.0, 1.0]  # target
    ]], dtype=np.float64, order="F").transpose()
    vecs_selection_mask = np.array([1], dtype=np.int32)  # Select the single pair

    signed_angles = clock_hand_angles_for_shift_vectors(fields, vecs_selection_mask, ANY_AXES)

    expected = PI/2
    assert np.isclose(signed_angles, expected), "expected PI/2 angle"


def test_shift_vectors_multiple_pairs():
    """Test multiple pairs of shift vectors"""

    # Three different rotations
    n_dims = 2
    n_vecs = 3

    fields = np.array([
        [
            [1.0, 0.0],  # origin
            [0.0, 1.0]  # target
        ],
        [
            [1.0, 0.0],  # origin
            [-1.0, 0.0]  # target
        ],
        [
            [1.0, 0.0],  # origin
            [0.0, -1.0]  # target
        ]
    ], dtype=np.float64, order="F").transpose()
    vecs_selection_mask = np.array([1, 1, 1], dtype=np.int32)  # Select all pairs

    signed_angles = clock_hand_angles_for_shift_vectors(fields, vecs_selection_mask, ANY_AXES)

    # Expected mathematically correct results
    expected = [PI/2, PI, -PI/2]  # 90° CCW, 180°, 90° CW
    assert all(np.isclose(signed_angles, expected)), "multiple pairs mismatch"


def test_shift_vectors_with_selection_mask():
    """Test shift vectors with selection mask"""

    # Four vectors, but only select 2nd and 4th
    n_dims = 2
    n_vecs = 4

    fields = np.array([
        [
            [1.0, 0.0],  # origin
            [0.0, 1.0]  # target
        ],
        [
            [1.0, 0.0],  # origin
            [-1.0, 0.0]  # target
        ],
        [
            [1.0, 0.0],  # origin
            [0.0, -1.0]  # target
        ],
        [
            [1.0, 0.0],  # origin
            [math.sqrt(2)/2, math.sqrt(2)/2]  # target
        ]
    ], dtype=np.float64, order="F").transpose()

    vecs_selection_mask = np.array([0, 1, 0, 1], dtype=np.int32)  # Select 2nd and 4th (FALSE, TRUE, FALSE, TRUE)

    signed_angles = clock_hand_angles_for_shift_vectors(fields, vecs_selection_mask, ANY_AXES)

    # Expected results based on R test: 2nd pair→180°, 4th pair→45°
    expected = [PI, PI/4]  # 180° and 45°
    assert all(np.isclose(signed_angles, expected)), "with selection mask mismatch"


def test_edge_cases():
    """Test edge cases and precision"""

    # Test 1: Nearly identical vectors
    epsilon = 1e-15
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, epsilon], dtype=np.float64)

    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)

    assert np.isclose(result, 0), "expected zero angle for nearly identical vectors"

    # Test 2: Denormalized vectors (large magnitude)
    v1 = np.array([100.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 50.0], dtype=np.float64)

    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)

    expected = PI/2
    assert np.isclose(result, expected), "expected PI/2 angle"

    # Test 3: Very small vectors
    tiny = 1e-14
    v1 = np.array([tiny, 0.0], dtype=np.float64)
    v2 = np.array([0.0, tiny], dtype=np.float64)

    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)
    assert np.isclose(result, expected), "expected PI/2 angle even for short vectors"


def test_consistency_between_functions():
    """Test consistency between single and batch functions"""

    # Test same calculation with both functions
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0], dtype=np.float64)

    # Single function
    single_result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)

    # Batch function with single pair - create like R
    fields = np.array([[
        [1.0, 0.0],  # origin
        [0.0, 1.0]  # target
    ]], dtype=np.float64, order="F").transpose()
    vecs_selection_mask = np.array([1], dtype=np.int32)
    batch_results = clock_hand_angles_for_shift_vectors(fields, vecs_selection_mask, ANY_AXES)

    assert np.isclose(single_result, batch_results[0]), "single result should match batch result"


def test_mathematical_properties():
    """Test mathematical properties (anti-commutativity)"""

    # Test anti-commutativity: angle(v1,v2) = -angle(v2,v1)
    v1 = np.array([1.0, 2.0], dtype=np.float64)
    v2 = np.array([3.0, 1.0], dtype=np.float64)

    # Normalize vectors
    v1 = v1 / np.sqrt(np.sum(v1**2))
    v2 = v2 / np.sqrt(np.sum(v2**2))


    result_12 = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)
    result_21 = clock_hand_angle_between_vectors(v2, v1, ANY_AXES)

    assert np.isclose(result_12, -result_21), "anti commutativity failed"


def test_performance():
    """Performance test with large-scale data"""

    # Large-scale test
    n_dims = 50
    n_vecs = 1000

    # Generate random-like data
    np.random.seed(12345)
    fields = np.random.randn(n_dims * 2 * n_vecs).astype(np.float64).reshape((n_dims, 2, n_vecs), order="F")

    # Normalize vectors (for more meaningful angles)
    for i in range(n_vecs):
        origin_vec = fields[:, 0, i].view()
        target_vec = fields[:, 1, i].view()

        origin_norm = np.sqrt(np.sum(origin_vec**2))
        target_norm = np.sqrt(np.sum(target_vec**2))

        if origin_norm > 0:
            fields[:, 0, i] /= origin_norm
        if target_norm > 0:
            fields[:, 1, i] /= target_norm

    vecs_selection_mask = np.ones(n_vecs, dtype=np.int32)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)

    # Time the operation
    start_time = time.time()

    signed_angles = clock_hand_angles_for_shift_vectors(fields, vecs_selection_mask, selected_axes)

    end_time = time.time()
    elapsed = end_time - start_time

    assert elapsed < 1, "performance very bad, should be faster than a second"


def test_invalid_selected_axes():
    """Test invalid selected axes (same indices)"""
    v1 = np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 1, 1, 0, 0], dtype=np.int32)
    assert_error(lambda: clock_hand_angle_between_vectors(v1, v2, selected_axes), "Expected Exception")


def test_out_of_bounds_selected_axes():
    """Test out-of-bounds selected axes"""
    v1 = np.array([1.0, 0.0, 0.0, 2.3], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0, 2.3], dtype=np.float64)
    selected_axes = np.array([1, 2, 5, 6], dtype=np.int32)
    assert_error(lambda: clock_hand_angle_between_vectors(v1, v2, selected_axes), "Expected Exception")


def test_zero_vectors():
    """Test zero vectors (should not produce NaN)"""
    v1 = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)
    result = clock_hand_angle_between_vectors(v1, v2, selected_axes)

    # TODO: Will never be the case, as there is no division inside, but instead we get a valid angle. What to do?
    assert not np.isnan(result), "Expected non-NaN angle for zero vectors"


def test_denormalized_vectors():
    """Test denormalized vectors (large magnitude)"""
    v1 = np.array([100.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 50.0], dtype=np.float64)
    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)
    expected = PI/2
    assert np.isclose(result, expected), "Expected PI/2 angle"


def test_tiny_vectors_precision():
    """Test tiny vectors near machine precision"""
    tiny = 1e-14
    v1 = np.array([tiny, 0.0], dtype=np.float64)
    v2 = np.array([0.0, tiny], dtype=np.float64)
    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)
    expected = PI/2
    assert np.isclose(result, expected), "expected PI/2 angle"


def test_huge_vectors_precision():
    """Test huge vectors near overflow"""
    huge_val = 1e14
    v1 = np.array([huge_val, 0.0], dtype=np.float64)
    v2 = np.array([0.0, huge_val], dtype=np.float64)
    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)
    expected = PI/2
    assert np.isclose(result, expected), "expected PI/2 angle"


def test_nearly_identical_vectors():
    """Test nearly identical vectors (precision boundary)"""
    epsilon = 1e-15
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, epsilon], dtype=np.float64)
    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)
    assert np.isclose(result, 0), "expected zero angle"


def test_nearly_opposite_vectors():
    """Test nearly opposite vectors (precision boundary)"""
    epsilon = 1e-15
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([-1.0, epsilon], dtype=np.float64)
    result = clock_hand_angle_between_vectors(v1, v2, ANY_AXES)
    assert np.isclose(result, PI), "expected close PI angle"


def test_mixed_positive_negative():
    """Test vectors with mixed positive/negative components"""
    v1 = np.array([1.0, -2.0, 3.0, 0.0], dtype=np.float64)
    v2 = np.array([-2.0, 1.0, -3.0, 0.0], dtype=np.float64)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)
    v1 = v1 / np.sqrt(np.sum(v1**2))
    v2 = v2 / np.sqrt(np.sum(v2**2))
    c1 = v1[1]*v2[2] - v1[2]*v2[1]
    c2 = v1[2]*v2[0] - v1[0]*v2[2]
    c3 = v1[0]*v2[1] - v1[1]*v2[0]

    result = clock_hand_angle_between_vectors(v1, v2, selected_axes)
    sign = np.sign(c1 + c2 + c3)
    assert np.isclose(result, sign * np.acos(v1 @ v2)), "Angle mismatch"


if __name__ == "__main__":
    run_all_tests(globals().values())
