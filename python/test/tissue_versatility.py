#!/usr/bin/env python3
"""
Comprehensive Python test suite for tissue versatility (mirrors Fortran and R unit tests)
Uses the modular tensoromics_functions module
"""

import numpy as np
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))
from tensoromics_functions import tox_calculate_tissue_versatility

# 1. Uniform expression (should yield TV=0)
def test_uniform_expression():
    """
    Test tissue versatility for uniform expression (should yield TV=0, angle=0).

    This test checks if the tissue versatility and angle are both effectively zero
    when the input expression is uniform across all dimensions.
    """
    expr = np.full((3, 1), 2.0)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    print("test_uniform_expression passed")

# 2. Single axis expression (should yield TV=1)
def test_single_axis_expression():
    """
    Test tissue versatility for single axis expression (should yield TV=1).

    This test verifies that the tissue versatility is 1 and the angle is positive
    when the input expression varies along a single axis.
    """
    expr = np.array([[0],[0],[5]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True, True])
    assert abs(res['tissue_versatilities'][0] - 1) < 1e-12
    assert res['tissue_angles_deg'][0] > 0
    print("test_single_axis_expression passed")

# 3. Null vector (should yield TV=1, angle=90)
def test_null_vector():
    """
    Test tissue versatility for a null vector (should yield TV=1, angle=90).

    This test checks the case where the input expression is a null vector,
    expecting the tissue versatility to be 1 and the angle to be 90 degrees.
    """
    expr = np.zeros((3,1), dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True, True])
    assert abs(res['tissue_versatilities'][0] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0] - 90) < 1e-12
    print("test_null_vector passed")

# 4. Partial axis selection (subspace)
def test_partial_axis_selection():
    """
    Test tissue versatility with partial axis selection (subspace).

    This test evaluates the tissue versatility and angle when only a subset of
    the axes are selected, effectively testing the function's behavior in a subspace.
    """
    expr = np.array([[1],[2],[3]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, False, True])
    assert 0 <= res['tissue_versatilities'][0] <= 1
    assert 0 <= res['tissue_angles_deg'][0] <= 90
    print("test_partial_axis_selection passed")

# 5. Mixed vectors (uniform, single axis, null)
def test_mixed_vectors():
    """
    Test tissue versatility for mixed vectors (uniform, single axis, null).

    This test combines different types of vectors (uniform, single-axis, null)
    to verify that the tissue versatility and angles are computed correctly for each.
    """
    expr = np.array([[1,0,0],[1,0,0],[1,2,0]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True, True, True], [True, True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1] - 1) < 1e-12
    assert abs(res['tissue_versatilities'][2] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert res['tissue_angles_deg'][1] > 0
    assert abs(res['tissue_angles_deg'][2] - 90) < 1e-12
    print("test_mixed_vectors passed")

# 6. Angle output in degrees for a known case (should be 45)
def test_angle_degrees():
    """
    Test angle output in degrees for a known case (should be 45).

    This test ensures that the angle is correctly calculated as 45 degrees
    for a standard basis vector in a 2D space.
    """
    expr = np.array([[1],[0]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True])
    assert abs(res['tissue_angles_deg'][0] - 45) < 1e-12
    print("test_angle_degrees passed")

# 7. Multiple vectors selection
def test_multiple_vectors_selection():
    """
    Test tissue versatility for multiple vector selection.

    This test checks the tissue versatility and angles when multiple vectors
    are selected, ensuring correct computation for each selected vector.
    """
    expr = np.array([[1,0,0],[1,2,0]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True, False, True], [True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][1] - 90) < 1e-12
    print("test_multiple_vectors_selection passed")

# 8. High-dimensional vectors (4D, 5D)
def test_high_dimensional_vectors():
    """
    Test tissue versatility for high-dimensional vectors (4D, 5D).

    This test evaluates the function's performance and correctness when
    handling vectors in higher dimensions, specifically 4D and 5D spaces.
    """
    expr4 = np.full((4,1), 1.0)
    expr5 = np.full((5,1), 2.0)
    res4 = tox_calculate_tissue_versatility(expr4, [True], [True, True, True, True])
    res5 = tox_calculate_tissue_versatility(expr5, [True], [True, True, True, True, True])
    assert abs(res4['tissue_versatilities'][0]) < 1e-12
    assert abs(res4['tissue_angles_deg'][0]) < 1e-12
    assert abs(res5['tissue_versatilities'][0]) < 1e-12
    assert abs(res5['tissue_angles_deg'][0]) < 1e-12
    print("test_high_dimensional_vectors passed")

# 9. Randomized vectors and axes
def test_randomized_vectors_axes():
    """
    Test tissue versatility for randomized vectors and axes.

    This test uses randomized input for vectors and axis selection to ensure
    the function's robustness and correctness across a wide range of inputs.
    """
    np.random.seed(42)
    n_axes = 5
    n_vecs = 4
    expr = np.random.rand(n_axes, n_vecs)
    res = tox_calculate_tissue_versatility(expr, [True]*n_vecs, [True, False, True, False, True])
    assert np.all((res['tissue_versatilities'] >= 0) & (res['tissue_versatilities'] <= 1))
    assert np.all((res['tissue_angles_deg'] >= 0) & (res['tissue_angles_deg'] <= 90))
    print("test_randomized_vectors_axes passed")

# 10. Numerical stability (very large/small values)
def test_numerical_stability():
    """
    Test numerical stability for very large and small values.

    This test checks the function's ability to handle extreme values in the input
    without losing precision or encountering errors.
    """
    expr = np.array([[1e15,1e-4],[1e15,1e-4],[1e15,1e-4]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True, True], [True, True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1]) < 1e-12
    assert abs(res['tissue_angles_deg'][1]) < 1e-12
    print("test_numerical_stability passed")

# 11. Invalid input: no axes selected (should raise RuntimeError)
def test_invalid_input_no_axes():
    """
    Test error handling for no axes selected (should raise RuntimeError).

    This test ensures that the function raises an appropriate error when
    no axes are selected in the input, preventing further computation.
    """
    expr = np.array([[1],[2],[3]], dtype=np.float64)
    error_raised = False
    try:
        tox_calculate_tissue_versatility(expr, [True], [False, False, False])
    except RuntimeError as e:
        error_raised = True
        # Check that the error message contains the expected text
        assert "Empty input arrays provided." in str(e)
    assert error_raised, "Expected RuntimeError was not raised"
    print("test_invalid_input_no_axes passed")

# 12. Multiple selection, partial axes
def test_multiple_selection_partial_axes():
    """
    Test tissue versatility for multiple selection with partial axes.

    This test evaluates the function's behavior and output when multiple selections
    are made with only partial axes being active, ensuring correct computation in
    such scenarios.
    """
    expr = np.array([[1,3,5],[2,4,6]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True, False, True], [True, False])
    assert len(res['tissue_versatilities']) == 2
    assert np.all((res['tissue_versatilities'] >= 0) & (res['tissue_versatilities'] <= 1))
    print("test_multiple_selection_partial_axes passed")


def main():
    """
    Run all tissue versatility test cases.

    This function executes all the individual test cases for the tissue versatility
    calculations, reporting the success or failure of each test.
    """
    print("=================================================")
    print("    TISSUE VERSATILITY FULL PYTHON INTERFACE TESTS")
    print("=================================================\n")
    
    test_uniform_expression()
    test_single_axis_expression()
    test_null_vector()
    test_partial_axis_selection()
    test_mixed_vectors()
    test_angle_degrees()
    test_multiple_vectors_selection()
    test_high_dimensional_vectors()
    test_randomized_vectors_axes()
    test_numerical_stability()
    test_invalid_input_no_axes()
    test_multiple_selection_partial_axes()
    
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all tissue versatility Python interface tests passed! ✓")

if __name__ == "__main__":
    main()
