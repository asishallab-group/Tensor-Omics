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
from test_helpers import run_all_tests, assert_error


# 1. Uniform expression (should yield TV=0)
def test_uniform_expression():
    expr = np.full((3, 1), 2.0)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12


# 2. Single axis expression (should yield TV=1)
def test_single_axis_expression():
    expr = np.array([[0],[0],[5]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True, True])
    assert abs(res['tissue_versatilities'][0] - 1) < 1e-12
    assert res['tissue_angles_deg'][0] > 0


# 3. Null vector (should yield TV=1, angle=90)
def test_null_vector():
    expr = np.zeros((3,1), dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True, True])
    assert abs(res['tissue_versatilities'][0] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0] - 90) < 1e-12


# 4. Partial axis selection (subspace)
def test_partial_axis_selection():
    expr = np.array([[1],[2],[3]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, False, True])
    assert 0 <= res['tissue_versatilities'][0] <= 1
    assert 0 <= res['tissue_angles_deg'][0] <= 90


# 5. Mixed vectors (uniform, single axis, null)
def test_mixed_vectors():
    expr = np.array([[1,0,0],[1,0,0],[1,2,0]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True, True, True], [True, True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1] - 1) < 1e-12
    assert abs(res['tissue_versatilities'][2] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert res['tissue_angles_deg'][1] > 0
    assert abs(res['tissue_angles_deg'][2] - 90) < 1e-12


# 6. Angle output in degrees for a known case (should be 45)
def test_angle_degrees():
    expr = np.array([[1],[0]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True], [True, True])
    assert abs(res['tissue_angles_deg'][0] - 45) < 1e-12


# 7. Multiple vectors selection
def test_multiple_vectors_selection():
    expr = np.array([[1,0,0],[1,2,0]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True, False, True], [True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-5
    assert abs(res['tissue_angles_deg'][1] - 90) < 1e-12


# 8. High-dimensional vectors (4D, 5D)
def test_high_dimensional_vectors():
    expr4 = np.full((4,1), 1.0)
    expr5 = np.full((5,1), 2.0)
    res4 = tox_calculate_tissue_versatility(expr4, [True], [True, True, True, True])
    res5 = tox_calculate_tissue_versatility(expr5, [True], [True, True, True, True, True])
    assert abs(res4['tissue_versatilities'][0]) < 1e-12
    assert abs(res4['tissue_angles_deg'][0]) < 1e-12
    assert abs(res5['tissue_versatilities'][0]) < 1e-12
    assert abs(res5['tissue_angles_deg'][0]) < 1e-5


# 9. Randomized vectors and axes
def test_randomized_vectors_axes():
    np.random.seed(42)
    n_axes = 5
    n_vecs = 4
    expr = np.random.rand(n_axes, n_vecs)
    res = tox_calculate_tissue_versatility(expr, [True]*n_vecs, [True, False, True, False, True])
    assert np.all((res['tissue_versatilities'] >= 0) & (res['tissue_versatilities'] <= 1))
    assert np.all((res['tissue_angles_deg'] >= 0) & (res['tissue_angles_deg'] <= 90))


# 10. Numerical stability (very large/small values)
def test_numerical_stability():
    expr = np.array([[1e15,1e-4],[1e15,1e-4],[1e15,1e-4]], dtype=np.float64)
    res = tox_calculate_tissue_versatility(expr, [True, True], [True, True, True])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1]) < 1e-12
    assert abs(res['tissue_angles_deg'][1]) < 1e-12


# 11. Invalid input: no axes selected (should raise RuntimeError)
def test_invalid_input_no_axes():
    expr = np.array([[1],[2],[3]], dtype=np.float64)
    assert_error(lambda: tox_calculate_tissue_versatility(expr, [True], [False, False, False]), "Expected error for no selected axes")


# 12. Multiple selection, partial axes
def test_multiple_selection_partial_axes():
    expr = np.array([[1, 3, 5],[2, 4, 6]], dtype=np.float64, order="F")

    res = tox_calculate_tissue_versatility(expr, [True, False, True], [True, False])
    assert len(res['tissue_versatilities']) == 2
    assert np.all((res['tissue_versatilities'] >= 0) & (res['tissue_versatilities'] <= 1))


if __name__ == "__main__":
    run_all_tests(globals().values())
