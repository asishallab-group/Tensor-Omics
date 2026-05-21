"""
Unit tests for relative axis contributions Python wrappers
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensoromics_functions import (
    relative_axes_changes_from_shift_vector,
    relative_axes_expression_from_expression_vector
)
from test_helpers import run_all_tests, assert_error

TOL = 1e-12


def assert_vec_sum(vec, expected_sum, tol=TOL, msg=""):
    actual_sum = np.sum(vec)
    assert np.isclose(actual_sum, expected_sum, atol=tol), f"{msg} Expected sum {expected_sum}, got {actual_sum}"


def assert_all_finite(vec, msg=""):
    assert np.all(np.isfinite(vec)), f"{msg} Non-finite values detected."


def assert_in_range(vec, min_val=0.0, max_val=1.0, msg=""):
    assert not np.any((vec < min_val) | (vec > max_val)), f"{msg} Values out of range [{min_val}, {max_val}]"


# Tests for relative_axes_changes_from_shift_vector
def test_shift_positive_vector():
    vec = np.array([1, 2, 3], dtype=np.float64)
    contrib = relative_axes_changes_from_shift_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="shift positive vector")
    assert_all_finite(contrib, "shift positive vector")
    assert_in_range(contrib, msg="shift positive vector")


def test_shift_negative_vector():
    vec = np.array([-1, -2, -3], dtype=np.float64)
    contrib = relative_axes_changes_from_shift_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="shift negative vector")
    assert_all_finite(contrib, "shift negative vector")
    assert_in_range(contrib, msg="shift negative vector")


def test_shift_mixed_vector():
    vec = np.array([2, -2, 4], dtype=np.float64)
    contrib = relative_axes_changes_from_shift_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="shift mixed vector")
    assert_all_finite(contrib, "shift mixed vector")
    assert_in_range(contrib, msg="shift mixed vector")


def test_shift_zero_vector():
    vec = np.array([0, 0, 0], dtype=np.float64)
    assert_error(lambda: relative_axes_changes_from_shift_vector(vec), "Expected exception for zero vector")


def test_shift_one_nonzero_axis():
    vec = np.array([0, 5, 0], dtype=np.float64)
    contrib = relative_axes_changes_from_shift_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="shift one nonzero axis")
    assert_all_finite(contrib, "shift one nonzero axis")
    assert_in_range(contrib, msg="shift one nonzero axis")


def test_shift_all_equal():
    vec = np.array([2, 2, 2], dtype=np.float64)
    contrib = relative_axes_changes_from_shift_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="shift all equal")
    assert_all_finite(contrib, "shift all equal")
    assert_in_range(contrib, msg="shift all equal")


def test_shift_large_vector():
    vec = np.ones(100, dtype=np.float64)
    contrib = relative_axes_changes_from_shift_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="shift large vector")
    assert_all_finite(contrib, "shift large vector")
    assert_in_range(contrib, msg="shift large vector")


# Tests for relative_axes_expression_from_expression_vector
def test_expr_positive_vector():
    vec = np.array([1, 2, 3], dtype=np.float64)
    contrib = relative_axes_expression_from_expression_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="expr positive vector")
    assert_all_finite(contrib, "expr positive vector")
    assert_in_range(contrib, msg="expr positive vector")


def test_expr_negative_vector():
    vec = np.array([-1, -2, -3], dtype=np.float64)
    contrib = relative_axes_expression_from_expression_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="expr negative vector")
    assert_all_finite(contrib, "expr negative vector")
    assert_in_range(contrib, msg="expr negative vector")


def test_expr_mixed_vector():
    vec = np.array([2, -2, 4], dtype=np.float64)
    contrib = relative_axes_expression_from_expression_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="expr mixed vector")
    assert_all_finite(contrib, "expr mixed vector")
    assert_in_range(contrib, msg="expr mixed vector")


def test_expr_zero_vector():
    vec = np.array([0, 0, 0], dtype=np.float64)
    assert_error(lambda: relative_axes_expression_from_expression_vector(vec), "Expected exception for zero vector")


def test_expr_one_nonzero_axis():
    vec = np.array([0, 5, 0], dtype=np.float64)
    contrib = relative_axes_expression_from_expression_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="expr one nonzero axis")
    assert_all_finite(contrib, "expr one nonzero axis")
    assert_in_range(contrib, msg="expr one nonzero axis")


def test_expr_all_equal():
    vec = np.array([2, 2, 2], dtype=np.float64)
    contrib = relative_axes_expression_from_expression_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="expr all equal")
    assert_all_finite(contrib, "expr all equal")
    assert_in_range(contrib, msg="expr all equal")


def test_expr_large_vector():
    vec = np.ones(100, dtype=np.float64)
    contrib = relative_axes_expression_from_expression_vector(vec)
    assert_vec_sum(contrib, 1.0, msg="expr large vector")
    assert_all_finite(contrib, "expr large vector")
    assert_in_range(contrib, msg="expr large vector")


if __name__ == "__main__":
    run_all_tests(globals().values())
