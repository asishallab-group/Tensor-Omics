"""
The `tox_relative_axis_plane_tools` suite: each published procedure can be called from Python,
returns the documented type and shape, and raises the documented error. Whether the numbers are
right is the Fortran suite's job (test/mod_test_tox_relative_axis_plane_tools.F90), as the coding
guide's "Where to Test What" asks. The one hand-derived value per procedure below proves only
that the binding passes its data through unchanged: argument order, layout, selection.

The binding derives every `n_selected_*` from its selection mask, so the Fortran check that a
count agrees with its mask cannot be reached from here; an empty selection can.
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import (
    clock_hand_angle_between_vectors,
    clock_hand_angles_for_shift_vectors,
    compute_relative_axis_contributions,
    omics_field_RAP_projection,
    omics_vector_RAP_projection,
    relative_axes_changes_from_shift_vector,
    relative_axes_expression_from_expression_vector,
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import (
    ERR_DIVISION_BY_ZERO,
    ERR_EMPTY_INPUT,
    ERR_INVALID_INPUT,
    ERR_NAN_INF,
)

N_AXES, N_VECS = 4, 5
AXES_SELECTION_MASK = np.array([True, False, True, True])
VECS_SELECTION_MASK = np.array([True, True, False, True, False])
N_SELECTED_AXES, N_SELECTED_VECS = 3, 3

#: the three contribution procedures share one signature and one contract
CONTRIBUTION_PROCEDURES = (compute_relative_axis_contributions, relative_axes_changes_from_shift_vector,
                           relative_axes_expression_from_expression_vector)

#: what the Fortran rejects in every real input
NON_FINITE_VALUES = (np.nan, np.inf)


def _vecs():
    """(axes, vectors), finite and all different, so nothing degenerates."""
    return np.fromfunction(lambda i_axis, i_vec: (i_axis + 1) * (i_vec + 2.0) + i_axis ** 2, (N_AXES, N_VECS))


def _fields():
    """(axes, origin/shift, fields): origins from `_vecs`, and shifts off them."""
    fields = np.empty((N_AXES, 2, N_VECS), order='F')
    fields[:, 0, :] = _vecs()
    fields[:, 1, :] = _vecs() ** 1.5
    return fields


def _with(array, value):
    """`array` with its first element replaced by `value`."""
    array = np.array(array, dtype=np.float64, order='F')
    array.flat[0] = value
    return array


def _assert_result(result, shape, name, order='F'):
    assert isinstance(result, np.ndarray), f"{name}: expected a numpy array, got {type(result).__name__}"
    assert result.dtype == np.float64, f"{name}: expected float64, got {result.dtype}"
    assert result.shape == shape, f"{name}: expected shape {shape}, got {result.shape}"
    assert not result.flags.writeable, f"{name}: a result is a value and must be read-only"
    contiguous = result.flags.f_contiguous if order == 'F' else result.flags.c_contiguous
    assert contiguous, f"{name}: expected memory order '{order}'"
    assert np.all(np.isfinite(result)), f"{name}: non-finite values in the result"


# ---- omics_vector_RAP_projection -------------------------------------------------------------

def test_omics_vector_RAP_projection():
    result = omics_vector_RAP_projection(_vecs(), VECS_SELECTION_MASK, AXES_SELECTION_MASK)
    _assert_result(result, (N_SELECTED_AXES, N_SELECTED_VECS), "omics_vector_RAP_projection")


def test_omics_vector_RAP_projection_keeps_axes_as_rows():
    # pass-through only: the projection subtracts the selected column's mean, 30, from
    # (10, 20, 60); a transposed matrix or a misapplied mask would give other numbers
    vecs = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 60.0]])
    result = omics_vector_RAP_projection(vecs, [False, True], [True, True, True])
    assert np.allclose(result, [[-20.0], [-10.0], [30.0]]), f"got {result.tolist()}"


def test_omics_vector_RAP_projection_rejects_a_mask_of_the_wrong_length():
    # the binding's own shape check, before the library is called
    assert_error(lambda: omics_vector_RAP_projection(_vecs(), VECS_SELECTION_MASK, AXES_SELECTION_MASK[:-1]),
                 "vecs has four axes, the mask three")


def test_omics_vector_RAP_projection_rejects_non_finite_vecs():
    for value in NON_FINITE_VALUES:
        assert_error(lambda: omics_vector_RAP_projection(_with(_vecs(), value), VECS_SELECTION_MASK,
                                                         AXES_SELECTION_MASK),
                     f"vecs holds {value}", ERR_NAN_INF)


def test_omics_vector_RAP_projection_rejects_an_empty_selection():
    nothing_selected = np.zeros(N_AXES, dtype=np.bool_)
    assert_error(lambda: omics_vector_RAP_projection(_vecs(), VECS_SELECTION_MASK, nothing_selected),
                 "no axis selected", ERR_EMPTY_INPUT)


# ---- omics_field_RAP_projection --------------------------------------------------------------

def test_omics_field_RAP_projection():
    result = omics_field_RAP_projection(_fields(), VECS_SELECTION_MASK, AXES_SELECTION_MASK)
    _assert_result(result, (N_SELECTED_AXES, N_SELECTED_VECS), "omics_field_RAP_projection")


def test_omics_field_RAP_projection_projects_the_selected_shift():
    # pass-through only: each field is [origin, shift], and the shift is what is projected. The
    # selected second field's shift (1, 2, 6) is mean-free (-2, -1, 3), whatever its origin (100
    # here); the first field is never selected
    fields = np.zeros((3, 2, 2), order='F')
    fields[:, :, 0] = 100.0
    fields[:, 0, 1] = 100.0
    fields[:, 1, 1] = [1.0, 2.0, 6.0]
    result = omics_field_RAP_projection(fields, [False, True], [True, True, True])
    expected = np.array([[-2.0], [-1.0], [3.0]])
    assert np.allclose(result, expected), f"expected {expected.tolist()}, got {result.tolist()}"


def test_omics_field_RAP_projection_rejects_a_mask_of_the_wrong_length():
    # the binding's own shape check, before the library is called
    assert_error(lambda: omics_field_RAP_projection(_fields(), VECS_SELECTION_MASK[:-1], AXES_SELECTION_MASK),
                 "fields has five fields, the mask four")


def test_omics_field_RAP_projection_rejects_non_finite_fields():
    for value in NON_FINITE_VALUES:
        assert_error(lambda: omics_field_RAP_projection(_with(_fields(), value), VECS_SELECTION_MASK,
                                                        AXES_SELECTION_MASK),
                     f"fields holds {value}", ERR_NAN_INF)


def test_omics_field_RAP_projection_rejects_an_empty_selection():
    nothing_selected = np.zeros(N_VECS, dtype=np.bool_)
    assert_error(lambda: omics_field_RAP_projection(_fields(), nothing_selected, AXES_SELECTION_MASK),
                 "no field selected", ERR_EMPTY_INPUT)


# ---- clock_hand_angle_between_vectors --------------------------------------------------------

def test_clock_hand_angle_between_vectors():
    result = clock_hand_angle_between_vectors([1.0, 0.0, 0.0], [0.6, 0.8, 0.0], [0.0, 1.0, 0.0])
    assert isinstance(result, float), f"expected a float, got {type(result).__name__}"
    assert math.isfinite(result), f"expected a finite angle, got {result}"


def test_clock_hand_angle_between_vectors_keeps_the_argument_order():
    # pass-through only: the turn from e1 to e2 moves along e2, which the reference (-1, 1)
    # counts positive, so +pi/2; swapped vectors turn along e1, which it counts negative
    result = clock_hand_angle_between_vectors([1.0, 0.0], [0.0, 1.0], [-1.0, 1.0])
    assert math.isclose(result, math.pi / 2), f"expected pi/2, got {result}"


def test_clock_hand_angle_between_vectors_rejects_vectors_of_different_length():
    # the binding's own shape check, before the library is called
    assert_error(lambda: clock_hand_angle_between_vectors([1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0]),
                 "v1 has two dimensions, v2 three")


def test_clock_hand_angle_between_vectors_rejects_a_non_finite_reference():
    for value in NON_FINITE_VALUES:
        assert_error(lambda: clock_hand_angle_between_vectors([1.0, 0.0], [0.0, 1.0], [value, 1.0]),
                     f"the reference holds {value}", ERR_NAN_INF)


def test_clock_hand_angle_between_vectors_rejects_empty_vectors():
    empty = np.empty(0)
    assert_error(lambda: clock_hand_angle_between_vectors(empty, empty, empty), "no dimensions", ERR_EMPTY_INPUT)


def test_clock_hand_angle_between_vectors_rejects_a_reference_that_orients_nothing():
    # the documented ERR_INVALID_INPUT: the turn is in the (1, 2) plane, and axis 3 cannot sign it
    assert_error(lambda: clock_hand_angle_between_vectors([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
                 "a reference orthogonal to the rotation", ERR_INVALID_INPUT)


# ---- clock_hand_angles_for_shift_vectors -----------------------------------------------------

def test_clock_hand_angles_for_shift_vectors():
    result = clock_hand_angles_for_shift_vectors(_fields(), VECS_SELECTION_MASK, np.arange(1.0, N_AXES + 1))
    _assert_result(result, (N_SELECTED_VECS,), "clock_hand_angles_for_shift_vectors", order='C')


def test_clock_hand_angles_for_shift_vectors_turns_from_origin_to_origin_plus_shift():
    # pass-through only, the rule of the single-pair test above: each field is [origin, shift] and
    # turns from o to o + s. The selected second field turns from e1 to e1 + (-1, 1) = e2, +pi/2;
    # reading the shift as the target would give 3pi/4, and the unselected first field turns to
    # e1 + (-1, -1) = -e2, -pi/2
    fields = np.zeros((2, 2, 2), order='F')
    fields[:, 0, :] = [[1.0, 1.0], [0.0, 0.0]]
    fields[:, 1, 0] = [-1.0, -1.0]
    fields[:, 1, 1] = [-1.0, 1.0]
    result = clock_hand_angles_for_shift_vectors(fields, [False, True], [-1.0, 1.0])
    assert np.allclose(result, [math.pi / 2]), f"expected [pi/2], got {result.tolist()}"


def test_clock_hand_angles_for_shift_vectors_rejects_a_mask_of_the_wrong_length():
    # the binding's own shape check, before the library is called
    assert_error(lambda: clock_hand_angles_for_shift_vectors(_fields(), VECS_SELECTION_MASK[:-1],
                                                             np.ones(N_AXES)),
                 "fields has five fields, the mask four")


def test_clock_hand_angles_for_shift_vectors_rejects_non_finite_fields():
    for value in NON_FINITE_VALUES:
        assert_error(lambda: clock_hand_angles_for_shift_vectors(_with(_fields(), value), VECS_SELECTION_MASK,
                                                                 np.ones(N_AXES)),
                     f"fields holds {value}", ERR_NAN_INF)


def test_clock_hand_angles_for_shift_vectors_rejects_an_empty_selection():
    nothing_selected = np.zeros(N_VECS, dtype=np.bool_)
    assert_error(lambda: clock_hand_angles_for_shift_vectors(_fields(), nothing_selected, np.ones(N_AXES)),
                 "no field selected", ERR_EMPTY_INPUT)


def test_clock_hand_angles_for_shift_vectors_rejects_a_reference_that_orients_nothing():
    # the documented ERR_INVALID_INPUT, as for a single pair
    fields = np.zeros((3, 2, 1), order='F')
    fields[0, 0, 0] = 1.0
    fields[1, 1, 0] = 1.0
    assert_error(lambda: clock_hand_angles_for_shift_vectors(fields, [True], [0.0, 0.0, 1.0]),
                 "a reference orthogonal to the rotation", ERR_INVALID_INPUT)


# ---- the three contribution procedures -------------------------------------------------------

def test_contributions():
    for procedure in CONTRIBUTION_PROCEDURES:
        _assert_result(procedure(np.arange(1.0, N_AXES + 1)), (N_AXES,), procedure.__name__, order='C')


def test_contributions_keep_the_axis_order():
    # pass-through only: |1| and |-3| are a quarter and three quarters of 4; reversed, the other way
    for procedure in CONTRIBUTION_PROCEDURES:
        result = procedure([1.0, -3.0])
        assert np.allclose(result, [0.25, 0.75]), f"{procedure.__name__}: got {result.tolist()}"


def test_contributions_reject_a_zero_vector():
    for procedure in CONTRIBUTION_PROCEDURES:
        assert_error(lambda: procedure(np.zeros(3)), f"{procedure.__name__}: nothing to apportion",
                     ERR_DIVISION_BY_ZERO)


def test_contributions_reject_non_finite_values():
    for procedure in CONTRIBUTION_PROCEDURES:
        for value in NON_FINITE_VALUES:
            assert_error(lambda: procedure([1.0, value]), f"{procedure.__name__}: vec holds {value}", ERR_NAN_INF)


def test_contributions_reject_an_empty_vector():
    for procedure in CONTRIBUTION_PROCEDURES:
        assert_error(lambda: procedure(np.empty(0)), f"{procedure.__name__}: no axes", ERR_EMPTY_INPUT)


if __name__ == '__main__':
    run_all_tests(globals().values())
