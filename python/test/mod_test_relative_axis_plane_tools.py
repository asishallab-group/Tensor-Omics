import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensor_omics import omics_vector_RAP_projection, omics_field_RAP_projection, normalize_unit_length
from tensor_omics.error_handling import ERR_DIVISION_BY_ZERO
from test_helpers import run_all_tests, assert_error

# (name, value, n_axes): vectors whose coordinates are all equal. Many of them used to project
# to a residue of rounding noise instead of zeros, [0.7]*6 to about -1.1e-16 per coordinate.
_SMALLEST_SUBNORMAL = np.nextafter(0.0, 1.0)
_EQUAL_COORDINATE_CASES = [
    ("[0.7]*6", 0.7, 6),
    ("[0.1]*3", 0.1, 3),
    ("[0.1]*7", 0.1, 7),
    ("[3.3]*13", 3.3, 13),
    ("[1e300]*3", 1e300, 3),
    ("[-2.5]*4", -2.5, 4),
    ("[tiny]*5", np.finfo(np.float64).tiny, 5),
    ("[tiny/4]*3, subnormal", np.finfo(np.float64).tiny / 4, 3),
    ("[smallest subnormal]*7", _SMALLEST_SUBNORMAL, 7),
    ("[5]*1, a single axis", 5.0, 1),
]


def _assert_exact_zeros_that_normalization_refuses(projection, case_name):
    assert np.all(projection == 0.0), f"{case_name}: projection {projection!r} is not exact zeros"
    column = np.array(projection, dtype=np.float64)
    assert_error(lambda: normalize_unit_length(column),
                 f"{case_name}: normalize_unit_length must refuse the zero vector",
                 code=ERR_DIVISION_BY_ZERO)


def test_omics_vector_RAP_projection_call():
    n_axes = 5
    n_selected_axes = 2
    n_vecs = 10
    n_selected_vecs = 5

    vecs = np.random.rand(n_axes, n_vecs)
    vecs_selection_mask = np.zeros(n_vecs, dtype=np.int32)
    vecs_selection_mask[np.random.choice(n_vecs, n_selected_vecs, replace=False)] = 1
    axes_selection_mask = np.zeros(n_axes, dtype=np.int32)
    axes_selection_mask[np.random.choice(n_axes, n_selected_axes, replace=False)] = 1

    projections = omics_vector_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

    for i_vec in range(projections.shape[1]):
        col = projections[:, i_vec]
        assert np.isclose(col.sum(), 0), f"{i_vec}. column not correct: Sum is {col.sum()}, should be zero"


def test_omics_field_RAP_projection_call():
    n_axes = 5
    n_selected_axes = 2
    n_vecs = 10
    n_selected_vecs = 5

    # fields are (n_axes, 2, n_fields): [:, 0, :] origins, [:, 1, :] targets
    vecs = np.asfortranarray(np.random.rand(n_axes, 2, n_vecs))
    vecs_selection_mask = np.zeros(n_vecs, dtype=np.int32)
    vecs_selection_mask[np.random.choice(n_vecs, n_selected_vecs, replace=False)] = 1
    axes_selection_mask = np.zeros(n_axes, dtype=np.int32)
    axes_selection_mask[np.random.choice(n_axes, n_selected_axes, replace=False)] = 1

    projections = omics_field_RAP_projection(vecs, vecs_selection_mask, axes_selection_mask)

    for i_vec in range(projections.shape[1]):
        col = projections[:, i_vec]
        assert np.isclose(col.sum(), 0), f"{i_vec}. column not correct: Sum is {col.sum()}, should be zero"


def test_omics_field_RAP_projection_refuses_fields_without_two_points():
    """`fields` is (n_axes, 2, n_fields); a middle axis of 1 used to be read past (#221)."""
    fields = np.asfortranarray(np.random.rand(3, 1, 2))
    mask = np.ones(2, dtype=bool)

    assert_error(lambda: omics_field_RAP_projection(fields, mask, np.ones(3, dtype=bool)),
                 "a field without its target must be refused")



def test_omics_vector_RAP_projection_equal_coordinates_give_exact_zeros():
    for case_name, value, n_axes in _EQUAL_COORDINATE_CASES:
        vecs = np.full((n_axes, 1), value, order="F")
        projections = np.asarray(omics_vector_RAP_projection(
            vecs, np.ones(1, dtype=bool), np.ones(n_axes, dtype=bool)))
        _assert_exact_zeros_that_normalization_refuses(projections[:, 0], case_name)


def test_omics_field_RAP_projection_equal_shifts_give_exact_zeros():
    for case_name, value, n_axes in _EQUAL_COORDINATE_CASES:
        fields = np.zeros((n_axes, 2, 1), order="F")
        fields[:, 0, 0] = value  # the origin; the target stays at 0, so every shift is `value`
        projections = np.asarray(omics_field_RAP_projection(
            fields, np.ones(1, dtype=bool), np.ones(n_axes, dtype=bool)))
        _assert_exact_zeros_that_normalization_refuses(projections[:, 0], "shift " + case_name)


if __name__ == "__main__":
    run_all_tests(globals().values())
