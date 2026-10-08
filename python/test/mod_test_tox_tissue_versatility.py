"""
The `tox_tissue_versatility` suite: each published procedure can be called from Python, returns the
documented type and shape, and raises the documented error. Whether the numbers are right is the
Fortran suite's job (test/mod_test_tox_tissue_versatility.F90), as the coding guide's "Where to Test
What" asks.

The one value needs no arithmetic beyond the trivial: over four selected axes a uniform vector
scores 0 at an angle of 0 and a single-axis vector scores 1, both exactly. It only shows that the
vectors mask selects and keeps the order, and that the axes mask reaches the library the right way
round: the fifth axis, which the mask leaves out, would break every uniform vector.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import compute_tissue_versatility
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import (
    ERR_EMPTY_INPUT,
    ERR_INVALID_INPUT,
    ERR_NAN_INF,
)

N_AXES, N_VECTORS = 5, 3

#: every axis but the last
AXES_SELECTION_MASK = np.array([True, True, True, True, False])


def _assert_result_array(result, shape, dtype, name, fortran_order):
    assert isinstance(result, np.ndarray), f"{name}: expected a numpy array, got {type(result).__name__}"
    assert result.dtype == dtype, f"{name}: expected {np.dtype(dtype)}, got {result.dtype}"
    assert result.shape == shape, f"{name}: expected shape {shape}, got {result.shape}"
    assert not result.flags.writeable, f"{name}: a result is a value, it should be read-only"
    if fortran_order:
        assert result.flags.f_contiguous, f"{name}: expected column-major (order='F')"
    else:
        assert result.flags.c_contiguous, f"{name}: expected a contiguous array"


def _expression_vectors():
    """(axes, vectors) in C order, so the binding has to convert it. Over the first four axes,
    vectors 1 and 3 are uniform and vector 2 sits on axis 4 alone; the fifth axis breaks the
    uniform ones, and a reversed axes mask would read it."""
    return np.array([[1.0, 0.0, 4.0],
                     [1.0, 0.0, 4.0],
                     [1.0, 0.0, 4.0],
                     [1.0, 3.0, 4.0],
                     [9.0, 0.0, 1.0]])


# -----------------------------------------------------------------------------------------------
# compute_tissue_versatility
# -----------------------------------------------------------------------------------------------

def test_compute_tissue_versatility():
    result = compute_tissue_versatility(_expression_vectors(), [False, True, True], AXES_SELECTION_MASK)
    assert isinstance(result, dict) and list(result) == ["tissue_versatilities", "tissue_angles_deg"], \
        f"expected a dict of tissue_versatilities and tissue_angles_deg, got {result!r}"
    tissue_versatilities = result["tissue_versatilities"]
    tissue_angles_deg = result["tissue_angles_deg"]
    _assert_result_array(tissue_versatilities, (2,), np.float64, "tissue_versatilities", fortran_order=False)
    _assert_result_array(tissue_angles_deg, (2,), np.float64, "tissue_angles_deg", fortran_order=False)
    # the one value: vector 2 (on one axis) then vector 3 (uniform). A reversed vectors mask would
    # give [0, 1], a reversed axes mask or none at all a nonzero second entry
    assert np.array_equal(tissue_versatilities, [1.0, 0.0]), f"tissue_versatilities: got {tissue_versatilities}"
    assert tissue_angles_deg[1] == 0.0, f"tissue_angles_deg: got {tissue_angles_deg}"


def test_compute_tissue_versatility_rejects_no_selected_vectors():
    # an empty selection is an error, not an empty result; the error names the output
    assert_error(lambda: compute_tissue_versatility(_expression_vectors(), np.zeros(N_VECTORS, dtype=bool),
                                                    AXES_SELECTION_MASK),
                 "no selected vectors", ERR_EMPTY_INPUT)


def test_compute_tissue_versatility_rejects_no_selected_axes():
    assert_error(lambda: compute_tissue_versatility(_expression_vectors(), np.ones(N_VECTORS, dtype=bool),
                                                    np.zeros(N_AXES, dtype=bool)),
                 "no selected axes", ERR_INVALID_INPUT)


def test_compute_tissue_versatility_rejects_no_axes():
    assert_error(lambda: compute_tissue_versatility(np.empty((0, N_VECTORS)), np.ones(N_VECTORS, dtype=bool),
                                                    np.empty(0, dtype=bool)),
                 "no axes", ERR_EMPTY_INPUT)


def test_compute_tissue_versatility_rejects_no_vectors():
    assert_error(lambda: compute_tissue_versatility(np.empty((N_AXES, 0)), np.empty(0, dtype=bool),
                                                    AXES_SELECTION_MASK),
                 "no vectors", ERR_EMPTY_INPUT)


def test_compute_tissue_versatility_rejects_nan():
    # the NaN sits in vector 1, which the mask leaves out: the whole matrix is checked
    expression_vectors = _expression_vectors()
    expression_vectors[2, 0] = np.nan
    assert_error(lambda: compute_tissue_versatility(expression_vectors, [False, True, True], AXES_SELECTION_MASK),
                 "a NaN expression", ERR_NAN_INF)


def test_compute_tissue_versatility_rejects_infinity():
    # the infinity sits on axis 5, which the mask leaves out: the whole matrix is checked
    expression_vectors = _expression_vectors()
    expression_vectors[4, 2] = np.inf
    assert_error(lambda: compute_tissue_versatility(expression_vectors, [False, True, True], AXES_SELECTION_MASK),
                 "an infinite expression", ERR_NAN_INF)


def test_compute_tissue_versatility_rejects_a_short_vectors_selection_mask():
    # the binding's own extent check, before the library is called
    assert_error(lambda: compute_tissue_versatility(_expression_vectors(), np.ones(N_VECTORS - 1, dtype=bool),
                                                    AXES_SELECTION_MASK),
                 "vectors_selection_mask shorter than n_vectors")


def test_compute_tissue_versatility_rejects_a_long_axes_selection_mask():
    # the binding's own extent check, before the library is called
    assert_error(lambda: compute_tissue_versatility(_expression_vectors(), np.ones(N_VECTORS, dtype=bool),
                                                    np.ones(N_AXES + 1, dtype=bool)),
                 "axes_selection_mask longer than n_axes")


def test_compute_tissue_versatility_rejects_a_vector():
    # the binding's own rank check, before the library is called
    assert_error(lambda: compute_tissue_versatility(np.ones(3), [True], [True, True, True]),
                 "expression_vectors must be 2-D")


if __name__ == '__main__':
    run_all_tests(globals().values())
