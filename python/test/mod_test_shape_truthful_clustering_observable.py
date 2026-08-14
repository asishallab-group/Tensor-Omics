"""
Python test suite for normal_error / tangent_scales / observable
(tox_shape_truthful_clustering_observable), mirroring
test/mod_test_shape_truthful_clustering_observable.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import normal_error, tangent_scales, observable, ensemble_final_observable
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT, ERR_EMPTY_INPUT


# =====================
# normal_error
# =====================
def test_normal_error_basic():
    # D=3, d=1 -- normal_error is the sum of the two smallest (normal-space) eigenvalues.
    eigenvalues = np.array([9.0, 4.0, 1.0])
    value = normal_error(1, eigenvalues)
    assert abs(value - 5.0) < 1e-12, f"expected 5.0, got {value}"


def test_normal_error_zero_tangent_dims():
    # d=0: no tangent directions, everything is "normal" -- the full sum.
    eigenvalues = np.array([9.0, 4.0, 1.0])
    value = normal_error(0, eigenvalues)
    assert abs(value - 14.0) < 1e-12, f"expected 14.0, got {value}"


def test_normal_error_all_tangent_dims():
    # d=D: every direction is tangent, nothing left over -- sum over the empty range is zero.
    eigenvalues = np.array([9.0, 4.0, 1.0])
    value = normal_error(3, eigenvalues)
    assert abs(value - 0.0) < 1e-12, f"expected 0.0, got {value}"


def test_normal_error_d_out_of_range():
    eigenvalues = np.array([9.0, 4.0, 1.0])
    assert_error(lambda: normal_error(4, eigenvalues), "Expected error for d > n_dimensions", ERR_INVALID_INPUT)


def test_normal_error_negative_eigenvalue():
    eigenvalues = np.array([9.0, -1.0, 1.0])
    assert_error(lambda: normal_error(1, eigenvalues), "Expected error for a negative eigenvalue", ERR_INVALID_INPUT)


def test_normal_error_zero_dimensions():
    eigenvalues = np.array([])
    assert_error(lambda: normal_error(0, eigenvalues), "Expected error for n_dimensions=0", ERR_EMPTY_INPUT)


# =====================
# tangent_scales
# =====================
def test_tangent_scales_basic():
    # D=3, d=2 -- tangent_scales is the square root of the two largest eigenvalues.
    eigenvalues = np.array([9.0, 4.0, 1.0])
    value = tangent_scales(2, eigenvalues)
    expected = np.array([3.0, 2.0])
    assert np.allclose(value, expected, atol=1e-12), f"expected {expected}, got {value}"


def test_tangent_scales_all_dims():
    eigenvalues = np.array([9.0, 4.0, 1.0])
    value = tangent_scales(3, eigenvalues)
    expected = np.array([3.0, 2.0, 1.0])
    assert np.allclose(value, expected, atol=1e-12), f"expected {expected}, got {value}"


def test_tangent_scales_zero_dims():
    # d=0: no tangent directions -- must return a well-defined, empty array.
    eigenvalues = np.array([9.0, 4.0, 1.0])
    value = tangent_scales(0, eigenvalues)
    assert len(value) == 0, f"expected an empty array, got {value}"


def test_tangent_scales_d_out_of_range():
    eigenvalues = np.array([9.0, 4.0, 1.0])
    assert_error(lambda: tangent_scales(4, eigenvalues), "Expected error for d > n_dimensions", ERR_INVALID_INPUT)


def test_tangent_scales_negative_eigenvalue():
    eigenvalues = np.array([9.0, -4.0, 1.0])
    assert_error(lambda: tangent_scales(2, eigenvalues), "Expected error for a negative eigenvalue", ERR_INVALID_INPUT)


# =====================
# observable
# =====================
def test_observable_full_rank_rectangle():
    # A rectangle in the z=0 plane, embedded in 3D: full economy-mode rank (rank=min(D,k)=3=D,
    # no zero-padding), with distinct, hand-computable eigenvalues.
    vectors = np.array([
        [0.0, 2.0, 0.0, 2.0],
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0],
    ], order='F')
    mask = np.ones(4, dtype=np.bool_)

    result = observable(vectors, mask)
    assert np.allclose(result['mu'], [1.0, 0.5, 0.0], atol=1e-6), f"got {result['mu']}"
    assert result['d'] == 2, f"expected d=2, got {result['d']}"
    assert abs(result['eigenvalues'][0] - 4.0 / 3.0) < 1e-6
    assert abs(result['eigenvalues'][1] - 1.0 / 3.0) < 1e-6
    assert abs(result['eigenvalues'][2] - 0.0) < 1e-6
    assert abs(result['normal_error_value'] - 0.0) < 1e-6
    assert abs(result['tangent_scales_value'][0] - np.sqrt(4.0 / 3.0)) < 1e-6
    assert abs(result['tangent_scales_value'][1] - np.sqrt(1.0 / 3.0)) < 1e-6
    assert result['G'] > 1e10, f"expected a huge spectral gap, got {result['G']}"
    # U columns are the standard basis up to sign (diagonal covariance, distinct eigenvalues).
    U = result['U']
    assert np.allclose(np.abs(U[:, 0]), [1.0, 0.0, 0.0], atol=1e-6)
    assert np.allclose(np.abs(U[:, 1]), [0.0, 1.0, 0.0], atol=1e-6)
    assert np.allclose(np.abs(U[:, 2]), [0.0, 0.0, 1.0], atol=1e-6)


def test_observable_low_rank_padding():
    # Three collinear points (intrinsic rank 1) embedded in a 5D ambient space: economy-mode
    # rank = min(D,k) = 3 < D = 5, so U columns 4-5 and eigenvalues 4-5 must be observable's
    # own zero-padding, not LAPACK output.
    vectors = np.zeros((5, 3), order='F')
    vectors[0, :] = [0.0, 1.0, 2.0]
    mask = np.ones(3, dtype=np.bool_)

    result = observable(vectors, mask)
    assert result['d'] == 1, f"expected d=1, got {result['d']}"
    assert abs(result['eigenvalues'][0] - 1.0) < 1e-8
    assert np.allclose(result['eigenvalues'][3:5], [0.0, 0.0], atol=1e-8)
    assert np.allclose(result['U'][:, 3], np.zeros(5), atol=1e-8)
    assert np.allclose(result['U'][:, 4], np.zeros(5), atol=1e-8)
    assert abs(result['tangent_scales_value'][0] - 1.0) < 1e-8
    assert np.allclose(result['tangent_scales_value'][1:5], [0.0, 0.0, 0.0, 0.0], atol=1e-8)
    assert result['G'] > 1e10, f"expected a huge spectral gap, got {result['G']}"


def test_observable_too_few_members():
    vectors = np.array([
        [0.0, 2.0, 0.0, 2.0],
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 0.0, 0.0],
    ], order='F')
    mask = np.zeros(4, dtype=np.bool_)
    mask[0] = True

    assert_error(lambda: observable(vectors, mask),
                 "Expected error for an ensemble with fewer than 2 members", ERR_INVALID_INPUT)


def test_observable_dimension_too_small():
    vectors = np.array([[0.0, 1.0, 2.0]], order='F')
    mask = np.ones(3, dtype=np.bool_)

    assert_error(lambda: observable(vectors, mask),
                 "Expected error for n_dimensions=1", ERR_INVALID_INPUT)


# =====================
# ensemble_final_observable
# =====================
def test_ensemble_final_observable_trailing_rejected_column():
    # D=2, o=2, one ensemble: column 1 is the true accepted state (d=1, G=1.0, mu=[0.5,0.0]);
    # column 2 is a rejected candidate with deliberately different values (d=0, G=99.0,
    # mu=[9,9], accepted=False). The extraction must land on column 1, never column 2.
    U = np.zeros((2, 2, 2, 1), order='F')
    d = np.array([[1], [0]], dtype=np.int32, order='F')
    S = np.zeros((2, 2, 1), order='F')
    mu = np.zeros((2, 2, 1), order='F')
    G = np.array([[1.0], [99.0]], order='F')
    k = np.array([[2], [3]], dtype=np.int32, order='F')
    accepted = np.array([[True], [False]], dtype=np.bool_, order='F')

    S[:, 0, 0] = [0.5, 2.0]
    S[:, 1, 0] = [7.0, 7.0]
    mu[:, 0, 0] = [0.5, 0.0]
    mu[:, 1, 0] = [9.0, 9.0]
    U[:, 0, 0, 0] = [1.0, 0.0]
    U[:, 1, 0, 0] = [0.0, 1.0]

    result = ensemble_final_observable(U, d, S, mu, G, k, accepted)

    assert result['ensemble_has_final'][0]
    assert result['ensemble_final_index'][0] == 1
    assert result['ensemble_d_final'][0] == 1
    assert result['ensemble_k_final'][0] == 2
    assert abs(result['ensemble_G_final'][0] - 1.0) < 1e-12
    assert np.allclose(result['ensemble_mu_final'][:, 0], [0.5, 0.0], atol=1e-12)
    assert np.allclose(result['ensemble_S_final'][:, 0], [0.5, 2.0], atol=1e-12)


def test_ensemble_final_observable_no_rejection():
    # Both history columns accepted: extraction must land on the last (most recent) column.
    U = np.zeros((2, 2, 2, 1), order='F')
    d = np.array([[0], [1]], dtype=np.int32, order='F')
    S = np.zeros((2, 2, 1), order='F')
    mu = np.zeros((2, 2, 1), order='F')
    G = np.array([[2.0], [1.5]], order='F')
    k = np.array([[2], [3]], dtype=np.int32, order='F')
    accepted = np.array([[True], [True]], dtype=np.bool_, order='F')

    S[:, 0, 0] = [2.0, 0.0]
    S[:, 1, 0] = [0.5, 3.0]
    mu[:, 0, 0] = [0.5, 0.0]
    mu[:, 1, 0] = [1.0, 0.0]

    result = ensemble_final_observable(U, d, S, mu, G, k, accepted)

    assert result['ensemble_has_final'][0]
    assert result['ensemble_final_index'][0] == 2
    assert result['ensemble_d_final'][0] == 1
    assert result['ensemble_k_final'][0] == 3
    assert abs(result['ensemble_G_final'][0] - 1.5) < 1e-12
    assert np.allclose(result['ensemble_mu_final'][:, 0], [1.0, 0.0], atol=1e-12)


def test_ensemble_final_observable_has_final_false():
    # Every history column has k=0 (unpopulated) -- has_final must be false and every
    # _final output must come back exactly zero.
    U = np.zeros((2, 2, 2, 1), order='F')
    d = np.zeros((2, 1), dtype=np.int32, order='F')
    S = np.zeros((2, 2, 1), order='F')
    mu = np.zeros((2, 2, 1), order='F')
    G = np.zeros((2, 1), order='F')
    k = np.zeros((2, 1), dtype=np.int32, order='F')
    accepted = np.zeros((2, 1), dtype=np.bool_, order='F')

    result = ensemble_final_observable(U, d, S, mu, G, k, accepted)

    assert not result['ensemble_has_final'][0]
    assert result['ensemble_final_index'][0] == 0
    assert result['ensemble_d_final'][0] == 0
    assert result['ensemble_k_final'][0] == 0
    assert abs(result['ensemble_G_final'][0] - 0.0) < 1e-12
    assert np.allclose(result['ensemble_mu_final'][:, 0], [0.0, 0.0], atol=1e-12)


def test_ensemble_final_observable_multi_ensemble_independence():
    # Two ensembles with independent histories: ensemble 1 has a trailing rejected column
    # (resolves to its own column 1), ensemble 2 has no rejection (resolves to its column 2).
    U = np.zeros((2, 2, 2, 2), order='F')
    d = np.array([[1, 0], [0, 1]], dtype=np.int32, order='F')
    S = np.zeros((2, 2, 2), order='F')
    mu = np.zeros((2, 2, 2), order='F')
    G = np.array([[1.0, 3.0], [99.0, 2.5]], order='F')
    k = np.array([[2, 4], [3, 5]], dtype=np.int32, order='F')
    accepted = np.array([[True, True], [False, True]], dtype=np.bool_, order='F')

    mu[:, 0, 0] = [0.5, 0.0]
    mu[:, 1, 0] = [9.0, 9.0]
    mu[:, 0, 1] = [2.0, 2.0]
    mu[:, 1, 1] = [3.0, 3.0]

    result = ensemble_final_observable(U, d, S, mu, G, k, accepted)

    assert result['ensemble_has_final'][0]
    assert result['ensemble_final_index'][0] == 1
    assert abs(result['ensemble_G_final'][0] - 1.0) < 1e-12

    assert result['ensemble_has_final'][1]
    assert result['ensemble_final_index'][1] == 2
    assert abs(result['ensemble_G_final'][1] - 2.5) < 1e-12
    assert result['ensemble_k_final'][1] == 5
    assert np.allclose(result['ensemble_mu_final'][:, 1], [3.0, 3.0], atol=1e-12)


def test_ensemble_final_observable_small_o_evicts_accepted():
    # o=1: the single available column holds a rejected candidate. No earlier accepted column
    # exists to fall back to -- has_final must be false.
    U = np.zeros((2, 2, 1, 1), order='F')
    d = np.zeros((1, 1), dtype=np.int32, order='F')
    S = np.zeros((2, 1, 1), order='F')
    mu = np.zeros((2, 1, 1), order='F')
    G = np.array([[99.0]], order='F')
    k = np.array([[3]], dtype=np.int32, order='F')
    accepted = np.array([[False]], dtype=np.bool_, order='F')

    mu[:, 0, 0] = [9.0, 9.0]

    result = ensemble_final_observable(U, d, S, mu, G, k, accepted)

    assert not result['ensemble_has_final'][0]
    assert result['ensemble_final_index'][0] == 0
    assert abs(result['ensemble_G_final'][0] - 0.0) < 1e-12


if __name__ == "__main__":
    run_all_tests(globals().values())
