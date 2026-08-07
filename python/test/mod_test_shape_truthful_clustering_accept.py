"""
Python test suite for accept_ensemble (tox_shape_truthful_clustering_accept), mirroring
test/mod_test_shape_truthful_clustering_accept.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import accept_ensemble
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT


def test_accept_ensemble_identical():
    # Identical basis, identical d, identical G: all three criteria trivially satisfied,
    # even at zero tolerance for d and G.
    U_t = np.eye(3, order='F')
    U_tp1 = U_t.copy(order='F')

    is_accepted = accept_ensemble(U_t, 2, 5.0, U_tp1, 2, 5.0, 0.1, 0, 0.0)
    assert is_accepted, "identical observables must be accepted"


def test_accept_ensemble_angle_exceeds_max():
    # A 60-degree rotation of a 1D tangent basis, rejected by an alpha_max of 30 degrees.
    theta = np.pi / 3.0
    alpha_max = np.pi / 6.0

    U_t = np.eye(2, order='F')
    U_tp1 = np.array([[np.cos(theta), 0.0], [np.sin(theta), 1.0]], order='F')

    is_accepted = accept_ensemble(U_t, 1, 1.0, U_tp1, 1, 1.0, alpha_max, 0, 1e10)
    assert not is_accepted, "a 60-degree tangent rotation must be rejected at alpha_max=30deg"


def test_accept_ensemble_angle_within_max():
    theta = np.pi / 3.0
    alpha_max = 7.0 * np.pi / 18.0  # 70 degrees

    U_t = np.eye(2, order='F')
    U_tp1 = np.array([[np.cos(theta), 0.0], [np.sin(theta), 1.0]], order='F')

    is_accepted = accept_ensemble(U_t, 1, 1.0, U_tp1, 1, 1.0, alpha_max, 0, 1e10)
    assert is_accepted, "a 60-degree tangent rotation must be accepted at alpha_max=70deg"


def test_accept_ensemble_d_mismatch_within_dmax():
    # d changed by 1 (2 -> 1): the angle criterion is vacuously satisfied (no common dimension
    # to compare), and d_max=1 tolerates the change.
    U_t = np.eye(3, order='F')
    U_tp1 = U_t.copy(order='F')

    is_accepted = accept_ensemble(U_t, 2, 1.0, U_tp1, 1, 1.0, 0.1, 1, 1e10)
    assert is_accepted, "a change in d of 1 must be accepted when d_max=1"


def test_accept_ensemble_d_mismatch_exceeds_dmax():
    U_t = np.eye(3, order='F')
    U_tp1 = U_t.copy(order='F')

    is_accepted = accept_ensemble(U_t, 2, 1.0, U_tp1, 1, 1.0, 0.1, 0, 1e10)
    assert not is_accepted, "a change in d of 1 must be rejected when d_max=0"


def test_accept_ensemble_g_ratio_exceeds_max():
    # G changes by a factor of 10 (ln(10) ~ 2.303): rejected at G_max=1.0.
    U_t = np.eye(3, order='F')
    U_tp1 = U_t.copy(order='F')

    is_accepted = accept_ensemble(U_t, 2, 1.0, U_tp1, 2, 10.0, 0.1, 0, 1.0)
    assert not is_accepted, "a 10x change in G must be rejected at G_max=1.0"


def test_accept_ensemble_nonpositive_g():
    U_t = np.eye(3, order='F')
    U_tp1 = U_t.copy(order='F')

    assert_error(lambda: accept_ensemble(U_t, 2, 0.0, U_tp1, 2, 1.0, 0.1, 0, 1e10),
                 "Expected error for G_t <= 0", ERR_INVALID_INPUT)


def test_accept_ensemble_d_out_of_range():
    U_t = np.eye(3, order='F')
    U_tp1 = U_t.copy(order='F')

    assert_error(lambda: accept_ensemble(U_t, 4, 1.0, U_tp1, 2, 1.0, 0.1, 0, 1e10),
                 "Expected error for d_t > n_dimensions", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
