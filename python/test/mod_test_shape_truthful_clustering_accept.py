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


def _identity_history(d_dim, o, d_value):
    U = np.eye(d_dim, order='F')
    U_history = np.zeros((d_dim, d_dim, o), dtype=np.float64, order='F')
    U_history[:, :, 0] = U
    d_history = np.zeros(o, dtype=np.int32)
    d_history[0] = d_value
    return U, U_history, d_history


def test_accept_ensemble_identical():
    # Identical bases at every retained iteration (U_first, U_history, U_tp1), identical d, G,
    # normal_error: all four criteria trivially satisfied, even at zero tolerance for d, G, RMSE.
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 2, U_history, d_history, 1, 5.0, 1.0, U, 2, 5.0, 1.0,
                                  0.1, 0, 0.0, 0.0)
    assert is_accepted, "identical observables must be accepted"


def test_accept_ensemble_chordal_exceeds_max():
    # A 60-degree rotation of a 1D tangent basis (U_first = U_history(1), isolating the
    # chordal-distance formula/threshold from the cumulative-drift machinery): chordal distance
    # = sin(60deg) ~ 0.866, rejected at a fraction-of-range of 0.5.
    theta = np.pi / 3.0
    U_t = np.eye(2, order='F')
    U_tp1 = np.array([[np.cos(theta), 0.0], [np.sin(theta), 1.0]], order='F')
    U_history = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    U_history[:, :, 0] = U_t
    d_history = np.array([1], dtype=np.int32)

    is_accepted = accept_ensemble(U_t, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.5, 0, 1e10, 1e10)
    assert not is_accepted, \
        "a 60-degree tangent rotation must be rejected at chordal_dist_max_as_prcnt_of_range=0.5"


def test_accept_ensemble_chordal_within_max():
    # The same 60-degree rotation, accepted once the fraction-of-range is raised to 0.9
    # (sin(60deg) ~ 0.866 <= 0.9).
    theta = np.pi / 3.0
    U_t = np.eye(2, order='F')
    U_tp1 = np.array([[np.cos(theta), 0.0], [np.sin(theta), 1.0]], order='F')
    U_history = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    U_history[:, :, 0] = U_t
    d_history = np.array([1], dtype=np.int32)

    is_accepted = accept_ensemble(U_t, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.9, 0, 1e10, 1e10)
    assert is_accepted, \
        "a 60-degree tangent rotation must be accepted at chordal_dist_max_as_prcnt_of_range=0.9"


def test_accept_ensemble_rejects_cumulative_drift_from_first():
    # The P5 regression this whole redesign fixes: a candidate only 5 degrees from the most
    # recently accepted state (U_history(1), at 75deg) -- which a step-to-step-only check would
    # accept -- but 80 degrees from the ensemble's own bootstrap state (U_first, at 0deg).
    U_first = np.eye(2, order='F')
    U_history = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    theta75 = 75.0 * np.pi / 180.0
    U_history[:, :, 0] = [[np.cos(theta75), 0.0], [np.sin(theta75), 1.0]]
    d_history = np.array([1], dtype=np.int32)

    theta80 = 80.0 * np.pi / 180.0
    U_tp1 = np.array([[np.cos(theta80), 0.0], [np.sin(theta80), 1.0]], order='F')

    is_accepted = accept_ensemble(U_first, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.5, 0, 1e10, 1e10)
    assert not is_accepted, \
        "a candidate close to the previous state but far from U_first must be rejected"


def test_accept_ensemble_accepts_small_drift_from_both():
    # Companion: small drift from both U_first (6deg) and U_history(1) (3deg) is accepted.
    U_first = np.eye(2, order='F')
    U_history = np.zeros((2, 2, 1), dtype=np.float64, order='F')
    theta3 = 3.0 * np.pi / 180.0
    U_history[:, :, 0] = [[np.cos(theta3), 0.0], [np.sin(theta3), 1.0]]
    d_history = np.array([1], dtype=np.int32)

    theta6 = 6.0 * np.pi / 180.0
    U_tp1 = np.array([[np.cos(theta6), 0.0], [np.sin(theta6), 1.0]], order='F')

    is_accepted = accept_ensemble(U_first, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.5, 0, 1e10, 1e10)
    assert is_accepted, "small drift from both U_first and U_history(1) must be accepted"


def test_accept_ensemble_d_to_first_exceeds_dmax():
    # d_to_last = |d_tp1 - d_history(1)| = 0 (fine on its own), but d_to_first =
    # |d_tp1 - d_first| = 2 > d_max=1: a d_to_last-only check would wrongly accept this.
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 0, U_history, d_history, 1, 1.0, 1.0, U, 2, 1.0, 1.0,
                                  0.9, 1, 1e10, 1e10)
    assert not is_accepted, \
        "max(d_to_first, d_to_last) must reject even when d_to_last alone is fine"


def test_accept_ensemble_d_two_fold_within_dmax():
    # Same setup, accepted once d_max is raised to 2 (>= max(d_to_first, d_to_last) = 2).
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 0, U_history, d_history, 1, 1.0, 1.0, U, 2, 1.0, 1.0,
                                  0.9, 2, 1e10, 1e10)
    assert is_accepted, "d_max=2 must tolerate max(d_to_first, d_to_last)=2"


def test_accept_ensemble_g_ratio_exceeds_max():
    # G changes by a factor of 10 (ln(10) ~ 2.303): rejected at G_max=1.0.
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 2, U_history, d_history, 1, 1.0, 1.0, U, 2, 10.0, 1.0,
                                  0.9, 0, 1.0, 1e10)
    assert not is_accepted, "a 10x change in G must be rejected at G_max=1.0"


def test_accept_ensemble_g_ratio_within_max():
    # Same 10x change in G, accepted once G_max is raised past ln(10) ~ 2.303.
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 2, U_history, d_history, 1, 1.0, 1.0, U, 2, 10.0, 1.0,
                                  0.9, 0, 3.0, 1e10)
    assert is_accepted, "a 10x change in G must be accepted at G_max=3.0"


def test_accept_ensemble_rmse_ratio_exceeds_max():
    # normal_error changes by a factor of 10 (RMSE by sqrt(10) ~ 3.162, log ~ 1.151): rejected
    # at RMSE_change_max=1.0.
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 2, U_history, d_history, 1, 1.0, 1.0, U, 2, 1.0, 10.0,
                                  0.9, 0, 1e10, 1.0)
    assert not is_accepted, "a 10x change in normal_error must be rejected at RMSE_change_max=1.0"


def test_accept_ensemble_rmse_ratio_within_max():
    # Same 10x change in normal_error, accepted once RMSE_change_max is raised past
    # |log(sqrt(10))| ~ 1.151.
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 2, U_history, d_history, 1, 1.0, 1.0, U, 2, 1.0, 10.0,
                                  0.9, 0, 1e10, 1.2)
    assert is_accepted, "a 10x change in normal_error must be accepted at RMSE_change_max=1.2"


def test_accept_ensemble_nonpositive_g():
    U, U_history, d_history = _identity_history(3, 1, 2)

    assert_error(lambda: accept_ensemble(U, 2, U_history, d_history, 1, 0.0, 1.0, U, 2, 1.0, 1.0,
                                         0.9, 0, 1e10, 1e10),
                 "Expected error for G_t <= 0", ERR_INVALID_INPUT)


def test_accept_ensemble_nonpositive_normal_error():
    # Negative normal_error is physically impossible (a sum of eigenvalues) and must be
    # rejected by validation -- unlike G_t, zero itself is valid (see the next test).
    U, U_history, d_history = _identity_history(3, 1, 2)

    assert_error(lambda: accept_ensemble(U, 2, U_history, d_history, 1, 1.0, -1.0, U, 2, 1.0, 1.0,
                                         0.9, 0, 1e10, 1e10),
                 "Expected error for normal_error_t < 0", ERR_INVALID_INPUT)


def test_accept_ensemble_zero_normal_error_is_accepted():
    # Zero normal_error at both t and t+1 (a perfectly flat/collinear ensemble) must not crash
    # or spuriously reject via a log(0/0) -- the +epsilon guard inside the RMSE ratio keeps it
    # well-defined and accepted.
    U, U_history, d_history = _identity_history(3, 1, 2)

    is_accepted = accept_ensemble(U, 2, U_history, d_history, 1, 1.0, 0.0, U, 2, 1.0, 0.0,
                                  0.9, 0, 1e10, 1e10)
    assert is_accepted, "zero normal_error at both t and t+1 must be accepted, not NaN-rejected"


def test_accept_ensemble_d_first_out_of_range():
    U, U_history, d_history = _identity_history(3, 1, 2)

    assert_error(lambda: accept_ensemble(U, 4, U_history, d_history, 1, 1.0, 1.0, U, 2, 1.0, 1.0,
                                         0.9, 0, 1e10, 1e10),
                 "Expected error for d_first > n_dimensions", ERR_INVALID_INPUT)


def test_accept_ensemble_chordal_frac_out_of_range():
    U, U_history, d_history = _identity_history(3, 1, 2)

    assert_error(lambda: accept_ensemble(U, 2, U_history, d_history, 1, 1.0, 1.0, U, 2, 1.0, 1.0,
                                         1.5, 0, 1e10, 1e10),
                 "Expected error for chordal_dist_max_as_prcnt_of_range > 1", ERR_INVALID_INPUT)


def test_accept_ensemble_history_len_out_of_range():
    U, U_history, d_history = _identity_history(3, 1, 2)

    assert_error(lambda: accept_ensemble(U, 2, U_history, d_history, 2, 1.0, 1.0, U, 2, 1.0, 1.0,
                                         0.9, 0, 1e10, 1e10),
                 "Expected error for history_len > o", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
