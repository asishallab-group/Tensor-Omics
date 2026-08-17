"""
Python test suite for ensemble_reconciliation / merge_to_super_ensembles
(tox_shape_truthful_clustering_reconciliation), mirroring
test/mod_test_shape_truthful_clustering_reconciliation.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import ensemble_reconciliation, merge_to_super_ensembles
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT, ERR_SIZE_MISMATCH

# ensemble_stop_reason values, 1-indexed -- tox_shape_truthful_clustering_impl's own
# STOP_REASON_MAX_SIZE(1)/STOP_REASON_REJECTED_AFTER_STABLE(2)/STOP_REASON_REJECTED_IMMEDIATELY(3)/
# STOP_REASON_FIXED_POINT(4). Not exported as named constants by any binding (see the Fortran
# kernel's own doc comment on why), so used as plain integers here, matching every other
# consumer of this array.
STOP_REASON_MAX_SIZE = 1
STOP_REASON_REJECTED_AFTER_STABLE = 2
STOP_REASON_REJECTED_IMMEDIATELY = 3
STOP_REASON_FIXED_POINT = 4


def _fixture():
    """N=14 vectors, 6 ensembles. E1={1,2,3,4}, E2={3,4,5,6}, E3={5,6,7,8}: a chain, E1-E2
    and E2-E3 each intersect at 2 members, both ensembles size 4, so OC = 2/min(4,4) = 0.5.
    E1-E3 do not intersect at all. E4={9,10}: isolated. E5={11,12,13}, E6={12,13,14}: a
    separate pair, intersecting at 2 members, both ensembles size 3, so OC = 2/min(3,3) = 2/3."""
    m = np.zeros((14, 6), dtype=np.bool_, order='F')
    m[0:4, 0] = True
    m[2:6, 1] = True
    m[4:8, 2] = True
    m[8:10, 3] = True
    m[10:13, 4] = True
    m[11:14, 5] = True
    return m


def _stop_reasons():
    """Every ensemble STOP_REASON_FIXED_POINT by default -- a neutral choice that changes
    nothing about any pre-existing test below when no allowed_stop_reasons filter is applied."""
    return np.full(6, STOP_REASON_FIXED_POINT, dtype=np.int32)


def _history():
    """Minimal, uniform D=2/o=1 history for 6 ensembles -- ensemble_reconciliation's new
    required history arguments, needed only to drive the dimension/variance-explained filters,
    which none of the tests below (other than the dedicated d_min/d_max/var_explained_min ones)
    actually exercise. Returns (U, d, S, mu, G, k, accepted)."""
    n_e = 6
    U = np.zeros((2, 2, 1, n_e), dtype=np.float64, order='F')
    d = np.zeros((1, n_e), dtype=np.int32, order='F')
    S = np.zeros((2, 1, n_e), dtype=np.float64, order='F')
    mu = np.zeros((2, 1, n_e), dtype=np.float64, order='F')
    G = np.zeros((1, n_e), dtype=np.float64, order='F')
    k = np.full((1, n_e), 2, dtype=np.int32, order='F')
    accepted = np.ones((1, n_e), dtype=np.bool_, order='F')
    return U, d, S, mu, G, k, accepted


def _reconcile(ensemble_masks, ensemble_stop_reason, max_group_size, history=None, **kwargs):
    """ensemble_reconciliation, with the history-array positional block filled in from
    `_history()` unless the caller supplies its own via `history=(U, d, S, mu, G, k, accepted)`."""
    U, d, S, mu, G, k, accepted = history if history is not None else _history()
    return ensemble_reconciliation(ensemble_masks, ensemble_stop_reason, U, d, S, mu, G, k, accepted,
                                   max_group_size, **kwargs)


def test_report_mode_no_transitive_merge():
    result = _reconcile(_fixture(), _stop_reasons(), 2, mode='report', report_overlap_coefficient=True)

    assert result['n_super_ensembles'] == 3
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2])
    assert np.array_equal(result['super_ensembles'][:, 1], [2, 3])
    assert np.array_equal(result['super_ensembles'][:, 2], [5, 6])
    assert abs(result['super_ensembles_overlap_coefficient'][0, 0] - 0.5) < 1e-9
    assert abs(result['super_ensembles_overlap_coefficient'][0, 1] - 0.5) < 1e-9
    assert abs(result['super_ensembles_overlap_coefficient'][0, 2] - 2.0 / 3.0) < 1e-9
    assert np.all(result['eligible'])


def test_merge_any_transitive():
    result = _reconcile(_fixture(), _stop_reasons(), 3, mode='merge_any', report_overlap_coefficient=True)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])
    assert abs(result['super_ensembles_overlap_coefficient'][0, 0] - 0.5) < 1e-9
    assert abs(result['super_ensembles_overlap_coefficient'][1, 0] - 0.5) < 1e-9
    assert np.array_equal(result['super_ensembles'][:, 1], [5, 6, 0])
    assert abs(result['super_ensembles_overlap_coefficient'][0, 1] - 2.0 / 3.0) < 1e-9


def test_merge_overlap_coefficient_threshold_excludes_weak_chain():
    # min_overlap_coefficient=0.6 excludes the 1-2-3 chain (OC = 0.5) but not the 5-6 pair (OC = 2/3).
    result = _reconcile(_fixture(), _stop_reasons(), 3, mode='merge_overlap_coefficient',
                        min_overlap_coefficient=0.6)

    assert result['n_super_ensembles'] == 1
    assert np.array_equal(result['super_ensembles'][:, 0], [5, 6, 0])


def test_merge_overlap_coefficient_threshold_includes_all():
    result = _reconcile(_fixture(), _stop_reasons(), 3, mode='merge_overlap_coefficient',
                        min_overlap_coefficient=0.4)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])


def test_overlap_coefficient_not_computed_unless_requested():
    result = _reconcile(_fixture(), _stop_reasons(), 3, mode='merge_any')

    assert result['n_super_ensembles'] == 2
    assert np.all(result['super_ensembles_overlap_coefficient'] == 0.0)


def test_group_exceeds_max_group_size():
    assert_error(lambda: _reconcile(_fixture(), _stop_reasons(), 2, mode='merge_any'),
                 "a 3-member group must not fit in max_group_size=2", ERR_SIZE_MISMATCH)


def test_n_ensembles_too_small():
    m = np.zeros((14, 1), dtype=np.bool_, order='F')
    sr = np.full(1, STOP_REASON_FIXED_POINT, dtype=np.int32)
    n_e = 1
    U = np.zeros((2, 2, 1, n_e), dtype=np.float64, order='F')
    d = np.zeros((1, n_e), dtype=np.int32, order='F')
    S = np.zeros((2, 1, n_e), dtype=np.float64, order='F')
    mu = np.zeros((2, 1, n_e), dtype=np.float64, order='F')
    G = np.zeros((1, n_e), dtype=np.float64, order='F')
    k = np.full((1, n_e), 2, dtype=np.int32, order='F')
    accepted = np.ones((1, n_e), dtype=np.bool_, order='F')
    assert_error(lambda: ensemble_reconciliation(m, sr, U, d, S, mu, G, k, accepted, 2),
                 "Expected error for n_ensembles=1", ERR_INVALID_INPUT)


def test_allowed_stop_reasons_excludes_pair_report_mode():
    """Ensemble 5 is STOP_REASON_REJECTED_IMMEDIATELY (rest STOP_REASON_FIXED_POINT); excluding
    STOP_REASON_REJECTED_IMMEDIATELY must drop the (5,6) pair from report mode's output, leaving
    only the untouched 1-2-3 chain's two pairs."""
    stop_reasons = _stop_reasons()
    stop_reasons[4] = STOP_REASON_REJECTED_IMMEDIATELY
    allowed = np.array([True, True, False, True], dtype=np.bool_)

    result = _reconcile(_fixture(), stop_reasons, 2, mode='report', allowed_stop_reasons=allowed)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2])
    assert np.array_equal(result['super_ensembles'][:, 1], [2, 3])
    assert not result['eligible'][4]
    assert not result['eligible_by_stop_condition'][4]


def test_allowed_stop_reasons_absent_matches_all_true():
    """Omitting allowed_stop_reasons (None), even with a genuinely mixed set of Stop Conditions
    across ensembles, must behave identically to no filtering at all."""
    stop_reasons = _stop_reasons()
    stop_reasons[0] = STOP_REASON_MAX_SIZE
    stop_reasons[3] = STOP_REASON_REJECTED_AFTER_STABLE
    stop_reasons[5] = STOP_REASON_REJECTED_IMMEDIATELY

    result = _reconcile(_fixture(), stop_reasons, 3, mode='merge_any', report_overlap_coefficient=True)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])
    assert np.array_equal(result['super_ensembles'][:, 1], [5, 6, 0])
    assert np.all(result['eligible'])


def test_allowed_stop_reasons_breaks_transitive_chain():
    """Ensemble 2 (the sole bridge of the 1-2-3 chain) is STOP_REASON_REJECTED_IMMEDIATELY;
    excluding it must fully break the chain (1 and 3 do not intersect directly), leaving only
    the untouched {5,6} group."""
    stop_reasons = _stop_reasons()
    stop_reasons[1] = STOP_REASON_REJECTED_IMMEDIATELY
    allowed = np.array([True, True, False, True], dtype=np.bool_)

    result = _reconcile(_fixture(), stop_reasons, 3, mode='merge_any', allowed_stop_reasons=allowed)

    assert result['n_super_ensembles'] == 1
    assert np.array_equal(result['super_ensembles'][:, 0], [5, 6, 0])


def test_allowed_stop_reasons_noop_when_no_ensemble_matches():
    """Every ensemble is STOP_REASON_FIXED_POINT; excluding STOP_REASON_REJECTED_AFTER_STABLE
    (matching none of them) must be a true no-op."""
    allowed = np.array([True, False, True, True], dtype=np.bool_)

    result = _reconcile(_fixture(), _stop_reasons(), 3, mode='merge_any', allowed_stop_reasons=allowed)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])
    assert np.array_equal(result['super_ensembles'][:, 1], [5, 6, 0])


def test_dimension_filter_excludes_pair():
    """Ensemble 5's final intrinsic dimension (2) exceeds d_max=1 (the rest are d=1); excluding
    it must drop the (5,6) pair from report mode's output, mirroring the stop-condition-filter
    test above but through the dimension filter."""
    U, d, S, mu, G, k, accepted = _history()
    d[0, :] = 1
    d[0, 4] = 2  # ensemble 5 (0-indexed 4)

    result = _reconcile(_fixture(), _stop_reasons(), 2, history=(U, d, S, mu, G, k, accepted),
                        mode='report', d_max=1)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2])
    assert np.array_equal(result['super_ensembles'][:, 1], [2, 3])
    assert not result['eligible_by_dimension'][4]
    assert result['eligible_by_dimension'][0]
    assert not result['eligible'][4]


def test_var_explained_filter_excludes_pair():
    """Ensemble 5's final variance explained (eigenvalues [1,100], d=1 -> 1/101 ~ 0.0099) falls
    far short of var_explained_min=0.5; the rest ([100,1] -> 100/101 ~ 0.99) comfortably clear
    it. Excluding ensemble 5 must drop the (5,6) pair."""
    U, d, S, mu, G, k, accepted = _history()
    d[0, :] = 1
    for e in range(6):
        S[:, 0, e] = [10.0, 1.0]  # eigenvalues [100,1], ve = 100/101
    S[:, 0, 4] = [1.0, 10.0]  # eigenvalues [1,100], ve = 1/101

    result = _reconcile(_fixture(), _stop_reasons(), 2, history=(U, d, S, mu, G, k, accepted),
                        mode='report', var_explained_min=0.5)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2])
    assert np.array_equal(result['super_ensembles'][:, 1], [2, 3])
    assert not result['eligible_by_var_explained'][4]
    assert result['eligible_by_var_explained'][0]


def test_merge_to_super_ensembles_all_eligible_matches_merge_any():
    """Direct test of the newly-independent merge_to_super_ensembles: an all-True eligible mask
    must reproduce test_merge_any_transitive's own result, with no history array involved."""
    eligible = np.ones(6, dtype=np.bool_)

    result = merge_to_super_ensembles(_fixture(), eligible, 3, mode='merge_any', report_overlap_coefficient=True)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])
    assert np.array_equal(result['super_ensembles'][:, 1], [5, 6, 0])


def test_merge_to_super_ensembles_excludes_ineligible_ensemble():
    """Ensemble 5 marked ineligible via a hand-constructed eligible mask: the (5,6) edge must
    never be considered, leaving only the untouched 1-2-3 chain."""
    eligible = np.ones(6, dtype=np.bool_)
    eligible[4] = False

    result = merge_to_super_ensembles(_fixture(), eligible, 3, mode='merge_any', report_overlap_coefficient=True)

    assert result['n_super_ensembles'] == 1
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])


if __name__ == "__main__":
    run_all_tests(globals().values())
