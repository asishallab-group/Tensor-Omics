"""
Python test suite for ensemble_reconciliation (tox_shape_truthful_clustering_reconciliation),
mirroring test/mod_test_shape_truthful_clustering_reconciliation.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import ensemble_reconciliation
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT, ERR_SIZE_MISMATCH


def _fixture():
    """N=14 vectors, 6 ensembles. E1={1,2,3,4}, E2={3,4,5,6}, E3={5,6,7,8}: a chain, E1-E2
    and E2-E3 each intersect at 2 members (JSI = 1/3), E1-E3 do not intersect at all.
    E4={9,10}: isolated. E5={11,12,13}, E6={12,13,14}: a separate pair (JSI = 1/2)."""
    m = np.zeros((14, 6), dtype=np.bool_, order='F')
    m[0:4, 0] = True
    m[2:6, 1] = True
    m[4:8, 2] = True
    m[8:10, 3] = True
    m[10:13, 4] = True
    m[11:14, 5] = True
    return m


def test_report_mode_no_transitive_merge():
    result = ensemble_reconciliation(_fixture(), 2, mode='report', report_jsi=True)

    assert result['n_super_ensembles'] == 3
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2])
    assert np.array_equal(result['super_ensembles'][:, 1], [2, 3])
    assert np.array_equal(result['super_ensembles'][:, 2], [5, 6])
    assert abs(result['super_ensembles_JSI'][0, 0] - 1.0 / 3.0) < 1e-9
    assert abs(result['super_ensembles_JSI'][0, 1] - 1.0 / 3.0) < 1e-9
    assert abs(result['super_ensembles_JSI'][0, 2] - 0.5) < 1e-9


def test_merge_any_transitive():
    result = ensemble_reconciliation(_fixture(), 3, mode='merge_any', report_jsi=True)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])
    assert abs(result['super_ensembles_JSI'][0, 0] - 1.0 / 3.0) < 1e-9
    assert abs(result['super_ensembles_JSI'][1, 0] - 1.0 / 3.0) < 1e-9
    assert np.array_equal(result['super_ensembles'][:, 1], [5, 6, 0])
    assert abs(result['super_ensembles_JSI'][0, 1] - 0.5) < 1e-9


def test_merge_jsi_threshold_excludes_weak_chain():
    # min_jsi=0.4 excludes the 1-2-3 chain (JSI ~= 0.333) but not the 5-6 pair (JSI = 0.5).
    result = ensemble_reconciliation(_fixture(), 3, mode='merge_jsi', min_jsi=0.4)

    assert result['n_super_ensembles'] == 1
    assert np.array_equal(result['super_ensembles'][:, 0], [5, 6, 0])


def test_merge_jsi_threshold_includes_all():
    result = ensemble_reconciliation(_fixture(), 3, mode='merge_jsi', min_jsi=0.3)

    assert result['n_super_ensembles'] == 2
    assert np.array_equal(result['super_ensembles'][:, 0], [1, 2, 3])


def test_jsi_not_computed_unless_requested():
    result = ensemble_reconciliation(_fixture(), 3, mode='merge_any')

    assert result['n_super_ensembles'] == 2
    assert np.all(result['super_ensembles_JSI'] == 0.0)


def test_group_exceeds_max_group_size():
    assert_error(lambda: ensemble_reconciliation(_fixture(), 2, mode='merge_any'),
                 "a 3-member group must not fit in max_group_size=2", ERR_SIZE_MISMATCH)


def test_n_ensembles_too_small():
    m = np.zeros((14, 1), dtype=np.bool_, order='F')
    assert_error(lambda: ensemble_reconciliation(m, 2),
                 "Expected error for n_ensembles=1", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
