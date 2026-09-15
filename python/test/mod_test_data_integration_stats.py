"""
Python test suite for data integration stats functions (the K-study GSL permutation test) in
tensoromics. Uses tensor_omics.py wrapper function (mirrors Fortran test suite).
"""
import numpy as np
import sys
import os

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from test_helpers import run_all_tests
from tensor_omics import gjct_permutation_test

TOL = 1e-12


def test_gjct_permutation_test_basic():
    """p_values must land in [0, 1] for a simple 2-study, 2-point scenario."""
    n_bins, n_points, n_studies, n_permutations = 3, 2, 2, 20

    mean_pmf_counts = np.zeros((n_bins, n_points), dtype=np.int32, order="F")
    mean_pmf_counts[:, 0] = [5, 3, 2]
    mean_pmf_counts[:, 1] = [2, 2, 6]
    mean_pmf_included_n_reps = np.array([10, 10], dtype=np.int32)
    included_n_reps = np.zeros((n_points, n_studies), dtype=np.int32, order="F")
    included_n_reps[:, 0] = [4, 4]
    included_n_reps[:, 1] = [6, 6]
    mean_pmf = mean_pmf_counts.astype(np.float64) / 10.0
    global_jsd_observed = np.array([0.3, 0.3], dtype=np.float64)

    p_values = gjct_permutation_test(
        n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
        included_n_reps, global_jsd_observed, random_seed=123
    )

    assert len(p_values) == n_studies
    assert np.all(p_values >= 0.0) and np.all(p_values <= 1.0)


def test_gjct_permutation_test_seeded_reproducibility():
    """Same random_seed twice must give identical p_values -- mirrors the Fortran test of the
    same name (test/mod_test_data_integration_js_comp_test.F90)."""
    n_bins, n_points, n_studies, n_permutations = 3, 2, 2, 20

    mean_pmf_counts = np.zeros((n_bins, n_points), dtype=np.int32, order="F")
    mean_pmf_counts[:, 0] = [5, 3, 2]
    mean_pmf_counts[:, 1] = [2, 2, 6]
    mean_pmf_included_n_reps = np.array([10, 10], dtype=np.int32)
    included_n_reps = np.zeros((n_points, n_studies), dtype=np.int32, order="F")
    included_n_reps[:, 0] = [4, 4]
    included_n_reps[:, 1] = [6, 6]
    mean_pmf = mean_pmf_counts.astype(np.float64) / 10.0
    global_jsd_observed = np.array([0.3, 0.3], dtype=np.float64)

    p_first = gjct_permutation_test(
        n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
        included_n_reps, global_jsd_observed, random_seed=123
    )
    p_second = gjct_permutation_test(
        n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
        included_n_reps, global_jsd_observed, random_seed=123
    )

    np.testing.assert_allclose(p_first, p_second, atol=TOL)


if __name__ == '__main__':
    run_all_tests(globals().values())
