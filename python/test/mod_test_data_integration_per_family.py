"""
Comprehensive Python test suite for data integration per family functions in tensoromics.
Uses tensoromics_functions.py wrapper function (mirrors Fortran test suite)
"""
import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from test_helpers import run_all_tests
from tensoromics_functions import (
    fjct_compute_jsd,
    fjct_compute_jsd_expert,
    fjct_compute_contribution_scores,
)

TOL = 1e-12


def test_fjct():
    # Parameters
    n_reps_S1 = 3
    n_reps_S2 = 4
    n_neighbors = 5
    n_points = 3
    k_families = 2

    n_bins = 4
    n_genes_S1 = 10
    n_genes_S2 = 10

    included_families = np.array([1, 2], dtype=np.int32)

    gene_to_family_S1 = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 2], dtype=np.int32)
    gene_to_family_S2 = np.array([1, 0, 0, 0, 0, 0, 0, 0, 3, 0], dtype=np.int32)

    # Allocate arrays
    neighborhood_genes_S1 = np.zeros((n_neighbors, n_points), dtype=np.int32)
    neighborhood_genes_S2 = np.zeros((n_neighbors, n_points), dtype=np.int32)

    neighborhood_residuals_S1 = np.zeros((n_reps_S1, n_neighbors, n_points), dtype=np.float64, order="F")
    neighborhood_residuals_S2 = np.zeros((n_reps_S2, n_neighbors, n_points), dtype=np.float64, order="F")

    # ============================================================
    # Test 1 — Simple symmetric case, no NaNs
    # ============================================================

    shared_residual_range = 2.0

    # S1 residuals
    neighborhood_residuals_S1[:, 1-1, 1-1] = np.array([-2.0, -0.5, 0.2])
    neighborhood_residuals_S1[:, 1-1, 3-1] = np.array([2.5, -3.0, 1.2])  # clamped to [2, -2, 1.2]

    neighborhood_genes_S1[:, 1-1] = np.array([1, 2, 3, 4, 5])
    neighborhood_genes_S1[:, 2-1] = np.array([2, 3, 4, 5, 6])
    neighborhood_genes_S1[:, 3-1] = np.array([10, 2, 3, 4, 5])

    # S2 residuals
    neighborhood_residuals_S2[:, 1-1, 1-1] = np.array([-2.0, np.nan, -0.5, 0.2])
    neighborhood_residuals_S2[:, 3-1, 3-1] = np.array([0.4, -0.1, 0.0, np.nan])

    neighborhood_genes_S2[:, 1-1] = np.array([1, 2, 3, 4, 5])
    neighborhood_genes_S2[:, 2-1] = np.array([2, 3, 4, 5, 6])
    neighborhood_genes_S2[:, 3-1] = np.array([9, 2, 3, 4, 5])

    # ============================================================
    # Family 1
    # ============================================================

    family_idx = 1

    mask_S1 = np.full((n_neighbors, n_points), False, order="F")
    mask_S2 = np.full((n_neighbors, n_points), False, order="F")
    mask_S1[0, 0] = True
    mask_S2[0, 0] = True

    expected_included_n_reps_S1 = np.array([3, 0, 0], dtype=np.int32)
    expected_included_n_reps_S2 = np.array([3, 0, 0], dtype=np.int32)
    expected_total_included_n_reps = 6
    expected_js_divergences = np.zeros(n_points, dtype=np.float64)
    expected_weights = np.array([1.0, 0.0, 0.0], dtype=np.float64)

    res = fjct_compute_jsd(
        family_idx,
        gene_to_family_S1,
        gene_to_family_S2,
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighborhood_genes_S1,
        neighborhood_genes_S2,
        n_bins,
        shared_residual_range
    )

    assert res["ierr"] == 0
    assert res["total_included_n_reps"] == expected_total_included_n_reps
    np.testing.assert_array_equal(res["included_n_reps_S1"], expected_included_n_reps_S1)
    np.testing.assert_array_equal(res["included_n_reps_S2"], expected_included_n_reps_S2)
    np.testing.assert_allclose(res["js_divergences"], expected_js_divergences, atol=TOL)
    np.testing.assert_allclose(res["weights"], expected_weights, atol=TOL)

    expected_global_jsd = expected_weights[2] * expected_js_divergences[2]
    assert abs(res["global_js_divergence"] - expected_global_jsd) < TOL

    res_expert = fjct_compute_jsd_expert(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        mask_S1,
        mask_S2,
        n_bins,
        shared_residual_range
    )

    assert res["ierr"] == 0
    assert res["total_included_n_reps"] == expected_total_included_n_reps
    np.testing.assert_array_equal(res["included_n_reps_S1"], expected_included_n_reps_S1)
    np.testing.assert_array_equal(res["included_n_reps_S2"], expected_included_n_reps_S2)
    np.testing.assert_allclose(res["js_divergences"], expected_js_divergences, atol=TOL)
    np.testing.assert_allclose(res["weights"], expected_weights, atol=TOL)

    expected_global_jsd = expected_weights[2] * expected_js_divergences[2]
    assert abs(res["global_js_divergence"] - expected_global_jsd) < TOL

    # ============================================================
    # Family 2
    # ============================================================

    family_idx = 2

    mask_S1 = np.full((n_neighbors, n_points), False, order="F")
    mask_S2 = np.full((n_neighbors, n_points), False, order="F")
    mask_S1[0, 2] = True

    expected_included_n_reps_S1 = np.array([0, 0, 3], dtype=np.int32)
    expected_included_n_reps_S2 = np.array([0, 0, 0], dtype=np.int32)
    expected_total_included_n_reps = 3

    expected_js_divergences = np.array([
        0.0,
        0.0,
        0.5 * np.log(2.0)
    ], dtype=np.float64)

    expected_weights = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    res = fjct_compute_jsd(
        family_idx,
        gene_to_family_S1,
        gene_to_family_S2,
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighborhood_genes_S1,
        neighborhood_genes_S2,
        n_bins,
        shared_residual_range
    )

    assert res["ierr"] == 0
    assert res["total_included_n_reps"] == expected_total_included_n_reps
    np.testing.assert_array_equal(res["included_n_reps_S1"], expected_included_n_reps_S1)
    np.testing.assert_array_equal(res["included_n_reps_S2"], expected_included_n_reps_S2)
    np.testing.assert_allclose(res["js_divergences"], expected_js_divergences, atol=TOL)
    np.testing.assert_allclose(res["weights"], expected_weights, atol=TOL)

    expected_global_jsd = expected_weights[2] * expected_js_divergences[2]
    assert abs(res["global_js_divergence"] - expected_global_jsd) < TOL

    res_expert = fjct_compute_jsd_expert(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        mask_S1,
        mask_S2,
        n_bins,
        shared_residual_range
    )

    assert res["ierr"] == 0
    assert res["total_included_n_reps"] == expected_total_included_n_reps
    np.testing.assert_array_equal(res["included_n_reps_S1"], expected_included_n_reps_S1)
    np.testing.assert_array_equal(res["included_n_reps_S2"], expected_included_n_reps_S2)
    np.testing.assert_allclose(res["js_divergences"], expected_js_divergences, atol=TOL)
    np.testing.assert_allclose(res["weights"], expected_weights, atol=TOL)

    expected_global_jsd = expected_weights[2] * expected_js_divergences[2]
    assert abs(res["global_js_divergence"] - expected_global_jsd) < TOL

    # ============================================================
    # Contribution scores
    # ============================================================

    total_included_n_reps = np.array([6, 3], dtype=np.int32)
    global_js_divergence = np.array([
        0.0,
        0.5 * np.log(2.0)
    ], dtype=np.float64)

    expected_support_weights = np.array([
        6.0 / 9.0,
        3.0 / 9.0
    ], dtype=np.float64)

    expected_contribution_scores = expected_support_weights * global_js_divergence

    res2 = fjct_compute_contribution_scores(
        global_js_divergence,
        total_included_n_reps
    )

    assert res2["ierr"] == 0
    np.testing.assert_allclose(res2["support_weights"], expected_support_weights, atol=TOL)
    np.testing.assert_allclose(res2["contribution_scores"], expected_contribution_scores, atol=TOL)


if __name__ == '__main__':
    run_all_tests(globals().values())
