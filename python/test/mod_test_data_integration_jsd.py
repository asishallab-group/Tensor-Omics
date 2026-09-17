"""
Comprehensive Python test suite for data integration and Jensen–Shannon divergence functions in tensoromics.
Uses tensor_omics.py wrapper function (mirrors Fortran test suite)
"""
import numpy as np
import sys
import os

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from test_helpers import run_all_tests, assert_error
from tensor_omics import (
    determine_shared_residual_range,
    determine_study_shared_residual_range,
    determine_shared_residual_range_expert,
    build_residual_histograms,
    compute_divergence_per_reference_point,
    compute_weighted_global_divergence,
    calc_pmf,
    determine_all_studies_shared_residual_range,
)
from tensor_omics.error_handling import ERR_INVALID_INPUT

# The filtered variant is gone: filtering is now the optional `neighbor_mask` argument of the
# base routine. gjct_permutation_test itself was replaced by a K-study, consensus-based version
# (see tox_data_integration_js_comp_test) -- the old 2-study test that lived in this file was
# removed rather than adapted, since the signature is unrelated. Its own direct test coverage now
# lives in mod_test_data_integration_stats.py.
build_residual_histograms_filtered = build_residual_histograms


TOL = 1e-12


def test_tox_determine_shared_residual_range():
    """Test determine_shared_residual_range (the plain, auto-sorting tier) directly, by cross-
    checking it against the already-tested _expert tier given the same pool pre-sorted."""
    pool = np.array([1.0, 5.0, 3.0, 8.0, 2.0, 7.0, 4.0], dtype=np.float64)
    perm = (np.argsort(pool, kind="mergesort").astype(np.int32) + 1)

    R_plain = determine_shared_residual_range(pool, 0.95)
    R_expert = determine_shared_residual_range_expert(pool, perm, 0.95)
    assert abs(R_plain - R_expert) < TOL, f"plain vs expert mismatch: {R_plain} vs {R_expert}"

    # Default residual_range_quantile (0.95)
    R_plain_default = determine_shared_residual_range(pool)
    assert abs(R_plain_default - R_expert) < TOL, "plain tier's default quantile should also be 0.95"


def test_tox_determine_study_shared_residual_range():

    # ============================================================
    # Test 1 — Basic correctness with simple values
    # ============================================================
    S1 = np.zeros((4, 2, 2), dtype=np.float64, order="F")
    S1[:, 0, 0] = [1,  2,  3, 4]
    S1[:, 1, 0] = [5,  6, -7, 8]
    S1[:, 0, 1] = [9, 10, 11, 12]
    S1[:, 1, 1] = [1, 1, 1, 1]

    S2 = np.zeros((3, 2, 2), dtype=np.float64, order="F")
    S2[:, 0, 0] = [2, -4,  6]
    S2[:, 1, 0] = [8,  1,  3]
    S2[:, 0, 1] = [5,  7,  9]
    S2[:, 1, 1] = [0,  1,  2]

    R = determine_study_shared_residual_range(S1, S2, 0.95)
    assert abs(R - 10.65) < TOL, f"Test 1 failed: expected 10.65, got {R}"

    # ============================================================
    # Test 2 — Custom quantile (median, q=0.5)
    # ============================================================
    R = determine_study_shared_residual_range(S1, S2, 0.5)
    assert abs(R - 4.0) < TOL, f"Test 2 failed: expected 4.0, got {R}"

    # ============================================================
    # Test 3 — Quantile < 0 → error
    # ============================================================
    assert_error(lambda: determine_study_shared_residual_range(S1, S2, -1.0), "Test 3 failed: expected ERR_INVALID_INPUT", ERR_INVALID_INPUT)

    # ============================================================
    # Test 4 — Quantile > 1 → error
    # ============================================================
    assert_error(lambda: determine_study_shared_residual_range(S1, S2, 1.5), "Test 4 failed: expected ERR_INVALID_INPUT", ERR_INVALID_INPUT)

    # ============================================================
    # Test 5 — NaNs must be ignored
    # ============================================================
    S1 = np.zeros((4, 2, 2), dtype=np.float64, order="F")
    S1[:, 0, 0] = [-1,  2,   3, 4]
    S1[:, 1, 0] = [5, 6, -7, 8]
    S1[:, 0, 1] = [9, 10, -11, 12]
    S1[:, 1, 1] = [1, 1, 1, 1]

    S2 = np.zeros((3, 2, 2), dtype=np.float64, order="F")
    S2[:, 0, 0] = [-1,  2, 3]
    S2[:, 1, 0] = [4, 5, 6]
    S2[:, 0, 1] = [-7, 8, 9]
    S2[:, 1, 1] = [10, -11, 12]

    S1[0, 0, 0] = np.nan
    S2[2, 1, 1] = np.nan

    R = determine_study_shared_residual_range(S1, S2, 0.95)
    assert abs(R - 11.0) < TOL, f"Test 5 failed: expected 11.0, got {R}"

    # ============================================================
    # Test 6 — All zeros
    # ============================================================
    S1[:, :, :] = 0
    S2[:, :, :] = 0

    R = determine_study_shared_residual_range(S1, S2, 0.95)
    assert abs(R - 0.0) < TOL, f"Test 6 failed: expected 0.0, got {R}"

    # ============================================================
    # Test 7 — Single residual (1×1)
    # ============================================================
    S1 = np.array([[[3.0]]], dtype=np.float64, order="F")
    S2 = np.array([[[-4.0]]], dtype=np.float64, order="F")

    R = determine_study_shared_residual_range(S1, S2, 0.95)
    assert abs(R - 3.95) < TOL, f"Test 7 failed: expected 3.95, got {R}"


def test_tox_determine_shared_residual_range_expert():
    # Helper to build pool + permutation
    def make_pool(S1, S2):
        pool = np.abs(np.concatenate([S1.ravel(order="F"), S2.ravel(order="F")]))
        perm = np.argsort(pool) + 1
        return pool, perm

    # ============================================================
    # Test 1 — Basic correctness with simple values
    # ============================================================
    S1 = np.array([
        [ 1,  2,  3],
        [ 4,  5,  6],
        [-7,  8,  9],
        [10, 11, 12],
    ], dtype=np.float64, order="F")

    S2 = np.array([
        [ 2, -4,  6],
        [ 8,  1,  3],
        [ 5,  7,  9],
        [ 0,  1,  2],
    ], dtype=np.float64, order="F")

    pool, perm = make_pool(S1, S2)
    R = determine_shared_residual_range_expert(pool, perm)
    assert abs(R - 10.85) < TOL, f"Test 1 failed: expected 10.85, got {R}"

    # ============================================================
    # Test 2 — Custom quantile (median, q=0.5)
    # ============================================================
    R = determine_shared_residual_range_expert(pool, perm, 0.5)
    assert abs(R - 5.0) < TOL, f"Test 2 failed: expected 5.0, got {R}"

    # ============================================================
    # Test 3 — Quantile < 0 → error
    # ============================================================
    assert_error(lambda: determine_shared_residual_range_expert(pool, perm, -1.0), "Test 3 failed: expected ERR_INVALID_INPUT", ERR_INVALID_INPUT)

    # ============================================================
    # Test 4 — Quantile > 1 → error
    # ============================================================
    assert_error(lambda: determine_shared_residual_range_expert(pool, perm, 1.5), "Test 4 failed: expected ERR_INVALID_INPUT", ERR_INVALID_INPUT)

    # ============================================================
    # Test 5 — NaNs must be ignored
    # ============================================================
    S1 = np.array([
        [np.nan,  2,   3],
        [ 4,      5,   6],
        [-7,      8, -11],
        [ 9,     10,  12],
    ], dtype=np.float64, order="F")

    S2 = S1.copy(order="F")
    S2[3, 2] = np.nan

    pool, perm = make_pool(S1, S2)
    R = determine_shared_residual_range_expert(pool, perm, 0.95)
    assert abs(R - 11.0) < TOL, f"Test 5 failed: expected 11.0, got {R}"

    # ============================================================
    # Test 6 — All zeros
    # ============================================================
    S1 = np.zeros((4, 3), dtype=np.float64, order="F")
    S2 = np.zeros((4, 3), dtype=np.float64, order="F")

    pool, perm = make_pool(S1, S2)
    R = determine_shared_residual_range_expert(pool, perm, 0.95)
    assert abs(R - 0.0) < TOL, f"Test 6 failed: expected 0.0, got {R}"

    # ============================================================
    # Test 7 — Single residual (1×1)
    # ============================================================
    S1 = np.array([[3.0]], dtype=np.float64, order="F")
    S2 = np.array([[-4.0]], dtype=np.float64, order="F")

    pool, perm = make_pool(S1, S2)
    R = determine_shared_residual_range_expert(pool, perm, 0.95)
    assert abs(R - 3.95) < TOL, f"Test 7 failed: expected 3.95, got {R}"


def test_tox_build_residual_histograms():
    n_reps = 3
    n_neighbors = 2
    n_points = 3
    n_bins = 4
    n_bins_per_point = np.full(n_points, n_bins, dtype=np.int32)
    R = 2.0

    # ============================================================
    # Test 1 — Simple symmetric case, no NaNs
    # ============================================================
    E = np.zeros((n_reps, n_neighbors, n_points), dtype=np.float64, order="F")
    E[:, 0, 0] = [-2.0, -0.5, 0.2]
    E[:, 1, 0] = [1.7, 0.9, -1.2]
    E[:, :, 1] = 0.0
    E[:, 0, 2] = [2.5, -3.0, 1.2]
    E[:, 1, 2] = [0.4, -0.1, 0.0]

    counts, pmf, included = build_residual_histograms_filtered(E, R, n_bins, n_bins_per_point, neighbor_mask=np.full((n_neighbors, n_points), False, order="F")).values()

    assert np.all(counts == 0), "All counts should be zero"
    assert np.allclose(pmf, 0.0, atol=TOL), "All pmfs should be zero"
    assert np.all(included == 0), "All included should be zero"

    # the mask is per (neighbor, point), not per residual
    def filtered(E, R, n_bins, n_bins_per_point): return build_residual_histograms_filtered(E, R, n_bins, n_bins_per_point, neighbor_mask=np.full(E.shape[1:], True, order="F"))

    for func in (build_residual_histograms, filtered):
        print(f"... test {func.__name__.replace("filtered", "build_residual_histograms_filtered")}")

        E = np.zeros((n_reps, n_neighbors, n_points), dtype=np.float64, order="F")
        E[:, 0, 0] = [-2.0, -0.5, 0.2]
        E[:, 1, 0] = [1.7, 0.9, -1.2]
        E[:, :, 1] = 0.0
        E[:, 0, 2] = [2.5, -3.0, 1.2]
        E[:, 1, 2] = [0.4, -0.1, 0.0]

        out = func(E, R, n_bins, n_bins_per_point)
        counts = out["counts"]
        pmf = out["pmf"]
        included = out["included_n_reps"]

        expected_counts = np.array([
            [2, 1, 2, 1],
            [0, 0, 6, 0],
            [1, 1, 2, 2],
        ], dtype=np.int32, order="F")

        expected_pmf = np.array([
            [2/6, 1/6, 2/6, 1/6],
            [0/6, 0/6, 6/6, 0/6],
            [1/6, 1/6, 2/6, 2/6],
        ], dtype=np.float64, order="F")

        assert np.array_equal(counts, expected_counts), "Test 1: counts mismatch"
        assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 1: pmf mismatch"
        assert np.array_equal(included, [6, 6, 6]), "Test 1: included mismatch"

        # ============================================================
        # Test 2 — NaNs must be ignored
        # ============================================================
        E[:, :, :] = 0
        E[0, 0, 0] = np.nan
        E[2, 0, 0] = np.nan
        E[1, 0, 1] = np.nan
        E[2, 1, 2] = np.nan

        out = func(E, R, n_bins, n_bins_per_point)
        counts = out["counts"]
        pmf = out["pmf"]
        included = out["included_n_reps"]

        expected_counts = np.array([
            [0, 0, 4, 0],
            [0, 0, 5, 0],
            [0, 0, 5, 0],
        ], dtype=np.int32, order="F")

        expected_pmf = np.array([
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ], dtype=np.float64, order="F")

        assert np.array_equal(counts, expected_counts), "Test 2: counts mismatch"
        assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 2: pmf mismatch"
        assert included[0] == 4
        assert included[1] == 5
        assert included[2] == 5

        # ============================================================
        # Test 3 — All NaN → pmf = 0, counts = 0, included = 0
        # ============================================================
        E[:, :, :] = np.nan

        out = func(E, R, n_bins, n_bins_per_point)
        counts = out["counts"]
        pmf = out["pmf"]
        included = out["included_n_reps"]

        assert np.all(counts == 0), "Test 3: counts mismatch"
        assert np.all(pmf == 0), "Test 3: pmf mismatch"
        assert np.all(included == 0), "Test 3: included mismatch"

        # ============================================================
        # Test 4 — Residuals exactly on boundaries
        # ============================================================
        E = np.array([
            [
                [-2.0, -2.0, -2.0],
                [-1.0, -1.0, -1.0],
                [0.0,  0.0,  0.0]
            ],
            [
                [1.0,  1.0,  1.0],
                [2.0,  2.0,  2.0],
                [0.0,  0.0,  0.0]
            ]
        ], dtype=np.float64, order="F")

        out = func(E, R, n_bins, n_bins_per_point)
        counts = out["counts"]
        pmf = out["pmf"]
        included = out["included_n_reps"]

        expected_counts = np.array([
            [1, 1, 2, 2],
            [1, 1, 2, 2],
            [1, 1, 2, 2],
        ], dtype=np.int32, order="F")

        expected_pmf = np.array([
            [1/6, 1/6, 1/3, 1/3],
            [1/6, 1/6, 1/3, 1/3],
            [1/6, 1/6, 1/3, 1/3],
        ], dtype=np.float64, order="F")

        assert np.array_equal(counts, expected_counts), "Test 4: counts mismatch"
        assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 4: pmf mismatch"
        assert np.all(included == 6), "Test 4: included mismatch"


def test_tox_compute_divergence_per_reference_point():
    n_points = 3
    n_bins = 4

    # ============================================================
    # Test 1 — Identical PMFs → JSD = 0
    # ============================================================
    p = np.array([
        [0.1, 0.2, 0.3, 0.4],
        [0.25, 0.25, 0.25, 0.25],
        [1.0, 0.0, 0.0, 0.0],
    ], dtype=np.float64, order="F")

    q = p.copy(order="F")

    jsd = compute_divergence_per_reference_point(p, q)
    expected = np.zeros(n_points, dtype=np.float64)

    assert np.allclose(jsd, expected, atol=TOL), "Test 1 failed: identical PMFs → JSD must be zero"

    # ============================================================
    # Test 2 — Completely disjoint PMFs → JSD = log(2)
    # ============================================================
    p = np.zeros((n_points, n_bins), dtype=np.float64, order="F")
    q = np.zeros((n_points, n_bins), dtype=np.float64, order="F")

    p[0, 0] = 1.0
    q[0, 1] = 1.0

    jsd = compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_points, dtype=np.float64)
    expected[0] = np.log(2.0) / np.log(2.0)  # rescaled onto 0..1 by dividing through LOG_2 = log(2): log(2)/log(2) = 1.0

    assert np.allclose(jsd, expected, atol=TOL), "Test 2 failed: disjoint PMFs → JSD=log(2)/log(2)"
    assert jsd[1] == 0.0, "Test 2 failed: row 2 must be zero"
    assert jsd[2] == 0.0, "Test 2 failed: row 3 must be zero"

    # ============================================================
    # Test 3 — Partially overlapping PMFs (analytic)
    # ============================================================
    p = np.zeros((n_points, n_bins), dtype=np.float64, order="F")
    q = np.zeros((n_points, n_bins), dtype=np.float64, order="F")

    p[0, :] = [0.5, 0.5, 0.0, 0.0]
    q[0, :] = [0.0, 1.0, 0.0, 0.0]

    jsd = compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_points, dtype=np.float64)
    expected[0] = (0.5 * (
        0.5 * np.log(2.0) +
        0.5 * np.log(2.0 / 3.0) +
        np.log(1.0 / 0.75)
    )) / np.log(2.0)  # rescaled onto 0..1 by dividing through LOG_2 = log(2)

    assert np.allclose(jsd, expected, atol=TOL), "Test 3 failed: analytic partial-overlap JSD mismatch"

    # ============================================================
    # Test 4 — Zero-probability bins handled correctly
    # ============================================================
    p = np.zeros((n_points, n_bins), dtype=np.float64, order="F")
    q = np.zeros((n_points, n_bins), dtype=np.float64, order="F")

    p[0, 0] = 1.0
    q[0, 2] = 1.0

    jsd = compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_points, dtype=np.float64)
    expected[0] = np.log(2.0) / np.log(2.0)  # rescaled onto 0..1 by dividing through LOG_2 = log(2): log(2)/log(2) = 1.0

    assert np.allclose(jsd, expected, atol=TOL), "Test 4 failed: zero-probability bins not handled correctly"

    # ============================================================
    # Test 5 — Multiple neighbors, mixed patterns
    # ============================================================
    p = np.array([
        [0.2, 0.0, 0.0, 0.25],
        [0.3, 1.0, 0.0, 0.25],
        [0.5, 0.0, 0.25, 0.25],
    ], dtype=np.float64, order="F")

    q = np.array([
        [0.2, 0.0, 0.0, 0.25],
        [0.3, 0.0, 0.0, 0.25],
        [0.5, 1.0, 0.25, 0.25],
    ], dtype=np.float64, order="F")

    jsd = compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_points, dtype=np.float64)
    expected[1] = (0.5 * (1.0 * np.log(1.0 / 0.5))) / np.log(2.0)  # rescaled onto 0..1: 0.5*log(2)/log(2) = 0.5
    expected[2] = (0.5 * (1.0 * np.log(1.0 / 0.5))) / np.log(2.0)

    assert np.allclose(jsd, expected, atol=TOL), "Test 5 failed: mixed patterns JSD mismatch"


def test_tox_compute_weighted_global_divergence():
    # ============================================================
    # Test 1 — Simple case: equal sample counts → uniform weights
    # ============================================================
    jsd = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64)
    n1  = np.array([5, 5, 5, 5], dtype=np.int32)
    n2  = np.array([5, 5, 5, 5], dtype=np.int32)

    global_jsd, w = compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)

    assert np.allclose(w, expected_weights, atol=TOL), "Test 1 failed: uniform weights mismatch"
    assert abs(global_jsd - 0.25) < TOL, "Test 1 failed: global JSD mismatch"

    # ============================================================
    # Test 2 — Unequal sample counts → weighted average
    # ============================================================
    jsd = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    n1  = np.array([10, 20, 30, 40], dtype=np.int32)
    n2  = np.array([ 0, 10, 10, 10], dtype=np.int32)

    global_jsd, w = compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([10, 30, 40, 50], dtype=np.float64) / 130.0

    assert np.allclose(w, expected_weights, atol=TOL), "Test 2 failed: weights mismatch"

    expected = (
        1.0 * (10/130) +
        2.0 * (30/130) +
        3.0 * (40/130) +
        4.0 * (50/130)
    )

    assert abs(global_jsd - expected) < TOL, "Test 2 failed: weighted global JSD mismatch"

    # ============================================================
    # Test 3 — Some neighborhoods have zero samples → weight = 0
    # ============================================================
    jsd = np.array([0.5, 1.0, 2.0, 4.0], dtype=np.float64)
    n1  = np.array([0, 10, 0, 5], dtype=np.int32)
    n2  = np.array([0,  0, 0, 5], dtype=np.int32)

    global_jsd, w = compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([0.0, 0.5, 0.0, 0.5], dtype=np.float64)

    assert np.allclose(w, expected_weights, atol=TOL), "Test 3 failed: zero-sample weights mismatch"
    assert abs(global_jsd - 2.5) < TOL, "Test 3 failed: global JSD mismatch"

    # ============================================================
    # Test 4 — All neighborhoods have zero samples → weights=0, global JSD=0
    # ============================================================
    jsd = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    n1  = np.array([0, 0, 0, 0], dtype=np.int32)
    n2  = np.array([0, 0, 0, 0], dtype=np.int32)

    global_jsd, w = compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.zeros(4, dtype=np.float64)

    assert np.allclose(w, expected_weights, atol=TOL), "Test 4 failed: all-zero weights mismatch"
    assert abs(global_jsd - 0.0) < TOL, "Test 4 failed: global JSD mismatch"

    # ============================================================
    # Test 5 — Mixed jsd values, mixed sample counts
    # ============================================================
    jsd = np.array([0.0, 0.5, 1.0, 2.0], dtype=np.float64)
    n1  = np.array([5, 0, 10, 5], dtype=np.int32)
    n2  = np.array([5, 5,  0, 5], dtype=np.int32)

    global_jsd, w = compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([10, 5, 10, 10], dtype=np.float64) / 35.0

    assert np.allclose(w, expected_weights, atol=TOL), "Test 5 failed: mixed weights mismatch"

    expected = (
        0.0 * (10/35) +
        0.5 * ( 5/35) +
        1.0 * (10/35) +
        2.0 * (10/35)
    )

    assert abs(global_jsd - expected) < TOL, "Test 5 failed: weighted global JSD mismatch"


def test_calc_pmf():
    """Test calc_pmf: the counts-to-pmf half of build_residual_histograms, factored out."""

    # ============================================================
    # Test 1 — matches build_residual_histograms's own basic fixture (see
    # test_tox_build_residual_histograms Test 1 above): feeding it the exact counts/
    # included_n_reps that routine produces must reproduce that fixture's expected pmf exactly.
    # ============================================================
    counts = np.array([
        [2, 1, 2, 1],
        [0, 0, 6, 0],
        [1, 1, 2, 2],
    ], dtype=np.int32, order="F")
    included_n_reps = np.array([6, 6, 6], dtype=np.int32, order="F")

    expected_pmf = np.array([
        [2/6, 1/6, 2/6, 1/6],
        [0/6, 0/6, 6/6, 0/6],
        [1/6, 1/6, 2/6, 2/6],
    ], dtype=np.float64, order="F")

    pmf = calc_pmf(counts, included_n_reps)
    assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 1 failed: pmf mismatch"

    # ============================================================
    # Test 2 — a reference point with zero included replicates (every residual there was NaN)
    # must report an all-zero pmf row rather than dividing by zero.
    # ============================================================
    counts = np.array([
        [3, 5],
        [0, 0],
    ], dtype=np.int32, order="F")
    included_n_reps = np.array([8, 0], dtype=np.int32, order="F")

    expected_pmf = np.array([
        [0.375, 0.625],
        [0.0, 0.0],
    ], dtype=np.float64, order="F")

    pmf = calc_pmf(counts, included_n_reps)
    assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 2 failed: pmf mismatch"

    # ============================================================
    # Test 3 — counts is documented as non-negative; a negative entry must be rejected.
    # ============================================================
    counts = np.array([[-1]], dtype=np.int32, order="F")
    included_n_reps = np.array([1], dtype=np.int32, order="F")

    assert_error(lambda: calc_pmf(counts, included_n_reps), "Test 3 failed: expected ERR_INVALID_INPUT", ERR_INVALID_INPUT)


def test_determine_all_studies_shared_residual_range():

    # ============================================================
    # Test 1 — Three single-replicate studies, hand-computed: pooled absolute residuals sorted
    # are [3, 4, 5]; the default 95% quantile has rank 0.95*(3-1)+1 = 2.9, so
    # R = sorted(2) + 0.9*(sorted(3)-sorted(2)) = 4 + 0.9*1 = 4.9
    # ============================================================
    all_studies = np.full((1, 1, 1, 3), np.nan, dtype=np.float64, order="F")
    all_studies[0, 0, 0, 0] = 3.0
    all_studies[0, 0, 0, 1] = -4.0
    all_studies[0, 0, 0, 2] = 5.0

    R = determine_all_studies_shared_residual_range(all_studies)
    assert abs(R - 4.9) < TOL, f"Test 1 failed: expected 4.9, got {R}"

    # ============================================================
    # Test 2 — Regression-safety cross-check: feeding the same two studies both through
    # determine_study_shared_residual_range (with S2 as-is) and through the N-study routine (S2
    # padded with NaN up to S1's replicate count, n_studies=2) must give the exact same range.
    # ============================================================
    S1 = np.zeros((4, 2, 2), dtype=np.float64, order="F")
    S1[:, 0, 0] = [1,  2,  3, 4]
    S1[:, 1, 0] = [5,  6, -7, 8]
    S1[:, 0, 1] = [9, 10, 11, 12]
    S1[:, 1, 1] = [1, 1, 1, 1]

    S2 = np.zeros((3, 2, 2), dtype=np.float64, order="F")
    S2[:, 0, 0] = [2, -4,  6]
    S2[:, 1, 0] = [8,  1,  3]
    S2[:, 0, 1] = [5,  7,  9]
    S2[:, 1, 1] = [0,  1,  2]

    R_two_study = determine_study_shared_residual_range(S1, S2, 0.95)

    all_studies = np.full((4, 2, 2, 2), np.nan, dtype=np.float64, order="F")
    all_studies[:, :, :, 0] = S1
    all_studies[:3, :, :, 1] = S2  # all_studies[3, :, :, 1] stays NaN -- padding S2 up to max_n_reps

    R_all_studies = determine_all_studies_shared_residual_range(all_studies)
    assert np.isclose(R_all_studies, R_two_study, atol=TOL), \
        f"Test 2 failed: N-study range {R_all_studies} != two-study range {R_two_study}"
    assert abs(R_all_studies - 10.65) < TOL, f"Test 2 failed: expected 10.65, got {R_all_studies}"


def main():
    run_all_tests(globals().values())


if __name__ == "__main__":
    main()
