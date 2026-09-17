"""
Comprehensive Python test suite for the gJCT/fJCT JSD-Comp-Test parameter-search building
blocks in tensoromics (Issue #126).
Uses tensor_omics.py wrapper functions (mirrors Fortran test suite
mod_test_data_integration_js_comp_test.F90 -- fixtures below are translated from that suite's
already hand-computed/verified numbers rather than re-derived).
"""
import numpy as np
import sys
import os

# Add parent directory to path to import tensor_omics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from test_helpers import run_all_tests, assert_error
from tensor_omics import (
    estimate_bin_count,
    estimate_bin_count_expert,
    determine_bin_count_occupancy,
    determine_bin_count_occupancy_expert,
    generate_js_comp_test_candidates,
    generate_js_comp_test_candidates_expert,
    check_neighborhood_overlaps,
    check_mean_pmf_min_counts,
    check_plateau_condition,
    create_mean_pmf,
    create_mean_pmf_only,
    bootstrap_histogram,
    calc_js_comp_test_candidate_bounds,
    calc_js_comp_test_n_top_k_jsds,
    run_js_comp_test,
    run_js_comp_test_parameter_search,
)
from tensor_omics.error_handling import ERR_INVALID_INPUT

TOL = 1e-12


def test_estimate_bin_count():
    # NOTE: estimate_bin_count/estimate_bin_count_expert now return a dict with keys n_bins,
    # sturges_bins and fd_bins (Issue #187 checkpoint 1). Per this project's testing philosophy
    # (Fortran_Coding_Guides.pdf Sec 17.1), this suite only checks call-ability, return
    # type/shape, and the n_bins == max(sturges_bins, fd_bins) structural invariant --
    # individual numerical correctness of sturges_bins/fd_bins is the Fortran suite's job
    # (mod_test_data_integration_js_comp_test.F90).

    # ============================================================
    # Test 1 -- basic hand-computed (test_estimate_bin_count_basic_hand_computed):
    # residuals already ascending, no ties needed. Sturges gives 1 + nint(log(1)/LOG_2) = 1.
    # 25th/75th percentiles (rank 3.25 and 7.75) interpolate to -1.75 and 2.75, so the
    # Freedman-Diaconis half-width is (2.75 - (-1.75)) / 1^(1/3) = 4.5; with
    # shared_residual_range=9.0, nint(9.0/4.5) = 2, which beats Sturges.
    # ============================================================
    residuals = np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)

    result = estimate_bin_count(residuals, max_n_reps_all_studies=1, n_neighbors=1, shared_residual_range=9.0)
    assert isinstance(result, dict), f"Test 1 failed: expected a dict, got {type(result)}"
    for key in ("n_bins", "sturges_bins", "fd_bins"):
        assert key in result, f"Test 1 failed: missing expected output key '{key}'"
    assert result["n_bins"] == 2, f"Test 1 failed: expected n_bins=2, got {result['n_bins']}"
    assert result["n_bins"] == max(result["sturges_bins"], result["fd_bins"]), \
        "Test 1 failed: n_bins must equal max(sturges_bins, fd_bins)"

    # The expert entry point, given the same sorting permutation, must agree.
    residuals_perm = (np.argsort(residuals, kind="mergesort") + 1).astype(np.int32)
    result_expert = estimate_bin_count_expert(residuals, residuals_perm, max_n_reps_all_studies=1, n_neighbors=1,
                                               shared_residual_range=9.0)
    assert result_expert["n_bins"] == 2, f"Test 1 (expert) failed: expected n_bins=2, got {result_expert['n_bins']}"
    assert result_expert["n_bins"] == max(result_expert["sturges_bins"], result_expert["fd_bins"]), \
        "Test 1 (expert) failed: n_bins must equal max(sturges_bins, fd_bins)"

    # ============================================================
    # Test 2 -- all-NaN pool falls back to a single bin (test_estimate_bin_count_all_nan_gives_one_bin)
    # ============================================================
    residuals_nan = np.full(4, np.nan, dtype=np.float64)
    result_nan = estimate_bin_count(residuals_nan, max_n_reps_all_studies=1, n_neighbors=1, shared_residual_range=9.0)
    assert result_nan["n_bins"] == 1, f"Test 2 failed: expected n_bins=1, got {result_nan['n_bins']}"
    assert result_nan["sturges_bins"] == 1, f"Test 2 failed: expected sturges_bins=1, got {result_nan['sturges_bins']}"
    assert result_nan["fd_bins"] == 1, f"Test 2 failed: expected fd_bins=1, got {result_nan['fd_bins']}"

    # ============================================================
    # Test 3 -- clamped to MAX_N_BINS=256 (test_estimate_bin_count_clamped_to_max_n_bins):
    # same residuals, but shared_residual_range=5000.0 -> raw estimate nint(5000.0/4.5)=1111,
    # far above MAX_N_BINS.
    # ============================================================
    result_clamped = estimate_bin_count(residuals, max_n_reps_all_studies=1, n_neighbors=1,
                                         shared_residual_range=5000.0)
    assert result_clamped["n_bins"] == 256, f"Test 3 failed: expected n_bins=256, got {result_clamped['n_bins']}"
    assert result_clamped["n_bins"] == max(result_clamped["sturges_bins"], result_clamped["fd_bins"]), \
        "Test 3 failed: n_bins must equal max(sturges_bins, fd_bins)"


def test_determine_bin_count_occupancy():
    # Call-ability and return type/shape only, per this project's testing philosophy
    # (Fortran_Coding_Guides.pdf Sec 17.1) -- numerical branch coverage of Issue #187's
    # geometric-search-then-refinement algorithm is the Fortran suite's job
    # (mod_test_data_integration_js_comp_test.F90's own 9 enumerated
    # test_determine_bin_count_occupancy_* tests).
    residuals = np.arange(-60, 60, dtype=np.float64)  # 120 values, matches a Fortran fixture

    result = determine_bin_count_occupancy(residuals, max_n_reps_all_studies=1, n_neighbors=1,
                                            shared_residual_range=60.0)
    assert isinstance(result, dict), f"expected a dict, got {type(result)}"
    for key in ("selected_n_bins", "occupancy_failed", "n_pooled_residuals", "min_bin_occupancy",
                "mean_bin_occupancy", "max_bin_occupancy", "sturges_bins", "fd_bins"):
        assert key in result, f"missing expected output key '{key}'"
    assert isinstance(result["occupancy_failed"], (bool, np.bool_)), \
        f"expected occupancy_failed to be a bool, got {type(result['occupancy_failed'])}"
    assert result["n_pooled_residuals"] == 120, \
        f"expected n_pooled_residuals=120, got {result['n_pooled_residuals']}"
    assert result["selected_n_bins"] == 12, f"expected selected_n_bins=12, got {result['selected_n_bins']}"
    assert not result["occupancy_failed"]

    # The expert entry point, given the same sorting permutation, must agree.
    residuals_perm = (np.argsort(residuals, kind="mergesort") + 1).astype(np.int32)
    result_expert = determine_bin_count_occupancy_expert(residuals, residuals_perm, max_n_reps_all_studies=1,
                                                           n_neighbors=1, shared_residual_range=60.0)
    assert result_expert["selected_n_bins"] == result["selected_n_bins"], \
        "expert entry point should agree with the plain one given the same sorted permutation"

    # All-NaN pool: occupancy_failed must be True, and every occupancy diagnostic must be 0.
    residuals_nan = np.full(4, np.nan, dtype=np.float64)
    result_nan = determine_bin_count_occupancy(residuals_nan, max_n_reps_all_studies=1, n_neighbors=1,
                                                 shared_residual_range=9.0)
    assert result_nan["occupancy_failed"]
    assert result_nan["n_pooled_residuals"] == 0
    assert result_nan["min_bin_occupancy"] == 0
    assert result_nan["max_bin_occupancy"] == 0


def test_generate_js_comp_test_candidates():
    residuals = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)

    # ============================================================
    # Test 1 -- small-N candidate-grid collapse at its exact threshold
    # (test_generate_js_comp_test_candidates_collapses_at_8742): with
    # max_n_genes_all_studies=8742, n_points_high = clamp(ceil(4*sqrt(8742)), 300, 1500) = 374
    # and n_points_low = max(300, ceil(0.2*374)) = 300; after one grid iteration n_points_high
    # becomes 374*0.8=299.2 < 300, so every candidate the grid returns shares n_points=374.
    # ============================================================
    out = generate_js_comp_test_candidates(8742, residuals, max_n_reps_all_studies=1, shared_residual_range=1.0)
    candidates = out["candidates_n_points_n_neighbors"]
    n_candidates = candidates.shape[1]
    assert n_candidates >= 1, "Test 1 failed: expected at least one candidate"
    assert np.all(candidates[0, :] == candidates[0, 0]), \
        "Test 1 failed: every candidate must share the same n_points"
    assert candidates[0, 0] == 374, f"Test 1 failed: n_points should be 374, got {candidates[0, 0]}"

    # ============================================================
    # Test 2 -- the other side of the bracket (test_generate_js_comp_test_candidates_has_two_distinct_at_8743):
    # with max_n_genes_all_studies=8743, n_points_high=ceil(4*sqrt(8743))=375 and
    # n_points_low=300; after one iteration n_points_high becomes 375*0.8=300.0, which is NOT
    # less than n_points_low=300, so a second, distinct n_points=300 candidate appears.
    # ============================================================
    out2 = generate_js_comp_test_candidates(8743, residuals, max_n_reps_all_studies=1, shared_residual_range=1.0)
    candidates2 = out2["candidates_n_points_n_neighbors"]
    n_candidates2 = candidates2.shape[1]
    n_distinct = len(set(candidates2[0, :].tolist()))
    assert n_distinct >= 2, "Test 2 failed: expected at least two distinct n_points values"
    assert candidates2[0, 0] == 375, f"Test 2 failed: first n_points should be 375, got {candidates2[0, 0]}"
    assert candidates2[0, n_candidates2 - 1] == 300, \
        f"Test 2 failed: last n_points should collapse to 300, got {candidates2[0, n_candidates2 - 1]}"

    # ============================================================
    # Test 3 -- validation (test_generate_js_comp_test_candidates_validation):
    # max_n_genes_all_studies must be positive; the expert tier accepts a valid permutation.
    # ============================================================
    assert_error(lambda: generate_js_comp_test_candidates(0, residuals, max_n_reps_all_studies=1,
                                                           shared_residual_range=1.0),
                 "Test 3 failed: expected ERR_INVALID_INPUT for max_n_genes_all_studies=0", ERR_INVALID_INPUT)

    residuals_perm = np.array([1, 2, 3, 4, 5], dtype=np.int32)
    out3 = generate_js_comp_test_candidates_expert(100, residuals, residuals_perm, max_n_reps_all_studies=1,
                                                    shared_residual_range=1.0)
    assert out3["candidates_n_points_n_neighbors"].shape[1] >= 1, \
        "Test 3 failed: expert tier should accept a valid permutation"


def test_check_neighborhood_overlaps():
    # ============================================================
    # Test 1 -- identical ranges always pass, even at the strictest threshold
    # (test_neighborhood_overlaps_all_pass)
    # ============================================================
    neighborhood_range = np.array([[1, 1, 1], [10, 10, 10]], dtype=np.int32, order='F')
    all_pass = check_neighborhood_overlaps(neighborhood_range, min_neighbor_overlap=1.0)
    assert all_pass, "Test 1 failed: identical ranges should always pass"

    # ============================================================
    # Test 2 -- overlap exactly at threshold passes (test_neighborhood_overlaps_exactly_at_threshold_passes):
    # ranges [1,10] and [6,10]: left_max=min(10,10)=10, right_min=max(1,6)=6, overlap=(10-6)/(10-1)=4/9.
    # The gate is overlap >= min_neighbor_overlap, so setting the threshold to that exact
    # fraction must still pass.
    # ============================================================
    neighborhood_range_2 = np.array([[1, 6], [10, 10]], dtype=np.int32, order='F')
    all_pass_2 = check_neighborhood_overlaps(neighborhood_range_2, min_neighbor_overlap=4.0 / 9.0)
    assert all_pass_2, "Test 2 failed: overlap exactly at threshold must pass"

    # ============================================================
    # Test 3 -- same 4/9 overlap, threshold 0.5 > 4/9, must fail
    # (test_neighborhood_overlaps_below_threshold_fails)
    # ============================================================
    all_fail = check_neighborhood_overlaps(neighborhood_range_2, min_neighbor_overlap=0.5)
    assert not all_fail, "Test 3 failed: overlap below threshold must fail"


def test_check_mean_pmf_min_counts():
    # ============================================================
    # Test 1 -- every bin comfortably above the minimum (test_mean_pmf_min_counts_all_pass)
    # n_bins_per_point uniformly equals n_bins (call-ability/shape only, per this file's convention)
    # ============================================================
    mean_pmf_counts = np.array([[10, 10], [10, 10]], dtype=np.int32, order='F')
    n_bins_per_point = np.array([2, 2], dtype=np.int32)
    assert check_mean_pmf_min_counts(mean_pmf_counts, n_bins_per_point, min_count=5), \
        "Test 1 failed: all bins above minimum should pass"

    # ============================================================
    # Test 2 -- every bin equals the minimum exactly: gate is >=, must pass
    # (test_mean_pmf_min_counts_exactly_at_threshold_passes)
    # ============================================================
    mean_pmf_counts_2 = np.array([[3, 3], [3, 3]], dtype=np.int32, order='F')
    n_bins_per_point_2 = np.array([2, 2], dtype=np.int32)
    assert check_mean_pmf_min_counts(mean_pmf_counts_2, n_bins_per_point_2, min_count=3), \
        "Test 2 failed: count exactly at minimum must pass"

    # ============================================================
    # Test 3 -- one bin one below the minimum: gate must fail
    # (test_mean_pmf_min_counts_below_threshold_fails)
    # ============================================================
    mean_pmf_counts_3 = np.array([[3, 3], [3, 2]], dtype=np.int32, order='F')
    n_bins_per_point_3 = np.array([2, 2], dtype=np.int32)
    assert not check_mean_pmf_min_counts(mean_pmf_counts_3, n_bins_per_point_3, min_count=3), \
        "Test 3 failed: one bin below minimum must fail the whole gate"


def test_check_plateau_condition():
    # ============================================================
    # Test 1 -- METHOD_JOIN_MIN requires every study to pass (test_check_plateau_condition_join_min_requires_all).
    # 2-of-3 studies passing must fail; 3-of-3 must succeed.
    # ============================================================
    confidence_interval = np.array([[0.4, 0.4, -2.0], [0.6, 0.6, -1.0]], dtype=np.float64, order='F')
    best_ci = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64, order='F')

    result = check_plateau_condition(confidence_interval, best_ci, best_candidate_index=0,
                                      best_exceeded_ci_overlap_count=0, candidate_index=1,
                                      join_method='join_min', succeeding_ci_overlap=0.9)
    assert not result["plateau_found"], "Test 1 failed: 2-of-3 studies must not plateau"

    confidence_interval[:, 2] = [0.4, 0.6]
    best_ci_2 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64, order='F')
    result2 = check_plateau_condition(confidence_interval, best_ci_2, best_candidate_index=0,
                                       best_exceeded_ci_overlap_count=0, candidate_index=1,
                                       join_method='join_min', succeeding_ci_overlap=0.9)
    assert result2["plateau_found"], "Test 1 failed: 3-of-3 studies must plateau"

    # ============================================================
    # Test 2 -- METHOD_JOIN_MAX succeeds once any one study passes
    # (test_check_plateau_condition_join_max_requires_any). 0-of-3 must fail, 1-of-3 must succeed.
    # ============================================================
    confidence_interval_max = np.array([[-2.0, -2.0, -2.0], [-1.0, -1.0, -1.0]], dtype=np.float64, order='F')
    best_ci_max = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64, order='F')
    result_max_0 = check_plateau_condition(confidence_interval_max, best_ci_max, best_candidate_index=0,
                                            best_exceeded_ci_overlap_count=0, candidate_index=1,
                                            join_method='join_max', succeeding_ci_overlap=0.9)
    assert not result_max_0["plateau_found"], "Test 2 failed: 0-of-3 studies must not plateau"

    confidence_interval_max[:, 0] = [0.4, 0.6]
    best_ci_max_2 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64, order='F')
    result_max_1 = check_plateau_condition(confidence_interval_max, best_ci_max_2, best_candidate_index=0,
                                            best_exceeded_ci_overlap_count=0, candidate_index=1,
                                            join_method='join_max', succeeding_ci_overlap=0.9)
    assert result_max_1["plateau_found"], "Test 2 failed: 1-of-3 studies must plateau"

    # ============================================================
    # Test 3 -- METHOD_JOIN_MEDIAN succeeds once a majority pass (n_studies=3 -> threshold 2)
    # (test_check_plateau_condition_join_median_requires_majority). 1-of-3 must fail, 2-of-3 must succeed.
    # ============================================================
    confidence_interval_med = np.array([[0.4, -2.0, -2.0], [0.6, -1.0, -1.0]], dtype=np.float64, order='F')
    best_ci_med = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64, order='F')
    result_med_1 = check_plateau_condition(confidence_interval_med, best_ci_med, best_candidate_index=0,
                                            best_exceeded_ci_overlap_count=0, candidate_index=1,
                                            join_method='join_median', succeeding_ci_overlap=0.9)
    assert not result_med_1["plateau_found"], "Test 3 failed: 1-of-3 must not plateau"

    confidence_interval_med[:, 1] = [0.4, 0.6]
    best_ci_med_2 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64, order='F')
    result_med_2 = check_plateau_condition(confidence_interval_med, best_ci_med_2, best_candidate_index=0,
                                            best_exceeded_ci_overlap_count=0, candidate_index=1,
                                            join_method='join_median', succeeding_ci_overlap=0.9)
    assert result_med_2["plateau_found"], "Test 3 failed: 2-of-3 must plateau"

    # ============================================================
    # Test 4 -- a candidate worse than the previous best plateaus immediately, regardless of
    # join_method, and none of the best_* state is overwritten
    # (test_check_plateau_condition_worse_than_previous_short_circuits).
    # ============================================================
    confidence_interval_worse = np.array([[-5.0, -5.0], [-4.0, -4.0]], dtype=np.float64, order='F')
    best_ci_worse = np.array([[0.1, 0.3], [0.2, 0.4]], dtype=np.float64, order='F')
    expected_best_ci_worse = best_ci_worse.copy()

    result_worse = check_plateau_condition(confidence_interval_worse, best_ci_worse, best_candidate_index=5,
                                            best_exceeded_ci_overlap_count=2, candidate_index=9,
                                            join_method='join_min', succeeding_ci_overlap=0.5)
    assert result_worse["plateau_found"], "Test 4 failed: a worse candidate must plateau immediately"
    assert result_worse["best_candidate_index"] == 5, "Test 4 failed: best_candidate_index must not be overwritten"
    assert result_worse["best_exceeded_ci_overlap_count"] == 2, \
        "Test 4 failed: best_exceeded_ci_overlap_count must not be overwritten"
    np.testing.assert_array_almost_equal(best_ci_worse, expected_best_ci_worse, decimal=12)


def test_create_mean_pmf():
    # ============================================================
    # Test 1 -- two studies, hand-computed (test_create_mean_pmf_two_studies_hand_computed):
    # S1 pmf=[0.3,0.7] counts=[3,7] included=10; S2 pmf=[0.5,0.5] counts=[5,5] included=10.
    # mean_pmf = ([0.3,0.7]+[0.5,0.5])/2 = [0.4,0.6], mean_pmf_counts=[8,12], included=20.
    # ============================================================
    pmfs = np.zeros((2, 1, 2), dtype=np.float64, order='F')
    pmfs[:, 0, 0] = [0.3, 0.7]
    pmfs[:, 0, 1] = [0.5, 0.5]
    counts = np.zeros((2, 1, 2), dtype=np.int32, order='F')
    counts[:, 0, 0] = [3, 7]
    counts[:, 0, 1] = [5, 5]
    included_n_reps = np.array([[10, 10]], dtype=np.int32, order='F')

    result = create_mean_pmf(pmfs, counts, included_n_reps)
    np.testing.assert_array_almost_equal(result["mean_pmf"][:, 0], [0.4, 0.6], decimal=12)
    assert list(result["mean_pmf_counts"][:, 0]) == [8, 12], "Test 1 failed: mean_pmf_counts mismatch"
    assert result["mean_pmf_included_n_reps"][0] == 20, "Test 1 failed: mean_pmf_included_n_reps mismatch"

    # ============================================================
    # Test 2 -- three studies, DOCUMENTS THE KNOWN LIMITATION: averages over ALL studies
    # including self, not a true leave-one-out background (test_create_mean_pmf_three_studies_includes_self).
    # S1 pmf=[0.2,0.8], S2 pmf=[0.4,0.6], S3 pmf=[0.9,0.1] -> mean_pmf=[0.5,0.5] over all three.
    # ============================================================
    pmfs3 = np.zeros((2, 1, 3), dtype=np.float64, order='F')
    pmfs3[:, 0, 0] = [0.2, 0.8]
    pmfs3[:, 0, 1] = [0.4, 0.6]
    pmfs3[:, 0, 2] = [0.9, 0.1]
    counts3 = np.zeros((2, 1, 3), dtype=np.int32, order='F')
    counts3[:, 0, 0] = [2, 8]
    counts3[:, 0, 1] = [4, 6]
    counts3[:, 0, 2] = [9, 1]
    included_n_reps3 = np.array([[10, 10, 10]], dtype=np.int32, order='F')

    result3 = create_mean_pmf(pmfs3, counts3, included_n_reps3)
    np.testing.assert_array_almost_equal(result3["mean_pmf"][:, 0], [0.5, 0.5], decimal=12)
    assert list(result3["mean_pmf_counts"][:, 0]) == [15, 15], "Test 2 failed: mean_pmf_counts mismatch"
    assert result3["mean_pmf_included_n_reps"][0] == 30, "Test 2 failed: mean_pmf_included_n_reps mismatch"

    # What a true leave-one-out background for study 1 would be, for contrast -- NOT what
    # mean_pmf above equals, and not asserted against it.
    leave_one_out_study1 = (np.array([0.4, 0.6]) + np.array([0.9, 0.1])) / 2.0
    assert abs(result3["mean_pmf"][0, 0] - leave_one_out_study1[0]) > TOL, \
        "Test 2 failed: mean_pmf must differ from the leave-one-out background -- known limitation was fixed " \
        "without updating this test"

    # ============================================================
    # Test 3 -- validation: pmfs must be in [0,1], counts must be non-negative
    # (test_create_mean_pmf_validation)
    # ============================================================
    pmfs_bad = np.array([[[1.5]]], dtype=np.float64, order='F')
    counts_bad = np.array([[[1]]], dtype=np.int32, order='F')
    included_n_reps_bad = np.array([[1]], dtype=np.int32, order='F')
    assert_error(lambda: create_mean_pmf(pmfs_bad, counts_bad, included_n_reps_bad),
                 "Test 3 failed: expected ERR_INVALID_INPUT for pmf above 1.0", ERR_INVALID_INPUT)

    pmfs_ok = np.array([[[0.5]]], dtype=np.float64, order='F')
    counts_neg = np.array([[[-1]]], dtype=np.int32, order='F')
    assert_error(lambda: create_mean_pmf(pmfs_ok, counts_neg, included_n_reps_bad),
                 "Test 3 failed: expected ERR_INVALID_INPUT for negative count", ERR_INVALID_INPUT)


def test_create_mean_pmf_only():
    # `create_mean_pmf_only` must produce exactly the same mean_pmf as `create_mean_pmf`
    # (test_create_mean_pmf_only_matches_create_mean_pmf).
    pmfs = np.zeros((2, 1, 2), dtype=np.float64, order='F')
    pmfs[:, 0, 0] = [0.3, 0.7]
    pmfs[:, 0, 1] = [0.5, 0.5]
    counts = np.zeros((2, 1, 2), dtype=np.int32, order='F')
    counts[:, 0, 0] = [3, 7]
    counts[:, 0, 1] = [5, 5]
    included_n_reps = np.array([[10, 10]], dtype=np.int32, order='F')

    full = create_mean_pmf(pmfs, counts, included_n_reps)
    only = create_mean_pmf_only(pmfs)

    np.testing.assert_array_almost_equal(only, full["mean_pmf"], decimal=12)


def test_bootstrap_histogram():
    # ============================================================
    # Test 1 -- same random_seed twice must give bit-for-bit identical confidence_interval
    # (test_bootstrap_histogram_seeded_reproducibility).
    # ============================================================
    mean_pmf_counts = np.array([[5, 2], [3, 2], [2, 6]], dtype=np.int32, order='F')
    mean_pmf_included_n_reps = np.array([10, 10], dtype=np.int32)
    included_n_reps = np.array([[4, 6], [4, 6]], dtype=np.int32, order='F')

    ci_first = np.array([[0.2, 0.3], [0.2, 0.3]], dtype=np.float64, order='F')
    bootstrap_histogram(15, mean_pmf_counts, mean_pmf_included_n_reps, included_n_reps, ci_first, random_seed=123)

    ci_second = np.array([[0.2, 0.3], [0.2, 0.3]], dtype=np.float64, order='F')
    bootstrap_histogram(15, mean_pmf_counts, mean_pmf_included_n_reps, included_n_reps, ci_second, random_seed=123)

    np.testing.assert_array_almost_equal(ci_first, ci_second, decimal=12)

    # ============================================================
    # Test 2 -- degenerate single-bin histogram: JSD is always exactly 0.0, so seeding the
    # confidence interval at [0.0, 0.0] must leave it collapsed at [0.0, 0.0]
    # (test_bootstrap_histogram_degenerate_single_bin_zero_ci).
    # ============================================================
    mean_pmf_counts_deg = np.array([[10]], dtype=np.int32, order='F')
    mean_pmf_included_n_reps_deg = np.array([10], dtype=np.int32)
    included_n_reps_deg = np.array([[5, 5]], dtype=np.int32, order='F')
    confidence_interval_deg = np.zeros((2, 2), dtype=np.float64, order='F')

    bootstrap_histogram(5, mean_pmf_counts_deg, mean_pmf_included_n_reps_deg, included_n_reps_deg,
                         confidence_interval_deg, random_seed=7)

    np.testing.assert_array_almost_equal(confidence_interval_deg, np.zeros((2, 2)), decimal=12)


def test_calc_js_comp_test_candidate_bounds():
    # Cross-check against the already-verified candidate-grid collapse fixture: with
    # max_n_genes_all_studies=8742 the grid's first (largest) n_points candidate is 374
    # (see test_generate_js_comp_test_candidates), and `max_n_points_candidate` is documented
    # as being exactly that value.
    bounds = calc_js_comp_test_candidate_bounds(8742)
    assert bounds["max_n_points_candidate"] == 374, \
        f"expected max_n_points_candidate=374, got {bounds['max_n_points_candidate']}"

    out = generate_js_comp_test_candidates(8742, np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64),
                                            max_n_reps_all_studies=1, shared_residual_range=1.0)
    candidates = out["candidates_n_points_n_neighbors"]
    assert bounds["max_n_neighbors_candidate"] >= int(np.max(candidates[1, :])), \
        "max_n_neighbors_candidate must be a safe upper bound on the grid's n_neighbors"


def test_calc_js_comp_test_n_top_k_jsds():
    # ============================================================
    # Test 1 -- default 2.5% significance level (test_calc_js_comp_test_n_top_k_jsds_default_hand_computed):
    # n_bootstraps=1000 -> n_top_k = max(1, floor(0.025*1000)) = 25.
    # ============================================================
    assert calc_js_comp_test_n_top_k_jsds(1000) == 25, "Test 1 failed: expected n_top_k=25"

    # ============================================================
    # Test 2 -- explicit significance level (test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level):
    # n_bootstraps=200, significance=5.0 -> n_top_k = max(1, floor(0.05*200)) = 10.
    # ============================================================
    assert calc_js_comp_test_n_top_k_jsds(200, two_sided_bootstrapping_significance_level=5.0) == 10, \
        "Test 2 failed: expected n_top_k=10"

    # ============================================================
    # Test 3 -- clamped to at least one (test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one):
    # n_bootstraps=10 -> floor(0.025*10)=0, must clamp up to 1.
    # ============================================================
    assert calc_js_comp_test_n_top_k_jsds(10) == 1, "Test 3 failed: expected n_top_k=1"


def test_run_js_comp_test_two_studies_hand_traceable():
    # End-to-end, fully closed-form 2-study case (test_run_js_comp_test_two_studies_hand_traceable),
    # with n_permutations=0 so the permutation loop never executes and the whole pipeline is
    # deterministic and hand-traceable. See the Fortran fixture's own long-form derivation of
    # EXPECTED_JSD = 1.5 - 0.75*log2(3).
    log2_3 = np.log(3.0) / np.log(2.0)
    expected_jsd = 1.5 - 0.75 * log2_3

    gene_means = np.array([[1.0, 1.0], [5.0, 5.0]], dtype=np.float64, order='F')
    gene_means_perms = np.array([[1, 1], [2, 2]], dtype=np.int32, order='F')

    residuals = np.full((2, 2, 2), np.nan, dtype=np.float64, order='F')
    residuals[:, 0, 0] = [-1.0, 1.0]
    residuals[:, 0, 1] = [-3.0, 3.0]

    x_star = np.array([1.0], dtype=np.float64)

    result = run_js_comp_test(n_neighbors=1, n_bins=4, shared_residual_range=4.0, gene_means=gene_means,
                               gene_means_perms=gene_means_perms, residuals=residuals, x_star=x_star,
                               n_permutations=0, random_seed=1)

    assert list(result["neighborhood_indices"][:, 0, 0]) == [1], "study 1 neighbor should be gene 1"
    assert list(result["neighborhood_indices"][:, 0, 1]) == [1], "study 2 neighbor should be gene 1"

    assert list(result["counts"][:, 0, 0]) == [0, 1, 1, 0], "study 1 counts mismatch"
    assert list(result["counts"][:, 0, 1]) == [1, 0, 0, 1], "study 2 counts mismatch"
    assert result["included_n_reps"][0, 0] == 2
    assert result["included_n_reps"][0, 1] == 2

    np.testing.assert_array_almost_equal(result["mean_pmf"][:, 0], [0.25, 0.25, 0.25, 0.25], decimal=12)
    assert list(result["mean_pmf_counts"][:, 0]) == [1, 1, 1, 1]
    assert result["mean_pmf_included_n_reps"][0] == 4

    assert abs(result["global_js_divergence"][0] - expected_jsd) < 1e-9, \
        f"study 1 global JSD mismatch: expected {expected_jsd}, got {result['global_js_divergence'][0]}"
    assert abs(result["global_js_divergence"][1] - expected_jsd) < 1e-9, \
        f"study 2 global JSD mismatch: expected {expected_jsd}, got {result['global_js_divergence'][1]}"

    assert abs(result["weights"][0, 0] - 1.0) < TOL
    assert abs(result["weights"][0, 1] - 1.0) < TOL

    assert abs(result["p_values"][0] - 0.0) < TOL, "n_permutations=0 -> p_values stay 0.0"
    assert abs(result["p_values"][1] - 0.0) < TOL, "n_permutations=0 -> p_values stay 0.0"


def test_run_js_comp_test_parameter_search():
    max_n_genes_all_studies = 2000
    max_n_reps_all_studies = 3
    n_studies = 2

    gene_means = np.empty((max_n_genes_all_studies, n_studies), dtype=np.float64, order='F')
    residuals = np.empty((max_n_reps_all_studies, max_n_genes_all_studies, n_studies), dtype=np.float64, order='F')
    for i_study in range(n_studies):
        gene_means[:, i_study] = np.arange(1, max_n_genes_all_studies + 1, dtype=np.float64)
        residuals[:, :, i_study] = np.array([-0.5, 0.0, 0.5]).reshape(3, 1)

    # ============================================================
    # Test 1 -- no candidate ever plateaus (min_residuals_per_bin impossibly high), so the
    # search falls back to the FIRST (finest-resolution) candidate and resets
    # best_candidate_pair_confidence_interval to -1.0 throughout
    # (test_param_search_no_plateau_falls_back_to_finest). Both n_points=300 and n_neighbors=26
    # are derived purely from the GAMMA-decay constants and max_n_genes_all_studies=2000.
    # ============================================================
    result = run_js_comp_test_parameter_search(gene_means, residuals, shared_residual_range=1.0, n_bootstraps=10,
                                                join_method='join_min', min_residuals_per_bin=1000000,
                                                random_seed=1)

    assert result["n_points"] == 300, f"expected n_points=300, got {result['n_points']}"
    assert result["n_neighbors"] == 26, f"expected n_neighbors=26, got {result['n_neighbors']}"
    assert result["plateau_established"] == False, "no candidate ever plateaus here"
    np.testing.assert_array_almost_equal(result["best_candidate_pair_confidence_interval"][:, 0], [-1.0, -1.0],
                                          decimal=12)
    np.testing.assert_array_almost_equal(result["best_candidate_pair_confidence_interval"][:, 1], [-1.0, -1.0],
                                          decimal=12)

    # ============================================================
    # Test 2 -- the GAMMA-decay grid collapses to exactly ONE candidate at
    # max_n_genes_all_studies=100 ((n_points, n_neighbors)=(300, 1)); with both admissibility
    # gates relaxed to their most permissive settings the sole candidate is returned regardless
    # of plateau, and the bootstrap actually runs, giving a real (not -1.0) confidence interval
    # (test_param_search_single_candidate_bypasses_plateau).
    # ============================================================
    max_n_genes_2 = 100
    gene_means_2 = np.empty((max_n_genes_2, n_studies), dtype=np.float64, order='F')
    residuals_2 = np.empty((max_n_reps_all_studies, max_n_genes_2, n_studies), dtype=np.float64, order='F')
    for i_study in range(n_studies):
        gene_means_2[:, i_study] = np.arange(1, max_n_genes_2 + 1, dtype=np.float64)
        residuals_2[:, :, i_study] = np.array([-0.5, 0.0, 0.5]).reshape(3, 1)

    result2 = run_js_comp_test_parameter_search(gene_means_2, residuals_2, shared_residual_range=1.0,
                                                 n_bootstraps=5, join_method='join_min',
                                                 min_residuals_per_bin=0, min_neighbor_overlap=0.0,
                                                 random_seed=1)

    assert result2["n_points"] == 300, f"expected n_points=300, got {result2['n_points']}"
    assert result2["n_neighbors"] == 1, f"expected n_neighbors=1, got {result2['n_neighbors']}"
    assert result2["plateau_established"] == True, "the single-candidate bypass counts as established"
    ci = result2["best_candidate_pair_confidence_interval"]
    assert ci[0, 0] >= 0.0 and ci[1, 0] <= 1.0, \
        "relaxed gates should let bootstrap actually run, giving a real (not -1.0) CI"

    # ============================================================
    # Test 3 -- callability/shape coverage only (numerical correctness of the effect-size
    # criterion itself is Fortran-only, per this project's testing convention) for the new
    # Issue #178 diagnostics: plateau_mode as a mode string, the new threshold optionals, and the
    # trace_* outputs' presence, dtype and shape. Reuses Test 2's single-candidate setup so the
    # search actually reaches and bootstraps a candidate.
    # ============================================================
    result3 = run_js_comp_test_parameter_search(gene_means_2, residuals_2, shared_residual_range=1.0,
                                                 n_bootstraps=5, join_method='join_min',
                                                 min_residuals_per_bin=0, min_neighbor_overlap=0.0,
                                                 plateau_mode='plateau_effect_size', delta_median_threshold=0.05,
                                                 delta_max_threshold=0.10, delta_epsilon=1e-10,
                                                 delta_min_consecutive_transitions=2, random_seed=1)

    # n_admissible_evaluated is DM_RESULT_SIZE_IS's own count argument, dropped from the Python
    # return since every trace_* array already comes back trimmed to exactly that length.
    for key in ("trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower",
                "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_delta",
                "trace_delta_median", "trace_delta_max", "plateau_established"):
        assert key in result3, f"missing expected output key '{key}'"
    assert isinstance(result3["plateau_established"], (bool, np.bool_)), \
        f"expected plateau_established to be a bool, got {type(result3['plateau_established'])}"
    n_admissible = result3["trace_n_points"].shape[-1]
    assert n_admissible >= 1, f"expected at least one admissible candidate, got {n_admissible}"
    assert result3["trace_global_js_divergence"].shape[-1] == n_admissible
    assert result3["trace_ci_width"].shape[-1] == n_admissible
    assert result3["trace_ci_width_relative"].shape == result3["trace_ci_width"].shape
    assert result3["trace_ci_width"].dtype == np.float64
    assert result3["trace_delta"].dtype == np.float64
    # The first (and here, only) admissible candidate has no predecessor to diff against.
    np.testing.assert_array_almost_equal(result3["trace_delta"][:, 0], [-1.0, -1.0], decimal=12)


def main():
    run_all_tests(globals().values())


if __name__ == "__main__":
    main()
