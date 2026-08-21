"""
Python test suite for sample_estimator_anchors / grow_estimator_anchor_clouds /
estimate_stc_parameters (tox_shape_truthful_clustering_parameter_estimation), mirroring
test/mod_test_shape_truthful_clustering_parameter_estimation.F90
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensor_omics import sample_estimator_anchors, grow_estimator_anchor_clouds, estimate_stc_parameters, build_kd_index
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import ERR_INVALID_INPUT, ERR_INTERNAL


# =====================
# sample_estimator_anchors
# =====================
# 11 points, density labels equal to point index (1..11, already ascending -- the sort
# permutation is the identity). n_anchors=5 gives percentiles 20/40/60/80/100, ranks
# 3/5/7/9/11 exactly (no interpolation rounding needed): anchor_indices = [3,5,7,9,11]
# (1-indexed, matching the Fortran kernel's own convention).
def test_sample_anchors_hand_computed():
    density_labels = np.arange(1, 12, dtype=np.float64)
    anchors = sample_estimator_anchors(density_labels, n_anchors=5)
    expected = np.array([3, 5, 7, 9, 11], dtype=np.int32)
    assert np.array_equal(anchors, expected), f"expected {expected}, got {anchors}"


# 5 points, n_anchors=5 (every point its own percentile mark): ranks 2/3/3/4/5 -- rank 3 is
# hit twice, so anchor_indices necessarily repeats an index. Documented, not a bug.
def test_sample_anchors_duplicates_possible():
    density_labels = np.arange(1, 6, dtype=np.float64)
    anchors = sample_estimator_anchors(density_labels, n_anchors=5)
    expected = np.array([2, 3, 3, 4, 5], dtype=np.int32)
    assert np.array_equal(anchors, expected), f"expected {expected}, got {anchors}"


def test_sample_anchors_invalid_n_anchors_zero():
    density_labels = np.arange(1, 12, dtype=np.float64)
    assert_error(lambda: sample_estimator_anchors(density_labels, n_anchors=0),
                 "Expected error for n_anchors < 1", ERR_INVALID_INPUT)


def test_sample_anchors_invalid_n_anchors_too_large():
    density_labels = np.arange(1, 12, dtype=np.float64)
    assert_error(lambda: sample_estimator_anchors(density_labels, n_anchors=12),
                 "Expected error for n_anchors > n_vectors", ERR_INVALID_INPUT)


# =====================
# grow_estimator_anchor_clouds
# =====================
# D=2, N=7 points on a line, (0,0)..(6,0). Two anchors at the opposite ends, point 1 (x=0,
# 1-indexed) and point 7 (x=6).
def _line_fixture_7():
    vectors = np.zeros((2, 7), dtype=np.float64, order='F')
    vectors[0, :] = np.arange(7, dtype=np.float64)
    return vectors


# seed_max_set_size=100 (grow until every point is claimed). Every round is an exact distance
# tie between the two clouds' own nearest-unclaimed candidate (both always 1.0 apart on this
# evenly-spaced line) -- ties are broken by whichever EA is scanned first (anchor 1), so
# anchor 1 wins every single round. Anchor 2's cloud never grows past its own single point.
def test_grow_clouds_symmetric_line_ties_favor_lower_index():
    vectors = _line_fixture_7()
    anchor_indices = np.array([1, 7], dtype=np.int32)
    result = grow_estimator_anchor_clouds(vectors, anchor_indices, seed_max_set_size=100.0)

    expected_sizes = np.array([6, 1], dtype=np.int32)
    assert np.array_equal(result["cloud_sizes"], expected_sizes), \
        f"expected {expected_sizes}, got {result['cloud_sizes']}"

    expected_cloud_1 = np.zeros(7, dtype=np.bool_)
    expected_cloud_1[0:6] = True
    assert np.array_equal(result["cloud_masks"][:, 0], expected_cloud_1)
    assert result["cloud_masks"][6, 1] and int(np.sum(result["cloud_masks"][:, 1])) == 1


# Same fixture, seed_max_set_size=50 -> ceiling(0.5*7)=4 total claims: 2 anchors already
# present plus 2 more rounds, both won by EA 1.
def test_grow_clouds_seed_max_set_size_stops_early():
    vectors = _line_fixture_7()
    anchor_indices = np.array([1, 7], dtype=np.int32)
    result = grow_estimator_anchor_clouds(vectors, anchor_indices, seed_max_set_size=50.0)

    expected_sizes = np.array([3, 1], dtype=np.int32)
    assert np.array_equal(result["cloud_sizes"], expected_sizes), \
        f"expected {expected_sizes}, got {result['cloud_sizes']}"
    assert int(np.sum(result["cloud_masks"])) == 4


# Default seed_max_set_size (5.0): ceiling(0.05*7)=1, clamped up to n_anchors=2 itself --
# actual_max_claims never exceeds the anchor count, so no growth happens at all.
def test_grow_clouds_default_seed_max_set_size_can_yield_no_growth():
    vectors = _line_fixture_7()
    anchor_indices = np.array([1, 7], dtype=np.int32)
    result = grow_estimator_anchor_clouds(vectors, anchor_indices)

    expected_sizes = np.array([1, 1], dtype=np.int32)
    assert np.array_equal(result["cloud_sizes"], expected_sizes), \
        f"expected {expected_sizes}, got {result['cloud_sizes']}"


def test_grow_clouds_invalid_seed_max_set_size():
    vectors = _line_fixture_7()
    anchor_indices = np.array([1, 7], dtype=np.int32)
    assert_error(lambda: grow_estimator_anchor_clouds(vectors, anchor_indices, seed_max_set_size=150.0),
                 "Expected error for seed_max_set_size > 100", ERR_INVALID_INPUT)


# =====================
# estimate_stc_parameters
# =====================
# D=2, N=21, a perfectly collinear, evenly-spaced line (0,0)..(20,0). Every estimator
# anchor's grown cloud is itself a sub-interval of the same line, so every EA agrees exactly
# on d=1 and on tangent direction -- chordal_dist_max_as_prcnt_of_range and d_max must both
# come out at (or, for the former, an SVD-residual hair above) 0. This fixture is exactly
# symmetric, so sample_estimator_anchors_impl's own tie-break (ties resolved by ascending point
# index) is what pins the anchor set to [3,5,7,9,11] and, through it, k_min/density_quantile
# below to a single deterministic outcome. k_min/k_density/density_quantile/G_max are
# cross-checked against this exact, already-verified, fully deterministic kernel's own real
# output (no randomness anywhere here).
def _collinear_line_21():
    vectors = np.zeros((2, 21), dtype=np.float64, order='F')
    vectors[0, :] = np.arange(21, dtype=np.float64)
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    return vectors, kd_indices, dimension_order


def test_estimate_parameters_collinear_line():
    vectors, kd_indices, dimension_order = _collinear_line_21()
    result = estimate_stc_parameters(vectors, kd_indices, dimension_order, seed_max_set_size=50.0)

    assert result["estimated_chordal_dist_max_as_prcnt_of_range"] < 1e-6
    assert abs(result["estimated_d_max"] - 0.0) < 1e-9
    assert abs(result["estimated_k_min"] - 2.0) < 1e-9
    assert abs(result["estimated_k_density"] - result["estimated_k_min"]) < 1e-9
    assert abs(result["estimated_density_quantile"] - 1.0) < 1e-9
    assert abs(result["estimated_G_max"] - 0.0) < 1e-9


# seed_max_set_size=0: every EA's cloud stays size 1 -- zero clouds ever reach the size >= 2
# a genuine observable/SVD needs. Fewer than 2 usable EAs is a genuine, data-dependent
# runtime failure, not a validation error.
def test_estimate_parameters_too_few_valid_eas():
    vectors, kd_indices, dimension_order = _collinear_line_21()
    assert_error(lambda: estimate_stc_parameters(vectors, kd_indices, dimension_order, seed_max_set_size=0.0),
                 "Expected error when fewer than 2 EAs ever grow past size 1", ERR_INTERNAL)


def test_estimate_parameters_invalid_n_anchors():
    vectors, kd_indices, dimension_order = _collinear_line_21()
    assert_error(lambda: estimate_stc_parameters(vectors, kd_indices, dimension_order, n_anchors=50),
                 "Expected error for n_anchors > n_vectors", ERR_INVALID_INPUT)


def test_estimate_parameters_invalid_seed_max_set_size():
    vectors, kd_indices, dimension_order = _collinear_line_21()
    assert_error(lambda: estimate_stc_parameters(vectors, kd_indices, dimension_order, seed_max_set_size=-1.0),
                 "Expected error for seed_max_set_size < 0", ERR_INVALID_INPUT)


# The Fortran suite's test_estimate_parameters_omitted_n_anchors_is_clamped is a regression
# test for a crash that could only happen with n_anchors truly *absent* at the Fortran ABI
# boundary -- not reproducible here, since the Python binding always resolves and passes
# n_anchors=5 explicitly (see misc/code_gen_footgun.md's third entry). What is worth covering
# from Python: that this always-explicit default of 5 still gets validated normally (a clean,
# typed error, not a crash) on a dataset smaller than it.
def test_estimate_parameters_default_n_anchors_too_large_for_dataset():
    vectors = np.zeros((2, 3), dtype=np.float64, order='F')
    vectors[0, :] = np.arange(3, dtype=np.float64)
    dimension_order = np.array([1, 2], dtype=np.int32)
    kd_indices = build_kd_index(vectors, dimension_order)
    assert_error(lambda: estimate_stc_parameters(vectors, kd_indices, dimension_order),
                 "Expected error for the default n_anchors=5 on a 3-point dataset", ERR_INVALID_INPUT)


if __name__ == "__main__":
    run_all_tests(globals().values())
