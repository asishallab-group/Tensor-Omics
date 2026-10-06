#!/usr/bin/env python3
"""Binding contract of tox_get_outliers_by_angle, the angle-based outlier detection.

Only what the Python layer adds is tested here: every published procedure can be called,
returns the documented type, shape and keys, passes one value through unchanged, and raises
the documented error for each reachable bad argument. The numbers themselves are pinned by
the Fortran suite, test/mod_test_tox_get_outliers_by_angle.F90.
"""

from math import pi as PI, sqrt, log
from pathlib import Path
import sys

import numpy as np

sys.path.append(str(Path(__file__).parent.parent))
from tensor_omics import (
    compute_family_direction,
    compute_angular_deviations,
    compute_family_direction_rap,
    compute_angular_deviations_rap,
    compute_relative_angular_deviations,
    compute_angle_outlier_threshold,
    flag_angle_outliers,
    detect_angle_outliers,
    detect_angle_outliers_rap,
)
from test_helpers import run_all_tests, assert_error
from tensor_omics.error_handling import (
    ERR_EMPTY_INPUT, ERR_INVALID_INPUT, ERR_NAN_INF,
    STAT_NO_STABLE_DIRECTION, STAT_NO_ANGULAR_VARIATION, STAT_TOO_FEW_MEMBERS, STAT_NO_FAMILY, STAT_ZERO_VECTOR,
)

TOL = 1e-12
FAMILY_KEYS = {"family_directions", "angular_dispersions", "member_counts", "status"}
FAMILY_KEYS_RAP = {"family_mean_angles", "angular_dispersions", "member_counts", "status"}
PIPELINE_KEYS = FAMILY_KEYS | {"relative_angular_deviations", "threshold", "is_outlier", "gene_status"}
PIPELINE_KEYS_RAP = FAMILY_KEYS_RAP | {"relative_angular_deviations", "threshold", "is_outlier", "gene_status"}

# One family of four genes whose unit vectors are e1, e1, e2, e2: direction (1,1)/sqrt(2),
# dispersion sqrt(ln 2); a zero gene, and a gene without a family. Columns are genes.
EXPRESSION = np.array([[1.0, 3.0, 0.0, 0.0, 0.0, 5.0],
                       [0.0, 0.0, 2.0, 0.5, 0.0, 5.0]])
GENE_TO_FAM = np.array([1, 1, 1, 1, 1, 0], dtype=np.int32)
# signed angles 0, pi/2, 0, pi/2 of one family, plus a gene without a family
SIGNED_ANGLES = np.array([0.0, PI / 2, 0.0, PI / 2, 1.0])
GENE_TO_FAM_RAP = np.array([1, 1, 1, 1, 0], dtype=np.int32)


def _assert_read_only(array, name):
    assert not array.flags.writeable, f"{name} should be read-only"


# ---------------------------------------------------------------------------------------------
# the steps
# ---------------------------------------------------------------------------------------------

def test_compute_family_direction_contract():
    # n_families is the caller's: a trailing family without genes is still reported
    result = compute_family_direction(2, EXPRESSION, GENE_TO_FAM)
    assert isinstance(result, dict) and set(result) == FAMILY_KEYS, f"keys: {set(result)}"
    assert result["family_directions"].shape == (2, 2) and result["family_directions"].dtype == np.float64
    assert result["angular_dispersions"].shape == (2,) and result["angular_dispersions"].dtype == np.float64
    assert result["member_counts"].dtype == np.int32 and result["status"].dtype == np.int32
    for name in FAMILY_KEYS:
        _assert_read_only(result[name], name)
    np.testing.assert_array_equal(result["member_counts"], [4, 0])  # the zero gene is not counted
    np.testing.assert_array_equal(result["status"], [0, STAT_TOO_FEW_MEMBERS])
    np.testing.assert_allclose(result["family_directions"][:, 0], [1 / sqrt(2), 1 / sqrt(2)], atol=TOL)
    np.testing.assert_allclose(result["angular_dispersions"][0], sqrt(log(2)), atol=TOL)


def test_compute_family_direction_bounds_pass_through():
    # a maximum below sqrt(ln 2) = 0.83 reaches Fortran: the family is then unstable
    status = compute_family_direction(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion=0.8)["status"]
    np.testing.assert_array_equal(status, [STAT_NO_STABLE_DIRECTION])
    status = compute_family_direction(1, EXPRESSION, GENE_TO_FAM, min_angular_dispersion=0.9)["status"]
    np.testing.assert_array_equal(status, [STAT_NO_ANGULAR_VARIATION])


def test_compute_family_direction_errors():
    assert_error(lambda: compute_family_direction(0, EXPRESSION, GENE_TO_FAM), "n_families 0", ERR_EMPTY_INPUT)
    assert_error(lambda: compute_family_direction(1, EXPRESSION, [1, 1, 1, 1, 2, 0]), "gene_to_fam above n_families",
                 ERR_INVALID_INPUT)
    assert_error(lambda: compute_family_direction(1, EXPRESSION, [1, 1, 1, 1, -1, 0]), "gene_to_fam below 1",
                 ERR_INVALID_INPUT)
    bad = EXPRESSION.copy()
    bad[0, 0] = np.nan
    assert_error(lambda: compute_family_direction(1, bad, GENE_TO_FAM), "NaN expression", ERR_NAN_INF)
    bad[0, 0] = np.inf
    assert_error(lambda: compute_family_direction(1, bad, GENE_TO_FAM), "infinite expression", ERR_NAN_INF)
    assert_error(lambda: compute_family_direction(1, EXPRESSION, GENE_TO_FAM, min_angular_dispersion=-0.1),
                 "negative min", ERR_INVALID_INPUT)
    assert_error(lambda: compute_family_direction(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion=-0.1),
                 "negative max", ERR_INVALID_INPUT)
    assert_error(lambda: compute_family_direction(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion=5.5),
                 "max above the limit 5", ERR_INVALID_INPUT)
    assert_error(lambda: compute_family_direction(1, np.empty((2, 0)), np.empty(0, dtype=np.int32)), "no genes",
                 ERR_EMPTY_INPUT)
    assert_error(lambda: compute_family_direction(1, EXPRESSION, GENE_TO_FAM[:3]), "gene_to_fam shorter than the genes")


def test_compute_angular_deviations_contract():
    directions = np.array([[1.0], [0.0]])
    deviations = compute_angular_deviations(EXPRESSION, directions, GENE_TO_FAM)
    assert isinstance(deviations, np.ndarray) and deviations.shape == (6,) and deviations.dtype == np.float64
    _assert_read_only(deviations, "angular_deviations")
    # (0, 2) against e1 is orthogonal; the zero gene and the gene without a family get the sentinel
    np.testing.assert_allclose(deviations[2], PI / 2, atol=TOL)
    np.testing.assert_array_equal(deviations[4:], [-1.0, -1.0])

    assert_error(lambda: compute_angular_deviations(EXPRESSION, np.array([[1.5], [0.0]]), GENE_TO_FAM),
                 "direction above 1", ERR_INVALID_INPUT)
    assert_error(lambda: compute_angular_deviations(EXPRESSION, np.array([[np.nan], [0.0]]), GENE_TO_FAM),
                 "NaN direction", ERR_NAN_INF)
    assert_error(lambda: compute_angular_deviations(EXPRESSION, directions, [1, 1, 1, 1, 2, 0]),
                 "gene_to_fam above the families in family_directions", ERR_INVALID_INPUT)


def test_compute_family_direction_rap_contract():
    result = compute_family_direction_rap(2, SIGNED_ANGLES, GENE_TO_FAM_RAP)
    assert isinstance(result, dict) and set(result) == FAMILY_KEYS_RAP, f"keys: {set(result)}"
    assert result["family_mean_angles"].shape == (2,) and result["family_mean_angles"].dtype == np.float64
    for name in FAMILY_KEYS_RAP:
        _assert_read_only(result[name], name)
    np.testing.assert_array_equal(result["member_counts"], [4, 0])
    np.testing.assert_array_equal(result["status"], [0, STAT_TOO_FEW_MEMBERS])
    np.testing.assert_allclose(result["family_mean_angles"][0], PI / 4, atol=TOL)
    assert result["family_mean_angles"][1] == -4.0, "a family without a mean gets the sentinel"

    assert_error(lambda: compute_family_direction_rap(0, SIGNED_ANGLES, GENE_TO_FAM_RAP), "n_families 0",
                 ERR_EMPTY_INPUT)
    # the range is (-pi, pi]: -pi itself is rejected, also for a gene without a family
    assert_error(lambda: compute_family_direction_rap(1, [0.0, 0.1, 0.2, 0.3, -PI], GENE_TO_FAM_RAP), "angle -pi",
                 ERR_INVALID_INPUT)
    assert_error(lambda: compute_family_direction_rap(1, [0.0, 0.1, 0.2, 4.0, 0.0], GENE_TO_FAM_RAP),
                 "angle above pi", ERR_INVALID_INPUT)
    assert_error(lambda: compute_family_direction_rap(1, [0.0, 0.1, np.nan, 0.3, 0.0], GENE_TO_FAM_RAP), "NaN angle",
                 ERR_NAN_INF)
    assert_error(lambda: compute_family_direction_rap(1, SIGNED_ANGLES, [1, 1, 1, 2, 0]), "gene_to_fam above",
                 ERR_INVALID_INPUT)
    assert_error(lambda: compute_family_direction_rap(1, SIGNED_ANGLES, GENE_TO_FAM_RAP, min_angular_dispersion=-1.0),
                 "negative min", ERR_INVALID_INPUT)


def test_compute_angular_deviations_rap_contract():
    deviations = compute_angular_deviations_rap(SIGNED_ANGLES, [PI, -4.0], [1, 1, 1, 1, 2])
    assert deviations.shape == (5,) and deviations.dtype == np.float64
    _assert_read_only(deviations, "angular_deviations")
    # 0 against pi is half a turn; a family without a mean gives the sentinel
    np.testing.assert_allclose(deviations[0], PI, atol=TOL)
    assert deviations[4] == -1.0

    assert_error(lambda: compute_angular_deviations_rap(SIGNED_ANGLES, [-PI], GENE_TO_FAM_RAP), "mean -pi",
                 ERR_INVALID_INPUT)
    assert_error(lambda: compute_angular_deviations_rap(SIGNED_ANGLES, [np.inf], GENE_TO_FAM_RAP), "infinite mean",
                 ERR_NAN_INF)
    assert_error(lambda: compute_angular_deviations_rap([0.0, 0.0, 0.0, 0.0, 5.0], [0.0], GENE_TO_FAM_RAP),
                 "angle above pi", ERR_INVALID_INPUT)
    assert_error(lambda: compute_angular_deviations_rap(SIGNED_ANGLES, [0.0], [1, 1, 1, 3, 0]), "gene_to_fam above",
                 ERR_INVALID_INPUT)


def test_compute_relative_angular_deviations_contract():
    relative = compute_relative_angular_deviations([1.0, 0.75, -1.0], [0.5, 0.25], [1, 2, 1])
    assert relative.shape == (3,) and relative.dtype == np.float64
    _assert_read_only(relative, "relative_angular_deviations")
    np.testing.assert_array_equal(relative, [2.0, 3.0, -1.0])

    assert_error(lambda: compute_relative_angular_deviations([1.0, 3.5, 0.0], [0.5, 0.25], [1, 2, 1]),
                 "deviation above pi", ERR_INVALID_INPUT)
    assert_error(lambda: compute_relative_angular_deviations([1.0, -0.5, 0.0], [0.5, 0.25], [1, 2, 1]),
                 "negative deviation", ERR_INVALID_INPUT)
    assert_error(lambda: compute_relative_angular_deviations([1.0, 0.5, 0.0], [0.5, -0.25], [1, 2, 1]),
                 "negative dispersion", ERR_INVALID_INPUT)
    assert_error(lambda: compute_relative_angular_deviations([1.0, 0.5, 0.0], [0.5, np.nan], [1, 2, 1]),
                 "NaN dispersion", ERR_NAN_INF)
    assert_error(lambda: compute_relative_angular_deviations([1.0, 0.5, 0.0], [0.5, 0.25], [1, 3, 1]),
                 "gene_to_fam above", ERR_INVALID_INPUT)


def test_compute_angle_outlier_threshold_contract():
    # type-7 quantile of 0.5, 1, 2, 4 at 0.5 is 1.5; the sentinel does not count
    threshold = compute_angle_outlier_threshold([0.5, 4.0, -1.0, 1.0, 2.0], quantile_level=0.5)
    assert isinstance(threshold, float) and threshold == 1.5, f"threshold: {threshold!r}"
    # the default level is 0.95
    assert compute_angle_outlier_threshold([0.5, 4.0, -1.0, 1.0, 2.0]) == \
        compute_angle_outlier_threshold([0.5, 4.0, -1.0, 1.0, 2.0], quantile_level=0.95)

    # a percentage is not a level
    assert_error(lambda: compute_angle_outlier_threshold([0.5, 1.0], quantile_level=95.0), "level 95",
                 ERR_INVALID_INPUT)
    assert_error(lambda: compute_angle_outlier_threshold([0.5, 1.0], quantile_level=-0.1), "level below 0",
                 ERR_INVALID_INPUT)
    assert_error(lambda: compute_angle_outlier_threshold([0.5, -0.5]), "negative value", ERR_INVALID_INPUT)
    assert_error(lambda: compute_angle_outlier_threshold([0.5, np.nan]), "NaN value", ERR_NAN_INF)
    assert_error(lambda: compute_angle_outlier_threshold([]), "no values", ERR_EMPTY_INPUT)


def test_flag_angle_outliers_contract():
    is_outlier = flag_angle_outliers([-1.0, 0.0, 0.5, 3.0], -5.0)
    assert isinstance(is_outlier, np.ndarray) and is_outlier.dtype == np.bool_ and is_outlier.shape == (4,)
    _assert_read_only(is_outlier, "is_outlier")
    # never the sentinel or a zero, whatever the threshold
    np.testing.assert_array_equal(is_outlier, [False, False, True, True])

    assert_error(lambda: flag_angle_outliers([0.5, 1.0], np.nan), "NaN threshold", ERR_NAN_INF)
    assert_error(lambda: flag_angle_outliers([0.5, -2.0], 1.0), "negative value", ERR_INVALID_INPUT)


# ---------------------------------------------------------------------------------------------
# the pipelines
# ---------------------------------------------------------------------------------------------

def test_detect_angle_outliers_contract():
    result = detect_angle_outliers(2, EXPRESSION, GENE_TO_FAM)
    assert isinstance(result, dict) and set(result) == PIPELINE_KEYS, f"keys: {set(result)}"
    assert isinstance(result["threshold"], float)
    assert result["relative_angular_deviations"].shape == (6,)
    assert result["is_outlier"].dtype == np.bool_ and result["is_outlier"].shape == (6,)
    assert result["gene_status"].dtype == np.int32 and result["gene_status"].shape == (6,)
    assert result["family_directions"].shape == (2, 2)
    for name in PIPELINE_KEYS - {"threshold"}:
        _assert_read_only(result[name], name)
    np.testing.assert_array_equal(result["gene_status"], [0, 0, 0, 0, STAT_ZERO_VECTOR, STAT_NO_FAMILY])
    np.testing.assert_array_equal(result["status"], [0, STAT_TOO_FEW_MEMBERS])

    # quantile_level reaches Fortran: at level 0 the threshold is the smallest value
    level_0 = detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, quantile_level=0.0)
    valid = level_0["relative_angular_deviations"][level_0["gene_status"] == 0]
    assert level_0["threshold"] == valid.min()

    assert_error(lambda: detect_angle_outliers(0, EXPRESSION, GENE_TO_FAM), "n_families 0", ERR_EMPTY_INPUT)
    assert_error(lambda: detect_angle_outliers(1, EXPRESSION, [1, 1, 1, 1, 2, 0]), "gene_to_fam above",
                 ERR_INVALID_INPUT)
    bad = EXPRESSION.copy()
    bad[1, 2] = np.nan
    assert_error(lambda: detect_angle_outliers(1, bad, GENE_TO_FAM), "NaN expression", ERR_NAN_INF)
    assert_error(lambda: detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, quantile_level=95.0), "level 95",
                 ERR_INVALID_INPUT)
    assert_error(lambda: detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, min_angular_dispersion=-1.0),
                 "negative min", ERR_INVALID_INPUT)
    assert_error(lambda: detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion=np.nan),
                 "NaN max", ERR_NAN_INF)


def test_detect_angle_outliers_rap_contract():
    result = detect_angle_outliers_rap(2, SIGNED_ANGLES, GENE_TO_FAM_RAP)
    assert isinstance(result, dict) and set(result) == PIPELINE_KEYS_RAP, f"keys: {set(result)}"
    assert isinstance(result["threshold"], float)
    assert result["family_mean_angles"].shape == (2,)
    assert result["is_outlier"].dtype == np.bool_ and result["is_outlier"].shape == (5,)
    for name in PIPELINE_KEYS_RAP - {"threshold"}:
        _assert_read_only(result[name], name)
    np.testing.assert_array_equal(result["gene_status"], [0, 0, 0, 0, STAT_NO_FAMILY])
    np.testing.assert_array_equal(result["member_counts"], [4, 0])

    assert_error(lambda: detect_angle_outliers_rap(0, SIGNED_ANGLES, GENE_TO_FAM_RAP), "n_families 0",
                 ERR_EMPTY_INPUT)
    assert_error(lambda: detect_angle_outliers_rap(1, [0.0, 0.1, 0.2, -PI, 0.0], GENE_TO_FAM_RAP), "angle -pi",
                 ERR_INVALID_INPUT)
    assert_error(lambda: detect_angle_outliers_rap(1, SIGNED_ANGLES, [1, 1, 1, 1, 2]), "gene_to_fam above",
                 ERR_INVALID_INPUT)
    assert_error(lambda: detect_angle_outliers_rap(1, SIGNED_ANGLES, GENE_TO_FAM_RAP, quantile_level=1.5),
                 "level above 1", ERR_INVALID_INPUT)
    assert_error(lambda: detect_angle_outliers_rap(1, SIGNED_ANGLES, GENE_TO_FAM_RAP, max_angular_dispersion=-1.0),
                 "negative max", ERR_INVALID_INPUT)


def test_no_expert_tier_published():
    # the expert tiers only hand over work arrays, so Python has none
    import tensor_omics
    for name in ("detect_angle_outliers_expert", "detect_angle_outliers_rap_expert",
                 "compute_angle_outlier_threshold_expert"):
        assert not hasattr(tensor_omics, name), f"{name} should not be published to Python"


if __name__ == "__main__":
    run_all_tests(globals().values())
