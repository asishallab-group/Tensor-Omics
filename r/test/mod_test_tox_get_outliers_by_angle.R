# Binding contract of tox_get_outliers_by_angle, the angle-based outlier detection.
#
# Only what the R layer adds is tested here: every published procedure can be called, returns
# the documented type, shape and names, passes one value through unchanged, and signals the
# documented error for each reachable bad argument. The numbers themselves are pinned by the
# Fortran suite, test/mod_test_tox_get_outliers_by_angle.F90.
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12
FAMILY_NAMES <- c("family_directions", "angular_dispersions", "member_counts", "status")
FAMILY_NAMES_RAP <- c("family_mean_angles", "angular_dispersions", "member_counts", "status")
PIPELINE_TAIL <- c("relative_angular_deviations", "threshold", "is_outlier", "gene_status")

# One family of four genes whose unit vectors are e1, e1, e2, e2: direction (1,1)/sqrt(2),
# dispersion sqrt(ln 2); a zero gene, and a gene without a family. Columns are genes.
EXPRESSION <- matrix(c(1, 0, 3, 0, 0, 2, 0, 0.5, 0, 0, 5, 5), nrow = 2)
GENE_TO_FAM <- c(1L, 1L, 1L, 1L, 1L, 0L)
# signed angles 0, pi/2, 0, pi/2 of one family, plus a gene without a family
SIGNED_ANGLES <- c(0, pi / 2, 0, pi / 2, 1)
GENE_TO_FAM_RAP <- c(1L, 1L, 1L, 1L, 0L)

# ---------------------------------------------------------------------------------------------
# the steps
# ---------------------------------------------------------------------------------------------

test_compute_family_direction_contract <- function() {
  # n_families is the caller's: a trailing family without genes is still reported
  res <- compute_family_direction(2, EXPRESSION, GENE_TO_FAM)
  assert_true(is.list(res) && setequal(names(res), FAMILY_NAMES), "family direction: names")
  assert_true(is.matrix(res$family_directions) && all(dim(res$family_directions) == c(2, 2)),
              "family direction: directions are a 2 x 2 matrix")
  assert_true(is.double(res$angular_dispersions) && length(res$angular_dispersions) == 2,
              "family direction: one dispersion per family")
  assert_equal_int(res$member_counts, c(4L, 0L), "family direction: the zero gene is not counted")
  assert_equal_int(res$status, c(0L, STAT_TOO_FEW_MEMBERS), "family direction: status")
  assert_equal_numeric(res$family_directions[, 1], c(1, 1) / sqrt(2), TOL, "family direction: direction")
  assert_equal_numeric(res$angular_dispersions[1], sqrt(log(2)), TOL, "family direction: dispersion")

  # the bounds reach Fortran: a maximum below sqrt(ln 2) = 0.83 makes the family unstable
  res <- compute_family_direction(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion = 0.8)
  assert_equal_int(res$status, STAT_NO_STABLE_DIRECTION, "family direction: max passed through")
  res <- compute_family_direction(1, EXPRESSION, GENE_TO_FAM, min_angular_dispersion = 0.9)
  assert_equal_int(res$status, STAT_NO_ANGULAR_VARIATION, "family direction: min passed through")
}

test_compute_family_direction_errors <- function() {
  assert_error(compute_family_direction(0, EXPRESSION, GENE_TO_FAM), "n_families 0", ERR_EMPTY_INPUT)
  assert_error(compute_family_direction(1, EXPRESSION, c(1, 1, 1, 1, 2, 0)), "gene_to_fam above n_families",
               ERR_INVALID_INPUT)
  assert_error(compute_family_direction(1, EXPRESSION, c(1, 1, 1, 1, -1, 0)), "gene_to_fam below 1",
               ERR_INVALID_INPUT)
  bad <- EXPRESSION
  bad[1, 1] <- NaN
  assert_error(compute_family_direction(1, bad, GENE_TO_FAM), "NaN expression", ERR_NAN_INF)
  bad[1, 1] <- Inf
  assert_error(compute_family_direction(1, bad, GENE_TO_FAM), "infinite expression", ERR_NAN_INF)
  assert_error(compute_family_direction(1, EXPRESSION, GENE_TO_FAM, min_angular_dispersion = -0.1), "negative min",
               ERR_INVALID_INPUT)
  assert_error(compute_family_direction(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion = -0.1), "negative max",
               ERR_INVALID_INPUT)
  assert_error(compute_family_direction(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion = 5.5),
               "max above the limit 5", ERR_INVALID_INPUT)
  assert_error(compute_family_direction(1, EXPRESSION, GENE_TO_FAM[1:3]), "gene_to_fam shorter than the genes")
}

test_compute_angular_deviations_contract <- function() {
  directions <- matrix(c(1, 0), nrow = 2)
  deviations <- compute_angular_deviations(EXPRESSION, directions, GENE_TO_FAM)
  assert_true(is.double(deviations) && length(deviations) == 6, "angular deviations: one per gene")
  # (0, 2) against e1 is orthogonal; the zero gene and the gene without a family get the sentinel
  assert_equal_numeric(deviations[3], pi / 2, TOL, "angular deviations: orthogonal")
  assert_true(identical(deviations[5:6], c(-1, -1)), "angular deviations: sentinels")

  assert_error(compute_angular_deviations(EXPRESSION, matrix(c(1.5, 0), nrow = 2), GENE_TO_FAM),
               "direction above 1", ERR_INVALID_INPUT)
  assert_error(compute_angular_deviations(EXPRESSION, matrix(c(NaN, 0), nrow = 2), GENE_TO_FAM), "NaN direction",
               ERR_NAN_INF)
  assert_error(compute_angular_deviations(EXPRESSION, directions, c(1, 1, 1, 1, 2, 0)),
               "gene_to_fam above the families in family_directions", ERR_INVALID_INPUT)
}

test_compute_family_direction_rap_contract <- function() {
  res <- compute_family_direction_rap(2, SIGNED_ANGLES, GENE_TO_FAM_RAP)
  assert_true(is.list(res) && setequal(names(res), FAMILY_NAMES_RAP), "rap family: names")
  assert_equal_int(res$member_counts, c(4L, 0L), "rap family: member counts")
  assert_equal_int(res$status, c(0L, STAT_TOO_FEW_MEMBERS), "rap family: status")
  assert_equal_numeric(res$family_mean_angles[1], pi / 4, TOL, "rap family: mean angle")
  assert_true(identical(res$family_mean_angles[2], -4), "rap family: a family without a mean gets the sentinel")

  assert_error(compute_family_direction_rap(0, SIGNED_ANGLES, GENE_TO_FAM_RAP), "n_families 0", ERR_EMPTY_INPUT)
  # the range is (-pi, pi]: -pi itself is rejected, also for a gene without a family
  assert_error(compute_family_direction_rap(1, c(0, 0.1, 0.2, 0.3, -pi), GENE_TO_FAM_RAP), "angle -pi",
               ERR_INVALID_INPUT)
  assert_error(compute_family_direction_rap(1, c(0, 0.1, 0.2, 4, 0), GENE_TO_FAM_RAP), "angle above pi",
               ERR_INVALID_INPUT)
  assert_error(compute_family_direction_rap(1, c(0, 0.1, NaN, 0.3, 0), GENE_TO_FAM_RAP), "NaN angle", ERR_NAN_INF)
  assert_error(compute_family_direction_rap(1, SIGNED_ANGLES, c(1, 1, 1, 2, 0)), "gene_to_fam above",
               ERR_INVALID_INPUT)
  assert_error(compute_family_direction_rap(1, SIGNED_ANGLES, GENE_TO_FAM_RAP, min_angular_dispersion = -1),
               "negative min", ERR_INVALID_INPUT)
}

test_compute_angular_deviations_rap_contract <- function() {
  deviations <- compute_angular_deviations_rap(SIGNED_ANGLES, c(pi, -4), c(1, 1, 1, 1, 2))
  assert_true(is.double(deviations) && length(deviations) == 5, "rap deviations: one per gene")
  # 0 against pi is half a turn; a family without a mean gives the sentinel
  assert_equal_numeric(deviations[1], pi, TOL, "rap deviations: half a turn")
  assert_true(identical(deviations[5], -1), "rap deviations: sentinel")

  assert_error(compute_angular_deviations_rap(SIGNED_ANGLES, -pi, GENE_TO_FAM_RAP), "mean -pi", ERR_INVALID_INPUT)
  assert_error(compute_angular_deviations_rap(SIGNED_ANGLES, Inf, GENE_TO_FAM_RAP), "infinite mean", ERR_NAN_INF)
  assert_error(compute_angular_deviations_rap(c(0, 0, 0, 0, 5), 0, GENE_TO_FAM_RAP), "angle above pi",
               ERR_INVALID_INPUT)
  assert_error(compute_angular_deviations_rap(SIGNED_ANGLES, 0, c(1, 1, 1, 3, 0)), "gene_to_fam above",
               ERR_INVALID_INPUT)
}

test_compute_relative_angular_deviations_contract <- function() {
  relative <- compute_relative_angular_deviations(c(1, 0.75, -1), c(0.5, 0.25), c(1, 2, 1))
  assert_true(identical(relative, c(2, 3, -1)), "relative: exact quotients and the sentinel")

  assert_error(compute_relative_angular_deviations(c(1, 3.5, 0), c(0.5, 0.25), c(1, 2, 1)), "deviation above pi",
               ERR_INVALID_INPUT)
  assert_error(compute_relative_angular_deviations(c(1, -0.5, 0), c(0.5, 0.25), c(1, 2, 1)), "negative deviation",
               ERR_INVALID_INPUT)
  assert_error(compute_relative_angular_deviations(c(1, 0.5, 0), c(0.5, -0.25), c(1, 2, 1)), "negative dispersion",
               ERR_INVALID_INPUT)
  assert_error(compute_relative_angular_deviations(c(1, 0.5, 0), c(0.5, NaN), c(1, 2, 1)), "NaN dispersion",
               ERR_NAN_INF)
  assert_error(compute_relative_angular_deviations(c(1, 0.5, 0), c(0.5, 0.25), c(1, 3, 1)), "gene_to_fam above",
               ERR_INVALID_INPUT)
}

test_compute_angle_outlier_threshold_contract <- function() {
  # type-7 quantile of 0.5, 1, 2, 4 at 0.5 is 1.5; the sentinel does not count
  threshold <- compute_angle_outlier_threshold(c(0.5, 4, -1, 1, 2), quantile_level = 0.5)
  assert_true(is.double(threshold) && length(threshold) == 1 && threshold == 1.5, "threshold: 1.5")
  # the default level is 0.95
  assert_true(identical(compute_angle_outlier_threshold(c(0.5, 4, -1, 1, 2)),
                        compute_angle_outlier_threshold(c(0.5, 4, -1, 1, 2), quantile_level = 0.95)),
              "threshold: default level 0.95")

  # a percentage is not a level
  assert_error(compute_angle_outlier_threshold(c(0.5, 1), quantile_level = 95), "level 95", ERR_INVALID_INPUT)
  assert_error(compute_angle_outlier_threshold(c(0.5, 1), quantile_level = -0.1), "level below 0", ERR_INVALID_INPUT)
  assert_error(compute_angle_outlier_threshold(c(0.5, -0.5)), "negative value", ERR_INVALID_INPUT)
  assert_error(compute_angle_outlier_threshold(c(0.5, NaN)), "NaN value", ERR_NAN_INF)
  assert_error(compute_angle_outlier_threshold(numeric(0)), "no values", ERR_EMPTY_INPUT)
}

test_flag_angle_outliers_contract <- function() {
  is_outlier <- flag_angle_outliers(c(-1, 0, 0.5, 3), -5)
  # never the sentinel or a zero, whatever the threshold
  assert_true(identical(is_outlier, c(FALSE, FALSE, TRUE, TRUE)), "flag: logical vector, the flag rule")

  assert_error(flag_angle_outliers(c(0.5, 1), NaN), "NaN threshold", ERR_NAN_INF)
  assert_error(flag_angle_outliers(c(0.5, -2), 1), "negative value", ERR_INVALID_INPUT)
}

# ---------------------------------------------------------------------------------------------
# the pipelines
# ---------------------------------------------------------------------------------------------

test_detect_angle_outliers_contract <- function() {
  res <- detect_angle_outliers(2, EXPRESSION, GENE_TO_FAM)
  assert_true(is.list(res) && setequal(names(res), c(FAMILY_NAMES, PIPELINE_TAIL)), "detect: names")
  assert_true(is.double(res$threshold) && length(res$threshold) == 1, "detect: threshold is a scalar")
  assert_true(is.logical(res$is_outlier) && length(res$is_outlier) == 6, "detect: one flag per gene")
  assert_true(all(dim(res$family_directions) == c(2, 2)), "detect: directions 2 x 2")
  assert_equal_int(res$gene_status, c(0L, 0L, 0L, 0L, STAT_ZERO_VECTOR, STAT_NO_FAMILY), "detect: gene status")
  assert_equal_int(res$status, c(0L, STAT_TOO_FEW_MEMBERS), "detect: status")

  # quantile_level reaches Fortran: at level 0 the threshold is the smallest value
  level_0 <- detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, quantile_level = 0)
  valid <- level_0$relative_angular_deviations[level_0$gene_status == 0L]
  assert_true(level_0$threshold == min(valid), "detect: level 0 passed through")

  assert_error(detect_angle_outliers(0, EXPRESSION, GENE_TO_FAM), "n_families 0", ERR_EMPTY_INPUT)
  assert_error(detect_angle_outliers(1, EXPRESSION, c(1, 1, 1, 1, 2, 0)), "gene_to_fam above", ERR_INVALID_INPUT)
  bad <- EXPRESSION
  bad[2, 3] <- NaN
  assert_error(detect_angle_outliers(1, bad, GENE_TO_FAM), "NaN expression", ERR_NAN_INF)
  assert_error(detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, quantile_level = 95), "level 95", ERR_INVALID_INPUT)
  assert_error(detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, min_angular_dispersion = -1), "negative min",
               ERR_INVALID_INPUT)
  assert_error(detect_angle_outliers(1, EXPRESSION, GENE_TO_FAM, max_angular_dispersion = NaN), "NaN max",
               ERR_NAN_INF)
}

test_detect_angle_outliers_rap_contract <- function() {
  res <- detect_angle_outliers_rap(2, SIGNED_ANGLES, GENE_TO_FAM_RAP)
  assert_true(is.list(res) && setequal(names(res), c(FAMILY_NAMES_RAP, PIPELINE_TAIL)), "detect rap: names")
  assert_true(is.double(res$threshold) && length(res$threshold) == 1, "detect rap: threshold is a scalar")
  assert_true(is.logical(res$is_outlier) && length(res$is_outlier) == 5, "detect rap: one flag per gene")
  assert_equal_int(res$gene_status, c(0L, 0L, 0L, 0L, STAT_NO_FAMILY), "detect rap: gene status")
  assert_equal_int(res$member_counts, c(4L, 0L), "detect rap: member counts")

  assert_error(detect_angle_outliers_rap(0, SIGNED_ANGLES, GENE_TO_FAM_RAP), "n_families 0", ERR_EMPTY_INPUT)
  assert_error(detect_angle_outliers_rap(1, c(0, 0.1, 0.2, -pi, 0), GENE_TO_FAM_RAP), "angle -pi", ERR_INVALID_INPUT)
  assert_error(detect_angle_outliers_rap(1, SIGNED_ANGLES, c(1, 1, 1, 1, 2)), "gene_to_fam above", ERR_INVALID_INPUT)
  assert_error(detect_angle_outliers_rap(1, SIGNED_ANGLES, GENE_TO_FAM_RAP, quantile_level = 1.5), "level above 1",
               ERR_INVALID_INPUT)
  assert_error(detect_angle_outliers_rap(1, SIGNED_ANGLES, GENE_TO_FAM_RAP, max_angular_dispersion = -1),
               "negative max", ERR_INVALID_INPUT)
}

test_no_expert_tier_published <- function() {
  # the expert tiers only hand over work arrays, so R has none
  for (name in c("detect_angle_outliers_expert", "detect_angle_outliers_rap_expert",
                 "compute_angle_outlier_threshold_expert"))
    assert_false(exists(name), paste(name, "should not be published to R"))
}

run_all_tests()
