source("rcpp/test_helpers.R")
source("rcpp/tensoromics_functions.R")

## R tests for tox angle / RAP pipeline
## Uses custom assertion helpers provided by the user.

EPS <- 1e-12
PI  <- base::pi

ERR_OK = 0L
ERR_INVALID_INPUT = 201L
ERR_NO_STABLE_DIRECTION = 401L
ERR_NO_ANGULAR_VARIATION = 402L

## ------------------------------------------------------------------
## 0. Unit-length normalization tests
## ------------------------------------------------------------------

test_normalize_vectors_unit_length_basic <- function() {
  n_samples <- 3L
  n_genes   <- 3L

  expression_vectors <- matrix(c(
    1.0, 2.0, 2.0,   # Gene 1: norm = 3.0
    0.0, 1.0, 0.0,   # Gene 2: norm = 1.0
    1.0, 0.0, 0.0    # Gene 3: norm = 1.0
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  res <- tox_normalize_vectors_unit_length(expression_vectors)

  assert_equal_int(res$ierr, ERR_OK, "normalize_vectors_unit_length_basic: ierr mismatch")

  unit_vectors <- res$unit_vectors
  for (i in seq_len(n_genes)) {
    norm_i <- sqrt(sum(unit_vectors[, i]^2))
    assert_true(abs(norm_i - 1.0) < EPS,
               sprintf("normalize_vectors_unit_length_basic: column %d not normalized (norm=%f)", i, norm_i))
  }
}

test_normalize_vectors_unit_length_zero_norm <- function() {
  n_samples <- 3L
  n_genes   <- 2L

  expression_vectors <- matrix(c(
    1.0, 2.0, 2.0,   # Gene 1: norm = 3.0
    0.0, 0.0, 0.0    # Gene 2: ZERO VECTOR
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  res <- tryCatch({
    tox_normalize_vectors_unit_length(expression_vectors)
  }, error = function(e) e)

  assert_true(inherits(res, "error"),
             "normalize_vectors_unit_length_zero_norm: expected an error for a zero-norm vector")
}

## ------------------------------------------------------------------
## 1. Spherical pipeline tests
## ------------------------------------------------------------------

test_compute_family_direction_basic <- function() {
  n_samples  <- 3L
  n_genes    <- 6L
  n_families <- 2L

  unit_vectors <- matrix(c(
    1.0, 0.0, 0.0,   # Gene 1
    0.9, 0.1, 0.0,   # Gene 2
    0.8, 0.2, 0.0,   # Gene 3
    0.0, 1.0, 0.0,   # Gene 4
    0.1, 0.9, 0.0,   # Gene 5
    0.2, 0.8, 0.0    # Gene 6
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  gene_to_fam <- c(1L, 1L, 1L, 2L, 2L, 2L)

  res <- tox_compute_family_direction(
    unit_vectors = unit_vectors,
    gene_to_fam  = gene_to_fam,
    min_angular_dispersion = 0,
    max_angular_dispersion = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "compute_family_direction_basic: ierr mismatch")
  assert_equal_int(res$status, ERR_OK, "compute_family_direction_basic: status mismatch")

  fam_dir <- res$family_directions
  ang_disp <- res$angular_dispersions

  norm1 <- sqrt(sum(fam_dir[, 1]^2))
  norm2 <- sqrt(sum(fam_dir[, 2]^2))

  assert_equal_numeric(norm1, 1.0, tol = EPS, msg = "compute_family_direction_basic: family 1 not unit")
  assert_equal_numeric(norm2, 1.0, tol = EPS, msg = "compute_family_direction_basic: family 2 not unit")

  assert_true(ang_disp[1] > 0 && ang_disp[1] < PI,
              "compute_family_direction_basic: family 1 dispersion out of range")
  assert_true(ang_disp[2] > 0 && ang_disp[2] < PI,
              "compute_family_direction_basic: family 2 dispersion out of range")

  assert_true(fam_dir[1, 1] > 0.9, "compute_family_direction_basic: family 1 direction incorrect")
  assert_true(abs(fam_dir[2, 1]) < 0.2, "compute_family_direction_basic: family 1 direction incorrect")
}

test_compute_family_direction_single_family <- function() {
  n_samples  <- 3L
  n_genes    <- 4L
  n_families <- 1L

  unit_vectors <- matrix(c(
    1.0, 0.0, 0.0,
    0.9, 0.1, 0.0,
    0.8, 0.2, 0.0,
    0.7, 0.3, 0.0
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  gene_to_fam <- rep(1L, n_genes)

  res <- tox_compute_family_direction(
    unit_vectors = unit_vectors,
    gene_to_fam  = gene_to_fam,
    min_angular_dispersion = 0,
    max_angular_dispersion = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "compute_family_direction_single_family: ierr mismatch")
  assert_true(res$angular_dispersions[1] > 0,
              "compute_family_direction_single_family: dispersion should be positive")
}

test_compute_family_direction_no_stable_direction <- function() {
  n_samples  <- 3L
  n_genes    <- 4L
  n_families <- 1L

  unit_vectors <- matrix(c(
    1.0,  0.0, 0.0,
   -1.0,  0.0, 0.0,
    0.0,  1.0, 0.0,
    0.0, -1.0, 0.0
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  gene_to_fam <- rep(1L, n_genes)

  res <- tox_compute_family_direction(
    unit_vectors = unit_vectors,
    gene_to_fam  = gene_to_fam,
    min_angular_dispersion = 0,
    max_angular_dispersion = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "compute_family_direction_no_stable_direction: ierr mismatch")
  assert_equal_numeric(res$angular_dispersions[1], -1.0, tol = EPS,
                       msg = "compute_family_direction_no_stable_direction: dispersion should be -1")
  assert_equal_int(res$status, ERR_NO_STABLE_DIRECTION,
                   "compute_family_direction_no_stable_direction: status mismatch")
}

test_compute_angles_to_direction_basic <- function() {
  n_samples  <- 3L
  n_genes    <- 4L
  n_families <- 2L

  unit_vectors <- matrix(c(
    1.0, 0.0, 0.0,
    0.8, 0.6, 0.0,
    0.0, 1.0, 0.0,
    0.6, 0.8, 0.0
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  family_directions <- matrix(c(
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0
  ), nrow = n_samples, ncol = n_families, byrow = FALSE)

  gene_to_fam <- c(1L, 1L, 2L, 2L)

  result <- tox_compute_angles_to_direction(
    unit_vectors      = unit_vectors,
    family_directions = family_directions,
    gene_to_fam       = gene_to_fam
  )

  ierr <- result$ierr
  assert_equal_int(ierr, ERR_OK, "compute_angles_to_direction_basic: ierr mismatch")

  assert_equal_numeric(result$angles[1], 0.0, tol = EPS,
                       msg = "compute_angles_to_direction_basic: gene 1 angle mismatch")
  assert_equal_numeric(result$angles[2], acos(0.8), tol = EPS,
                       msg = "compute_angles_to_direction_basic: gene 2 angle mismatch")
  assert_equal_numeric(result$angles[3], 0.0, tol = EPS,
                       msg = "compute_angles_to_direction_basic: gene 3 angle mismatch")
  assert_equal_numeric(result$angles[4], acos(0.8), tol = EPS,
                       msg = "compute_angles_to_direction_basic: gene 4 angle mismatch")
}

test_compute_angles_to_direction_invalid_family <- function() {
  n_samples  <- 3L
  n_genes    <- 4L
  n_families <- 2L

  unit_vectors <- matrix(c(
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
    1.0 / sqrt(3), 1.0 / sqrt(3), 1.0 / sqrt(3)
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  family_directions <- matrix(c(
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0
  ), nrow = n_samples, ncol = n_families, byrow = FALSE)

  gene_to_fam <- c(1L, 0L, -3L, 2L)

  assert_error(
    tox_compute_angles_to_direction(
      unit_vectors      = unit_vectors,
      family_directions = family_directions,
      gene_to_fam       = gene_to_fam
    )
  , "compute_angles_to_direction_invalid_family: Expected error")
}

test_z_scores_by_dispersion_basic <- function() {
  n_genes    <- 5L
  n_families <- 2L

  angles <- c(0.2, 0.4, -1.0, 0.3, 0.6)
  angular_dispersions <- c(0.1, 0.2)
  gene_to_fam <- c(1L, 1L, 0L, 2L, 2L)

  result <- tox_z_scores_by_dispersion(
    angles             = angles,
    gene_to_fam        = gene_to_fam,
    angular_dispersions = angular_dispersions
  )

  ierr <- result$ierr
  assert_equal_int(ierr, ERR_OK, "z_scores_by_dispersion_basic: ierr mismatch")

  assert_equal_numeric(result$z_scores[1], 2.0, tol = EPS,
                       msg = "z_scores_by_dispersion_basic: gene 1 mismatch")
  assert_equal_numeric(result$z_scores[2], 4.0, tol = EPS,
                       msg = "z_scores_by_dispersion_basic: gene 2 mismatch")
  assert_equal_numeric(result$z_scores[3], -1.0, tol = EPS,
                       msg = "z_scores_by_dispersion_basic: gene 3 mismatch")
  assert_equal_numeric(result$z_scores[4], 1.5, tol = EPS,
                       msg = "z_scores_by_dispersion_basic: gene 4 mismatch")
  assert_equal_numeric(result$z_scores[5], 3.0, tol = EPS,
                       msg = "z_scores_by_dispersion_basic: gene 5 mismatch")
}

test_z_scores_by_dispersion_invalid <- function() {
  n_genes    <- 4L
  n_families <- 2L

  angles <- c(-1.0, 0.5, 0.3, 0.4)
  angular_dispersions <- c(-1.0, 0.0)
  gene_to_fam <- c(1L, 1L, -3L, 2L)

  assert_error(
    tox_z_scores_by_dispersion(
      angles             = angles,
      gene_to_fam        = gene_to_fam,
      angular_dispersions = angular_dispersions
    )
  , "z_scores_by_dispersion_invalid: Expected error")
}

test_angle_outliers_basic <- function() {
  z_scores <- c(1.0, 1.5, 2.0, 2.5, 3.0, 2.5, 2.0, -1.0)
  percentile    <- 0.8

  res <- tox_angle_outliers(
    z_scores = z_scores,
    percentile    = percentile
  )

  assert_equal_int(res$ierr, ERR_OK, "angle_outliers_basic: ierr mismatch")

  threshold <- res$threshold
  is_outlier <- as.logical(res$is_outlier)

  assert_equal_numeric(threshold, 2.5, tol = EPS,
                       msg = "angle_outliers_basic: threshold mismatch")

  assert_false(is_outlier[1], "angle_outliers_basic: gene 1 should not be outlier")
  assert_false(is_outlier[2], "angle_outliers_basic: gene 2 should not be outlier")
  assert_false(is_outlier[3], "angle_outliers_basic: gene 3 should not be outlier")
  assert_true(is_outlier[4],  "angle_outliers_basic: gene 4 should be outlier")
  assert_true(is_outlier[5],  "angle_outliers_basic: gene 5 should be outlier")
  assert_true(is_outlier[6],  "angle_outliers_basic: gene 6 should be outlier")
  assert_false(is_outlier[7], "angle_outliers_basic: gene 7 should not be outlier")
  assert_false(is_outlier[8], "angle_outliers_basic: gene 8 should not be outlier")
}

test_angle_outliers_no_valid <- function() {
  z_scores <- c(-1.0, -1.0, -1.0)
  percentile    <- 0.9

  res <- tox_angle_outliers(
    z_scores = z_scores,
    percentile    = percentile
  )

  assert_equal_int(res$ierr, ERR_OK, "angle_outliers_no_valid: ierr mismatch")
  is_outlier <- as.logical(res$is_outlier)
  assert_false(any(is_outlier), "angle_outliers_no_valid: no outliers expected")
}

test_angle_outliers_all_outliers <- function() {
  z_scores <- c(1.0, 1.1, 1.2, 1.3, 1.4)

  res0 <- tox_angle_outliers(
    z_scores = z_scores,
    percentile    = 0.0
  )
  assert_equal_int(res0$ierr, ERR_OK, "angle_outliers_all_outliers: ierr mismatch (0%)")
  assert_true(all(as.logical(res0$is_outlier)),
              "angle_outliers_all_outliers: all should be outliers at 0th percentile")

  res100 <- tox_angle_outliers(
    z_scores = z_scores,
    percentile    = 1.0
  )
  assert_equal_int(res100$ierr, ERR_OK, "angle_outliers_all_outliers: ierr mismatch (100%)")
  assert_false(all(as.logical(res100$is_outlier)),
               "angle_outliers_all_outliers: not all should be outliers at 100th percentile")
}

test_detect_angle_outliers_basic <- function() {
  n_samples  <- 3L
  n_genes    <- 10L
  n_families <- 2L

  expression_vectors <- matrix(c(
    1.0, 0.1, 0.0,
    0.9, 0.2, 0.0,
    0.8, 0.3, 0.0,
    0.7, 0.4, 0.0,
    0.0, 1.0, 0.0,
    0.1, 1.0, 0.1,
    0.2, 0.9, 0.1,
    0.3, 0.8, 0.1,
    0.4, 0.7, 0.1,
    1.0, 0.0, 0.0
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  gene_to_fam <- c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L)

  res <- tox_detect_angle_outliers_pipeline(
    expression_vectors      = expression_vectors,
    gene_to_fam             = gene_to_fam,
    percentile              = 0.8,
    min_angular_dispersion  = 0,
    max_angular_dispersion  = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "detect_angle_outliers_basic: ierr mismatch")

  is_outlier <- as.logical(res$is_outlier)

  assert_true(is_outlier[5],  "detect_angle_outliers_basic: gene 5 should be outlier")
  assert_true(is_outlier[10], "detect_angle_outliers_basic: gene 10 should be outlier")
  assert_false(any(is_outlier[1:4]), "detect_angle_outliers_basic: genes 1-4 should not be outliers")
  assert_false(any(is_outlier[6:9]), "detect_angle_outliers_basic: genes 6-9 should not be outliers")
}

test_detect_angle_outliers_single_family <- function() {
  n_samples  <- 3L
  n_genes    <- 8L
  n_families <- 1L

  expression_vectors <- matrix(c(
    1.0, 0.1, 0.0,
    0.9, 0.2, 0.0,
    0.8, 0.3, 0.0,
    0.7, 0.4, 0.0,
    0.6, 0.5, 0.0,
    0.5, 0.6, 0.0,
    0.4, 0.7, 0.0,
    0.0, 1.0, 0.0
  ), nrow = n_samples, ncol = n_genes, byrow = FALSE)

  gene_to_fam <- rep(1L, n_genes)

  res <- tox_detect_angle_outliers_pipeline(
    expression_vectors      = expression_vectors,
    gene_to_fam             = gene_to_fam,
    percentile              = 0.9,
    min_angular_dispersion  = 0,
    max_angular_dispersion  = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "detect_angle_outliers_single_family: ierr mismatch")

  is_outlier <- as.logical(res$is_outlier)
  assert_true(is_outlier[8], "detect_angle_outliers_single_family: gene 8 should be outlier")
  assert_false(any(is_outlier[1:7]),
               "detect_angle_outliers_single_family: genes 1-7 should not be outliers")
}

## ------------------------------------------------------------------
## 2. RAP pipeline tests
## ------------------------------------------------------------------

test_compute_family_direction_rap_basic <- function() {
  n_genes    <- 6L
  n_families <- 2L

  rap_angles <- c(
    0.5,
    0.6,
    0.4,
    PI / 2,
    PI / 2 + 0.1,
    PI / 2 - 0.1
  )

  gene_to_fam <- c(1L, 1L, 1L, 2L, 2L, 2L)

  res <- tox_compute_family_direction_rap(
    rap_angles            = rap_angles,
    gene_to_fam           = gene_to_fam,
    min_family_dispersion = 0,
    max_family_dispersion = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "compute_family_direction_rap_basic: ierr mismatch")

  fam_mean <- res$family_mean_angles
  fam_disp <- res$family_dispersions

  assert_equal_numeric(fam_mean[1], 0.5, tol = EPS,
                       msg = "compute_family_direction_rap_basic: family 1 mean mismatch")
  assert_equal_numeric(fam_mean[2], PI / 2, tol = EPS,
                       msg = "compute_family_direction_rap_basic: family 2 mean mismatch")

  assert_true(fam_disp[1] > 0 && fam_disp[1] < PI,
              "compute_family_direction_rap_basic: family 1 dispersion out of range")
  assert_true(fam_disp[2] > 0 && fam_disp[2] < PI,
              "compute_family_direction_rap_basic: family 2 dispersion out of range")
}

test_compute_family_direction_rap_single_family <- function() {
  n_genes    <- 5L
  n_families <- 1L

  rap_angles <- c(
    PI / 4,
    PI / 4 + 0.05,
    PI / 4 - 0.05,
    PI / 4 + 0.03,
    PI / 4 - 0.03
  )

  gene_to_fam <- rep(1L, n_genes)

  res <- tox_compute_family_direction_rap(
    rap_angles            = rap_angles,
    gene_to_fam           = gene_to_fam,
    min_family_dispersion = 0,
    max_family_dispersion = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "compute_family_direction_rap_single_family: ierr mismatch")
  assert_equal_numeric(res$family_mean_angles[1], PI / 4, tol = EPS,
                       msg = "compute_family_direction_rap_single_family: mean mismatch")
  assert_true(res$family_dispersions[1] > 0,
              "compute_family_direction_rap_single_family: dispersion should be positive")
}

test_compute_family_direction_rap_no_stable_direction <- function() {
  n_genes    <- 4L
  n_families <- 1L

  rap_angles <- c(
    0.0,
    (175 / 180) * PI,
    0.0,
    (175 / 180) * PI
  )

  gene_to_fam <- rep(1L, n_genes)

  res <- tox_compute_family_direction_rap(
    rap_angles            = rap_angles,
    gene_to_fam           = gene_to_fam,
    min_family_dispersion = 0.01,
    max_family_dispersion = sqrt(-2.0 * log(0.5))
  )

  assert_equal_int(res$ierr, ERR_OK, "compute_family_direction_rap_no_stable_direction: ierr mismatch")
  assert_equal_numeric(res$family_dispersions[1], -1.0, tol = EPS,
                       msg = "compute_family_direction_rap_no_stable_direction: dispersion should be -1")
  assert_equal_int(res$status, ERR_NO_STABLE_DIRECTION,
                   "compute_family_direction_rap_no_stable_direction: status mismatch")
}

test_compute_family_direction_rap_all_same_direction <- function() {
  n_genes    <- 5L
  n_families <- 1L

  rap_angles <- rep(PI / 3, n_genes)
  gene_to_fam <- rep(1L, n_genes)

  res <- tox_compute_family_direction_rap(
    rap_angles            = rap_angles,
    gene_to_fam           = gene_to_fam,
    min_family_dispersion = 0.01,
    max_family_dispersion = sqrt(-2.0 * log(0.5))
  )

  assert_equal_int(res$ierr, ERR_OK, "compute_family_direction_rap_all_same_direction: ierr mismatch")
  assert_equal_numeric(res$family_mean_angles[1], PI / 3, tol = EPS,
                       msg = "compute_family_direction_rap_all_same_direction: mean mismatch")
  assert_equal_numeric(res$family_dispersions[1], -1.0, tol = EPS,
                       msg = "compute_family_direction_rap_all_same_direction: dispersion should be -1")
  assert_equal_int(res$status, ERR_NO_ANGULAR_VARIATION,
                   "compute_family_direction_rap_all_same_direction: status mismatch")
}

test_compute_angular_deviations_rap_basic <- function() {
  n_genes    <- 5L
  n_families <- 2L

  rap_angles <- c(0.0, 0.2, 1.5, 1.7, -1.0)
  family_mean_angles <- c(0.1, 1.6)
  gene_to_fam <- c(1L, 1L, 2L, 2L, 0L)

  result <- tox_compute_angular_deviations_rap(
    rap_angles         = rap_angles,
    family_mean_angles = family_mean_angles,
    gene_to_fam        = gene_to_fam
  )

  ierr <- result$ierr
  assert_equal_int(ierr, ERR_OK, "compute_angular_deviations_rap_basic: ierr mismatch")

  assert_equal_numeric(result$angular_deviations[1], 0.1, tol = EPS,
                       msg = "compute_angular_deviations_rap_basic: gene 1 mismatch")
  assert_equal_numeric(result$angular_deviations[2], 0.1, tol = EPS,
                       msg = "compute_angular_deviations_rap_basic: gene 2 mismatch")
  assert_equal_numeric(result$angular_deviations[3], 0.1, tol = EPS,
                       msg = "compute_angular_deviations_rap_basic: gene 3 mismatch")
  assert_equal_numeric(result$angular_deviations[4], 0.1, tol = EPS,
                       msg = "compute_angular_deviations_rap_basic: gene 4 mismatch")
  assert_equal_numeric(result$angular_deviations[5], -1.0, tol = EPS,
                       msg = "compute_angular_deviations_rap_basic: gene 5 mismatch")
}

test_compute_angular_deviations_rap_invalid_family <- function() {
  n_genes    <- 4L
  n_families <- 2L

  rap_angles <- c(0.1, 0.2, 0.3, 0.4)
  family_mean_angles <- c(0.15, 0.35)
  gene_to_fam <- c(1L, 1L, -3L, 2L)

  assert_error(
    tox_compute_angular_deviations_rap(
      rap_angles         = rap_angles,
      family_mean_angles = family_mean_angles,
      gene_to_fam        = gene_to_fam
    )
  , "compute_angular_deviations_rap_invalid_family: Expected error")
}

test_z_scores_by_dispersion_rap_basic <- function() {
  n_genes    <- 5L
  n_families <- 2L

  angular_deviations <- c(0.2, 0.4, -1.0, 0.3, 0.6)
  family_dispersions <- c(0.1, 0.2)
  gene_to_fam <- c(1L, 1L, 0L, 2L, 2L)

  result <- tox_z_scores_by_dispersion_rap(
    angular_deviations = angular_deviations,
    family_dispersions = family_dispersions,
    gene_to_fam        = gene_to_fam
  )

  ierr <- result$ierr
  assert_equal_int(ierr, ERR_OK, "z_scores_by_dispersion_rap_basic: ierr mismatch")

  assert_equal_numeric(result$z_scores[1], 2.0, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_basic: gene 1 mismatch")
  assert_equal_numeric(result$z_scores[2], 4.0, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_basic: gene 2 mismatch")
  assert_equal_numeric(result$z_scores[3], -1.0, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_basic: gene 3 mismatch")
  assert_equal_numeric(result$z_scores[4], 1.5, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_basic: gene 4 mismatch")
  assert_equal_numeric(result$z_scores[5], 3.0, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_basic: gene 5 mismatch")
}

test_z_scores_by_dispersion_rap_invalid <- function() {
  n_genes    <- 4L
  n_families <- 2L

  angular_deviations <- c(-1.0, 0.5, 0.3, 0.4)
  family_dispersions <- c(0.0, 0.1)
  gene_to_fam <- c(1L, 1L, -3L, 2L)

  assert_error(
    tox_z_scores_by_dispersion_rap(
      angular_deviations = angular_deviations,
      family_dispersions = family_dispersions,
      gene_to_fam        = gene_to_fam
    )
  , "z_scores_by_dispersion_rap_invalid: Expected error")
}

test_z_scores_by_dispersion_rap_zero_dispersion <- function() {
  n_genes    <- 3L
  n_families <- 1L

  angular_deviations <- c(0.1, 0.2, 0.3)
  family_dispersions <- c(0.0)
  gene_to_fam <- rep(1L, n_genes)

  result <- tox_z_scores_by_dispersion_rap(
    angular_deviations = angular_deviations,
    family_dispersions = family_dispersions,
    gene_to_fam        = gene_to_fam
  )

  ierr <- result$ierr
  assert_equal_int(ierr, ERR_OK, "z_scores_by_dispersion_rap_zero_dispersion: ierr mismatch")

  assert_equal_numeric(result$z_scores[1], -1.0, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_zero_dispersion: gene 1 mismatch")
  assert_equal_numeric(result$z_scores[2], -1.0, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_zero_dispersion: gene 2 mismatch")
  assert_equal_numeric(result$z_scores[3], -1.0, tol = EPS,
                       msg = "z_scores_by_dispersion_rap_zero_dispersion: gene 3 mismatch")
}

test_angle_outliers_rap_basic <- function() {
  z_scores <- c(1.0, 1.5, 2.0, 2.5, 3.0, 2.5, 2.0, -1.0)
  percentile    <- 0.8

  res <- tox_angle_outliers_rap(
    z_scores = z_scores,
    percentile    = percentile
  )

  assert_equal_int(res$ierr, ERR_OK, "angle_outliers_rap_basic: ierr mismatch")

  threshold <- res$threshold
  is_outlier <- as.logical(res$is_outlier)

  assert_equal_numeric(threshold, 2.5, tol = EPS,
                       msg = "angle_outliers_rap_basic: threshold mismatch")

  assert_false(is_outlier[1], "angle_outliers_rap_basic: gene 1 should not be outlier")
  assert_false(is_outlier[2], "angle_outliers_rap_basic: gene 2 should not be outlier")
  assert_false(is_outlier[3], "angle_outliers_rap_basic: gene 3 should not be outlier")
  assert_true(is_outlier[4],  "angle_outliers_rap_basic: gene 4 should be outlier")
  assert_true(is_outlier[5],  "angle_outliers_rap_basic: gene 5 should be outlier")
  assert_true(is_outlier[6],  "angle_outliers_rap_basic: gene 6 should be outlier")
  assert_false(is_outlier[7], "angle_outliers_rap_basic: gene 7 should not be outlier")
  assert_false(is_outlier[8], "angle_outliers_rap_basic: gene 8 should not be outlier")
}

test_angle_outliers_rap_no_valid <- function() {
  z_scores <- c(-1.0, -1.0, -1.0)
  percentile    <- 0.9

  res <- tox_angle_outliers_rap(
    z_scores = z_scores,
    percentile    = percentile
  )

  assert_equal_int(res$ierr, ERR_OK, "angle_outliers_rap_no_valid: ierr mismatch")

  is_outlier <- as.logical(res$is_outlier)
  assert_false(any(is_outlier), "angle_outliers_rap_no_valid: no outliers expected")
}

test_angle_outliers_rap_all_outliers <- function() {
  z_scores <- c(1.0, 1.1, 1.2, 1.3, 1.4)

  res0 <- tox_angle_outliers_rap(
    z_scores = z_scores,
    percentile    = 0.0
  )
  assert_equal_int(res0$ierr, ERR_OK, "angle_outliers_rap_all_outliers: ierr mismatch (0%)")
  assert_true(all(as.logical(res0$is_outlier)),
              "angle_outliers_rap_all_outliers: all should be outliers at 0th percentile")

  res100 <- tox_angle_outliers_rap(
    z_scores = z_scores,
    percentile    = 1.0
  )
  assert_equal_int(res100$ierr, ERR_OK, "angle_outliers_rap_all_outliers: ierr mismatch (100%)")
  assert_false(all(as.logical(res100$is_outlier)),
               "angle_outliers_rap_all_outliers: not all should be outliers at 100th percentile")
}

test_detect_angle_outliers_pipeline_rap_basic <- function() {
  n_genes    <- 10L
  n_families <- 2L

  rap_angles <- c(
    0.5, 0.6, 0.4, 0.7, 1.5,
    1.6, 1.7, 1.5, 1.6, 0.0
  )

  gene_to_fam <- c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L)

  res <- tox_detect_angle_outliers_pipeline_rap(
    rap_angles            = rap_angles,
    gene_to_fam           = gene_to_fam,
    percentile            = 0.85,
    min_family_dispersion = 0,
    max_family_dispersion = PI
  )

  assert_equal_int(res$ierr, ERR_OK, "detect_angle_outliers_pipeline_rap_basic: ierr mismatch")

  is_outlier <- as.logical(res$is_outlier)
  # We don't have the exact Fortran pattern here (truncated), but we can at least
  # assert that some outliers exist and not all are outliers.
  assert_true(any(is_outlier), "detect_angle_outliers_pipeline_rap_basic: at least one outlier expected")
  assert_false(all(is_outlier), "detect_angle_outliers_pipeline_rap_basic: not all should be outliers")
}

run_all_tests()
# force(test_compute_family_direction_rap_no_stable_direction())