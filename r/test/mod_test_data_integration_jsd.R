source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12

# Test determine_shared_residual_range (the plain, auto-sorting tier) directly, by cross-checking
# it against the already-tested _expert tier given the same pool pre-sorted.
test_determine_shared_residual_range <- function() {
  pool <- c(1.0, 5.0, 3.0, 8.0, 2.0, 7.0, 4.0)
  perm <- order(pool)

  R_plain <- determine_shared_residual_range(pool, 0.95)
  R_expert <- determine_shared_residual_range_expert(pool, perm, 0.95)
  assert_true(abs(R_plain - R_expert) < TOL)

  # Default residual_range_quantile (0.95)
  R_plain_default <- determine_shared_residual_range(pool)
  assert_true(abs(R_plain_default - R_expert) < TOL)
}

test_determine_study_shared_residual_range <- function() {

  # Helper
  approx_equal <- function(a, b) abs(a - b) < TOL

  # Test 1 — Basic correctness
  S1 <- array(c(
    1,2,3,4,
    5,6,-7,8,
    9,10,11,12,
    1,1,1,1
  ), dim = c(4, 2, 2))

  S2 <- array(c(
    2,-4,6,8,
    1,3,5,7,
    9,0,1,2
  ), dim = c(3, 2, 2))

  R <- determine_study_shared_residual_range(S1, S2, 0.95)
  assert_true(approx_equal(R, 10.65), "Test 1 failed: expected ~10.65")

  # Test 2 — Custom quantile
  R <- determine_study_shared_residual_range(S1, S2, 0.5)
  assert_true(approx_equal(R, 4.0), "Test 2 failed: expected ~4.0")

  # Test 3 — Quantile < 0 → error
  assert_error(
    determine_study_shared_residual_range(S1, S2, -1),
    "Test 3 failed: expected error for negative quantile", ERR_INVALID_INPUT
  )

  # Test 4 — Quantile > 1 → error
  assert_error(
    determine_study_shared_residual_range(S1, S2, 1.5),
    "Test 4 failed: expected error for quantile > 1", ERR_INVALID_INPUT
  )

  # Test 5 — NaNs ignored
  S1 <- array(c(
    NA_real_,2,3,
    4,5,6,
    -7,8,-11,
    9,10,12,
    1,1,1,1
  ), dim = c(4,2,2))
  S2 <- array(c(
    -1,2,3,
    4,5,6,
    -7,8,-11,
    9,10,NA_real_
  ), dim = c(3,2,2))

  R <- determine_study_shared_residual_range(S1, S2, 0.95)
  assert_true(approx_equal(R, 11.0), "Test 5 failed: expected ~11.0")

  # Test 6 — All zeros
  S1 <- array(0, dim = c(4,2,2))
  S2 <- array(0, dim = c(3,2,2))
  R <- determine_study_shared_residual_range(S1, S2, 0.95)
  assert_true(approx_equal(R, 0.0), "Test 6 failed: expected 0")

  # Test 7 — Single residual
  S1 <- array(3, dim = c(1, 1, 1))
  S2 <- array(-4, dim = c(1, 1, 1))
  R <- determine_study_shared_residual_range(S1, S2, 0.95)
  assert_true(approx_equal(R, 3.95), "Test 7 failed: expected ~3.95")
}

test_determine_shared_residual_range_expert <- function() {
  TOL <- 1e-12
  approx_equal <- function(a, b) abs(a - b) < TOL

  make_pool <- function(S1, S2) {
    pool <- abs(c(S1, S2))
     perm <- order(pool)
    list(pool = pool, perm = perm)
  }

  S1 <- array(c(
    1,2,3,4,
    5,6,-7,8,
    9,10,11,12,
    1,1,1,1
  ), dim = c(4, 2, 2))

  S2 <- array(c(
    2,-4,6,8,
    1,3,5,7,
    9,0,1,2
  ), dim = c(3, 2, 2))

  pp <- make_pool(S1, S2)

  # Test 1
  R <- determine_shared_residual_range_expert(pp$pool, pp$perm, 0.95)
  assert_true(approx_equal(R, 10.65), "Test 1 failed")

  # Test 2
  R <- determine_shared_residual_range_expert(pp$pool, pp$perm, 0.5)
  assert_true(approx_equal(R, 4.0), "Test 2 failed")

  # Test 3
  assert_error(
    determine_shared_residual_range_expert(pp$pool, pp$perm, -1),
    "Test 3 failed", ERR_INVALID_INPUT
  )

  # Test 4
  assert_error(
    determine_shared_residual_range_expert(pp$pool, pp$perm, 1.5),
    "Test 4 failed", ERR_INVALID_INPUT
  )

  # Test 5 — NaNs ignored
  S1 <- array(c(
    NA_real_,2,3,
    4,5,6,
    -7,8,-11,
    9,10,12,
    1,1,1,1
  ), dim = c(4,2,2))
  S2 <- array(c(
    -1,2,3,
    4,5,6,
    -7,8,-11,
    9,10,NA_real_
  ), dim = c(3,2,2))

  pp <- make_pool(S1, S2)
  R <- determine_shared_residual_range_expert(pp$pool, pp$perm, 0.95)
  assert_true(approx_equal(R, 11.0), "Test 5 failed")

  # Test 6 — All zeros
  S1 <- array(0, dim = c(4,2,2))
  S2 <- array(0, dim = c(3,2,2))
  pp <- make_pool(S1, S2)
  R <- determine_shared_residual_range_expert(pp$pool, pp$perm, 0.95)
  assert_true(approx_equal(R, 0.0), "Test 6 failed")

  # Test 7 — Single residual
  S1 <- array(3, dim = c(1, 1, 1))
  S2 <- array(-4, dim = c(1, 1, 1))
  pp <- make_pool(S1, S2)
  R <- determine_shared_residual_range_expert(pp$pool, pp$perm, 0.95)
  assert_true(approx_equal(R, 3.95), "Test 7 failed")
}

test_tox_build_residual_histograms <- function() {
    n_reps      <- 3
    n_neighbors <- 2
    n_points    <- 3
    n_bins      <- 4
    n_bins_per_point <- rep(as.integer(n_bins), n_points)
    Rval        <- 2.0
    R_low       <- rep(-Rval, n_points)
    R_high      <- rep(Rval, n_points)

    # ============================================================
    # Test 1 — Simple symmetric case, no NaNs
    # ============================================================
    E <- array(0, dim = c(n_reps, n_neighbors, n_points))

    E[,1,1] <- c(-2.0, -0.5, 0.2)
    E[,2,1] <- c( 1.7,  0.9, -1.2)
    E[,,2]  <- 0.0
    E[,1,3] <- c( 2.5, -3.0, 1.2)
    E[,2,3] <- c( 0.4, -0.1, 0.0)

    neighbor_mask <- array(FALSE, dim = c(n_neighbors, n_points))

    out <- build_residual_histograms(
        E, R_low, R_high, n_bins, n_bins_per_point, neighbor_mask = neighbor_mask
    )

    counts   <- out$counts
    pmf      <- out$pmf
    included <- out$included_n_residuals

    assert_true(all(counts == 0), "All counts should be zero")
    assert_true(all(abs(pmf - 0.0) < 1e-12), "All pmfs should be zero")
    assert_true(all(included == 0), "All included should be zero")

    filtered <- function(E, R_low, R_high, n_bins, n_bins_per_point) {
        build_residual_histograms(
            E, R_low, R_high, n_bins, n_bins_per_point,
            neighbor_mask = matrix(TRUE, nrow = dim(E)[3], ncol=dim(E)[2], byrow=TRUE)
        )
    }

    for (func in list(build_residual_histograms, filtered)) {

        E <- array(0, dim = c(n_reps, n_neighbors, n_points))
        E[,1,1] <- c(-2.0, -0.5, 0.2)
        E[,2,1] <- c( 1.7,  0.9, -1.2)
        E[,,2]  <- 0.0
        E[,1,3] <- c( 2.5, -3.0, 1.2)
        E[,2,3] <- c( 0.4, -0.1, 0.0)

        out <- func(E, R_low, R_high, n_bins, n_bins_per_point)
        counts   <- out$counts
        pmf      <- out$pmf
        included <- out$included_n_residuals

        expected_counts <- matrix(
            c(2L,1L,2L,1L,
              0L,0L,6L,0L,
              1L,1L,2L,2L),
            nrow = n_points, byrow = TRUE
        )

        expected_pmf <- expected_counts / 6

        assert_equal_int(out$counts, expected_counts, "Test 1: counts mismatch")
        assert_equal_numeric(out$pmf, expected_pmf, "Test 1: pmf mismatch")
        assert_equal_int(out$included, as.integer(c(6,6,6)), "Test 1: included mismatch")

        # ============================================================
        # Test 2 — NaNs must be ignored
        # ============================================================
        E[,,] <- 0
        E[1,1,1] <- NaN
        E[3,1,1] <- NaN
        E[2,1,2] <- NaN
        E[3,2,3] <- NaN

        out <- func(E, R_low, R_high, n_bins, n_bins_per_point)
        counts   <- out$counts
        pmf      <- out$pmf
        included <- out$included_n_residuals

        expected_counts <- matrix(
            as.integer(
            c(0,0,4,0,
              0,0,5,0,
              0,0,5,0)
            ),
            nrow = 3, byrow = TRUE
        )

        expected_pmf <- matrix(
            c(0,0,1,0,
              0,0,1,0,
              0,0,1,0),
            nrow = 3, byrow = TRUE
        )
        assert_equal_int(out$counts, expected_counts, "Test 2: counts mismatch")
        assert_equal_int(out$pmf, expected_pmf, "Test 2: pmf mismatch")
        assert_equal_int(out$included, as.integer(c(4,5,5)), "Test 2: included mismatch")

        # ============================================================
        # Test 3 — All NaN → pmf = 0, counts = 0, included = 0
        # ============================================================
        E[,,] <- NaN

        out <- func(E, R_low, R_high, n_bins, n_bins_per_point)
        assert_true(all(out$counts == 0), "test 3: All counts should be zero")
        assert_true(all(abs(out$pmf - 0.0) < 1e-12), "test 3: All pmfs should be zero")
        assert_true(all(out$included == 0), "test 3: All included should be zero")

        # ============================================================
        # Test 4 — Residuals exactly on boundaries
        # ============================================================
        E <- array(c(
          -2,-1,0,1,2,0,
          -2,-1,0,1,2,0,
          -2,-1,0,1,2,0
        ), dim = c(3,2,3))

        out <- func(E, R_low, R_high, n_bins, n_bins_per_point)
        counts   <- out$counts
        pmf      <- out$pmf
        included <- out$included_n_residuals

        expected_counts <- matrix(
            as.integer(
            c(1,1,2,2,
              1,1,2,2,
              1,1,2,2)
            ),
            nrow = 3, byrow = TRUE
        )

        expected_pmf <- expected_counts / 6

        assert_equal_int(out$counts, expected_counts, "Test 4: counts mismatch")
        assert_equal_int(out$pmf, expected_pmf, "Test 4: pmf mismatch")
        assert_equal_int(out$included, as.integer(c(6,6,6)), "Test 4: included mismatch")
    }
}

test_compute_divergence_per_reference_point <- function() {
  TOL <- 1e-12
  approx_equal <- function(a, b) abs(a - b) < TOL

  # Test 1 — identical PMFs
  p <- matrix(c(
    0.1,0.2,0.3,0.4,
    0.25,0.25,0.25,0.25,
    1,0,0,0
  ), 3, 4, byrow = TRUE)

  q <- p

  jsd <- compute_divergence_per_reference_point(p, q)
  assert_true(all(abs(jsd) < TOL), "Test 1 failed")

  # Test 2 — disjoint PMFs
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,1] <- 1
  q[1,2] <- 1

  jsd <- compute_divergence_per_reference_point(p, q)
  assert_true(approx_equal(jsd[1], log(2) / log(2)), "Test 2 failed")  # rescaled onto 0..1: log(2)/log(2) = 1.0
  assert_true(jsd[2] == 0, "Test 2 row 2 failed")
  assert_true(jsd[3] == 0, "Test 2 row 3 failed")

  # Test 3 — partial overlap
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,] <- c(0.5,0.5,0,0)
  q[1,] <- c(0,1,0,0)

  jsd <- compute_divergence_per_reference_point(p, q)

  expected <- (0.5 * (
    0.5 * log(2) +
    0.5 * log(2/3) +
    log(1/0.75)
  )) / log(2)  # rescaled onto 0..1 by dividing through LOG_2 = log(2)

  assert_true(approx_equal(jsd[1], expected), "Test 3 failed")

  # Test 4 — zero-probability bins
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,1] <- 1
  q[1,3] <- 1

  jsd <- compute_divergence_per_reference_point(p, q)
  assert_true(approx_equal(jsd[1], log(2) / log(2)), "Test 4 failed")  # rescaled onto 0..1: log(2)/log(2) = 1.0

  # Test 5 — mixed patterns
  p <- matrix(
    c(
      0.2, 0.3, 0.5,
      0.0, 1.0, 0.0,
      0.0, 0.0, 0.25,
      0.25, 0.25, 0.25
    ),
    nrow = 3,
    ncol = 4
  )

  q <- matrix(
    c(
      0.2, 0.3, 0.5,
      0.0, 0.0, 1.0,
      0.0, 0.0, 0.25,
      0.25, 0.25, 0.25
    ),
    nrow = 3,
    ncol = 4
  )

  jsd <- compute_divergence_per_reference_point(p, q)

  assert_true(approx_equal(jsd[2], (0.5 * log(2)) / log(2)), "Test 5 row 2 failed")  # rescaled: 0.5*log(2)/log(2) = 0.5
  assert_true(approx_equal(jsd[3], (0.5 * log(2)) / log(2)), "Test 5 row 3 failed")
}

test_compute_weighted_global_divergence <- function() {
  TOL <- 1e-12
  approx_equal <- function(a, b) abs(a - b) < TOL

  # Test 1 — uniform weights
  jsd <- c(0.1,0.2,0.3,0.4)
  n1 <- c(5L,5L,5L,5L)
  n2 <- c(5L,5L,5L,5L)
  out <- compute_weighted_global_divergence(jsd, n1, n2)

  assert_true(all(out$weights == 0.25), "Test 1 weights mismatch")
  assert_true(approx_equal(out$global_js_divergence, 0.25), "Test 1 global JSD mismatch")

  # Test 2 — unequal sample counts
  jsd <- c(1,2,3,4)
  n1 <- c(10L,20L,30L,40L)
  n2 <- c(0L,10L,10L,10L)

  out <- compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(10,30,40,50) / 130
  assert_true(all(abs(out$weights - expected_w) < TOL), "Test 2 weights mismatch")

  expected <- sum(jsd * expected_w)
  assert_true(approx_equal(out$global_js_divergence, expected), "Test 2 global JSD mismatch")

  # Test 3 — zero-sample neighborhoods
  jsd <- c(0.5,1.0,2.0,4.0)
  n1 <- c(0L,10L,0L,5L)
  n2 <- c(0L,0L,0L,5L)

  out <- compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(0,0.5,0,0.5)
  assert_true(all(abs(out$weights - expected_w) < TOL), "Test 3 weights mismatch")
  assert_true(approx_equal(out$global_js_divergence, 2.5), "Test 3 global JSD mismatch")

  # Test 4 — all zero samples
  jsd <- c(1,2,3,4)
  n1 <- c(0L,0L,0L,0L)
  n2 <- c(0L,0L,0L,0L)

  out <- compute_weighted_global_divergence(jsd, n1, n2)

  assert_true(all(out$weights == 0), "Test 4 weights mismatch")
  assert_true(out$global_js_divergence == 0, "Test 4 global JSD mismatch")

  # Test 5 — mixed jsd values, mixed sample counts
  jsd <- c(0.0, 0.5, 1.0, 2.0)
  n1  <- c(5L, 0L, 10L, 5L)
  n2  <- c(5L, 5L,  0L, 5L)

  out <- compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(10, 5, 10, 10) / 35
  assert_true(all(abs(out$weights - expected_w) < TOL),
             "Test 5 failed: weights mismatch")

  expected <- sum(jsd * expected_w)
  assert_true(abs(out$global_js_divergence - expected) < TOL,
             "Test 5 failed: global JSD mismatch")
}

test_calc_pmf <- function() {
  TOL <- 1e-12

  # ============================================================
  # Test 1 — matches build_residual_histograms's own basic fixture (see
  # test_tox_build_residual_histograms Test 1 above): feeding it the exact counts/
  # included_n_reps that routine produces must reproduce that fixture's expected pmf exactly.
  # ============================================================
  counts <- matrix(
      as.integer(
      c(2,1,2,1,
        0,0,6,0,
        1,1,2,2)
      ),
      nrow = 3, byrow = TRUE
  )
  included_n_reps <- as.integer(c(6, 6, 6))

  expected_pmf <- counts / matrix(included_n_reps, nrow = 3, ncol = 4)

  pmf <- calc_pmf(counts, included_n_reps)
  assert_equal_numeric(pmf, expected_pmf, tol = TOL, msg = "Test 1: pmf mismatch")

  # ============================================================
  # Test 2 — a reference point with zero included replicates (every residual there was NaN)
  # must report an all-zero pmf row rather than dividing by zero.
  # ============================================================
  counts <- matrix(as.integer(c(3, 5, 0, 0)), nrow = 2, byrow = TRUE)
  included_n_reps <- as.integer(c(8, 0))

  expected_pmf <- matrix(c(0.375, 0.625, 0.0, 0.0), nrow = 2, byrow = TRUE)

  pmf <- calc_pmf(counts, included_n_reps)
  assert_equal_numeric(pmf, expected_pmf, tol = TOL, msg = "Test 2: pmf mismatch")

  # ============================================================
  # Test 3 — counts is documented as non-negative; a negative entry must be rejected.
  # ============================================================
  counts <- matrix(as.integer(-1), nrow = 1, ncol = 1)
  included_n_reps <- as.integer(1)

  assert_error(
    calc_pmf(counts, included_n_reps),
    "Test 3 failed: expected ERR_INVALID_INPUT", ERR_INVALID_INPUT
  )
}

test_determine_all_studies_shared_residual_range <- function() {
  TOL <- 1e-12

  # ============================================================
  # Test 1 — Three single-replicate studies, hand-computed: pooled absolute residuals sorted
  # are [3, 4, 5]; the default 95% quantile has rank 0.95*(3-1)+1 = 2.9, so
  # R = sorted(2) + 0.9*(sorted(3)-sorted(2)) = 4 + 0.9*1 = 4.9
  # ============================================================
  all_studies <- array(NA_real_, dim = c(1, 1, 1, 3))
  all_studies[1, 1, 1, 1] <- 3.0
  all_studies[1, 1, 1, 2] <- -4.0
  all_studies[1, 1, 1, 3] <- 5.0

  R <- determine_all_studies_shared_residual_range(all_studies)
  assert_true(abs(R - 4.9) < TOL, "Test 1 failed: expected 4.9")

  # ============================================================
  # Test 2 — Regression-safety cross-check: feeding the same two studies both through
  # determine_study_shared_residual_range (with S2 as-is) and through the N-study routine (S2
  # padded with NA up to S1's replicate count, n_studies=2) must give the exact same range.
  # ============================================================
  S1 <- array(c(
    1,2,3,4,
    5,6,-7,8,
    9,10,11,12,
    1,1,1,1
  ), dim = c(4, 2, 2))

  S2 <- array(c(
    2,-4,6,8,
    1,3,5,7,
    9,0,1,2
  ), dim = c(3, 2, 2))

  R_two_study <- determine_study_shared_residual_range(S1, S2, 0.95)

  all_studies <- array(NA_real_, dim = c(4, 2, 2, 2))
  all_studies[, , , 1] <- S1
  all_studies[1:3, , , 2] <- S2  # all_studies[4, , , 2] stays NA -- padding S2 up to max_n_reps

  R_all_studies <- determine_all_studies_shared_residual_range(all_studies)
  assert_true(abs(R_all_studies - R_two_study) < TOL,
             "Test 2 failed: N-study range does not match two-study range")
  assert_true(abs(R_all_studies - 10.65) < TOL, "Test 2 failed: expected 10.65")
}

# test_gjct_permutation_test was removed here: gjct_permutation_test was replaced by a K-study,
# consensus-based version (see tox_data_integration_js_comp_test) -- the old 2-study test that
# lived in this file was removed rather than adapted, since the signature is unrelated. Its own
# direct test coverage now lives in mod_test_data_integration_stats.R.

run_all_tests()