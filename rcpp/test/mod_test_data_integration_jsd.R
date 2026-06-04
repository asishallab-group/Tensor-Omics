source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

TOL <- 1e-12

test_determine_shared_residual_range <- function() {

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

  R <- tox_determine_shared_residual_range(S1, S2, 95)
  assert_true(approx_equal(R, 10.65), "Test 1 failed: expected ~10.65")

  # Test 2 — Custom quantile
  R <- tox_determine_shared_residual_range(S1, S2, 50)
  assert_true(approx_equal(R, 4.0), "Test 2 failed: expected ~4.0")

  # Test 3 — Quantile < 0 → error
  assert_error(
    tox_determine_shared_residual_range(S1, S2, -1),
    "Test 3 failed: expected error for negative quantile"
  )

  # Test 4 — Quantile > 100 → error
  assert_error(
    tox_determine_shared_residual_range(S1, S2, 150),
    "Test 4 failed: expected error for quantile > 100"
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

  R <- tox_determine_shared_residual_range(S1, S2, 95)
  assert_true(approx_equal(R, 11.0), "Test 5 failed: expected ~11.0")

  # Test 6 — All zeros
  S1 <- array(0, dim = c(4,2,2))
  S2 <- array(0, dim = c(3,2,2))
  R <- tox_determine_shared_residual_range(S1, S2, 95)
  assert_true(approx_equal(R, 0.0), "Test 6 failed: expected 0")

  # Test 7 — Single residual
  S1 <- array(3, dim = c(1, 1, 1))
  S2 <- array(-4, dim = c(1, 1, 1))
  R <- tox_determine_shared_residual_range(S1, S2, 95)
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
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assert_true(approx_equal(R, 10.65), "Test 1 failed")

  # Test 2
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 50)
  assert_true(approx_equal(R, 4.0), "Test 2 failed")

  # Test 3
  assert_error(
    tox_determine_shared_residual_range_expert(pp$pool, pp$perm, -1),
    "Test 3 failed"
  )

  # Test 4
  assert_error(
    tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 150),
    "Test 4 failed"
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
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assert_true(approx_equal(R, 11.0), "Test 5 failed")

  # Test 6 — All zeros
  S1 <- array(0, dim = c(4,2,2))
  S2 <- array(0, dim = c(3,2,2))
  pp <- make_pool(S1, S2)
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assert_true(approx_equal(R, 0.0), "Test 6 failed")

  # Test 7 — Single residual
  S1 <- array(3, dim = c(1, 1, 1))
  S2 <- array(-4, dim = c(1, 1, 1))
  pp <- make_pool(S1, S2)
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assert_true(approx_equal(R, 3.95), "Test 7 failed")
}

test_tox_build_residual_histograms <- function() {
    n_reps      <- 3
    n_neighbors <- 2
    n_points    <- 3
    n_bins      <- 4
    Rval        <- 2.0

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

    out <- tox_build_residual_histograms_filtered(
        E, Rval, n_bins, neighbor_mask = neighbor_mask
    )

    counts   <- out$counts
    pmf      <- out$pmf
    included <- out$included_n_residuals

    assert_true(all(counts == 0), "All counts should be zero")
    assert_true(all(abs(pmf - 0.0) < 1e-12), "All pmfs should be zero")
    assert_true(all(included == 0), "All included should be zero")

    filtered <- function(E, Rval, n_bins) {
        tox_build_residual_histograms_filtered(
            E, Rval, n_bins,
            neighbor_mask = matrix(TRUE, nrow = dim(E)[3], ncol=dim(E)[2], byrow=TRUE)
        )
    }

    for (func in list(tox_build_residual_histograms, filtered)) {

        E <- array(0, dim = c(n_reps, n_neighbors, n_points))
        E[,1,1] <- c(-2.0, -0.5, 0.2)
        E[,2,1] <- c( 1.7,  0.9, -1.2)
        E[,,2]  <- 0.0
        E[,1,3] <- c( 2.5, -3.0, 1.2)
        E[,2,3] <- c( 0.4, -0.1, 0.0)

        out <- func(E, Rval, n_bins)
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

        out <- func(E, Rval, n_bins)
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

        out <- func(E, Rval, n_bins)
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

        out <- func(E, Rval, n_bins)
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

  jsd <- tox_compute_divergence_per_reference_point(p, q)
  assert_true(all(abs(jsd) < TOL), "Test 1 failed")

  # Test 2 — disjoint PMFs
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,1] <- 1
  q[1,2] <- 1

  jsd <- tox_compute_divergence_per_reference_point(p, q)
  assert_true(approx_equal(jsd[1], log(2)), "Test 2 failed")
  assert_true(jsd[2] == 0, "Test 2 row 2 failed")
  assert_true(jsd[3] == 0, "Test 2 row 3 failed")

  # Test 3 — partial overlap
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,] <- c(0.5,0.5,0,0)
  q[1,] <- c(0,1,0,0)

  jsd <- tox_compute_divergence_per_reference_point(p, q)

  expected <- 0.5 * (
    0.5 * log(2) +
    0.5 * log(2/3) +
    log(1/0.75)
  )

  assert_true(approx_equal(jsd[1], expected), "Test 3 failed")

  # Test 4 — zero-probability bins
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,1] <- 1
  q[1,3] <- 1

  jsd <- tox_compute_divergence_per_reference_point(p, q)
  assert_true(approx_equal(jsd[1], log(2)), "Test 4 failed")

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

  jsd <- tox_compute_divergence_per_reference_point(p, q)

  assert_true(approx_equal(jsd[2], 0.5 * log(2)), "Test 5 row 2 failed")
  assert_true(approx_equal(jsd[3], 0.5 * log(2)), "Test 5 row 3 failed")
}

test_compute_weighted_global_divergence <- function() {
  TOL <- 1e-12
  approx_equal <- function(a, b) abs(a - b) < TOL

  # Test 1 — uniform weights
  jsd <- c(0.1,0.2,0.3,0.4)
  n1 <- c(5L,5L,5L,5L)
  n2 <- c(5L,5L,5L,5L)
  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  assert_true(all(out$weights == 0.25), "Test 1 weights mismatch")
  assert_true(approx_equal(out$global_js_divergence, 0.25), "Test 1 global JSD mismatch")

  # Test 2 — unequal sample counts
  jsd <- c(1,2,3,4)
  n1 <- c(10L,20L,30L,40L)
  n2 <- c(0L,10L,10L,10L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(10,30,40,50) / 130
  assert_true(all(abs(out$weights - expected_w) < TOL), "Test 2 weights mismatch")

  expected <- sum(jsd * expected_w)
  assert_true(approx_equal(out$global_js_divergence, expected), "Test 2 global JSD mismatch")

  # Test 3 — zero-sample neighborhoods
  jsd <- c(0.5,1.0,2.0,4.0)
  n1 <- c(0L,10L,0L,5L)
  n2 <- c(0L,0L,0L,5L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(0,0.5,0,0.5)
  assert_true(all(abs(out$weights - expected_w) < TOL), "Test 3 weights mismatch")
  assert_true(approx_equal(out$global_js_divergence, 2.5), "Test 3 global JSD mismatch")

  # Test 4 — all zero samples
  jsd <- c(1,2,3,4)
  n1 <- c(0L,0L,0L,0L)
  n2 <- c(0L,0L,0L,0L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  assert_true(all(out$weights == 0), "Test 4 weights mismatch")
  assert_true(out$global_js_divergence == 0, "Test 4 global JSD mismatch")

  # Test 5 — mixed jsd values, mixed sample counts
  jsd <- c(0.0, 0.5, 1.0, 2.0)
  n1  <- c(5L, 0L, 10L, 5L)
  n2  <- c(5L, 5L,  0L, 5L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(10, 5, 10, 10) / 35
  assert_true(all(abs(out$weights - expected_w) < TOL),
             "Test 5 failed: weights mismatch")

  expected <- sum(jsd * expected_w)
  assert_true(abs(out$global_js_divergence - expected) < TOL,
             "Test 5 failed: global JSD mismatch")
}

test_gjct_permutation_test <- function() {
  n_reps_S1 <- 4L
  n_reps_S2 <- 3L
  n_neighbors <- 1L
  n_points <- 2L
  n_permutations <- 2L
  n_bins <- 4L
  random_seed <- 666L
  huge <- .Machine$double.xmax

  # Base vector (same as Fortran)
  S_12 <- c(
    1, 2, 3, 4,
    5, 6, -7, 8,
    2, -4, 6, 8,
    1, 3
  )

  # Split into S1 and S2 arrays
  S1_vec <- S_12[seq_len(n_reps_S1 * n_neighbors * n_points)]
  S2_vec <- S_12[(n_reps_S1 * n_neighbors * n_points + 1L):length(S_12)]

  S1_arr <- array(S1_vec, dim = c(n_reps_S1, n_neighbors, n_points))
  S2_arr <- array(S2_vec, dim = c(n_reps_S2, n_neighbors, n_points))

  # Case A: all null >= observed → p = 1
  res_p1 <- tox_gjct_permutation_test_filtered(
    S1_arr, S2_arr,
    global_jsd_observed = huge,
    n_bins = n_bins,
    shared_residual_range = 10,
    n_permutations = n_permutations,
    random_seed = random_seed,
    neighbor_mask_S1 = matrix(FALSE, nrow = n_neighbors, ncol=n_points),
    neighbor_mask_S2 = matrix(FALSE, nrow = n_neighbors, ncol=n_points)
  )

  assert_true(abs(res_p1$p_value - 1/3) < 1e-12, "Test 3A: p-value should be 1 when for huge observed JSD and without included residuals")
  print(args(tox_gjct_permutation_test_filtered))
  filtered <- function(S1_arr, S2_arr, global_jsd_observed, n_bins, shared_residual_range, n_permutations, random_seed) {
      tox_gjct_permutation_test_filtered(
          S1_arr,
          S2_arr,
          global_jsd_observed=global_jsd_observed,
          n_bins=n_bins,
          shared_residual_range=shared_residual_range,
          n_permutations=n_permutations,
          random_seed=random_seed,
          neighbor_mask_S1 = matrix(TRUE, nrow = dim(S1_arr)[3], ncol=dim(S1_arr)[2], byrow=TRUE),
          neighbor_mask_S2 = matrix(TRUE, nrow = dim(S2_arr)[3], ncol=dim(S2_arr)[2], byrow=TRUE)
      )
  }

  for (func in list(tox_gjct_permutation_test, filtered)) {

    # Base vector (same as Fortran)
    S_12 <- c(
      1, 2, 3, 4,
      5, 6, -7, 8,
      2, -4, 6, 8,
      1, 3
    )

    # Split into S1 and S2 arrays
    S1_vec <- S_12[seq_len(n_reps_S1 * n_neighbors * n_points)]
    S2_vec <- S_12[(n_reps_S1 * n_neighbors * n_points + 1L):length(S_12)]

    S1_arr <- array(S1_vec, dim = c(n_reps_S1, n_neighbors, n_points))
    S2_arr <- array(S2_vec, dim = c(n_reps_S2, n_neighbors, n_points))

    # Case A: all null >= observed → p = 1
    res_p1 <- func(
      S1_arr, S2_arr,
      global_jsd_observed = 0,
      n_bins = n_bins,
      shared_residual_range = 10,
      n_permutations = n_permutations,
      random_seed = random_seed
    )
    check_err_code(res_p1$ierr)

    assert_true(abs(res_p1$p_value - 1.0) < 1e-12,
               "Test 3A: p-value should be 1 when observed JSD = 0")

  #   # Case B: none null >= observed → p = 1/(n_permutations+1)
  #   res_p2 <- func(
  #     S1_arr, S2_arr,
  #     global_jsd_observed = huge,
  #     n_bins = n_bins,
  #     shared_residual_range = 10,
  #     n_permutations = n_permutations,
  #     random_seed = random_seed
  #   )
  #   check_err_code(res_p2$ierr)

  #   expected <- 1.0 / (n_permutations + 1.0)
  #   assert_true(abs(res_p2$p_value - expected) < 1e-12,
  #              "Test 3B: p-value should be 1/(n_permutations+1) for huge observed JSD")

  }

  invisible(TRUE)
}

run_all_tests()