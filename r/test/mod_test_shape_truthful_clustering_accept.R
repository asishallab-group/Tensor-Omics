source("r/load_tensor_omics.R")
source("r/test_helpers.R")

.identity_history <- function(d_dim, o, d_value) {
  U <- diag(d_dim)
  U_history <- array(0.0, dim = c(d_dim, d_dim, o))
  U_history[, , 1] <- U
  d_history <- integer(o)
  d_history[1] <- d_value
  list(U = U, U_history = U_history, d_history = d_history)
}

# Identical bases at every retained iteration (U_first, U_history, U_tp1), identical d, G,
# normal_error: all four criteria trivially satisfied, even at zero tolerance for d, G, RMSE.
test_accept_ensemble_identical <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 5.0, 1.0, f$U, 2, 5.0, 1.0,
                                  0.1, 0, 0.0, 0.0)
  assert_true(is_accepted)
}

# A 60-degree rotation of a 1D tangent basis (U_first = U_history(1), isolating the
# chordal-distance formula/threshold from the cumulative-drift machinery): chordal distance
# = sin(60deg) ~ 0.866, rejected at a fraction-of-range of 0.5.
test_accept_ensemble_chordal_exceeds_max <- function() {
  theta <- pi / 3.0
  U_t <- diag(2)
  U_tp1 <- rbind(c(cos(theta), 0.0), c(sin(theta), 1.0))
  U_history <- array(0.0, dim = c(2, 2, 1))
  U_history[, , 1] <- U_t
  d_history <- c(1L)

  is_accepted <- accept_ensemble(U_t, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.5, 0, 1e10, 1e10)
  assert_true(!is_accepted)
}

# The same 60-degree rotation, accepted once the fraction-of-range is raised to 0.9
# (sin(60deg) ~ 0.866 <= 0.9).
test_accept_ensemble_chordal_within_max <- function() {
  theta <- pi / 3.0
  U_t <- diag(2)
  U_tp1 <- rbind(c(cos(theta), 0.0), c(sin(theta), 1.0))
  U_history <- array(0.0, dim = c(2, 2, 1))
  U_history[, , 1] <- U_t
  d_history <- c(1L)

  is_accepted <- accept_ensemble(U_t, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.9, 0, 1e10, 1e10)
  assert_true(is_accepted)
}

# The P5 regression this whole redesign fixes: a candidate only 5 degrees from the most
# recently accepted state (U_history(1), at 75deg) -- which a step-to-step-only check would
# accept -- but 80 degrees from the ensemble's own bootstrap state (U_first, at 0deg).
test_accept_ensemble_rejects_cumulative_drift_from_first <- function() {
  U_first <- diag(2)
  theta75 <- 75.0 * pi / 180.0
  U_history <- array(0.0, dim = c(2, 2, 1))
  U_history[, , 1] <- rbind(c(cos(theta75), 0.0), c(sin(theta75), 1.0))
  d_history <- c(1L)

  theta80 <- 80.0 * pi / 180.0
  U_tp1 <- rbind(c(cos(theta80), 0.0), c(sin(theta80), 1.0))

  is_accepted <- accept_ensemble(U_first, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.5, 0, 1e10, 1e10)
  assert_true(!is_accepted)
}

# Companion: small drift from both U_first (6deg) and U_history(1) (3deg) is accepted.
test_accept_ensemble_accepts_small_drift_from_both <- function() {
  U_first <- diag(2)
  theta3 <- 3.0 * pi / 180.0
  U_history <- array(0.0, dim = c(2, 2, 1))
  U_history[, , 1] <- rbind(c(cos(theta3), 0.0), c(sin(theta3), 1.0))
  d_history <- c(1L)

  theta6 <- 6.0 * pi / 180.0
  U_tp1 <- rbind(c(cos(theta6), 0.0), c(sin(theta6), 1.0))

  is_accepted <- accept_ensemble(U_first, 1, U_history, d_history, 1, 1.0, 1.0, U_tp1, 1, 1.0, 1.0,
                                  0.5, 0, 1e10, 1e10)
  assert_true(is_accepted)
}

# d_to_last = |d_tp1 - d_history(1)| = 0 (fine on its own), but d_to_first =
# |d_tp1 - d_first| = 2 > d_max=1: a d_to_last-only check would wrongly accept this.
test_accept_ensemble_d_to_first_exceeds_dmax <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 0, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 1.0, 1.0,
                                  0.9, 1, 1e10, 1e10)
  assert_true(!is_accepted)
}

# Same setup, accepted once d_max is raised to 2 (>= max(d_to_first, d_to_last) = 2).
test_accept_ensemble_d_two_fold_within_dmax <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 0, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 1.0, 1.0,
                                  0.9, 2, 1e10, 1e10)
  assert_true(is_accepted)
}

# G changes by a factor of 10 (ln(10) ~ 2.303): rejected at G_max=1.0.
test_accept_ensemble_g_ratio_exceeds_max <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 10.0, 1.0,
                                  0.9, 0, 1.0, 1e10)
  assert_true(!is_accepted)
}

# Same 10x change in G, accepted once G_max is raised past ln(10) ~ 2.303.
test_accept_ensemble_g_ratio_within_max <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 10.0, 1.0,
                                  0.9, 0, 3.0, 1e10)
  assert_true(is_accepted)
}

# normal_error changes by a factor of 10 (RMSE by sqrt(10) ~ 3.162, log ~ 1.151): rejected at
# RMSE_change_max=1.0.
test_accept_ensemble_rmse_ratio_exceeds_max <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 1.0, 10.0,
                                  0.9, 0, 1e10, 1.0)
  assert_true(!is_accepted)
}

# Same 10x change in normal_error, accepted once RMSE_change_max is raised past
# |log(sqrt(10))| ~ 1.151.
test_accept_ensemble_rmse_ratio_within_max <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 1.0, 10.0,
                                  0.9, 0, 1e10, 1.2)
  assert_true(is_accepted)
}

test_accept_ensemble_nonpositive_g <- function() {
  f <- .identity_history(3, 1, 2)

  assert_error(accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 0.0, 1.0, f$U, 2, 1.0, 1.0,
                               0.9, 0, 1e10, 1e10),
               "Expected error for G_t <= 0", ERR_INVALID_INPUT)
}

# Negative normal_error is physically impossible (a sum of eigenvalues) and must be rejected
# by validation -- unlike G_t, zero itself is valid (see the next test).
test_accept_ensemble_nonpositive_normal_error <- function() {
  f <- .identity_history(3, 1, 2)

  assert_error(accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 1.0, -1.0, f$U, 2, 1.0, 1.0,
                               0.9, 0, 1e10, 1e10),
               "Expected error for normal_error_t < 0", ERR_INVALID_INPUT)
}

# Zero normal_error at both t and t+1 (a perfectly flat/collinear ensemble) must not crash or
# spuriously reject via a log(0/0) -- the +epsilon guard inside the RMSE ratio keeps it
# well-defined and accepted.
test_accept_ensemble_zero_normal_error_is_accepted <- function() {
  f <- .identity_history(3, 1, 2)

  is_accepted <- accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 1.0, 0.0, f$U, 2, 1.0, 0.0,
                                  0.9, 0, 1e10, 1e10)
  assert_true(is_accepted)
}

test_accept_ensemble_d_first_out_of_range <- function() {
  f <- .identity_history(3, 1, 2)

  assert_error(accept_ensemble(f$U, 4, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 1.0, 1.0,
                               0.9, 0, 1e10, 1e10),
               "Expected error for d_first > n_dimensions", ERR_INVALID_INPUT)
}

test_accept_ensemble_chordal_frac_out_of_range <- function() {
  f <- .identity_history(3, 1, 2)

  assert_error(accept_ensemble(f$U, 2, f$U_history, f$d_history, 1, 1.0, 1.0, f$U, 2, 1.0, 1.0,
                               1.5, 0, 1e10, 1e10),
               "Expected error for chordal_dist_max_as_prcnt_of_range > 1", ERR_INVALID_INPUT)
}

test_accept_ensemble_history_len_out_of_range <- function() {
  f <- .identity_history(3, 1, 2)

  assert_error(accept_ensemble(f$U, 2, f$U_history, f$d_history, 2, 1.0, 1.0, f$U, 2, 1.0, 1.0,
                               0.9, 0, 1e10, 1e10),
               "Expected error for history_len > o", ERR_INVALID_INPUT)
}

run_all_tests()
