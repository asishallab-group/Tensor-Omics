source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# Identical basis, identical d, identical G: all three criteria trivially satisfied, even at
# zero tolerance for d and G.
test_accept_ensemble_identical <- function() {
  U_t <- diag(3)
  U_tp1 <- U_t

  is_accepted <- accept_ensemble(U_t, 2, 5.0, U_tp1, 2, 5.0, 0.1, 0, 0.0)
  assert_true(is_accepted)
}

# A 60-degree rotation of a 1D tangent basis, rejected by an alpha_max of 30 degrees.
test_accept_ensemble_angle_exceeds_max <- function() {
  theta <- pi / 3.0
  alpha_max <- pi / 6.0

  U_t <- diag(2)
  U_tp1 <- rbind(c(cos(theta), 0.0), c(sin(theta), 1.0))

  is_accepted <- accept_ensemble(U_t, 1, 1.0, U_tp1, 1, 1.0, alpha_max, 0, 1e10)
  assert_true(!is_accepted)
}

test_accept_ensemble_angle_within_max <- function() {
  theta <- pi / 3.0
  alpha_max <- 7.0 * pi / 18.0  # 70 degrees

  U_t <- diag(2)
  U_tp1 <- rbind(c(cos(theta), 0.0), c(sin(theta), 1.0))

  is_accepted <- accept_ensemble(U_t, 1, 1.0, U_tp1, 1, 1.0, alpha_max, 0, 1e10)
  assert_true(is_accepted)
}

# d changed by 1 (2 -> 1): the angle criterion is vacuously satisfied (no common dimension to
# compare), and d_max=1 tolerates the change.
test_accept_ensemble_d_mismatch_within_dmax <- function() {
  U_t <- diag(3)
  U_tp1 <- U_t

  is_accepted <- accept_ensemble(U_t, 2, 1.0, U_tp1, 1, 1.0, 0.1, 1, 1e10)
  assert_true(is_accepted)
}

test_accept_ensemble_d_mismatch_exceeds_dmax <- function() {
  U_t <- diag(3)
  U_tp1 <- U_t

  is_accepted <- accept_ensemble(U_t, 2, 1.0, U_tp1, 1, 1.0, 0.1, 0, 1e10)
  assert_true(!is_accepted)
}

# G changes by a factor of 10 (ln(10) ~ 2.303): rejected at G_max=1.0.
test_accept_ensemble_g_ratio_exceeds_max <- function() {
  U_t <- diag(3)
  U_tp1 <- U_t

  is_accepted <- accept_ensemble(U_t, 2, 1.0, U_tp1, 2, 10.0, 0.1, 0, 1.0)
  assert_true(!is_accepted)
}

test_accept_ensemble_nonpositive_g <- function() {
  U_t <- diag(3)
  U_tp1 <- U_t

  assert_error(accept_ensemble(U_t, 2, 0.0, U_tp1, 2, 1.0, 0.1, 0, 1e10),
               "Expected error for G_t <= 0", ERR_INVALID_INPUT)
}

test_accept_ensemble_d_out_of_range <- function() {
  U_t <- diag(3)
  U_tp1 <- U_t

  assert_error(accept_ensemble(U_t, 4, 1.0, U_tp1, 2, 1.0, 0.1, 0, 1e10),
               "Expected error for d_t > n_dimensions", ERR_INVALID_INPUT)
}

run_all_tests()
