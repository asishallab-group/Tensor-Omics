source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12

# =====================
# normal_error
# =====================
test_normal_error_basic <- function() {
  # D=3, d=1 -- normal_error is the sum of the two smallest (normal-space) eigenvalues.
  eigenvalues <- c(9.0, 4.0, 1.0)
  value <- normal_error(1, eigenvalues)
  assert_true(abs(value - 5.0) < TOL)
}

test_normal_error_zero_tangent_dims <- function() {
  # d=0: no tangent directions, everything is "normal" -- the full sum.
  eigenvalues <- c(9.0, 4.0, 1.0)
  value <- normal_error(0, eigenvalues)
  assert_true(abs(value - 14.0) < TOL)
}

test_normal_error_all_tangent_dims <- function() {
  # d=D: every direction is tangent, nothing left over -- sum over the empty range is zero.
  eigenvalues <- c(9.0, 4.0, 1.0)
  value <- normal_error(3, eigenvalues)
  assert_true(abs(value - 0.0) < TOL)
}

test_normal_error_d_out_of_range <- function() {
  eigenvalues <- c(9.0, 4.0, 1.0)
  assert_error(normal_error(4, eigenvalues), "Expected error for d > n_dimensions", ERR_INVALID_INPUT)
}

test_normal_error_negative_eigenvalue <- function() {
  eigenvalues <- c(9.0, -1.0, 1.0)
  assert_error(normal_error(1, eigenvalues), "Expected error for a negative eigenvalue", ERR_INVALID_INPUT)
}

test_normal_error_zero_dimensions <- function() {
  eigenvalues <- numeric(0)
  assert_error(normal_error(0, eigenvalues), "Expected error for n_dimensions=0", ERR_EMPTY_INPUT)
}

# =====================
# tangent_scales
# =====================
test_tangent_scales_basic <- function() {
  # D=3, d=2 -- tangent_scales is the square root of the two largest eigenvalues.
  eigenvalues <- c(9.0, 4.0, 1.0)
  value <- tangent_scales(2, eigenvalues)
  expected <- c(3.0, 2.0)
  assert_true(all(abs(value - expected) < TOL))
}

test_tangent_scales_all_dims <- function() {
  eigenvalues <- c(9.0, 4.0, 1.0)
  value <- tangent_scales(3, eigenvalues)
  expected <- c(3.0, 2.0, 1.0)
  assert_true(all(abs(value - expected) < TOL))
}

test_tangent_scales_zero_dims <- function() {
  # d=0: no tangent directions -- must return a well-defined, empty array.
  eigenvalues <- c(9.0, 4.0, 1.0)
  value <- tangent_scales(0, eigenvalues)
  assert_true(length(value) == 0)
}

test_tangent_scales_d_out_of_range <- function() {
  eigenvalues <- c(9.0, 4.0, 1.0)
  assert_error(tangent_scales(4, eigenvalues), "Expected error for d > n_dimensions", ERR_INVALID_INPUT)
}

test_tangent_scales_negative_eigenvalue <- function() {
  eigenvalues <- c(9.0, -4.0, 1.0)
  assert_error(tangent_scales(2, eigenvalues), "Expected error for a negative eigenvalue", ERR_INVALID_INPUT)
}

# =====================
# observable
# =====================
test_observable_full_rank_rectangle <- function() {
  # A rectangle in the z=0 plane, embedded in 3D: full economy-mode rank (rank=min(D,k)=3=D,
  # no zero-padding), with distinct, hand-computable eigenvalues.
  vectors <- rbind(
    c(0.0, 2.0, 0.0, 2.0),
    c(0.0, 0.0, 1.0, 1.0),
    c(0.0, 0.0, 0.0, 0.0)
  )
  mask <- rep(TRUE, 4)

  result <- observable(vectors, mask)
  assert_true(all(abs(result$mu - c(1.0, 0.5, 0.0)) < 1e-6))
  assert_true(result$d == 2)
  assert_true(abs(result$eigenvalues[1] - 4.0 / 3.0) < 1e-6)
  assert_true(abs(result$eigenvalues[2] - 1.0 / 3.0) < 1e-6)
  assert_true(abs(result$eigenvalues[3] - 0.0) < 1e-6)
  assert_true(abs(result$normal_error_value - 0.0) < 1e-6)
  assert_true(abs(result$tangent_scales_value[1] - sqrt(4.0 / 3.0)) < 1e-6)
  assert_true(abs(result$tangent_scales_value[2] - sqrt(1.0 / 3.0)) < 1e-6)
  assert_true(result$G > 1e10)
  # U columns are the standard basis up to sign (diagonal covariance, distinct eigenvalues).
  assert_true(all(abs(abs(result$U[, 1]) - c(1.0, 0.0, 0.0)) < 1e-6))
  assert_true(all(abs(abs(result$U[, 2]) - c(0.0, 1.0, 0.0)) < 1e-6))
  assert_true(all(abs(abs(result$U[, 3]) - c(0.0, 0.0, 1.0)) < 1e-6))
}

test_observable_low_rank_padding <- function() {
  # Three collinear points (intrinsic rank 1) embedded in a 5D ambient space: economy-mode
  # rank = min(D,k) = 3 < D = 5, so U columns 4-5 and eigenvalues 4-5 must be observable's
  # own zero-padding, not LAPACK output.
  vectors <- matrix(0.0, nrow = 5, ncol = 3)
  vectors[1, ] <- c(0.0, 1.0, 2.0)
  mask <- rep(TRUE, 3)

  result <- observable(vectors, mask)
  assert_true(result$d == 1)
  assert_true(abs(result$eigenvalues[1] - 1.0) < 1e-8)
  assert_true(all(abs(result$eigenvalues[4:5]) < 1e-8))
  assert_true(all(abs(result$U[, 4]) < 1e-8))
  assert_true(all(abs(result$U[, 5]) < 1e-8))
  assert_true(abs(result$tangent_scales_value[1] - 1.0) < 1e-8)
  assert_true(all(abs(result$tangent_scales_value[2:5]) < 1e-8))
  assert_true(result$G > 1e10)
}

test_observable_too_few_members <- function() {
  vectors <- rbind(
    c(0.0, 2.0, 0.0, 2.0),
    c(0.0, 0.0, 1.0, 1.0),
    c(0.0, 0.0, 0.0, 0.0)
  )
  mask <- rep(FALSE, 4)
  mask[1] <- TRUE

  assert_error(observable(vectors, mask),
               "Expected error for an ensemble with fewer than 2 members", ERR_INVALID_INPUT)
}

test_observable_dimension_too_small <- function() {
  vectors <- rbind(c(0.0, 1.0, 2.0))
  mask <- rep(TRUE, 3)

  assert_error(observable(vectors, mask), "Expected error for n_dimensions=1", ERR_INVALID_INPUT)
}

# =====================
# ensemble_final_observable
# =====================
test_ensemble_final_observable_trailing_rejected_column <- function() {
  # D=2, o=2, one ensemble: column 1 is the true accepted state (d=1, G=1.0, mu=[0.5,0.0]);
  # column 2 is a rejected candidate with deliberately different values (d=0, G=99.0,
  # mu=[9,9], accepted=FALSE). The extraction must land on column 1, never column 2.
  U <- array(0.0, dim = c(2, 2, 2, 1))
  d <- matrix(c(1L, 0L), nrow = 2, ncol = 1)
  S <- array(0.0, dim = c(2, 2, 1))
  mu <- array(0.0, dim = c(2, 2, 1))
  G <- matrix(c(1.0, 99.0), nrow = 2, ncol = 1)
  k <- matrix(c(2L, 3L), nrow = 2, ncol = 1)
  accepted <- matrix(c(TRUE, FALSE), nrow = 2, ncol = 1)

  S[, 1, 1] <- c(0.5, 2.0)
  S[, 2, 1] <- c(7.0, 7.0)
  mu[, 1, 1] <- c(0.5, 0.0)
  mu[, 2, 1] <- c(9.0, 9.0)
  U[, 1, 1, 1] <- c(1.0, 0.0)
  U[, 2, 1, 1] <- c(0.0, 1.0)

  result <- ensemble_final_observable(U, d, S, mu, G, k, accepted)

  assert_true(result$ensemble_has_final[1])
  assert_true(result$ensemble_final_index[1] == 1)
  assert_true(result$ensemble_d_final[1] == 1)
  assert_true(result$ensemble_k_final[1] == 2)
  assert_true(abs(result$ensemble_G_final[1] - 1.0) < TOL)
  assert_true(all(abs(result$ensemble_mu_final[, 1] - c(0.5, 0.0)) < TOL))
  assert_true(all(abs(result$ensemble_S_final[, 1] - c(0.5, 2.0)) < TOL))
}

test_ensemble_final_observable_no_rejection <- function() {
  # Both history columns accepted: extraction must land on the last (most recent) column.
  U <- array(0.0, dim = c(2, 2, 2, 1))
  d <- matrix(c(0L, 1L), nrow = 2, ncol = 1)
  S <- array(0.0, dim = c(2, 2, 1))
  mu <- array(0.0, dim = c(2, 2, 1))
  G <- matrix(c(2.0, 1.5), nrow = 2, ncol = 1)
  k <- matrix(c(2L, 3L), nrow = 2, ncol = 1)
  accepted <- matrix(c(TRUE, TRUE), nrow = 2, ncol = 1)

  S[, 1, 1] <- c(2.0, 0.0)
  S[, 2, 1] <- c(0.5, 3.0)
  mu[, 1, 1] <- c(0.5, 0.0)
  mu[, 2, 1] <- c(1.0, 0.0)

  result <- ensemble_final_observable(U, d, S, mu, G, k, accepted)

  assert_true(result$ensemble_has_final[1])
  assert_true(result$ensemble_final_index[1] == 2)
  assert_true(result$ensemble_d_final[1] == 1)
  assert_true(result$ensemble_k_final[1] == 3)
  assert_true(abs(result$ensemble_G_final[1] - 1.5) < TOL)
  assert_true(all(abs(result$ensemble_mu_final[, 1] - c(1.0, 0.0)) < TOL))
}

test_ensemble_final_observable_has_final_false <- function() {
  # Every history column has k=0 (unpopulated) -- has_final must be FALSE and every _final
  # output must come back exactly zero.
  U <- array(0.0, dim = c(2, 2, 2, 1))
  d <- array(0L, dim = c(2, 1))
  S <- array(0.0, dim = c(2, 2, 1))
  mu <- array(0.0, dim = c(2, 2, 1))
  G <- array(0.0, dim = c(2, 1))
  k <- array(0L, dim = c(2, 1))
  accepted <- array(FALSE, dim = c(2, 1))

  result <- ensemble_final_observable(U, d, S, mu, G, k, accepted)

  assert_true(!result$ensemble_has_final[1])
  assert_true(result$ensemble_final_index[1] == 0)
  assert_true(result$ensemble_d_final[1] == 0)
  assert_true(result$ensemble_k_final[1] == 0)
  assert_true(abs(result$ensemble_G_final[1] - 0.0) < TOL)
  assert_true(all(abs(result$ensemble_mu_final[, 1] - c(0.0, 0.0)) < TOL))
}

test_ensemble_final_observable_multi_ensemble_independence <- function() {
  # Two ensembles with independent histories: ensemble 1 has a trailing rejected column
  # (resolves to its own column 1), ensemble 2 has no rejection (resolves to its column 2).
  U <- array(0.0, dim = c(2, 2, 2, 2))
  d <- matrix(c(1L, 0L, 0L, 1L), nrow = 2, ncol = 2)
  S <- array(0.0, dim = c(2, 2, 2))
  mu <- array(0.0, dim = c(2, 2, 2))
  G <- matrix(c(1.0, 99.0, 3.0, 2.5), nrow = 2, ncol = 2)
  k <- matrix(c(2L, 3L, 4L, 5L), nrow = 2, ncol = 2)
  accepted <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2, ncol = 2)

  mu[, 1, 1] <- c(0.5, 0.0)
  mu[, 2, 1] <- c(9.0, 9.0)
  mu[, 1, 2] <- c(2.0, 2.0)
  mu[, 2, 2] <- c(3.0, 3.0)

  result <- ensemble_final_observable(U, d, S, mu, G, k, accepted)

  assert_true(result$ensemble_has_final[1])
  assert_true(result$ensemble_final_index[1] == 1)
  assert_true(abs(result$ensemble_G_final[1] - 1.0) < TOL)

  assert_true(result$ensemble_has_final[2])
  assert_true(result$ensemble_final_index[2] == 2)
  assert_true(abs(result$ensemble_G_final[2] - 2.5) < TOL)
  assert_true(result$ensemble_k_final[2] == 5)
  assert_true(all(abs(result$ensemble_mu_final[, 2] - c(3.0, 3.0)) < TOL))
}

test_ensemble_final_observable_small_o_evicts_accepted <- function() {
  # o=1: the single available column holds a rejected candidate. No earlier accepted column
  # exists to fall back to -- has_final must be FALSE.
  U <- array(0.0, dim = c(2, 2, 1, 1))
  d <- array(0L, dim = c(1, 1))
  S <- array(0.0, dim = c(2, 1, 1))
  mu <- array(0.0, dim = c(2, 1, 1))
  G <- matrix(99.0, nrow = 1, ncol = 1)
  k <- matrix(3L, nrow = 1, ncol = 1)
  accepted <- matrix(FALSE, nrow = 1, ncol = 1)

  mu[, 1, 1] <- c(9.0, 9.0)

  result <- ensemble_final_observable(U, d, S, mu, G, k, accepted)

  assert_true(!result$ensemble_has_final[1])
  assert_true(result$ensemble_final_index[1] == 0)
  assert_true(abs(result$ensemble_G_final[1] - 0.0) < TOL)
}

run_all_tests()
