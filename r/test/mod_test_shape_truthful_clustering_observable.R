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

run_all_tests()
