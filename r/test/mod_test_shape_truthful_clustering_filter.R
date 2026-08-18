source("r/load_tensor_omics.R")
source("r/test_helpers.R")

STOP_REASON_MAX_SIZE <- 1L
STOP_REASON_REJECTED_AFTER_STABLE <- 2L
STOP_REASON_REJECTED_IMMEDIATELY <- 3L
STOP_REASON_FIXED_POINT <- 4L

# =====================
# filter_ensembles_by_stop_condition
# =====================
test_filter_stop_condition_exact_exclusion <- function() {
  stop_reason <- c(STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE,
                   STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT)
  allowed <- c(TRUE, TRUE, FALSE, TRUE)

  eligible <- filter_ensembles_by_stop_condition(stop_reason, allowed)

  assert_true(eligible[1])
  assert_true(eligible[2])
  assert_true(!eligible[3])
  assert_true(eligible[4])
}

test_filter_stop_condition_absent_is_noop <- function() {
  stop_reason <- c(STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE,
                   STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT)

  eligible <- filter_ensembles_by_stop_condition(stop_reason)

  assert_true(all(eligible))
}

test_filter_stop_condition_all_four_values <- function() {
  for (r in 1:4) {
    stop_reason <- c(r)
    allowed <- rep(TRUE, 4)
    allowed[r] <- FALSE

    eligible <- filter_ensembles_by_stop_condition(stop_reason, allowed)

    assert_true(!eligible[1])
  }
}

test_filter_stop_condition_invalid_value <- function() {
  stop_reason <- c(5L)
  assert_error(filter_ensembles_by_stop_condition(stop_reason),
               "an out-of-range ensemble_stop_reason value must be rejected", ERR_INVALID_INPUT)
}

test_filter_stop_condition_zero_ensembles <- function() {
  stop_reason <- integer(0)
  eligible <- filter_ensembles_by_stop_condition(stop_reason)
  assert_true(length(eligible) == 0)
}

# =====================
# filter_ensembles_by_dimension
# =====================
test_filter_dimension_d_min_only <- function() {
  d_final <- c(0L, 1L, 2L, 3L)
  has_final <- rep(TRUE, 4)

  eligible <- filter_ensembles_by_dimension(3L, d_final, has_final, filter_dim_min = 2L)

  assert_true(!eligible[1])
  assert_true(!eligible[2])
  assert_true(eligible[3])
  assert_true(eligible[4])
}

test_filter_dimension_d_max_only <- function() {
  d_final <- c(0L, 1L, 2L, 3L)
  has_final <- rep(TRUE, 4)

  eligible <- filter_ensembles_by_dimension(3L, d_final, has_final, filter_dim_max = 1L)

  assert_true(eligible[1])
  assert_true(eligible[2])
  assert_true(!eligible[3])
  assert_true(!eligible[4])
}

test_filter_dimension_both_bounds <- function() {
  d_final <- c(0L, 1L, 2L, 3L)
  has_final <- rep(TRUE, 4)

  eligible <- filter_ensembles_by_dimension(3L, d_final, has_final, filter_dim_min = 1L, filter_dim_max = 2L)

  assert_true(!eligible[1])
  assert_true(eligible[2])
  assert_true(eligible[3])
  assert_true(!eligible[4])
}

test_filter_dimension_both_absent_is_noop <- function() {
  d_final <- c(0L, 3L)
  has_final <- c(TRUE, FALSE)

  eligible <- filter_ensembles_by_dimension(3L, d_final, has_final)

  assert_true(all(eligible))
}

test_filter_dimension_no_final_excluded_once_bound_present <- function() {
  d_final <- c(1L)
  has_final <- c(FALSE)

  eligible <- filter_ensembles_by_dimension(3L, d_final, has_final, filter_dim_min = 0L, filter_dim_max = 3L)

  assert_true(!eligible[1])
}

test_filter_dimension_invalid_n_dimensions <- function() {
  d_final <- c(0L)
  has_final <- c(TRUE)
  assert_error(filter_ensembles_by_dimension(1L, d_final, has_final, filter_dim_min = 0L),
               "n_dimensions=1 must be rejected (minimum is 2)", ERR_INVALID_INPUT)
}

test_filter_dimension_d_min_exceeds_d_max_still_computes <- function() {
  # filter_dim_min > filter_dim_max is not itself validated -- simply unsatisfiable, every ensemble ineligible.
  d_final <- c(0L, 1L, 2L)
  has_final <- rep(TRUE, 3)

  eligible <- filter_ensembles_by_dimension(2L, d_final, has_final, filter_dim_min = 2L, filter_dim_max = 1L)

  assert_true(!any(eligible))
}

# =====================
# filter_ensembles_by_var_explained
# =====================
test_filter_var_explained_clean_fixture <- function() {
  # S=[10,1], k=2, d=1 -> eigenvalues [100,1] -> ve = 100/101 ~ 0.9901.
  S_final <- matrix(c(10.0, 1.0), nrow = 2, ncol = 1)
  d_final <- c(1L)
  k_final <- c(2L)
  has_final <- c(TRUE)

  eligible <- filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min = 0.9)

  assert_true(eligible[1])
}

test_filter_var_explained_threshold_at_boundary <- function() {
  S_final <- matrix(c(10.0, 1.0), nrow = 2, ncol = 1)
  d_final <- c(1L)
  k_final <- c(2L)
  has_final <- c(TRUE)
  ve <- 100.0 / 101.0

  assert_true(filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final,
                                                var_explained_min = ve - 1e-9)[1])
  assert_true(filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min = ve)[1])
  assert_true(!filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final,
                                                 var_explained_min = ve + 1e-9)[1])
}

test_filter_var_explained_k_le_one_guard <- function() {
  S_final <- matrix(c(10.0, 1.0, 10.0, 1.0), nrow = 2, ncol = 2)
  d_final <- c(1L, 1L)
  k_final <- c(1L, 0L)
  has_final <- c(TRUE, TRUE)

  eligible <- filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min = 0.0)

  assert_true(!eligible[1])
  assert_true(!eligible[2])
}

test_filter_var_explained_absent_is_noop <- function() {
  S_final <- matrix(c(10.0, 1.0), nrow = 2, ncol = 1)
  d_final <- c(1L)
  k_final <- c(1L) # degenerate k, irrelevant when absent
  has_final <- c(TRUE)

  eligible <- filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final)

  assert_true(eligible[1])
}

test_filter_var_explained_invalid_threshold <- function() {
  S_final <- matrix(c(10.0, 1.0), nrow = 2, ncol = 1)
  d_final <- c(1L)
  k_final <- c(2L)
  has_final <- c(TRUE)

  assert_error(filter_ensembles_by_var_explained(S_final, d_final, k_final, has_final, var_explained_min = 1.5),
               "var_explained_min > 1.0 must be rejected", ERR_INVALID_INPUT)
}

# =====================
# filter_ensembles (combined orchestrator)
# =====================
test_filter_ensembles_combined_different_criteria <- function() {
  # D=2, o=1, 4 ensembles, each failing a different single criterion (or none).
  U <- array(0.0, dim = c(2, 2, 1, 4))
  mu <- array(0.0, dim = c(2, 1, 4))
  G <- array(0.0, dim = c(1, 4))
  k <- array(2L, dim = c(1, 4))
  accepted <- array(TRUE, dim = c(1, 4))
  d <- array(1L, dim = c(1, 4))
  d[1, 2] <- 2L # ensemble 2 fails filter_dim_max=1
  S <- array(0.0, dim = c(2, 1, 4))
  S[, 1, 1] <- c(10.0, 1.0)
  S[, 1, 2] <- c(10.0, 1.0)
  S[, 1, 3] <- c(1.0, 10.0) # ensemble 3 fails var_explained_min=0.5
  S[, 1, 4] <- c(10.0, 1.0)

  stop_reason <- rep(STOP_REASON_FIXED_POINT, 4)
  stop_reason[1] <- STOP_REASON_REJECTED_IMMEDIATELY # ensemble 1 fails stop condition
  allowed <- c(TRUE, TRUE, FALSE, TRUE)

  result <- filter_ensembles(U, d, S, mu, G, k, accepted, stop_reason,
                             allowed_stop_reasons = allowed, filter_dim_max = 1L, var_explained_min = 0.5)

  assert_true(!result$eligible_by_stop_condition[1])
  assert_true(result$eligible_by_dimension[1])
  assert_true(result$eligible_by_var_explained[1])
  assert_true(!result$eligible[1])

  assert_true(result$eligible_by_stop_condition[2])
  assert_true(!result$eligible_by_dimension[2])
  assert_true(result$eligible_by_var_explained[2])
  assert_true(!result$eligible[2])

  assert_true(result$eligible_by_stop_condition[3])
  assert_true(result$eligible_by_dimension[3])
  assert_true(!result$eligible_by_var_explained[3])
  assert_true(!result$eligible[3])

  assert_true(result$eligible_by_stop_condition[4])
  assert_true(result$eligible_by_dimension[4])
  assert_true(result$eligible_by_var_explained[4])
  assert_true(result$eligible[4])
}

test_filter_ensembles_all_omitted_is_noop <- function() {
  U <- array(0.0, dim = c(2, 2, 1, 2))
  d <- array(0L, dim = c(1, 2))
  S <- array(0.0, dim = c(2, 1, 2))
  mu <- array(0.0, dim = c(2, 1, 2))
  G <- array(0.0, dim = c(1, 2))
  k <- array(0L, dim = c(1, 2))
  accepted <- array(FALSE, dim = c(1, 2))
  stop_reason <- rep(STOP_REASON_FIXED_POINT, 2)

  result <- filter_ensembles(U, d, S, mu, G, k, accepted, stop_reason)

  assert_true(all(result$eligible))
  assert_true(all(result$eligible_by_stop_condition))
  assert_true(all(result$eligible_by_dimension))
  assert_true(all(result$eligible_by_var_explained))
}

run_all_tests()
