source("r/load_tensor_omics.R")
source("r/test_helpers.R")

# ensemble_stop_reason values, 1-indexed -- tox_shape_truthful_clustering_impl's own
# STOP_REASON_MAX_SIZE(1)/STOP_REASON_REJECTED_AFTER_STABLE(2)/STOP_REASON_REJECTED_IMMEDIATELY(3)/
# STOP_REASON_FIXED_POINT(4). Not exported as named constants by any binding (see the Fortran
# kernel's own doc comment on why), so used as plain integers here, matching every other
# consumer of this array.
STOP_REASON_MAX_SIZE <- 1L
STOP_REASON_REJECTED_AFTER_STABLE <- 2L
STOP_REASON_REJECTED_IMMEDIATELY <- 3L
STOP_REASON_FIXED_POINT <- 4L

# N=14 vectors, 6 ensembles. E1={1,2,3,4}, E2={3,4,5,6}, E3={5,6,7,8}: a chain, E1-E2 and
# E2-E3 each intersect at 2 members, both ensembles size 4, so OC = 2/min(4,4) = 0.5. E1-E3 do
# not intersect at all. E4={9,10}: isolated. E5={11,12,13}, E6={12,13,14}: a separate pair,
# intersecting at 2 members, both ensembles size 3, so OC = 2/min(3,3) = 2/3.
fixture <- function() {
  m <- matrix(FALSE, nrow = 14, ncol = 6)
  m[1:4, 1] <- TRUE
  m[3:6, 2] <- TRUE
  m[5:8, 3] <- TRUE
  m[9:10, 4] <- TRUE
  m[11:13, 5] <- TRUE
  m[12:14, 6] <- TRUE
  m
}

# Every ensemble STOP_REASON_FIXED_POINT by default -- a neutral choice that changes nothing
# about any pre-existing test below when no allowed_stop_reasons filter is applied.
stop_reasons <- function() rep(STOP_REASON_FIXED_POINT, 6)

# Minimal, uniform D=2/o=1 history for 6 ensembles -- ensemble_reconciliation's new required
# history arguments, needed only to drive the dimension/variance-explained filters, which none
# of the tests below (other than the dedicated d_min/d_max/var_explained_min ones) exercise.
history_fixture <- function() {
  n_e <- 6
  list(
    ensemble_U_history = array(0.0, dim = c(2, 2, 1, n_e)),
    ensemble_d_history = array(0L, dim = c(1, n_e)),
    ensemble_S_history = array(0.0, dim = c(2, 1, n_e)),
    ensemble_mu_history = array(0.0, dim = c(2, 1, n_e)),
    ensemble_G_history = array(0.0, dim = c(1, n_e)),
    ensemble_k_history = array(2L, dim = c(1, n_e)),
    ensemble_accepted_history = array(TRUE, dim = c(1, n_e))
  )
}

# ensemble_reconciliation, with the history-array block filled in from history_fixture() unless
# the caller supplies its own via hist = list(...).
reconcile <- function(ensemble_masks, ensemble_stop_reason, hist = NULL, ...) {
  h <- if (is.null(hist)) history_fixture() else hist
  do.call(ensemble_reconciliation, c(list(
    ensemble_masks = ensemble_masks, ensemble_stop_reason = ensemble_stop_reason,
    ensemble_U_history = h$ensemble_U_history, ensemble_d_history = h$ensemble_d_history,
    ensemble_S_history = h$ensemble_S_history, ensemble_mu_history = h$ensemble_mu_history,
    ensemble_G_history = h$ensemble_G_history, ensemble_k_history = h$ensemble_k_history,
    ensemble_accepted_history = h$ensemble_accepted_history
  ), list(...)))
}

test_report_mode_no_transitive_merge <- function() {
  res <- reconcile(fixture(), stop_reasons(), mode = "report", report_overlap_coefficient = TRUE,
                   max_group_size = 2)

  assert_true(res$n_super_ensembles == 3)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2)))
  assert_true(all(res$super_ensembles[, 2] == c(2, 3)))
  assert_true(all(res$super_ensembles[, 3] == c(5, 6)))
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 1] - 0.5) < 1e-9)
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 2] - 0.5) < 1e-9)
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 3] - 2 / 3) < 1e-9)
  assert_true(all(res$eligible))
}

test_merge_any_transitive <- function() {
  res <- reconcile(fixture(), stop_reasons(), mode = "merge_any", report_overlap_coefficient = TRUE,
                   max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 1] - 0.5) < 1e-9)
  assert_true(abs(res$super_ensembles_overlap_coefficient[2, 1] - 0.5) < 1e-9)
  assert_true(all(res$super_ensembles[, 2] == c(5, 6, 0)))
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 2] - 2 / 3) < 1e-9)
}

test_merge_overlap_coefficient_threshold_excludes_weak_chain <- function() {
  # min_overlap_coefficient=0.6 excludes the 1-2-3 chain (OC = 0.5) but not the 5-6 pair (OC = 2/3).
  res <- reconcile(fixture(), stop_reasons(), mode = "merge_overlap_coefficient",
                   min_overlap_coefficient = 0.6, max_group_size = 3)

  assert_true(res$n_super_ensembles == 1)
  assert_true(all(res$super_ensembles[, 1] == c(5, 6, 0)))
}

test_merge_overlap_coefficient_threshold_includes_all <- function() {
  res <- reconcile(fixture(), stop_reasons(), mode = "merge_overlap_coefficient",
                   min_overlap_coefficient = 0.4, max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
}

test_overlap_coefficient_not_computed_unless_requested <- function() {
  res <- reconcile(fixture(), stop_reasons(), mode = "merge_any", max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles_overlap_coefficient == 0))
}

test_group_exceeds_max_group_size <- function() {
  assert_error(reconcile(fixture(), stop_reasons(), mode = "merge_any", max_group_size = 2),
               "a 3-member group must not fit in max_group_size=2", ERR_SIZE_MISMATCH)
}

test_n_ensembles_too_small <- function() {
  m <- matrix(FALSE, nrow = 14, ncol = 1)
  sr <- rep(STOP_REASON_FIXED_POINT, 1)
  n_e <- 1
  h <- list(
    ensemble_U_history = array(0.0, dim = c(2, 2, 1, n_e)),
    ensemble_d_history = array(0L, dim = c(1, n_e)),
    ensemble_S_history = array(0.0, dim = c(2, 1, n_e)),
    ensemble_mu_history = array(0.0, dim = c(2, 1, n_e)),
    ensemble_G_history = array(0.0, dim = c(1, n_e)),
    ensemble_k_history = array(2L, dim = c(1, n_e)),
    ensemble_accepted_history = array(TRUE, dim = c(1, n_e))
  )
  assert_error(reconcile(m, sr, hist = h, max_group_size = 2),
               "Expected error for n_ensembles=1", ERR_INVALID_INPUT)
}

test_allowed_stop_reasons_excludes_pair_report_mode <- function() {
  # Ensemble 5 is STOP_REASON_REJECTED_IMMEDIATELY (rest STOP_REASON_FIXED_POINT); excluding
  # STOP_REASON_REJECTED_IMMEDIATELY must drop the (5,6) pair from report mode's output,
  # leaving only the untouched 1-2-3 chain's two pairs.
  sr <- stop_reasons()
  sr[5] <- STOP_REASON_REJECTED_IMMEDIATELY
  allowed <- c(TRUE, TRUE, FALSE, TRUE)

  res <- reconcile(fixture(), sr, mode = "report", allowed_stop_reasons = allowed, max_group_size = 2)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2)))
  assert_true(all(res$super_ensembles[, 2] == c(2, 3)))
  assert_true(!res$eligible[5])
  assert_true(!res$eligible_by_stop_condition[5])
}

test_allowed_stop_reasons_absent_matches_all_true <- function() {
  # Omitting allowed_stop_reasons (NULL), even with a genuinely mixed set of Stop Conditions
  # across ensembles, must behave identically to no filtering at all.
  sr <- stop_reasons()
  sr[1] <- STOP_REASON_MAX_SIZE
  sr[4] <- STOP_REASON_REJECTED_AFTER_STABLE
  sr[6] <- STOP_REASON_REJECTED_IMMEDIATELY

  res <- reconcile(fixture(), sr, mode = "merge_any", report_overlap_coefficient = TRUE,
                   max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
  assert_true(all(res$super_ensembles[, 2] == c(5, 6, 0)))
  assert_true(all(res$eligible))
}

test_allowed_stop_reasons_breaks_transitive_chain <- function() {
  # Ensemble 2 (the sole bridge of the 1-2-3 chain) is STOP_REASON_REJECTED_IMMEDIATELY;
  # excluding it must fully break the chain (1 and 3 do not intersect directly), leaving only
  # the untouched {5,6} group.
  sr <- stop_reasons()
  sr[2] <- STOP_REASON_REJECTED_IMMEDIATELY
  allowed <- c(TRUE, TRUE, FALSE, TRUE)

  res <- reconcile(fixture(), sr, mode = "merge_any", allowed_stop_reasons = allowed, max_group_size = 3)

  assert_true(res$n_super_ensembles == 1)
  assert_true(all(res$super_ensembles[, 1] == c(5, 6, 0)))
}

test_allowed_stop_reasons_noop_when_no_ensemble_matches <- function() {
  # Every ensemble is STOP_REASON_FIXED_POINT; excluding STOP_REASON_REJECTED_AFTER_STABLE
  # (matching none of them) must be a true no-op.
  allowed <- c(TRUE, FALSE, TRUE, TRUE)

  res <- reconcile(fixture(), stop_reasons(), mode = "merge_any", allowed_stop_reasons = allowed,
                   max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
  assert_true(all(res$super_ensembles[, 2] == c(5, 6, 0)))
}

test_dimension_filter_excludes_pair <- function() {
  # Ensemble 5's final intrinsic dimension (2) exceeds d_max=1 (the rest are d=1); excluding it
  # must drop the (5,6) pair from report mode's output.
  h <- history_fixture()
  h$ensemble_d_history[1, ] <- 1L
  h$ensemble_d_history[1, 5] <- 2L

  res <- reconcile(fixture(), stop_reasons(), hist = h, mode = "report", d_max = 1L, max_group_size = 2)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2)))
  assert_true(all(res$super_ensembles[, 2] == c(2, 3)))
  assert_true(!res$eligible_by_dimension[5])
  assert_true(res$eligible_by_dimension[1])
  assert_true(!res$eligible[5])
}

test_var_explained_filter_excludes_pair <- function() {
  # Ensemble 5's final variance explained (eigenvalues [1,100], d=1 -> 1/101 ~ 0.0099) falls
  # far short of var_explained_min=0.5; the rest ([100,1] -> 100/101 ~ 0.99) comfortably clear
  # it. Excluding ensemble 5 must drop the (5,6) pair.
  h <- history_fixture()
  h$ensemble_d_history[1, ] <- 1L
  for (e in 1:6) h$ensemble_S_history[, 1, e] <- c(10.0, 1.0)
  h$ensemble_S_history[, 1, 5] <- c(1.0, 10.0)

  res <- reconcile(fixture(), stop_reasons(), hist = h, mode = "report", var_explained_min = 0.5, max_group_size = 2)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2)))
  assert_true(all(res$super_ensembles[, 2] == c(2, 3)))
  assert_true(!res$eligible_by_var_explained[5])
  assert_true(res$eligible_by_var_explained[1])
}

test_merge_to_super_ensembles_all_eligible_matches_merge_any <- function() {
  # Direct test of the newly-independent merge_to_super_ensembles: an all-TRUE eligible mask
  # must reproduce test_merge_any_transitive's own result, with no history array involved.
  eligible <- rep(TRUE, 6)

  res <- merge_to_super_ensembles(fixture(), eligible, mode = "merge_any", report_overlap_coefficient = TRUE,
                                  max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
  assert_true(all(res$super_ensembles[, 2] == c(5, 6, 0)))
}

test_merge_to_super_ensembles_excludes_ineligible_ensemble <- function() {
  # Ensemble 5 marked ineligible via a hand-constructed eligible mask: the (5,6) edge must
  # never be considered, leaving only the untouched 1-2-3 chain.
  eligible <- rep(TRUE, 6)
  eligible[5] <- FALSE

  res <- merge_to_super_ensembles(fixture(), eligible, mode = "merge_any", report_overlap_coefficient = TRUE,
                                  max_group_size = 3)

  assert_true(res$n_super_ensembles == 1)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
}

run_all_tests()
