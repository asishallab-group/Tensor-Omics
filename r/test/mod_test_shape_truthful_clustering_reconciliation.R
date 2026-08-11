source("r/load_tensor_omics.R")
source("r/test_helpers.R")

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

test_report_mode_no_transitive_merge <- function() {
  res <- ensemble_reconciliation(fixture(), mode = "report", report_overlap_coefficient = TRUE, max_group_size = 2)

  assert_true(res$n_super_ensembles == 3)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2)))
  assert_true(all(res$super_ensembles[, 2] == c(2, 3)))
  assert_true(all(res$super_ensembles[, 3] == c(5, 6)))
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 1] - 0.5) < 1e-9)
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 2] - 0.5) < 1e-9)
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 3] - 2 / 3) < 1e-9)
}

test_merge_any_transitive <- function() {
  res <- ensemble_reconciliation(fixture(), mode = "merge_any", report_overlap_coefficient = TRUE, max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 1] - 0.5) < 1e-9)
  assert_true(abs(res$super_ensembles_overlap_coefficient[2, 1] - 0.5) < 1e-9)
  assert_true(all(res$super_ensembles[, 2] == c(5, 6, 0)))
  assert_true(abs(res$super_ensembles_overlap_coefficient[1, 2] - 2 / 3) < 1e-9)
}

test_merge_overlap_coefficient_threshold_excludes_weak_chain <- function() {
  # min_overlap_coefficient=0.6 excludes the 1-2-3 chain (OC = 0.5) but not the 5-6 pair (OC = 2/3).
  res <- ensemble_reconciliation(fixture(), mode = "merge_overlap_coefficient", min_overlap_coefficient = 0.6, max_group_size = 3)

  assert_true(res$n_super_ensembles == 1)
  assert_true(all(res$super_ensembles[, 1] == c(5, 6, 0)))
}

test_merge_overlap_coefficient_threshold_includes_all <- function() {
  res <- ensemble_reconciliation(fixture(), mode = "merge_overlap_coefficient", min_overlap_coefficient = 0.4, max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles[, 1] == c(1, 2, 3)))
}

test_overlap_coefficient_not_computed_unless_requested <- function() {
  res <- ensemble_reconciliation(fixture(), mode = "merge_any", max_group_size = 3)

  assert_true(res$n_super_ensembles == 2)
  assert_true(all(res$super_ensembles_overlap_coefficient == 0))
}

test_group_exceeds_max_group_size <- function() {
  assert_error(ensemble_reconciliation(fixture(), mode = "merge_any", max_group_size = 2),
               "a 3-member group must not fit in max_group_size=2", ERR_SIZE_MISMATCH)
}

test_n_ensembles_too_small <- function() {
  m <- matrix(FALSE, nrow = 14, ncol = 1)
  assert_error(ensemble_reconciliation(m, max_group_size = 2),
               "Expected error for n_ensembles=1", ERR_INVALID_INPUT)
}

run_all_tests()
