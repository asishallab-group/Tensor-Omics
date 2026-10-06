source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12

test_gjct_permutation_test_basic <- function() {
  # Both p-value vectors must have one entry per study and land in (0, 1] for a simple 2-study,
  # 2-point scenario (the Laplace correction keeps them above 0).
  n_bins <- 3; n_points <- 2; n_studies <- 2; n_permutations <- 20

  mean_pmf_counts <- matrix(0L, n_bins, n_points)
  mean_pmf_counts[, 1] <- c(5, 3, 2)
  mean_pmf_counts[, 2] <- c(2, 2, 6)
  mean_pmf_included_n_reps <- c(10L, 10L)
  included_n_reps <- matrix(0L, n_points, n_studies)
  included_n_reps[, 1] <- c(4, 4)
  included_n_reps[, 2] <- c(6, 6)
  mean_pmf <- mean_pmf_counts / 10.0
  global_jsd_observed <- c(0.3, 0.3)

  result <- gjct_permutation_test(n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
                                  included_n_reps, global_jsd_observed, random_seed = 123L)

  for (key in c("p_values", "p_values_observed_consensus")) {
    p_values <- result[[key]]
    assert_true(length(p_values) == n_studies)
    assert_true(all(p_values > 0.0 & p_values <= 1.0))
  }
}

# Same random_seed twice must give identical p_values and p_values_observed_consensus -- mirrors
# the Fortran test of the same name (test/mod_test_data_integration_js_comp_test.F90).
test_gjct_permutation_test_seeded_reproducibility <- function() {
  n_bins <- 3; n_points <- 2; n_studies <- 2; n_permutations <- 20

  mean_pmf_counts <- matrix(0L, n_bins, n_points)
  mean_pmf_counts[, 1] <- c(5, 3, 2)
  mean_pmf_counts[, 2] <- c(2, 2, 6)
  mean_pmf_included_n_reps <- c(10L, 10L)
  included_n_reps <- matrix(0L, n_points, n_studies)
  included_n_reps[, 1] <- c(4, 4)
  included_n_reps[, 2] <- c(6, 6)
  mean_pmf <- mean_pmf_counts / 10.0
  global_jsd_observed <- c(0.3, 0.3)

  p_first <- gjct_permutation_test(n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
                                   included_n_reps, global_jsd_observed, random_seed = 123L)
  p_second <- gjct_permutation_test(n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
                                    included_n_reps, global_jsd_observed, random_seed = 123L)

  for (key in c("p_values", "p_values_observed_consensus")) {
    assert_equal_numeric(p_first[[key]], p_second[[key]], TOL, paste("same random_seed twice must give identical", key))
  }
}

# Per-study draws that do not fit the pooled counts (pool c(2, 2), each of two studies drawing 4)
# fail with invalid input, not tied to any one argument -- mirrors the Fortran test of the same
# name (test/mod_test_data_integration_js_comp_test.F90).
test_gjct_permutation_test_rejects_draws_larger_than_pool <- function() {
  mean_pmf_counts <- matrix(c(2L, 2L), nrow = 2)
  mean_pmf <- matrix(c(0.5, 0.5), nrow = 2)
  mean_pmf_included_n_reps <- c(8L)
  included_n_reps <- matrix(c(4L, 4L), nrow = 1)
  global_jsd_observed <- c(0.1, 0.1)

  call <- function() {
    gjct_permutation_test(3L, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
                          included_n_reps, global_jsd_observed, random_seed = 7L)
  }
  assert_error(call(), "expected ERR_INVALID_INPUT for draws larger than the pooled counts", ERR_INVALID_INPUT)
  err <- tryCatch({ call(); NULL }, error = function(e) e)
  assert_true(is.na(err$argument), "expected no argument named")
}

run_all_tests()
