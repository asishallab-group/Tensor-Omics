source("r/load_tensor_omics.R")
source("r/test_helpers.R")

TOL <- 1e-12

test_gjct_permutation_test_basic <- function() {
  # p_values must land in [0, 1] for a simple 2-study, 2-point scenario.
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

  p_values <- gjct_permutation_test(n_permutations, mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps,
                                    included_n_reps, global_jsd_observed, random_seed = 123L)

  assert_true(length(p_values) == n_studies)
  assert_true(all(p_values >= 0.0 & p_values <= 1.0))
}

# Same random_seed twice must give identical p_values -- mirrors the Fortran test of the same name
# (test/mod_test_data_integration_js_comp_test.F90).
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

  assert_equal_numeric(p_first, p_second, TOL, "same random_seed twice must give identical p_values")
}

run_all_tests()
