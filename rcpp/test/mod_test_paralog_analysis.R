# Comprehensive R test suite for tox_paralog_analysis (mirrors python/test/mod_test_tox_paralog_analysis.py)

source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

normalize_unit <- function(v) {
  n <- sqrt(sum(v * v))
  if (n == 0) return(v)
  v / n
}

test_paralog_functions <- function() {

  # Testing mask logic
  n_paralogs <- 5L
  i_paralog <- 2L

  chunk_count <- tox_mask_chunk_count(n_paralogs)
  assert_true(chunk_count == 1L)

  bit_mask <- integer(chunk_count)
  bit_mask[(i_paralog %/% 32L) + 1L] <- bitwShiftL(1L, i_paralog %% 32L)
  state <- tox_mask_check_state(bit_mask, i_paralog + 1L)
  assert_true(isTRUE(state))

  # Testing pattern filtering
  n_families <- 1L
  gene_to_fam <- rep(1L, n_paralogs)
  angles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  threshold <- 0.6

  dosage_mask <- tox_filter_paralogs_by_pattern_dosage_effect(angles, threshold, gene_to_fam, n_families)
  assert_true(all(dim(dosage_mask) == c(1L, 1L)))
  assert_true(dosage_mask[1, 1] == 7L)

  subfunc_mask <- tox_filter_paralogs_by_pattern_subfunctionalization(angles, threshold, gene_to_fam, n_families)
  assert_true(all(dim(subfunc_mask) == c(1L, 1L)))
  assert_true(subfunc_mask[1, 1] == 24L)

  # Testing work-array size calculation
  max_subset_size <- 3L
  work_info <- tox_calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, as.integer(dosage_mask[, 1]))
  assert_true(work_info$work_array_size == 3L)
  assert_true(work_info$actual_max_subset_size == 3L)

  # Testing dosage-effect detection
  ancestor <- c(1.0, 1.0)
  paralogs <- matrix(c(
    1.1, 1.2, 1.3, 1.4, 1.5,
    0.9, 0.8, 0.7, 0.6, 0.5
  ), nrow = 2, byrow = TRUE)

  dosage_result <- tox_detect_dosage_effect(
    ancestor = ancestor,
    genes = paralogs,
    filtered_paralogs_mask = as.integer(dosage_mask[, 1]),
    max_subset_size = n_paralogs,
    gain_gamma = 0.1,
    max_angle = pi
  )

  assert_true(dosage_result$n_results == 3L)
  assert_true(all(dosage_result$results[1, ] == c(6L, 3L, 5L)))

  # Testing subfunctionalization detection
  norms <- sqrt(colSums(paralogs ^ 2))
  sorted_perm <- as.integer(order(norms))

  subfunc_result <- tox_detect_subfunctionalization(
    ancestor = ancestor,
    genes = paralogs,
    rdi_threshold = 0.5,
    filtered_paralogs_mask = as.integer(subfunc_mask[, 1]),
    max_subset_size = n_paralogs,
    paralog_norms = norms,
    sorted_paralog_norms_perm = sorted_perm
  )

  assert_true(subfunc_result$n_results == 0L)
  assert_true(ncol(subfunc_result$results) == 0L)

  # Edge cases
  error_caught <- FALSE
  tryCatch({
    tox_mask_chunk_count(0L)
  }, error = function(e) {
    error_caught <<- TRUE
  })
  assert_true(error_caught)

  single_mask <- tox_mask_chunk_count(1L)
  assert_true(single_mask == 1L)

}


test_detect_neofunctionalization <- function() {

  # Case 1: differences below threshold -> all FALSE
  ancestors <- matrix(c(5, 3,
                        2, 1), nrow = 2, ncol = 2)
  ancestors <- apply(ancestors, 2, normalize_unit)

  gene_to_fam <- as.integer(c(1, 2, 1))
  thresholds <- c(0.05, 0.05)

  genes <- matrix(0.0, nrow = 2, ncol = 3)
  for (i in seq_len(3)) {
    genes[, i] <- ancestors[, gene_to_fam[i]]
  }

  neofunc <- tox_detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
  expected <- matrix(FALSE, nrow = 3, ncol = 2)
  assert_true(identical(neofunc, expected))

  # Case 2: differences above threshold -> some TRUE
  thresholds <- c(0.2, 0.2)
  genes2 <- matrix(0.0, nrow = 2, ncol = 3)
  for (i in seq_len(3)) {
    genes2[, i] <- ancestors[, gene_to_fam[i]] - thresholds * gene_to_fam[i]
  }

  neofunc2 <- tox_detect_neofunctionalization(ancestors, genes2, gene_to_fam, thresholds)
  expected2 <- matrix(c(FALSE, TRUE, FALSE,
                        FALSE, TRUE, FALSE), nrow = 3, ncol = 2)
  assert_true(identical(neofunc2, expected2))

}

run_all_tests()
