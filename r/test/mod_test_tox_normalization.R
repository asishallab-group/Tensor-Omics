# The `tox_normalization` suite: each published procedure can be called from R, returns the
# documented type and shape, and raises the documented error. Whether the numbers are right is the
# Fortran suite's job (test/mod_test_tox_normalization.F90), as the coding guide's "Where to Test
# What" asks.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

N_REPLICATES <- 6L
N_GENES <- 10L
# two tissues of three replicates each
REPS_PER_TISSUE <- c(3L, 3L)

# (replicates, genes), every gene's spread proportional to its level: the mean-sd trend the LOESS
# step fits, positive throughout, so `log2(x + 1)` accepts every normalized value.
.expr <- function() {
  outer(seq_len(N_REPLICATES), seq_len(N_GENES), function(i_replicate, i_gene) i_gene * (10 + 0.5 * i_replicate))
}

.assert_matrix <- function(result, dims, name) {
  assert_true(is.matrix(result), paste0(name, ": expected a matrix"))
  assert_true(is.double(result), paste0(name, ": expected a double matrix"))
  assert_true(identical(dim(result), as.integer(dims)),
              sprintf("%s: expected dim %s, got %s", name, toString(dims), toString(dim(result))))
  assert_true(all(is.finite(result)), paste0(name, ": non-finite values in the result"))
}

test_normalize_unit_length <- function() {
  result <- normalize_unit_length(c(3, 4))
  assert_true(is.double(result) && length(result) == 2L, "expected a double vector of length 2")
  assert_true(all(is.finite(result)), "non-finite values in the result")
}

test_normalize_unit_length_rejects_zero_vector <- function() {
  assert_error(normalize_unit_length(c(0, 0, 0)), "a zero vector has no direction", ERR_DIVISION_BY_ZERO)
}

test_normalization_pipeline <- function() {
  for (use_quantile in c(FALSE, TRUE)) {
    result <- normalization_pipeline(.expr(), REPS_PER_TISSUE, span = 0.75, degree = 2L, use_quantile = use_quantile)
    .assert_matrix(result, c(length(REPS_PER_TISSUE), N_GENES), paste0("use_quantile=", use_quantile))
  }
}

test_normalization_pipeline_rejects_reps_not_summing_to_the_replicates <- function() {
  assert_error(normalization_pipeline(.expr(), c(3L, 2L)),
               "reps_per_tissue sums to 5, but expr has 6 replicates", ERR_SIZE_MISMATCH)
}

test_normalize_by_std_dev <- function() {
  .assert_matrix(normalize_by_std_dev(.expr()), c(N_REPLICATES, N_GENES), "defaults")
  .assert_matrix(normalize_by_std_dev(.expr(), span = 0.75, degree = 1L), c(N_REPLICATES, N_GENES), "span and degree")
}

test_normalize_by_std_dev_rejects_too_few_varying_genes <- function() {
  # the LOESS fit needs five genes that vary across replicates; four is one short
  assert_error(normalize_by_std_dev(.expr()[, 1:4]), "four genes are too few for the fit", ERR_INVALID_INPUT)
}

test_root_mean_sq_normalization <- function() {
  .assert_matrix(root_mean_sq_normalization(.expr()), c(N_REPLICATES, N_GENES), "root_mean_sq_normalization")
}

test_quantile_normalization <- function() {
  result <- quantile_normalization(.expr())
  assert_true(is.list(result), "expected a list")
  assert_true(identical(sort(names(result)), c("normalized_expr", "rank_means")),
              paste("unexpected names", toString(names(result))))
  .assert_matrix(result$normalized_expr, c(N_REPLICATES, N_GENES), "normalized_expr")
  assert_true(is.double(result$rank_means) && length(result$rank_means) == N_GENES,
              "rank_means: expected a double vector, one per gene")
}

test_log2_transformation <- function() {
  tissue_averages <- calc_tiss_avg(REPS_PER_TISSUE, .expr())
  .assert_matrix(log2_transformation(tissue_averages), c(length(REPS_PER_TISSUE), N_GENES), "log2_transformation")
}

test_log2_transformation_rejects_values_at_minus_one <- function() {
  # log2(x + 1) is only defined above -1
  assert_error(log2_transformation(matrix(c(1, -1), nrow = 1)), "log2(0) is undefined", ERR_INVALID_INPUT)
}

test_every_matrix_procedure_rejects_an_empty_expr <- function() {
  empty <- matrix(numeric(0), nrow = N_REPLICATES, ncol = 0)
  procedures <- list(normalize_by_std_dev = normalize_by_std_dev,
                     root_mean_sq_normalization = root_mean_sq_normalization,
                     quantile_normalization = quantile_normalization,
                     log2_transformation = log2_transformation)
  for (name in names(procedures)) {
    assert_error(procedures[[name]](empty), paste0(name, ": no genes"), ERR_EMPTY_INPUT)
  }
}

test_calc_tiss_avg <- function() {
  .assert_matrix(calc_tiss_avg(REPS_PER_TISSUE, .expr()), c(length(REPS_PER_TISSUE), N_GENES), "calc_tiss_avg")
}

test_calc_tiss_avg_rejects_a_tissue_without_replicates <- function() {
  assert_error(calc_tiss_avg(c(3L, 0L, 3L), .expr()), "a tissue with no replicates has no average", ERR_INVALID_INPUT)
}

test_calc_tiss_avg_rejects_reps_not_summing_to_the_replicates <- function() {
  # expr has 6 replicates: 5 would silently shift every gene after the first onto the wrong
  # rows, 7 would read past the end of the matrix
  for (reps_per_tissue in list(c(3L, 2L), c(3L, 4L))) {
    assert_error(calc_tiss_avg(reps_per_tissue, .expr()),
                 paste0("reps_per_tissue ", toString(reps_per_tissue), ", but expr has 6 replicates"), ERR_SIZE_MISMATCH)
  }
}

test_calc_fchange <- function() {
  tissue_averages <- calc_tiss_avg(c(2L, 2L, 2L), .expr())
  # one control against two conditions
  .assert_matrix(calc_fchange(c(1L, 1L), c(2L, 3L), tissue_averages), c(2L, N_GENES), "calc_fchange")
}

test_calc_fchange_rejects_a_tissue_past_the_last <- function() {
  tissue_averages <- calc_tiss_avg(c(2L, 2L, 2L), .expr())
  assert_error(calc_fchange(1L, 4L, tissue_averages), "there are only three tissues", ERR_INVALID_INPUT)
}

test_calc_fchange_rejects_unpaired_tissues <- function() {
  tissue_averages <- calc_tiss_avg(c(2L, 2L, 2L), .expr())
  # the binding's own shape check, before the library is called
  assert_error(calc_fchange(c(1L, 1L), 2L, tissue_averages), "two controls, one condition")
}

run_all_tests()
