# Test script for TensorOmics normalization pipeline with mocked data
# This script generates a small synthetic dataset and runs the full pipeline,
# including fold change calculation for multiple conditions.

source("r/load_tensor_omics.R")
source("r/test_helpers.R")

test_basic_calling <- function() {
  # --- Mocked gene expression data ---
  set.seed(42)
  gene_ids <- paste0("gene", 1:10)
  tissues <- c("brain", "muscle")
  diets <- c("dietM", "dietP", "dietQ")
  reps <- 1:2

  # Build column names: e.g., brain_dietM_1, brain_dietM_2, ...
  row_names <- as.vector(sapply(tissues, function(tis) sapply(diets, function(diet) paste0(tis, "_", diet, "_", reps))))
  row_names <- as.vector(row_names)

  # Create a 10 genes x (2 tissues x 3 diets x 2 reps = 12) matrix with stable variance
  mock_matrix <- matrix(0, nrow=12, ncol=10)
  for (j in 1:nrow(mock_matrix)) {
    for (i in 1:ncol(mock_matrix)) {
      mock_matrix[j, i] <- i * 2 + j * 0.5
    }
  }
  dimnames(mock_matrix) <- list(row_names, gene_ids)
  force(dimnames)

  # --- Run normalization pipeline ---
  normalized_matrix_std <- normalize_by_std_dev(mock_matrix, span = 0.75, degree = 2)
  # quantile_normalization returns the rank means alongside the matrix
  normalized_matrix_qtl <- quantile_normalization(normalized_matrix_std)$normalized_expr
  normalized_matrix_log <- log2_transformation(normalized_matrix_qtl)
  dimnames(normalized_matrix_log) <- list(row_names, gene_ids)

  # 2 tissues x 3 diets, 2 replicates each: the routine takes the replicate counts
  # rather than parsing them back out of the row names
  reps_per_tissue <- rep(2L, 6L)
  averaged_mat <- calc_tiss_avg(reps_per_tissue, normalized_matrix_log)

  # --- Fold changes, control paired with condition one-to-one ---
  # averaged tissues are brain_dietM/P/Q then muscle_dietM/P/Q, so dietM is 1 and 4
  fc_df <- calc_fchange(
    control_tissues = c(1L, 4L),
    condition_tissues = c(2L, 5L),
    expr = averaged_mat
  )
  # --- Print results ---
  invisible(force(mock_matrix))
  invisible(force(normalized_matrix_std))
  invisible(force(normalized_matrix_qtl))
  invisible(force(normalized_matrix_log))
  invisible(force(averaged_mat))
  invisible(force(fc_df))
}

run_all_tests()
