# Test script for TensorOmics normalization pipeline with mocked data
# This script generates a small synthetic dataset and runs the full pipeline,
# including fold change calculation for multiple conditions.

source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

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
  normalized_matrix_std <- tox_normalize_by_std_dev(mock_matrix, span = 0.75, degree = 2)
  normalized_matrix_qtl <- tox_quantile_normalization(normalized_matrix_std)
  normalized_matrix_log <- tox_log2_transformation(normalized_matrix_qtl)
  dimnames(normalized_matrix_log) <- list(row_names, gene_ids)
  averaged_df <- tox_calculate_tissue_averages(normalized_matrix_log)
  # `tox_calculate_tissue_averages` returns a data.frame; downstream
  # functions expect a numeric matrix. Convert here to avoid validation errors.
  averaged_mat <- as.matrix(averaged_df)
  # --- Calculate fold changes for multiple diets ---
  fc_df <- tox_calculate_fold_changes(
    expr = averaged_mat,
    control_pattern = "dietM",
    condition_patterns = c("dietP", "dietQ")
  )
  # --- Print results ---
  invisible(force(mock_matrix))
  invisible(force(normalized_matrix_std))
  invisible(force(normalized_matrix_qtl))
  invisible(force(normalized_matrix_log))
  invisible(force(averaged_df))
  invisible(force(fc_df))
}

run_all_tests()
