# Test script for TensorOmics normalization pipeline with mocked data
# This script generates a small synthetic dataset and runs the full pipeline,
# including fold change calculation for multiple conditions.

source("r/tensoromics_functions.R")

# --- Mocked gene expression data ---
set.seed(42)
gene_ids <- paste0("gene", 1:5)
tissues <- c("brain", "muscle")
diets <- c("dietM", "dietP", "dietQ")
reps <- 1:2

# Build column names: e.g., brain_dietM_1, brain_dietM_2, ...
col_names <- as.vector(sapply(tissues, function(tis) sapply(diets, function(diet) paste0(tis, "_", diet, "_", reps))))
col_names <- as.vector(col_names)

# Create a 5 genes x (2 tissues x 3 diets x 2 reps = 12) matrix with random values
mock_matrix <- matrix(runif(5 * 12, min=10, max=100), nrow=5, ncol=12)
dimnames(mock_matrix) <- list(gene_ids, col_names)

# --- Run normalization pipeline ---
normalized_matrix_std <- tox_normalize_by_std_dev(mock_matrix)
normalized_matrix_qtl <- tox_quantile_normalization(normalized_matrix_std)
normalized_matrix_log <- tox_log2_transformation(normalized_matrix_qtl)
averaged_df <- tox_calculate_tissue_averages(normalized_matrix_log)

# --- Calculate fold changes for multiple diets ---
fc_df <- tox_calculate_fc_by_patterns(
  df = averaged_df,
  control_pattern = "dietM",
  condition_patterns = c("dietP", "dietQ")
)

# --- Print results ---
cat("\nMocked input matrix:\n")
print(mock_matrix)
cat("\nAveraged by tissue group:\n")
print(averaged_df)
cat("\nFold changes (dietP, dietQ vs dietM):\n")
print(fc_df)
