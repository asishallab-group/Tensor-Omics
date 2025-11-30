# Comprehensive R test suite for tox_normalization_pipeline (mirrors Fortran unit tests)
source("r/tensoromics_functions.R")

# Helper: compare two numeric vectors/matrices with tolerance
assert_equal <- function(x, y, tol=1e-12, msg="") {
  stopifnot(all(abs(x - y) < tol))
  if (msg != "") cat(msg, "✓\n")
}


basic_test <- function() {
  n_genes <- 2; n_tissues <- 2; n_grps <- 2
  input_matrix <- matrix(c(1,2,3,4), nrow=n_genes, ncol=n_tissues)
  group_s <- c(1L,2L); group_c <- c(1L,1L)
  cat("[DEBUG] basic_test: n_genes=", n_genes, " n_tissues=", n_tissues, " n_grps=", n_grps, "\n")
  cat("[DEBUG] input_matrix=", paste(input_matrix, collapse=","), "\n")
  cat("[DEBUG] group_s=", paste(group_s, collapse=","), " group_c=", paste(group_c, collapse=","), "\n")
  result <- tox_normalization_pipeline(input_matrix, group_s, group_c)
  cat("[DEBUG] buf_log=", paste(result, collapse=","), "\n")
  stopifnot(all(!is.na(result)))
  stopifnot(all(result >= 0))
  stopifnot(dim(result)[1] == n_genes && dim(result)[2] == n_grps)
  cat("basic_test passed\n")
}

edge_case_test <- function() {
  n_genes <- 2; n_tissues <- 2; n_grps <- 2
  input_matrix <- matrix(0, nrow=n_genes, ncol=n_tissues)
  group_s <- c(1L,2L); group_c <- c(1L,1L)
  cat("[DEBUG] edge_case_test: n_genes=", n_genes, " n_tissues=", n_tissues, " n_grps=", n_grps, "\n")
  cat("[DEBUG] input_matrix=", paste(input_matrix, collapse=","), "\n")
  cat("[DEBUG] group_s=", paste(group_s, collapse=","), " group_c=", paste(group_c, collapse=","), "\n")
  result <- tox_normalization_pipeline(input_matrix, group_s, group_c)
  cat("[DEBUG] buf_log=", paste(result, collapse=","), "\n")
  stopifnot(all(result == 0))
  cat("edge_case_test passed\n")
}

pipeline_vs_manual <- function() {
  n_genes <- 2; n_tissues <- 2; n_grps <- 2
  col_names <- c("muscle_dietM_1", "muscle_dietP_1")
  input_matrix <- matrix(c(1,2,3,4), nrow=n_genes, ncol=n_tissues)
  dimnames(input_matrix) <- list(paste0("gene", 1:n_genes), col_names)
  group_s <- c(1L,2L); group_c <- c(1L,1L)
  cat("[DEBUG] pipeline_vs_manual: n_genes=", n_genes, " n_tissues=", n_tissues, " n_grps=", n_grps, "\n")
  cat("[DEBUG] input_matrix=", paste(input_matrix, collapse=","), "\n")
  cat("[DEBUG] colnames=", paste(colnames(input_matrix), collapse=","), "\n")
  cat("[DEBUG] group_s=", paste(group_s, collapse=","), " group_c=", paste(group_c, collapse=","), "\n")
  buf_stddev <- tox_normalize_by_std_dev(input_matrix)
  buf_quant <- tox_quantile_normalization(buf_stddev)
  buf_avg <- tox_calculate_tissue_averages(buf_quant)
  cat("[DEBUG] buf_avg=", paste(as.vector(as.matrix(buf_avg)), collapse=","), " colnames=", paste(colnames(buf_avg), collapse=","), "\n")
  manual_out <- tox_log2_transformation(as.matrix(buf_avg))
  cat("[DEBUG] manual_out=", paste(manual_out, collapse=","), "\n")
  result <- tox_normalization_pipeline(input_matrix, group_s, group_c)
  cat("[DEBUG] buf_log=", paste(result, collapse=","), "\n")
  assert_equal(as.vector(result), as.vector(manual_out), msg="pipeline_vs_manual")
}

# Run all tests
cat("=================================================\n")
cat("    NORMALIZATION PIPELINE FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

basic_test()
edge_case_test()
pipeline_vs_manual()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all normalization pipeline R interface tests passed! ✓\n")
