# Comprehensive R test suite for normalization_pipeline (mirrors Fortran unit tests)
source("r/load_tensor_omics.R")
source("r/test_helpers.R")

test_basic <- function() {
  n_genes <- 10; n_tissues <- 6; n_grps <- 2
  input_matrix <- matrix(0, ncol=n_genes, nrow=n_tissues)
  for (i in 1:n_genes) {
    for (j in 1:n_tissues) {
      input_matrix[j, i] <- i + j * 0.25
    }
  }
  rownames(input_matrix) <- c(
    "muscle_dietM_1", "muscle_dietM_2", "muscle_dietM_3",
    "muscle_dietP_1", "muscle_dietP_2", "muscle_dietP_3"
  )
  colnames(input_matrix) <- paste0("gene", 1:n_genes)
  group_s <- c(1L,4L); group_c <- c(3L,3L)
  result <- normalization_pipeline(input_matrix, group_c, span = 0.75, degree = 2, use_quantile = TRUE)
  assert_true(all(!is.na(result)))
  assert_true(all(result >= 0))
  assert_true(dim(result)[2] == n_genes && dim(result)[1] == n_grps)
}

test_edge_case <- function() {
  n_genes <- 10; n_tissues <- 6; n_grps <- 2
  input_matrix <- matrix(0, ncol=n_genes, nrow=n_tissues)
  rownames(input_matrix) <- c(
    "muscle_dietM_1", "muscle_dietM_2", "muscle_dietM_3",
    "muscle_dietP_1", "muscle_dietP_2", "muscle_dietP_3"
  )
  colnames(input_matrix) <- paste0("gene", 1:n_genes)
  group_s <- c(1L,4L); group_c <- c(3L,3L)

  err <- tryCatch({
    normalization_pipeline(input_matrix, group_c, span = 0.75, degree = 2, use_quantile = TRUE)
    NULL
  }, error = function(e) e)

  assert_true(!is.null(err))
}

test_pipeline_vs_manual <- function() {
  n_genes <- 10; n_tissues <- 6; n_grps <- 2
  input_matrix <- matrix(0, ncol=n_genes, nrow=n_tissues)
  for (i in 1:n_genes) {
    for (j in 1:n_tissues) {
      input_matrix[j, i] <- i * 2 + j * 0.5
    }
  }
  rownames(input_matrix) <- c(
    "muscle_dietM_1", "muscle_dietM_2", "muscle_dietM_3",
    "muscle_dietP_1", "muscle_dietP_2", "muscle_dietP_3"
  )
  colnames(input_matrix) <- paste0("gene", 1:n_genes)
  group_s <- c(1L,4L); group_c <- c(3L,3L)

  buf_stddev <- normalize_by_std_dev(input_matrix)
  # quantile_normalization returns the rank means alongside the matrix
  buf_quant <- quantile_normalization(buf_stddev)$normalized_expr
  dimnames(buf_quant) <- dimnames(input_matrix)
  buf_avg <- calc_tiss_avg(group_c, buf_quant)
  manual_out <- log2_transformation(as.matrix(buf_avg))

  dimnames(buf_stddev) <- dimnames(input_matrix)
  buf_avg_no_quant <- calc_tiss_avg(group_c, buf_stddev)
  manual_out_no_quant <- log2_transformation(as.matrix(buf_avg_no_quant))

  result <- normalization_pipeline(input_matrix, group_c, span = 0.75, degree = 2, use_quantile = TRUE)
  assert_equal_numeric(as.vector(result), as.vector(manual_out), msg="test_pipeline_vs_manual")

  result_no_quant <- normalization_pipeline(input_matrix, group_c, span = 0.75, degree = 2, use_quantile = FALSE)
  assert_equal_numeric(as.vector(result_no_quant), as.vector(manual_out_no_quant), msg="test_pipeline_vs_manual_no_quant")
}

run_all_tests()