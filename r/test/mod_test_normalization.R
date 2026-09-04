# Load all TensorOmics helper functions
source("r/load_tensor_omics.R")
source("r/test_helpers.R")
# tox_diagnose_data_quality / tox_clean_data_for_normalization are pure R -- no Fortran
# behind them -- so they live outside the generated binding
source("r/data_prep.R")

test_basic_calling <- function() {
  # === Example of full normalization pipeline ===

  # Load raw expression data
  input_file <- "material/kallisto_sex_data.tsv"  # Update the input file path if necessary
  output_file <- "results/normalization_sexdata.tsv"  # Output file path for normalized data

  df <- read.table(input_file, header = TRUE, sep = "\t")
  # Prepare matrix for processing (removing the gene ID column)
  gene_ids <- df[[1]]              # Save gene identifiers
  col_names <- colnames(df)[-1]    
  df_matrix <- as.matrix(df[,-1])   # Convert the expression values into a matrix

  # === Diagnose data quality before normalization ===
  diagnostics <- tox_diagnose_data_quality(df_matrix, show_details=FALSE)

  # Clean data if there are problems
  if (diagnostics$problems$na_count > 0 || diagnostics$problems$inf_count > 0 || diagnostics$problems$nan_count > 0) {
    # Clean the data using our cleaning function
    df_matrix_clean <- tox_clean_data_for_normalization(
      df_matrix,
      remove_all_zero_genes = TRUE,
      na_strategy = "impute_mean",  # Impute NA with gene means instead of removing genes
      min_expression_threshold = 0.0,  # Don't set a threshold for TPM data
      convert_small_to_zero = FALSE    # Preserve small TPM values
    )
    
    # Update gene_ids to match cleaned matrix
    gene_ids <- gene_ids[1:nrow(df_matrix_clean)]
  } else {
    df_matrix_clean <- df_matrix
  }
  # === Apply normalization steps sequentially ===
  normalized_matrix_std <- normalize_by_std_dev(df_matrix_clean)    # Normalize by standard deviation
  # quantile_normalization returns the rank means alongside the matrix
  normalized_matrix_qtl <- quantile_normalization(normalized_matrix_std)$normalized_expr  # Apply quantile normalization
  normalized_matrix_log <- log2_transformation(normalized_matrix_qtl)     # Log2(x+1) transformation
  # one replicate per row: the routine takes the counts rather than parsing row names
  averaged_df <- calc_tiss_avg(rep(1L, nrow(normalized_matrix_log)), normalized_matrix_log)

  # force(df_matrix_clean)
  invisible(force(normalized_matrix_std))
  invisible(force(normalized_matrix_qtl))
  invisible(force(normalized_matrix_log))
  invisible(force(averaged_df))

  # Create final normalized dataframe with averages
  normalized_df <- data.frame(gene_id = gene_ids, averaged_df)
  if (!dir.exists("results")) dir.create("results", recursive = TRUE)
  write.table(normalized_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

run_all_tests()
