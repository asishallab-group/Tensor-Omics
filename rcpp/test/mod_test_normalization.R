# Load all TensorOmics helper functions
source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

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
  normalized_matrix_std <- tox_normalize_by_std_dev(df_matrix_clean)    # Normalize by standard deviation
  normalized_matrix_qtl <- tox_quantile_normalization(normalized_matrix_std)  # Apply quantile normalization
  normalized_matrix_log <- tox_log2_transformation(normalized_matrix_qtl)     # Log2(x+1) transformation
  averaged_df <- tox_calculate_tissue_averages(normalized_matrix_log)         # Average replicates by tissue

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
