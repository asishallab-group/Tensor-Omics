# Load all TensorOmics helper functions
source("r/tensoromics_functions.R")

# === Example of full normalization pipeline ===

# Load raw expression data
input_file <- "material/kallisto_sex_data.tsv"  # Update the input file path if necessary
output_file <- "results/normalization_sexdata.tsv"  # Output file path for normalized data

df <- read.table(input_file, header = TRUE, sep = "\t")
print(head(df))  # Display the first few rows of the data
# Prepare matrix for processing (removing the gene ID column)
gene_ids <- df[[1]]              # Save gene identifiers
col_names <- colnames(df)[-1]    
df_matrix <- as.matrix(df[,-1])   # Convert the expression values into a matrix

# === Diagnose data quality before normalization ===
cat("Analyzing data quality before normalization...\n")
diagnostics <- tox_diagnose_data_quality(df_matrix)

# Clean data if there are problems
if (diagnostics$problems$na_count > 0 || diagnostics$problems$inf_count > 0 || diagnostics$problems$nan_count > 0) {
  cat("Data contains problematic values. Cleaning data...\n")
  
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
  
  cat("Data cleaning completed. Proceeding with normalization...\n")
} else {
  df_matrix_clean <- df_matrix
  cat("Data is clean. Proceeding with normalization...\n")
}

# === Apply normalization steps sequentially ===
normalized_matrix_std <- tox_normalize_by_std_dev(df_matrix_clean)    # Normalize by standard deviation
normalized_matrix_qtl <- tox_quantile_normalization(normalized_matrix_std)  # Apply quantile normalization
normalized_matrix_log <- tox_log2_transformation(normalized_matrix_qtl)     # Log2(x+1) transformation
averaged_df <- tox_calculate_tissue_averages(normalized_matrix_log)         # Average replicates by tissue

# # === Calculate fold changes between specified groups ===
# fc_df <- calculate_fc_by_patterns(
#   df = averaged_df,
#   control_pattern = "dietM",           # Define control pattern
#   condition_patterns = c("dietP")       # Define conditions pattern(s)
# )

# === Save the resulting normalized and transformed data ===
print("Sample of original matrix:")
print(head(df_matrix_clean))
print("Sample of normalized by std dev:")
print(head(normalized_matrix_std))
print("Sample of quantile normalized:")
print(head(normalized_matrix_qtl))
print("Sample of log2 transformed:")
print(head(normalized_matrix_log))
print("Sample of tissue averages:")
print(head(averaged_df))

# Create final normalized dataframe with averages
normalized_df <- data.frame(gene_id = gene_ids, averaged_df)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
write.table(normalized_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Normalization pipeline completed successfully!\n")
cat("Output saved to:", output_file, "\n")
cat("Final dimensions:", nrow(normalized_df), "genes x", ncol(normalized_df)-1, "conditions\n")
