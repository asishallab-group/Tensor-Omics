# === Load the shared library ===
# dyn.load("build/libtensor-omics.so")
dyn.load("/mnt/datos/Documentos/tensor_omics/check_code/tensor-omics/build/libtensor-omics.so")
#' Check error code and throw informative error if needed
#' 
#' @param ierr Error code from Fortran routine
tox_errors <- function(ierr) {
  msg <- switch(
    as.character(ierr),
    "0" = NULL,
    "101" = "File could not be opened",
    "102" = "Could not read magic number", 
    "103" = "Could not read array type code",
    "104" = "Could not read array dimension number",
    "105" = "Could not read array dimensions",
    "106" = "Could not read character length",
    "107" = "Could not read array data",
    "200" = "Invalid file format (magic number mismatch)",
    "201" = "Invalid input parameters",
    "202" = "No axes selected (empty input)",
    "5002" = "File not open or unit not connected",
    "9999" = "Unknown error",
    paste("Unknown Fortran error code:", ierr)
  )
  
  if (!is.null(msg)) {
    stop(msg)
  }
}

#' Normalize gene expression values by standard deviation
#'
#' This function wraps the Fortran subroutine `normalize_by_std_dev`
#' to normalize each gene's expression across tissues by the standard deviation.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A normalized numeric matrix with the same dimensions and names as the input.
#' @details
#' - The input matrix is flattened into a column-major vector.
#' - Calls a Fortran routine that normalizes each gene across tissues.
#' - Restores the original row and column names after normalization.
#'
#' @examples
#' normalized_matrix <- tox_normalize_by_std_dev(input_matrix)
tox_normalize_by_std_dev <- function(input_matrix) {
  n_genes <- nrow(input_matrix)  # Number of genes (rows)
  n_tissues <- ncol(input_matrix)  # Number of tissues (columns)

  # print(head(input_matrix))

  # Prepare the input vector (flatten matrix column-major) and allocate output space
  input_vector <- as.numeric(as.vector(input_matrix))
  output_vector <- numeric(n_genes * n_tissues)

  # Validate input data before calling Fortran
  if (any(is.na(input_vector))) {
    stop("Error: Input matrix contains NA values. Found ", sum(is.na(input_vector)), " NA values.")
  }
  if (any(is.infinite(input_vector))) {
    stop("Error: Input matrix contains infinite values. Found ", sum(is.infinite(input_vector)), " infinite values.")
  }
  if (any(is.nan(input_vector))) {
    stop("Error: Input matrix contains NaN values. Found ", sum(is.nan(input_vector)), " NaN values.")
  }

  # Call the Fortran subroutine
  result <- .Fortran("normalize_by_std_dev_r",
               as.integer(n_genes),
               as.integer(n_tissues),
               input_vector,
               output_vector)

  matrix(result[[4]], nrow = n_genes, ncol = n_tissues,
         dimnames = dimnames(input_matrix))

}

#' Quantile normalization of gene expression values
#'
#' This function wraps the Fortran subroutine `quantile_normalization`
#' to apply quantile normalization across all tissue columns.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A quantile-normalized numeric matrix with the same dimensions and names as the input.
#' @details
#' - The input matrix is flattened into a column-major vector.
#' - The Fortran subroutine sorts and aligns the distributions across tissues.
#' - After normalization, the original row and column names are restored.
#'
#' @examples
#' normalized_matrix <- tox_quantile_normalization(input_matrix)
tox_quantile_normalization <- function(input_matrix) {
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)

  # Flatten input matrix (column-major) and preallocate vectors
  input_vector <- as.numeric(as.vector(input_matrix))
  output_vector <- numeric(n_genes * n_tissues)
  temp_col <- numeric(n_genes)
  rank_means <- numeric(n_genes)
  perm <- integer(n_genes)

  # Estimar tamaño máximo para la pila (según pseudocódigo: log2(n) + 10)
  max_stack <- as.integer(ceiling(log2(n_genes)) + 10)
  stack_left <- integer(max_stack)
  stack_right <- integer(max_stack)

  # Fortran interop: asegurar tipos
  storage.mode(input_vector) <- "double"
  storage.mode(output_vector) <- "double"
  storage.mode(temp_col) <- "double"
  storage.mode(rank_means) <- "double"
  storage.mode(perm) <- "integer"
  storage.mode(stack_left) <- "integer"
  storage.mode(stack_right) <- "integer"

  # Initialize first stack entry manually
  # stack_left[1] <- 1L
  # stack_right[1] <- n_genes

  result <- .Fortran("quantile_normalization_r",
    as.integer(n_genes),
    as.integer(n_tissues),
    input_vector,
    output_vector,
    temp_col,
    rank_means,
    perm,
    stack_left,
    stack_right,
    as.integer(max_stack)
  )

  matrix(result[[4]], nrow = n_genes, ncol = n_tissues,
         dimnames = dimnames(input_matrix))
}




#' Apply log2(x + 1) transformation to gene expression values
#'
#' This function wraps the Fortran subroutine `log2_transformation`
#' to apply a log2(x + 1) transformation to each element in the input matrix.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A numeric matrix with log2-transformed expression values, 
#' preserving the same dimensions and names as the input.
#' @details
#' - The input matrix is flattened into a column-major vector.
#' - The Fortran subroutine applies a log2(x+1) transformation to each value.
#' - After transformation, the original row and column names are restored.
#'
#' @examples
#' log_matrix <- tox_log2_transformation(input_matrix)
tox_log2_transformation <- function(input_matrix) {
  n_genes <- nrow(input_matrix)  # Number of genes (rows)
  n_tissues <- ncol(input_matrix)  # Number of tissues (columns)

  # # Save column and row names for later restoration
  # col_names <- colnames(input_matrix)
  # row_names <- rownames(input_matrix)

  # Prepare the input vector (flatten matrix column-major) and allocate output space
  input_vector <- as.numeric(as.vector(input_matrix))
  output_vector <- numeric(n_genes * n_tissues)

  # Call the Fortran subroutine
  result <- .Fortran("log2_transformation_r",
               as.integer(n_genes),
               as.integer(n_tissues),
               input_vector,
               output_vector)

  # Reconstruct the transformed matrix
  matrix(result[[4]], nrow = n_genes, ncol = n_tissues,
  dimnames = dimnames(input_matrix))

  # # Restore row and column names
  # colnames(normalized_matrix) <- col_names
  # rownames(normalized_matrix) <- row_names

  # return(normalized_matrix)
}

#' Parse tissue group name from column name
#'
#' Helper function to extract the tissue group from a column name.
#' Handles various replicate naming patterns:
#' - "muscle_dietM_1" -> "muscle_dietM" 
#' - "Adipose_rep1" -> "Adipose"
#' - "Brain_rep2" -> "Brain"
#' - "tissue_condition_rep3" -> "tissue_condition"
#'
#' @param colname A string with the column name to parse.
#' @return A string representing the parsed tissue group name.
#' @examples
#' tox_parse_tissue_group("muscle_dietM_1") # returns "muscle_dietM"
#' tox_parse_tissue_group("Adipose_rep1")   # returns "Adipose"
#' tox_parse_tissue_group("brain_dietP")    # returns "brain_dietP"
tox_parse_tissue_group <- function(colname) {
  parts <- strsplit(colname, "_")[[1]]
  
  # Pattern 1: ends with just a number (e.g., "muscle_dietM_1")
  if (length(parts) >= 2 && grepl("^[0-9]+$", parts[length(parts)])) {
    return(paste(parts[1:(length(parts)-1)], collapse = "_"))
  }
  
  # Pattern 2: ends with "rep" followed by number (e.g., "Adipose_rep1")
  if (length(parts) >= 2 && grepl("^rep[0-9]+$", parts[length(parts)])) {
    return(paste(parts[1:(length(parts)-1)], collapse = "_"))
  }
  
  # Pattern 3: ends with "replicate" followed by number (e.g., "Tissue_replicate1")
  if (length(parts) >= 2 && grepl("^replicate[0-9]+$", parts[length(parts)])) {
    return(paste(parts[1:(length(parts)-1)], collapse = "_"))
  }
  
  # If no pattern matches, return full name
  return(colname)
}


#' Calculate average expression across replicates for each tissue group
#'
#' This function wraps the Fortran subroutine `calc_tiss_avg`
#' to compute the mean expression value for replicates grouped by tissue.
#'
#' @param df A data frame or matrix with genes as rows and tissue replicates as columns.
#' @return A data frame with genes as rows and averaged tissues as columns.
#' @details
#' - Replicate columns are grouped based on their parsed tissue group names.
#' - The input matrix is sorted according to groups before being processed.
#' - Calls a Fortran subroutine that calculates averages within each group.
#' - The result restores the original gene IDs as row names.
#'
#' @examples
#' averaged_df <- tox_calculate_tissue_averages(df)
tox_calculate_tissue_averages <- function(df) {
  n_genes <- nrow(df)      # Number of genes (rows)
  n_columns <- ncol(df)    # Number of columns (tissues)

  # --- Parse all column names to find their corresponding tissue group ---
  tissue_groups <- sapply(colnames(df), tox_parse_tissue_group)

  # --- Identify unique tissue groups ---
  unique_groups <- unique(tissue_groups)
  n_groups <- length(unique_groups)

  # --- Initialize mapping for groups ---
  group_starts <- integer(n_groups)
  group_counts <- integer(n_groups)

  # --- Sort the dataframe by tissue group name ---
  df_sorted <- df[, order(tissue_groups)]
  sorted_tissue_groups <- tissue_groups[order(tissue_groups)]

  current_group <- sorted_tissue_groups[1]
  group_starts[1] <- 1
  group_counts[1] <- 1
  group_idx <- 1

  # --- Build group_starts and group_counts arrays ---
  for (i in 2:length(sorted_tissue_groups)) {
    if (sorted_tissue_groups[i] == current_group) {
      group_counts[group_idx] <- group_counts[group_idx] + 1
    } else
      {
      group_idx <- group_idx + 1
      group_starts[group_idx] <- i
      group_counts[group_idx] <- 1
      current_group <- sorted_tissue_groups[i]
    }
  }

  # --- Prepare input vector and allocate output space ---
  input_vector <- as.numeric(as.vector(as.matrix(df_sorted)))
  output_vector <- numeric(n_genes * n_groups)

  # --- Call the Fortran subroutine to calculate averages ---
  result <- .Fortran("calc_tiss_avg_r",
               as.integer(n_genes),
               as.integer(n_groups),
               as.integer(group_starts),
               as.integer(group_counts),
               as.numeric(input_vector),
               as.numeric(output_vector))

  # --- Reconstruct output matrix ---
  output_matrix <- matrix(result[[6]], nrow = n_genes, ncol = n_groups)

  # --- Restore column and row names ---
  colnames(output_matrix) <- unique_groups
  rownames(output_matrix) <- rownames(df)

  return(as.data.frame(output_matrix))
}

#' Prepare control and condition column indices based on naming patterns
#'
#' This helper function searches for columns in the input dataframe that match
#' specified control and condition patterns. It builds the mapping necessary 
#' to calculate fold changes.
#'
#' @param df A data frame with expression data, genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern used to identify control columns.
#' @param condition_patterns A character vector with patterns used to identify condition columns.
#' @return A list containing:
#' \describe{
#'   \item{control_cols}{Indices of control columns}
#'   \item{condition_cols}{Indices of condition columns}
#'   \item{condition_labels}{Column names for the resulting fold changes}
#' }
#'
#' @examples
#' indices_info <- tox_prepare_indices_by_patterns(df, control_pattern = "dietM", condition_patterns = c("dietP"))
tox_prepare_indices_by_patterns <- function(df, control_pattern, condition_patterns) {
  colnames_df <- colnames(df)

  control_cols <- integer(0)         # Will store control column indices
  condition_cols <- integer(0)       # Will store condition column indices
  condition_labels <- character(0)   # Labels for resulting logFC columns

  # --- Find all control columns matching control_pattern ---
  control_candidates <- grep(control_pattern, colnames_df, value = TRUE)

  if (length(control_candidates) == 0) {
    stop("No control columns found matching pattern: ", control_pattern)
  }

  # --- For each control, find associated condition columns ---
  for (control in control_candidates) {
    tissue_prefix <- sub(paste0("_", control_pattern), "", control)

    for (cond_pattern in condition_patterns) {
      pattern_to_match <- paste0("^", tissue_prefix, "_", cond_pattern)
      matching_conditions <- grep(pattern_to_match, colnames_df, value = TRUE)

      for (condition in matching_conditions) {
        control_cols <- c(control_cols, which(colnames_df == control))
        condition_cols <- c(condition_cols, which(colnames_df == condition))
        condition_labels <- c(condition_labels, paste0(condition, "_logFC"))
      }
    }
  }

  if (length(control_cols) == 0) {
    stop("No valid control-condition pairs found!")
  }

  return(list(control_cols = control_cols, condition_cols = condition_cols, condition_labels = condition_labels))
}

#' Calculate log2 fold changes based on control and condition patterns
#'
#' This function wraps the Fortran subroutine `calc_fchange`
#' to calculate the fold changes between control and condition columns.
#'
#' @param df A data frame with genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern to detect control columns.
#' @param condition_patterns A character vector with patterns to detect condition columns.
#' @return A data frame with genes as rows and log2 fold change values as columns.
#'
#' @examples
#' fc_df <- tox_calculate_fc_by_patterns(df, control_pattern = "dietM", condition_patterns = c("dietP"))
tox_calculate_fc_by_patterns <- function(df, control_pattern, condition_patterns) {
  n_genes <- nrow(df)       # Number of genes
  n_columns <- ncol(df)     # Number of columns (conditions)

  # --- Identify control and condition columns ---
  indices_info <- tox_prepare_indices_by_patterns(df, control_pattern, condition_patterns)
  control_cols <- indices_info$control_cols
  condition_cols <- indices_info$condition_cols
  condition_labels <- indices_info$condition_labels
  
  n_pairs <- length(control_cols)

  # --- Prepare input and output vectors ---
  input_vector <- as.numeric(as.vector(as.matrix(df)))
  output_vector <- numeric(n_genes * n_pairs)

  print(control_cols)
  print(condition_cols)
  print(head(input_vector))
  # --- Call Fortran subroutine to calculate fold changes ---
  result <- .Fortran("calc_fchange_r",
               as.integer(n_genes),
               as.integer(n_columns),   # Pass n_cols as required by Fortran
               as.integer(n_pairs),
               as.integer(control_cols),
               as.integer(condition_cols),
               input_vector,
               output_vector)

  # --- Reconstruct the fold change matrix ---
  output_matrix <- matrix(result[[7]], nrow = n_genes, ncol = n_pairs)
  colnames(output_matrix) <- condition_labels
  rownames(output_matrix) <- rownames(df)

  return(as.data.frame(output_matrix))
}

#' Diagnose data quality issues in gene expression matrix
#'
#' This function examines the input matrix for common data quality issues
#' that could cause problems in downstream analysis.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @param show_details Logical indicating whether to show detailed information.
#' @return A list with diagnostic information about the data quality.
#' @examples
#' diagnostics <- tox_diagnose_data_quality(input_matrix)
tox_diagnose_data_quality <- function(input_matrix, show_details = TRUE) {
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)
  total_values <- n_genes * n_tissues
  
  # Check for different types of problematic values
  na_count <- sum(is.na(input_matrix))
  inf_count <- sum(is.infinite(input_matrix))
  nan_count <- sum(is.nan(input_matrix))
  zero_count <- sum(input_matrix == 0, na.rm = TRUE)
  negative_count <- sum(input_matrix < 0, na.rm = TRUE)
  
  # Find problematic genes (rows with issues)
  genes_with_na <- which(apply(input_matrix, 1, function(x) any(is.na(x))))
  genes_with_inf <- which(apply(input_matrix, 1, function(x) any(is.infinite(x))))
  genes_with_nan <- which(apply(input_matrix, 1, function(x) any(is.nan(x))))
  genes_all_zero <- which(apply(input_matrix, 1, function(x) all(x == 0, na.rm = TRUE)))
  
  # Summary statistics
  if (na_count == 0 && inf_count == 0 && nan_count == 0) {
    min_val <- min(input_matrix, na.rm = TRUE)
    max_val <- max(input_matrix, na.rm = TRUE)
    mean_val <- mean(input_matrix, na.rm = TRUE)
  } else {
    min_val <- NA
    max_val <- NA
    mean_val <- NA
  }
  
  diagnostics <- list(
    dimensions = c(genes = n_genes, tissues = n_tissues, total_values = total_values),
    problems = list(
      na_count = na_count,
      inf_count = inf_count,
      nan_count = nan_count,
      zero_count = zero_count,
      negative_count = negative_count
    ),
    problematic_genes = list(
      genes_with_na = genes_with_na,
      genes_with_inf = genes_with_inf,
      genes_with_nan = genes_with_nan,
      genes_all_zero = genes_all_zero
    ),
    statistics = list(
      min_val = min_val,
      max_val = max_val,
      mean_val = mean_val
    )
  )
  
  if (show_details) {
    cat("=== DATA QUALITY DIAGNOSTICS ===\n")
    cat("Matrix dimensions:", n_genes, "genes x", n_tissues, "tissues (", total_values, "total values)\n\n")
    
    cat("Problem summary:\n")
    cat("  - NA values:", na_count, "(", round(100*na_count/total_values, 2), "%)\n")
    cat("  - Infinite values:", inf_count, "(", round(100*inf_count/total_values, 2), "%)\n")
    cat("  - NaN values:", nan_count, "(", round(100*nan_count/total_values, 2), "%)\n")
    cat("  - Zero values:", zero_count, "(", round(100*zero_count/total_values, 2), "%)\n")
    cat("  - Negative values:", negative_count, "(", round(100*negative_count/total_values, 2), "%)\n\n")
    
    cat("Problematic genes:\n")
    cat("  - Genes with NA:", length(genes_with_na), "\n")
    cat("  - Genes with Inf:", length(genes_with_inf), "\n") 
    cat("  - Genes with NaN:", length(genes_with_nan), "\n")
    cat("  - Genes all zero:", length(genes_all_zero), "\n\n")
    
    if (na_count == 0 && inf_count == 0 && nan_count == 0) {
      cat("Data range:\n")
      cat("  - Min value:", min_val, "\n")
      cat("  - Max value:", max_val, "\n")
      cat("  - Mean value:", mean_val, "\n\n")
    }
    
    # Show examples of problematic genes
    if (length(genes_with_na) > 0) {
      cat("First few genes with NA values:\n")
      print(head(genes_with_na, 5))
      cat("\n")
    }
    
    if (length(genes_with_inf) > 0) {
      cat("First few genes with infinite values:\n")
      print(head(genes_with_inf, 5))
      cat("\n")
    }
  }
  
  return(invisible(diagnostics))
}

#' Clean data by removing or imputing problematic values
#'
#' This function handles NA, NaN, Inf values and genes that are all zeros
#' to prepare data for Fortran normalization routines.
#'
#' @param df_matrix A numeric matrix with genes as rows and tissues as columns
#' @param remove_all_zero_genes Logical, whether to remove genes that are all zeros
#' @param na_strategy Strategy for handling NA values: "remove_genes", "remove_samples", "impute_zero", "impute_mean"
#' @param min_expression_threshold Minimum expression value to consider (values below this become 0)
#' @return A cleaned matrix ready for normalization
tox_clean_data_for_normalization <- function(df_matrix, 
                                        remove_all_zero_genes = TRUE,
                                        na_strategy = "remove_genes",
                                        min_expression_threshold = 0.0,  # Changed default to 0.0
                                        convert_small_to_zero = FALSE) {   # New parameter to control this behavior
  
  cat("=== CLEANING DATA FOR NORMALIZATION ===\n")
  original_dims <- dim(df_matrix)
  cat("Original dimensions:", original_dims[1], "genes x", original_dims[2], "tissues\n")
  
  # Step 1: Handle very small values (only if explicitly requested)
  if (convert_small_to_zero && min_expression_threshold > 0.0) {
    small_values <- df_matrix > 0 & df_matrix < min_expression_threshold
    if (sum(small_values, na.rm = TRUE) > 0) {
      cat("Converting", sum(small_values, na.rm = TRUE), "values <", min_expression_threshold, "to zero\n")
      df_matrix[small_values] <- 0
    }
  } else {
    cat("Preserving all small values (convert_small_to_zero = FALSE)\n")
  }
  
  # Step 2: Handle infinite values
  inf_values <- is.infinite(df_matrix)
  if (sum(inf_values, na.rm = TRUE) > 0) {
    cat("WARNING: Converting", sum(inf_values, na.rm = TRUE), "infinite values to NA\n")
    df_matrix[inf_values] <- NA
  }
  
  # Step 3: Handle NaN values
  nan_values <- is.nan(df_matrix)
  if (sum(nan_values, na.rm = TRUE) > 0) {
    cat("WARNING: Converting", sum(nan_values, na.rm = TRUE), "NaN values to NA\n")
    df_matrix[nan_values] <- NA
  }
  
  # Step 4: Handle NA values according to strategy
  na_count <- sum(is.na(df_matrix))
  if (na_count > 0) {
    cat("Handling", na_count, "NA values using strategy:", na_strategy, "\n")
    
    if (na_strategy == "remove_genes") {
      # Remove genes with any NA values
      genes_with_na <- apply(df_matrix, 1, function(x) any(is.na(x)))
      df_matrix <- df_matrix[!genes_with_na, , drop = FALSE]
      cat("Removed", sum(genes_with_na), "genes with NA values\n")
      
    } else if (na_strategy == "remove_samples") {
      # Remove samples/tissues with any NA values
      samples_with_na <- apply(df_matrix, 2, function(x) any(is.na(x)))
      df_matrix <- df_matrix[, !samples_with_na, drop = FALSE]
      cat("Removed", sum(samples_with_na), "samples with NA values\n")
      
    } else if (na_strategy == "impute_zero") {
      # Replace NA with 0
      df_matrix[is.na(df_matrix)] <- 0
      cat("Imputed", na_count, "NA values with zero\n")
      
    } else if (na_strategy == "impute_mean") {
      # Replace NA with gene mean (row-wise)
      for (i in 1:nrow(df_matrix)) {
        na_positions <- is.na(df_matrix[i, ])
        if (any(na_positions)) {
          gene_mean <- mean(df_matrix[i, ], na.rm = TRUE)
          if (is.finite(gene_mean)) {
            df_matrix[i, na_positions] <- gene_mean
          } else {
            df_matrix[i, na_positions] <- 0  # If all values are NA, use 0
          }
        }
      }
      cat("Imputed", na_count, "NA values with gene means\n")
      
    } else if (na_strategy == "smart_impute") {
      # More sophisticated strategy: remove genes with >50% NA, impute the rest
      na_threshold <- 0.5  # Remove genes with more than 50% NA values
      
      genes_with_many_na <- apply(df_matrix, 1, function(x) {
        sum(is.na(x)) / length(x) > na_threshold
      })
      
      if (sum(genes_with_many_na) > 0) {
        df_matrix <- df_matrix[!genes_with_many_na, , drop = FALSE]
        cat("Removed", sum(genes_with_many_na), "genes with >", na_threshold*100, "% NA values\n")
      }
      
      # Impute remaining NA values with gene means
      remaining_na <- sum(is.na(df_matrix))
      if (remaining_na > 0) {
        for (i in 1:nrow(df_matrix)) {
          na_positions <- is.na(df_matrix[i, ])
          if (any(na_positions)) {
            gene_mean <- mean(df_matrix[i, ], na.rm = TRUE)
            if (is.finite(gene_mean)) {
              df_matrix[i, na_positions] <- gene_mean
            } else {
              df_matrix[i, na_positions] <- 0
            }
          }
        }
        cat("Imputed", remaining_na, "remaining NA values with gene means\n")
      }
    }
  }
  
  # Step 5: Handle all-zero genes
  if (remove_all_zero_genes) {
    all_zero_genes <- apply(df_matrix, 1, function(x) all(x == 0, na.rm = TRUE))
    if (sum(all_zero_genes) > 0) {
      df_matrix <- df_matrix[!all_zero_genes, , drop = FALSE]
      cat("Removed", sum(all_zero_genes), "genes with all zero values\n")
    }
  }
  
  # Final validation
  final_dims <- dim(df_matrix)
  cat("Final dimensions:", final_dims[1], "genes x", final_dims[2], "tissues\n")
  cat("Genes removed:", original_dims[1] - final_dims[1], "\n")
  cat("Samples removed:", original_dims[2] - final_dims[2], "\n")
  
  # Check for remaining problematic values
  remaining_na <- sum(is.na(df_matrix))
  remaining_inf <- sum(is.infinite(df_matrix))
  remaining_nan <- sum(is.nan(df_matrix))
  
  if (remaining_na > 0 || remaining_inf > 0 || remaining_nan > 0) {
    cat("WARNING: Still have problematic values:\n")
    cat("  NA:", remaining_na, "\n")
    cat("  Inf:", remaining_inf, "\n") 
    cat("  NaN:", remaining_nan, "\n")
  } else {
    cat("✓ Data is clean and ready for Fortran normalization\n")
  }
  
  return(df_matrix)
}

# ===================================================================
# TISSUE VERSATILITY FUNCTIONS
# ===================================================================

#' Calculate Tissue Versatility
#' 
#' Computes normalized tissue versatility for selected expression vectors.
#' The metric is based on the angle between each gene expression vector and the space diagonal.
#' Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
#' This function automatically checks for errors and throws informative exceptions.
#' 
#' @param expression_vectors Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param vector_selection Logical vector indicating which vectors to process (length n_vectors)
#' @param axis_selection Logical vector indicating which axes to include in calculation (length n_axes)
#' 
#' @return List containing:
#'   \item{tissue_versatilities}{Normalized tissue versatility values [0,1] for selected vectors}
#'   \item{tissue_angles_deg}{Angles in degrees [0,90] for selected vectors}
#'   \item{n_selected_vectors}{Number of vectors processed}
#'   \item{n_selected_axes}{Number of axes used in calculation}
#' 
tox_calculate_tissue_versatility <- function(expression_vectors, vector_selection, axis_selection) {
  # Input validation
  if (!is.matrix(expression_vectors)) {
    stop("expression_vectors must be a matrix")
  }
  if (!is.logical(vector_selection) && !is.numeric(vector_selection)) {
    stop("vector_selection must be logical or numeric")
  }
  if (!is.logical(axis_selection) && !is.numeric(axis_selection)) {
    stop("axis_selection must be logical or numeric")
  }
  
  # Convert to logical if numeric (0/1)
  if (is.numeric(vector_selection)) {
    vector_selection <- as.logical(vector_selection)
  }
  if (is.numeric(axis_selection)) {
    axis_selection <- as.logical(axis_selection)
  }
  
  # Dimensions and counts
  n_axes <- nrow(expression_vectors)
  n_vectors <- ncol(expression_vectors)
  n_selected_vectors <- sum(vector_selection)
  n_selected_axes <- sum(axis_selection)
  
  # Validate dimensions
  if (length(vector_selection) != n_vectors) {
    stop("vector_selection length must match number of columns in expression_vectors")
  }
  if (length(axis_selection) != n_axes) {
    stop("axis_selection length must match number of rows in expression_vectors")
  }
  
  # Prepare output arrays
  tissue_versatilities <- rep(0.0, n_selected_vectors)
  tissue_angles_deg <- rep(0.0, n_selected_vectors)
  ierr <- as.integer(0)
  
  # Call Fortran wrapper
  result <- .Fortran("compute_tissue_versatility_r",
                     n_axes = as.integer(n_axes),
                     n_vectors = as.integer(n_vectors),
                     expression_vectors = as.double(expression_vectors),
                     exp_vecs_selection_index = as.logical(vector_selection),
                     n_selected_vectors = as.integer(n_selected_vectors),
                     axes_selection = as.logical(axis_selection),
                     n_selected_axes = as.integer(n_selected_axes),
                     tissue_versatilities = as.double(tissue_versatilities),
                     tissue_angles_deg = as.double(tissue_angles_deg),
                     ierr = ierr)
  
  # Check for errors and throw informative messages
  tox_errors(result$ierr)
  
  # Return structured result (no ierr since we checked for errors)
  return(list(
    tissue_versatilities = result$tissue_versatilities,
    tissue_angles_deg = result$tissue_angles_deg,
    n_selected_vectors = n_selected_vectors,
    n_selected_axes = n_selected_axes
  ))
}


#' Compute family scaling factors using LOESS smoothing
#'
#' This function calculates scaling factors for gene families based on distance distributions.
#' It uses LOESS smoothing to estimate the relationship between family median distances
#' and their standard deviations.
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param n_families Integer number of families
#' @return List with components:
#'   - dscale: Scaling factors for each family
#'   - loess_x: Family median distances 
#'   - loess_y: Family standard deviations
#'   - indices_used: Number of genes used per family
tox_compute_family_scaling <- function(distances, gene_to_fam, n_families) {
  n_genes <- length(distances)
  
 
  # Prepare output arrays
  dscale <- double(n_families)
  loess_x <- double(n_families)
  loess_y <- double(n_families)
  indices_used <- integer(n_families)
  error_code <- integer(1)
  
  # Call the main interface Fortran wrapper
  result <- .Fortran("compute_family_scaling_r",
    n_genes = as.integer(n_genes),
    n_families = as.integer(n_families),
    distances = as.double(distances),
    gene_to_fam = as.integer(gene_to_fam),
    dscale = dscale,
    loess_x = loess_x,
    loess_y = loess_y,
    indices_used = indices_used,
    error_code = error_code
  )
  
  # Check for errors and throw informative messages
  tox_errors(result$error_code)
  
  return(list(
    dscale = result$dscale,
    loess_x = result$loess_x,
    loess_y = result$loess_y,
    indices_used = result$indices_used
  ))
}


#' Compute family scaling factors using LOESS smoothing (Expert Version)
#'
#' Expert version of compute_family_scaling with user-provided work arrays.
#' This version requires pre-allocated work arrays for maximum performance and control.
#' Use this when you need fine-grained control over memory allocation or are calling
#' this function many times in a tight loop.
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param n_families Integer number of families
#' @param perm_tmp Pre-allocated permutation array for sorting (n_genes)
#' @param stack_left_tmp Pre-allocated stack array for sorting (n_genes)
#' @param stack_right_tmp Pre-allocated stack array for sorting (n_genes)
#' @param family_distances Pre-allocated work array for family distances (n_genes)
#' @return List with components:
#'   - dscale: Scaling factors for each family
#'   - loess_x: Family median distances 
#'   - loess_y: Family standard deviations
#'   - indices_used: Number of genes used per family
#'   - perm_tmp: Final state of permutation array
#'   - stack_left_tmp: Final state of left stack array
#'   - stack_right_tmp: Final state of right stack array
#'   - family_distances: Final state of family distances array
tox_compute_family_scaling_expert <- function(distances, gene_to_fam, n_families, 
                                        perm_tmp, stack_left_tmp, stack_right_tmp, family_distances) {
  n_genes <- length(distances)
  
  # Validate inputs
  if (length(gene_to_fam) != n_genes) {
    stop("Length of gene_to_fam must equal length of distances")
  }
  if (any(gene_to_fam < 1 | gene_to_fam > n_families)) {
    stop("gene_to_fam indices must be between 1 and n_families")
  }
  if (length(perm_tmp) != n_genes) {
    stop("perm_tmp must have same length as distances")
  }
  if (length(stack_left_tmp) != n_genes) {
    stop("stack_left_tmp must have same length as distances")
  }
  if (length(stack_right_tmp) != n_genes) {
    stop("stack_right_tmp must have same length as distances")
  }
  if (length(family_distances) != n_genes) {
    stop("family_distances must have same length as distances")
  }
  
  # Prepare output arrays
  dscale <- numeric(n_families)
  loess_x <- numeric(n_families)
  loess_y <- numeric(n_families)
  indices_used <- integer(n_families)
  error_code <- integer(1)
  
  # Call Fortran expert routine
  result <- .Fortran("compute_family_scaling_expert_r",
    n_genes = as.integer(n_genes),
    n_families = as.integer(n_families),
    distances = as.double(distances),
    gene_to_fam = as.integer(gene_to_fam),
    dscale = dscale,
    loess_x = loess_x,
    loess_y = loess_y,
    indices_used = indices_used,
    perm_tmp = as.integer(perm_tmp),
    stack_left_tmp = as.integer(stack_left_tmp),
    stack_right_tmp = as.integer(stack_right_tmp),
    family_distances = as.double(family_distances),
    error_code = error_code
  )
  
  # Check for errors and throw informative messages
  tox_errors(result$error_code)
  
  return(list(
    dscale = result$dscale,
    loess_x = result$loess_x,
    loess_y = result$loess_y,
    indices_used = result$indices_used,
    perm_tmp = result$perm_tmp,
    stack_left_tmp = result$stack_left_tmp,
    stack_right_tmp = result$stack_right_tmp,
    family_distances = result$family_distances
  ))
}


#' Compute Relative Distance Index (RDI) for genes
#'
#' This function calculates the Relative Distance Index for each gene,
#' which is the absolute distance divided by the family scaling factor.
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param dscale Numeric vector of scaling factors for each family
#' @return List with components:
#'   - rdi: Relative Distance Index for each gene
#'   - sorted_rdi: RDI values sorted in ascending order
tox_compute_rdi <- function(distances, gene_to_fam, dscale) {
  n_genes <- length(distances)
  n_families <- length(dscale)
  
  # Validate inputs
  if (length(gene_to_fam) != n_genes) {
    stop("Length of gene_to_fam must equal length of distances")
  }
  if (any(gene_to_fam < 1 | gene_to_fam > n_families)) {
    stop("gene_to_fam indices must be between 1 and n_families")
  }
  
  # Prepare output arrays and work arrays
  rdi <- double(n_genes)
  sorted_rdi <- double(n_genes)
  perm <- as.integer(seq_along(distances))
  stack_left <- integer(n_genes)
  stack_right <- integer(n_genes)
  
  # Call the Fortran wrapper
  result <- .Fortran("compute_rdi_r",
    n_genes = as.integer(n_genes),
    n_families = as.integer(n_families),
    distances = as.double(distances),
    gene_to_fam = as.integer(gene_to_fam),
    dscale = as.double(dscale),
    rdi = rdi,
    sorted_rdi = sorted_rdi,
    perm = perm,
    stack_left = stack_left,
    stack_right = stack_right
  )
  
  return(list(
    rdi = result$rdi,
    sorted_rdi = result$sorted_rdi
  ))
}


#' Identify outliers based on RDI percentiles
#'
#' This function identifies outliers by comparing each gene's RDI value
#' against a percentile threshold of the sorted RDI distribution.
#'
#' @param rdi Numeric vector of RDI values
#' @param percentile Percentile threshold (default: 95.0)
#' @return List with components:
#'   - is_outlier: Logical vector indicating outliers
#'   - threshold: The RDI threshold value used
tox_identify_outliers <- function(rdi, percentile = 95.0) {
  n_genes <- length(rdi)
  

  # Prepare sorted RDI array (only non-negative values)
  sorted_rdi <- sort(rdi[rdi >= 0])
  if (length(sorted_rdi) < n_genes) {
    sorted_rdi <- c(sorted_rdi, rep(0, n_genes - length(sorted_rdi)))
  }
  
  # Prepare output arrays
  is_outlier <- logical(n_genes)
  threshold <- double(1)
  
  # Call the Fortran wrapper
  result <- .Fortran("identify_outliers_r",
    n_genes = as.integer(n_genes),
    rdi = as.double(rdi),
    sorted_rdi = as.double(sorted_rdi),
    is_outlier = is_outlier,
    threshold = threshold,
    percentile = as.double(percentile)
  )
  
  return(list(
    is_outlier = result$is_outlier,
    threshold = result$threshold
  ))
}


#' Complete outlier detection workflow
#'
#' This function performs the complete outlier detection workflow:
#' 1. Computes family scaling factors using LOESS
#' 2. Calculates RDI values
#' 3. Identifies outliers based on percentile threshold
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param n_families Integer number of families
#' @param percentile Percentile threshold for outlier detection (default: 95.0)
#' @return List with components:
#'   - is_outlier: Logical vector indicating outliers
#'   - loess_x: Family median distances
#'   - loess_y: Family standard deviations
#'   - loess_n: Number of genes used per family
tox_detect_outliers <- function(distances, gene_to_fam, n_families, percentile = 95.0) {
  n_genes <- length(distances)
  
  
  # Prepare output arrays and work arrays
  work_array <- double(n_genes)
  perm <- seq_len(n_genes)
  stack_left <- integer(n_genes)
  stack_right <- integer(n_genes)
  is_outlier <- logical(n_genes)
  loess_x <- double(n_families)
  loess_y <- double(n_families)
  loess_n <- integer(n_families)
  error_code <- integer(1)
  
  # Call the comprehensive Fortran wrapper
  result <- .Fortran("detect_outliers_r",
    n_genes = as.integer(n_genes),
    n_families = as.integer(n_families),
    distances = as.double(distances),
    gene_to_fam = as.integer(gene_to_fam),
    work_array = work_array,
    perm = perm,
    stack_left = stack_left,
    stack_right = stack_right,
    is_outlier = is_outlier,
    loess_x = loess_x,
    loess_y = loess_y,
    loess_n = loess_n,
    error_code = error_code,
    percentile = as.double(percentile)
  )
  
  # Check for errors and throw informative messages
  tox_errors(result$error_code)
  
  return(list(
    is_outlier = result$is_outlier,
    loess_x = result$loess_x,
    loess_y = result$loess_y,
    loess_n = result$loess_n
  ))
}

#' Perform 2D LOESS smoothing
#'
#' This function wraps the Fortran subroutine `loess_smooth_2d_r`
#' to perform 2D LOESS smoothing with kernel weighting.
#' Error handling is performed by Fortran and reported via error codes.
#'
#' @param x_ref Reference x coordinates (numeric vector)
#' @param y_ref Reference y values (numeric matrix, 1 row)
#' @param x_query Query x coordinates where smoothed values are needed (numeric vector)
#' @param indices_used Indices of reference points to use (integer vector, 1-based)
#' @param kernel_sigma Kernel bandwidth parameter (numeric)
#' @param kernel_cutoff Kernel cutoff distance (numeric)
#' @return A list containing the smoothed y values at query points
#' @details
#' - Uses a Gaussian kernel for distance weighting
#' - Only uses reference points specified by indices_used
#' - Returns smoothed values at all query points
#' - Fortran handles all validation and error checking
#'
tox_loess_smooth_2d <- function(x_ref, y_ref, x_query, indices_used = NULL, 
                           kernel_sigma, kernel_cutoff) {
  
  # Ensure y_ref is a matrix with proper dimensions
  if (!is.matrix(y_ref)) {
    y_ref <- matrix(y_ref, nrow = 1)
  }
  
  n_total <- as.integer(length(x_ref))
  n_target <- as.integer(length(x_query))
  
  # If indices_used not provided, use all points
  if (is.null(indices_used)) {
    indices_used <- 1:n_total
  }
  n_used <- as.integer(length(indices_used))
  
  # Prepare output matrix and error code
  y_out <- matrix(0.0, nrow = 1, ncol = n_target)
  ierr <- as.integer(0)
  
  # Ensure all inputs are proper types for Fortran
  x_ref <- as.double(x_ref)
  y_ref <- matrix(as.double(y_ref), nrow = 1)
  x_query <- as.double(x_query)
  indices_used <- as.integer(indices_used)
  kernel_sigma <- as.double(kernel_sigma)
  kernel_cutoff <- as.double(kernel_cutoff)
  
  # Call Fortran subroutine - let Fortran handle all validation
  result <- .Fortran("loess_smooth_2d_r",
    n_total = n_total,
    n_target = n_target,
    x_ref = x_ref,
    y_ref = y_ref,
    indices_used = indices_used,
    n_used = n_used,
    x_query = x_query,
    kernel_sigma = kernel_sigma,
    kernel_cutoff = kernel_cutoff,
    y_out = y_out,
    ierr = ierr
  )
  
  # Check for errors and throw informative messages
  tox_errors(result$ierr)
  
  return(list(
    y_out = result$y_out,
    smoothed_values = as.vector(result$y_out)
  ))
}

# ===================================================================
# EUCLIDEAN DISTANCE FUNCTIONS
# ===================================================================

#' Calculate Euclidean distance between two vectors
#' 
#' Computes the Euclidean distance between two vectors of the same dimension.
#' This function automatically checks for errors and throws informative exceptions.
#' 
#' @param vec1 First vector (numeric)
#' @param vec2 Second vector (numeric, same length as vec1)
#' 
#' @return Numeric value representing the Euclidean distance between the vectors
#' 
tox_euclidean_distance <- function(vec1, vec2) {
  # Input validation
  if (!is.numeric(vec1) || !is.numeric(vec2)) {
    stop("Both vectors must be numeric")
  }
  if (length(vec1) != length(vec2)) {
    stop("Vectors must have the same length")
  }
  if (length(vec1) == 0) {
    stop("Vectors cannot be empty")
  }
  
  d <- as.integer(length(vec1))
  result <- 0.0
  
  # Call Fortran wrapper
  fortran_result <- .Fortran("euclidean_distance_r",
                            vec1 = as.double(vec1),
                            vec2 = as.double(vec2),
                            d = d,
                            result = as.double(result))
  
  return(fortran_result$result)
}

#' Calculate distances from genes to their family centroids
#' 
#' Computes the Euclidean distance from each gene to its corresponding family centroid.
#' This function automatically checks for errors and throws informative exceptions.
#' 
#' @param genes Matrix of gene expression data (genes as columns, dimensions as rows)
#' @param centroids Matrix of family centroids (families as columns, dimensions as rows) 
#' @param gene_to_fam Integer vector mapping each gene to its family index (1-based)
#' @param d Integer number of dimensions
#' 
#' @return Numeric vector of distances from each gene to its family centroid
#' 
tox_distance_to_centroid <- function(genes, centroids, gene_to_fam, d) {
  # Input validation
  if (!is.numeric(genes) || !is.numeric(centroids)) {
    stop("genes and centroids must be numeric")
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    stop("gene_to_fam must be numeric or integer")
  }
  if (!is.numeric(d) && !is.integer(d)) {
    stop("d must be numeric or integer")
  }
  
  # Convert to appropriate types
  genes <- as.double(genes)
  centroids <- as.double(centroids)
  gene_to_fam <- as.integer(gene_to_fam)
  d <- as.integer(d)
  
  # Calculate dimensions
  n_genes <- as.integer(length(genes) / d)
  n_families <- as.integer(length(centroids) / d)
  
  # Validate dimensions
  if (length(genes) %% d != 0) {
    stop("Length of genes must be divisible by d")
  }
  if (length(centroids) %% d != 0) {
    stop("Length of centroids must be divisible by d")
  }
  if (length(gene_to_fam) != n_genes) {
    stop("Length of gene_to_fam must equal number of genes")
  }
  if (any(gene_to_fam < 0)) {
    stop("gene_to_fam indices must be between 0 and n_families (0 = no family assignment)")
  }
  
  # Prepare output array
  distances <- rep(0.0, n_genes)
  
  # Call Fortran wrapper
  result <- .Fortran("distance_to_centroid_r",
                    n_genes = n_genes,
                    n_families = n_families,
                    genes = genes,
                    centroids = centroids,
                    gene_to_fam = gene_to_fam,
                    distances = as.double(distances),
                    d = d)
  
  # Fortran returns -1 for genes without valid family assignment
  
  return(result$distances)
}

# ===================================================================
# SHIFT VECTOR FIELD FUNCTIONS
# ===================================================================
#' Calculate Shift Vector Field 

#' Computes the shift vector field for each gene expression vector based on its family centroid.
#' The shift vector is defined as the difference between the gene expression vector and its corresponding family centroid,
#' starting at the expression vector and pointing to its family centroid.
#' This function automatically checks for errors and throws informative exceptions.
#'
#' @param expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param family_centroids: Matrix where each column is a family centroid vector (n_axes x n_families)
#' @param gene_to_centroid: Array mapping each gene to its corresponding family centroid ID in family_centroids (length n_vectors)
#'
#' @return List containing:
#'   \item{shift_vectors}{The computed shift vectors for each gene expression vector}
#'

tox_compute_shift_vector_field <- function(expression_vectors, family_centroids, gene_to_centroid) {
  # Input validation
  if (!is.matrix(expression_vectors)) {
    stop("expression_vectors must be a matrix")
  }

  if (!is.matrix(family_centroids)) {
    stop("family_centroids must be a matrix")
  }
  
  # Dimensions and counts
  n_axes_genes <- nrow(expression_vectors)
  n_vectors <- ncol(expression_vectors)
  n_axes_centroids <- nrow(family_centroids)
  n_families <- ncol(family_centroids)
  
  # Validate length of gene_to_centroid
  if (n_vectors != length(gene_to_centroid)) {
    stop("number of expression_vectors must be equal to length of gene_to_centroid")
  }

  # Validate dimensions
  if (n_axes_genes != n_axes_centroids) {
    stop("family_centroids must have the same number of axes as expression_vectors")
  }
  
  # Prepare output arrays
  shift_vectors <- matrix(0.0, nrow = 2 *n_axes_genes, ncol = n_vectors)
  ierr <- as.integer(0)
  
  # Call Fortran wrapper
  result <- .Fortran("compute_shift_vector_field_r",
                     n_axes_genes = as.integer(n_axes_genes),
                     n_vectors = as.integer(n_vectors),
                     n_families = as.integer(n_families),
                     expression_vectors = as.double(expression_vectors),
                     family_centroids = as.double(family_centroids),
                     gene_to_centroid = as.integer(gene_to_centroid),
                     shift_vectors = as.double(shift_vectors),
                     ierr = ierr)
  
  # Check for errors and throw informative messages
  tox_errors(result$ierr)
  
  # Return structured result (no ierr since we checked for errors)
  return(list(
    shift_vectors = result$shift_vectors
  ))
}

# ===================================================================
# GENE CENTROIDS FUNCTIONS
# ===================================================================
#' Calculate Gene Centroids

#' Computes the centroids for each gene family based on the expression vectors of its member genes.
#' This function automatically checks for errors and throws informative exceptions.
#'
#' @param expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param gene_to_family: Array mapping each gene to its corresponding family ID (length n_vectors)
#' @param n_families: Total number of gene families
#' @param ortholog_set: Logical array indicating if a gene is part of a specific subset (e.g., orthologs)
#' @param mode: Character string indicating the mode of operation ('all' or 'ortho')
#'
#' @return List containing:
#'   \item{centroid_matrix}{The computed centroids for each gene family}
#'

tox_group_centroid <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  
  # 1) Validate inputs
  if (!is.matrix(expression_vectors) || !is.numeric(expression_vectors)) {
    stop("`expression_vectors` must be a numeric matrix.")
  }
  n_axes <- nrow(expression_vectors)
  n_genes <- ncol(expression_vectors)

  if (!is.integer(gene_to_family) || length(gene_to_family) != n_genes) {
    stop("`gene_to_family` must be an integer vector of length n_genes.")
  }
  if (!is.logical(ortholog_set) || length(ortholog_set) != n_genes) {
    stop("`ortholog_set` must be a logical vector of length n_genes.")
  }
  if (!mode %in% c('all', 'ortho')) {
    stop("`mode` must be either 'all' or 'ortho'.")
  }

  # 2) Prepare inputs/outputs for Fortran
  use_all_mode <- (mode == 'all')
  centroid_matrix_out <- matrix(0.0, nrow = n_axes, ncol = n_families)
  selected_indices_ws <- integer(n_genes) # Workspace buffer
  ierr <- as.integer(0)
  # 3) Call Fortran
  result <- .Fortran("group_centroid_r",
                     expression_vectors = as.double(expression_vectors),
                     n_axes = as.integer(n_axes),
                     n_genes = as.integer(n_genes),
                     gene_to_family = as.integer(gene_to_family),
                     num_families = as.integer(n_families),
                     centroid_matrix = centroid_matrix_out,
                     use_all_mode = as.logical(use_all_mode),
                     ortholog_set = as.logical(ortholog_set),
                     selected_indices = selected_indices_ws,
                     selected_indices_len = as.integer(n_genes),
                     ierr = ierr)
  
  # Check for errors and throw informative messages
  tox_errors(result$ierr)

  # 4) Return the populated output matrix (no ierr since we checked for errors)
  return(result$centroid_matrix)
}

#' Compute the element-wise mean for a given set of gene expression vectors
#'
#' This function wraps the Fortran subroutine `mean_vector_r`
#' to compute the centroid (mean vector) for a selected set of genes.
#'
#' @param expression_vectors Numeric matrix (n_axes x n_genes) of gene expression vectors
#' @param gene_indices Integer vector of column indices of selected genes (1-based)
#'
#' @return Numeric vector of length n_axes representing the computed centroid
#'

tox_mean_vector <- function(expression_vectors, gene_indices) {
  # Validate inputs
  if (!is.matrix(expression_vectors) || !is.numeric(expression_vectors)) {
    stop("`expression_vectors` must be a numeric matrix.")
  }
  n_axes <- nrow(expression_vectors)
  n_genes <- ncol(expression_vectors)
  n_selected_genes <- length(gene_indices)
  if (!is.integer(gene_indices) || any(gene_indices < 1) || any(gene_indices > n_genes)) {
    stop("`gene_indices` must be integer indices between 1 and n_genes.")
  }

  centroid_col <- numeric(n_axes)
  ierr <- as.integer(0)
  
  # Call Fortran wrapper
  result <- .Fortran("mean_vector_r",
                     expression_vectors = as.double(expression_vectors),
                     n_axes = as.integer(n_axes),
                     n_genes = as.integer(n_genes),
                     gene_indices = as.integer(gene_indices),
                     n_selected_genes = as.integer(n_selected_genes),
                     centroid_col = centroid_col,
                     ierr = ierr)
  
  # Check for errors and throw informative messages
  tox_errors(result$ierr)
  
  return(result$centroid_col)
}