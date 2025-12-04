library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library

lib_path <- shQuote(normalizePath("build"))

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)

cat("✓ TensorOmics Rcpp functions loaded successfully\n")

source("rcpp/error_handling.R")

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
  validate_numeric_vector(vec1)
  validate_numeric_vector(vec2)
  validate_same_length(vec1, vec2)
  validate_nonempty(vec1)

  # Call Rcpp wrapper 
  return(tox_euclidean_distance_rcpp(as.numeric(vec1), as.numeric(vec2)))
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
  # R-layer validation in rcpp/ (kept here because r/ must not be changed)
  validate_numeric_vector(genes)
  validate_numeric_vector(centroids)
  validate_positive_integer_scalar(d)

#  # Validate flattened lengths are compatible with d
  validate_divisible_length(genes, d)
  validate_divisible_length(centroids, d)
  # Calculate dimensions
  n_genes <- as.integer(length(genes) / d)
  n_families <- as.integer(length(centroids) / d)
  validate_gene_to_family(gene_to_fam, n_genes, n_families)
  validate_length_equals_n(gene_to_fam, n_genes)

  return(tox_distance_to_centroid_rcpp(genes, centroids, gene_to_fam, d))
}


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
  # R-layer validation (kept in rcpp/ to avoid touching r/)
  validate_numeric_matrix(expression_vectors)
  n_axes <- nrow(as.matrix(expression_vectors))
  n_vectors <- ncol(as.matrix(expression_vectors))

  # Ensure selection vectors have correct lengths and types
  #Input validation
  validate_numeric_matrix(expression_vectors)

  # Ensure selectors have expected lengths
  validate_logical_or_index_vector(vector_selection)
  validate_logical_or_index_vector(axis_selection)

  #Convert to appropriate types for Rcpp
  if (is.numeric(vector_selection)) {
    vector_selection <- as.integer(as.logical(vector_selection))
  } else {
    vector_selection <- as.integer(vector_selection)
  }
  if (is.numeric(axis_selection)) {
    axis_selection <- as.integer(as.logical(axis_selection))
  } else {
    axis_selection <- as.integer(axis_selection)
  }

  # Call Rcpp wrapper
  result <- tox_calculate_tissue_versatility_rcpp(expression_vectors, vector_selection, axis_selection)
  
  # Check for errors
  check_err_code(result$ierr)
  
  # Return structured result 
  return(list(
    tissue_versatilities = result$tissue_versatilities,
    tissue_angles_deg = result$tissue_angles_deg,
    n_selected_vectors = result$n_selected_vectors,
    n_selected_axes = result$n_selected_axes
  ))

}
# ===================================================================
# OUTLIER DETECTION FUNCTIONS 
# ===================================================================
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
  # Input validation
  validate_numeric_vector(distances)
  n_genes <- as.integer(length(distances))

  # Call Rcpp wrapper
  result <- tox_detect_outliers_rcpp(distances, gene_to_fam, n_families, percentile)

  # Check for error (ierr returned by the Rcpp wrapper)
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    is_outlier = result$is_outlier,
    loess_x = result$loess_x,
    loess_y = result$loess_y,
    loess_n = result$loess_n
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
  # Input validation
  validate_numeric_vector(distances)
  n_genes <- as.integer(length(distances))
  validate_length_equals_n(gene_to_fam, n_genes)
  validate_index_bounds(gene_to_fam, low = 1, high = n_families)

  # Call the Rcpp forwarder.
  result <- tox_compute_family_scaling_rcpp(distances, gene_to_fam, n_families)

  # Check for error
  check_err_code(result$ierr)
  
  # Return structured result
  return(list(
    dscale = as.numeric(result$dscale),
    loess_x = as.numeric(result$loess_x),
    loess_y = as.numeric(result$loess_y),
    indices_used = as.integer(result$indices_used)
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
                                              perm_tmp, stack_left_tmp, stack_right_tmp,
                                              family_distances) {
# Input validation
  validate_numeric_vector(distances)
  n_genes <- as.integer(length(distances))
  validate_length_equals_n(gene_to_fam, n_genes)
  validate_index_bounds(gene_to_fam, low = 1, high = n_families)
  validate_length_equals_n(perm_tmp, n_genes)
  validate_length_equals_n(stack_left_tmp, n_genes)
  validate_length_equals_n(stack_right_tmp, n_genes)
  validate_length_equals_n(family_distances, n_genes)

  # Call the Rcpp forwarder.
  result <- tox_compute_family_scaling_expert_rcpp(n_families, distances, gene_to_fam, perm_tmp, stack_left_tmp, stack_right_tmp,
    family_distances
  )

   # Check for error
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    dscale = result$dscale,
    loess_x = result$loess_x,
    loess_y = result$loess_y,
    indices_used = result$indices_used,
    perm_tmp     = result$perm_tmp,
    stack_left_tmp   = result$stack_left_tmp,
    stack_right_tmp  = result$stack_right_tmp,
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
# Input validation
  validate_numeric_vector(distances)
  n_genes <- as.integer(length(distances))
  validate_length_equals_n(gene_to_fam, n_genes)
  n_families <- as.integer(length(dscale))
  validate_index_bounds(gene_to_fam, low = 1, high = n_families)
  # Call Rcpp forwarder
  result <- tox_compute_rdi_rcpp(distances, gene_to_fam, dscale)

  # Return 
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
  #Calculate dimensions
  n_genes <- as.integer(length(rdi))
  # Call the Rcpp forwarder
  result <- tox_identify_outliers_rcpp(rdi, percentile)
  # Return structured result
  return(list(
    is_outlier = result$is_outlier,
    threshold = result$threshold
  ))

}

# NORMALIZATION FUNCTIONS
# ===================================================================
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

  
  result <- tox_normalize_by_std_dev_rcpp(input_matrix)  
  return(matrix(result$output_vector, nrow = nrow(input_matrix), ncol = ncol(input_matrix), dimnames = dimnames(input_matrix)))
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
  validate_matrix(input_matrix)
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)
  result <- tox_quantile_normalization_rcpp(input_matrix)

  check_err_code(result$ierr)
  
  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
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
  validate_matrix(input_matrix)
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)
  result <- tox_log2_transformation_rcpp(input_matrix)
  
  check_err_code(result$ierr)
  
  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
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
  validate_matrix(df)

  tissue_groups <- as.character(sapply(colnames(df), tox_parse_tissue_group))
  unique_groups <- unique(tissue_groups)
  n_groups <- length(unique_groups)
  group_starts <- integer(n_groups)
  group_counts <- integer(n_groups)
  df_sorted <- df[, order(tissue_groups)]
  sorted_tissue_groups <- tissue_groups[order(tissue_groups)]
  current_group <- as.character(sorted_tissue_groups[1])
  group_starts[1] <- 1
  group_counts[1] <- 1
  group_idx <- 1
  for (i in 2:length(sorted_tissue_groups)) {
    if (!is.na(sorted_tissue_groups[i]) && !is.na(current_group) && as.character(sorted_tissue_groups[i]) == as.character(current_group)) {
      group_counts[group_idx] <- group_counts[group_idx] + 1
    } else {
      group_idx <- group_idx + 1
      group_starts[group_idx] <- i
      group_counts[group_idx] <- 1
      current_group <- as.character(sorted_tissue_groups[i])
    }
  }

  result <- tox_calc_tiss_avg_rcpp(df, group_starts, group_counts)
  
  check_err_code(result$ierr)

  n_genes <- nrow(df)
  n_cols <- length(result$output_vector) / n_genes
  output_matrix <- matrix(result$output_vector, nrow = n_genes, ncol = n_cols)

  colnames(output_matrix) <- unique_groups
  rownames(output_matrix) <- rownames(df)
  return(as.data.frame(output_matrix))
}


#' Calculate log2 fold changes based on control and condition patterns
#' @param df A data frame with genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern to detect control columns.
#' @param condition_patterns A character vector with patterns to detect condition columns.
tox_calculate_fc_by_patterns <- function(df, control_pattern, condition_patterns) {
  validate_matrix(df)
  validate_string_scalar(control_pattern)
  validate_character_vector(condition_patterns)

  # --- Identify control and condition columns ---
  indices_info <- tox_prepare_indices_by_patterns(df, control_pattern, condition_patterns)
  control_cols <- indices_info$control_cols
  condition_cols <- indices_info$condition_cols
  condition_labels <- indices_info$condition_labels

  n_pairs <- length(control_cols)

  result <- tox_calc_fchange_rcpp(df, control_cols, condition_cols)

  check_err_code(result$ierr)
  
  n_genes <- nrow(df)
  output_matrix <- matrix(result$output_vector, nrow = n_genes, ncol = n_pairs)

  colnames(output_matrix) <- condition_labels
  rownames(output_matrix) <- rownames(df)

  return(as.data.frame(output_matrix))
}


#' Complete normalization pipeline for gene expression data (up to log2(x+1))
#' @param input_matrix Numeric matrix (genes x tissues)
#' @param group_s Integer vector: start column index for each replicate group (1-based)
#' @param group_c Integer vector: number of columns per replicate group
tox_normalization_pipeline <- function(input_matrix, group_s, group_c) {
  validate_matrix(input_matrix)
  group_s <- as.integer(group_s)
  group_c <- as.integer(group_c)
  validate_group_vectors(group_s, group_c, ncol(input_matrix))

  result <- tox_normalization_pipeline_rcpp(input_matrix, group_s, group_c)
  
  check_err_code(result$ierr)
  

  return(matrix(result$buf_log, nrow = nrow(input_matrix), ncol = length(group_s)))
}

##################################################
## Helper Functions for normalization
##################################################


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
  # R-layer validation (kept in rcpp/)
  validate_numeric_matrix(expression_vectors)
  validate_numeric_matrix(family_centroids)
  # Ensure matching axes (rows)
  validate_matching_rows(expression_vectors, family_centroids)

  # gene_to_centroid should be integer vector with length equal to number of vectors
  gene_to_centroid <- as.integer(gene_to_centroid)
  n_vectors <- ncol(as.matrix(expression_vectors))
  validate_length_equals_n(gene_to_centroid, n_vectors)




  result <- tox_compute_shift_vector_field_rcpp(expression_vectors, family_centroids, gene_to_centroid)
  check_err_code(result$ierr)
  return(list(shift_vectors = result$shift_vectors))
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
  # R-layer validation (kept in rcpp/)
  validate_group_centroid_inputs(expression_vectors, gene_to_family, n_families, ortholog_set, mode)

  result <- tox_group_centroid_rcpp(expression_vectors, gene_to_family, n_families, ortholog_set, mode)
  check_err_code(result$ierr)
  return(result)
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
  # R-layer validation (kept in rcpp/)
  validate_mean_vector_inputs(expression_vectors, gene_indices)

  result <- tox_mean_vector_rcpp(expression_vectors, gene_indices)
  check_err_code(result$ierr)
  return(result)
}
