
#> f42_helper-import_libs: Import necessary packages
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

#> tox_euclidean_distance:euclidean_distance_c: Calculate Euclidean distance between two vectors
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
  validate_is_numeric_vector(vec1)
  validate_is_numeric_vector(vec2)
  validate_is_same_length(vec1, vec2)
  validate_is_nonempty(vec1)

  # Call Rcpp wrapper
  return(tox_euclidean_distance_rcpp(as.numeric(vec1), as.numeric(vec2)))
}


#> tox_euclidean_distance:distance_to_centroid_c: Calculate distance from each gene to its family centroid
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
  validate_is_numeric_vector(genes)
  validate_is_numeric_vector(centroids)
  validate_is_positive_integer_scalar(d)

#  # Validate flattened lengths are compatible with d
  validate_is_divisible_length(genes, d)
  validate_is_divisible_length(centroids, d)
  # Calculate dimensions
  n_genes <- as.integer(length(genes) / d)
  n_families <- as.integer(length(centroids) / d)
  validate_is_gene_to_family(gene_to_fam, n_genes, n_families)
  validate_is_length_equals(gene_to_fam, n_genes)

  return(tox_distance_to_centroid_rcpp(genes, centroids, gene_to_fam, d))
}


#> tox_tissue_versatility:compute_tissue_versatility_c: Computes normalized tissue versatility for selected expression vectors
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
  validate_is_numeric_matrix(expression_vectors)
  n_axes <- nrow(as.matrix(expression_vectors))
  n_vectors <- ncol(as.matrix(expression_vectors))

  # Ensure selection vectors have correct lengths and types
  #Input validation
  validate_is_numeric_matrix(expression_vectors)

  # Ensure selectors have expected lengths
  validate_is_logical_or_index_vector(vector_selection)
  validate_is_logical_or_index_vector(axis_selection)

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

#> tox_get_outliers:detect_outliers_c: Complete outlier detection pipeline
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

#> tox_get_outliers:compute_family_scaling_c: Compute family scaling factors for outlier detection
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
  validate_is_numeric_vector(distances)
  n_genes <- as.integer(length(distances))
  validate_is_length_equals(gene_to_fam, n_genes)
  validate_is_in_range(gene_to_fam, min = 1, max = n_families)

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

#> tox_get_outliers:compute_family_scaling_expert_c: Compute family scaling factors using LOESS smoothing (Expert Version)
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
  validate_is_numeric_vector(distances)
  n_genes <- as.integer(length(distances))
  validate_is_length_equals(gene_to_fam, n_genes)
  validate_is_in_range(gene_to_fam, min = 1, max = n_families)
  validate_is_length_equals(perm_tmp, n_genes)
  validate_is_length_equals(stack_left_tmp, n_genes)
  validate_is_length_equals(stack_right_tmp, n_genes)
  validate_is_length_equals(family_distances, n_genes)

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

#> tox_get_outliers:compute_rdi_c: Compute Relative Distance Index (RDI) for outlier detection
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

#> tox_get_outliers:identify_outliers_c: Identify outliers based on RDI percentile or threshold
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

# ===================================================================
# NORMALIZATION FUNCTIONS
# ===================================================================

#> tox_normalization:normalize_by_std_dev_c: Normalize gene expression values by standard deviation
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

#> tox_normalization:quantile_normalization_c: Quantile normalization of gene expression values
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
  validate_numeric_matrix(input_matrix)
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)
  result <- tox_quantile_normalization_rcpp(input_matrix)

  check_err_code(result$ierr)

  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
}

#> tox_normalization:log2_transformation_c: Apply log2(x + 1) transformation to gene expression values
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
  validate_numeric_matrix(input_matrix)
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)
  result <- tox_log2_transformation_rcpp(input_matrix)

  check_err_code(result$ierr)

  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
}

#> tox_normalization:calc_tiss_avg_c: Calculate average expression across replicates for each tissue group
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
  validate_numeric_matrix(df)

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


#> tox_normalization:calc_fchange_c: Calculate log2 fold changes between control and condition columns
#' Calculate log2 fold changes based on control and condition patterns
#' @param df A data frame with genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern to detect control columns.
#' @param condition_patterns A character vector with patterns to detect condition columns.
tox_calculate_fold_changes <- function(df, control_pattern, condition_patterns) {
  validate_numeric_matrix(df)
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


#> tox_normalization:normalization_pipeline_c: Complete normalization pipeline for gene expression data (up to log2(x+1))
#' Complete normalization pipeline for gene expression data (up to log2(x+1))
#' @param input_matrix Numeric matrix (genes x tissues)
#' @param group_s Integer vector: start column index for each replicate group (1-based)
#' @param group_c Integer vector: number of columns per replicate group
tox_normalization_pipeline <- function(input_matrix, group_s, group_c) {
  validate_numeric_matrix(input_matrix)
  validate_group_vectors(group_s, group_c, ncol(input_matrix))

  result <- tox_normalization_pipeline_rcpp(
    input_matrix,
    as.integer(group_s),
    as.integer(group_c)
  )

  check_err_code(result$ierr)


  return(matrix(result$buf_log, nrow = nrow(input_matrix), ncol = length(group_s)))
}

##################################################
## Helper Functions for normalization
##################################################

#> f42_helper: Parse tissue group name from column name
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


#> f42_helper: Diagnose data quality issues in gene expression matrix
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

#> f42_helper: Clean data by removing or imputing problematic values
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

#> f42_helper: Prepare control and condition column indices based on naming patterns
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

#> tox_shift_vectors:compute_shift_vector_field_c: Computes the shift vector field for each gene expression vector based on its family centroid
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
  validate_matching_rows(expression_vectors, family_centroids)

  # Validate gene_to_centroid before casting
  n_vectors <- ncol(as.matrix(expression_vectors))
  validate_length_equals_n(gene_to_centroid, n_vectors)

  result <- tox_compute_shift_vector_field_rcpp(
    expression_vectors,
    family_centroids,
    as.integer(gene_to_centroid)
  )
  check_err_code(result$ierr)
  return(list(shift_vectors = result$shift_vectors))
}

# ===================================================================
# GENE CENTROIDS FUNCTIONS
# ===================================================================


#> tox_gene_centroids:group_centroid_c: Computes expression centroids for groups of genes
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

#> tox_gene_centroids:mean_vector_c: Compute the element-wise mean for a given set of gene expression vectors
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

# ============================================================
#  1) Array metadata wrapper
# ============================================================

#> f42_array_utils:get_array_metadata_C: Get array metadata from serialized file
tox_get_array_metadata <- function(filename, max_dims = 5L, with_clen = FALSE) {
  # Coerce to expected base types
  filename  <- as.character(filename)
  max_dims  <- as.integer(max_dims)
  with_clen <- as.logical(with_clen)

 
  # Expected to return a list with names dims, ndim, and optionally clen
  result <- tox_get_array_metadata_rcpp(filename, max_dims, with_clen)

  check_err_code(result$ierr)

  # Keep the public interface consistent with the old R function
  if (isTRUE(with_clen)) {
    return(list(
      dims = result$dims,
      ndim = result$ndim,
      clen = result$clen
    ))
  } else {
    return(list(
      dims = result$dims,
      ndim = result$ndim
    ))
}
}



# ============================================================
#  2) Deserialization (int / real / char)
# ============================================================

#> f42_deserialize_int:deserialize_int_nd_C: Deserialize integer array from binary file
tox_deserialize_int_array <- function(filename, max_dims = 5L) {
  # validate inputs
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)

  meta <- tox_get_array_metadata(filename, max_dims)
  total_size <- prod(meta$dims)

  filename_ascii <- utf8ToInt(filename)
 # Call Rcpp wrapper
  result <- tox_deserialize_int_array_rcpp(filename, max_dims)
  
  check_err_code(result$ierr)

  array(result$values, dim = result$dims[1:result$ndim])
}


#> f42_deserialize_real:deserialize_real_nd_C: Deserialize real array from binary file
tox_deserialize_real_array <- function(filename, max_dims = 5L) {
#validate inputs
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)

  meta <- tox_get_array_metadata(filename, max_dims)
  total_size <- prod(meta$dims)

  # Coerce to base types
  filename <- as.character(filename)
  max_dims <- as.integer(max_dims)

 
  # Call Rcpp wrapper
  result <- tox_deserialize_real_array_rcpp(filename, max_dims)

  # Check error code
  check_err_code(result$ierr)

  # Shape flat vector into array
  array(result$values, dim = result$dims[1:result$ndim])
  }

#> f42_deserialize_char:deserialize_char_nd_C: Deserialize character array from binary file
tox_deserialize_char_array <- function(filename, max_dims = 5L) {
  #validate inputs
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)
  # Load metadata dimensions + clen
  meta <- tox_get_array_metadata(filename, max_dims, with_clen = TRUE)
  
  # Coerce to base types
  filename <- as.character(filename)
  max_dims <- as.integer(max_dims)

  # Call Rcpp wrapper
  result <- tox_deserialize_char_array_rcpp(filename, max_dims)

  # Check error code
  check_err_code(result$ierr)

  # Shape flat vector into array
  array(result$values, dim = result$dims[1:result$ndim])
}


# ============================================================
#  3) Serialization (int / real / char)
# ============================================================

#> f42_serialize_int:serialize_int_nd_C: Serialize integer array to binary file
tox_serialize_int_array <- function(arr, filename) {
  #validate inputs
  validate_array_or_vector(arr)
  validate_filename(filename)
  
  # Coerce to base types
  filename <- as.character(filename)
  arr <- as.array(arr)


  result <- tox_serialize_int_array_rcpp(arr, filename)

  check_err_code(result$ierr)
  invisible(NULL)
}


#> f42_serialize_real:serialize_real_nd_C: Serialize real array to binary file
tox_serialize_real_array <- function(arr, filename) {
  #validate inputs
  validate_array_or_vector(arr)
  validate_filename(filename)
  

  # Coerce to base types
    filename <- as.character(filename)
    arr <- as.array(arr)

  # Call Rcpp wrapper
  result <- tox_serialize_real_array_rcpp(arr, filename)

  check_err_code(result$ierr)
  invisible(NULL)
}


#> f42_serialize_char:serialize_char_nd_C: Serialize character array to binary file
tox_serialize_char_array <- function(arr, filename) {
  #validate inputs
  validate_array_or_vector(arr)
  validate_filename(filename)

 
  # Coerce to base types
  filename <- as.character(filename)
  arr <- as.array(arr)

   
  # Call Rcpp wrapper
  result <- tox_serialize_char_array_rcpp(arr, filename)

  check_err_code(result$ierr)
  invisible(NULL)
}

# ============================================================
#  Deserialization (logical / complex)
# ============================================================

#> f42_deserialize_logical:deserialize_logical_nd_C: Deserialize logical array from binary file
tox_deserialize_logical_array <- function(filename, max_dims = 5L) {
  # validate inputs
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)

  meta <- tox_get_array_metadata(filename, max_dims, with_clen = FALSE)
  
  # Coerce to base types
  filename <- as.character(filename)
  max_dims <- as.integer(max_dims)

  # Call Rcpp wrapper
  result <- tox_deserialize_logical_array_rcpp(filename, max_dims)

  # Check error code
  check_err_code(result$ierr)

  # Shape flat vector into array
  array(result$values, dim = result$dims[1:result$ndim])
}

#> f42_deserialize_complex:deserialize_complex_nd_C: Deserialize complex array from binary file
tox_deserialize_complex_array <- function(filename, max_dims = 5L) {
  # validate inputs
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)

  meta <- tox_get_array_metadata(filename, max_dims, with_clen = FALSE)
  
  # Coerce to base types
  filename <- as.character(filename)
  max_dims <- as.integer(max_dims)

  # Call Rcpp wrapper
  result <- tox_deserialize_complex_array_rcpp(filename, max_dims)

  # Check error code
  check_err_code(result$ierr)

  # Shape flat vector into array
  array(result$values, dim = result$dims[1:result$ndim])
}

# ============================================================
#  Serialization (logical / complex)
# ============================================================

#> f42_serialize_logical:serialize_logical_nd_C: Serialize logical array to binary file
tox_serialize_logical_array <- function(arr, filename) {
  # validate inputs
  validate_array_or_vector(arr)
  validate_filename(filename)
  
  # Coerce to base types
  filename <- as.character(filename)
  arr <- as.array(arr)

  # Call Rcpp wrapper
  result <- tox_serialize_logical_array_rcpp(arr, filename)

  check_err_code(result$ierr)
  invisible(NULL)
}

#> f42_serialize_complex:serialize_complex_nd_C: Serialize complex array to binary file
tox_serialize_complex_array <- function(arr, filename) {
  # validate inputs
  validate_array_or_vector(arr)
  validate_filename(filename)
  
  # Coerce to base types
  filename <- as.character(filename)
  arr <- as.array(arr)

  # Call Rcpp wrapper
  result <- tox_serialize_complex_array_rcpp(arr, filename)

  check_err_code(result$ierr)
  invisible(NULL)
}

# ============================================================
#  Trajectory Contribution Analysis (TCA) Functions
# ============================================================

#> tox_trajectory_contribution_analysis:compute_contributions_c: Compute local and total contributions for factor-dependent pairs
#' Compute local and total contributions for factor-dependent pairs
#'
#' Computes contribution metrics for a single factor-dependent timeseries pair
#' using the specified baseline mode.
#'
#' @param factor NumericVector of factor values across timepoints
#' @param dependent NumericVector of dependent values across timepoints (same length as factor)
#' @param mode Character string specifying baseline mode: "raw", "min", or "mean"
#'
#' @return List with components:
#'   - local_contributions: NumericVector of per-timepoint contributions
#'   - total_contribution: Numeric scalar (sum of contributions)
#'   - ierr: Integer error code (0 = success)
#'
tox_compute_contributions <- function(factor, dependent, mode = "raw") {
  # Input validation
  validate_numeric_vector(factor)
  validate_numeric_vector(dependent)
  validate_same_length(factor, dependent)
  validate_nonempty(factor)
  
  if (!mode %in% c("raw", "min", "mean")) {
    stop("mode must be one of: 'raw', 'min', 'mean'")
  }
  
  # Coerce to base types
  factor <- as.numeric(factor)
  dependent <- as.numeric(dependent)
  mode <- as.character(mode)
  
  # Call Rcpp wrapper
  result <- tox_compute_contributions_rcpp(factor, dependent, mode)
  
  check_err_code(result$ierr)
  return(result)
}

#> tox_trajectory_contribution_analysis:compute_all_contributions_c: Compute contributions for all factor-dependent pairs in trajectories
#' Compute contributions for all factor-dependent pairs in trajectories
#'
#' Computes contributions across all selected factor-dependent pairs in 3D trajectory data
#' (n_factors × n_samples × n_timepoints).
#'
#' @param trajectories NumericVector of flattened 3D trajectory data
#' @param dims IntegerVector of length 3: (n_factors, n_samples, n_timepoints)
#' @param factor_indices IntegerVector of 1-based factor indices to analyze
#' @param dependent_indices IntegerVector of 1-based dependent indices to analyze
#' @param mode Character string specifying baseline mode: "raw", "min", or "mean"
#'
#' @return List with components:
#'   - local_contributions: 4D array of per-timepoint contributions
#'   - total_contributions: 3D array of total contributions per factor-dependent-sample
#'   - temp_factors: NumericMatrix of extracted factor trajectories
#'   - temp_dependent: NumericVector of extracted dependent trajectory
#'   - ierr: Integer error code (0 = success)
#'
tox_compute_all_contributions <- function(trajectories, dims, factor_indices, dependent_indices, mode = "raw") {
  # Input validation
  validate_numeric_vector(trajectories)
  validate_integer_vector(dims)
  validate_integer_vector(factor_indices)
  validate_integer_vector(dependent_indices)
  
  if (!mode %in% c("raw", "min", "mean")) {
    stop("mode must be one of: 'raw', 'min', 'mean'")
  }
  
  # Coerce to base types
  trajectories <- as.numeric(trajectories)
  dims <- as.integer(dims)
  factor_indices <- as.integer(factor_indices)
  dependent_indices <- as.integer(dependent_indices)
  mode <- as.character(mode)
  
  # Call Rcpp wrapper
  result <- tox_compute_all_contributions_rcpp(trajectories, dims, factor_indices, dependent_indices, mode)
  
  check_err_code(result$ierr)
  return(result)
}

#> tox_trajectory_contribution_analysis:compute_baselines_factor_dependent_c: Compute baseline values for factor and dependent variables
#' Compute baseline values for factor and dependent variables
#'
#' Computes baseline values for a factor-dependent pair using the specified mode
#' (raw = no baseline, min = minimum value, mean = mean value).
#'
#' @param factor NumericVector of factor values across timepoints
#' @param dependent NumericVector of dependent values across timepoints (same length as factor)
#' @param mode Character string specifying baseline mode: "raw", "min", or "mean"
#'
#' @return List with components:
#'   - factor_baseline: Numeric scalar baseline for factor
#'   - dependent_baseline: Numeric scalar baseline for dependent
#'   - ierr: Integer error code (0 = success)
#'
tox_compute_baselines_factor_dependent <- function(factor, dependent, mode = "raw") {
  # Input validation
  validate_numeric_vector(factor)
  validate_numeric_vector(dependent)
  validate_same_length(factor, dependent)
  validate_nonempty(factor)
  
  if (!mode %in% c("raw", "min", "mean")) {
    stop("mode must be one of: 'raw', 'min', 'mean'")
  }
  
  # Coerce to base types
  factor <- as.numeric(factor)
  dependent <- as.numeric(dependent)
  mode <- as.character(mode)
  
  # Call Rcpp wrapper
  result <- tox_compute_baselines_factor_dependent_rcpp(factor, dependent, mode)
  
  check_err_code(result$ierr)
  return(result)
}

#> tox_trajectory_contribution_analysis:perform_permutation_test_c: Perform permutation test for contribution significance
#' Perform permutation test for contribution significance
#'
#' Tests the significance of a contribution by comparing observed contributions
#' to a permutation distribution.
#'
#' @param trajectories NumericVector of flattened 3D trajectory data
#' @param dims IntegerVector of length 3: (n_factors, n_samples, n_timepoints)
#' @param factor_idx Integer 1-based index of the factor
#' @param dependent_idx Integer 1-based index of the dependent variable
#' @param sample_idx Integer 1-based index of the sample
#' @param mode Character string specifying baseline mode: "raw", "min", or "mean"
#' @param n_permutations Integer number of permutations to perform
#' @param random_seed Integer random seed for reproducibility
#'
#' @return List with components:
#'   - local_contributions: 2D array of permutation local contributions
#'   - total_contributions: NumericVector of permutation total contributions
#'   - temp_factor: NumericVector of extracted factor trajectory
#'   - temp_dependent: NumericVector of extracted dependent trajectory
#'   - ierr: Integer error code (0 = success)
#'
tox_perform_permutation_test <- function(trajectories, dims, factor_idx, dependent_idx, sample_idx,
                                         mode = "raw", n_permutations = 1000, random_seed = 12345) {
  # Input validation
  validate_numeric_vector(trajectories)
  validate_nonempty(trajectories)
  validate_integer_vector(dims)
  validate_positive_integer_scalar(factor_idx, "factor_idx")
  validate_positive_integer_scalar(dependent_idx, "dependent_idx")
  validate_positive_integer_scalar(sample_idx, "sample_idx")
  validate_positive_integer_scalar(n_permutations, "n_permutations")
  validate_positive_integer_scalar(random_seed, "random_seed")
  
  if (!mode %in% c("raw", "min", "mean")) {
    stop("mode must be one of: 'raw', 'min', 'mean'")
  }
  
  # Coerce to base types
  trajectories <- as.numeric(trajectories)
  dims <- as.integer(dims)
  factor_idx <- as.integer(factor_idx)
  dependent_idx <- as.integer(dependent_idx)
  sample_idx <- as.integer(sample_idx)
  mode <- as.character(mode)
  n_permutations <- as.integer(n_permutations)
  random_seed <- as.integer(random_seed)
  
  # Call Rcpp wrapper
  result <- tox_perform_permutation_test_rcpp(trajectories, dims, factor_idx, dependent_idx, sample_idx, 
                                              mode, n_permutations, random_seed)
  
  check_err_code(result$ierr)
  return(result)
}

#> tox_trajectory_contribution_analysis:compute_p_values_c: Compute p-values for contribution significance
#' Compute p-values for contribution significance
#'
#' Computes p-values by comparing observed contributions to permutation distributions.
#' P-values represent the fraction of permutations where the contribution was >= observed.
#'
#' @param local_contributions_observed NumericVector of observed per-timepoint contributions (length n_timepoints)
#' @param total_contribution_observed Numeric scalar of observed total contribution
#' @param local_contributions_perm NumericVector of flattened permutation local contributions (length n_timepoints * n_permutations)
#' @param total_contributions_perm NumericVector of permutation total contributions (length n_permutations)
#' @param dims IntegerVector of length 2: (n_timepoints, n_permutations)
#'
#' @return List with components:
#'   - local_p_values: NumericVector of p-values for each timepoint
#'   - total_p_value: Numeric scalar p-value for total contribution
#'   - ierr: Integer error code (0 = success)
#'
tox_compute_p_values <- function(local_contributions_observed, total_contribution_observed, 
                                  local_contributions_perm, total_contributions_perm, dims) {
  # Input validation
  validate_numeric_vector(local_contributions_observed)
  validate_nonempty(local_contributions_observed)
  validate_positive_numeric_scalar(total_contribution_observed, "total_contribution_observed")
  validate_numeric_vector(local_contributions_perm)
  validate_numeric_vector(total_contributions_perm)
  validate_integer_vector(dims)
  
  if (length(dims) != 2) {
    stop("dims must have length 2: (n_timepoints, n_permutations)")
  }
  
  # Call Rcpp wrapper with type conversion
  result <- tox_compute_p_values_rcpp(
    as.numeric(local_contributions_observed),
    as.numeric(total_contribution_observed),
    as.numeric(local_contributions_perm),
    as.numeric(total_contributions_perm),
    as.integer(dims)
  )
  
  check_err_code(result$ierr)
  return(result)
}

# ============================================================
#  Trajectory Normalization Functions
# ============================================================

#> tox_trajectory_normalization:normalize_variable_timeseries_c: Normalize a variable timeseries using min-max scaling
#' Normalize a variable timeseries using min-max scaling
#'
#' Performs min-max normalization (scaling to [0,1] range) on a single variable timeseries.
#'
#' @param v NumericVector of timeseries values to normalize
#'
#' @return List with components:
#'   - v_norm: NumericVector of normalized values in [0,1] range
#'   - ierr: Integer error code (0 = success)
#'   - status: Integer status code (additional warnings/info)
#'
tox_normalize_variable_timeseries <- function(v) {
  # Input validation
  validate_numeric_vector(v)
  validate_nonempty(v)
  
  # Coerce to base types
  v <- as.numeric(v)
  
  # Call Rcpp wrapper
  result <- tox_normalize_variable_timeseries_rcpp(v)
  
  check_err_code(result$ierr)
  return(result)
}

#> tox_trajectory_normalization:normalize_single_trajectory_c: Normalize a single trajectory (multiple factors across timepoints)
#' Normalize a single trajectory (multiple factors across timepoints)
#'
#' Performs min-max normalization independently for each factor in a trajectory
#' (n_timepoints × n_factors matrix).
#'
#' @param trajectory NumericMatrix with dimensions (n_timepoints × n_factors)
#'   Each column represents a single factor normalized across timepoints
#'
#' @return List with components:
#'   - trajectory_norm: NumericMatrix with same dimensions as input, normalized per factor
#'   - ierr: Integer error code (0 = success)
#'   - status: Integer status code (additional warnings/info)
#'
tox_normalize_single_trajectory <- function(trajectory) {
  # Input validation
  validate_numeric_matrix(trajectory, "trajectory")
  
  if (nrow(trajectory) == 0 || ncol(trajectory) == 0) {
    stop("trajectory must have non-zero dimensions")
  }
  
  # Coerce to base types
  trajectory <- as.matrix(trajectory)
  mode(trajectory) <- "numeric"
  
  # Call Rcpp wrapper
  result <- tox_normalize_single_trajectory_rcpp(trajectory)
  
  check_err_code(result$ierr)
  return(result)
}

#> tox_trajectory_normalization:normalize_all_trajectories_c: Normalize all trajectories in 3D trajectory data
#' Normalize all trajectories in 3D trajectory data
#'
#' Performs min-max normalization independently for each factor across all samples and timepoints.
#' Input: 3D array of dimensions (n_factors × n_samples × n_timepoints)
#'
#' @param trajectories NumericVector of flattened 3D trajectory data
#' @param dims IntegerVector of length 3: (n_factors, n_samples, n_timepoints)
#'
#' @return List with components:
#'   - trajectories_norm: NumericVector of normalized flattened data with 3D dimensions restored
#'   - ierr: Integer error code (0 = success)
#'   - status: Integer status code (additional warnings/info)
#'
tox_normalize_all_trajectories <- function(trajectories, dims) {
  # Input validation
  validate_numeric_vector(trajectories)
  validate_nonempty(trajectories)
  validate_integer_vector(dims)
  
  # Coerce to base types
  trajectories <- as.numeric(trajectories)
  dims <- as.integer(dims)
  
  # Call Rcpp wrapper
  result <- tox_normalize_all_trajectories_rcpp(trajectories, dims)
  
  check_err_code(result$ierr)
  return(result)
}



# ============================================================
#  4) KD-tree index (multidimensional) + spherical KD
# ============================================================

#> f42_kd_tree:build_kd_index_C: Build KD-tree index for multidimensional points
build_kd_index <- function(X, dim_order = NULL) {
  # Input validation using standardized validation functions
  validate_numeric_matrix(X, "X")
  
  d <- as.integer(nrow(X))
  n <- as.integer(ncol(X))
  
  if (is.null(dim_order)) {
    dim_order <- 1:d  # Default: use dimensions in order
  }
  
  # Convert and validate dim_order
  dim_order <- as.integer(dim_order)
  
  # Ensure X is a numeric (double) matrix with expected dimensions
  X <- matrix(as.numeric(X), nrow = d, ncol = n)
  
  # Call low-level Rcpp wrapper with explicit dimensions (num_dimensions, num_points)
  result <- tox_build_kd_index_rcpp(X, dim_order)
  
  # Check backend error code
  check_err_code(result$ierr)
  
  # Return KD index vector
  result$kd_ix
}

#> f42_kd_tree:build_spherical_kd_index_C: Build Spherical KD-Tree index (R-level wrapper)
#' Build Spherical KD-Tree index (R-level wrapper)
build_spherical_kd <- function(V, dim_order = NULL) {
  # R-layer validation
  validate_numeric_matrix(V, "V")
  
  d <- as.integer(nrow(V))
  n <- as.integer(ncol(V))
  
  if (is.null(dim_order)) {
    dim_order <- 1:d  # Default: use dimensions in order
  }
  
  # Validate dim_order
  validate_logical_or_index_vector(dim_order, NULL, "dim_order")
  dim_order <- as.integer(dim_order)
  
  # Ensure V is a numeric (double) matrix with expected dimensions
  V <- matrix(as.numeric(V), nrow = d, ncol = n)
  
  # Call Rcpp wrapper (which wraps the Fortran call)
  # Underlying Fortran computes and returns its own `dimension_order`.
  # Do not pass `dim_order` to the low-level routine; it returns its own ordering.
  result <- tox_build_spherical_kd_rcpp(V,dim_order)
  
  # Check backend error code
  check_err_code(result$ierr)
  
  # Return sphere index vector
  result$sphere_ix
}
#' Get a point from the KD index (R-level helper)
get_kd_point <- function(X, kd_ix, position) {
  # Input validation using standardized validation functions
  validate_numeric_matrix(X, "X")
  validate_array_or_vector(kd_ix, "kd_ix")
  
  # Extract the point at the specified position
  # Note: X is organized as (dimensions, points), so we need to extract column kd_ix[position]
  X[, kd_ix[position], drop = FALSE]
}
  



# ============================================================
#  5) BST index + range query
# ============================================================

#> f42_binary_search_tree:build_bst_index_C: Build binary search tree index
#  BST index
build_bst_index <- function(x) {
  validate_numeric_vector(x, "x")

  # Call Rcpp wrapper with type conversion
  result <- build_bst_index_rcpp(as.numeric(x))

  # Return index vector
  return(result)
}
 
#' Get sorted value by BST index (R helper)
get_sorted_value <- function(x, ix, position) {
  validate_numeric_vector(x, "x")
  validate_array_or_vector(ix, "ix")
  
  # Cast once and use
  ix_int <- as.integer(ix)
  position_int <- as.integer(position)
  return(as.numeric(x)[ix_int[position_int]])
}

# ============================================================
# BST range query

#> f42_binary_search_tree:bst_range_query_C: Query BST index for values in range [lo, hi]
bst_range_query <- function(x, ix, lo, hi) {

  validate_numeric_vector(x, "x")
  validate_integer_vector(ix, "ix")
  validate_equal_length(x, ix, "x", "ix")

  # Call Rcpp wrapper with type conversion
  result <- bst_range_query_rcpp(
    as.numeric(x),
    as.integer(ix),
    as.numeric(lo),
    as.numeric(hi)
  )

  # Check for errors
  check_err_code(result$ierr)

  # Keep the same public return structure as old version
  list(
    indices = result$indices,
    count   = result$count
  )
}


# ============================================================
#  6) which() helper via C backend
# ============================================================

#> f42_utils:which_C: Extract indices where mask is non-zero
tox_which <- function(mask, m_max = length(mask)) {
  validate_integer_vector(mask, "mask")

  result <- tox_which_rcpp(
    mask = as.integer(mask),
    m_max = as.integer(m_max)
  )

  check_err_code(result$ierr)

  result$idx_out[seq_len(result$m_out)]
}


# ============================================================
#  7) LOESS smoothing (2D)
# ============================================================

#> tox_get_outliers:loess_smooth_2d_c: Perform 2D LOESS smoothing with kernel weighting
tox_loess_smooth_2d <- function(
  x_ref,
  y_ref,
  x_query,
  indices_used = NULL,
  kernel_sigma,
  kernel_cutoff
) {




  validate_numeric_vector(x_ref, "x_ref")
  validate_numeric_vector(x_query, "x_query")

  # y_ref can be vector or matrix; convert first, then validate
  if (!is.matrix(y_ref)) {
    y_ref <- matrix(y_ref, nrow = 1L)
  }
  validate_numeric_matrix(y_ref, "y_ref")
 
  # Indices used: default to all, then validate
  # total number of reference points
  n_total <- length(x_ref)
  if (is.null(indices_used)) {
    indices_used <- seq_len(n_total)
  }
  validate_index_vector(indices_used, n_total, "indices_used")

  # Kernel parameters
  validate_positive_numeric_scalar(kernel_sigma, "kernel_sigma")
  validate_positive_numeric_scalar(kernel_cutoff, "kernel_cutoff")

  ## 2) Convert to types suitable for Rcpp / Fortran

  x_ref        <- as.numeric(x_ref)
  x_query      <- as.numeric(x_query)
  y_ref        <- matrix(as.numeric(y_ref), nrow = 1L)
  indices_used <- as.integer(indices_used)
  kernel_sigma <- as.numeric(kernel_sigma)
  kernel_cutoff<- as.numeric(kernel_cutoff)

  ## 3) Call Rcpp wrapper

  result <- tox_loess_smooth_2d_rcpp(
    n_total      = as.integer(n_total),
    n_target     = as.integer(length(x_query)),
    x_ref        = x_ref,
    y_ref        = y_ref,
    indices_used = indices_used,
    n_used       = as.integer(length(indices_used)),
    x_query      = x_query,
    kernel_sigma = kernel_sigma,
    kernel_cutoff= kernel_cutoff
  )

  ## 4) Error handling

  check_err_code(result$ierr)

  ## 5) Structured return

  y_out <- result$y_out
  list(
    y_out           = y_out,
    smoothed_values = as.vector(y_out)
  )
}

#> tox_relative_axis_plane_tools:omics_vector_RAP_projection_c: Project selected vectors onto a Relative Axis Plane (RAP)
#' Project selected vectors onto a Relative Axis Plane (RAP)
#' 
#' This function validates inputs and calls the C/Rcpp RAP projection wrapper.
#' @param vecs Numeric matrix (n_axes x n_vecs)
#' @param vecs_selection_mask Logical or integer vector (length n_vecs)
#' @param axes_selection_mask Logical or integer vector (length n_axes)
#' @return Matrix of projections (n_selected_axes x n_selected_vecs)
#' @export
tox_omics_vector_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
  validate_numeric_matrix(vecs, "vecs")
  n_axes <- nrow(vecs)
  n_vecs <- ncol(vecs)
  validate_logical_or_index_vector(vecs_selection_mask, expected_length = n_vecs, name = "vecs_selection_mask")
  validate_logical_or_index_vector(axes_selection_mask, expected_length = n_axes, name = "axes_selection_mask")

  res <- tox_omics_vector_RAP_projection_rcpp(
    vecs,
    as.integer(vecs_selection_mask),
    as.integer(axes_selection_mask)
  )
  check_err_code(res$ierr)
  return(res$projections)
}

#> tox_relative_axis_plane_tools:omics_field_RAP_projection_c: Project selected vector fields onto a Relative Axis Plane (RAP)
#' Project selected vector fields onto a Relative Axis Plane (RAP)
#' 
#' This function validates inputs and calls the C/Rcpp RAP field projection wrapper.
#' @param vecs Numeric matrix (2*n_axes x n_vecs)
#' @param vecs_selection_mask Logical or integer vector (length n_vecs)
#' @param axes_selection_mask Logical or integer vector (length n_axes)
#' @return Matrix of projections (n_selected_axes x n_selected_vecs)
#' @export
tox_omics_field_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
  validate_numeric_matrix(vecs, "vecs")
  n_axes <- nrow(vecs) / 2
  n_vecs <- ncol(vecs)
  validate_logical_or_index_vector(vecs_selection_mask, expected_length = n_vecs, name = "vecs_selection_mask")
  validate_logical_or_index_vector(axes_selection_mask, expected_length = n_axes, name = "axes_selection_mask")

  res <- tox_omics_field_RAP_projection_rcpp(
    vecs,
    as.integer(vecs_selection_mask),
    as.integer(axes_selection_mask)
  )
  check_err_code(res$ierr)
  return(res$projections)
}

#> tox_relative_axis_plane_tools:relative_axes_changes_from_shift_vector_c: Compute relative axis changes from a RAP-projected and normalized shift vector
#' Compute relative axis changes from a RAP-projected and normalized shift vector
#' 
#' @param vec Numeric vector (RAP-projected and normalized shift vector)
#' @return Numeric vector of fractional axis contributions (sums to 1)
#' @export
tox_relative_axes_changes_from_shift_vector <- function(vec) {
  validate_numeric_vector(vec, "vec")
  res <- tox_relative_axes_changes_from_shift_vector_rcpp(vec)
  return(res)
}

#> tox_relative_axis_plane_tools:relative_axes_expression_from_expression_vector_c: Compute relative axis contributions from a RAP-projected and normalized expression vector
#' Compute relative axis contributions from a RAP-projected and normalized expression vector
#' 
#' @param vec Numeric vector (RAP-projected and normalized expression vector)
#' @return Numeric vector of fractional axis contributions (sums to 1)
#' @export
tox_relative_axes_expression_from_expression_vector <- function(vec) {
  validate_numeric_vector(vec, "vec")
  res <- tox_relative_axes_expression_from_expression_vector_rcpp(vec)
  return(res)
}

#> tox_relative_axis_plane_tools:clock_hand_angle_between_vectors_c: Compute signed clock hand angle between two RAP-projected and normalized vectors
#' Compute signed clock hand angle between two RAP-projected and normalized vectors
#' 
#' @param v1 Numeric vector (first normalized vector in RAP space)
#' @param v2 Numeric vector (second normalized vector in RAP space)
#' @param selected_axes_for_signed Integer vector of length 3 (axes for directionality, ignored if length(v1) <= 3)
#' @return Signed angle in radians [-pi, pi]
#' @export
tox_clock_hand_angle_between_vectors <- function(v1, v2, selected_axes_for_signed) {
  validate_numeric_vector(v1, "v1")
  validate_numeric_vector(v2, "v2")
  validate_same_length(v1, v2, "v1", "v2")
  validate_integer_vector(selected_axes_for_signed, "selected_axes_for_signed", expected_length = 3)
  res <- tox_clock_hand_angle_between_vectors_rcpp(v1, v2, selected_axes_for_signed)
  return(res)
}

#> tox_relative_axis_plane_tools:clock_hand_angles_for_shift_vectors_c: Compute signed rotation angles for pairs of RAP-projected and normalized vectors
#' Compute signed rotation angles for pairs of RAP-projected and normalized vectors
#' 
#' @param origins Numeric matrix (n_dims x n_vecs), first set of vectors
#' @param targets Numeric matrix (n_dims x n_vecs), second set of vectors
#' @param vecs_selection_mask Integer or logical vector (length n_vecs)
#' @param selected_axes_for_signed Integer vector of length 3 (axes for directionality)
#' @return Numeric vector of signed angles (radians)
#' @export
tox_clock_hand_angles_for_shift_vectors <- function(origins, targets, vecs_selection_mask, selected_axes_for_signed) {
  validate_numeric_matrix(origins, "origins")
  validate_numeric_matrix(targets, "targets")
  validate_same_length(vecs_selection_mask, rep(0, ncol(origins)), "vecs_selection_mask", "origins columns")
  validate_integer_vector(selected_axes_for_signed, "selected_axes_for_signed", expected_length = 3)
  res <- tox_clock_hand_angles_for_shift_vectors_rcpp(origins, targets, as.integer(vecs_selection_mask), selected_axes_for_signed)
  return(res)
}
#> tox_clustering:cluster_factor_trajectories_k_means_c: K-means clustering on factor trajectories
#' K-means clustering on factor trajectories
#'
#' Performs k-means clustering on factor trajectories (factor evolution over time).
#' All validation is performed in R. Returns centroids, labels, label_counts, and ierr.
#' @param n_clusters Number of clusters (integer)
#' @param trajectories Numeric vector (flattened, n_factors * n_samples * n_timepoints)
#' @param n_factors Number of factors (integer)
#' @param n_samples Number of samples (integer)
#' @param n_timepoints Number of timepoints (integer)
#' @param centroids Numeric matrix (n_factors x n_clusters)
#' @param max_iterations Maximum number of iterations (integer)
#' @return List with centroids, labels, label_counts, ierr
#' @export
tox_cluster_factor_trajectories_k_means <- function(n_clusters, trajectories, n_factors, n_samples, n_timepoints, centroids, max_iterations = 300) {
  validate_positive_integer_scalar(n_clusters, "n_clusters")
  validate_numeric_vector(trajectories, "trajectories")
  validate_positive_integer_scalar(n_factors, "n_factors")
  validate_positive_integer_scalar(n_samples, "n_samples")
  validate_positive_integer_scalar(n_timepoints, "n_timepoints")
  validate_numeric_matrix(centroids, "centroids")
  validate_positive_integer_scalar(max_iterations, "max_iterations")
  # Check dimensions
  if (length(trajectories) != n_factors * n_samples * n_timepoints) stop("trajectories length mismatch")
  validate_matrix_shape_factor_centroids(centroids, n_factors, n_clusters)
  res <- tox_cluster_factor_trajectories_k_means_rcpp(n_clusters, as.numeric(trajectories), n_factors, n_samples, n_timepoints, centroids, max_iterations)
  check_err_code(res$ierr)
  return(res)
}

#> tox_clustering:k_means_clustering_c: K-means clustering (general)
#' K-means clustering (general)
#'
#' Performs k-means clustering on general data points.
#' All validation is performed in R. Returns centroids, labels, label_counts, and ierr.
#' @param n_clusters Number of clusters (integer)
#' @param data_points Numeric matrix (n_dims x n_points)
#' @param n_points Number of points (integer)
#' @param n_dims Number of dimensions (integer)
#' @param centroids Numeric matrix (n_dims x n_clusters)
#' @param max_iterations Maximum number of iterations (integer)
#' @return List with centroids, labels, label_counts, ierr
#' @export
tox_k_means_clustering <- function(n_clusters, data_points, n_points, n_dims, centroids, max_iterations = 300) {
  validate_positive_integer_scalar(n_clusters, "n_clusters")
  validate_numeric_matrix(data_points, "data_points")
  validate_positive_integer_scalar(n_points, "n_points")
  validate_positive_integer_scalar(n_dims, "n_dims")
  validate_numeric_matrix(centroids, "centroids")
  validate_positive_integer_scalar(max_iterations, "max_iterations")
  validate_matrix_shape_data_points(data_points, n_dims, n_points)
  validate_matrix_shape_centroids(centroids, n_dims, n_clusters)
  res <- tox_k_means_clustering_rcpp(n_clusters, data_points, n_points, n_dims, centroids, max_iterations)
  check_err_code(res$ierr)
  return(res)
}
