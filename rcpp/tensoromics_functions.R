
#> f42_helper-import_libs: Import necessary packages
library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library

lib_path <- shQuote(normalizePath("build"))

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv, cacheDir = "rcpp/rcpp_cache")

cat("✓ TensorOmics Rcpp functions loaded successfully\n")

source("rcpp/error_handling.R")

# ===================================================================
# EUCLIDEAN DISTANCE FUNCTIONS
# ===================================================================

#> tox_euclidean_distance:euclidean_distance_c: Calculate Euclidean distance between two vectors
#' Computes the Euclidean distance between two vectors of the same dimension.
#' 
#' @param vec1 First vector (numeric)
#' @param vec2 Second vector (numeric, same length as vec1)
#' 
#' @return Numeric value representing the Euclidean distance between the vectors
#' 
tox_euclidean_distance <- function(vec1, vec2) {
  # Input Validation
  validate_numeric_vector(vec1)
  validate_numeric_vector(vec2)
  validate_same_length(vec1, vec2)
  validate_nonempty(vec1)

  # Call the Rcpp forwarder
  result <- tox_euclidean_distance_rcpp(vec1, vec2)

  # Return distance
  return(result)
}

#> tox_euclidean_distance:distance_to_centroid_c: Calculate distance from each gene to its family centroid
#' Calculate distances from genes to their family centroids
#' 
#' Computes the Euclidean distance from each gene to its corresponding family centroid.
#' 
#' @param genes Matrix of gene expression data (genes as columns, dimensions as rows)
#' @param centroids Matrix of family centroids (families as columns, dimensions as rows) 
#' @param gene_to_fam Integer vector mapping each gene to its family index (1-based)
#' @param d Integer number of dimensions
#' 
#' @return Numeric vector of distances from each gene to its family centroid
#' 
tox_distance_to_centroid <- function(genes, centroids, gene_to_fam, d) {
  # Input Validation
  validate_numeric_vector(genes)
  validate_numeric_vector(centroids)
  validate_positive_integer_scalar(d)

  # Validate flattened lengths are compatible with d
  validate_divisible_length(genes, d)
  validate_divisible_length(centroids, d)

  # Calculate dimensions
  n_genes <- as.integer(length(genes) / d)
  n_families <- as.integer(length(centroids) / d)
  validate_gene_to_family_centroid(gene_to_fam, n_genes, n_families)
  validate_length_equals_n(gene_to_fam, n_genes)

  # Call the Rcpp forwarder
  result <- tox_distance_to_centroid_rcpp(genes,
                                          centroids,
                                          gene_to_fam,
                                          d
  )

  # Return distance vector
  return(result)
}


# ===================================================================
# TISSUE VERSATILITY FUNCTIONS
# ===================================================================

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
#' @return A list with:
#' \describe{
#'   \item{tissue_versatilities}{Normalized tissue versatility values in [0, 1] for selected vectors}
#'   \item{tissue_angles_deg}{Angles in degrees [0, 90] for selected vectors}
#'   \item{n_selected_vectors}{Integer number of vectors processed}
#'   \item{n_selected_axes}{Integer number of axes used in calculation}
#'   \item{ierr}{Integer status code from backend routine}
#' }

tox_calculate_tissue_versatility <- function(expression_vectors, vector_selection, axis_selection) {
  # Input Validation
  validate_numeric_matrix(expression_vectors)
  n_axes <- nrow(expression_vectors)
  n_vectors <- ncol(expression_vectors)

  # Ensure selectors have expected lengths
  validate_logical_or_index_vector(vector_selection, expected_length = n_vectors, name = "vector_selection")
  validate_logical_or_index_vector(axis_selection, expected_length = n_axes, name = "axis_selection")

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

  # Call the Rcpp forwarder
  result <- tox_calculate_tissue_versatility_rcpp(expression_vectors,
                                                  vector_selection,
                                                  axis_selection
  )
  
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
#' @return A list with:
#' \describe{
#'   \item{is_outlier}{Logical vector indicating outliers}
#'   \item{loess_x}{Numeric vector of family median distances}
#'   \item{loess_y}{Numeric vector of family standard deviations}
#'   \item{loess_n}{Integer vector of number of genes used per family}
#'   \item{ierr}{Integer status code from backend routine}
#' }
tox_detect_outliers <- function(distances, gene_to_fam, n_families, percentile = 95.0) {

  # Input Validation
  validate_numeric_vector(distances)
  n_genes <- length(distances)
  validate_logical_or_index_vector(gene_to_fam, expected_length = n_genes, name = "gene_to_fam")
  
  
  validate_length_equals_n(gene_to_fam, n_genes)
  validate_positive_integer_scalar(n_families)
  validate_numeric_scalar(percentile)
  validate_gene_to_family_outliers(gene_to_fam, n_genes, n_families)

  # Call the Rcpp forwarder
  result <- tox_detect_outliers_rcpp(distances,
                                     gene_to_fam,
                                     n_families,
                                     percentile
  )

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
#' @return A list with:
#' \describe{
#'   \item{dscale}{Numeric vector of scaling factors for each family}
#'   \item{loess_x}{Numeric vector of family median distances}
#'   \item{loess_y}{Numeric vector of family standard deviations}
#'   \item{indices_used}{Integer vector of number of genes used per family}
#'   \item{ierr}{Integer status code from backend routine}
#' }
tox_compute_family_scaling <- function(distances, gene_to_fam, n_families) {

  # Input Validation
  validate_numeric_vector(distances)
  n_genes <- length(distances)
  validate_length_equals_n(gene_to_fam, n_genes)

  # Call the Rcpp forwarder
  result <- tox_compute_family_scaling_rcpp(distances,
                                            gene_to_fam,
                                            n_families
  )

  # Check for errors
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
#' @return A list with:
#' \describe{
#'   \item{dscale}{Numeric vector of scaling factors for each family}
#'   \item{loess_x}{Numeric vector of family median distances}
#'   \item{loess_y}{Numeric vector of family standard deviations}
#'   \item{indices_used}{Integer vector of number of genes used per family}
#'   \item{perm_tmp}{Final state of permutation array}
#'   \item{stack_left_tmp}{Final state of left stack array}
#'   \item{stack_right_tmp}{Final state of right stack array}
#'   \item{family_distances}{Final state of family distances array}
tox_compute_family_scaling_expert <- function(distances, gene_to_fam, n_families,
                                              perm_tmp, stack_left_tmp, stack_right_tmp,
                                              family_distances) {
  # Input Validation
  validate_numeric_vector(distances)
  n_genes <- length(distances)
  validate_length_equals_n(gene_to_fam, n_genes)
  validate_length_equals_n(perm_tmp, n_genes)
  validate_length_equals_n(stack_left_tmp, n_genes)
  validate_length_equals_n(stack_right_tmp, n_genes)
  validate_length_equals_n(family_distances, n_genes)


  # Call the Rcpp forwarder
  result <- tox_compute_family_scaling_expert_rcpp(n_families,
                                                   distances,
                                                   gene_to_fam,
                                                   perm_tmp,
                                                   stack_left_tmp,
                                                   stack_right_tmp,
                                                   family_distances
  )

  # Check for errors
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
#' @return A list with:
#' \describe{
#'   \item{rdi}{Numeric vector of Relative Distance Index values for each gene}
#'   \item{sorted_rdi}{Numeric vector of RDI values sorted in ascending order}
#' }
tox_compute_rdi <- function(distances, gene_to_fam, dscale) {

  # Input Validation
  validate_numeric_vector(distances)
  n_genes <- length(distances)
  validate_length_equals_n(gene_to_fam, n_genes)
  n_families <- length(dscale)
  validate_index_bounds(gene_to_fam, low = 1, high = n_families)

  # Call the Rcpp forwarder
  result <- tox_compute_rdi_rcpp(distances,
                                 gene_to_fam,
                                 dscale
  )

  # Return structured result
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
#' @return A list with:
#' \describe{
#'   \item{is_outlier}{Logical vector indicating outliers}
#'   \item{threshold}{Numeric RDI threshold value used}
#' }
tox_identify_outliers <- function(rdi, percentile = 95.0) {

  # Input Validation
  validate_numeric_vector(rdi)
  validate_numeric_scalar(percentile)

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
#'
#' 
#' This function wraps the Fortran subroutine `normalize_by_std_dev`
#' to normalize each gene's expression across tissues by the standard deviation.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A normalized numeric matrix with the same dimensions and names as the input.
 
tox_normalize_by_std_dev <- function(input_matrix) {

  # Input Validation
  validate_numeric_matrix_values(input_matrix)

  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)
  
  # Call the Rcpp forwarder
  result <- tox_normalize_by_std_dev_rcpp(input_matrix)

  # Check for errors
  check_err_code(result$ierr)

  # Return the normalized matrix with original dimensions and names
  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
}

#> tox_normalization:quantile_normalization_c: Quantile normalization of gene expression values
#' Quantile normalization of gene expression values
#'
#' This function wraps the Fortran subroutine `quantile_normalization`
#' to apply quantile normalization across all tissue columns.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A quantile-normalized numeric matrix with the same dimensions and names as the input.
#'
#' 
tox_quantile_normalization <- function(input_matrix) {

  # Input Validation
  validate_matrix(input_matrix)

  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)

  # Call the Rcpp forwarder
  result <- tox_quantile_normalization_rcpp(input_matrix)
 
  # Check for errors
  check_err_code(result$ierr)
  
  # Return the normalized matrix with original dimensions and names
  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
}

#> tox_normalization:log2_transformation_c: Apply log2(x + 1) transformation to gene expression values
#' Apply log2(x + 1) transformation to gene expression values
#'
#' This function wraps the Fortran subroutine `log2_transformation`
#' to apply a log2(x + 1) transformation to each element in the input matrix.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A numeric matrix with log2-transformed expression values, with the same dimensions and names as the input.
#'
tox_log2_transformation <- function(input_matrix) {

  # Input Validation
  validate_matrix(input_matrix)

  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)

  # Call the Rcpp forwarder
  result <- tox_log2_transformation_rcpp(input_matrix)
  
  # Check for errors 
  check_err_code(result$ierr)
  
  # Return the log2-transformed matrix with original dimensions and names
  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
}

#> tox_normalization:normalize_unit_length_c: Normalize a vector to unit length
#' Normalize a numeric vector to unit length
#'
#' @param vector Numeric vector
#' @return Numeric vector with unit norm
#' 
tox_normalize_unit_length <- function(vector) {
  
  # Input Validation
  validate_numeric_vector(vector)
  validate_nonempty(vector)

  # Call the Rcpp forwarder
  result <- tox_normalize_unit_length_rcpp(vector)
  
  # Check for errors
  check_err_code(result$ierr)

  # Return the normalized vector
  return(result$vector)
}

#> tox_normalization:calc_tiss_avg_c: Calculate average expression across replicates for each tissue group
#' Calculate average expression across replicates for each tissue group
#'
#' This function wraps the Fortran subroutine `calc_tiss_avg`
#' to compute the mean expression value for replicates grouped by tissue.
#'
#' @param df A data frame or matrix with genes as rows and tissue replicates as columns.
#' @return A data frame with genes as rows and averaged tissues as columns.
#' 
tox_calculate_tissue_averages <- function(df) {
  
  # Input Validation
  validate_matrix(df)

  tissue_groups <- vapply(colnames(df), tox_parse_tissue_group, FUN.VALUE = character(1), USE.NAMES = FALSE)
  unique_groups <- unique(tissue_groups)

  n_groups <- length(unique_groups)
  group_starts <- integer(n_groups)
  group_counts <- integer(n_groups)
  df_sorted <- df[, order(tissue_groups)]
  sorted_tissue_groups <- tissue_groups[order(tissue_groups)]
  current_group <- sorted_tissue_groups[1]
  group_starts[1] <- 1
  group_counts[1] <- 1
  group_idx <- 1
  for (i in 2:length(sorted_tissue_groups)) {
    if (!is.na(sorted_tissue_groups[i]) && !is.na(current_group) && sorted_tissue_groups[i] == current_group) {
      group_counts[group_idx] <- group_counts[group_idx] + 1
    } else {
      group_idx <- group_idx + 1
      group_starts[group_idx] <- i
      group_counts[group_idx] <- 1
      current_group <- sorted_tissue_groups[i]
    }
  }
  # Call the Rcpp forwarder
  result <- tox_calc_tiss_avg_rcpp(df_sorted,
                                   group_starts,
                                   group_counts
  )

  # Check for errors
  check_err_code(result$ierr)

  n_genes <- nrow(df_sorted)
  n_cols <- length(result$output_vector) / n_genes
  output_matrix <- matrix(result$output_vector, nrow = n_genes, ncol = n_cols)

  colnames(output_matrix) <- unique_groups
  rownames(output_matrix) <- rownames(df_sorted)

  # Return the averaged matrix as a data frame
  return(as.data.frame(output_matrix))
}


#> tox_normalization:calc_fchange_c: Calculate log2 fold changes between control and condition columns
#' Calculate log2 fold changes based on control and condition patterns
#' @param df A data frame with genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern to detect control columns.
#' @param condition_patterns A character vector with patterns to detect condition columns.
#' @return A data frame with log2 fold changes for each gene and condition, with genes as rows and conditions as columns.
tox_calculate_fold_changes <- function(df, control_pattern, condition_patterns) {
  # Input Validation
  
  validate_matrix(df)
  validate_string_scalar(control_pattern)
  validate_character_vector(condition_patterns) 
 

  # --- Identify control and condition columns ---
  indices_info <- tox_prepare_indices_by_patterns(df, control_pattern, condition_patterns)
  control_cols <- indices_info$control_cols
  condition_cols <- indices_info$condition_cols
  condition_labels <- indices_info$condition_labels

  n_pairs <- length(control_cols)
  
  # Call the Rcpp forwarder
  result <- tox_calc_fchange_rcpp(df,
                                  control_cols,
                                  condition_cols
  )
    
  # Check for errors
  check_err_code(result$ierr)
  
  n_genes <- nrow(df)
  output_matrix <- matrix(result$output_vector, nrow = n_genes, ncol = n_pairs)

  colnames(output_matrix) <- condition_labels
  rownames(output_matrix) <- rownames(df)
  
  # Return the fold change matrix as a data frame
  return(as.data.frame(output_matrix))
}


#> tox_normalization:normalization_pipeline_c: Complete normalization pipeline for gene expression data (up to log2(x+1))
#' Complete normalization pipeline for gene expression data (up to log2(x+1))
#' @param input_matrix Numeric matrix (genes x tissues)
#' @param group_s Integer vector: start column index for each replicate group (1-based)
#' @param group_c Integer vector: number of columns per replicate group
#' 
#' @return Numeric matrix with log2-transformed average expression values for each tissue group (genes x groups).
tox_normalization_pipeline <- function(input_matrix, group_s, group_c) {
  # Input Validation
  validate_matrix(input_matrix)
  
  validate_group_vectors(group_s, group_c, ncol(input_matrix))

  # Call the Rcpp forwarder
  result <- tox_normalization_pipeline_rcpp(input_matrix,
                                            group_s,
                                            group_c
  )
  
  # Check for errors
  check_err_code(result$ierr)
  
  # Return the normalized matrix with original row names and new column names based on group_s and group_c
  return(matrix(result$buf_log, nrow = nrow(input_matrix), ncol = length(group_s)))
}

# ===================================================================
# TRAJECTORY NORMALIZATION FUNCTIONS
# ===================================================================

#> tox_trajectory_normalization:normalize_variable_timeseries_c: Normalize a single variable across time using min-max scaling.
#' Normalize a single variable across time using min-max scaling
#'
#' @param v Numeric vector time series (length = n_points)
#' @return A list with:
#' \describe{
#'   \item{v_norm}{Numeric vector of normalized values in [0, 1]}
#'   \item{status}{Integer status code returned by the backend routine}
#' }
#'
tox_normalize_variable_timeseries <- function(v) {

  # Input Validation
  validate_numeric_vector(v)
  validate_nonempty(v)

  # Call the Rcpp forwarder
  result <- tox_normalize_variable_timeseries_rcpp(v)

  # Check for errors
  check_err_code(result$ierr)
  
  # Return structured result
  return(list(
    v_norm = as.numeric(result$v_norm),
    status = as.integer(result$status)
  ))
}

#> tox_trajectory_normalization:normalize_single_trajectory_c: Normalize all factors in a single trajectory independently across time.
#' Normalize all factors in a single trajectory independently across time
#'
#' @param trajectory Numeric matrix (n_timepoints x n_factors)
#' @return A list with:
#' \describe{
#'   \item{traj_norm}{Numeric matrix (n_timepoints x n_factors) with normalized values}
#'   \item{status}{Integer status code returned by the backend routine}
#' }
#'
tox_normalize_single_trajectory <- function(trajectory) {

  # Input Validation
  validate_numeric_matrix(trajectory)
  n_timepoints <- nrow(trajectory)
  n_factors <- ncol(trajectory)
  traj_dimnames <- dimnames(trajectory)

  
  # Call the Rcpp forwarder
  result <- tox_normalize_single_trajectory_rcpp(trajectory)

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    traj_norm = matrix(result$traj_norm,
      nrow = n_timepoints,
      ncol = n_factors,
      dimnames = traj_dimnames),
    status = as.integer(result$status)
  ))
}

#> tox_trajectory_normalization:normalize_all_trajectories_c: Normalize all trajectories across multiple entities.
#' Normalize all trajectories across multiple entities
#'
#' @param trajectories 3D numeric array (n_factors x n_samples x n_timepoints)
#' @return A list with:
#' \describe{
#'   \item{traj_norm}{3D numeric array (n_factors x n_samples x n_timepoints) of normalized trajectories}
#'   \item{status}{Integer status code returned by the backend routine}
#' }
#'
tox_normalize_all_trajectories <- function(trajectories) {
  # Input Validation
  validate_numeric_array(trajectories)

  dims <- dim(trajectories)
  validate_length_equals_n(dims, 3)

  n_factors <- dims[1]
  n_samples <- dims[2]
  n_timepoints <- dims[3]
  traj_dimnames <- dimnames(trajectories)


  # Call the Rcpp forwarder
  result <- tox_normalize_all_trajectories_rcpp(trajectories,
                                                n_factors,
                                                n_samples,
                                                n_timepoints
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    traj_norm = array(as.numeric(result$traj_norm), dim = dims, dimnames = traj_dimnames),
    status = as.integer(result$status)
  ))
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

##################################################
## Helper functions for KD-tree
##################################################


#> f42_helper: Get a point from the KD index
#' Get a point from the KD index (R-level helper)
#' @param X A numeric matrix of shape (dimensions x points) representing the original data
#' @param kd_ix An integer vector of indices representing the KD-tree order of the points
#' @param position An integer scalar specifying the position in the KD index to retrieve
#' @return A numeric vector representing the point at the specified position in the KD index
get_kd_point <- function(X, kd_ix, position) {
  # Input validation using standardized validation functions
  validate_numeric_matrix(X, "X")
  validate_array_or_vector(kd_ix, "kd_ix")
  validate_numeric_scalar(position, "position")
  if (position != as.integer(position)) {
    stop("position must be an integer-valued scalar")
  }
  position <- as.integer(position)
  validate_index_bounds(position, low = 1, high = length(kd_ix), name = "position")
  
  # Extract the point at the specified position
  # Note: X is organized as (dimensions, points), so we need to extract column kd_ix[position]
  X[, kd_ix[position], drop = FALSE]
}
  
# ===================================================================
# SHIFT VECTOR FIELD FUNCTIONS
# ===================================================================
#> tox_shift_vectors:compute_shift_vector_field_c: Calculate shift vector field for gene expression
#' Calculate Shift Vector Field
#' Computes the shift vector field for each gene expression vector based on its family centroid.
#' The shift vector is defined as the difference between the gene expression vector and its corresponding family centroid,
#' starting at the expression vector and pointing to its family centroid.
#' This function automatically checks for errors and throws informative exceptions.
#'
#' @param expression_vectors Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param family_centroids Matrix where each column is a family centroid vector (n_axes x n_families)
#' @param gene_to_centroid Integer vector mapping each gene to its corresponding family centroid ID in `family_centroids` (length n_vectors)
#'
#' @return Numeric vector containing flattened shift-vector field output
#'
tox_compute_shift_vector_field <- function(expression_vectors, family_centroids, gene_to_centroid) {
  # Input Validation
  validate_numeric_matrix(expression_vectors)
  validate_numeric_matrix(family_centroids)
  

  # Ensure matching axes (rows)
  validate_matching_rows(expression_vectors, family_centroids)

  # Ensure mapping length matches number of vectors
  n_vectors <- ncol(expression_vectors)
  validate_logical_or_index_vector(gene_to_centroid, expected_length = n_vectors, name = "gene_to_centroid")
  gene_to_centroid <- as.integer(gene_to_centroid)

  # Call the Rcpp forwarder
  result <- tox_compute_shift_vector_field_rcpp(expression_vectors,
                                                family_centroids,
                                                gene_to_centroid
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return shift vectors
  return(result$shift_vectors)
}



# ===================================================================
# GENE CENTROIDS FUNCTIONS
# ===================================================================
#> tox_gene_centroids:group_centroid_c: Calculate centroids for gene families
#' Calculate Gene Centroids
#' Computes the centroids for each gene family based on the expression vectors of its member genes.
#' This function automatically checks for errors and throws informative exceptions.
#'
#' @param expression_vectors Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param gene_to_family Integer vector mapping each gene to its corresponding family ID (length n_vectors)
#' @param n_families Integer total number of gene families
#' @param ortholog_set Logical vector indicating if a gene is part of a specific subset (e.g., orthologs)
#' @param mode Character string indicating the mode of operation (`all` or `ortho`)
#'
#' @return Numeric matrix of shape (n_axes x n_families) with computed family centroids
#'
tox_group_centroid <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  # Input Validation
  validate_group_centroid_inputs(expression_vectors, gene_to_family, n_families, ortholog_set, mode)

  # Call the Rcpp forwarder
  result <- tox_group_centroid_rcpp(expression_vectors,
                                    gene_to_family,
                                    n_families,
                                    ortholog_set,
                                    mode
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return only the centroid matrix
  return(result$centroid_matrix)
}

#> tox_gene_centroids:mean_vector_c: Compute element-wise mean of gene expression vectors
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

  # Input Validation
  validate_mean_vector_inputs(expression_vectors, gene_indices)

  # Call the Rcpp forwarder
  result <- tox_mean_vector_rcpp(expression_vectors, gene_indices)

  # Check for errors
  check_err_code(result$ierr)

  # Return centroid vector
  return(result$centroid_col)
}


# ===================================================================
# DATA INTEGRATION FUNCTIONS
# ===================================================================

#> tox_data_integration:determine_shared_residual_range_expert_c: Compute shared residual range R from a precomputed residual pool
#' Compute the shared residual range R from a precomputed residual pool
#'
#' This function wraps the Fortran subroutine `determine_shared_residual_range_expert_c`
#' which computes the residual range from a flattened residual pool and its sorting permutation.
#'
#' @param residual_pool Numeric vector of absolute residuals (NaNs removed)
#' @param residual_pool_perm Integer vector giving the permutation that sorts `residual_pool`
#' @param residual_range_quantile Numeric scalar (default 95.0)
#'
#' @return Numeric scalar: the shared residual range R
#'
tox_determine_shared_residual_range_expert <- function(residual_pool, residual_pool_perm, residual_range_quantile = 95.0) {
  # Input Validation
  validate_numeric_vector(residual_pool)
  validate_integer_vector(residual_pool_perm)
  validate_nonempty_vector(residual_pool)
  validate_equal_length(residual_pool, residual_pool_perm)
  validate_numeric_scalar(residual_range_quantile)

  # Call the Rcpp forwarder
  result <- tox_determine_shared_residual_range_expert_rcpp(residual_pool,
                                                            residual_pool_perm,
                                                            residual_range_quantile
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return shared residual range
  return(result$shared_R)
}

#> tox_data_integration:determine_shared_residual_range_c: Compute shared residual range R from two residual matrices
#' Compute the shared residual range R from two residual matrices
#'
#' This function wraps the Fortran subroutine
#' `determine_shared_residual_range_c`, which computes the shared
#' residual range from two residual matrices (S1 and S2).
#'
#' @param neighborhood_residuals_S1 Numeric array (n_reps x n_neighbors x n_points)
#' @param neighborhood_residuals_S2 Numeric array (n_reps x n_neighbors x n_points)
#' @param residual_range_quantile Numeric scalar (default 95.0)
#'
#' @return Numeric scalar: the shared residual range R
#'
tox_determine_shared_residual_range <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, residual_range_quantile = 95.0) {
  # Input Validation
  validate_numeric_array(neighborhood_residuals_S1)
  validate_numeric_array(neighborhood_residuals_S2)
  validate_numeric_scalar(residual_range_quantile)

  # Call the Rcpp forwarder
  result <- tox_determine_shared_residual_range_rcpp(neighborhood_residuals_S1,
                                                     neighborhood_residuals_S2,
                                                     residual_range_quantile
  )
  
  # Check for errors
  check_err_code(result$ierr)

  # Return shared residual range
  return(result$shared_R)
}

#> tox_data_integration:build_residual_histograms_c: Build residual histograms and PMFs
#' Build residual histograms and probability mass functions (PMFs)
#'
#' This function wraps the Fortran subroutine
#' `build_residual_histograms_c`, which computes histogram bin counts
#' and normalized PMFs for each neighbor.
#'
#' @param neighborhood_residuals Numeric array (n_reps x n_neighbors x n_points)
#' @param shared_residual_range Numeric scalar R
#' @param n_bins Integer number of histogram bins
#'
#' @return A list with:
#' \describe{
#'   \item{counts}{Integer matrix (n_neighbors x n_bins)}
#'   \item{pmf}{Numeric matrix (n_neighbors x n_bins)}
#'   \item{included_n_residuals}{Integer vector (length n_points)}
#' }
#'
tox_build_residual_histograms <- function(neighborhood_residuals, shared_residual_range, n_bins) {

  # Input Validation
  validate_numeric_array(neighborhood_residuals)
  validate_length_equals_n(dim(neighborhood_residuals), 3)
  validate_numeric_scalar(shared_residual_range)
  validate_positive_integer_scalar(n_bins)

  # Call the Rcpp forwarder
  result <- tox_build_residual_histograms_rcpp(neighborhood_residuals,
                                               shared_residual_range,
                                               n_bins
  )
 
  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    counts = result$counts,
    pmf = result$pmf,
    included_n_residuals = result$included_n_residuals
  ))
}

#> tox_data_integration:build_residual_histograms_filtered_c: Build residual histograms and PMFs
#' Build residual histograms and probability mass functions (PMFs)
#'
#' This function wraps the Fortran subroutine
#' `build_residual_histograms_filtered_c`, which computes histogram bin counts
#' and normalized PMFs for each neighbor.
#'
#' @param neighborhood_residuals Numeric array (n_reps x n_neighbors x n_points)
#' @param shared_residual_range Numeric scalar R
#' @param n_bins Integer number of histogram bins
#' @param neighbor_mask Logical or integer matrix (n_neighbors x n_points), non-zero is TRUE.
#'
#' @return A list with:
#'   \describe{
#'     \item{counts}{Integer matrix (n_neighbors x n_bins)}
#'     \item{pmf}{Numeric matrix (n_neighbors x n_bins)}
#'     \item{included_n_residuals}{Integer vector (length n_points)}
#'   }
#'
tox_build_residual_histograms_filtered <- function(neighborhood_residuals, shared_residual_range, n_bins, neighbor_mask) {

  # Input Validation
  validate_numeric_array(neighborhood_residuals)
  dims <- dim(neighborhood_residuals)
  validate_length_equals_n(dims, 3)
  validate_numeric_scalar(shared_residual_range)
  validate_positive_integer_scalar(n_bins)

  # Convert mask to integer representation expected by backend
  mask_int <- neighbor_mask * 1L

  # Call the Rcpp forwarder
  result <- tox_build_residual_histograms_filtered_rcpp(neighborhood_residuals,
                                                        shared_residual_range,
                                                        n_bins,
                                                        mask_int
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    counts = result$counts,
    pmf = result$pmf,
    included_n_residuals = result$included_n_residuals
  ))
}

#> tox_data_integration:compute_divergence_per_reference_point_c: Compute per-neighbor Jensen–Shannon divergences
#' Compute per-neighbor Jensen–Shannon divergences
#'
#' This function wraps the Fortran subroutine
#' `compute_divergence_per_reference_point_c`, which computes the
#' Jensen–Shannon divergence for each neighbor based on two PMF matrices.
#'
#' @param pmf_S1 Numeric matrix (n_neighbors × n_bins)
#' @param pmf_S2 Numeric matrix (n_neighbors × n_bins)
#'
#' @return Numeric vector of length n_neighbors containing the JSD values
#'
tox_compute_divergence_per_reference_point <- function(pmf_S1, pmf_S2) {
  # Input Validation
  validate_numeric_matrix(pmf_S1)
  validate_numeric_matrix(pmf_S2)
  validate_matching_rows(pmf_S1, pmf_S2)

  # Call the Rcpp forwarder
  result <- tox_compute_divergence_per_reference_point_rcpp(pmf_S1, pmf_S2)

  # Check for errors
  check_err_code(result$ierr)

  # Return JSD values
  return(result$js_divergences)
}

#> tox_data_integration:compute_weighted_global_divergence_c: Compute weighted global Jensen–Shannon divergence
#' Compute weighted global Jensen–Shannon divergence
#'
#' This function wraps the Fortran subroutine
#' `compute_weighted_global_divergence_c`, which computes the weighted
#' global Jensen–Shannon divergence from per-neighbor divergences and
#' sample counts.
#'
#' @param js_divergences Numeric vector (length n_neighbors)
#' @param included_n_residuals_S1 Integer vector (length n_neighbors)
#' @param included_n_residuals_S2 Integer vector (length n_neighbors)
#'
#' @return A list with:
#'   \describe{
#'     \item{global_js_divergence}{Numeric scalar}
#'     \item{weights}{Numeric vector (length n_points)}
#'   }
#'
tox_compute_weighted_global_divergence <- function(js_divergences, included_n_residuals_S1, included_n_residuals_S2) {
  # Input Validation
  validate_numeric_vector(js_divergences)
  validate_integer_vector(included_n_residuals_S1)
  validate_integer_vector(included_n_residuals_S2)
  validate_equal_length(js_divergences, included_n_residuals_S1)
  validate_equal_length(js_divergences, included_n_residuals_S2)

  # Call the Rcpp forwarder
  result <- tox_compute_weighted_global_divergence_rcpp(js_divergences,
                                                        included_n_residuals_S1,
                                                        included_n_residuals_S2
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    global_js_divergence = result$global_js_divergence,
    weights = result$weights
  ))
}

#> tox_data_integration:gjct_permutation_test_c: Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
#' Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
#'
#' This function wraps the Fortran subroutine
#' `gjct_permutation_test_c`, which computes the global JSD values for permutations of the residuals of both studies.
#' To create a permutation, the residuals of a reference point will be concatenated, shuffled and reassigned in the shuffled order.
#' @param neighborhood_residuals_S1 Numeric array (n_reps x n_neighbors x n_points)
#' @param neighborhood_residuals_S2 Numeric array (n_reps x n_neighbors x n_points)
#' @param global_jsd_observed Numeric scalar of the observed weighted global JSD value
#' @param n_bins Integer number of bins used to compute `global_jsd_observed`
#' @param shared_residual_range Numeric scalar of shared residual range used to compute `global_jsd_observed`
#' @param n_permutations Integer number of permutations to perform
#' @param random_seed Optional integer random seed
#'
#' @return A list with:
#'   \describe{
#'     \item{jsd_null}{Numeric vector (length n_permutations)}
#'     \item{p_value}{Numeric scalar}
#'   }
#'
tox_gjct_permutation_test <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, global_jsd_observed, n_bins, shared_residual_range, n_permutations,
  random_seed) {
  # Input Validation
  validate_numeric_array(neighborhood_residuals_S1)
  validate_numeric_array(neighborhood_residuals_S2)
  validate_numeric_scalar(global_jsd_observed)
  validate_positive_integer_scalar(n_bins)
  validate_numeric_scalar(shared_residual_range)
  validate_positive_integer_scalar(n_permutations)
  

  # Call the Rcpp forwarder
  result <- tox_gjct_permutation_test_rcpp(neighborhood_residuals_S1,
                                           neighborhood_residuals_S2,
                                           global_jsd_observed,
                                           n_bins,
                                           shared_residual_range,
                                           n_permutations,
                                           random_seed
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    jsd_null = result$jsd_null,
    p_value = result$p_value,
    ierr = result$ierr
  ))
}

#> tox_data_integration:gjct_permutation_test_filtered_c: Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
#' Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
#'
#' This function wraps the Fortran subroutine
#' `gjct_permutation_test_filtered_c`, which computes the global JSD values for permutations of the residuals of both studies.
#' To create a permutation, the residuals of a reference point will be concatenated, shuffled and reassigned in the shuffled order.
#'
#' @param neighborhood_residuals_S1 Numeric array (n_reps x n_neighbors x n_points)
#' @param neighborhood_residuals_S2 Numeric array (n_reps x n_neighbors x n_points)
#' @param global_jsd_observed Numeric scalar of the observed weighted global JSD value
#' @param n_bins Integer number of bins used to compute `global_jsd_observed`
#' @param shared_residual_range Numeric scalar of shared residual range used to compute `global_jsd_observed`
#' @param n_permutations Integer number of permutations to perform
#' @param neighbor_mask_S1 Logical or integer matrix (n_neighbors x n_points), non-zero is TRUE
#' @param neighbor_mask_S2 Logical or integer matrix (n_neighbors x n_points), non-zero is TRUE
#' @param random_seed Optional integer random seed
#'
#' @return A list with:
#'   \describe{
#'     \item{jsd_null}{Numeric vector (length n_permutations)}
#'     \item{p_value}{Numeric scalar}
#'   }
#'
tox_gjct_permutation_test_filtered <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, global_jsd_observed, n_bins, shared_residual_range,
  n_permutations,
  neighbor_mask_S1,
  neighbor_mask_S2,
  random_seed
) {
  # Input Validation
  validate_numeric_array(neighborhood_residuals_S1)
  validate_numeric_array(neighborhood_residuals_S2)
  validate_numeric_scalar(global_jsd_observed)
  validate_positive_integer_scalar(n_bins)
  validate_numeric_scalar(shared_residual_range)
  validate_positive_integer_scalar(n_permutations)
  

  # Convert masks to integer representation expected by backend
  mask_int_S1 <- neighbor_mask_S1 * 1L
  mask_int_S2 <- neighbor_mask_S2 * 1L

  # Call the Rcpp forwarder
  result <- tox_gjct_permutation_test_filtered_rcpp(neighborhood_residuals_S1,
                                                    neighborhood_residuals_S2,
                                                    global_jsd_observed,
                                                    n_bins,
                                                    shared_residual_range,
                                                    n_permutations,
                                                    mask_int_S1,
                                                    mask_int_S2,
                                                    random_seed
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    jsd_null = result$jsd_null,
    p_value = result$p_value,
    ierr = result$ierr
  ))
}

#> tox_data_integration:compute_gene_means_c: Compute per-gene mean expression values
#' Compute per-gene mean expression values
#'
#' This function computes the mean expression for each gene across replicates,
#' handling NaN values appropriately.
#'
#' @param expr Numeric matrix of expression data (replicates × genes)
#' @return Numeric vector of per-gene mean expression values

tox_compute_gene_means <- function(expr) {
  # Input Validation
  validate_numeric_matrix(expr)
  
  # Call the Rcpp forwarder
  result <- tox_compute_gene_means_rcpp(expr)

  # Check for errors
  check_err_code(result$ierr)
  
  # Return means
  return(result$means)
}

#> tox_data_integration:compute_residuals_c: Compute signed residuals
#' Compute signed residuals
#'
#' This function computes signed residuals for each gene and replicate.
#'
#' @param expr Numeric matrix of expression data (replicates × genes)
#' @param means Numeric vector of per-gene mean expression values
#' @return Numeric matrix of signed residuals (replicates x genes)

tox_compute_residuals <- function(expr, means) {
  # Input Validation
  validate_numeric_matrix(expr)
  validate_numeric_vector(means)
  validate_length_equals_n(means, ncol(expr))
  
  # Call the Rcpp forwarder
  result <- tox_compute_residuals_rcpp(expr, means)

  # Check for errors
  check_err_code(result$ierr)
  
  # Return residual matrix
  return(result$resid)
}

#> tox_data_integration:pool_means_c: Pool mean expression values across studies
#' Pool mean expression values across studies
#'
#' This function pools per-gene mean expression values from two studies
#' and computes reference points for neighborhood construction.
#'
#' @param mean_S1 Numeric vector of per-gene mean expression values for study 1
#' @param mean_S2 Numeric vector of per-gene mean expression values for study 2
#' @param n_points Number of reference points to define
#' @return A list with:
#' \describe{
#'   \item{n_pool}{Integer count of valid (non-NA) pooled mean-expression values}
#'   \item{x_star}{Numeric vector of mean-expression reference points}
#' }

tox_pool_means <- function(mean_S1, mean_S2, n_points) {
  # Input Validation
  validate_numeric_vector(mean_S1)
  validate_numeric_vector(mean_S2)
  validate_positive_integer_scalar(n_points)
  
  # Call the Rcpp forwarder
  result <- tox_pool_means_rcpp(mean_S1,
                                mean_S2,
                                n_points
  )

  # Check for errors
  check_err_code(result$ierr)
  
  # Return structured result
  return(list(
    n_pool = result$n_pool,
    x_star = result$x_star
  ))
}

#> tox_data_integration:pool_means_expert_c: Pool mean expression values across studies
#' Pool mean expression values across studies
#'
#' This function pools per-gene mean expression values from two studies
#' and computes reference points for neighborhood construction.
#'
#' @param pooled_means Numeric vector of merged per-gene mean expression values from both studies
#' @param pooled_perm Integer vector of indices that would sort pooled_means (NaN values should be last values in sorting)
#' @param n_points Number of reference points to define
#' @return A list with:
#' \describe{
#'   \item{n_pool}{Integer count of valid (non-NA) pooled mean-expression values}
#'   \item{x_star}{Numeric vector of mean-expression reference points}
#' }

tox_pool_means_expert <- function(pooled_means, pooled_perm, n_points) {
  # Input Validation
  validate_numeric_vector(pooled_means)
  validate_logical_or_index_vector(pooled_perm)
  validate_same_length(pooled_means, pooled_perm)
  validate_positive_integer_scalar(n_points)
  
  # Call the Rcpp forwarder
  result <- tox_pool_means_expert_rcpp(pooled_means,
                                       pooled_perm,
                                       n_points
  )

  # Check for errors
  check_err_code(result$ierr)
  
  # Return structured result
  return(list(
    n_pool = result$n_pool,
    x_star = result$x_star
  ))
}

#> tox_data_integration:construct_neighborhoods_c: Construct neighborhood-based residual sets (kNN)
#' Construct neighborhood-based residual sets (kNN)
#'
#' This function constructs neighborhoods around reference points and
#' collects residuals from genes in each neighborhood.
#'
#' @param x_star Numeric vector of mean-expression reference points
#' @param n_pool Integer number of pooled means used to build neighborhoods
#' @param mean_S Numeric vector of per-gene mean expression values for the study
#' @param resid_S Numeric matrix of signed residuals for the study (replicates × genes)
#' @param desired_n_neighbors Desired size of the neighborhood, might be lower in the end (default=0, means automatic detection)
#' @return A list with:
#' \describe{
#'   \item{neighborhood_residuals}{Numeric array of residual vectors for each neighborhood}
#'   \item{neighborhood_indices}{Integer matrix of selected neighborhood gene indices (1-based)}
#' }
#' 
tox_construct_neighborhoods <- function(x_star, n_pool, mean_S, resid_S,
                                        desired_n_neighbors = 0) {
  # Input Validation
  validate_numeric_vector(x_star)
  validate_numeric_vector(mean_S)
  validate_numeric_matrix(resid_S)
  validate_positive_integer_scalar(n_pool)

  
  n_points <- length(x_star)
  
  # Call the Rcpp forwarder
  result <- tox_construct_neighborhoods_rcpp(x_star,
                                             n_pool,
                                             mean_S,
                                             resid_S,
                                             desired_n_neighbors
  )
  
  # Check for errors
  check_err_code(result$ierr)
  
  # Return structured result
  return(list(
    neighborhood_residuals = result$neighborhood_residuals,
    neighborhood_indices = result$neighborhood_indices
  ))
}

#> tox_data_integration:fjct_compute_jsd_c: Compute family-level JSD
#' Compute family-level JSD
#'
#' @param family_idx Integer scalar.
#' @param gene_to_family_S1 Integer vector (n_genes_S1).
#' @param gene_to_family_S2 Integer vector (n_genes_S2).
#' @param neighborhood_residuals_S1 Numeric array (n_reps_S1 × n_neighbors × n_points).
#' @param neighborhood_residuals_S2 Numeric array (n_reps_S2 × n_neighbors × n_points).
#' @param neighborhood_genes_S1 Integer matrix (n_neighbors × n_points).
#' @param neighborhood_genes_S2 Integer matrix (n_neighbors × n_points).
#' @param n_bins Integer scalar.
#' @param shared_residual_range Numeric scalar.
#'
#' @return A list with:
#' \describe{
#'   \item{js_divergences}{Numeric vector of Jensen-Shannon divergences per reference point}
#'   \item{included_n_reps_S1}{Integer vector of included replicate counts for study 1}
#'   \item{included_n_reps_S2}{Integer vector of included replicate counts for study 2}
#'   \item{total_included_n_reps}{Integer vector of total included replicate counts}
#'   \item{global_js_divergence}{Numeric scalar weighted global divergence}
#'   \item{weights}{Numeric vector of per-reference weights}
#' }
tox_fjct_compute_jsd <- function(
  family_idx,
  gene_to_family_S1,
  gene_to_family_S2,
  neighborhood_residuals_S1,
  neighborhood_residuals_S2,
  neighborhood_genes_S1,
  neighborhood_genes_S2,
  n_bins,
  shared_residual_range
) {
  # Input Validation
  validate_numeric_array(neighborhood_residuals_S1)
  validate_numeric_array(neighborhood_residuals_S2)
  validate_positive_integer_scalar(family_idx)
  validate_integer_vector(gene_to_family_S1)
  validate_integer_vector(gene_to_family_S2)
  validate_numeric_matrix(neighborhood_genes_S1)
  validate_numeric_matrix(neighborhood_genes_S2)

  # Call the Rcpp forwarder
  result <- tox_fjct_compute_jsd_alloc_rcpp(family_idx,
                                            gene_to_family_S1,
                                            gene_to_family_S2,
                                            neighborhood_residuals_S1,
                                            neighborhood_residuals_S2,
                                            neighborhood_genes_S1,
                                            neighborhood_genes_S2,
                                            n_bins,
                                            shared_residual_range
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    js_divergences = result$js_divergences,
    included_n_reps_S1 = result$included_n_reps_S1,
    included_n_reps_S2 = result$included_n_reps_S2,
    total_included_n_reps = result$total_included_n_reps,
    global_js_divergence = result$global_js_divergence,
    weights = result$weights
  ))
}

#> tox_data_integration:fjct_compute_jsd_expert_c: Compute family-level JSD (expert variant with masks)
#' Compute family-level JSD (expert variant with masks)
#'
#' @param neighborhood_residuals_S1 Numeric array (n_reps_S1 × n_neighbors × n_points).
#' @param neighborhood_residuals_S2 Numeric array (n_reps_S2 × n_neighbors × n_points).
#' @param neighbor_mask_S1 Logical or integer matrix (n_neighbors x n_points), non-zero is TRUE.
#' @param neighbor_mask_S2 Logical or integer matrix (n_neighbors x n_points), non-zero is TRUE.
#' @param n_bins Integer scalar.
#' @param shared_residual_range Numeric scalar.
#'
#' @return A list with:
#' \describe{
#'   \item{js_divergences}{Numeric vector of Jensen-Shannon divergences per reference point}
#'   \item{included_n_reps_S1}{Integer vector of included replicate counts for study 1}
#'   \item{included_n_reps_S2}{Integer vector of included replicate counts for study 2}
#'   \item{total_included_n_reps}{Integer vector of total included replicate counts}
#'   \item{global_js_divergence}{Numeric scalar weighted global divergence}
#'   \item{weights}{Numeric vector of per-reference weights}
#'   \item{pmf_S1}{Numeric matrix or array with PMFs for study 1}
#'   \item{pmf_S2}{Numeric matrix or array with PMFs for study 2}
#'   \item{tmp_counts}{Integer matrix or array with intermediate histogram counts}
#' }
tox_fjct_compute_jsd_expert <- function(
  neighborhood_residuals_S1,
  neighborhood_residuals_S2,
  neighbor_mask_S1,
  neighbor_mask_S2,
  n_bins,
  shared_residual_range
) {
  # Input Validation
  validate_numeric_array(neighborhood_residuals_S1)
  validate_numeric_array(neighborhood_residuals_S2)
  validate_positive_integer_scalar(n_bins)
  validate_numeric_scalar(shared_residual_range)
  


  # Convert masks to integer representation expected by backend
  mask_int_S1 <- neighbor_mask_S1 * 1L
  mask_int_S2 <- neighbor_mask_S2 * 1L

  # Call the Rcpp forwarder
  result <- tox_fjct_compute_jsd_expert_rcpp(neighborhood_residuals_S1,
                                             neighborhood_residuals_S2,
                                             mask_int_S1,
                                             mask_int_S2,
                                             n_bins,
                                             shared_residual_range
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    js_divergences = result$js_divergences,
    included_n_reps_S1 = result$included_n_reps_S1,
    included_n_reps_S2 = result$included_n_reps_S2,
    total_included_n_reps = result$total_included_n_reps,
    global_js_divergence = result$global_js_divergence,
    weights = result$weights,
    pmf_S1 = result$pmf_S1,
    pmf_S2 = result$pmf_S2,
    tmp_counts = result$tmp_counts
  ))
}

#> tox_data_integration:fjct_compute_contribution_scores_c: Compute per-family contribution scores
#' Compute per-family contribution scores
#'
#' @param global_js_divergences Numeric vector (k_families).
#' @param total_included_n_reps_per_f Integer vector (k_families).
#'
#' @return A list with:
#' \describe{
#'   \item{support_weights}{Numeric vector of support weights per family}
#'   \item{contribution_scores}{Numeric vector of normalized contribution scores per family}
#' }
tox_fjct_compute_contribution_scores <- function(
  global_js_divergences,
  total_included_n_reps_per_f
) {
  # Input Validation
  validate_numeric_vector(global_js_divergences)
  validate_integer_vector(total_included_n_reps_per_f)
  validate_same_length(global_js_divergences, total_included_n_reps_per_f)

  # Call the Rcpp forwarder
  result <- tox_fjct_compute_contribution_scores_rcpp(global_js_divergences, total_included_n_reps_per_f)

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    support_weights = result$support_weights,
    contribution_scores = result$contribution_scores
  ))
}

# ===================================================================
# PARALOG ANALYSIS FUNCTIONS
# ===================================================================

#> tox_paralog_analysis:detect_neofunctionalization_c: Identify neofunctionalization for genes
#' Identify neofunctionalization for genes
#'
#' Checks whether expression differences between genes and their mapped ancestors
#' exceed per-axis thresholds.
#'
#' @param ancestors Numeric matrix (n_axes x n_families) of ancestor vectors.
#' @param genes Numeric matrix (n_axes x n_genes) of gene vectors.
#' @param gene_to_fam Integer vector of length n_genes mapping genes to family indices (1-based).
#' @param thresholds Numeric vector of length n_axes with per-axis thresholds.
#' @return Logical matrix of shape (n_genes x n_axes), where TRUE indicates neofunctionalization.
tox_detect_neofunctionalization <- function(ancestors, genes, gene_to_fam, thresholds) {
  # Input Validation
  validate_numeric_matrix(ancestors)
  validate_numeric_matrix(genes)

  n_axes <- nrow(ancestors)
  n_families <- ncol(ancestors)
  n_genes <- ncol(genes)

 

  validate_length_equals_n(gene_to_fam, n_genes)
  validate_index_bounds(gene_to_fam, low = 1, high = n_families)
  validate_length_equals_n(thresholds, n_axes)

  # Call the Rcpp forwarder
  result <- tox_detect_neofunctionalization_rcpp(ancestors,
                                                 genes,
                                                 gene_to_fam,
                                                 thresholds
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return neofunctionalization mask
  return(result$neofunc)
}

#> tox_paralog_analysis:mask_check_state_c: Check the state of a specific gene in a bit mask
#' Check the state of a specific gene in a bit mask.
#'
#' @param bit_mask Integer vector representing a bit mask (chunks of 32 bits).
#' @param i_gene Integer index of the gene to check (1-based).
#' @return Logical value indicating inactive (FALSE) or active (TRUE).
tox_mask_check_state <- function(bit_mask, i_gene) {
  # Input Validation
  
  validate_integer_vector(bit_mask)
  validate_positive_integer_scalar(i_gene)
  validate_index_bounds(i_gene, low = 1, high = length(bit_mask) * 32)

  # Call the Rcpp forwarder
  result <- tox_mask_check_state_rcpp(bit_mask, i_gene)

  # Check for errors
  check_err_code(result$ierr)

  # Return bit state
  return(as.logical(result$state))
}

#> tox_paralog_analysis:mask_chunk_count_c: Compute the number of 32-bit chunks needed to encode a given number of genes
#' Compute number of 32-bit chunks needed for a gene mask.
#'
#' @param n_genes Number of genes to encode.
#' @return Integer count of 32-bit chunks.
tox_mask_chunk_count <- function(n_genes) {
  # Input Validation

  validate_positive_integer_scalar(n_genes)

  # Call the Rcpp forwarder
  result <- tox_mask_chunk_count_rcpp(n_genes)

  # Check for errors
  check_err_code(result$ierr)

  # Return chunk count
  return(as.integer(result$count))
}

#> tox_paralog_analysis:calc_work_arr_paralog_subsets_size_c: Calculate the required work array size for paralog subset analysis
#' Calculate required work array size for paralog subset analysis.
#'
#' @param max_subset_size Desired maximum subset size.
#' @param n_genes Total number of genes.
#' @param filtered_paralogs_mask Integer bit mask (chunks of 32 bits).
#' @return A list with:
#' \describe{
#'   \item{actual_max_subset_size}{Integer adjusted maximum subset size used by backend}
#'   \item{work_array_size}{Integer required work-array size}
#' }
tox_calc_work_arr_paralog_subsets_size <- function(max_subset_size, n_genes, filtered_paralogs_mask) {
  # Input Validation


  validate_positive_integer_scalar(max_subset_size)
  validate_positive_integer_scalar(n_genes)
  validate_integer_vector(filtered_paralogs_mask)

  # Call the Rcpp forwarder
  result <- tox_calc_work_arr_paralog_subsets_size_rcpp(max_subset_size,
                                                        n_genes,
                                                        filtered_paralogs_mask
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    actual_max_subset_size = as.integer(result$actual_max_subset_size),
    work_array_size = as.integer(result$work_array_size)
  ))
}

#> tox_paralog_analysis:filter_paralogs_by_pattern_dosage_effect_c: Filter paralogs by dosage effect using angle threshold, within the grouped slice
#' Filter paralogs by dosage-effect pattern.
#'
#' @param gene_angles Numeric vector of gene angles.
#' @param threshold Numeric filtering threshold.
#' @param gene_to_fam Integer vector mapping genes to family indices.
#' @param n_families Number of families.
#' @return Integer matrix mask (n_mask_chunks x n_families).
tox_filter_paralogs_by_pattern_dosage_effect <- function(gene_angles, threshold, gene_to_fam, n_families) {
  # Input Validation

  validate_numeric_vector(gene_angles)
  validate_numeric_scalar(threshold)
  validate_length_equals_n(gene_to_fam, length(gene_angles))
  validate_positive_integer_scalar(n_families)
  validate_index_bounds(gene_to_fam, low = 1, high = n_families)

  # Call the Rcpp forwarder
  result <- tox_filter_paralogs_by_pattern_dosage_effect_rcpp(gene_angles,
                                                              threshold,
                                                              gene_to_fam,
                                                              n_families
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return family-wise masks
  return(result$masks)
}

#> tox_paralog_analysis:filter_paralogs_by_pattern_subfunctionalization_c: Filter paralogs by subfunctionalization pattern using angle threshold, within the grouped slice
#' Filter paralogs by subfunctionalization pattern.
#'
#' @param gene_angles Numeric vector of gene angles.
#' @param threshold Numeric filtering threshold.
#' @param gene_to_fam Integer vector mapping genes to family indices.
#' @param n_families Number of families.
#' @return Integer matrix mask (n_mask_chunks x n_families).
tox_filter_paralogs_by_pattern_subfunctionalization <- function(gene_angles, threshold, gene_to_fam, n_families) {
  # Input Validation


  validate_numeric_vector(gene_angles)
  validate_numeric_scalar(threshold)
  validate_length_equals_n(gene_to_fam, length(gene_angles))
  validate_positive_integer_scalar(n_families)
  validate_index_bounds(gene_to_fam, low = 1, high = n_families)

  # Call the Rcpp forwarder
  result <- tox_filter_paralogs_by_pattern_subfunctionalization_rcpp(gene_angles,
                                                                     threshold,
                                                                     gene_to_fam,
                                                                     n_families
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return family-wise masks
  return(result$masks)
}

#> tox_paralog_analysis:detect_subfunctionalization_c: Detect subfunctionalization among paralogs based on residual distance and pruning
#' Detect subfunctionalization among paralogs.
#'
#' @param ancestor Numeric vector (n_dims).
#' @param genes Numeric matrix (n_dims x n_genes).
#' @param rdi_threshold Numeric residual distance threshold.
#' @param filtered_paralogs_mask Integer bit mask (chunks of 32 bits).
#' @param max_subset_size Desired maximum subset size.
#' @param paralog_norms Numeric vector of paralog norms (length n_genes).
#' @param sorted_paralog_norms_perm Integer permutation vector (length n_genes).
#' @return A list with:
#' \describe{
#'   \item{n_results}{Integer number of accepted subsets}
#'   \item{results}{Integer matrix containing subset bitmasks}
#' }
tox_detect_subfunctionalization <- function(ancestor, genes, rdi_threshold,
                                            filtered_paralogs_mask, max_subset_size,
                                            paralog_norms, sorted_paralog_norms_perm) {
  # Input Validation


  validate_numeric_vector(ancestor)
  validate_numeric_matrix(genes)
  validate_length_equals_n(ancestor, nrow(genes))
  validate_numeric_scalar(rdi_threshold)
  validate_integer_vector(filtered_paralogs_mask)
  validate_positive_integer_scalar(max_subset_size)

  n_genes <- ncol(genes)
  validate_length_equals_n(paralog_norms, n_genes)
  validate_length_equals_n(sorted_paralog_norms_perm, n_genes)
  validate_index_bounds(sorted_paralog_norms_perm, low = 1, high = n_genes)

  # Call the Rcpp forwarder
  result <- tox_detect_subfunctionalization_rcpp(ancestor,
                                                 genes,
                                                 rdi_threshold,
                                                 filtered_paralogs_mask,
                                                 max_subset_size,
                                                 paralog_norms,
                                                 sorted_paralog_norms_perm
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    n_results = as.integer(result$n_results),
    results = result$results
  ))
}

#> tox_paralog_analysis:detect_dosage_effect_c: Detect dosage effect among paralogs using gain and angle thresholds
#' Detect dosage effect among paralogs.
#'
#' @param ancestor Numeric vector (n_dims).
#' @param genes Numeric matrix (n_dims x n_genes).
#' @param filtered_paralogs_mask Integer bit mask (chunks of 32 bits).
#' @param max_subset_size Desired maximum subset size.
#' @param gain_gamma Numeric gain threshold.
#' @param max_angle Numeric maximum angle threshold (radians).
#' @return A list with:
#' \describe{
#'   \item{n_results}{Integer number of accepted subsets}
#'   \item{results}{Integer matrix containing subset bitmasks}
#' }
tox_detect_dosage_effect <- function(ancestor, genes,
                                     filtered_paralogs_mask, max_subset_size,
                                     gain_gamma = 0.1, max_angle = pi) {
  # Input Validation
 
  validate_numeric_vector(ancestor)
  validate_numeric_matrix(genes)
  validate_length_equals_n(ancestor, nrow(genes))
  validate_integer_vector(filtered_paralogs_mask)
  validate_positive_integer_scalar(max_subset_size)
  validate_positive_numeric_scalar(gain_gamma)
  validate_positive_numeric_scalar(max_angle)

  # Call the Rcpp forwarder
  result <- tox_detect_dosage_effect_rcpp(ancestor,
                                          genes,
                                          filtered_paralogs_mask,
                                          max_subset_size,
                                          gain_gamma,
                                          max_angle
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    n_results = as.integer(result$n_results),
    results = result$results
  ))
}

# ============================================================
# ARRAY UTILITIES FUNCTIONS
# ============================================================

#> f42_array_utils:get_array_metadata_c: Get array metadata from serialized file
#' Get array metadata from serialized file
#' @param filename Path to the serialized array file
#' @param max_dims Maximum number of dimensions to read from the file (default: 5)
#' @param with_clen Logical indicating whether to also read character length metadata (default: FALSE)
#' @return A list containing:
#' \describe{
#'   \item{dims}{Integer vector of array dimensions}
#'   \item{ndim}{Integer scalar of number of dimensions}
#'   \item{clen}{Integer scalar of character length (only if with_clen = TRUE)}
#' }
tox_get_array_metadata <- function(filename, max_dims = 5L, with_clen = FALSE) {
  # Input Validation
  validate_filename(filename)
  validate_max_dims(max_dims)

  
  # Call the Rcpp forwarder
  result <- tox_get_array_metadata_rcpp(filename,
                                        max_dims,
                                        with_clen
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return metadata shape
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
#  DESERIALIZATION FUNCTIONS 
# ============================================================
#> f42_deserialize_int:deserialize_int_nd_C: Deserialize integer array from file
#' @param filename Path to the serialized integer array file
#' @param max_dims Maximum number of dimensions to read from the file (default: 5)
#' @return An integer array with dimensions as specified in the file
tox_deserialize_int_array <- function(filename, max_dims = 5L) {
  # Input Validation
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)

 
  # Call the Rcpp forwarder
  result <- tox_deserialize_int_array_rcpp(filename, max_dims)

  # Check for errors
  check_err_code(result$ierr)

  # Return deserialized array
  return(array(result$values, dim = result$dims[1:result$ndim]))
}

#> f42_deserialize_real:deserialize_real_nd_C: Deserialize real/double array from file
#' @param filename Path to the serialized real array file
#' @param max_dims Maximum number of dimensions to read from the file (default: 5)
#' @return A numeric array with dimensions as specified in the file
tox_deserialize_real_array <- function(filename, max_dims = 5L) {
  # Input Validation
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)



  # Call the Rcpp forwarder
  result <- tox_deserialize_real_array_rcpp(filename, max_dims)

  # Check for errors
  check_err_code(result$ierr)

  # Return deserialized array
  return(array(result$values, dim = result$dims[1:result$ndim]))
}

#> f42_deserialize_char:deserialize_char_nd_C: Deserialize character array from file
#' @param filename Path to the serialized character array file
#' @param max_dims Maximum number of dimensions to read from the file (default: 5)
#' @return A character array with dimensions as specified in the file
tox_deserialize_char_array <- function(filename, max_dims = 5L) {
  # Input Validation
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)

  # Call the Rcpp forwarder
  result <- tox_deserialize_char_array_rcpp(filename, max_dims)

  # Check for errors
  check_err_code(result$ierr)

  # Return deserialized array
  return(array(result$values, dim = result$dims[1:result$ndim]))
}
#> f42_deserialize_logical:deserialize_logical_nd_C: Deserialize logical array from file
#' @param filename Path to the serialized logical array file
#' @param max_dims Maximum number of dimensions to read from the file (default: 5)
#' @return A logical array with dimensions as specified in the file
tox_deserialize_logical_array <- function(filename, max_dims = 5L) {
  # Input Validation
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)
  

  # Call the Rcpp forwarder
  result <- tox_deserialize_logical_array_rcpp(filename, max_dims)

  # Check for errors
  check_err_code(result$ierr)

  # Return deserialized array
  return(array(result$values, dim = result$dims[1:result$ndim]))
}

#> f42_deserialize_complex:deserialize_complex_nd_C: Deserialize complex array from file
#' @param filename Path to the serialized complex array file
#' @param max_dims Maximum number of dimensions to read from the file (default: 5)
#' @return A complex array with dimensions as specified in the file
tox_deserialize_complex_array <- function(filename, max_dims = 5L) {
  # Input Validation
  validate_filename(filename)
  validate_max_dims(max_dims)
  validate_file_exists(filename)
  

  # Call the Rcpp forwarder
  result <- tox_deserialize_complex_array_rcpp(filename, max_dims)

  # Check for errors
  check_err_code(result$ierr)

  # Return deserialized array
  return(array(result$values, dim = result$dims[1:result$ndim]))
}

# ============================================================
#  SERIALIZATION FUNCTIONS 
# ============================================================

#> f42_serialize_int:serialize_int_nd_C: Serialize integer array to file
#' @param arr An integer array to serialize
#' @param filename Path to the output file where the array will be serialized
#' @return NULL (invisible) - the function is called for its side effect of writing to a file
tox_serialize_int_array <- function(arr, filename) {
  # Input Validation
  validate_array_or_vector(arr)
  validate_filename(filename)
  
  
  # Call the Rcpp forwarder
  ierr <- tox_serialize_int_array_rcpp(arr, filename)

  # Check for errors
  check_err_code(ierr)

  # Return invisibly on success
  return(invisible(NULL))
}

#> f42_serialize_real:serialize_real_nd_C: Serialize real/double array to file
#' @param arr A numeric array to serialize
#' @param filename Path to the output file where the array will be serialized
#' @return NULL (invisible) - the function is called for its side effect of writing to a file
tox_serialize_real_array <- function(arr, filename) {
  # Input Validation
  validate_array_or_vector(arr)
  validate_filename(filename)



  # Call the Rcpp forwarder
  ierr <- tox_serialize_real_array_rcpp(arr, filename)

  # Check for errors
  check_err_code(ierr)

  # Return invisibly on success
  return(invisible(NULL))
}

#> f42_serialize_char:serialize_char_nd_C: Serialize character array to file
#' Serialize character array to file
#'
#' @param arr A character array to serialize
#' @param filename Path to the output file where the array will be serialized
#' @return NULL (invisible) - the function is called for its side effect of writing to a file
tox_serialize_char_array <- function(arr, filename) {
  # Input Validation
  validate_array_or_vector(arr)
  validate_filename(filename)



  # Call the Rcpp forwarder
  ierr <- tox_serialize_char_array_rcpp(arr, filename)

  # Check for errors
  check_err_code(ierr)

  # Return invisibly on success
  return(invisible(NULL))
}

#> f42_serialize_logical:serialize_logical_nd_C: Serialize logical array to file
#' @param arr A logical array to serialize
#' @param filename Path to the output file where the array will be serialized
#' @return NULL (invisible) - the function is called for its side effect of writing to a file
tox_serialize_logical_array <- function(arr, filename) {
  # Input Validation
  validate_array_or_vector(arr)
  validate_filename(filename)
  


  # Call the Rcpp forwarder
  ierr <- tox_serialize_logical_array_rcpp(arr, filename)

  # Check for errors
  check_err_code(ierr)

  # Return invisibly on success
  return(invisible(NULL))
}

#> f42_serialize_complex:serialize_complex_nd_C: Serialize complex array to file
#' @param arr A complex array to serialize
#' @param filename Path to the output file where the array will be serialized
#' @return NULL (invisible) - the function is called for its side effect of writing to a file
tox_serialize_complex_array <- function(arr, filename) {
  # Input Validation
  validate_array_or_vector(arr)
  validate_filename(filename)


  # Call the Rcpp forwarder
  ierr <- tox_serialize_complex_array_rcpp(arr, filename)

  # Check for errors
  check_err_code(ierr)

  # Return invisibly on success
  return(invisible(NULL))
}



# ============================================================
# KD TREE FUNCTIONS
# ============================================================

#> f42_kd_tree:build_kd_index_c: Build KD-tree index for multidimensional data
#' Build KD-tree index for multidimensional data
#'
#' @param X A numeric matrix of shape (dimensions x points) representing the data to index
#' @param dim_order Integer vector specifying the order of dimensions used for splitting
#'   (default: `1:nrow(X)`).
#' @return An integer vector of indices representing the KD-tree order of the points
build_kd_index <- function(X, dim_order = NULL) {
  # Input Validation
  validate_numeric_matrix(X)

  d <- nrow(X)
  
  
  if (is.null(dim_order)) {
    dim_order <- seq_len(d)
  }
  
  # Validate and coerce split dimension order
  validate_logical_or_index_vector(dim_order, expected_length = d, name = "dim_order")
  
  validate_length_equals_n(dim_order, d)
  
  
  # Call the Rcpp forwarder
  result <- tox_build_kd_index_rcpp(X, dim_order)
  
  # Check for errors
  check_err_code(result$ierr)
  
  # Return KD index vector
  return(result$kd_ix)
}

#> f42_kd_tree:build_spherical_kd_c: Build spherical KD-tree index
#' Build spherical KD-tree index for multidimensional data
#'
#' @param V A numeric matrix of shape (dimensions x points) representing the data to index
#' @param dim_order Integer vector specifying the order of dimensions used for splitting
#'   (default: `1:nrow(V)`).
#' @return An integer vector of indices representing the spherical KD-tree order of the points
build_spherical_kd <- function(V, dim_order = NULL) {
  # Input Validation
  validate_numeric_matrix(V)

  d <- nrow(V)
  n <- ncol(V)


  if (is.null(dim_order)) {
    dim_order <- seq_len(d)
  }
  
  # Validate and coerce split dimension order
  validate_logical_or_index_vector(dim_order, expected_length = d, name = "dim_order")
  
  validate_length_equals_n(dim_order, d)
  
  # Call the Rcpp forwarder
  result <- tox_build_spherical_kd_rcpp(V, dim_order)
  
  # Check for errors
  check_err_code(result$ierr)
  
  # Return sphere index vector
  return(result$sphere_ix)
}


# ============================================================
#  BINARY SEARCH TREE FUNCTIONS
# ============================================================

#> f42_binary_search_tree:build_bst_index_c: Build binary search tree index for a numeric vector
#' Build binary search tree index for a numeric vector
#'
#' @param x A numeric vector for which to build the BST index.
#' @return Integer vector of indices representing the BST order of the input vector.

build_bst_index <- function(x) {
  # Input Validation
  validate_numeric_vector(x)



  # Call the Rcpp forwarder
  result <- tox_build_bst_index_rcpp(x)

  # Return index vector
  return(result)
}
 
#> f42_helper: Get sorted value by BST index
#' Get sorted value by BST index (R helper)
get_sorted_value <- function(x, ix, position) {
  # Input Validation
  validate_numeric_vector(x)
  validate_array_or_vector(ix)
  validate_numeric_scalar(position, "position")
  if (position != as.integer(position)) {
    stop("position must be an integer-valued scalar")
  }
  position <- as.integer(position)
  validate_index_bounds(position, low = 1, high = length(ix))

  # Return sorted value at selected position
  return(as.numeric(x)[ix[position]])
}


#> f42_binary_search_tree:bst_range_query_c: Query BST index for values in range
#' Query BST index for values in range
#'
#' @param x A numeric vector of values corresponding to the BST index.
#' @param ix An integer vector of indices representing the BST order of the input vector
#' @param lo A numeric scalar representing the lower bound of the query range
#' @param hi A numeric scalar representing the upper bound of the query range
#' @return A list containing:
#' \describe{
#'   \item{indices}{Integer vector of original indices of values in the range [lo, hi]}
#'   \item{count}{Integer scalar of the number of values in the range}
#' }
#' 
bst_range_query <- function(x, ix, lo, hi) {
  # Input Validation
  validate_numeric_vector(x)
  validate_integer_vector(ix)
  validate_equal_length(x, ix)
  validate_numeric_scalar(lo)
  validate_numeric_scalar(hi)

 

  # Call the Rcpp forwarder
  result <- tox_bst_range_query_rcpp(x,
                                     ix,
                                     lo,
                                     hi
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return standardized public structure
  return(list(
    indices = result$output_indices,
    count = result$num_matches
  ))
}



# ============================================================
#  RELATIVE AXIS PLANE TOOLS FUNCTIONS
# ============================================================

#> tox_relative_axis_plane_tools:omics_vector_RAP_projection_c: Project vectors onto Relative Axis Plane
#' Project selected vectors onto a Relative Axis Plane (RAP)
#' 
#' This function validates inputs and calls the C/Rcpp RAP projection wrapper.
#' @param vecs Numeric matrix (n_axes x n_vecs).
#' @param vecs_selection_mask Logical or integer vector (length n_vecs).
#' @param axes_selection_mask Logical or integer vector (length n_axes).
#' @return Numeric matrix of projections (n_selected_axes x n_selected_vecs).
#'
tox_omics_vector_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
  # Input Validation
  validate_numeric_matrix(vecs)
  n_axes <- nrow(vecs)
  n_vecs <- ncol(vecs)
  validate_logical_or_index_vector(vecs_selection_mask, expected_length = n_vecs, name = "vecs_selection_mask")
  validate_logical_or_index_vector(axes_selection_mask, expected_length = n_axes, name = "axes_selection_mask")

  # Convert selection masks to integer 0/1 representation expected by backend
  vecs_selection_mask <- as.integer(as.logical(vecs_selection_mask))
  axes_selection_mask <- as.integer(as.logical(axes_selection_mask))
  

  # Call the Rcpp forwarder
  result <- tox_omics_vector_RAP_projection_rcpp(vecs,
                                                 vecs_selection_mask,
                                                 axes_selection_mask
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return RAP projections
  return(result$projections)
}

#> tox_relative_axis_plane_tools:omics_field_RAP_projection_c: Project vector fields onto Relative Axis Plane
#' Project selected vector fields onto a Relative Axis Plane (RAP)
#' 
#' This function validates inputs and calls the C/Rcpp RAP field projection wrapper.
#' @param vecs Numeric matrix (2*n_axes x n_vecs).
#' @param vecs_selection_mask Logical or integer vector (length n_vecs).
#' @param axes_selection_mask Logical or integer vector (length n_axes).
#' @return Numeric matrix of projections (n_selected_axes x n_selected_vecs).
#' 
tox_omics_field_RAP_projection <- function(vecs, vecs_selection_mask, axes_selection_mask) {
  # Input Validation
  validate_numeric_matrix(vecs)
  
  n_rows <- nrow(vecs)
  n_axes <- as.integer(n_rows / 2L)
  n_vecs <- ncol(vecs)

  validate_logical_or_index_vector(vecs_selection_mask, expected_length = n_vecs, name = "vecs_selection_mask")
  validate_logical_or_index_vector(axes_selection_mask, expected_length = n_axes, name = "axes_selection_mask")

  # Convert selection masks to integer 0/1 representation expected by backend
  vecs_selection_mask <- as.integer(as.logical(vecs_selection_mask))
  axes_selection_mask <- as.integer(as.logical(axes_selection_mask))
  

  # Call the Rcpp forwarder
  result <- tox_omics_field_RAP_projection_rcpp(vecs,
                                                vecs_selection_mask,
                                                axes_selection_mask
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return RAP projections
  return(result$projections)
}

#> tox_relative_axis_plane_tools:relative_axes_changes_from_shift_vector_c: Compute relative axis changes from shift vector
#' Compute relative axis changes from a RAP-projected and normalized shift vector
#' 
#' @param vec Numeric vector (RAP-projected and normalized shift vector).
#' @return Numeric vector of fractional axis contributions (sums to 1).
#' 
tox_relative_axes_changes_from_shift_vector <- function(vec) {
  # Input Validation
  validate_numeric_vector(vec)

  # Preserve stable zero-vector behavior expected by tests/callers.
  if (length(vec) > 0L && all(vec == 0, na.rm = TRUE)) {
    return(rep(0, length(vec)))
  }
  
  # Call the Rcpp forwarder
  result <- tox_relative_axes_changes_from_shift_vector_rcpp(vec)

  # Check for errors
  check_err_code(result$ierr)

  # Return axis contributions
  return(result$contributions)
}


#> tox_relative_axis_plane_tools:clock_hand_angle_between_vectors_c: Compute signed angle between two vectors
#' Compute signed clock hand angle between two RAP-projected and normalized vectors
#' 
#' @param v1 Numeric vector (first normalized vector in RAP space).
#' @param v2 Numeric vector (second normalized vector in RAP space).
#' @param selected_axes_for_signed Integer vector of length 3 (axes for directionality,
#'   ignored if length(v1) <= 3). Default: `c(1, 2, 3)`.
#' @return Signed angle in radians in `[-pi, pi]`.
#' 
tox_clock_hand_angle_between_vectors <- function(v1, v2, selected_axes_for_signed = c(1L, 2L, 3L)) {
  # Input Validation
  validate_numeric_vector(v1)
  validate_numeric_vector(v2)
  validate_same_length(v1, v2)

  n_dims <- length(v1)

  if (missing(selected_axes_for_signed)) {
    if (n_dims <= 1L) {
      selected_axes_for_signed <- c(1L, 1L, 1L)
    } else if (n_dims == 2L) {
      selected_axes_for_signed <- c(1L, 2L, 1L)
    } else {
      selected_axes_for_signed <- c(1L, 2L, 3L)
    }
  }

  # Accept numeric index triplets and normalize to integer for backend.
  validate_logical_or_index_vector(selected_axes_for_signed, expected_length = 3, name = "selected_axes_for_signed")
  selected_axes_for_signed <- as.integer(selected_axes_for_signed)

  validate_length_equals_n(selected_axes_for_signed, 3)
  validate_index_bounds(selected_axes_for_signed, low = 1, high = n_dims)

  # Call the Rcpp forwarder
  result <- tox_clock_hand_angle_between_vectors_rcpp(v1,
                                                      v2,
                                                      selected_axes_for_signed
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return signed angle
  return(result$signed_angle)
}

#> tox_relative_axis_plane_tools:clock_hand_angles_for_shift_vectors_c: Compute signed rotation angles for vector pairs
#' Compute signed rotation angles for pairs of RAP-projected and normalized vectors
#'
#' @param origins Numeric matrix (n_dims x n_vecs), first set of vectors.
#' @param targets Numeric matrix (n_dims x n_vecs), second set of vectors.
#' @param vecs_selection_mask Integer or logical vector (length n_vecs).
#'   Default: `rep(1, ncol(origins))` (all vectors selected).
#' @param selected_axes_for_signed Integer vector of length 3 (axes for directionality).
#'   Default: `c(1, 2, 3)`.
#' @return Numeric vector of signed angles (radians).
#' 
tox_clock_hand_angles_for_shift_vectors <- function(origins, targets, vecs_selection_mask = rep(1L, ncol(origins)), selected_axes_for_signed = c(1L, 2L, 3L)) {
  # Input Validation
  validate_numeric_matrix(origins)
  validate_numeric_matrix(targets)
  validate_matching_rows(origins, targets)
  validate_matching_cols(origins, targets)

  n_vecs <- ncol(origins)
  validate_logical_or_index_vector(vecs_selection_mask, expected_length = n_vecs, name = "vecs_selection_mask")
  vecs_selection_mask <- as.integer(as.logical(vecs_selection_mask))

  n_dims <- nrow(origins)

  if (missing(selected_axes_for_signed)) {
    if (n_dims <= 1L) {
      selected_axes_for_signed <- c(1L, 1L, 1L)
    } else if (n_dims == 2L) {
      selected_axes_for_signed <- c(1L, 2L, 1L)
    } else {
      selected_axes_for_signed <- c(1L, 2L, 3L)
    }
  }

  # Accept numeric index triplets and normalize to integer for backend.
  validate_logical_or_index_vector(selected_axes_for_signed, expected_length = 3, name = "selected_axes_for_signed")
  selected_axes_for_signed <- as.integer(selected_axes_for_signed)

  validate_length_equals_n(selected_axes_for_signed, 3)
  validate_index_bounds(selected_axes_for_signed, low = 1, high = n_dims)

  # Call the Rcpp forwarder
  result <- tox_clock_hand_angles_for_shift_vectors_rcpp(origins,
                                                         targets,
                                                         vecs_selection_mask,
                                                         selected_axes_for_signed
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return signed angles
  return(result$signed_angles)
}

#> tox_relative_axis_plane_tools:relative_axes_expression_from_expression_vector_c: Compute relative axis contributions from expression vector
#' Compute relative axis contributions from a RAP-projected and normalized expression vector
#' 
#' @param vec Numeric vector (RAP-projected and normalized expression vector).
#' @return Numeric vector of fractional axis contributions (sums to 1).
#'
tox_relative_axes_expression_from_expression_vector <- function(vec) {
  # Input Validation
  validate_numeric_vector(vec)

  # Preserve stable zero-vector behavior expected by tests/callers.
  if (length(vec) > 0L && all(vec == 0, na.rm = TRUE)) {
    return(rep(0, length(vec)))
  }
  
  # Call the Rcpp forwarder
  result <- tox_relative_axes_expression_from_expression_vector_rcpp(vec)

  # Check for errors
  check_err_code(result$ierr)

  # Return axis contributions
  return(result$contributions)
}

# ============================================================
#  CLUSTERING FUNCTIONS
# ============================================================

#> tox_clustering:cluster_factor_trajectories_k_means_c: K-means clustering on factor trajectories
#'
#' Performs k-means clustering on factor trajectories (factor evolution over time).
#' All validation is performed in R.
#' @param n_clusters Number of clusters (integer)
#' @param trajectories Numeric vector (flattened, n_factors * n_samples * n_timepoints)
#' @param n_factors Number of factors (integer)
#' @param n_samples Number of samples (integer)
#' @param n_timepoints Number of timepoints (integer)
#' @param centroids Numeric matrix (n_factors x n_clusters)
#' @param max_iterations Maximum number of iterations (integer)
#' @return A list with:
#' \describe{
#'   \item{centroids}{Numeric matrix of updated centroids}
#'   \item{labels}{Integer vector of cluster labels}
#'   \item{label_counts}{Integer vector of cluster membership counts}
#' }
#' 
tox_cluster_factor_trajectories_k_means <- function(n_clusters, trajectories, n_factors, n_samples, n_timepoints, centroids, max_iterations = 300) {
  # Input Validation
  validate_positive_integer_scalar(n_clusters)
  validate_numeric_vector(trajectories)
  validate_positive_integer_scalar(n_factors)
  validate_positive_integer_scalar(n_samples)
  validate_positive_integer_scalar(n_timepoints)
  validate_numeric_matrix(centroids)
  validate_positive_integer_scalar(max_iterations)

  
  

  # Validate centroid shape
  validate_matrix_shape_factor_centroids(centroids, n_factors, n_clusters)

  # Call the Rcpp forwarder
  result <- tox_cluster_factor_trajectories_k_means_rcpp(trajectories = trajectories,
                                                         centroids = centroids,
                                                         n_clusters = n_clusters,
                                                         n_factors = n_factors,
                                                         n_samples = n_samples,
                                                         n_timepoints = n_timepoints,
                                                         max_iterations = max_iterations
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    centroids = result$centroids,
    labels = result$labels,
    label_counts = result$label_counts
  ))
}

#> tox_clustering:k_means_clustering_c: K-means clustering (general)
#'
#' Performs k-means clustering on general data points.
#' All validation is performed in R.
#' @param n_clusters Number of clusters (integer)
#' @param data_points Numeric matrix (n_dims x n_points)
#' @param n_points Number of points (integer)
#' @param n_dims Number of dimensions (integer)
#' @param centroids Numeric matrix (n_dims x n_clusters)
#' @param max_iterations Maximum number of iterations (integer)
#' @return A list with:
#' \describe{
#'   \item{centroids}{Numeric matrix of updated centroids}
#'   \item{labels}{Integer vector of cluster labels}
#'   \item{label_counts}{Integer vector of cluster membership counts}
#' }
#' 
tox_k_means_clustering <- function(n_clusters, data_points, n_points, n_dims, centroids, max_iterations = 300) {
  # Input Validation
  validate_positive_integer_scalar(n_clusters)
  validate_numeric_matrix(data_points)
  validate_positive_integer_scalar(n_points)
  validate_positive_integer_scalar(n_dims)
  validate_numeric_matrix(centroids)
  validate_positive_integer_scalar(max_iterations)

 

  # Validate matrix shapes
  validate_matrix_shape_data_points(data_points, n_dims, n_points)
  validate_matrix_shape_centroids(centroids, n_dims, n_clusters)

  # Call the Rcpp forwarder
  result <- tox_k_means_clustering_rcpp(n_clusters,
                                        data_points,
                                        n_points,
                                        n_dims,
                                        centroids,
                                        max_iterations
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    centroids = result$centroids,
    labels = result$labels,
    label_counts = result$label_counts
  ))
}

#> tox_clustering:linkage_clustering_c: Hierarchical linkage clustering
#' Perform hierarchical linkage clustering
#'
#' @param distances Numeric square distance matrix (n_points x n_points)
#' @param method Linkage method: "average", "weighted", or "ward"
#' @return A list with:
#' \describe{
#'   \item{merge_i}{Integer vector of first merge indices}
#'   \item{merge_j}{Integer vector of second merge indices}
#'   \item{heights}{Numeric vector of merge heights}
#'   \item{cluster_sizes}{Integer vector of cluster sizes after each merge}
#' }
tox_linkage_clustering <- function(distances, method = "average") {
  # Input Validation
  validate_numeric_matrix(distances)
 
 

  # Call the Rcpp forwarder
  result <- tox_linkage_clustering_rcpp(distances, method)

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    merge_i = result$merge_i,
    merge_j = result$merge_j,
    heights = result$heights,
    cluster_sizes = result$cluster_sizes
  ))
}

# ============================================================
#  TRAJECTORY ANALYSIS FUNCTIONS
# ============================================================
#> tox_trajectory_contribution_analysis:compute_baselines_factor_dependent_c: Compute scalar baselines for a factor and dependent variable
#' Compute scalar baselines for a factor and dependent variable
#'
#' @param factor Numeric vector (factor time series)
#' @param dependent Numeric vector (dependent time series)
#' @param mode Baseline computation mode: "raw", "min", or "mean"
#' @return A list with:
#' \describe{
#'   \item{factor_baseline}{Numeric scalar baseline for the factor}
#'   \item{dependent_baseline}{Numeric scalar baseline for the dependent variable}
#' }
tox_compute_baselines_factor_dependent <- function(factor, dependent, mode = "raw") {
  # Input Validation
  validate_numeric_vector(factor)
  validate_numeric_vector(dependent)
  validate_same_length(factor, dependent, "factor", "dependent")
  validate_nonempty(factor)
  validate_string_scalar(mode)
 
 
  # Call the Rcpp forwarder
  result <- tox_compute_baselines_factor_dependent_rcpp(factor,
                                                        dependent,
                                                        mode
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    factor_baseline = result$factor_baseline,
    dependent_baseline = result$dependent_baseline
  ))
}
#> tox_trajectory_contribution_analysis:perform_permutation_test_c: Perform permutation test for trajectory contributions
#' Perform permutation test for trajectory contributions
#'
#' @param trajectories 3D numeric array (factors x samples x timepoints)
#' @param factor_idx Integer index of factor
#' @param dependent_idx Integer index of dependent
#' @param sample_idx Integer index of sample
#' @param mode Baseline mode: "raw", "min", or "mean"
#' @param n_permutations Number of permutations
#' @param random_seed Optional integer seed
#' @return A list with:
#' \describe{
#'   \item{local_contributions}{Numeric matrix (n_timepoints x n_permutations)}
#'   \item{total_contributions}{Numeric vector (length n_permutations)}
#' }
#' @examples
#' res <- tox_perform_permutation_test(trajectories, 1, 2, 1, "raw", 100)
tox_perform_permutation_test <- function(trajectories, factor_idx, dependent_idx, sample_idx, mode, n_permutations, random_seed = NULL) {
  # Input Validation
  validate_numeric_array(trajectories)
  dims <- dim(trajectories)
  validate_length_equals_n(dims, 3)
  n_factors <- dims[1]
  n_samples <- dims[2]
  n_timepoints <- dims[3]

  validate_string_scalar(mode)
 
 
  validate_positive_integer_scalar(factor_idx)
  validate_positive_integer_scalar(dependent_idx)
  validate_positive_integer_scalar(sample_idx)
  validate_index_bounds(factor_idx, low = 1, high = n_factors)
  validate_index_bounds(dependent_idx, low = 1, high = n_factors)
  validate_index_bounds(sample_idx, low = 1, high = n_samples)
  validate_positive_integer_scalar(n_permutations)

  if (is.null(random_seed)) {
    random_seed <- 0L
  }
  validate_is_integer(random_seed, "random_seed")
  validate_length_equals(random_seed, 1, "random_seed")
 
 

  # Call the Rcpp forwarder
  result <- tox_perform_permutation_test_rcpp(trajectories,
                                              n_factors,
                                              n_samples,
                                              n_timepoints,
                                              factor_idx,
                                              dependent_idx,
                                              sample_idx,
                                              mode,
                                              n_permutations,
                                              random_seed
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    local_contributions = result$local_contributions,
    total_contributions = result$total_contributions
  ))
}

#> tox_trajectory_contribution_analysis:compute_p_values_c: Compute p-values for observed vs permutation contributions
#' Compute p-values for observed vs permutation contributions
#'
#' @param local_contributions_observed Numeric vector (n_timepoints)
#' @param total_contribution_observed Numeric scalar
#' @param local_contributions_perm Matrix (n_timepoints x n_permutations)
#' @param total_contributions_perm Numeric vector (n_permutations)
#' @return A list with:
#' \describe{
#'   \item{local_p_values}{Numeric vector of local p-values (length n_timepoints)}
#'   \item{total_p_value}{Numeric scalar total p-value}
#' }
#' @examples
#' res <- tox_compute_p_values(obs, obs_total, perm, perm_total)
tox_compute_p_values <- function(local_contributions_observed, total_contribution_observed, local_contributions_perm, total_contributions_perm) {
  # Input Validation
  validate_numeric_vector(local_contributions_observed)
  validate_numeric_scalar(total_contribution_observed)
  validate_numeric_matrix(local_contributions_perm)
  validate_numeric_vector(total_contributions_perm)

  n_timepoints <- length(local_contributions_observed)
  n_permutations <- ncol(local_contributions_perm)

  validate_length_equals_n(total_contributions_perm, n_permutations)
  validate_length_equals_n(local_contributions_observed, nrow(local_contributions_perm))

  # Call the Rcpp forwarder
  result <- tox_compute_p_values_rcpp(local_contributions_observed,
                                      total_contribution_observed,
                                      local_contributions_perm,
                                      total_contributions_perm,
                                      n_timepoints,
                                      n_permutations
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    local_p_values = result$local_p_values,
    total_p_value = result$total_p_value
  ))
}

#> tox_trajectory_contribution_analysis:compute_contributions_c: Compute contributions for a factor-dependent pair
#' Compute contributions for a factor-dependent pair
#'
#' @param factor Numeric vector (n_timepoints)
#' @param dependent Numeric vector (n_timepoints)
#' @param mode Baseline mode: "raw", "min", or "mean"
#' @return A list with:
#' \describe{
#'   \item{local_contributions}{Numeric vector of local contributions}
#'   \item{total_contribution}{Numeric scalar total contribution}
#' }
#' @examples
#' res <- tox_compute_contributions(factor, dependent, "raw")
tox_compute_contributions <- function(factor, dependent, mode) {
  # Input Validation
  validate_numeric_vector(factor)
  validate_numeric_vector(dependent)
  validate_same_length(factor, dependent, "factor", "dependent")
  validate_nonempty(factor)
  validate_string_scalar(mode)


  # Call the Rcpp forwarder
  result <- tox_compute_contributions_rcpp(factor,
                                           dependent,
                                           mode
  )
    
  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    local_contributions = result$local_contributions,
    total_contribution = result$total_contribution
  ))
}

#> tox_trajectory_contribution_analysis:compute_all_contributions_c: Compute all contributions for selected factor-dependent pairs
#' Compute all contributions for selected factor-dependent pairs
#'
#' @param trajectories 3D numeric array (factors x samples x timepoints)
#' @param factor_indices Integer vector of selected factor indices
#' @param dependent_indices Integer vector of selected dependent indices
#' @param mode Baseline mode: "raw", "min", or "mean"
#' @return A list with:
#' \describe{
#'   \item{local_contributions}{Numeric 4D array (n_timepoints x n_selected_factors x n_selected_dependents x n_samples)}
#'   \item{total_contributions}{Numeric 3D array (n_selected_factors x n_selected_dependents x n_samples)}
#' }
#' @examples
#' res <- tox_compute_all_contributions(trajectories, c(1,2), c(1,2), "raw")
tox_compute_all_contributions <- function(trajectories, factor_indices, dependent_indices, mode) {
  # Input Validation
  validate_numeric_array(trajectories)
  dims <- dim(trajectories)
  validate_length_equals_n(dims, 3)
  n_factors <- dims[1]
  n_samples <- dims[2]
  n_timepoints <- dims[3]
  validate_string_scalar(mode)
  
  
  validate_integer_vector(factor_indices)
  validate_integer_vector(dependent_indices)

  validate_index_bounds(factor_indices, low = 1, high = n_factors)
  validate_index_bounds(dependent_indices, low = 1, high = n_factors)

  n_selected_factors <- length(factor_indices)
  n_selected_dependents <- length(dependent_indices)

  # Call the Rcpp forwarder
  result <- tox_compute_all_contributions_rcpp(trajectories,
                                               n_factors,
                                               n_samples,
                                               n_timepoints,
                                               factor_indices,
                                               n_selected_factors,
                                               dependent_indices,
                                               n_selected_dependents,
                                               mode
  )

  # Check for errors
  check_err_code(result$ierr)
  
  # Reshape to match Fortran output: (n_timepoints, n_selected_factors, n_selected_dependents, n_samples)
  local_contributions <- array(result$local_contributions,
    dim = c(n_timepoints, n_selected_factors, n_selected_dependents, n_samples))
  total_contributions <- array(result$total_contributions,
    dim = c(n_selected_factors, n_selected_dependents, n_samples))
  # Return structured result
  return(list(
    local_contributions = local_contributions,
    total_contributions = total_contributions
  ))
}

#> tox_trajectory_contribution_analysis:compute_velocity_trajectories_c: Compute velocity for all trajectories
#' Compute velocity for all trajectories
#' 
#' @param trajectories 3D numeric array (factors x samples x timepoints)
#' @return Numeric array (n_timepoints-1 x n_factors x n_samples)
tox_compute_velocity_trajectories <- function(trajectories) {
    # Input Validation
    validate_numeric_array(trajectories)
    dims <- dim(trajectories)
    validate_length_equals_n(dims, 3)
    n_factors <- dims[1]
    n_samples <- dims[2]
    n_timepoints <- dims[3]
    # Call the Rcpp forwarder
      result <- tox_compute_velocity_trajectories_rcpp(trajectories,
        n_factors,
        n_samples,
        n_timepoints
      )
    # Check for errors
    check_err_code(result$ierr)

    # Return structured result
    return(array(result$velocity, dim = c(n_timepoints - 1, n_factors, n_samples)))
    
  
  }

#> tox_trajectory_contribution_analysis:compute_acceleration_from_velocity_c: Compute acceleration for all velocity trajectories
#' Compute acceleration for all velocity trajectories
#' 
#' @param velocity 3D numeric array of velocity trajectories (factors x samples x timepoints-1)
#' @return Numeric array (n_timepoints-2 x n_factors x n_samples)
tox_compute_acceleration_from_velocity <- function(velocity) {
    # Input Validation
    validate_numeric_array(velocity)
    dims <- dim(velocity)
    validate_length_equals_n(dims, 3)
    n_vel <- dims[1]
    n_factors <- dims[2]
    n_samples <- dims[3]
    n_timepoints <- n_vel + 1
    # Call the Rcpp forwarder
    result <- tox_compute_acceleration_from_velocity_rcpp(velocity,
      n_factors,
      n_samples,
      n_timepoints
    )
    # Check for errors
    check_err_code(result$ierr)

    # Return structured result
    return(array(result$acceleration, dim = c(n_timepoints - 2, n_factors, n_samples)))
  }

  

#> tox_trajectory_contribution_analysis:compute_velocity_trajectory_c: Compute velocity for a single trajectory
#' Compute velocity for a single trajectory
#' 
#' @param trajectory Numeric vector (length n_timepoints)
#' @return Numeric vector (length n_timepoints-1)
tox_compute_velocity_trajectory <- function(trajectory) {
    # Input Validation
    validate_numeric_vector(trajectory)
    n_timepoints <- length(trajectory)
    # Call the Rcpp forwarder
    result <- tox_compute_velocity_trajectory_rcpp(trajectory, n_timepoints)
    # Return velocity vector directly
    return(result)
  }

#> tox_trajectory_contribution_analysis:compute_acceleration_from_velocity_trajectory_c: Compute acceleration for a single velocity trajectory
#' Compute acceleration for a single velocity trajectory
#' 
#' @param velocity Numeric vector (length n_timepoints-1)
#' @return Numeric vector (length n_timepoints-2)
tox_compute_acceleration_from_velocity_trajectory <- function(velocity) {
    # Input Validation
    validate_numeric_vector(velocity)
    n_timepoints <- length(velocity) + 1
    # Call the Rcpp forwarder
    result <- tox_compute_acceleration_from_velocity_trajectory_rcpp(velocity, n_timepoints)
    # Return acceleration vector directly
    return(result)
  }

#> tox_trajectory_contribution_analysis:compute_velocity_acceleration_contributions_c: Compute velocity and acceleration contributions for a single factor-dependent pair
#' Compute velocity and acceleration contributions for a single factor-dependent pair
#' 
#' @param factor Numeric vector (n_timepoints)
#' @param dependent Numeric vector (n_timepoints)
#' @param mode Baseline mode: "raw", "min", or "mean"
#' @return A list with:
#' \describe{
#'  \item{contrib_velocity}{Numeric scalar total velocity contribution}
#' \item{velocity_contribution_series}{Numeric vector of velocity contributions over time}
#' \item{contrib_acceleration}{Numeric scalar total acceleration contribution}
#' \item{acceleration_contribution_series}{Numeric vector of acceleration contributions over time}
#' }
#' 
tox_compute_velocity_acceleration_contributions <- function(trajectories, mode) {
    # Input Validation
    validate_numeric_array(trajectories)
    dims <- dim(trajectories)
    validate_length_equals_n(dims, 3)
    n_factors <- dims[1]
    n_samples <- dims[2]
    n_timepoints <- dims[3]
    validate_string_scalar(mode)
    
    # Call the Rcpp forwarder
    result <- tox_compute_velocity_acceleration_contributions_rcpp(trajectories,
      n_factors,
      n_samples,
      n_timepoints,
      mode
    )
    # Check for errors
    check_err_code(result$ierr)

    # Return structured result
    return(list(
      contrib_velocity = result$contrib_velocity,
      velocity_contribution_series = result$velocity_contribution_series,
      contrib_acceleration = result$contrib_acceleration,
      acceleration_contribution_series = result$acceleration_contribution_series
    ))
  }

#> tox_trajectory_contribution_analysis:compute_velocity_acceleration_contributions_alloc_c: Compute velocity and acceleration contributions with pre-allocated output
#' Compute velocity and acceleration contributions with pre-allocated output
#' 
#' @param trajectories 3D numeric array (factors x samples x timepoints)
#' @param mode Baseline mode: "raw", "min", or "mean"
#' @return A list with:
#' \describe{
#' \item{contrib_velocity}{Numeric scalar total velocity contribution}
#' \item{velocity_contribution_series}{Numeric vector of velocity contributions over time}
#' \item{contrib_acceleration}{Numeric scalar total acceleration contribution}
#' \item{acceleration_contribution_series}{Numeric vector of acceleration contributions over time}
#' }
tox_compute_velocity_acceleration_contributions_alloc <- function(trajectories, mode) {
    # Input Validation
    validate_numeric_array(trajectories)
    dims <- dim(trajectories)
    validate_length_equals_n(dims, 3)
    n_factors <- dims[1]
    n_samples <- dims[2]
    n_timepoints <- dims[3]
    validate_string_scalar(mode)
    # Call the Rcpp forwarder
      result <- tox_compute_velocity_acceleration_contributions_alloc_rcpp(trajectories,
        n_factors,
        n_samples,
        n_timepoints,
        mode
      )
    # Check for errors
    check_err_code(result$ierr)

    # Return structured result
    return(list(
      contrib_velocity = result$contrib_velocity,
      velocity_contribution_series = result$velocity_contribution_series,
      contrib_acceleration = result$contrib_acceleration,
      acceleration_contribution_series = result$acceleration_contribution_series
    ))
  }

# ============================================================
#  UTILITY FUNCTION
# ============================================================

#> f42_utils:loess_smooth_2d_c: 2D LOESS smoothing for trajectory data
#' @param x_ref Numeric vector of reference x-coordinates (length n_ref)
#' @param y_ref Numeric vector or matrix of reference y-coordinates (length n_ref or n_ref x n_y)
#' @param x_query Numeric vector of query x-coordinates (length n_query)
#' @param indices_used Integer vector of indices (1-based) indicating which reference points to use.
#' @param kernel_sigma Numeric scalar for the Gaussian kernel bandwidth
#' @param kernel_cutoff Numeric scalar for the Gaussian kernel cutoff distance
#' @return A list containing:
#' \describe{
#'   \item{y_out}{Numeric vector of smoothed values at `x_query`}
#'   \item{smoothed_values}{Numeric alias of `y_out` for backward compatibility}
#' }
#' 
tox_loess_smooth_2d <- function(x_ref, y_ref, x_query, indices_used = NULL, kernel_sigma, kernel_cutoff) {
  # Input Validation
  validate_numeric_vector(x_ref, "x_ref")
  validate_nonempty(x_ref, "x_ref")

  if (is.matrix(y_ref)) {
    validate_numeric_matrix(y_ref, "y_ref")
    y_ref <- as.vector(y_ref)
  } else {
    validate_numeric_vector(y_ref, "y_ref")
  }
  validate_length_equals_n(y_ref, length(x_ref), "y_ref")

  validate_numeric_vector(x_query, "x_query")
  validate_nonempty(x_query, "x_query")

  if (!is.null(indices_used)) {
    validate_logical_or_index_vector(indices_used, name = "indices_used")
    validate_nonempty(indices_used, "indices_used")
  }

  validate_numeric_scalar(kernel_sigma, "kernel_sigma")
  validate_numeric_scalar(kernel_cutoff, "kernel_cutoff")

  # Indices used: default to all, then validate
  n_total <- length(x_ref)
  if (is.null(indices_used)) {
    indices_used <- seq_len(n_total)
  }
  
 

  

  # Call the Rcpp forwarder
  result <- tox_loess_smooth_2d_rcpp(n_total = n_total,
                                     n_target = length(x_query),
                                     x_ref = x_ref,
                                     y_ref = y_ref,
                                     indices_used = indices_used,
                                     n_used = length(indices_used),
                                     x_query = x_query,
                                     kernel_sigma = kernel_sigma,
                                     kernel_cutoff = kernel_cutoff
  )

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    y_out = result$y_out,
    smoothed_values = result$smoothed_values
  ))
}

#> f42_utils:compute_edf_c: Compute Empirical Distribution Function (EDF)
#' Compute Empirical Distribution Function (EDF)
#' @param values Numeric vector of observed values
#' @param perm Optional integer vector of permutation indices (sorted by values[perm])
#' @return A list with:
#' \describe{
#'   \item{unique_values}{Numeric vector of unique sorted values}
#'   \item{cdf_values}{Numeric vector of cumulative distribution values}
#'   \item{n_unique}{Integer number of unique values}
#' }
tox_compute_edf <- function(values, perm = NULL) {
  # Input Validation
  validate_numeric_vector(values)

  if (is.null(perm)) {
    # Call the Rcpp forwarder (auto permutation)
    result <- tox_compute_edf_rcpp(values)
  } else {
    # Call the Rcpp forwarder (expert permutation)

    validate_is_permutation(perm, n = length(values), name = "perm")
    result <- tox_compute_edf_expert_rcpp(values, perm)
  }

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    unique_values = result$unique_values,
    cdf_values = result$cdf_values,
    n_unique = as.integer(result$n_unique)
  ))
}

#> f42_utils:compute_edf_expert_c: Expert EDF with pre-sorted permutation
#' Compute Empirical Distribution Function (EDF) using expert permutation input
#' @param values Numeric vector of observed values
#' @param perm Integer vector of permutation indices (sorted by values[perm])
#' @return A list with:
#' \describe{
#'   \item{unique_values}{Numeric vector of unique sorted values}
#'   \item{cdf_values}{Numeric vector of cumulative distribution values}
#'   \item{n_unique}{Integer number of unique values}
#' }
tox_compute_edf_expert <- function(values, perm) {
  # Input Validation
  validate_numeric_vector(values)
  validate_is_permutation(perm, n = length(values), name = "perm")

  # Call the Rcpp forwarder
  result <- tox_compute_edf_expert_rcpp(values, perm)

  # Check for errors
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    unique_values = result$unique_values,
    cdf_values = result$cdf_values,
    n_unique = as.integer(result$n_unique)
  ))
}

#> f42_utils:which_c: Get indices of TRUE or non-zero elements in a mask
#' This function takes a logical or integer mask and returns the indices of the TRUE or non-zero elements, up to a specified maximum.
#' @param mask A logical or integer vector serving as the mask for which to find indices.
#' @param m_max An integer scalar specifying the maximum number of indices to return (default: length of mask).
#' @return An integer vector of indices corresponding to TRUE or non-zero elements in `mask`.
#' 
tox_which <- function(mask, m_max = length(mask)) {
  # Input Validation

  if (is.logical(mask)) {
    mask <- as.integer(mask)
  }
  validate_integer_vector(mask)
  validate_positive_integer_scalar(m_max)

  # Call the Rcpp forwarder
  result <- tox_which_rcpp(mask = mask, m_max = m_max)

  # Check for errors
  check_err_code(result$ierr)

  # Return indices of TRUE/non-zero elements (capped by m_max)
  m_out <- as.integer(result$m_out)
  n_take <- min(m_out, m_max, length(result$idx_out))
  if (n_take <= 0L) {
    return(integer(0))
  }

  return(result$idx_out[seq_len(n_take)])
}