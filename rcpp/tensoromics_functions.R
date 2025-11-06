library(Rcpp)

# Resolve script directory so paths work regardless of working directory
get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", cmd, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  of <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of)))
  # fallback to cwd
  return(getwd())
}
script_dir <- get_script_dir()

# Get absolute path to build directory if it exists
build_dir <- file.path(script_dir, "build")
if (!file.exists(build_dir)) build_dir <- file.path(dirname(script_dir), "build")
if (file.exists(build_dir)) {
  lib_path <- normalizePath(build_dir)
  Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))
} else {
  warning("build/ directory not found; skipping PKG_LIBS setting. If Fortran-linked routines fail, ensure build/ contains the lib.")
}

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
# Try to source the cpp file relative to script_dir or the project rcpp/ folder
cpp_candidates <- c(
  file.path(script_dir, "tensoromics_functions.cpp"),
  file.path(script_dir, "rcpp", "tensoromics_functions.cpp"),
  file.path(dirname(script_dir), "rcpp", "tensoromics_functions.cpp")
)
cpp_path <- NULL
for (p in cpp_candidates) if (file.exists(p)) { cpp_path <- p; break }
if (is.null(cpp_path)) stop("Could not find tensoromics_functions.cpp in expected locations. Tried:\n  ", paste(cpp_candidates, collapse = "\n  "))
sourceCpp(cpp_path, env = .GlobalEnv)

cat("✓ TensorOmics Rcpp functions loaded (sourceCpp ok)\n")

# Try to source error_handling.R from common locations
eh_candidates <- c(
  file.path(script_dir, "..", "r", "error_handling.R"),
  file.path(script_dir, "..", "rcpp", "error_handling.R"),
  file.path(script_dir, "error_handling.R"),
  file.path(script_dir, "r", "error_handling.R")
)
eh_found <- NULL
for (p in eh_candidates) if (file.exists(p)) { eh_found <- p; break }
if (!is.null(eh_found)) {
  source(eh_found)
} else {
  warning("error_handling.R not found; check error-handling helpers if you encounter failures.")
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
  genes <- as.numeric(genes)
  centroids <- as.numeric(centroids)
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
  
  # Call Rcpp wrapper
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
  
  # Convert to appropriate types for Rcpp
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
  
  # Validate dimensions
  if (length(vector_selection) != ncol(expression_vectors)) {
    stop("vector_selection length must match number of columns in expression_vectors")
  }
  if (length(axis_selection) != nrow(expression_vectors)) {
    stop("axis_selection length must match number of rows in expression_vectors")
  }
  
  # Call Rcpp wrapper
  result <- tox_calculate_tissue_versatility_rcpp(expression_vectors, vector_selection, axis_selection)
  
  result$tissue_versatilities <- pmin(1, pmax(0, result$tissue_versatilities))

  # Check for errors
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  # Return structured result 
  return(list(
    tissue_versatilities = result$tissue_versatilities,
    tissue_angles_deg = result$tissue_angles_deg,
    n_selected_vectors = result$n_selected_vectors,
    n_selected_axes = result$n_selected_axes
  ))
}

# ===================================================================
# NORMALIZATION FUNCTIONS
# ===================================================================

tox_normalize_by_std_dev <- function(input) {
  if (!is.matrix(input)) 
    stop("Input must be a matrix")

  result <- normalize_by_std_dev_rcpp(input)

  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  return(result$output)
}


tox_quantile_normalization <- function(input, max_stack = 10000) {

  # Validate input type
  if (!is.matrix(input)) 
    stop("Input must be a matrix")
  if (!is.numeric(max_stack)) 
    stop("max_stack must be numeric")

  # Call the Rcpp wrapper (Fortran runs underneath)
  result <- quantile_normalization_rcpp(input, as.integer(max_stack))

  # Handle errors using the centralized error handling system
  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  # Return only the normalized matrix (not the full result list)
  return(result$output)
}


tox_normalize_data <- function(input, group_s, group_c, max_stack = 10000) {

  # Validate input type
  if (!is.matrix(input))
    stop("Input must be a matrix")

  # Ensure grouping vectors are integer
  if (!is.integer(group_s)) group_s <- as.integer(group_s)
  if (!is.integer(group_c)) group_c <- as.integer(group_c)

  # Call the Rcpp wrapper
  result <- tox_normalize_data_rcpp(input, group_s, group_c, as.integer(max_stack))

  # Handle Fortran error codes
  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  # Return only processed buffers (not the error code)
  return(list(
    buf_stddev = result$buf_stddev,
    buf_quant  = result$buf_quant,
    buf_avg    = result$buf_avg,
    buf_log    = result$buf_log
  ))
}

tox_log2_transformation <- function(input) {

  # Validate input
  if (!is.matrix(input))
    stop("Input must be a matrix")

  # Call the Rcpp wrapper
  result <- log2_transformation_rcpp(input)

  # Check error code
  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  # Return transformed output
  return(result$output)
}


tox_calc_tiss_avg <- function(input, group_s, group_c) {

  # Validate input
  if (!is.matrix(input))
    stop("Input must be a matrix")

  # Call Rcpp wrapper
  result <- calc_tiss_avg_rcpp(input, as.integer(group_s), as.integer(group_c))

  # Handle Fortran errors
  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  return(result$output)
}

tox_calc_fchange <- function(input, control_cols, cond_cols) {

  # Validate input
  if (!is.matrix(input))
    stop("Input must be a matrix")

  # Call Rcpp wrapper
  result <- calc_fchange_rcpp(input, as.integer(control_cols), as.integer(cond_cols))

  # Check for errors
  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  return(result$output)
}
#####
cat("✓ Added normalization functions successfully\n")

