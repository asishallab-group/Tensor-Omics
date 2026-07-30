
library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library

lib_path <- shQuote(normalizePath("build"))

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv, cacheDir = "rcpp/rcpp_cache")
tox_validate_gene_to_family_mapping_rcpp

cat("✓ TensorOmics Rcpp Tox Data functions loaded successfully\n")

source("rcpp/error_handling.R")

#> tox_data_read_write:read_gene_ids_from_tsv_file_c: Read gene ids from a tsv file
#' Read gene IDs from a TSV file (R wrapper)
#' @param filename Path to TSV file
#' @param n_genes Number of genes (rows)
#' @param gene_ids_len Max length of gene ID string
#' @param n_header_rows Number of header rows
#' @param gene_col Column index for gene IDs
#' @return List with gene_ids (character vector) and ierr
read_gene_ids_from_tsv_file <- function(filename, n_genes, gene_ids_len, n_header_rows, gene_col) {
  res <- tox_read_gene_ids_from_tsv_file_rcpp(filename, n_genes, gene_ids_len, n_header_rows, gene_col)
  check_err_code(res$ierr)
  return(res)
}

#> tox_data_read_write:read_expression_vectors_tsv_c: Read expression vectors from given tabular (csv/tsv) files
#' Read expression vectors from TSV file(s) (R wrapper)
#' @param file_list Character vector of file paths (length 1 for single file)
#' @param gene_ids Character vector of gene IDs (order to match)
#' @param n_samples Number of samples (rows)
#' @param n_header_rows Number of header rows
#' @param gene_col Column index for gene IDs
#' @param value_cols Integer vector of value column indices
#' @param delimiter Delimiter (default tab)
#' @return List with expression_vectors (matrix) and ierr
read_expression_vectors_tsv <- function(file_list, gene_ids, n_samples, n_header_rows, gene_col, value_cols, delimiter = "\t") {
  # Convert file_list and gene_ids to raw matrices
  res <- tox_read_expression_vectors_tsv_rcpp(
    file_list,
    gene_ids,
    as.integer(value_cols),
    delimiter,
    as.integer(n_samples),
    as.integer(n_header_rows),
    as.integer(gene_col)
  )
  check_err_code(res$ierr)
  return(res)
}

#> tox_data_read_write:read_orthofinder_file_c: Read an orthofinder family file and map genes to families
#' Read gene family assignments from a file
#' @param filename Character string of the filename
#' @param gene_ids Character vector of gene IDs to match
#' @param n_families Number of unique gene families expected
#' @param family_len Maximum length of each family ID
#' @return A list with elements:
#'   - family_ids: Character vector of family IDs read from the file
#'   - gene_to_fam: Integer vector mapping each gene to its family index
#'   - ierr: Integer error code (0 if successful)
#' Note: If genes are not found in the gene IDs list a message will be printed but NO error code is thrown
read_orthofinder_file <- function(filename, gene_ids, n_families, family_len) {
  n_families <- as.integer(n_families)
  family_len <- as.integer(family_len)


  result <- tox_read_orthofinder_file_rcpp(filename, gene_ids, n_families, family_len)

  check_err_code(result$ierr)
  return(result)
}

#> tox_data_tools:filter_unassigned_genes_c: Filter out genes that are not assigned to any family (where gene_to_fam == 0).
#' Filter out genes without family assignments
#' @param gene_ids Character vector of gene IDs
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @return A list with elements:    
#'   - mask: Logical Mask indicating unassigned genes
#'   - n_genes_kept: Count of TRUE in mask
filter_unassigned_genes <- function(gene_ids, gene_to_fam) {
  
  res <- tox_filter_unassigned_genes_rcpp(gene_ids, gene_to_fam)
  res$mask <- res$mask != 0
  check_err_code(res$ierr)
  return(res)
}

#> tox_data_validation:validate_data_structure_c: Validate overall data structure consistency. Confirms sizes and dependencies as far as possible.
#' Validate overall data structure
#' @param n_genes Number of genes
#' @param n_families Number of gene families
#' @param n_samples Number of samples
#' @param d Dimensionality of expression vectors (should be 2 * n_samples)
#' @param gene_ids Character vector of gene IDs
#' @param gene_family_ids Character vector of gene family IDs
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param shift_vectors Numeric matrix of shift vectors (2*n_samples x n_genes)
validate_data_structure <- function(n_genes, n_families, n_samples, d,
                                    gene_ids, gene_family_ids,
                                    gene_to_fam, expression_vectors,
                                    family_centroids, shift_vectors) {
  
  n_genes    <- as.integer(n_genes)
  n_families <- as.integer(n_families)
  n_samples  <- as.integer(n_samples)
  d          <- as.integer(d)

  expr_matrix <- as.matrix(expression_vectors)
  fam_matrix <- as.matrix(family_centroids)
  shift_matrix <- as.matrix(shift_vectors)

  ierr <- tox_validate_data_structure_rcpp(n_genes, n_families, n_samples, d,
                                    gene_ids, gene_family_ids,
                                    gene_to_fam, expression_vectors,
                                    family_centroids, shift_vectors)

  check_err_code(ierr)
  return(ierr)
}

#> tox_data_validation:validate_gene_to_family_mapping_c: Validate gene to family mapping
#' Validate gene to family mapping
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param n_families Number of families
validate_gene_to_family_mapping <- function(gene_to_fam, n_families) {
  ierr <- tox_validate_gene_to_family_mapping_rcpp(gene_to_fam, n_families)
  check_err_code(ierr)
  return(ierr)
}

#> tox_data_validation:validate_expression_data_c: Validate expression data
#' Validate expression data
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param check_non_negative Logical flag to check for non-negative values
validate_expression_data <- function(expression_vectors, check_non_negative = TRUE) {
  expr_matrix <- as.matrix(expression_vectors)

  ierr <- tox_validate_expression_data_rcpp(expression_vectors, check_non_negative)

  check_err_code(ierr)
  return(ierr)
}

#> tox_data_validation:validate_family_centroids_c: Validate family centroids, checks for NaN/Inf
#' Validate family centroids
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
validate_family_centroids <- function(family_centroids) {
  fam_matrix <- as.matrix(family_centroids)
  ierr <- tox_validate_family_centroids_rcpp(family_centroids)
  check_err_code(ierr)
  return(ierr)
}

#> tox_data_validation:validate_shift_vectors_c: Validate shift vectors, checks if datatypes are correct and if the general structure matches
#' Validate shift vectors
#' @param shift_vectors Numeric matrix of shift vectors (2*n_samples x n_genes)
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
validate_shift_vectors <- function(shift_vectors, expression_vectors, family_centroids, gene_to_fam) {
  expr_matrix <- as.matrix(expression_vectors)
  fam_matrix <- as.matrix(family_centroids)
  shift_matrix <- as.matrix(shift_vectors)

  n_genes <- ncol(expr_matrix)
  n_families <- ncol(fam_matrix)

  ierr <- tox_validate_shift_vectors_rcpp(shift_vectors, expression_vectors, family_centroids, gene_to_fam)
  check_err_code(ierr)
  return(ierr)
}

#> tox_data_validation:validate_string_array_uniqueness_c: Validate uniqueness of strings
#' Validate uniqueness of string array
#' @param string_arr Character vector of gene IDs
#' Note: Uses hashset internally which may increase memory usage temporarily for large datasets
validate_string_array_uniqueness <- function(string_arr) {
  ierr <- tox_validate_string_array_uniqueness_rcpp(string_arr)
  check_err_code(ierr)
  return(ierr)
}

#> tox_data_validation:validate_all_data_c: Comprehensive validation of all data components. This function performs all individual validations in one go.
#' Validate all data components together
#' @param gene_ids Character vector of gene IDs
#' @param gene_family_ids Character vector of gene family IDs
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param shift_vectors Numeric matrix of shift vectors (2*n_samples x n_genes)
validate_all_data <- function(gene_ids, gene_family_ids,
                              gene_to_fam, expression_vectors,
                              family_centroids, shift_vectors) {
  expr_matrix <- as.matrix(expression_vectors)
  fam_matrix <- as.matrix(family_centroids)
  shift_matrix <- as.matrix(shift_vectors)

  ierr <- tox_validate_all_data_rcpp(gene_ids, gene_family_ids,
                              gene_to_fam, expression_vectors,
                              family_centroids, shift_vectors)

  check_err_code(ierr)
  return(ierr)
}


#> tox_data_archive:create_zip_archive_c: Low-level function to create zip archive from keys and filenames.
create_zip_archive <- function(zip_filename, keys, filenames) {
  # stop("Zip archive helpers have been removed. Use an external zip tool instead.")
}

#> f42_helper: Save tox data to zip archive
#' Save standard conform tox data directly to zip archive
#' @param zip_filename Name of the zip file to create
#' @param gene_ids Character vector of gene IDs
#' @param gene_ids_name Filename for gene IDs in the archive
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param expression_vectors_name Filename for expression vectors in the archive
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param gene_to_fam_name Filename for gene to family mapping in the archive
#' @param family_ids Character vector of family IDs
#' @param family_ids_name Filename for family IDs in the archive
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param family_centroids_name Filename for family centroids in the archive
#' @param shift_vectors Numeric matrix of shift vectors (2*n_samples x n_genes)
#' @param shift_vectors_name Filename for shift vectors in the archive
#' @param debug Specifying whether to print status messages or not
save_tox_data <- function(zip_filename,
                                 gene_ids = NULL, gene_ids_name = NULL,
                                 expression_vectors = NULL, expression_vectors_name = NULL,
                                 gene_to_fam = NULL, gene_to_fam_name = NULL,
                                 family_ids = NULL, family_ids_name = NULL,
                                 family_centroids = NULL, family_centroids_name = NULL,
                                 shift_vectors = NULL, shift_vectors_name = NULL, debug = TRUE) {
  
 
  # Validation moved to error_handling.R
  validate_non_empty_string(zip_filename)
    
    # Write files to temporary directory
  temp_files <- character(0)
  keys <- character(0)
  filenames <- character(0)
  
  if (!is.null(gene_ids)) {
    validate_character_vector(gene_ids)
    tox_serialize_char_array(gene_ids, gene_ids_name)
    if(debug) {message(paste("Wrote gene IDs to", gene_ids_name))}
    temp_files <- c(temp_files, gene_ids_name)
    keys <- c(keys, "gene_ids")
    filenames <- c(filenames, gene_ids_name)
  }

  if (!is.null(expression_vectors)) {
    validate_numeric_matrix(expression_vectors)
    tox_serialize_real_array(expression_vectors, expression_vectors_name)
    if(debug) {message(paste("Wrote expression vectors to", expression_vectors_name))}
    temp_files <- c(temp_files, expression_vectors_name)
    keys <- c(keys, "expression")
    filenames <- c(filenames, expression_vectors_name)
  }
  
  if (!is.null(gene_to_fam)) {
    validate_integer_vector(gene_to_fam)
    tox_serialize_int_array(gene_to_fam, gene_to_fam_name)
    if(debug) {message(paste("Wrote gene to family to", gene_to_fam_name))}
    temp_files <- c(temp_files, gene_to_fam_name)
    keys <- c(keys, "gene_to_family")
    filenames <- c(filenames, gene_to_fam_name)
  }
  
  if (!is.null(family_ids)) {
    validate_character_vector(family_ids)
    tox_serialize_char_array(family_ids, family_ids_name)
    if(debug) {message(paste("Wrote family IDs to", family_ids_name))}
    temp_files <- c(temp_files, family_ids_name)
    keys <- c(keys, "family_ids")
    filenames <- c(filenames, family_ids_name)
  }
  
  if (!is.null(family_centroids)) {
    validate_numeric_matrix(family_centroids)
    tox_serialize_real_array(family_centroids, family_centroids_name)
    if(debug){ message(paste("Wrote family centroids to", family_centroids_name))}
    temp_files <- c(temp_files, family_centroids_name)
    keys <- c(keys, "family_centroids")
    filenames <- c(filenames, family_centroids_name)
  }
  
  if (!is.null(shift_vectors)) {
    validate_numeric_matrix(shift_vectors)
    tox_serialize_real_array(shift_vectors, shift_vectors_name)
    if(debug) { message(paste("Wrote shift vectors to", shift_vectors_name))}
    temp_files <- c(temp_files, shift_vectors_name)
    keys <- c(keys, "shift_vectors")
    filenames <- c(filenames, shift_vectors_name)
  }

  # Call the low-level function
  if (length(keys) > 0) {
    create_zip_archive(zip_filename, keys, filenames)
  } else {
    warning("No valid data provided to save - skipping archive creation")
  }

  # Clean up temporary files
  for (temp_file in temp_files) {
    if (file.exists(temp_file)) {
      file.remove(temp_file)
    }
  }
}



#> f42_helper: Read tox data from zip archive
#' Load standard conform tox data directly from zip archive
#' @param zip_filename Name of the zip file to read from
#' @param gene_ids If not NULL, will attempt to read gene IDs from archive
#' @param expression_vectors If not NULL, will attempt to read expression vectors from archive
#' @param gene_to_fam If not NULL, will attempt to read gene to family mapping from archive
#' @param family_ids If not NULL, will attempt to read family IDs from archive
#' @param family_centroids If not NULL, will attempt to read family centroids from archive
#' @param shift_vectors If not NULL, will attempt to read shift vectors from archive
read_tox_data <- function(zip_filename,
                          gene_ids = NULL,
                          expression_vectors = NULL,
                          gene_to_fam = NULL,
                          family_ids = NULL,
                          family_centroids = NULL,
                          shift_vectors = NULL) {
  # stop("Zip archive helpers have been removed. Use an external zip tool instead.")
}

#> tox_data_archive:extract_zip_archive_c: Extract a zip archive created by create_zip_archive
extract_zip_archive <- function(zip_filename) {
  # stop("Zip archive helpers have been removed. Use an external zip tool instead.")
}



