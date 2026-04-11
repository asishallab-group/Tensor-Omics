
library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library

lib_path <- shQuote(normalizePath("build"))

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)


cat("✓ TensorOmics Rcpp functions loaded successfully\n")

source("rcpp/error_handling.R")



# Helper functions for string to raw conversions
strings_to_raw_matrix <- function(arr, clen) {
  n <- length(arr)
  # Create a matrix of raw bytes with dimensions clen x n
  mat <- matrix(raw(1), nrow = clen, ncol = n)
  for (i in seq_along(arr)) {
    # Convert string to raw bytes
    raw_bytes <- charToRaw(arr[i])
    length_bytes <- length(raw_bytes)
    if (length_bytes > 0) {
      # Copy raw bytes into matrix
      mat[1:min(length_bytes, clen), i] <- raw_bytes[1:min(length_bytes, clen)]
    }
    # Note: No null termination needed - Fortran handles this
  }
  mat  # Return the raw matrix
}

# Helper function for raw to string conversion
raw_matrix_to_strings <- function(raw_mat, clen) {
  # raw_mat: raw matrix with dimensions clen x n
  if (is.null(dim(raw_mat))) {
    # 1D vector: convert to matrix with one column
    raw_mat <- matrix(raw_mat, nrow = clen)
  }
  n <- ncol(raw_mat)
  strings <- character(n)
  for (i in seq_len(n)) {
    raw_vec <- raw_mat[, i]
    # Find first null byte (if any)
    null_pos <- which(raw_vec == as.raw(0))
    if (length(null_pos) > 0) {
      end_pos <- min(null_pos[1] - 1, length(raw_vec))
    } else {
      end_pos <- length(raw_vec)
    }
    if (end_pos > 0) {
      strings[i] <- rawToChar(raw_vec[1:end_pos])
    } else {
      strings[i] <- ""
    }
  }
  strings
}

# Helper function for string to raw conversion
string_to_raw <- function(s, len) {
  raw_bytes <- charToRaw(s)
  if (length(raw_bytes) < len) {
    # Pad with zeros if needed (Fortran will handle null termination)
    raw_bytes <- c(raw_bytes, raw(1))
    length(raw_bytes) <- len
    raw_bytes[is.na(raw_bytes)] <- as.raw(0)
  }
  raw_bytes[1:len]  # Ensure exact length
}

#' Read gene IDs from a TSV file (R wrapper)
#' @param filename Path to TSV file
#' @param n_genes Number of genes (rows)
#' @param gene_ids_len Max length of gene ID string
#' @param n_header_rows Number of header rows
#' @param gene_col Column index for gene IDs
#' @return List with gene_ids (character vector) and ierr
read_gene_ids_from_tsv_file <- function(filename, n_genes, gene_ids_len, n_header_rows, gene_col) {
  filename_raw <- charToRaw(filename)
  res <- tox_read_gene_ids_from_tsv_file_rcpp(filename_raw, n_genes, gene_ids_len, n_header_rows, gene_col)
  gene_ids <- raw_matrix_to_strings(res$gene_ids_raw, gene_ids_len)
  check_err_code(res$ierr)
  list(gene_ids = gene_ids, ierr = res$ierr)
}

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
  file_list_raw <- strings_to_raw_matrix(file_list, max(nchar(file_list)))
  gene_ids_raw <- strings_to_raw_matrix(gene_ids, max(nchar(gene_ids)))
  delimiter_raw <- charToRaw(delimiter)
  res <- tox_read_expression_vectors_tsv_rcpp(
    file_list_raw,
    gene_ids_raw,
    as.integer(value_cols),
    delimiter_raw,
    as.integer(n_samples),
    as.integer(n_header_rows),
    as.integer(gene_col)
  )
  check_err_code(res$ierr)
  list(expression_vectors = res$expression_vectors, ierr = res$ierr)
}

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
  list(
    family_ids = raw_matrix_to_strings(result$family_ids_raw, family_raw_len),
    gene_to_fam = result$gene_to_fam,
    ierr = result$ierr
  )
}

#' Filter out genes without family assignments
#' @param gene_ids Character vector of gene IDs
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @return A list with elements:    
#'   - gene_ids: Filtered character vector of gene IDs
#'   - expression_vectors: Filtered numeric matrix of expression values
#'   - gene_to_fam: Filtered integer vector mapping each gene to its family index
filter_unassigned_genes <- function(gene_ids, expression_vectors, gene_to_fam) {
  
  mask <- gene_to_fam != 0L
  filtered_gene_ids <- gene_ids[mask]
  filtered_expression_vectors <- expression_vectors[, mask, drop = FALSE]
  filtered_gene_to_fam <- gene_to_fam[mask]
  n_genes_kept <- sum(mask)
  list(
    gene_ids = filtered_gene_ids,
    expression_vectors = filtered_expression_vectors,
    gene_to_fam = filtered_gene_to_fam,
    n_genes_kept = n_genes_kept
  )
}


# R wrappers for validation routines (via C bindings)
# Uses raw conversion helpers: strings_to_raw_matrix, raw_matrix_to_strings

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
  list(ierr = ierr)
}

#' Validate gene to family mapping
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param n_genes Number of genes
#' @param n_families Number of gene families
validate_gene_to_family_mapping <- function(gene_to_fam, n_genes, n_families) {
  ierr <- tox_validate_gene_to_family_mapping_rcpp(gene_to_fam, n_genes, n_families)
  check_err_code(ierr)
  list(ierr = ierr)
}

#' Validate expression data
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param check_non_negative Logical flag to check for non-negative values
validate_expression_data <- function(expression_vectors, n_genes, n_samples, check_non_negative = TRUE) {
  expr_matrix <- as.matrix(expression_vectors)

  ierr <- tox_validate_expression_data_rcpp(expression_vectors, 
  n_genes, n_samples, check_non_negative)

  check_err_code(ierr)
  list(ierr = ierr)
}

#' Validate family centroids
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param n_families Number of gene families
#' @param n_samples Number of samples
validate_family_centroids <- function(family_centroids, n_families, n_samples) {
  fam_matrix <- as.matrix(family_centroids)
  ierr <- tox_validate_family_centroids_rcpp(
    family_centroids,
    n_families,
    n_samples
  )
  check_err_code(ierr)
  list(ierr = ierr)
}

#' Validate shift vectors
#' @param shift_vectors Numeric matrix of shift vectors (2*n_samples x n_genes)
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param n_samples Number of samples
validate_shift_vectors <- function(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_samples) {
  expr_matrix <- as.matrix(expression_vectors)
  fam_matrix <- as.matrix(family_centroids)
  shift_matrix <- as.matrix(shift_vectors)

  n_genes <- ncol(expr_matrix)
  n_families <- ncol(fam_matrix)

  ierr <- tox_validate_shift_vectors_rcpp(shift_vectors, 
  expression_vectors, family_centroids, gene_to_fam, n_samples)
  check_err_code(ierr)
  list(ierr = ierr)
}

#' Validate uniqueness of string array
#' @param string_arr Character vector of gene IDs
#' @param n_strings Number of genes
#' Note: Uses hashset internally which may increase memory usage temporarily for large datasets
validate_string_array_uniqueness <- function(string_arr, n_strings) {
  n_strings <- as.integer(n_strings)

  ierr <- tox_validate_string_array_uniqueness_rcpp(string_arr, n_strings)
  check_err_code(ierr)
  list(ierr = ierr)
}

#' Validate all data components together
#' @param n_genes Number of genes
#' @param n_families Number of gene families
#' @param n_samples Number of samples
#' @param gene_ids Character vector of gene IDs
#' @param gene_family_ids Character vector of gene family IDs
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param family_centroids Numeric matrix of family centroids (n_samples x n_families)
#' @param shift_vectors Numeric matrix of shift vectors (2*n_samples x n_genes)
validate_all_data <- function(n_genes, n_families, n_samples,
                              gene_ids, gene_family_ids,
                              gene_to_fam, expression_vectors,
                              family_centroids, shift_vectors) {
  n_genes    <- as.integer(n_genes)
  n_families <- as.integer(n_families)
  n_samples  <- as.integer(n_samples)

  expr_matrix <- as.matrix(expression_vectors)
  fam_matrix <- as.matrix(family_centroids)
  shift_matrix <- as.matrix(shift_vectors)

  ierr <- tox_validate_all_data_rcpp(n_genes, n_families, n_samples,
                              gene_ids, gene_family_ids,
                              gene_to_fam, expression_vectors,
                              family_centroids, shift_vectors)

  check_err_code(ierr)
  list(ierr = ierr)
}


create_zip_archive <- function(zip_filename, keys, filenames) {
  stop("Zip archive helpers have been removed. Use an external zip tool instead.")
}

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
save_tox_data <- function(zip_filename,
                                 gene_ids = NULL, gene_ids_name = NULL,
                                 expression_vectors = NULL, expression_vectors_name = NULL,
                                 gene_to_fam = NULL, gene_to_fam_name = NULL,
                                 family_ids = NULL, family_ids_name = NULL,
                                 family_centroids = NULL, family_centroids_name = NULL,
                                 shift_vectors = NULL, shift_vectors_name = NULL) {
  
 
  # Validation moved to error_handling.R
  validate_non_empty_string(zip_filename)
  validate_character_vector(gene_ids)
  validate_numeric_matrix(expression_vectors)
  validate_integer_vector(gene_to_fam)
  validate_character_vector(family_ids)
  validate_numeric_matrix(family_centroids)
  validate_numeric_matrix(shift_vectors)
    
    # Write files to temporary directory
  temp_files <- character(0)
  keys <- character(0)
  filenames <- character(0)
  
  if (gene_ids_array_valid && gene_ids_name_valid) {
    tox_serialize_char_array(gene_ids, gene_ids_name)
    if(debug) {message(paste("Wrote gene IDs to", gene_ids_name))}
    temp_files <- c(temp_files, gene_ids_name)
    keys <- c(keys, "gene_ids")
    filenames <- c(filenames, gene_ids_name)
  }

  if (expression_vectors_array_valid && expression_vectors_name_valid) {
    tox_serialize_real_array(expression_vectors, expression_vectors_name)
    if(debug) {message(paste("Wrote expression vectors to", expression_vectors_name))}
    temp_files <- c(temp_files, expression_vectors_name)
    keys <- c(keys, "expression")
    filenames <- c(filenames, expression_vectors_name)
  }
  
  if (gene_to_fam_array_valid && gene_to_fam_name_valid) {
    tox_serialize_int_array(gene_to_fam, gene_to_fam_name)
    if(debug) {message(paste("Wrote gene to family to", gene_to_fam_name))}
    temp_files <- c(temp_files, gene_to_fam_name)
    keys <- c(keys, "gene_to_family")
    filenames <- c(filenames, gene_to_fam_name)
  }
  
  if (family_ids_array_valid && family_ids_name_valid) {
    tox_serialize_char_array(family_ids, family_ids_name)
    if(debug) {message(paste("Wrote family IDs to", family_ids_name))}
    temp_files <- c(temp_files, family_ids_name)
    keys <- c(keys, "family_ids")
    filenames <- c(filenames, family_ids_name)
  }
  
  if (family_centroids_array_valid && family_centroids_name_valid) {
    tox_serialize_real_array(family_centroids, family_centroids_name)
    if(debug){ message(paste("Wrote family centroids to", family_centroids_name))}
    temp_files <- c(temp_files, family_centroids_name)
    keys <- c(keys, "family_centroids")
    filenames <- c(filenames, family_centroids_name)
  }
  
  if (shift_vectors_array_valid && shift_vectors_name_valid) {
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
  stop("Zip archive helpers have been removed. Use an external zip tool instead.")
}

extract_zip_archive <- function(zip_filename) {
  stop("Zip archive helpers have been removed. Use an external zip tool instead.")
}



